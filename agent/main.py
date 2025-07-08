import os
from agent.search import fetch_pubmed_articles, fetch_arxiv_articles, save_raw_results
from agent.summarize_and_prepare_email import main as prepare_email_html
from agent import send_digest_email

import re

def extract_lead_author(article):
    authors = article.get("authors", [])
    return authors[0] if authors else "Unknown"

def extract_email(text_blob):
    if not isinstance(text_blob, str):
        return ""
    matches = re.findall(r'\b[\w.-]+?@\w+?\.\w+?\b', text_blob)
    return matches[0] if matches else ""

def main():
    print("Fetching articles...")

    query = "echocardiography AND (deep learning OR artificial intelligence)"
    pubmed_articles = fetch_pubmed_articles(query, max_results=10)
    arxiv_articles = fetch_arxiv_articles(query, max_results=10)

    save_raw_results(pubmed_articles, arxiv_articles)

    all_articles = pubmed_articles + arxiv_articles
    print(f"Total articles fetched: {len(all_articles)}")

    # Enrich articles with lead author and email
    for article in all_articles:
        article["lead_author"] = extract_lead_author(article)
        # Try pulling email from summary or anywhere
        email_source = article.get("summary", "") + " " + article.get("title", "")
        article["author_email"] = extract_email(email_source)

    # Generate and send email
    print("Summarizing and sending email...")
    html_content = prepare_email_html()
    send_digest_email.send_email(html_content)

if __name__ == "__main__":
    main()
