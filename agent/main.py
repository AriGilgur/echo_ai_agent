import os
import sys
from agent.search import fetch_pubmed_articles, fetch_arxiv_articles, save_raw_results
from agent.summarize_and_prepare_email import main as prepare_email_html
from agent import send_digest_email

def main():
    print("Fetching articles...")

    query = "echocardiography AND (deep learning OR artificial intelligence)"
    pubmed_articles = fetch_pubmed_articles(query, max_results=10)
    arxiv_articles = fetch_arxiv_articles(query, max_results=10)

    all_articles = pubmed_articles + arxiv_articles
    print(f"Total articles fetched: {len(all_articles)}")

    save_raw_results(pubmed_articles, arxiv_articles)

    print("Summarizing and sending email...")
    html_content = prepare_email_html()  # No arguments — it reads from CSV and generates HTML
    send_digest_email.send_email(html_content)

if __name__ == "__main__":
    main()
