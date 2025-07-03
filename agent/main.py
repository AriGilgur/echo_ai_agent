import os
from agent.search import fetch_pubmed_articles, fetch_arxiv_articles, save_raw_results
from agent.summarize_and_prepare_email import main as prepare_email_html
from agent import send_digest_email

def main():
    print("Fetching articles...")

    query = "echocardiography AND (deep learning OR artificial intelligence)"
    pubmed_articles = fetch_pubmed_articles(query, max_results=10)
    arxiv_articles = fetch_arxiv_articles(query, max_results=10)

    save_raw_results(pubmed_articles, arxiv_articles)

    all_articles = pubmed_articles + arxiv_articles
    print(f"Total articles fetched: {len(all_articles)}")

    # Optionally call extract_emails() if implemented
    # for article in all_articles:
    #     extract_emails(article)

    # Simple email content with article titles and links
    titles_html = "<ul>" + "".join(
        f"<li><a href='{a['link']}'>{a['title']}</a></li>" for a in all_articles
    ) + "</ul>"
    html_content = f"<h1>Weekly Echo-AI Articles</h1>{titles_html}"

    print("Sending email...")
    send_digest_email.send_email(html_content)

if __name__ == "__main__":
    main()

