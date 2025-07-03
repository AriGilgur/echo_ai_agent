import os
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail
from agent.search import fetch_pubmed_articles, fetch_arxiv_articles, save_raw_results
# Uncomment and implement if you have these
# from agent.extract_authors import extract_emails
# from agent.summarize_and_prepare_email import generate_digest_html
from agent.summarize_and_prepare_email import main as prepare_email_html
from agent import send_digest_email

def main():
    print("Fetching articles...")
    pubmed_articles = fetch_pubmed_articles()
    arxiv_articles = fetch_arxiv_articles()

    save_raw_results(pubmed_articles, arxiv_articles)

    all_articles = pubmed_articles + arxiv_articles
    print(f"Total articles fetched: {len(all_articles)}")

    # Uncomment if implemented:
    # for article in all_articles:
    #     extract_emails(article)
    #
    # html_content = generate_digest_html(all_articles)

    # For now, let's send a simple HTML email listing hyperlinked titles
    titles_html = "<ul>" + "".join(f"<li><a href='{a['link']}'>{a['title']}</a></li>" for a in all_articles) + "</ul>"
    html_content = f"<h1>Weekly Echo-AI Articles</h1>{titles_html}"

    print("Sending email...")
    send_digest_email.send_email(html_content)

if __name__ == "__main__":
    main()
