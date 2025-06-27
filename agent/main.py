import os
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail
from agent.search import search_pubmed, search_arxiv, save_raw_results
from agent.extract_authors import extract_emails
from agent.summarize_and_prepare_email import generate_digest_html

def send_digest_email(html_content):
    message = Mail(
        from_email='gilgurari@gmail.com',
        to_emails=os.environ.get('RECIPIENTS'),
        subject='Weekly Digest: Echo-AI PubMed and arXiv Articles',
        html_content=html_content
    )
    try:
        sg = SendGridAPIClient(os.environ.get('SENDGRID_API_KEY'))
        response = sg.send(message)
        print(f"Email sent! Status code: {response.status_code}")
    except Exception as e:
        print(f"Error sending email: {e}")

def main():
    print("Fetching articles...")
    pubmed_results = search_pubmed()
    arxiv_results = search_arxiv()

    # Save raw results (optional but helpful for debugging)
    save_raw_results(pubmed_results, arxiv_results)

    # Extract papers list from results
    pubmed_articles = pubmed_results.get('papers', [])
    arxiv_articles = arxiv_results.get('papers', [])
    
    all_articles = pubmed_articles + arxiv_articles
    print(f"Total articles fetched: {len(all_articles)}")

    for article in all_articles:
        extract_emails(article)

    html_content = generate_digest_html(all_articles)
    send_digest_email(html_content)

if __name__ == "__main__":
    main()
