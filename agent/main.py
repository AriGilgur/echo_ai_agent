import os
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail
from agent.search import fetch_pubmed_articles, fetch_arxiv_articles, save_raw_results
# Assuming you have these modules, otherwise comment these out or implement accordingly:
# from agent.extract_authors import extract_emails
# from agent.summarize_and_prepare_email import generate_digest_html
from agent.summarize_and_prepare_email import main as prepare_email_html
from agent.send_email import send_digest_email

def main():
    print("Preparing email content...")
    html_content = prepare_email_html()
    print("Sending email...")
    send_digest_email(html_content)

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
    pubmed_articles = fetch_pubmed_articles()
    arxiv_articles = fetch_arxiv_articles()

    save_raw_results(pubmed_articles, arxiv_articles)

    all_articles = pubmed_articles + arxiv_articles
    print(f"Total articles fetched: {len(all_articles)}")

    # If you have extract_emails and generate_digest_html functions, uncomment these lines:
    # for article in all_articles:
    #     extract_emails(article)
    #
    # html_content = generate_digest_html(all_articles)

    # For now, let's send a simple HTML email listing titles
    titles_html = "<ul>" + "".join(f"<li>{a['title']}</li>" for a in all_articles) + "</ul>"
    html_content = f"<h1>Weekly Echo-AI Articles</h1>{titles_html}"

    send_digest_email(html_content)

if __name__ == "__main__":
    main()
