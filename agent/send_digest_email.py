import os
import pandas as pd
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

def load_digest():
    digest_path = "data/digest_ready.csv"
    if not os.path.exists(digest_path):
        raise FileNotFoundError("digest_ready.csv not found. Run summarization first.")
    return pd.read_csv(digest_path)

def format_html(df):
    html = "<h2>Weekly Digest: Echo-AI Articles</h2><ul>"
    for _, row in df.iterrows():
        html += f"<li><strong>{row['title']}</strong><br>"
        html += f"{row['summary']}<br>"
        html += f"<i>Lead Author: {row['lead_author']}, Email: {row['author_email']}</i><br>"
        html += f"<a href='{row['link']}'>Read Full Article</a></li><br><br>"
    html += "</ul>"
    return html

# This function is exported so main.py can import and use it
def send_digest_email(html_content):
    raw_recipients = os.environ.get("RECIPIENTS", "")
    recipient_list = [email.strip() for email in raw_recipients.split(",") if email.strip()]

    message = Mail(
        from_email='gilgurari@gmail.com',
        to_emails=recipient_list,
        subject="Weekly Digest: Echo-AI Articles",
        html_content=html_content
    )

    try:
        sg = SendGridAPIClient(os.environ.get("SENDGRID_API_KEY"))
        response = sg.send(message)
        print(f"Email sent! Status code: {response.status_code}")
    except Exception as e:
        print(f"Failed to send email: {e}")

def main():
    df = load_digest()
    html_content = format_html(df)
    send_digest_email(html_content)

if __name__ == "__main__":
    main()
