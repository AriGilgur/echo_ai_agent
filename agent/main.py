import os
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

def send_test_email():
    print("Preparing to send test email...")
    message = Mail(
        from_email='youremail@example.com',  # Replace with your verified sender email
        to_emails=os.environ.get('RECIPIENTS'),  # Recipients from GitHub Secrets
        subject='Test Email from GitHub Actions',
        plain_text_content='Hello! This is a test email from your GitHub workflow.'
    )
    try:
        sg = SendGridAPIClient(os.environ.get('SENDGRID_API_KEY'))
        response = sg.send(message)
        print(f"Email sent! Status code: {response.status_code}")
        print(f"Response body: {response.body}")
    except Exception as e:
        print(f"Error sending email: {e}")

def main():
    print("Starting agent...")
    send_test_email()
    print("Finished sending email.")

if __name__ == "__main__":
    main()
