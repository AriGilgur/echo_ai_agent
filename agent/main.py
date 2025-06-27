import os
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail
from Bio import Entrez
import arxiv

# Set your email
Entrez.email = "youremail@example.com"  # ← CHANGE THIS

def fetch_pubmed(term="echo AI", retmax=5):
    from Bio import Entrez
    Entrez.email = "youremail@example.com"
    handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
    record = Entrez.read(handle)
    ids = record['IdList']
    summaries = []
    if ids:
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        for article in records['PubmedArticle']:
            data = article['MedlineCitation']['Article']
            title = data.get('ArticleTitle', '')
            authors = ", ".join([f"{a['LastName']} {a['Initials']}" for a in data.get('AuthorList', []) if 'LastName' in a])
            summaries.append(f"**{title}**\nAuthors: {authors}\n")
    return summaries

def fetch_arxiv(query="echo AI", max_results=5):
    search = arxiv.Search(query=query, max_results=max_results, sort_by=arxiv.SortCriterion.SubmittedDate)
    summaries = []
    for result in search.results():
        authors = ", ".join([a.name for a in result.authors])
        summaries.append(f"**{result.title}**\nAuthors: {authors}\n{result.entry_id}\n")
    return summaries

def send_email(content):
message = Mail(
    from_email='gilgurari@gmail.com',
    to_emails=os.environ.get('RECIPIENTS'),
    subject='Test Email from GitHub Actions',
    plain_text_content='Hello! This is a test email from your GitHub workflow.'
)

    )
    try:
        sg = SendGridAPIClient(os.environ.get('SENDGRID_API_KEY'))
        response = sg.send(message)
        print(f"Email sent! Status code: {response.status_code}")
    except Exception as e:
        print(f"Error sending email: {e}")

def main():
    print("Fetching articles...")
    pubmed = fetch_pubmed()
    arxiv_papers = fetch_arxiv()
    body = "📚 Weekly Digest: echo AI articles\n\n"
    body += "🧬 PubMed:\n" + "\n".join(pubmed) + "\n\n"
    body += "🛰️ arXiv:\n" + "\n".join(arxiv_papers)
    send_email(body)

if __name__ == "__main__":
    main()
