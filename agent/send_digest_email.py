from Bio import Entrez
import arxiv
from datetime import datetime, timedelta
import pytz
import os
import json
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail
from dotenv import load_dotenv

Entrez.email = "anna@icardio.com"

def fetch_pubmed_articles(query="echocardiography AI", max_results=10):
    today = datetime.today()
    last_week = today - timedelta(days=7)

    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        mindate=last_week.strftime("%Y/%m/%d"),
        maxdate=today.strftime("%Y/%m/%d"),
        datetype="pdat",
        retmax=max_results
    )
    record = Entrez.read(handle)
    ids = record["IdList"]
    if not ids:
        return []

    summaries = Entrez.esummary(db="pubmed", id=",".join(ids))
    results = Entrez.read(summaries)

    articles = []
    for res in results:
        authors = res.get("AuthorList", [])
        lead_author = authors[0] if authors else "Unknown"
        articles.append({
            "title": res.get("Title"),
            "abstract": "",  # You may add fetching abstract here later if needed
            "lead_author": lead_author,
            "author_email": "Not available",
            "published": res.get("PubDate"),
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{res.get('Id')}/"
        })
    return articles

def fetch_arxiv_articles(query, max_results=25):
    search = arxiv.Search(
        query=query,
        max_results=max_results,
        sort_by=arxiv.SortCriterion.SubmittedDate,
        sort_order=arxiv.SortOrder.Descending,
    )
    cutoff_date = datetime.now(pytz.UTC) - timedelta(days=7)

    papers = []
    for result in search.results():
        pub_date = result.published
        if pub_date.tzinfo is None:
            pub_date = pub_date.replace(tzinfo=pytz.UTC)
        if pub_date < cutoff_date:
            continue

        lead_author = result.authors[0].name if result.authors else "Unknown"
        papers.append({
            "title": result.title,
            "abstract": result.summary,
            "lead_author": lead_author,
            "author_email": "Not available",
            "published": pub_date.strftime("%Y-%m-%d"),
            "link": result.entry_id,
        })
    return papers

def save_to_csv(pubmed, arxiv):
    import pandas as pd
    all_articles = pubmed + arxiv
    df = pd.DataFrame(all_articles)
    os.makedirs("data", exist_ok=True)
    df.to_csv("data/papers_with_authors.csv", index=False)
    print(f"Saved {len(all_articles)} articles to data/papers_with_authors.csv")

load_dotenv()

def send_email(html_content):
    """
    Send an HTML email using SendGrid.
    """
    message = Mail(
        from_email="gilgurari@gmail.com",
        to_emails=["anna@icardio.ai", "vlad@abcmilwaukee.com"],
        subject="🫀 Weekly Echo-AI Research Digest",
        html_content=html_content
    )

    try:
        sg = SendGridAPIClient(os.getenv("SENDGRID_API_KEY"))
        response = sg.send(message)
        print(f"✅ Email sent: {response.status_code}")
    except Exception as e:
        print(f"❌ Error sending email: {str(e)}")

if __name__ == "__main__":
    query = "echocardiography AI"
    pubmed = fetch_pubmed_articles(query=query, max_results=20)
    arxiv = fetch_arxiv_articles(query=query, max_results=20)
    save_to_csv(pubmed, arxiv)
