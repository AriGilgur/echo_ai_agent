import os
import json
import re
from datetime import datetime, timedelta
from Bio import Entrez
import arxiv
import pytz

# Set your email here for NCBI API rules
Entrez.email = "anna@icardio.com"


def extract_email_from_affiliation(affiliation):
    # Simple regex to find emails in text
    email_pattern = r'[\w\.-]+@[\w\.-]+\.\w+'
    matches = re.findall(email_pattern, affiliation or "")
    return matches[0] if matches else None

def fetch_pubmed_articles_with_emails(query="echocardiography AI", max_results=10):
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
        print("No PubMed results found for query.")
        return []

    # Fetch full details for each article to get affiliations
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(handle)

    articles = []
    for article in records['PubmedArticle']:
        article_title = article['MedlineCitation']['Article']['ArticleTitle']
        abstract = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [""])[0]
        pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
        pmid = article['MedlineCitation']['PMID']
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

        authors = []
        emails = []
        author_list = article['MedlineCitation']['Article'].get('AuthorList', [])
        for author in author_list:
            name = ""
            if 'LastName' in author and 'ForeName' in author:
                name = f"{author['ForeName']} {author['LastName']}"
            elif 'CollectiveName' in author:
                name = author['CollectiveName']
            authors.append(name)

            # Try extract emails from AffiliationInfo
            affs = author.get('AffiliationInfo', [])
            for aff in affs:
                email = extract_email_from_affiliation(aff.get('Affiliation', ''))
                if email:
                    emails.append(email)

        lead_author = authors[0] if authors else "Unknown"
        author_email = emails[0] if emails else "No email found"

        articles.append({
            "title": article_title,
            "summary": abstract,
            "authors": authors,
            "published": str(pub_date),
            "link": link,
            "lead_author": lead_author,
            "author_email": author_email,
        })

    return articles


def fetch_arxiv_articles(query, max_results=25):
    search = arxiv.Search(
        query=query,
        max_results=max_results,
        sort_by=arxiv.SortCriterion.SubmittedDate,
        sort_order=arxiv.SortOrder.Descending
    )

    papers = []
    cutoff_date = datetime.now(pytz.UTC) - timedelta(days=7)

    for result in search.results():
        published_aware = result.published
        if published_aware.tzinfo is None:
            published_aware = published_aware.replace(tzinfo=pytz.UTC)

        if published_aware < cutoff_date:
            continue

        title = result.title.lower()
        summary = result.summary.lower()

        if "echocardiography" in title + summary and ("ai" in title + summary or "artificial intelligence" in title + summary):
            lead_author = result.authors[0].name if result.authors else "Unknown"
            papers.append({
                "title": result.title,
                "summary": result.summary,
                "authors": [a.name for a in result.authors],
                "lead_author": lead_author,
                "author_email": "Not available",  # arXiv does not provide emails
                "published": published_aware.strftime("%Y-%m-%d"),
                "link": result.entry_id,
            })

    return papers

def save_raw_results(pubmed_results, arxiv_results):
    today_str = datetime.today().strftime("%Y%m%d")
    os.makedirs("data", exist_ok=True)
    path = f"data/raw_{today_str}.json"

    with open(path, "w", encoding="utf-8") as f:
        json.dump({"pubmed": pubmed_results, "arxiv": arxiv_results}, f, indent=2)

    print(f"Saved results to {path}")

if __name__ == "__main__":
    pubmed = fetch_pubmed_articles()
    arxiv = fetch_arxiv_articles()
    save_raw_results(pubmed, arxiv)





