import os
import json
import re
from datetime import datetime, timedelta
from Bio import Entrez
import arxiv
import pytz

# Set your email here for NCBI API rules
Entrez.email = "anna@icardio.com"

def extract_email(text):
    if not text:
        return None
    matches = re.findall(r'[\w\.-]+@[\w\.-]+\.\w+', text)
    return matches[0] if matches else None

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
        print("No PubMed results found for query.")
        return []

    summaries = Entrez.esummary(db="pubmed", id=",".join(ids))
    results = Entrez.read(summaries)

    articles = []
    for res in results:
        lead_author = res.get("AuthorList", [])[0] if res.get("AuthorList") else "Unknown"
        email = extract_email(res.get("Summary", ""))
        articles.append({
            "title": res.get("Title"),
            "summary": res.get("Summary", ""),
            "authors": res.get("AuthorList", []),
            "lead_author": lead_author,
            "author_email": email or "Not available",
            "published": res.get("PubDate"),
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{res.get('Id')}/"
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





