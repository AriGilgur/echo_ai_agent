# agent/search.py

import os
import json
from datetime import datetime, timedelta
from Bio import Entrez
import arxiv

# Set your email here for NCBI API rules
Entrez.email = "anna@icardio.com"

def search_pubmed():
    today = datetime.today()
    last_week = today - timedelta(days=7)

    query = '("echocardiography"[MeSH]) AND ("deep learning" OR "artificial intelligence")'

    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        mindate=last_week.strftime("%Y/%m/%d"),
        maxdate=today.strftime("%Y/%m/%d"),
        datetype="pdat",
        retmax=10
    )
    record = Entrez.read(handle)
    ids = record["IdList"]

    if not ids:
        print("No PubMed results found for query.")
        return {"source": "pubmed", "papers": []}

    summaries = Entrez.esummary(db="pubmed", id=",".join(ids))
    results = Entrez.read(summaries)

    return {"source": "pubmed", "papers": results}

def search_arxiv():
    search = arxiv.Search(
        query='echocardiography AND (deep learning OR artificial intelligence)',
        max_results=10,
        sort_by=arxiv.SortCriterion.SubmittedDate,
        sort_order=arxiv.SortOrder.Descending
    )

    papers = []
    for result in search.results():
        papers.append({
            "title": result.title,
            "summary": result.summary,
            "authors": [a.name for a in result.authors],
            "published": result.published.strftime("%Y-%m-%d"),
            "link": result.entry_id,
        })

    return {"source": "arxiv", "papers": papers}

def save_raw_results(pubmed_results, arxiv_results):
    today_str = datetime.today().strftime("%Y%m%d")
    os.makedirs("data", exist_ok=True)
    path = f"data/raw_{today_str}.json"

    with open(path, "w", encoding="utf-8") as f:
        json.dump([pubmed_results, arxiv_results], f, indent=2)

    print(f"Saved results to {path}")

if __name__ == "__main__":
    pubmed_data = search_pubmed()
    arxiv_data = search_arxiv()
    save_raw_results(pubmed_data, arxiv_data)



