import os
import json
from datetime import datetime, timedelta
from Bio import Entrez
import arxiv

# Set your email here for NCBI API rules
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
        print("No PubMed results found for query.")
        return []

    summaries = Entrez.esummary(db="pubmed", id=",".join(ids))
    results = Entrez.read(summaries)

    articles = []
    for res in results:
        articles.append({
            "title": res.get("Title"),
            "summary": res.get("Summary", ""),
            "authors": res.get("AuthorList", []),
            "published": res.get("PubDate"),
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{res.get('Id')}/"
        })

    return articles

def fetch_arxiv_articles(query="echocardiography AI", max_results=10):
    search = arxiv.Search(
        query=query,
        max_results=max_results,
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





