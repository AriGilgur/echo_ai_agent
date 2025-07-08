from agent.search import fetch_pubmed_articles, fetch_arxiv_articles, save_raw_results
from agent.summarize_and_prepare_email import main as prepare_email_html
from agent.send_digest_email import send_email

def main():
    query = "echocardiography AND (deep learning OR artificial intelligence)"
    print("Fetching articles...")
    pubmed_articles = fetch_pubmed_articles(query, max_results=10)
    arxiv_articles = fetch_arxiv_articles(query, max_results=10)

    save_raw_results(pubmed_articles, arxiv_articles)

    all_articles = pubmed_articles + arxiv_articles

    # Extract lead author and email for each article (assuming authors is a list of names)
    for art in all_articles:
        authors = art.get("authors", [])
        art["lead_author"] = authors[0] if authors else "Unknown"
        # For PubMed, author emails often not directly available; you may want to extend extraction if possible
        art["author_email"] = art.get("author_email", "No email")

    print(f"Total articles fetched: {len(all_articles)}")

    print("Summarizing and sending email...")
    html_content = prepare_email_html(all_articles)  # Pass articles with lead_author and author_email
    send_email(html_content)

if __name__ == "__main__":
    main()
