import re
import pandas as pd

MASTER_FILE = "data/papers_master.csv"
OUTPUT_FILE = "data/papers_with_authors.csv"

def extract_email(affiliation):
    if not affiliation:
        return None
    # Simple regex to find emails in affiliation string
    emails = re.findall(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", affiliation)
    return emails[0] if emails else None

def create_google_scholar_url(name, institution=None):
    base_url = "https://scholar.google.com/scholar?q="
    query = name
    if institution:
        query += f" {institution}"
    return base_url + query.replace(" ", "+")

def main():
    df = pd.read_csv(MASTER_FILE)

    # Add columns for lead_author and author_email (or fallback)
    lead_authors = []
    author_emails = []

    for _, row in df.iterrows():
        source = row["source"]
        # For PubMed, we’d ideally have affiliation data, but here we simulate
        # If you saved affiliations separately, load them here to extract emails
        lead_author = None
        email = None

        # Example: extract lead author from title/authors field
        # (In practice, you need to extract authors list from raw data)
        # For demo, just assign dummy data
        lead_author = "First Author"  # Replace with real extraction logic
        email = None  # Replace with real extraction logic

        # If no email found and source is arxiv, create Google Scholar fallback link
        if not email and source == "arxiv":
            email = create_google_scholar_url(lead_author)

        lead_authors.append(lead_author)
        author_emails.append(email)

    df["lead_author"] = lead_authors
    df["author_email"] = author_emails

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Saved papers with authors to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
