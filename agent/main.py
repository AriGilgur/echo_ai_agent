import os
import pandas as pd
from datetime import datetime
from dotenv import load_dotenv
from openai import OpenAI
from agent.search import fetch_pubmed_articles, fetch_arxiv_articles, save_to_csv
from agent.send_digest_email import send_digest_email

load_dotenv()

# Load OpenAI API key
openai_api_key = os.getenv("OPENAI_API_KEY")
if not openai_api_key:
    raise ValueError("OPENAI_API_KEY not set in environment variables")
client = OpenAI(api_key=openai_api_key)

MASTER_FILE = "data/papers_with_authors.csv"

def summarize_abstract(abstract_text):
    if not abstract_text:
        return "No abstract available."
    prompt = (
        "In two sentences, summarize the following abstract for a cardiology AI research team:\n"
        f"{abstract_text}\n"
    )
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": prompt},
        ],
        max_tokens=150,
        temperature=0.5,
    )
    return response.choices[0].message.content.strip()

def generate_digest_html(articles):
    html = """
    <html>
    <head>
      <style>
        body {font-family: Arial, sans-serif; background-color: #f8f9fa; color: #212529; padding: 20px;}
        .container {max-width: 700px; margin: auto; background-color: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);}
        h2 {text-align: center; color: #007BFF; margin-bottom: 30px;}
        .article {border-bottom: 1px solid #dee2e6; padding-bottom: 15px; margin-bottom: 15px;}
        .article:last-child {border-bottom: none; margin-bottom: 0; padding-bottom: 0;}
        .title {font-size: 18px; font-weight: bold; color: #0056b3; text-decoration: none;}
        .summary {margin: 8px 0; font-size: 14px; line-height: 1.4; color: #495057;}
        .author-info {font-size: 13px; color: #6c757d;}
        .read-more {display: inline-block; margin-top: 8px; padding: 6px 12px; background-color: #007BFF; color: white !important; text-decoration: none; border-radius: 4px; font-size: 13px;}
      </style>
    </head>
    <body>
      <div class="container">
        <h2>Weekly Echo-AI Articles Digest</h2>
    """
    for article in articles:
        title = article.get("title", "No Title")
        link = article.get("link", "#")
        summary = article.get("summary", "No summary available.")
        lead_author = article.get("lead_author", "Unknown Author")
        author_email = article.get("author_email", "No email provided")

        html += f"""
        <div class="article">
          <a href="{link}" class="title">{title}</a>
          <p class="summary">{summary}</p>
          <p class="author-info">Lead Author: {lead_author} | Email: {author_email}</p>
          <a href="{link}" class="read-more">Read Full Article</a>
        </div>
        """
    html += """
      </div>
    </body>
    </html>
    """
    return html

def main():
    print("Fetching articles...")
    query = "echocardiography AI"
    pubmed_articles = fetch_pubmed_articles(query=query, max_results=20)
    arxiv_articles = fetch_arxiv_articles(query=query, max_results=20)

    print(f"Fetched {len(pubmed_articles)} PubMed and {len(arxiv_articles)} arXiv articles.")

    save_to_csv(pubmed_articles, arxiv_articles)

    # Load saved CSV and summarize abstracts
    df = pd.read_csv(MASTER_FILE)
    summaries = [summarize_abstract(row.get("abstract", "")) for _, row in df.iterrows()]
    df["summary"] = summaries

    # Save updated CSV with summaries
    df.to_csv(MASTER_FILE, index=False)
    print(f"Summaries added and saved to {MASTER_FILE}")

    # Generate HTML email content
    articles_data = df.to_dict(orient="records")
    html_content = generate_digest_html(articles_data)

    # Send email (adjust recipient and sender inside send_email)
    print("Sending email...")
    send_email(html_content)
    print("Email sent.")

if __name__ == "__main__":
    main()

