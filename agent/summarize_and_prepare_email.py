import os
import pandas as pd
from datetime import datetime
from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()

openai_api_key = os.getenv("OPENAI_API_KEY")
if not openai_api_key:
    raise ValueError("OPENAI_API_KEY not set in environment variables")

client = OpenAI(api_key=openai_api_key)

MASTER_FILE = "data/papers_with_authors.csv"

def summarize_abstract(abstract_text):
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
    """
    Generate an organized HTML email body from a list of article dicts.
    Each dict should have: title, link, summary, lead_author, author_email.
    """
    html = """
    <html>
    <head>
      <style>
        body {
          font-family: Arial, sans-serif;
          background-color: #f8f9fa;
          color: #212529;
          padding: 20px;
        }
        .container {
          max-width: 700px;
          margin: auto;
          background-color: white;
          padding: 20px;
          border-radius: 8px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        h2 {
          text-align: center;
          color: #007BFF;
          margin-bottom: 30px;
        }
        .article {
          border-bottom: 1px solid #dee2e6;
          padding-bottom: 15px;
          margin-bottom: 15px;
        }
        .article:last-child {
          border-bottom: none;
          margin-bottom: 0;
          padding-bottom: 0;
        }
        .title {
          font-size: 18px;
          font-weight: bold;
          color: #0056b3;
          text-decoration: none;
        }
        .summary {
          margin: 8px 0;
          font-size: 14px;
          line-height: 1.4;
          color: #495057;
        }
        .author-info {
          font-size: 13px;
          color: #6c757d;
        }
        .read-more {
          display: inline-block;
          margin-top: 8px;
          padding: 6px 12px;
          background-color: #007BFF;
          color: white !important;
          text-decoration: none;
          border-radius: 4px;
          font-size: 13px;
        }
      </style>
    </head>
    <body>
      <div class="container">
        <h2>Weekly Echo-AI Articles Digest</h2>
    """

    for article in articles:
        title = article.get("title", "No Title")
        link = article.get("link", "#")
        summary = article.get("summary", "")
        lead_author = article.get("lead_author", "Unknown Author")
        author_email = article.get("author_email", "No email")

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
    df = pd.read_csv(MASTER_FILE)

    # Add summaries if not already present
    if "summary" not in df.columns or df["summary"].isnull().all():
        summaries = [summarize_abstract(row.get("abstract", "")) for _, row in df.iterrows()]
        df["summary"] = summaries

    digest_df = df[["title", "summary", "lead_author", "author_email", "link"]]
    digest_df.to_csv("data/digest_ready.csv", index=False)
    print("Summaries created and saved to data/digest_ready.csv")

    html_content = generate_digest_html(digest_df.to_dict(orient="records"))
    return html_content

if __name__ == "__main__":
    html = main()
    with open("data/email_digest.html", "w", encoding="utf-8") as f:
        f.write(html)
    print("HTML digest saved to data/email_digest.html")
