import os
import pandas as pd
from datetime import datetime
from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()  # loads .env file if you have one with OPENAI_API_KEY

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
    summary = response.choices[0].message.content.strip()
    return summary

def generate_digest_html(articles):
    """
    Generate an HTML email body from a list of article dicts.
    Each article dict must contain keys: 'title', 'link', 'summary', 'lead_author', 'author_email'.
    """
    html = """
    <html>
    <body>
        <h2>Weekly Echo-AI Articles Digest</h2>
        <ul>
    """
    for article in articles:
        title = article.get("title", "No Title")
        link = article.get("link", "#")
        summary = article.get("summary", "")
        lead_author = article.get("lead_author", "Unknown Author")
        author_email = article.get("author_email", "No email")

        html += f"""
        <li style="margin-bottom: 20px;">
            <a href="{link}" style="font-size: 18px; font-weight: bold; text-decoration: none; color: #2a6ebb;">{title}</a><br/>
            <small>Lead Author: {lead_author} ({author_email})</small><br/>
            <p style="max-width: 600px;">{summary}</p>
        </li>
        """

    html += """
        </ul>
        <p>--<br/>Echo-AI Weekly Digest</p>
    </body>
    </html>
    """
    return html

def main():
    # Load papers with authors info
    df = pd.read_csv(MASTER_FILE)

    # Generate summaries for each abstract
    summaries = []
    for idx, row in df.iterrows():
        abstract = row.get("abstract", "")
        summary = summarize_abstract(abstract)
        summaries.append(summary)
    df["summary"] = summaries

    # Save the digest CSV for records or other uses
    digest_df = df[["title", "summary", "lead_author", "author_email", "link"]]
    digest_df.to_csv("data/digest_ready.csv", index=False)
    print("Summaries created and saved to data/digest_ready.csv")

    # Convert to list of dicts for HTML generation
    articles = digest_df.to_dict(orient="records")
    html_content = generate_digest_html(articles)

    # Return the html_content so your main script can send it as email body
    return html_content


if __name__ == "__main__":
    # For testing, just generate and save the html digest to a file
    html = main()
    with open("data/email_digest.html", "w", encoding="utf-8") as f:
        f.write(html)
    print("HTML digest saved to data/email_digest.html")
