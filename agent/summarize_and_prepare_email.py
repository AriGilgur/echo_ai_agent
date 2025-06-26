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

def main():
    df = pd.read_csv(MASTER_FILE)

    summaries = []
    for idx, row in df.iterrows():
        abstract = row.get("abstract", "")
        summary = summarize_abstract(abstract)
        summaries.append(summary)

    df["summary"] = summaries

    # Keep only relevant columns for the digest
    digest_df = df[["title", "summary", "lead_author", "author_email", "link"]]

    # Save or prepare for email sending
    digest_df.to_csv("data/digest_ready.csv", index=False)
    print("Summaries created and saved to data/digest_ready.csv")

if __name__ == "__main__":
    main()
