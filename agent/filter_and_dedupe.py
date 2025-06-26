import os
import json
import pandas as pd
from datetime import datetime, timedelta

DATA_DIR = "data"
RAW_FILE = None
MASTER_FILE = os.path.join(DATA_DIR, "papers_master.csv")

def load_latest_raw():
    files = [f for f in os.listdir(DATA_DIR) if f.startswith("raw_") and f.endswith(".json")]
    if not files:
        print("No raw data files found.")
        return None
    latest_file = sorted(files)[-1]
    global RAW_FILE
    RAW_FILE = os.path.join(DATA_DIR, latest_file)
    print(f"Loading raw data from {RAW_FILE}")
    with open(RAW_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data

def raw_to_dataframe(raw_data):
    records = []
    for source_data in raw_data:
        source = source_data.get("source")
        papers = source_data.get("papers", [])
        for paper in papers:
            if source == "pubmed":
                record = {
                    "source": source,
                    "id": paper.get("Id"),
                    "title": paper.get("Title"),
                    "abstract": paper.get("Title") + " (Abstract unavailable)",
                    "published": paper.get("PubDate"),
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{paper.get('Id')}/"
                }
            elif source == "arxiv":
                record = {
                    "source": source,
                    "id": paper.get("link", paper.get("entry_id", "")),
                    "title": paper.get("title"),
                    "abstract": paper.get("summary"),
                    "published": paper.get("published"),
                    "link": paper.get("link")
                }
            else:
                continue
            records.append(record)
    df = pd.DataFrame(records)
    return df

def filter_by_date(df, days=7):
    cutoff = datetime.today() - timedelta(days=days)
    df["published"] = pd.to_datetime(df["published"], errors="coerce")
    filtered_df = df[df["published"] >= cutoff]
    return filtered_df

def load_master():
    if os.path.exists(MASTER_FILE):
        return pd.read_csv(MASTER_FILE)
    else:
        columns = ["source", "id", "title", "abstract", "published", "link"]
        empty_df = pd.DataFrame(columns=columns)
        empty_df.to_csv(MASTER_FILE, index=False)
        return empty_df

def deduplicate(df, master_df):
    merged = pd.merge(df, master_df, on=["source", "id"], how="left", indicator=True)
    new_papers = merged[merged["_merge"] == "left_only"].drop(columns=["_merge"])
    return new_papers

def save_new_papers(df):
    if df.empty:
        print("No new papers found after filtering and deduplication.")
        return
    df_to_save = df[["source", "id", "title", "abstract", "published", "link"]]
    df_to_save.to_csv(MASTER_FILE, mode="a", header=False, index=False)
    print(f"Appended {len(df)} new papers to {MASTER_FILE}")

if __name__ == "__main__":
    raw_data = load_latest_raw()
    if raw_data is None:
        exit()

    df = raw_to_dataframe(raw_data)
    df = filter_by_date(df, days=7)

    master_df = load_master()
    new_papers = deduplicate(df, master_df)

    save_new_papers(new_papers)

    print("Filtering and deduplication done. New papers ready for summary.")
