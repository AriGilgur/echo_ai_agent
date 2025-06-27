# agent/main.py

import os
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

def format_papers(papers, source):
    if not papers["papers"]:
        return f"No new {source} papers found this week.\n"

    output = [f"\n### {source.upper()} Results:\n"]
    for paper in papers["papers"]:
        title = paper.get("title", "No title")
        authors = ", ".join(paper.get("authors", [])) if "authors" in paper else paper.get("AuthorList", ["N/A"])
        date = paper.get("published", "No date")
        link = paper.get("link", "No link")
        output.append(f"- **{title}**\n  - Authors: {authors}\n  - Date: {date}\n  - Link: {link}\n")
    return "\n".join(output)

def send_email(pubmed, arxiv):
    pubmed_text = format_papers(pubmed, "PubMed")
    arxiv_text = format_papers(arxiv, "arXiv")
    body = f"""Hello,\n\nHere is this week's AI echo article digest:\n{pubmed_text}\n{arxiv_text}\n\n– Echo Agent"""

    message = Mail(
        from_email='gilgurari@gmail.com',
        to_emails=os.environ.get('RECIPIENTS'),
        subject='Weekly AI Echocardiography Digest',
        plain_text_content=body
    )

    try:
        sg = SendGridAPIClient(os.environ.get('SENDGRID_API_KEY'))
        response = sg.send(message)
        print(f"✅ Email sent! Status code: {response.status_code}")
    except Exception as e:
        print(f"❌ Error sending email: {e}")

def main():
    print("🔍 Searching PubMed and arXiv...")
    pubmed = search_pubmed()
    arxiv = search_arxiv()
    save_raw_results(pubmed, arxiv)
    print("📧 Sending email...")
    send_email(pubmed, arxiv)

if __name__ == "__main__":
    main()




