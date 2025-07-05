Echo-AI Email Agent

This GitHub repository contains a weekly email agent that automatically scrapes PubMed and arXiv for new echocardiography + AI research papers, summarizes the abstracts using GPT-4o, and sends a digest email every Monday.

✨ Features

Queries PubMed and arXiv for articles on "echocardiography" and "AI"

Extracts lead authors and emails (when available)

Summarizes abstracts with GPT-4o

Sends a weekly HTML email digest via SendGrid

Avoids resending the same articles using a master CSV

GitHub Actions automation (runs every Monday at 16:00 UTC)

📊 Architecture Diagram

graph TD
    A[GitHub Action Triggers Weekly] --> B[Run agent/main.py]
    B --> C[Fetch articles from PubMed and arXiv]
    C --> D[Summarize abstracts via OpenAI API]
    D --> E[Extract author names and emails]
    E --> F[Append to master CSV to prevent duplicates]
    F --> G[Generate HTML email digest]
    G --> H[Send via SendGrid]
    H --> I[Recipient: anna@icardio.ai]

📝 Setup Instructions

1. Clone the repo

git clone https://github.com/<your-username>/echo-ai-email-agent.git
cd echo-ai-email-agent

2. Install dependencies

pip install -r requirements.txt

3. Create .env file

OPENAI_API_KEY=your_openai_api_key
SENDGRID_API_KEY=your_sendgrid_api_key
RECIPIENTS=anna@icardio.ai

4. Set GitHub Secrets (for automated workflow)

OPENAI_API_KEY

SENDGRID_API_KEY

RECIPIENTS

5. Directory Structure

.
├── agent/
│   ├── main.py
│   ├── search.py
│   ├── extract_authors.py
│   └── summarize_and_prepare_email.py
├── data/
│   ├── papers_master.csv
│   └── raw_*.json
├── .github/workflows/agent.yml
└── requirements.txt

📷 Sample Email
![Screenshot 2025-07-05 153907](https://github.com/user-attachments/assets/a03ca75f-ea06-44e4-a640-732db7411237)



🚀 License

MIT License


