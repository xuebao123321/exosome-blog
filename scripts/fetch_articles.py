import os
import json
import openai
from datetime import datetime
from Bio import Entrez

Entrez.email = "your_email@example.com"  # 可替换为你自己的邮箱
openai.api_key = os.getenv("OPENAI_API_KEY")
keywords = ["exosome", "癌症", "外泌体治疗"]

def fetch_articles():
    all_results = []
    for keyword in keywords:
        handle = Entrez.esearch(db="pubmed", term=keyword, retmax=3, sort="relevance")
        record = Entrez.read(handle)
        ids = record["IdList"]
        if not ids:
            continue
        summaries = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
        abstracts = summaries.read().split("\n\n")
        for i, abstract in enumerate(abstracts):
            if not abstract.strip():
                continue
            title = f"{keyword} 相关文章 {i+1}"
            result = {
                "id": f"{keyword}-{i}",
                "keyword": keyword,
                "title": title,
                "abstract_en": abstract.strip(),
                "abstract_zh": translate(abstract.strip())
            }
            all_results.append(result)
    return all_results

def translate(text):
    prompt = f"请将以下英文摘要翻译成简洁专业的中文：\\n\\n{text}"
    try:
        response = openai.ChatCompletion.create(
            model="gpt-4",
            messages=[{"role": "user", "content": prompt}]
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        return f"[翻译失败] {str(e)}"

def main():
    articles = fetch_articles()
    os.makedirs("public/articles", exist_ok=True)
    for article in articles:
        with open(f"public/articles/{article['id']}.json", "w") as f:
            json.dump(article, f, ensure_ascii=False, indent=2)

if __name__ == "__main__":
    main()
