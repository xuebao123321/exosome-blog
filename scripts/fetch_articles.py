import os
import json
import openai
from datetime import datetime
from Bio import Entrez
from pathlib import Path

# 设置 Entrez 和 OpenAI
Entrez.email = "your_email@example.com"
openai.api_key = os.getenv("OPENAI_API_KEY")

KEYWORD = "exosome"
MAX_RESULTS = 5
OUTPUT_DIR = Path("public/articles")
MARKDOWN_DIR = Path("public/posts")
HTML_DIR = Path("public/html")

# 自动提取简短标题
def extract_short_title(text):
    prompt = f"从以下英文摘要生成一个不超过15个汉字的中文标题：\n\n{text}"
    try:
        response = openai.ChatCompletion.create(
            model="gpt-3",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.7,
            max_tokens=100
        )
        return response["choices"][0]["message"]["content"].strip()
    except Exception as e:
        print("标题提取失败：", e)
        return "未知标题"

# 翻译摘要
def translate_text(text):
    prompt = f"请将以下英文科学摘要翻译成通顺的中文：\n\n{text}"
    try:
        response = openai.ChatCompletion.create(
            model="gpt-3",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.7,
            max_tokens=1000
        )
        return response["choices"][0]["message"]["content"].strip()
    except Exception as e:
        print("翻译失败：", e)
        return "[翻译失败]"

# 抓取文章 ID
def fetch_pubmed_ids():
    with Entrez.esearch(db="pubmed", term=KEYWORD, retmax=MAX_RESULTS, sort="pub+date") as handle:
        record = Entrez.read(handle)
    return record.get("IdList", [])

# 获取文章详细信息
def fetch_article_details(pmid):
    with Entrez.efetch(db="pubmed", id=pmid, retmode="xml") as handle:
        records = Entrez.read(handle)
    articles = records.get("PubmedArticle", [])
    if articles:
        return articles[0]
    return {}

# 生成 Markdown 内容
def generate_markdown(title, abstract_en, abstract_zh, pmid):
    return f"""---
title: "{title}"
date: {datetime.now().strftime('%Y-%m-%d %H:%M')}
pmid: {pmid}
---

## 英文摘要
{abstract_en}

## 中文翻译
{abstract_zh}
"""

# 生成 HTML 内容（可选）
def generate_html(title, abstract_en, abstract_zh, pmid):
    return f"""<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>{title}</title></head>
<body>
  <h1>{title}</h1>
  <p><strong>PMID:</strong> {pmid}</p>
  <h2>英文摘要</h2>
  <p>{abstract_en}</p>
  <h2>中文翻译</h2>
  <p>{abstract_zh}</p>
</body>
</html>
"""

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    MARKDOWN_DIR.mkdir(parents=True, exist_ok=True)
    HTML_DIR.mkdir(parents=True, exist_ok=True)

    pmids = fetch_pubmed_ids()
    print("抓取到的 PubMed ID：", pmids)

    all_data = []

    for idx, pmid in enumerate(pmids):
        details = fetch_article_details(pmid)
        article_title = details.get("MedlineCitation", {}).get("Article", {}).get("ArticleTitle", "No Title")
        abstracts = details.get("MedlineCitation", {}).get("Article", {}).get("Abstract", {}).get("AbstractText", [])
        if isinstance(abstracts, list):
            abstract_en = " ".join(abstracts)
        else:
            abstract_en = abstracts or "无摘要内容"
        source = details.get("MedlineCitation", {}).get("MedlineJournalInfo", {}).get("MedlineTA", "")
        pub_date = details.get("MedlineCitation", {}).get("Article", {}).get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        pub_time = pub_date.get("Year", "Unknown")

        # 翻译 & 提取标题
        abstract_zh = translate_text(abstract_en)
        short_title = extract_short_title(abstract_en)

        data = {
            "id": f"exosome-{idx}",
            "keyword": KEYWORD,
            "title": short_title,
            "pmid": pmid,
            "source": source,
            "date": pub_time,
            "time": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "abstract_en": abstract_en,
            "abstract_zh": abstract_zh
        }

        all_data.append(data)

        # 写入 JSON 文件（单条）
        with open(OUTPUT_DIR / f"{data['id']}.json", "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

        # 写入 Markdown 文件
        with open(MARKDOWN_DIR / f"{data['id']}.md", "w", encoding="utf-8") as f:
            f.write(generate_markdown(short_title, abstract_en, abstract_zh, pmid))

        # 写入 HTML 页面（可选）
        with open(HTML_DIR / f"{data['id']}.html", "w", encoding="utf-8") as f:
            f.write(generate_html(short_title, abstract_en, abstract_zh, pmid))

    # 写入总索引 JSON 文件（用于前端展示）
    with open(OUTPUT_DIR / "index.json", "w", encoding="utf-8") as f:
        json.dump(all_data, f, indent=2, ensure_ascii=False)

    print(f"✅ 抓取完成，共写入 {len(pmids)} 篇文章。")

if __name__ == "__main__":
    main()
