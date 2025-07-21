import os
import requests
import openai
from datetime import datetime
from bs4 import BeautifulSoup

# 设置 OpenAI API Key
openai.api_key = os.getenv("OPENAI_API_KEY")

# PubMed 搜索关键词
KEYWORDS = ["exosome", "癌症", "外泌体治疗"]

# PubMed 搜索 API
def fetch_pubmed_articles(keywords):
    query = " OR ".join(keywords)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={query}&retmode=json&retmax=5"
    res = requests.get(url)
    ids = res.json().get("esearchresult", {}).get("idlist", [])
    return ids

# 获取文章详情
def fetch_article_details(pubmed_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pubmed_id}&retmode=json"
    res = requests.get(url)
    summary = res.json().get("result", {}).get(pubmed_id, {})
    title = summary.get("title", "No Title")
    authors = [a['name'] for a in summary.get("authors", [])]
    link = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
    return {"title": title, "authors": authors, "link": link}

# 翻译函数（使用 GPT）
def translate_text(text):
    try:
        response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            temperature=0.7,
            messages=[
                {"role": "system", "content": "你是一个优秀的生物医学翻译助手"},
                {"role": "user", "content": f"请将以下英文翻译成流畅自然的中文：\n\n{text}"}
            ],
            timeout=15,
        )
        return response["choices"][0]["message"]["content"].strip()
    except Exception as e:
        print("翻译失败：", str(e))
        return "❌ 翻译失败：OpenAI API 返回错误或没有响应"

# 保存为 Markdown 文件
def save_article(article, content_cn):
    slug = article["title"][:50].replace(" ", "-").replace("/", "-")
    filename = f"articles/{slug}.md"
    date = datetime.utcnow().strftime("%Y-%m-%d")

    md = f"""---
title: "{article['title']}"
date: {date}
authors: {', '.join(article['authors'])}
link: {article['link']}
---

{content_cn}
"""

    os.makedirs("articles", exist_ok=True)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(md)

# 主函数
def main():
    ids = fetch_pubmed_articles(KEYWORDS)
    print("获取文章 ID：", ids)

    for pmid in ids:
        details = fetch_article_details(pmid)
        print("获取文章：", details["title"])
        content_cn = translate_text(details["title"])
        save_article(details, content_cn)

if __name__ == "__main__":
    main()
