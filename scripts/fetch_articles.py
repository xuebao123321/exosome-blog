import os
import json
from datetime import datetime
from Bio import Entrez
import openai

# 配置参数
Entrez.email = "your_email@example.com"  # 替换为你的邮箱
openai.api_key = os.getenv("OPENAI_API_KEY")

# 关键词设置
KEYWORDS = ["exosome", "癌症", "外泌体治疗"]

# 根据关键词从 PubMed 搜索文章 ID（这里仅搜索最近5篇）
def fetch_pubmed_ids():
    query = " OR ".join(KEYWORDS)
    handle = Entrez.esearch(db="pubmed", term=query, retmode="xml", retmax="5")
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])

# 根据文章 ID 获取文章详情
def fetch_article_details(pmid):
    handle = Entrez.esummary(db="pubmed", id=pmid, retmode="xml")
    summary = Entrez.read(handle)
    handle.close()
    result = summary.get("result", {}).get(pmid, {})
    # 如果没有标题或详情，则赋个默认值
    title = result.get("title", "无标题")
    # 构造 PubMed 链接
    link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    # 这里也可提取作者等信息，如 result.get("authors")（可选）
    return {
        "id": f"exosome-{pmid}",
        "title": title,
        "link": link,
        "abstract_en": result.get("title", "")  # 此处简化为使用标题作为摘要示例；你可以改为其他摘要字段
    }

# 使用 OpenAI API 翻译文本（从英文翻译成中文）
def translate_text(text):
    try:
        response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "system", "content": "你是一个优秀的生物医学翻译助手，尽量用简洁和专业的中文翻译下文。"},
                {"role": "user", "content": f"请将下面的英文翻译成中文：\n\n{text}"}
            ],
            temperature=0.3,
            timeout=15,
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print("翻译失败：", str(e))
        return "[翻译失败]"

# 自动截断标题至不超过15个字符（超长则加省略号）
def shorten_title(title):
    if len(title) > 15:
        return title[:15] + "…"
    return title

# 主函数，整合流程
def main():
    pmids = fetch_pubmed_ids()
    print("抓取到的 PubMed ID：", pmids)
    # 确保保存文件的目录存在
    output_dir = os.path.join("public", "articles")
    os.makedirs(output_dir, exist_ok=True)
    
    for pmid in pmids:
        details = fetch_article_details(pmid)
        # 翻译文章标题作为示例（你也可以翻译其他摘要内容）
        translated = translate_text(details["title"])
        
        # 获取当前日期和时间
        now = datetime.utcnow()
        date_str = now.strftime("%Y-%m-%d")
        time_str = now.strftime("%H:%M UTC")
        
        # 组装 JSON 数据，注意：这里用原标题生成，但前端会自动截断显示
        article = {
            "id": details["id"],
            "title": details["title"],
            "abstract_en": details["abstract_en"],
            "abstract_zh": translated,
            "link": details["link"],
            "date": date_str,
            "time": time_str
        }
        
        # 保存 JSON 文件，文件名可以用 date + id 形式保证唯一性
        filename = f"{date_str}_{pmid}.json"
        filepath = os.path.join(output_dir, filename)
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(article, f, ensure_ascii=False, indent=2)
        print(f"保存文章：{filename}")

if __name__ == "__main__":
    main()
