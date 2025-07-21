import os, json, openai
from datetime import datetime
from Bio import Entrez
from pathlib import Path

# ===== 配置 =====
Entrez.email = "your_email@example.com"  # ⚠️ 请替换为你自己的邮箱
openai.api_key = os.getenv("OPENAI_API_KEY")  # 确保 GitHub Secrets 已设置

KEYWORD = "exosome"
MAX_RESULTS = 5
OUT = Path("public/articles")

# ===== 获取 PubMed ID =====
def fetch_ids():
    try:
        h = Entrez.esearch(db="pubmed", term=KEYWORD, retmax=MAX_RESULTS, sort="pub+date")
        return Entrez.read(h).get("IdList", [])
    except Exception as e:
        print("抓取 PubMed ID 失败:", e)
        return []

# ===== 获取文章详情 =====
def fetch_details(pmid):
    try:
        with Entrez.efetch(db="pubmed", id=pmid, retmode="xml") as h:
            rec = Entrez.read(h)
        arts = rec.get("PubmedArticle", [])
        if not arts:
            print(f"未找到文章：{pmid}")
            return {}
        art = arts[0]
        cite = art.get("MedlineCitation", {})
        artinfo = cite.get("Article", {})
        title_en = artinfo.get("ArticleTitle", "")
        abst = artinfo.get("Abstract", {}).get("AbstractText", [""])
        abstract_en = " ".join(abst) if isinstance(abst, list) else abst
        src = cite.get("MedlineJournalInfo", {}).get("MedlineTA", "")
        pub = artinfo.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        year = pub.get("Year", "")
        month = pub.get("Month", "01")
        return {
            "pmid": pmid,
            "title_en": title_en,
            "abstract_en": abstract_en,
            "source": src,
            "pubyear": year,
            "pubmonth": month,
        }
    except Exception as e:
        print(f"获取文章详情失败（PMID:{pmid}）:", e)
        return {}

# ===== ChatGPT 翻译函数 =====
def translate(text, prompt):
    if not openai.api_key:
        print("❌ OPENAI_API_KEY 未设置")
        return "[未设置 API KEY]"
    try:
        print(f"🚀 调用 GPT 翻译（提示：{prompt[:10]}...）")
        r = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "user", "content": prompt + "\n\n" + text}],
            temperature=0.3,
            max_tokens=500
        )
        return r.choices[0].message.content.strip()
    except Exception as e:
        print("翻译失败:", e)
        with open("translate_error.log", "a", encoding="utf-8") as f:
            f.write(f"翻译失败: {e}\n")
        return "[翻译失败]"

# ===== 主程序入口 =====
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data_list = []

    pmids = fetch_ids()
    print(f"抓取到的 PubMed ID：{pmids}")

    for pmid in pmids:
        d = fetch_details(pmid)
        if not d or not d.get("abstract_en"):
            continue

        # 翻译和标题生成
        zh = translate(d["abstract_en"], "翻译为中文：")
        title = translate(d["abstract_en"], "生成不超过15字的中文标题：")

        dt = datetime.utcnow()
        rec = {
            "id": pmid,
            "title": title[:15] if title else "[无标题]",
            "full_title": d["title_en"],
            "abstract_en": d["abstract_en"],
            "abstract_zh": zh,
            "source": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "date": f"{d['pubyear']}-{d.get('pubmonth', '01')}",
            "time": dt.strftime("%Y-%m-%d %H:%M UTC")
        }

        # 写入单篇 JSON
        with open(OUT / f"{pmid}.json", "w", encoding="utf-8") as f:
            json.dump(rec, f, ensure_ascii=False, indent=2)

        data_list.append(rec)

    # 写入索引页
    with open(OUT / "index.json", "w", encoding="utf-8") as f:
        json.dump(data_list, f, ensure_ascii=False, indent=2)

if __name__ == "__main__":
    main()
