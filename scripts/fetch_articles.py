import os, json, openai
from datetime import datetime
from Bio import Entrez
from pathlib import Path

# 配置
Entrez.email = "your_email@example.com"
openai.api_key = os.getenv("OPENAI_API_KEY")
KEYWORD = "exosome"
MAX_RESULTS = 5
OUT = Path("public/articles")

def fetch_ids():
    h = Entrez.esearch(db="pubmed", term=KEYWORD, retmax=MAX_RESULTS, sort="pub+date")
    return Entrez.read(h).get("IdList", [])

def fetch_details(pmid):
    with Entrez.efetch(db="pubmed", id=pmid, retmode="xml") as h:
        rec = Entrez.read(h)
    arts = rec.get("PubmedArticle", [])
    if not arts: return {}
    art = arts[0]
    cite = art.get("MedlineCitation", {})
    artinfo = cite.get("Article", {})
    title_en = artinfo.get("ArticleTitle", "")
    abst = artinfo.get("Abstract", {}).get("AbstractText", [""])
    abstract_en = " ".join(abst) if isinstance(abst, list) else abst
    src = cite.get("MedlineJournalInfo", {}).get("MedlineTA", "")
    pub = artinfo.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
    year = pub.get("Year", "")
    month = pub.get("Month", "")
    return {
        "pmid": pmid,
        "title_en": title_en,
        "abstract_en": abstract_en,
        "source": src,
        "pubyear": year,
        "pubmonth": month or "",
    }

def translate(text, prompt):
    if not openai.api_key:
        print("未设置 OPENAI_API_KEY！")
        return "[未设置 API KEY]"
    try:
        print(f"翻译中：{prompt[:10]}...")
        r = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "user", "content": prompt + "\n\n" + text}],
            temperature=0.3,
            max_tokens=500
        )
        return r.choices[0].message.content.strip()
    except Exception as e:
        print("翻译失败:", e)
        # ✅ 这行必须缩进！确保是 try 的 except 内部！
        with open("translate_error.log", "a", encoding="utf-8") as f:
            f.write(f"翻译失败: {e}\n")
        return "[翻译失败]"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data_list = []
    for pmid in fetch_ids():
        d = fetch_details(pmid)
        if not d: continue
        zh = translate(d["abstract_en"], "翻译为中文：")
        title = translate(d["abstract_en"], "生成不超15字中文标题：")
        dt = datetime.utcnow()
        rec = {
            "id": pmid,
            "title": title[:15],
            "full_title": d["title_en"],
            "abstract_en": d["abstract_en"],
            "abstract_zh": zh,
            "source": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "date": f"{d['pubyear']}-{d.get('pubmonth', '01') or '01'}",
            "time": dt.strftime("%Y-%m-%d %H:%M UTC")
        }
        data_list.append(rec)
        with open(OUT / f"{pmid}.json","w",encoding="utf-8") as f:
            json.dump(rec, f, ensure_ascii=False, indent=2)
    with open(OUT / "index.json","w",encoding="utf-8") as f:
        json.dump(data_list, f, ensure_ascii=False, indent=2)

if __name__=="__main__":
    main()
