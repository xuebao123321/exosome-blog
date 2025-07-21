import fs from "fs";
import path from "path";

export async function getStaticProps() {
  const dir = path.join(process.cwd(), "public/articles");
  const files = fs.existsSync(dir)
    ? fs.readdirSync(dir).filter(f => f.endsWith(".json"))
    : [];

  const articles = files.map(filename => {
    const filePath = path.join(dir, filename);
    const data = JSON.parse(fs.readFileSync(filePath, "utf-8"));

    // 自动标题截断（最多15个字）
    let shortTitle = data.title || "无标题";
    if (shortTitle.length > 15) {
      shortTitle = shortTitle.substring(0, 15) + "…";
    }

    const date = data.date || filename.split("_")[0]; // 从文件名提取日期
    const time = data.time || "未知时间";
    const source = data.link || "#";

    return {
      id: data.id || filename,
      title: shortTitle,
      fullTitle: data.title || "无标题",
      abstract_en: data.abstract_en || "",
      abstract_zh: data.abstract_zh || "",
      date,
      time,
      source
    };
  });

  return { props: { articles } };
}

export default function Articles({ articles }) {
  return (
    <div style={{ maxWidth: 800, margin: "auto", padding: "2rem" }}>
      <h1 style={{ fontSize: "2rem", marginBottom: "1.5rem" }}>🧬 外泌体文章更新</h1>
      {articles.length === 0 && <p>暂无文章，请等待采集。</p>}
      {articles.map(article => (
        <div key={article.id} style={{ marginBottom: "2rem", borderBottom: "1px solid #ccc", paddingBottom: "1rem" }}>
          <h2 style={{ fontSize: "1.2rem", color: "#333" }}>{article.title}</h2>
          <p style={{ fontSize: "0.9rem", color: "#666" }}>
            📅 {article.date} 🕘 {article.time}
          </p>
          <p style={{ fontSize: "0.9rem" }}>
            🔗 来源链接：<a href={article.source} target="_blank" rel="noopener noreferrer">{article.source}</a>
          </p>
          <details style={{ marginTop: "0.5rem" }}>
            <summary style={{ cursor: "pointer", fontWeight: "bold" }}>查看详情</summary>
            <p><strong>英文摘要：</strong> {article.abstract_en}</p>
            <p><strong>中文翻译：</strong> {article.abstract_zh}</p>
          </details>
        </div>
      ))}
    </div>
  );
}
