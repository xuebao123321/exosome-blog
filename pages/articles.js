import fs from 'fs';
import path from 'path';

export async function getStaticProps() {
  const dir = path.join(process.cwd(), 'public/articles');
  const files = fs.readdirSync(dir);
  const articles = files.map(filename => {
    const content = fs.readFileSync(path.join(dir, filename), 'utf-8');
    return JSON.parse(content);
  });
  return {
    props: { articles }
  };
}

export default function Articles({ articles }) {
  return (
    <div style={{ maxWidth: 700, margin: "auto", padding: 20 }}>
      <h1>外泌体相关文章</h1>
      {articles.length === 0 && <p>暂无文章，可等待定时任务首次执行。</p>}
      {articles.map(a => (
        <div key={a.id} style={{ marginBottom: 40 }}>
          <h2>{a.title}</h2>
          <p><b>英文：</b> {a.abstract_en}</p>
          <p><b>中文翻译：</b> {a.abstract_zh}</p>
          <hr />
        </div>
      ))}
    </div>
    <div>
  {articles.map(article => (
    <div key={article.slug} style={{ marginBottom: "2rem", paddingBottom: "1rem", borderBottom: "1px solid #ccc" }}>
      <h2 style={{ fontSize: "1.5rem", fontWeight: "bold" }}>{article.title}</h2>
      <p style={{ fontSize: "0.9rem", color: "#555" }}>
        🗓️ 日期：{article.date} ⏰ 时间：{article.time}
      </p>
      <p style={{ fontSize: "0.9rem", color: "#555" }}>
        ✍️ 作者：{article.authors}
      </p>
      <p style={{ fontSize: "0.9rem" }}>
        🔗 来源：<a href={article.source} target="_blank" rel="noopener noreferrer" style={{ color: "#0070f3" }}>
          {article.source}
        </a>
      </p>
      <div style={{ marginTop: "1rem", lineHeight: "1.6" }}>
        {article.content}
      </div>
    </div>
  ))}
</div>

  );
}
