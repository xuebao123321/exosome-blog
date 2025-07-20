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
    <div style={{ maxWidth: 700, margin: 'auto', padding: 20 }}>
      <h1>外泌体相关文章</h1>
      {articles.map((a) => (
        <div key={a.id} style={{ marginBottom: 40 }}>
          <h2>{a.title}</h2>
          <p><b>英文摘要：</b><br />{a.abstract_en}</p>
          <p><b>中文翻译：</b><br />{a.abstract_zh}</p>
          <hr />
        </div>
      ))}
    </div>
  );
}
