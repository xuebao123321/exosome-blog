"use client";

import React, { useEffect, useState } from "react";

type Article = {
  id: string;
  title: string;
  pmid: string;
  source: string;
  date: string;
  time: string;
  abstract_zh: string;
  abstract_en: string;
};

export default function ArticlesPage() {
  const [articles, setArticles] = useState<Article[]>([]);

  useEffect(() => {
    fetch("/articles/index.json")
      .then((res) => res.json())
      .then(setArticles)
      .catch((err) => console.error("加载失败：", err));
  }, []);

  return (
    <main className="max-w-3xl mx-auto px-4 py-8">
      <h1 className="text-3xl font-bold mb-6">最新外泌体研究</h1>
      {articles.length === 0 ? (
        <p>正在加载文章...</p>
      ) : (
        <ul className="space-y-6">
          {articles.map((article) => (
            <li
              key={article.id}
              className="border p-4 rounded-2xl shadow hover:shadow-lg transition"
            >
              <h2 className="text-xl font-semibold">{article.title}</h2>
              <p className="text-sm text-gray-500 mb-1">
                来源：{article.source} | 发表时间：{article.date} | 抓取：{article.time}
              </p>
              <p className="mt-2 text-gray-700 line-clamp-4">{article.abstract_zh}</p>
              <a
                href={`https://pubmed.ncbi.nlm.nih.gov/${article.pmid}`}
                target="_blank"
                className="inline-block mt-3 text-blue-600 hover:underline"
              >
                查看原文 &rarr;
              </a>
            </li>
          ))}
        </ul>
      )}
    </main>
  );
}
