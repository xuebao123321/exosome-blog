import React from 'react';

export default function Home() {
  return (
    <div style={{ maxWidth: 600, margin: 'auto', padding: 20 }}>
      <h1>外泌体科普站</h1>
      <p>欢迎来到外泌体科普网站！这里将通过 AI 实时更新最新科研资讯。</p>
      <p>文章列表：</p>
      <ul>
        <li>示例文章1：什么是外泌体？</li>
        <li>示例文章2：外泌体在癌症研究中的应用</li>
      </ul>
    </div>
  );
}