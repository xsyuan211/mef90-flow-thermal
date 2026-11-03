// 生成用于 GitHub API PUT 的 JS 代码片段
const fs = require('fs');
const path = require('path');

const filePath = path.join(__dirname, 'm_Flow', 'm_MEF90_FlowAssembly.F90');
const content = fs.readFileSync(filePath, 'utf8');
const b64 = Buffer.from(content).toString('base64');

// 输出可在浏览器 DevTools 执行的 JS
const js = `
const b64Content = "${b64.replace(/\\/g, '\\\\').replace(/"/g, '\\"')}";
const sha = "7fc6a4c909d9a70848cd5f76edeb462d8025b4d7";
fetch("https://api.github.com/repos/xsyuan211/mef90-flow-thermal/contents/m_Flow/m_MEF90_FlowAssembly.F90", {
  method: "PUT",
  headers: {
    "Accept": "application/vnd.github.v3+json",
    "Content-Type": "application/json"
  },
  body: JSON.stringify({
    message: "Fix SNES INTEGER type mismatch: extract Core subroutines to avoid PETSC_NULL_OBJECT passing",
    content: b64Content,
    sha: sha
  })
}).then(r => r.json()).then(j => {
  window.__gh_result = j;
  return j.commit?.sha || "NO_COMMIT_SHA";
}).catch(e => "ERR: " + e.message);
`;

fs.writeFileSync('_gh_put.js', js);
console.log('Generated _gh_put.js, b64 len:', b64.length);
