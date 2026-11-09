"use strict";
(async function () {
  const byId = id => document.getElementById(id);
  try {
    const response = await fetch("vocabulary-results.json?v=20261109");
    if (!response.ok) throw new Error(`Results HTTP ${response.status}`);
    const data = await response.json();
    let budget = "16000", method = "hybrid", step = "encode", mode = "text";
    const number = n => n.toLocaleString("en-US");
    function cell(row, text, tag = "td") {
      const node = document.createElement(tag); node.textContent = text; row.append(node); return node;
    }
    function results() {
      const r = data.summary.budgets[budget];
      byId("v-mean").textContent = `${r.whole_fragment.mean.toFixed(3)} → ${r.hybrid.mean.toFixed(3)}`;
      byId("v-reduction").textContent = `${(r.length_reduction * 100).toFixed(2)}% shorter`;
      const body = byId("v-stat-body"); body.replaceChildren();
      for (const [title, values] of [["Whole fragments + base fallback", r.whole_fragment], ["Whole fragments + BPE fallback", r.hybrid]]) {
        const tr = document.createElement("tr");
        [title, "44,329", number(values.total), values.mean.toFixed(6), number(values.median), number(values.p95)].forEach(v => cell(tr,v));
        body.append(tr);
      }
      byId("v-changes").textContent = `${number(r.shorter)} shorter · ${number(r.equal)} unchanged · ${number(r.longer)} longer. Both methods: 44,329 / 44,329 exact SAFE round trips.`;
      byId("v-allocation").textContent = budget === "16000" ? "16k hybrid = 2,136 base + 6,932 complete fragments + 6,932 substrings." : "32k hybrid = 2,136 base + 7,466 complete fragments + 22,398 substrings.";
    }
    function pseudocode() {
      byId("v-code-title").textContent = data.methods[method].title;
      byId("v-code").textContent = data.methods[method][step];
      document.querySelectorAll("[data-v-step]").forEach(button => button.setAttribute("aria-pressed",String(button.dataset.vStep===step)));
    }
    function tokens(node, values, ids, fragmentIndex) {
      values.forEach((value, i) => {
        const span = document.createElement("span"); span.className="v-token";
        span.textContent=mode==="ids" ? ids[i] : value;
        span.title=`${value} · ID ${ids[i]}`;
        if (fragmentIndex !== undefined) span.dataset.vFragment=fragmentIndex;
        node.append(span);
      });
    }
    function example() {
      const body=byId("v-example-body"), sequence=byId("v-sequence");
      body.replaceChildren(); sequence.replaceChildren();
      data.example.fragments.forEach((f,i)=> {
        const tr=document.createElement("tr");
        cell(tr,`${i+1}. ${f.text}`).className="v-fragment";
        tokens(cell(tr,""),f.old_tokens,f.old_ids);
        tokens(cell(tr,""),f.tokens,f.ids);
        cell(tr,`${f.old_ids.length} → ${f.ids.length}`);
        body.append(tr);
        tokens(sequence,f.tokens,f.ids,i);
        if(i<7){const dot=document.createElement("span");dot.className="v-token v-dot";dot.textContent=mode==="ids"?"11":".";dot.title="Fragment separator · ID 11";sequence.append(dot);}
      });
      const dots=document.createElement("tr");["Separating dots","7","7","7 → 7"].forEach(v=>cell(dots,v));body.append(dots);
      document.querySelectorAll("[data-v-mode]").forEach(button=>button.setAttribute("aria-pressed",String(button.dataset.vMode===mode)));
    }
    byId("v-budget").addEventListener("change",e=>{budget=e.target.value;results();});
    byId("v-method").addEventListener("change",e=>{method=e.target.value;pseudocode();});
    document.querySelectorAll("[data-v-step]").forEach(b=>b.addEventListener("click",()=>{step=b.dataset.vStep;pseudocode();}));
    document.querySelectorAll("[data-v-mode]").forEach(b=>b.addEventListener("click",()=>{mode=b.dataset.vMode;example();}));
    byId("v-safe").textContent=data.example.safe;
    results(); pseudocode(); example();
    byId("v-loading").hidden=true; byId("v-content").hidden=false;
  } catch (error) {
    byId("v-loading").textContent=`Could not load the experiment: ${error.message}. The README contains the full report.`;
  }
})();
