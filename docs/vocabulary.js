"use strict";
(async function () {
  const byId = id => document.getElementById(id);
  try {
    const response = await fetch("vocabulary-results.json?v=20261109-examples");
    if (!response.ok) throw new Error(`Results HTTP ${response.status}`);
    const data = await response.json();
    let budget = "16000", method = "hybrid", step = "encode", mode = "text", exampleIndex = 0;
    const examples = data.examples;
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
      const selected = examples[exampleIndex];
      const count = selected.fragments.length, dotsCount = count - 1;
      byId("v-example-label").textContent = `${selected.chembl_id} · TEST RECORD ${selected.record_index} · FIXED 16k EXAMPLE`;
      byId("v-example-title").textContent = `${selected.old_length} → ${selected.length} tokens. ${count} fragments.`;
      byId("v-safe").textContent = selected.safe;
      byId("v-sequence-title").textContent = `All ${selected.length} hybrid tokens in order`;
      byId("v-count-explanation").textContent = `Fragment content: ${selected.fragments.map(f=>f.ids.length).join(" + ")} = ${selected.length-dotsCount} tokens; plus ${dotsCount} dots = ${selected.length}. Both complete encodings restore the same SAFE and canonical isomeric molecule. Complete-fragment matches are locked to one token; only unmatched fragments use fallback.`;
      const body=byId("v-example-body"), sequence=byId("v-sequence");
      sequence.setAttribute("aria-label",`Complete ${selected.length}-token sequence`);
      body.replaceChildren(); sequence.replaceChildren();
      selected.fragments.forEach((f,i)=> {
        const tr=document.createElement("tr");
        cell(tr,`${i+1}. ${f.text}`).className="v-fragment";
        tokens(cell(tr,""),f.old_tokens,f.old_ids);
        const hybridCell = cell(tr,"");
        const badge = document.createElement("small");
        badge.className = "v-match-label";
        badge.textContent = f.encode_mode === "whole" ? "Complete match · 1 token" : f.encode_mode === "base" ? "Base vocabulary · 1 token" : "Substring fallback";
        hybridCell.append(badge);
        tokens(hybridCell,f.tokens,f.ids);
        cell(tr,`${f.old_ids.length} → ${f.ids.length}`);
        body.append(tr);
        tokens(sequence,f.tokens,f.ids,i);
        if(i<dotsCount){const dot=document.createElement("span");dot.className="v-token v-dot";dot.textContent=mode==="ids"?"11":".";dot.title="Fragment separator · ID 11";sequence.append(dot);}
      });
      const dots=document.createElement("tr");["Separating dots",dotsCount,dotsCount,`${dotsCount} → ${dotsCount}`].forEach(v=>cell(dots,v));body.append(dots);
      const total=document.createElement("tr");["Total",selected.old_length,selected.length,`${selected.old_length} → ${selected.length}`].forEach(v=>cell(total,v));body.append(total);
      document.querySelectorAll("[data-v-mode]").forEach(button=>button.setAttribute("aria-pressed",String(button.dataset.vMode===mode)));
    }
    byId("v-budget").addEventListener("change",e=>{budget=e.target.value;results();});
    byId("v-method").addEventListener("change",e=>{method=e.target.value;pseudocode();});
    document.querySelectorAll("[data-v-step]").forEach(b=>b.addEventListener("click",()=>{step=b.dataset.vStep;pseudocode();}));
    document.querySelectorAll("[data-v-mode]").forEach(b=>b.addEventListener("click",()=>{mode=b.dataset.vMode;example();}));
    examples.forEach((e,i)=>{const option=document.createElement("option");option.value=String(i);option.textContent=`${e.chembl_id} · ${e.fragments.length} fragments · ${e.old_length} → ${e.length} tokens`;byId("v-example-select").append(option);});
    byId("v-example-select").addEventListener("change",e=>{exampleIndex=Number(e.target.value);example();});
    results(); pseudocode(); example();
    byId("v-loading").hidden=true; byId("v-content").hidden=false;
  } catch (error) {
    byId("v-loading").textContent=`Could not load the experiment: ${error.message}. The README contains the full report.`;
  }
})();
