'use strict';
const $ = id => document.getElementById(id);
const esc = value => String(value).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt = n => Number(n).toLocaleString('en-US');
let data, cohort='shared', method='graph', step='encode', example=1, fragment=0;
function renderStats(){
  $('shared').setAttribute('aria-pressed',cohort==='shared');$('full').setAttribute('aria-pressed',cohort==='full');
  $('cohort-caption').textContent=cohort==='shared'?'44,328 identical molecules in every row, matched by record index, ChEMBL ID and original SMILES.':'Different inclusion criteria: each row uses its full available cohort. Fixed motifs retain only boundary-preserving records.';
  $('stats-table').querySelector('tbody').innerHTML=data.stats.filter(s=>s.cohort===cohort).map(s=>`<tr class="${s.unit.includes('excluded')?'motif':''}"><td><strong>${esc(s.method)}</strong><small>${esc(s.unit)}</small></td><td>${fmt(s.n)}</td><td>${fmt(s.total)}</td><td class="mean">${s.mean.toFixed(6)}</td><td>${fmt(s.median)}</td><td>${fmt(s.p95)}</td><td>${fmt(s.maximum)}</td></tr>`).join('');
}
function renderMethod(){
  const m=data.methods.find(m=>m.id===method);
  $('method-title').textContent=m.name;$('method-training').textContent=m.training;$('method-decode').textContent=m.decode;
  $('method-flow').innerHTML=m.flow.map((x,i)=>(i?'<span class="flow-arrow" aria-hidden="true">→</span>':'')+`<span class="flow-step">${esc(x)}</span>`).join('');
  $('pseudocode').textContent=m.code[step];document.querySelectorAll('[data-step]').forEach(b=>b.setAttribute('aria-pressed',b.dataset.step===step));
  const detailBox=$('graph-encode-details');
  detailBox.hidden=method!=='graph'||step!=='encode';
  if(!detailBox.hidden){
    detailBox.innerHTML='<h3 class="expanded-heading">Expanded encode steps</h3><p>Open each step for the implementation and source links. These are explanatory algorithms for the existing code; no vocabulary is retrained.</p>'+data.graph_encode_details.map(d=>`<details class="encode-detail"><summary>${esc(d.title)}</summary><div class="detail-body"><p>${esc(d.description)}</p><pre><code>${esc(d.code)}</code></pre><p>${esc(d.notes)}</p><div class="source-links">${d.sources.map(s=>`<a href="https://github.com/zeyuanyu195-design/molecular-tokenizer/blob/1ae2f0db56a0b1bd8d4140ac65da74b634321ae6/${s.path}#L${s.line}">${esc(s.label)} ↗</a>`).join('')}</div></div></details>`).join('')+
    '<aside class="note connection-note"><h3>Where does SAFE store attachment information?</h3><p><strong>SAFE serializes connectivity into the string.</strong> Dots separate fragments; matching labels identify the two endpoints. The atom carrying each label specifies the attachment position, and bond syntax specifies the bond type.</p><pre>CCO → cut C–O → CC[1*] . O[1*] → CC1.O1 → CCO</pre><p>This manually selected cut illustrates serialization, not a claim that GraphBPE always splits ethanol this way. Dummy-fragment spelling is schematic; the current converter round-trips <code>CC1.O1</code> to <code>CCO</code>.</p><p>In <code>CC1.O1</code>, the two “1” labels connect the second carbon to oxygen. The label is neither an atom index nor a token ID. Moving a label can change the molecule: <code>C1NCC.C1 → CCNCC</code>, whereas <code>CN1CC.C1 → CCN(C)C</code>. For a double-bond example, <code>CC=1.CC1 → CC=CC</code>; this example specifies bond order, not E/Z stereochemistry.</p><p><strong>No separate edge table is needed to decode complete SAFE.</strong> However, the current motif IDs omit this attachment syntax and cannot reconstruct the molecule alone. The reported motif counts exclude connectivity overhead.</p></aside>';
  }

}
function renderExample(){
  const e=data.examples[example];fragment=0;
  $('original-smiles').textContent=e.smiles;$('record-label').textContent=`${e.chembl_id} · Input record ${fmt(e.record_index)} · Boundaries preserved`;
  const card=(cls,title,n,unit,body,decoded)=>`<article class="method-card ${cls}"><h3>${title}</h3><div class="card-metric">${n}<span>${unit}</span></div>${body}<span class="label">DECODED SMILES</span><code>${esc(decoded)}</code><p class="check">Molecular round trip verified</p></article>`;
  $('example-comparison').innerHTML=card('brics','BRICS+SAFE',e.brics_ids.length,'SAFE text tokens',`<p>${e.brics_safe.split('.').length} SAFE fragment(s)</p><span class="label">COMPLETE SAFE</span><code>${esc(e.brics_safe)}</code>`,e.brics_decoded)+card('demo','Original DemoDiff',e.demo_ids.length,'node tokens · connectivity excluded',`<p>Original SMILES → GraphBPE</p><span class="label">NODE IDS</span><code>${esc(JSON.stringify(e.demo_ids))}</code><p>Separate edge/attachment information is required for decoding.</p>`,e.demo_decoded)+card('graph','GraphBPE+SAFE',e.graph_ids.length,'motif tokens · connectivity excluded',`<p>${e.graph_ids.length} SAFE fragments · ${e.graph_text_length} text tokens with the SAFE text tokenizer</p><span class="label">COMPLETE SAFE</span><code>${esc(e.graph_safe)}</code><span class="label">FIXED MOTIF IDS · SAFE ORDER</span><code>${esc(JSON.stringify(e.graph_ids))}</code>`,e.graph_decoded);
  $('fragment-buttons').innerHTML=e.graph_safe.split('.').map((f,i)=>`<button data-fragment="${i}" aria-pressed="${i===0}" aria-label="Fragment ${i+1}: ${esc(f)}">${i+1}. ${esc(f)}</button>`).join('');
  $('fragment-table').querySelector('tbody').innerHTML=e.graph_safe.split('.').map((f,i)=>`<tr data-row="${i}"><td>${i+1}</td><td>${esc(f)}</td><td>${esc(e.graph_keys[i])}</td><td>${e.graph_ids[i]}</td></tr>`).join('');
  $('connectivity').textContent=`node_ids = ${JSON.stringify(e.demo_ids)}\n\nNode index : motif key\n${e.demo_keys.map((k,i)=>`${i}: ${k}`).join('\n')}\n\nadj_edge_ids =\n${e.adj.map(r=>JSON.stringify(r)).join('\n')}\n\nbond_adj =\n${e.bond.map(r=>JSON.stringify(r)).join('\n')}\n\nposition_adj =\n${e.pos.map(r=>JSON.stringify(r)).join('\n')}\n\n-1 = no entry; matrix indices start at 0.`;
  $('text-tokens').textContent=`Tokens (${e.brics_tokens.length}):\n${JSON.stringify(e.brics_tokens)}\n\nIDs:\n${JSON.stringify(e.brics_ids)}`;
  renderFragment();
}
function renderFragment(){
  const e=data.examples[example];$('fragment-number').textContent=`SAFE FRAGMENT ${fragment+1}`;$('fragment-safe').textContent=e.graph_safe.split('.')[fragment];$('fragment-key').textContent=e.graph_keys[fragment];$('fragment-id').textContent=e.graph_ids[fragment];
  document.querySelectorAll('[data-fragment]').forEach(b=>b.setAttribute('aria-pressed',Number(b.dataset.fragment)===fragment));
  document.querySelectorAll('[data-row]').forEach(r=>r.classList.toggle('selected',Number(r.dataset.row)===fragment));
}
async function start(){
  try{
    const response=await fetch('results.json?v=20261103-pseudocode');if(!response.ok)throw new Error(`HTTP ${response.status}`);data=await response.json();
    $('shared').addEventListener('click',()=>{cohort='shared';renderStats();});$('full').addEventListener('click',()=>{cohort='full';renderStats();});
    $('method-select').addEventListener('change',e=>{method=e.target.value;renderMethod();});
    document.querySelectorAll('[data-step]').forEach(b=>b.addEventListener('click',()=>{step=b.dataset.step;renderMethod();}));
    $('example-select').addEventListener('change',e=>{example=Number(e.target.value);renderExample();});
    $('fragment-buttons').addEventListener('click',e=>{const b=e.target.closest('[data-fragment]');if(b){fragment=Number(b.dataset.fragment);renderFragment();}});
    $('fingerprints').innerHTML=data.sources.map(s=>`<p><strong>${esc(s.label)}</strong> · ${esc(s.file)}<code>SHA256 ${esc(s.sha256)}</code></p>`).join('');
    renderStats();renderMethod();renderExample();
  }catch(error){$('load-error').hidden=false;$('load-error').textContent='The result data could not be loaded. Please reload, or read the static report in the README.';console.error(error);}
}
start();
