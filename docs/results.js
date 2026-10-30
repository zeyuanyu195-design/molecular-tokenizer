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
    const response=await fetch('results.json');if(!response.ok)throw new Error(`HTTP ${response.status}`);data=await response.json();
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
