/* Real saved IDs from the 899,836-record experiment; chemistry checked offline. */
(() => {
  const $=id=>document.getElementById(id), names={whole:'Whole-fragment · 16k',A:'A · Complete syntax blocks',B:'B · Protected structures'};
  const node=(tag,text,cls)=>{const e=document.createElement(tag);e.textContent=text;if(cls)e.className=cls;return e;};
  let data,example,budget='20000',fragment=0;
  function chips(container,ids,model){
    for(const id of ids){const d=model.dictionary[id], b=node('button',d.text,'sc-token'+(d.kind==='marker'?' sc-begin':''));b.type='button';b.append(node('small',`ID ${id}`));b.onclick=()=>{$('lx-inspector').textContent=`ID ${id} → ${d.text}. ${d.kind==='marker'?'Logical-block marker: counts as one model token; emits no SAFE text.':'Exact dictionary payload. IDs must be interpreted using this method’s dictionary.'}`;};container.append(b);}
  }
  function fragments(){
    const keep=Number($('lx-fragment').value), all=$('lx-all').checked;
    $('lx-fragment').replaceChildren();
    example.fragments.forEach((f,i)=>{if(all||example.featured.includes(i)){const o=node('option',`${i+1}. ${f}`);o.value=i;$('lx-fragment').append(o);}});
    fragment=[...$('lx-fragment').options].some(o=>Number(o.value)===keep)?keep:example.featured[0];
    $('lx-fragment').value=fragment;render();
  }
  function render(){
    const models=example.budgets[budget];
    $('lx-summary').textContent=`${example.chembl_id} · TEST record ${example.record_index.toLocaleString()} · ${example.fragments.length} original SAFE fragments. All dictionaries were trained in the 899,836-molecule experiment. Each total includes ${example.fragments.length-1} dot tokens.`;
    $('lx-smiles').textContent=example.smiles;$('lx-safe').textContent=example.safe;
    $('lx-selected').textContent=`Original fragment ${fragment+1}: ${example.fragments[fragment]}`;
    $('lx-cards').replaceChildren();$('lx-sequences').replaceChildren();
    for(const mode of ['whole','A','B']){
      const m=models[mode],p=m.fragments[fragment],card=node('article','',`sc-card sc-card-${mode}`),h=node('div','','sc-card-heading');
      h.append(node('h4',names[mode]),node('strong',`${m.length} tokens · whole molecule`));card.append(h,node('p',`${p.ids.length} token${p.ids.length===1?'':'s'} · selected fragment`,'sc-subtle'));
      card.append(node('p',p.hit?'Frozen 16k whole-fragment hit: one unchanged ID.':mode==='whole'?'Whole-fragment miss → original base tokenizer fallback.':p.ids.length===1?'Complete entry found in the extended vocabulary.':'Whole-fragment miss → shortest legal span sequence, with framed lexical fallback when needed.','sc-subtle'));
      if(p.units.length){const d=node('details','');d.append(node('summary',`Allowed logical blocks (${p.units.length})`),node('code',p.units.join(' | '),'sc-block-list'));card.append(d);}
      const box=node('div','','sc-token-box');chips(box,p.ids,m);card.append(box);$('lx-cards').append(card);
      const d=node('details','');d.append(node('summary',`${names[mode]}: all ${m.length} model IDs`));const seq=node('div','','sc-token-box');chips(seq,m.ids,m);d.append(seq);$('lx-sequences').append(d);
    }
    $('lx-inspector').textContent='Click a token to inspect its actual ID. Whole-molecule totals and selected-fragment counts are shown separately. Markers and dots count as tokens.';
    $('lx-range').value=0;replay();
  }
  function replay(){const m=example.budgets[budget][$('lx-decode').value],s=$('lx-range');s.max=m.length;const n=Math.min(Number(s.value),m.length);s.value=n;const text=m.ids.slice(0,n).map(id=>m.dictionary[id].payload).join('');$('lx-restored').textContent=text||'(No SAFE characters emitted yet.)';$('lx-status').textContent=n===m.length?`${n}/${m.length} IDs: ${text===example.safe?'exact SAFE verified in this browser; molecular identity including stereochemistry verified offline.':'ERROR: mismatch'}`:`${n}/${m.length} IDs decoded. Partial strings need not be valid molecules.`;}
  async function init(){try{
    const response=await fetch('long-examples.json?v=20261118');if(!response.ok)throw Error(response.status);data=await response.json();
    for(const e of data.examples)for(const ms of Object.values(e.budgets))for(const m of Object.values(ms)){
      if(m.ids.length!==m.length||m.ids.map(id=>m.dictionary[id].payload).join('')!==e.safe)throw Error('Dictionary audit failed');
    }
    data.examples.forEach((e,i)=>{const o=node('option',`${e.chembl_id} · ${e.fragments.length} fragments · ${['whole','A','B'].map(k=>e.budgets['20000'][k].length).join(' / ')} tokens (16k / A+4k / B+4k)`);o.value=i;$('lx-molecule').append(o);});
    $('lx-molecule').onchange=()=>{example=data.examples[Number($('lx-molecule').value)];$('lx-fragment').replaceChildren();fragments();};
    $('lx-budget').onchange=()=>{budget=$('lx-budget').value;render();};$('lx-all').onchange=fragments;$('lx-fragment').onchange=()=>{fragment=Number($('lx-fragment').value);render();};
    $('lx-decode').onchange=()=>{$('lx-range').value=0;replay();};$('lx-range').oninput=replay;$('lx-next').onclick=()=>{$('lx-range').value=Number($('lx-range').value)+1;replay();};$('lx-complete').onclick=()=>{$('lx-range').value=$('lx-range').max;replay();};
    example=data.examples[0];fragments();$('lx-loading').hidden=true;$('lx-content').hidden=false;
  }catch(e){$('lx-loading').textContent='Unable to load verified examples: '+e.message;}}
  init();
})();
