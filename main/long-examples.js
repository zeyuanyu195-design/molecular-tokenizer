/* Real saved IDs from the 899,836-record experiment; chemistry checked offline. */
(() => {
  const $=id=>document.getElementById(id), names={whole:'Whole-fragment · 16k',A:'A · Complete syntax blocks',B:'B · Protected structures'};
  const node=(tag,text,cls)=>{const e=document.createElement(tag);e.textContent=text;if(cls)e.className=cls;return e;};
  let data,example,budget='20000',fragment=0;
  let graphReady, clickVersion=0;
  function occurrences(model){
    let start=0, fi=0, within=0;
    return model.ids.map((id,index)=>{
      const entry=model.dictionary[id], end=start+entry.payload.length;
      const token={id,index,fi,within,start,end,entry};
      start=end;
      if(entry.kind==='dot'){fi++;within=0;}else{within++;}
      return token;
    });
  }
  function chips(container,tokens,mode){
    for(const token of tokens){
      const d=token.entry, b=node('button',d.text,'sc-token'+(d.kind==='dot'?' sc-dot':d.kind==='marker'?' sc-begin':''));
      b.type='button';b.append(node('small',`ID ${token.id}`));
      if(mode!=='whole'){
        b.dataset.scheme=mode;b.dataset.occurrence=token.index;b.setAttribute('aria-pressed','false');
      }
      b.onclick=async()=>{
        if(mode==='whole'){$('lx-inspector').textContent=`ID ${token.id} → ${d.text}. Exact dictionary payload.`;return;}
        if(fragment>=0 && (token.fi!==fragment || d.kind==='dot')){
          fragment=-1;$('lx-fragment').value='-1';render();
        }
        const ticket=++clickVersion;
        await graphReady;
        if(ticket!==clickVersion)return;
        const message=window.highlightExampleToken(mode,token);
        if(!message)return;
        for(const chip of document.querySelectorAll(`#long-examples [data-scheme="${mode}"]`)){
          chip.setAttribute('aria-pressed',String(Number(chip.dataset.occurrence)===token.index));
        }
        $('lx-inspector').textContent=message;
        const graph=$('graph-'+mode), rect=graph.getBoundingClientRect();
        if(rect.top<0 || rect.top>innerHeight*.55)graph.scrollIntoView({block:'center'});
      };
      container.append(b);
    }
  }
  function fragments(){
    $('lx-fragment').replaceChildren();
    const full = node('option', 'Whole molecule · all SAFE fragments');
    full.value = '-1';
    $('lx-fragment').append(full);
    example.fragments.forEach((text, index) => {
      const featured = example.featured.includes(index) ? ' · contrasting splits' : '';
      const option = node('option', `${index + 1}. ${text}${featured}`);
      option.value = index;
      $('lx-fragment').append(option);
    });
    fragment = -1;
    render();
  }
  function render(){
    const models=example.budgets[budget];
    clickVersion++;
    graphReady=window.showExampleGraphs(example,budget,fragment);
    $('lx-summary').textContent=`${example.chembl_id} · TEST record ${example.record_index.toLocaleString()} · ${example.fragments.length} original SAFE fragments. All dictionaries were trained in the 899,836-molecule experiment. Each total includes ${example.fragments.length-1} dot tokens.`;
    $('lx-smiles').textContent=example.smiles;$('lx-safe').textContent=example.safe;
    const entire = fragment === -1;
    $('lx-selected').textContent = entire ? example.safe : `Original fragment ${fragment+1}: ${example.fragments[fragment]}`;
    $('lx-cards').replaceChildren();$('lx-sequences').replaceChildren();
    for(const mode of ['whole','A','B']){
      const m=models[mode],p=m.fragments[fragment],card=node('article','',`sc-card sc-card-${mode}`),h=node('div','','sc-card-heading');
      h.append(node('h4',names[mode]),node('strong',`${m.length} tokens · whole molecule`));card.append(h);
      if (entire) {
        card.append(node('p', `${example.fragments.length} fragments · ${example.fragments.length - 1} dot tokens (ID 11)`, 'sc-subtle'));
      } else {
        card.append(node('p',`${p.ids.length} token${p.ids.length===1?'':'s'} · selected fragment`,'sc-subtle'));
        card.append(node('p',p.hit?'Frozen 16k whole-fragment hit: one unchanged ID.':mode==='whole'?'Whole-fragment miss → original base tokenizer fallback.':p.ids.length===1?'Complete entry found in the extended vocabulary → one ID.':'Whole-fragment miss → shortest legal span sequence, with framed lexical fallback when needed.','sc-subtle'));
        if(p.units.length){const d=node('details','');d.append(node('summary',`Allowed logical blocks (${p.units.length})`),node('code',p.units.join(' | '),'sc-block-list'));card.append(d);}
      }
      const tokens=occurrences(m);
      const box=node('div','','sc-token-box');chips(box,entire ? tokens : tokens.filter(t=>t.fi===fragment && t.entry.kind!=='dot'),mode);card.append(box);$('lx-cards').append(card);
      const d=node('details','');d.append(node('summary',`${names[mode]}: all ${m.length} model IDs`));const seq=node('div','','sc-token-box');chips(seq,tokens,mode);d.append(seq);$('lx-sequences').append(d);
    }
    $('lx-inspector').textContent='Click a token to inspect its actual ID. Whole-molecule totals and selected-fragment counts are shown separately. Markers and dots count as tokens.';
    $('lx-range').value=0;replay();
  }
  function replay(){const m=example.budgets[budget][$('lx-decode').value],s=$('lx-range');s.max=m.length;const n=Math.min(Number(s.value),m.length);s.value=n;const text=m.ids.slice(0,n).map(id=>m.dictionary[id].payload).join('');$('lx-restored').textContent=text||'(No SAFE characters emitted yet.)';$('lx-status').textContent=n===m.length?`${n}/${m.length} IDs: ${text===example.safe?'exact SAFE verified in this browser; molecular identity including stereochemistry verified offline.':'ERROR: mismatch'}`:`${n}/${m.length} IDs decoded. Partial strings need not be valid molecules.`;}
  async function init(){try{
    const response=await fetch('long-examples.json?v=20261118-multi');if(!response.ok)throw Error(response.status);data=await response.json();
    for(const e of data.examples)for(const ms of Object.values(e.budgets))for(const m of Object.values(ms)){
      if(m.ids.length!==m.length||m.ids.map(id=>m.dictionary[id].payload).join('')!==e.safe)throw Error('Dictionary audit failed');
    }
    data.examples.forEach((e,i)=>{const o=node('option',`${e.chembl_id} · ${e.fragments.length} fragments · ${['whole','A','B'].map(k=>e.budgets['20000'][k].length).join(' / ')} tokens (16k / A+4k / B+4k)`);o.value=i;$('lx-molecule').append(o);});
    $('lx-molecule').onchange=()=>{example=data.examples[Number($('lx-molecule').value)];$('lx-fragment').replaceChildren();fragments();};
    $('lx-budget').onchange=()=>{budget=$('lx-budget').value;render();};$('lx-fragment').onchange=()=>{fragment=Number($('lx-fragment').value);render();};
    $('lx-decode').onchange=()=>{$('lx-range').value=0;replay();};$('lx-range').oninput=replay;$('lx-next').onclick=()=>{$('lx-range').value=Number($('lx-range').value)+1;replay();};$('lx-complete').onclick=()=>{$('lx-range').value=$('lx-range').max;replay();};
    example=data.examples[0];fragments();$('lx-loading').hidden=true;$('lx-content').hidden=false;
  }catch(e){$('lx-loading').textContent='Unable to load verified examples: '+e.message;}}
  init();
})();
