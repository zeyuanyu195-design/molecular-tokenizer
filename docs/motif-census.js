/* Verified descriptive census; never computes molecule reconstruction in-browser. */
(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const fmt = n => Number(n).toLocaleString('en-US');
  const pct = n => (100*n).toFixed(3)+'%';
  const names = {motif_graph:'Canonical motif graph',literal_SAFE:'Exact labelled SAFE string',normalized_label_pattern:'Label-normalized text pattern',canonical_port_graph:'Local connection graph'};
  const text = (tag,value) => {const e=document.createElement(tag);e.textContent=value;return e;};
  const row = values => {const r=document.createElement('tr');values.forEach(v=>r.append(text('td',v)));return r;};
  async function start(){
    const response=await fetch('motif-census.json?v=20261112');
    if(!response.ok)throw Error('Results unavailable ('+response.status+')');
    const d=await response.json();
    if(d.cohort_size!==899836||d.verified.status!=='verified'||d.motifs.length!==2740)throw Error('Unexpected census metadata');
    const s=d.summary,v=d.verified;
    const countBody=$('mc-count-body');
    Object.keys(names).forEach(key=>countBody.append(row([names[key],fmt(s.representations[key].unique_types),fmt(v.labelled_only_representations[key].unique_types)])));
    const variant=v.labelled_literal_variants_per_motif;
    $('mc-variants').textContent=`Across ${fmt(variant.observed_motif_keys)} canonical analysis groups with cross-fragment labels: mean ${variant.mean.toFixed(2)} literal variants; median ${variant.median}; maximum ${fmt(variant.max)}. ${fmt(variant.motifs_with_multiple_variants)} groups have multiple variants. These are group-level statistics, not per-original-ID statistics.`;
    function renderSharing(){
      const multi=$('mc-size').value==='multi';
      const prob=multi?'share_multiatom_labelled_probability':'share_any_labelled_probability';
      const ci=multi?'share_multiatom_labelled_95pct_ci':'share_any_labelled_95pct_ci';
      $('mc-share-body').replaceChildren();$('mc-bars').replaceChildren();
      for(const [key,name] of Object.entries(names)){
        const r=v.labelled_only_representations[key];
        $('mc-share-body').append(row([name,pct(r[prob]),r[ci].map(pct).join(' – ')]));
        const barRow=text('div','');barRow.className='mc-bar-row';barRow.append(text('span',name));
        const track=text('div','');track.className='mc-bar-track';
        const bar=text('div','');bar.className='mc-bar'+(key==='literal_SAFE'?' literal':'');bar.style.width=(r[prob]/.25*100)+'%';track.append(bar);
        barRow.append(track,text('strong',pct(r[prob])));$('mc-bars').append(barRow);
      }
      $('mc-size-note').textContent=multi?'Only fragments containing at least two real atoms and at least one external connection count as a match. All 899,836 molecules remain in the denominator.':'Includes single-atom NPE units, but only fragments with at least one external connection count as a match. All 899,836 molecules remain in the denominator.';
    }
    $('mc-size').addEventListener('change',renderSharing);renderSharing();
    const select=$('mc-motif');
    function inspect(){
      const m=d.motifs[Number(select.value)];if(!m)return;
      $('mc-key').textContent=m.key;
      $('mc-group-info').textContent=`Canonical analysis group · ${m.atoms} real atoms · ${fmt(m.molecules)} molecules · ${fmt(m.occurrences)} fragment occurrences. Counts below include variants without external connections.`;
      $('mc-id-list').replaceChildren(...m.ids.map(i=>text('code',`NPE ID ${i.id}: ${i.token} · ${fmt(i.occurrences)} occurrences`)));
      $('mc-literal').textContent=fmt(m.literal);$('mc-normalized').textContent=fmt(m.normalized);$('mc-ports').textContent=fmt(m.ports);
      $('mc-id-note').textContent=m.ids.length>1?'This group combines '+m.ids.length+' original IDs only for analysis. The trained vocabulary and IDs are unchanged; these are not variant counts for each individual ID.':'This canonical analysis group contains one observed original NPE ID.';
    }
    function filter(){
      const q=$('mc-search').value.trim().toLowerCase();
      const keep=select.value;select.replaceChildren();
      d.motifs.forEach((m,i)=>{if(!q||m.key.toLowerCase().includes(q)||m.ids.some(v=>String(v.id)===q||v.token.toLowerCase().includes(q))){const o=text('option',`${m.key} · ${fmt(m.literal)} literal variants`);o.value=i;select.append(o);}});
      $('mc-match-count').textContent=`${fmt(select.options.length)} / 2,740 canonical analysis groups. Search by motif SMILES or observed NPE ID.`;
      if([...select.options].some(o=>o.value===keep))select.value=keep;
      $('mc-explorer-detail').hidden=!select.options.length;
      inspect();
    }
    $('mc-search').addEventListener('input',filter);select.addEventListener('change',inspect);filter();
    const cco=d.motifs.findIndex(m=>m.key==='CCO');if(cco>=0){select.value=cco;inspect();}
    for(const r of v.corpus_growth)$('mc-growth-body').append(row([fmt(r.molecules),fmt(r.unique_labelled_literal),fmt(r.unique_labelled_ports)]));
    for(const r of v.top_shared_multiatom_labelled_fragments.slice(0,5))$('mc-shared-examples').append(row([r.safe_fragment,fmt(r.molecule_count),pct(r.molecule_prevalence),r.examples.join(', ')]));
    $('mc-content').hidden=false;$('mc-loading').hidden=true;
  }
  start().catch(e=>{$('mc-loading').textContent='The census could not be loaded. Please reload, or open the linked JSON results. '+e.message;$('mc-loading').className='mc-error';});
})();
