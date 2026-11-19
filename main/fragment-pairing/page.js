const el=id=>document.getElementById(id);
const schemes=['Original labels','Relative pairing','Relative pairing · reordered'];
const topics=['One open site','Close and open in one token','Choose the second pending site','Internal ring and external double bond','Two occurrences of the same methyl fragment'];
let data,molecule,steps,index=0;
const lex=/\{\+\}|\{-\d+\}|\[[^\]]+\]|Cl|Br|%\(\d+\)|%\d{2}|./g;
function replay(ids,dictionary,atomOrder,fragmentOrder){
  return replayFragments(ids.map(i=>dictionary[i]),atomOrder,fragmentOrder);
}
function replayFragments(parts,atomOrder,fragmentOrder){
  const reserved=new Set(parts.flatMap(s=>s.match(lex)).filter(s=>/^(\d|%)/.test(s)).map(s=>Number(s.replace(/\D/g,''))));
  const used=new Set(reserved),pending=[],output=[],result=[];
  let nextAtom=0;
  for(const [fi,part] of parts.entries()){
    let atom=null,text='',order='';
    const before=pending.slice(),atoms=[],actions=[],connections=[],branches=[],listSizes=[];
    for(const word of part.match(lex)){
      if(word==='{+}' || word.startsWith('{-')){
        let label;
        if(word==='{+}'){
          label=1;while(used.has(label))label++;
          used.add(label);pending.push({atom,fragment:fragmentOrder?.[fi] ?? fi+1,label,order});
          actions.push('Open a site at atom '+atom+' (rank '+pending.length+').');
        }else{
          const rank=Number(word.slice(2,-1));if(rank<1||rank>pending.length)throw Error('Invalid closing rank');
          const site=pending.splice(rank-1,1)[0];label=site.label;used.delete(label);
          const bond=order||site.order;
          connections.push([site.atom,atom,bond]);
          actions.push('Close rank '+rank+': join atom '+atom+' to atom '+site.atom+' by a '+({'=':'double','#':'triple'}[bond]||'single')+' bond.');
        }
        listSizes.push(pending.length);
        text+=label<10?label:label<100?'%'+label:'%('+label+')';order='';
      }else{
        text+=word;
        if(/^(\[|Cl$|Br$|[BCNOPSFIbcnops*]$)/.test(word)){
          atom=atomOrder?.[nextAtom] ?? nextAtom;nextAtom++;atoms.push(atom);order='';
        }
        else if(word==='(')branches.push(atom);
        else if(word===')')atom=branches.pop();
        else if(/^[=#:/\\-]$/.test(word))order=word;
        else order='';
      }
    }
    output.push(text);result.push({before,after:pending.slice(),atoms,actions,connections,listSizes,text:output.join('.')});
  }
  if(pending.length)throw Error('Unclosed connection');
  return result;
}
function shape(tag,attrs){
  const n=document.createElementNS('http://www.w3.org/2000/svg',tag);
  for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);
  return n;
}
function drawStep(svg,ex,s){
  svg.replaceChildren();
  const active=new Set(s.atoms);
  for(const [a,b] of ex.bonds){
    if(!active.has(a)||!active.has(b))continue;
    const [x1,y1]=ex.points[a],[x2,y2]=ex.points[b];
    svg.append(shape('line',{x1,y1,x2,y2,stroke:'#f6c460','stroke-width':12,'stroke-linecap':'round'}));
  }
  for(const a of s.atoms){
    const [cx,cy]=ex.points[a];svg.append(shape('circle',{cx,cy,r:13,fill:'#ffdb76',stroke:'#ab6c10','data-atom':a}));
  }
  for(const [a,b] of s.connections){
    const [x1,y1]=ex.points[a],[x2,y2]=ex.points[b];
    svg.append(shape('line',{x1,y1,x2,y2,stroke:'#d76145','stroke-width':9,'stroke-linecap':'round','data-bond':a+'-'+b}));
  }
  svg.append(shape('image',{href:ex.id+'.svg',x:0,y:0,width:540,height:290}));
}
function render(){
  const ex=molecule,s=steps[index];
  drawStep(el('diagram'),ex,s);
  for(const [i,b] of Array.from(el('tokens').children).entries())b.setAttribute('aria-pressed',String(i===index));
  el('step').textContent='Token '+(index+1)+' of '+steps.length;
  el('actions').textContent=s.actions.join(' ')||'Emit this fragment.';
  for(const key of ['before','after']){
    el(key).textContent=s[key].map((v,i)=>'#'+(i+1)+' '+ex.symbols[v.atom]+'('+v.atom+') · F'+v.fragment).join('   |   ')||'Empty';
  }
  el('emitted').textContent=s.text;
  el('outcome').textContent=index===steps.length-1?'All sites paired. Recovered canonical isomeric SMILES (verified with RDKit): '+ex.smiles:'';
  el('prev').disabled=index===0;el('next').disabled=index===steps.length-1;
}
function choose(){
  molecule=data.examples[Number(el('molecule').value)];
  steps=replay(molecule.ids,data.dictionary);index=0;
  if(steps.at(-1).text!==molecule.unified)throw Error('Replay mismatch');
  el('original').textContent=molecule.safe;el('tokens').replaceChildren();
  molecule.ids.forEach((id,i)=>{
    const b=document.createElement('button');b.textContent=data.dictionary[id];
    const sub=document.createElement('small');sub.textContent='F'+(i+1);b.append(sub);
    b.onclick=()=>{index=i;render();};el('tokens').append(b);
  });render();
}
function row(table,cells){
  const tr=table.insertRow();for(const v of cells)tr.insertCell().textContent=v;
}
function setupComparison(ex,pairData){
  const decks=[];
  el('comparison-decks').replaceChildren();
  el('event-molecule').textContent=ex.id;
  el('compare-source').textContent=ex.safe;
  el('compare-smiles').textContent=ex.smiles;
  ex.runs.forEach((run,side)=>{
    const pairs=side===1?pairData.sequences[ex.id]:null;
    const parts=pairs?pairs.map(([a,b])=>pairData.motifs[a].variants[b]):run.ids.map(id=>run.dictionary[id]);
    const frames=replayFragments(parts,run.atomOrder,run.order);
    const sizes=frames.flatMap(frame=>frame.listSizes);
    if(Math.max(0,...sizes)!==run.peak)throw Error('Pending peak mismatch');
    el(side===0?'original-events':'reordered-events').textContent=sizes.join(' → ');
    if(frames.at(-1).text!==run.unified)throw Error('Fragment order replay mismatch');
    const panel=el('order-template').content.firstElementChild.cloneNode(true);
    const at=role=>panel.querySelector('[data-role="'+role+'"]');
    at('title').textContent=side===0?'Original SAFE order':'Reordered · [a,b] pairs';
    at('lookup').hidden=!pairs;
    at('peak').textContent='Peak '+run.peak+' pending';
    at('final').textContent=frames.at(-1).text;
    parts.forEach((part,i)=>{
      const btn=document.createElement('button');
      btn.type='button';btn.textContent=part;
      const meta=document.createElement('small');
      meta.textContent='F'+run.order[i]+(pairs?' · '+JSON.stringify(pairs[i]):'');
      btn.append(meta);btn.onclick=()=>selectFragment(run.order[i]);at('chips').append(btn);
    });
    decks.push({run,frames,at,pairs,parts});
    el('comparison-decks').append(panel);
  });
  function selectFragment(fragment){
    for(const {run,frames,at,pairs,parts} of decks){
      const i=run.order.indexOf(fragment),frame=frames[i];
      for(const [j,button] of Array.from(at('chips').children).entries()){
        button.setAttribute('aria-pressed',String(j===i));
      }
      at('step').textContent='F'+fragment+' · step '+(i+1)+' of '+run.ids.length;
      if(pairs){
        const [a,b]=pairs[i],motif=pairData.motifs[a];
        at('raw').textContent=ex.safe.split('.')[fragment-1];
        at('motif-label').textContent='a = '+a+' · original NPE motif';
        at('motif').textContent=motif.motif;
        at('variant-label').textContent='b = '+b+' · lookup in '+motif.count+' local variants';
        at('variant').textContent=parts[i];
      }
      at('actions').textContent=frame.actions.join(' ');
      for(const key of ['before','after']){
        at(key).textContent=frame[key].map((site,j)=>'#'+(j+1)+' '+ex.symbols[site.atom]+'('+site.atom+') · F'+site.fragment).join(' | ')||'Empty';
      }
      at('partial').textContent=frame.text;
      if(pairs && i===frames.length-1){
        at('actions').textContent+=' All sites paired. This complete SAFE string was verified against the original canonical isomeric SMILES shown below.';
      }
      drawStep(at('graph'),ex,frame);
    }
  }
  const firstPeak=decks[0].frames.findIndex(frame=>frame.listSizes.includes(ex.runs[0].peak));
  selectFragment(ex.id==='CHEMBL43138'?3:ex.runs[0].order[Math.max(0,firstPeak)]);
} 
async function start(){
  const res=await fetch('data.json?v=exact45402');if(!res.ok)throw Error(res.status);
  data=await res.json();const c=data.comparison,fmt=n=>n.toLocaleString('en-US');
  schemes.forEach((name,i)=>{
    row(el('rows'),[name,fmt(c.vocab[i]),c.mean.toFixed(4),c.median,fmt(c.strict[i])+' / '+fmt(c.n)]);
    row(el('stats'),[name,c.stats[i][0].toFixed(2),c.stats[i][1],fmt(c.stats[i][2])]);
  });
  data.examples.forEach((ex,i)=>el('molecule').add(new Option(ex.id+' · '+topics[i],i)));
  el('molecule').onchange=choose;el('prev').onclick=()=>{index--;render();};el('next').onclick=()=>{index++;render();};choose();
  const response=await fetch('reorder.json?v=exact45402');
  if(!response.ok)throw Error(response.status);
  const reordered=await response.json();
  const examples=[reordered.example,...reordered.additionalExamples];
  const picker=el('order-molecule');
  for(const ex of examples){
    const [original,ordered]=ex.runs;
    picker.add(new Option(ex.id+' · '+original.ids.length+' fragments · peak '+original.peak+' → '+ordered.peak,ex.id));
  }
  const searchResponse=await fetch('search-traces.json');
  if(!searchResponse.ok)throw Error(searchResponse.status);
  const traces=await searchResponse.json();
  const pairResponse=await fetch('pair-examples.json');
  if(!pairResponse.ok)throw Error(pairResponse.status);
  const pairData=await pairResponse.json();
  picker.onchange=()=>{
    const ex=examples[picker.selectedIndex];
    setupComparison(ex,pairData);
    setupSearch(ex,traces[ex.id]);
  };
  picker.onchange();
}
start().catch(e=>el('error').textContent='Unable to load this experiment: '+e.message);
