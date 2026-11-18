const el=id=>document.getElementById(id);
const schemes=['raw','unified','relative'],names=['Original labels','Unified labels','Relative pairing'];
const topics=['One open site','Close and open in one token','Choose the second pending site','Internal ring and external double bond','Reuse one ID for two methyl groups'];
let data,molecule,steps,index=0;
const lex=/\{\+\}|\{-\d+\}|\[[^\]]+\]|Cl|Br|%\(\d+\)|%\d{2}|./g;
function replay(ids,dictionary){
  const parts=ids.map(i=>dictionary[i]);
  const reserved=new Set(parts.flatMap(s=>s.match(lex)).filter(s=>/^(\d|%)/.test(s)).map(s=>Number(s.replace(/\D/g,''))));
  const used=new Set(reserved),pending=[],output=[],result=[];
  let nextAtom=0;
  for(const [fi,part] of parts.entries()){
    let atom=null,text='',order='';
    const before=pending.slice(),atoms=[],actions=[],connections=[],branches=[];
    for(const word of part.match(lex)){
      if(word==='{+}' || word.startsWith('{-')){
        let label;
        if(word==='{+}'){
          label=1;while(used.has(label))label++;
          used.add(label);pending.push({atom,fragment:fi+1,label,order});
          actions.push('Open a site at atom '+atom+' (rank '+pending.length+').');
        }else{
          const rank=Number(word.slice(2,-1));if(rank<1||rank>pending.length)throw Error('Invalid closing rank');
          const site=pending.splice(rank-1,1)[0];label=site.label;used.delete(label);
          const bond=order||site.order;
          connections.push([site.atom,atom,bond]);
          actions.push('Close rank '+rank+': join atom '+atom+' to atom '+site.atom+' by a '+({'=':'double','#':'triple'}[bond]||'single')+' bond.');
        }
        text+=label<10?label:label<100?'%'+label:'%('+label+')';order='';
      }else{
        text+=word;
        if(/^(\[|Cl$|Br$|[BCNOPSFIbcnops*]$)/.test(word)){atom=nextAtom++;atoms.push(atom);order='';}
        else if(word==='(')branches.push(atom);
        else if(word===')')atom=branches.pop();
        else if(/^[=#:/\\-]$/.test(word))order=word;
        else order='';
      }
    }
    output.push(text);result.push({before,after:pending.slice(),atoms,actions,connections,text:output.join('.')});
  }
  if(pending.length)throw Error('Unclosed connection');
  return result;
}
function shape(tag,attrs){
  const n=document.createElementNS('http://www.w3.org/2000/svg',tag);
  for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);
  return n;
}
function render(){
  const ex=molecule,s=steps[index],svg=el('diagram');
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
  for(const [i,b] of Array.from(el('tokens').children).entries())b.setAttribute('aria-pressed',String(i===index));
  el('step').textContent='Token '+(index+1)+' of '+steps.length+' · ID '+ex.ids[index];
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
    const sub=document.createElement('small');sub.textContent='F'+(i+1)+' · ID '+id;b.append(sub);
    b.onclick=()=>{index=i;render();};el('tokens').append(b);
  });render();
}
function row(table,cells){
  const tr=table.insertRow();for(const v of cells)tr.insertCell().textContent=v;
}
async function start(){
  const res=await fetch('data.json');if(!res.ok)throw Error(res.status);
  data=await res.json();const c=data.comparison,fmt=n=>n.toLocaleString('en-US');
  schemes.forEach((scheme,i)=>{
    row(el('rows'),[names[i],fmt(c.vocabulary_sizes.all[scheme]),c.mean_tokens,c.median_tokens,fmt(c.strict_isomeric_molecule_roundtrips[scheme])+' / 50,000',fmt(c.original_SAFE_text_roundtrips[scheme])+' / 50,000']);
    row(el('stats'),[names[i],data.stats[i][0].toFixed(2),data.stats[i][1],fmt(data.stats[i][2])]);
  });
  data.examples.forEach((ex,i)=>el('molecule').add(new Option(ex.id+' · '+topics[i],i)));
  el('molecule').onchange=choose;el('prev').onclick=()=>{index--;render();};el('next').onclick=()=>{index++;render();};choose();
}
start().catch(e=>el('error').textContent='Unable to load this experiment: '+e.message);
