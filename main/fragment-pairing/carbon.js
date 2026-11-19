async function startCarbon(){
  const response=await fetch('carbon-examples.json');
  if(!response.ok)throw Error(response.status);
  const data=await response.json(),templates=data.templates,base=data.radix,N=templates.length;
  function unpack(b){
    const t=b%N,template=templates[t],ranks=[],remainders=[];
    let p=Math.floor(b/N);
    const word=template.replaceAll('{-k}',()=>{
      const k=p%base;
      if(!k)throw Error('Missing closing rank');
      ranks.push(k);remainders.push(p+' mod '+base+' = '+k);
      p=Math.floor(p/base);
      return '{-'+k+'}';
    });
    if(p)throw Error('Unused closing rank');
    return {word,t,template,ranks,remainders};
  }
  let current,stage=0;
  const choices=[];
  for(const mol of Object.values(data.molecules)){
    const words=mol.pairs.map(([a,b])=>a===2897?unpack(b).word:data.other[a+','+b]);
    const frames=replayFragments(words,mol.atomOrder,mol.order);
    if(frames.at(-1).text!==mol.decoded)throw Error('Carbon replay mismatch');
    for(const example of mol.examples){
      const i=example.index,[a,b]=mol.pairs[i],v=unpack(b);
      if(a!==2897||v.word!==example.word)throw Error('Variant mismatch');
      const button=document.createElement('button');
      const sub=document.createElement('small');
      sub.textContent=mol.id+' · F'+mol.order[i];button.textContent=v.word;button.append(sub);
      const choice={mol,i,b,v,frame:frames[i],old:example.old,button,safe:frames.at(-1).text};
      choices.push(choice);button.onclick=()=>choose(choice);el('carbon-options').append(button);
    }
  }
  function choose(choice){
    current=choice;stage=0;
    el('carbon-result').open=false;
    for(const c of choices)c.button.setAttribute('aria-pressed',String(c===choice));
    draw(); 
  }
  function draw(){
    const {mol,i,b,v,frame,old,safe}=current;
    const params='['+v.ranks.join(', ')+']';
    const terms=v.ranks.map((k,j)=>j===0?String(k):base**j+' × '+k);
    const p=v.ranks.reduce((sum,k,j)=>sum+k*base**j,0);
    if(v.t+N*p!==b)throw Error('Packed ID mismatch');
    const stages=[
      ['Read the complete carbon fragment',mol.safe.split('.')[mol.order[i]-1]+' → '+v.word,'In the reordered sequence, this is fragment F'+mol.order[i]+'. Its original local variant ID was '+old+'.'],
      ['Extract template and closing ranks','t = '+v.t+'   '+v.template+'   ranks = '+params,'Replace each closing rank with k, reading from left to right. Openings, bonds and stereochemistry keep their positions.'],
      ['Pack the parameters into b','p = '+(terms.length?terms.join(' + ')+' = ':'')+p+'\nb = '+v.t+' + '+N+' × '+p+' = '+b,'Output one position: [2897, '+b+']. The first rank occupies the lowest base-eight digit. With no closing ranks, p is zero.'],
      ['Recover the template and fill its slots','t = '+b+' mod '+N+' = '+v.t+'\np = floor('+b+' / '+N+') = '+p+'\n'+(v.remainders.length?v.remainders.join(' → ')+'\n':'')+v.template+'  →  '+v.word,'Read one base-eight digit per slot from right to left in the packed number, dividing p by eight after each digit. Insert the ranks into the template from left to right.'],
      ['Reconnect the highlighted atom',v.word+'\nPending sites: '+frame.before.length+' → '+frame.after.length,frame.actions.join(' ')+' Continuing through the complete sequence yields the decoded structure below.']
    ];
    const [title,math,note]=stages[stage];
    el('carbon-step').textContent='Step '+(stage+1)+' / '+stages.length;
    el('carbon-title').textContent=title;el('carbon-math').textContent=math;
    el('carbon-math').style.whiteSpace='pre-wrap';el('carbon-note').textContent=note;
    el('carbon-location').textContent=mol.id+' · atom '+frame.atoms.join(', ')+' · token '+(i+1)+' of '+mol.order.length;
    drawStep(el('carbon-graph'),mol,{...frame,connections:stage===4?frame.connections:[]});
    el('carbon-safe').textContent=safe;el('carbon-smiles').textContent=mol.smiles;
    el('carbon-result').open=stage===4;
    el('carbon-prev').disabled=stage===0;el('carbon-next').disabled=stage===4;
  }
  el('carbon-prev').onclick=()=>{stage--;draw();};el('carbon-next').onclick=()=>{stage++;draw();};
  choose(choices[1]);
}
startCarbon().catch(e=>el('carbon-error').textContent='Unable to load carbon examples: '+e.message);
