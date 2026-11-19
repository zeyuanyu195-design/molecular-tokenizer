function fragmentPositions(edges,order){
  // The example fragment graphs are trees. Lay out successive depth levels left to right.
  const children=new Map(order.map(f=>[f,[]]));
  for(const [a,b] of edges){children.get(a).push(b);children.get(b).push(a);}
  const visited=new Set(),positions=new Map();
  let leaf=0,depthMax=0;
  function place(f,depth){
    visited.add(f);
    depthMax=Math.max(depthMax,depth);
    const ys=[];
    for(const child of children.get(f)){
      if(!visited.has(child))ys.push(place(child,depth+1));
    }
    const y=ys.length?(ys[0]+ys.at(-1))/2:leaf++;
    positions.set(f,[depth,y]);
    return y;
  }
  for(const f of order)if(!visited.has(f))place(f,0);
  for(const [f,[d,y]] of positions){
    positions.set(f,[32+476*d/Math.max(1,depthMax),leaf===1?116:27+179*y/(leaf-1)]);
  }
  return positions;
}

function setupSearch(ex,data){
  const seed=el('search-seed'),map=el('search-map'),best=ex.runs[1].order;
  const positions=fragmentPositions(data.edges,best);
  const path=order=>order.map(f=>'F'+f).join(' → ');
  const bonds=ids=>ids.map(id=>data.edges[id].map(f=>'F'+f).join('–')).join(', ')||'∅';
  seed.replaceChildren();
  for(const trace of data.traces)seed.add(new Option('F'+trace.order[0],trace.order[0]));
  seed.value=best[0];
  let step=0;
  function paint(){
    const trace=data.traces[seed.selectedIndex],s=trace.steps[step];
    const prefix=trace.order.slice(0,step+1),current=prefix.at(-1);
    map.replaceChildren();
    for(const [id,[a,b]] of data.edges.entries()){
      const [x1,y1]=positions.get(a),[x2,y2]=positions.get(b);
      const pending=s.after.includes(id);
      map.append(shape('line',{x1,y1,x2,y2,stroke:pending?'#aa520a':'#a9b8b4','stroke-width':pending?4:2,'stroke-dasharray':pending?'8 4':'none'}));
    }
    for(const [f,[cx,cy]] of positions){
      const active=f===current,read=prefix.includes(f);
      map.append(shape('circle',{cx,cy,r:19,fill:active?'#ffe5a1':read?'#246b5c':'white',stroke:active?'#9c5909':'#548276','stroke-width':active?3:1.5}));
      const t=shape('text',{x:cx,y:cy+5,'text-anchor':'middle',fill:read&&!active?'white':'#163b35'});
      t.textContent='F'+f;map.append(t);
    }
    el('search-prev').disabled=step===0;
    el('search-next').disabled=step===trace.steps.length-1;
    el('search-step').textContent='Step '+(step+1)+' / '+trace.steps.length;
    el('search-path').textContent=path(prefix);
    el('search-queue').textContent='['+bonds(s.before)+'] → ['+bonds(s.after)+']';
    el('search-ops').textContent=s.ops.join('  ')||'No external connection';
    el('search-rule').textContent=step===0?'The first fragment is fixed by the selected start.':'When possible, only fragments touching a pending bond compete. The smallest score wins, comparing columns from left to right. Ties use the fragment text key, then the original fragment number.';
    const local=el('search-local');local.replaceChildren();
    for(const [f,...score] of s.choices)row(local,['F'+f+(f===current?' · chosen':''),...score]);
    const ranking=el('search-ranking');ranking.replaceChildren();
    for(const [seq,score] of data.ranking){
      const tags=[];
      if(path(seq)===path(best))tags.push('selected');
      if(path(seq)===path(ex.runs[0].order))tags.push('original');
      if(path(seq)===path(trace.order))tags.push('from this start');
      row(ranking,[path(seq)+(tags.length?' · '+tags.join(', '):''),...score]);
    }
    const winner=data.ranking.find(([seq])=>path(seq)===path(best));
    const base=data.ranking.find(([seq])=>path(seq)===path(ex.runs[0].order));
    el('search-result').textContent='Selected for '+ex.id+' · '+path(best)+'. Original score ('+base[1].join(', ')+') → selected score ('+winner[1].join(', ')+'). Decoding the selected IDs recovers the original canonical isomeric SMILES.';
  }
  seed.onchange=()=>{step=0;paint();};
  el('search-prev').onclick=()=>{step--;paint();};
  el('search-next').onclick=()=>{step++;paint();};
  paint();
}
