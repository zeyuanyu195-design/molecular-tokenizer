(() => {
  const ns = 'http://www.w3.org/2000/svg';
  const colors = ['#4da49d','#db977a','#9e91cb','#90ad61','#d08eaf','#6ea6cc','#d0af61','#83aa88'];
  const cache = new Map();
  const byId = id => document.getElementById(id);
  let current, request = 0;

  function svg(tag, attrs) {
    const element = document.createElementNS(ns, tag);
    for (const [key, value] of Object.entries(attrs)) element.setAttribute(key, value);
    return element;
  }

  function drawPanel(mode, graph, budget, fragment, view) {
    const panel = byId('graph-' + mode);
    const groups = graph.layouts[budget][mode].flatMap((record, fi) =>
      fragment < 0 || fi === fragment ? record[view].map(row => ({...row, fi})) : []);
    const owner = new Map();
    groups.forEach((group, i) => group.atoms.forEach(atom => owner.set(atom,i)));
    const fragmentOf = new Map();
    graph.fragments.forEach((atoms,i) => atoms.forEach(a => fragmentOf.set(a,i)));
    const figure = svg('svg', {viewBox:'0 0 720 440', role:'img','aria-label':'Scheme '+mode+' molecular graph and group boundaries'});
    let radius = Math.min(11, Math.max(4, Math.hypot(...graph.points[graph.bonds[0][0]].map((v,i)=>v-graph.points[graph.bonds[0][1]][i]))*.29));
    if(fragment >= 0){
      const points = graph.fragments[fragment].map(a=>graph.points[a]);
      const xs=points.map(p=>p[0]), ys=points.map(p=>p[1]);
      const pad=radius*4;
      figure.setAttribute('viewBox',[Math.min(...xs)-pad,Math.min(...ys)-pad,Math.max(...xs)-Math.min(...xs)+2*pad,Math.max(...ys)-Math.min(...ys)+2*pad].join(' '));
    }
    const halos = [];
    for(const [atom,group] of owner){
      const [cx,cy]=graph.points[atom];
      const halo=svg('circle',{cx,cy,r:radius,fill:colors[group%colors.length],opacity:.58});
      figure.append(halo); halos.push([halo,group]);
    }
    figure.append(svg('image',{x:0,y:0,width:720,height:440,href:'data:image/svg+xml;charset=utf-8,'+encodeURIComponent(graph.svg)}));
    const cuts=[],external=[];
    for(const [a,b] of graph.bonds){
      if(!owner.has(a)&&!owner.has(b))continue;
      if(owner.get(a)===owner.get(b))continue;
      const cross=fragmentOf.get(a)!==fragmentOf.get(b);
      (cross?external:cuts).push(a+'–'+b);
      const [x,y]=graph.points[a], [u,v]=graph.points[b];
      const midX=(x+u)/2, midY=(y+v)/2, len=Math.hypot(u-x,v-y);
      const dx=(v-y)/len*radius*.85,dy=(x-u)/len*radius*.85;
      const line=svg('line',{x1:midX-dx,y1:midY-dy,x2:midX+dx,y2:midY+dy,stroke:cross?'#555b66':'#c44a14','stroke-width':cross?2:3,'stroke-dasharray':cross?'3 2':'none'});
      figure.append(line);
    }
    const heading=document.createElement('h4'); heading.textContent='Scheme '+mode;
    const caption=document.createElement('p');
    const unit=view==='allowed'?'allowed block':'encoded group';
    const countIds=groups.reduce((total,row)=>total+row.ids.length,0);
    const ids=view==='encoded'?` · ${countIds} fragment token ID${countIds===1?'':'s'}`:'';
    caption.textContent=`${groups.length} ${unit}${groups.length===1?'':'s'}${ids} · ${cuts.length} within-fragment boundary bond${cuts.length===1?'':'s'}`;
    const chooser=document.createElement('select');
    chooser.setAttribute('aria-label', 'Highlight a group in scheme '+mode);
    chooser.add(new Option('Show all groups',''));
    groups.forEach((row,i)=>chooser.add(new Option(`F${row.fi+1} · group ${i+1} · ${row.text}`,String(i))));
    const detail=document.createElement('p');detail.className='gx-detail';
    function highlight(){
      const id=chooser.value===''?-1:Number(chooser.value);
      for(const [circle,group] of halos)circle.setAttribute('opacity',id<0||id===group?'.65':'.10');
      const row=groups[id];
      detail.textContent=row ? `${row.text} · atoms ${row.atoms.join(', ')}${row.ids.length?' · IDs '+row.ids.join(', '):''}` : '';
    }
    chooser.onchange=highlight;
    const legend=document.createElement('details');
    const summary=document.createElement('summary');summary.textContent='Atom groups and boundary bonds';
    const list=document.createElement('ol');
    groups.forEach((row,i)=>{
      const item=document.createElement('li');
      const button=document.createElement('button');button.type='button';button.style.borderLeft='5px solid '+colors[i%colors.length];
      button.textContent=`F${row.fi+1} · ${row.text} · atoms ${row.atoms.join(', ')}${row.ids.length?' · IDs '+row.ids.join(', '):''}`;
      button.onclick=()=>{chooser.value=String(i);highlight();};item.append(button);list.append(item);
    });
    const bonds=document.createElement('p');bonds.textContent='Within-fragment boundaries: '+(cuts.join(', ')||'none')+'. SAFE inter-fragment bonds: '+(external.join(', ')||'none')+'.';
    legend.append(summary,bonds,list);
    panel.replaceChildren(heading,caption,figure,chooser,detail,legend);
  }
  function render(){
    if(!current)return;
    const {graph,budget,fragment}=current;
    const view=byId('gx-view').value;
    for(const mode of ['A','B'])drawPanel(mode,graph,budget,fragment,view);
    const retained=graph.protected.filter(([,atoms])=>graph.fragments.some((f,i)=>(fragment<0||fragment===i)&&atoms.every(a=>f.includes(a))));
    byId('gx-protected').textContent='B protection inside the displayed fragments: '+(retained.map(([name,atoms])=>name+' ('+atoms.join(', ')+')').join(' · ')||'none')+'.';
  }
  byId('gx-view').onchange=render;
  window.showExampleGraphs=async function(example,budget,fragment){
    const version=++request;current=null;
    const status=byId('gx-status');status.textContent='Loading molecule…';
    byId('gx-panels').hidden=true;byId('gx-protected').textContent='';
    try{
      if(!cache.has(example.chembl_id)){
        cache.set(example.chembl_id,fetch('example-graphs/'+example.chembl_id+'.json').then(r=>{if(!r.ok)throw Error(r.status);return r.json();}));
      }
      const graph=await cache.get(example.chembl_id);
      if(version!==request)return;
      current={graph,budget,fragment};render();status.textContent='';byId('gx-panels').hidden=false;
    }catch(e){cache.delete(example.chembl_id);if(version===request)status.textContent='Could not load molecular diagram: '+e.message;}
  };
})();
