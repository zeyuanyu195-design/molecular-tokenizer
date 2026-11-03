'use strict';
(async function(){
  const byId=id=>document.getElementById(id);
  try {
    const response=await fetch('figures/token-distributions.json?v=20261103-simple');
    if(!response.ok)throw new Error('Histogram data HTTP '+response.status);
    const report=await response.json();
    const show=()=>{
      const kind=byId('distribution-unit').value;
      const base='figures/token-distributions-'+kind;
      byId('distribution-image').src=base+'.png?v=simple';
      byId('distribution-image').alt='BRICS versus GraphBPE: '+(kind==='fragments'?'original fragment counts':'SAFE text token lengths')+' on 44,328 matched molecules; overlaid histogram';
      byId('distribution-image-link').href=base+'.png?v=simple';
      byId('distribution-png').href=base+'.png?v=simple';byId('distribution-pdf').href=base+'.pdf?v=simple';
      byId('distribution-description').textContent=kind==='fragments'?'One original fragment = one counting token. BRICS mean: 6.102915; GraphBPE mean: 8.412832. Both medians are 5. These counts exclude connectivity overhead.':'Same SAFE text tokenizer for both representations. BRICS+SAFE mean: 60.222320; GraphBPE+SAFE mean: 68.847861. Connectivity syntax is included.';
      const body=byId('distribution-stats').querySelector('tbody');body.replaceChildren();
      for(const [name,stats] of Object.entries(report[kind])){
        const row=document.createElement('tr');
        for(const value of [name,stats.n.toLocaleString('en-US'),stats.mean.toFixed(6),stats.median,Number(stats.p95.toFixed(2)),Number(stats.p99.toFixed(2)),stats.maximum]){
          const cell=document.createElement('td');cell.textContent=value;row.appendChild(cell);
        }body.appendChild(row);
      }
    };
    byId('distribution-unit').addEventListener('change',show);show();
  }catch(error){byId('distribution-error').hidden=false;console.error(error);}
})();
