(async () => {
  const el = name => document.getElementById('threshold-' + name);
  const fmt = value => value.toLocaleString('en-US');
  try {
    const response = await fetch('threshold-results.json');
    if (!response.ok) throw new Error('Could not load the threshold results.');
    const {rows, examples} = await response.json();
    const baseline = rows.find(row => row.threshold === 1).dictionary_entries;
    for (const row of rows) {
      const saved = baseline - row.dictionary_entries;
      const tr = el('rows').insertRow();
      const cells = ['≥' + row.threshold, fmt(row.complete_entries), fmt(row.scaffold_entries),
        fmt(row.dictionary_entries), saved ? fmt(saved) + ' (' + (100 * saved / baseline).toFixed(1) + '%)' : 'Baseline',
        row.direct_hit_percent.toFixed(2) + '%', row.flat_mean.toFixed(3)];
      for (const text of cells) tr.insertCell().textContent = text;
      el('level').add(new Option('≥' + row.threshold + ' occurrences', row.threshold));
    }
    for (const e of examples) {
      el('case').add(new Option(e.cid + ' · variant seen ' + e.frequency + (e.frequency === 1 ? ' time' : ' times'), e.cid));
    }
    function render() {
      const e = examples.find(item => item.cid === el('case').value);
      const t = Number(el('level').value);
      const state = e.levels[t];
      const fallback = Array.isArray(state.b);
      el('location').textContent = e.cid + ' · fragment ' + e.position + ' of ' + e.motif_count;
      el('graph').src = 'threshold-' + e.cid + '.svg';
      el('graph').alt = 'Structure of ' + e.cid + ' with the selected fragment highlighted in gold';
      el('word').textContent = e.word;
      el('decision').textContent = 'Observed ' + e.frequency + (e.frequency === 1 ? ' time' : ' times') +
        ' under NPE node ' + e.a + '. ' + (fallback ? 'Below the threshold of ' : 'Meets the threshold of ') + t + '.';
      el('path').textContent = fallback ? 'Fallback to scaffold and connection actions' : 'Direct lookup of the complete variant';
      el('actions').replaceChildren();
      if (fallback) {
        const [index, actions] = state.b;
        el('payload').textContent = 'Look up scaffold ' + index + ' under node ' + e.a + ' to obtain ' + state.scaffold + '.';
        for (const [offset, op, rank] of actions) {
          const closing = op.endsWith('CLOSE');
          const bond = op.slice(0, closing ? -5 : -4);
          const marker = bond + (closing ? '{-' + rank + '}' : '{+}');
          const item = document.createElement('li');
          item.textContent = 'At offset ' + offset + ', insert ' + marker + ' to ' +
            (closing ? 'close pending site ' + rank : 'open a connection') + (bond ? ' with bond ' + bond : '') + '.';
          el('actions').append(item);
        }
      } else {
        el('payload').textContent = 'Look up [a, b] = [' + e.a + ', ' + state.b + '] to recover the entire highlighted fragment.';
      }
      el('ids').textContent = '[' + state.ids.join(', ') + ']';
      el('cost').textContent = 'Selected fragment: ' + state.ids.length + ' flat ' + (state.ids.length === 1 ? 'token' : 'tokens') +
        '. Whole molecule: ' + state.flat_length + ' flat tokens across ' + e.motif_count + ' fragment positions; ' +
        state.misses + (state.misses === 1 ? ' fragment uses' : ' fragments use') + ' fallback.';
      el('safe').textContent = e.safe;
      el('smiles').textContent = e.smiles;
      el('example').hidden = false;
    }
    el('level').value = '3';
    el('case').addEventListener('change', render);
    el('level').addEventListener('change', render);
    render();
  } catch (error) {
    el('error').textContent = 'Threshold comparison unavailable. ' + error.message;
  }
})();
