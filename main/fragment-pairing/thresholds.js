(async () => {
  const el = name => document.getElementById('threshold-' + name);
  const fmt = value => value.toLocaleString('en-US');
  try {
    const response = await fetch('threshold-results.json?v=fragment-pairs');
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
      el('case').add(new Option(e.cid, e.cid));
    }
    let index;
    function selectMolecule() {
      const e = examples.find(item => item.cid === el('case').value);
      index = e.position - 1;
      const buttons = e.fragments.map((fragment, i) => {
        const button = document.createElement('button');
        button.type = 'button';
        button.textContent = (i + 1) + ' · ' + fragment.word;
        button.addEventListener('click', () => { index = i; render(); });
        return button;
      });
      el('fragments').replaceChildren(...buttons);
      render();
    }
    function render() {
      const e = examples.find(item => item.cid === el('case').value);
      const t = Number(el('level').value);
      const fragment = e.fragments[index];
      const state = fragment.levels[t];
      const fallback = Array.isArray(state.b);
      Array.from(el('fragments').children).forEach((button, i) => button.setAttribute('aria-pressed', String(i === index)));
      el('location').textContent = e.cid + ' · fragment ' + (index + 1) + ' of ' + e.motif_count;
      el('graph').src = fragment.graph;
      el('graph').alt = 'Structure of ' + e.cid + ' with fragment ' + (index + 1) + ' highlighted in gold';
      el('word').textContent = fragment.word;
      el('decision').textContent = 'Observed ' + fragment.frequency + (fragment.frequency === 1 ? ' time' : ' times') +
        ' under NPE node ' + fragment.a + '. ' + (fallback ? 'Below the threshold of ' : 'Meets the threshold of ') + t + '.';
      el('pair').textContent = JSON.stringify([fragment.a, state.b]).replaceAll(',', ', ');
      el('node').textContent = 'a = ' + fragment.a + ' · original NPE node ID';
      el('path').textContent = fallback ? 'Fallback · b = [scaffold ID, connection actions]' : 'Direct hit · b is a complete variant ID';
      el('connections').hidden = !fallback;
      el('actions').replaceChildren();
      if (fallback) {
        const [index, actions] = state.b;
        el('payload').textContent = 'Local scaffold ID ' + index + ' selects ' + state.scaffold + '.';
        if (!actions.length) el('actions').append(Object.assign(document.createElement('li'), {textContent:'No external connections.'}));
        for (const [offset, op, rank] of actions) {
          const closing = op.endsWith('CLOSE');
          const bond = op.slice(0, closing ? -5 : -4);
          const marker = bond + (closing ? '{-' + rank + '}' : '{+}');
          const item = document.createElement('li');
          item.textContent = JSON.stringify(rank === undefined ? [offset, op] : [offset, op, rank]) + ' · At offset ' + offset + ', insert ' + marker + ' to ' +
            (closing ? 'close pending site ' + rank : 'open a connection') + (bond ? ' with bond ' + bond : '') + '.';
          el('actions').append(item);
        }
      } else {
        el('payload').textContent = 'b = ' + state.b + ' selects the complete fragment ' + fragment.word + '.';
      }
      el('ids').textContent = '[' + state.ids.join(', ') + ']';
      el('flat').open = false;
      el('flat-label').textContent = 'Expanded one-dimensional IDs · ' + state.ids.length + (state.ids.length === 1 ? ' token' : ' tokens');
      el('cost').textContent = 'Whole molecule · ' + e.motif_count + ' fragment positions · ' +
        e.levels[t].flat_length + ' tokens after expansion';
      el('safe').textContent = e.safe;
      el('smiles').textContent = e.smiles;
      el('example').hidden = false;
    }
    el('level').value = '3';
    el('case').addEventListener('change', selectMolecule);
    el('level').addEventListener('change', render);
    selectMolecule();
  } catch (error) {
    el('error').textContent = 'Threshold comparison unavailable. ' + error.message;
  }
})();
