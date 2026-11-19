(async () => {
  const el = name => document.getElementById('threshold-' + name);
  const fmt = value => value.toLocaleString('en-US');
  try {
    const response = await fetch('threshold-results.json?v=octal-actions-v1');
    if (!response.ok) throw new Error('Could not load the threshold results.');
    const {cohort, rows, examples, packing} = await response.json();
    el('bond-codes').textContent = packing.bonds.map((bond, i) => i+' = '+(bond || 'implicit bond')).join(' · ');
    el('length').textContent = 'All four thresholds · ' + fmt(cohort.tokens) + ' fragment tokens across ' +
      fmt(cohort.molecules) + ' molecules · mean ' + (cohort.tokens / cohort.molecules).toFixed(3) +
      ' tokens per molecule · median ' + cohort.median;
    const baseline = rows.find(row => row.threshold === 1).dictionary_entries;
    for (const row of rows) {
      const saved = baseline - row.dictionary_entries;
      const tr = el('rows').insertRow();
      const cells = ['≥' + row.threshold, fmt(row.complete_entries), fmt(row.scaffold_entries),
        fmt(row.dictionary_entries), saved ? fmt(saved) + ' (' + (100 * saved / baseline).toFixed(1) + '%)' : 'Baseline',
        row.direct_hit_percent.toFixed(2) + '%', fmt(row.strict_roundtrips) + ' / ' + fmt(cohort.molecules)];
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
      const fallback = BigInt(state.b) >= BigInt(state.D);
      const reader = document.getElementById('octal-reader');
      reader.hidden = !(fallback && e.cid === 'CHEMBL405225' && index === e.position - 1);
      if (!reader.hidden) showOctalReader(state, packing, fragment.word);
      Array.from(el('fragments').children).forEach((button, i) => button.setAttribute('aria-pressed', String(i === index)));
      el('location').textContent = e.cid + ' · fragment ' + (index + 1) + ' of ' + e.fragments.length;
      el('graph').src = fragment.graph;
      el('graph').alt = 'Structure of ' + e.cid + ' with fragment ' + (index + 1) + ' highlighted in gold';
      el('word').textContent = fragment.word;
      el('decision').textContent = 'Observed ' + fragment.frequency + (fragment.frequency === 1 ? ' time' : ' times') +
        ' under NPE node ' + fragment.a + '. ' + (fallback ? 'Below the threshold of ' : 'Meets the threshold of ') + t + '.';
      el('pair').textContent = '[' + fragment.a + ', ' + state.b + ']';
      el('node').textContent = 'a = ' + fragment.a + ' · original NPE node ID';
      el('path').textContent = fallback ? 'Fallback · calculate the integer b' : 'Direct hit · b is a complete variant ID';
      el('connections').hidden = !fallback;
      el('actions').replaceChildren();
      if (fallback) {
        const actions = state.actions;
        el('payload').textContent = 'Under node ' + fragment.a + ', D = ' + state.D + ' complete variants and S = ' +
          state.S + ' fallback scaffolds. Local scaffold s = ' + state.s + ' selects ' + state.scaffold + '.';
        if (!actions.length) el('actions').append(Object.assign(document.createElement('li'), {textContent:'No external connections.'}));
        for (const [i, [offset, op, rank]] of actions.entries()) {
          const closing = op.endsWith('CLOSE');
          const bond = op.slice(0, closing ? -5 : -4);
          const marker = bond + (closing ? '{-' + rank + '}' : '{+}');
          const item = document.createElement('li');
          item.textContent = JSON.stringify(rank === undefined ? [offset, op] : [offset, op, rank]) + ' · At offset ' + offset + ', insert ' + marker + ' to ' +
            (closing ? 'close pending site ' + rank : 'open a connection') + (bond ? ' with bond ' + bond : '') +
            '. Octal fields: ' + state.groups[i].join(' | ') + '.';
          el('actions').append(item);
        }
        const octal = '1' + state.groups.flat().join('');
        const q = BigInt('0o' + octal);
        const D = BigInt(state.D), S = BigInt(state.S), s = BigInt(state.s);
        if (D+s+S*q !== BigInt(state.b)) throw new Error('Packed value does not match the saved token.');
        el('pack').textContent = [
          'q = ' + octal + '₈ = ' + q + ' (decimal)',
          'b = D + s + S × q',
          '  = ' + D + ' + ' + s + ' + ' + S + ' × ' + q + ' = ' + state.b,
          'Decode: (' + state.b + ' − ' + D + ') ÷ ' + S + ' gives quotient ' + q + ' and remainder ' + s + '.'
        ].join('\n');
      } else {
        el('payload').textContent = 'The complete word ' + fragment.word + ' is entry ' + state.b +
          ' (counting from zero) among D = ' + state.D + ' retained variants of node ' + fragment.a +
          '. Lookup gives b = ' + state.b + '.';
      }
      const firstMiss = rows.find(row => fragment.levels[row.threshold].s !== undefined);
      el('show-fallback').hidden = fallback || !firstMiss;
      if (firstMiss) {
        el('show-fallback').textContent = 'See octal calculation at ≥' + firstMiss.threshold;
        el('show-fallback').onclick = () => { el('level').value = String(firstMiss.threshold); render(); };
      }
      el('cost').textContent = 'Whole molecule · ' + e.fragments.length + ' fragment tokens';
      el('safe').textContent = e.safe;
      el('smiles').textContent = e.smiles;
      el('example').hidden = false;
    }
    el('level').value = '3';
    el('case').addEventListener('change', selectMolecule);
    el('level').addEventListener('change', render);
    selectMolecule();
    if (location.hash === '#octal-reader') document.getElementById('octal-reader').scrollIntoView();
  } catch (error) {
    el('error').textContent = 'Threshold comparison unavailable. ' + error.message;
  }
})();
