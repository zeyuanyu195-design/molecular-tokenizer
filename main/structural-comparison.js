/* Audited static results. Browser replay uses only real ID dictionary excerpts. */
(() => {
  'use strict';
  const root = document.getElementById('structural-comparison');
  if (!root) return;
  const el = id => document.getElementById(id);
  const make = (tag, text, cls) => { const n = document.createElement(tag); if (text !== undefined) n.textContent = text; if (cls) n.className = cls; return n; };
  const names = {hybrid: 'Original hybrid', lexical: 'No new symbol-only combinations', A: 'A · Complete syntax blocks', B: 'B · Protected structures'};
  let data, budget = '16000', fragment = 0, algorithm = 'A', functionName = 'encode', decodeMethod = 'A';
  const code = {
    train: mode => `TRAIN(${mode}, K):
    use TRAIN only (35,518 records)
    keep the previous whole-fragment ranking
    base = previous 2,136 entries
         + 170 complete TRAIN ring-label lexemes
         + BLOCK_BEGIN + BLOCK_END  # M = 2,308

    for each original SAFE fragment:
        blocks = top_level_atoms_with_complete_branches()
        keep bracket atoms and ring labels indivisible
${mode === 'B' ? '        merge blocks across in-fragment rings, amides,\n        and carboxylic-acid / ester groups\n        verify connected blocks and edges between neighbors\n' : '        verify connected blocks and edges between neighbors\n'}        count legal contiguous spans in TRAIN
        # At most 16 blocks / 32 atoms / 128 lexemes

    reserve 6,932 whole entries at 16k; 7,466 at 32k
    fill remaining slots with legal spans
    compare frequency vs estimated-saving ranking
    select by VALIDATION total length; freeze before TEST`,
    encode: mode => `ENCODE(safe, ${mode}, vocabulary):
    parse the complete SAFE to locate atoms and blocks
    for each original dot-separated fragment:
        if the full fragment has a dictionary entry:
            emit its ID; do not split it further
            continue

        blocks = complete_syntax_blocks(fragment)
${mode === 'B' ? '        blocks = coalesce_protected_groups(blocks)\n        # Only groups wholly inside this original fragment\n' : ''}
        cost[end] = 0
        for i from last block to first:
            fallback = [BLOCK_BEGIN]
                     + complete_lexeme_IDs(blocks[i])
                     + [BLOCK_END]
            cost[i] = len(fallback) + cost[i + 1]
            for each legal dictionary span blocks[i:j]:
                consider 1 + cost[j]
            keep the minimum-cost continuation
        emit the chosen sequence

    insert the original dots; preserve all SAFE labels
    # A protected logical block may use multiple model IDs`,
    decode: () => `DECODE(ids, dictionary):
    for each ID:
        if it is BLOCK_BEGIN or BLOCK_END:
            track the logical block; emit no text
        else:
            emit the exact dictionary payload

    safe = concatenate_payloads_in_order()
    molecule = SAFE.decode(safe, fix=False)
    return molecule

VERIFY_OFFLINE:
    reconstructed SAFE == original SAFE
    canonical isomeric molecule == source molecule
    every logical cut respects the selected A/B blocks
    # No per-molecule edge table is needed
    # Browser verifies text; molecular chemistry was checked offline`
  };
  function tableRows(target, cohort) {
    target.replaceChildren();
    for (const key of ['hybrid','lexical','A','B']) {
      const tr = make('tr');
      const th = make('th', names[key]); th.scope = 'row'; tr.append(th);
      for (const k of ['16000','32000']) tr.append(make('td', data.budgets[k][key][cohort].mean.toFixed(3), 'sc-number'));
      target.append(tr);
    }
  }
  function renderAlgorithm() {
    const steps = [
      ['01 · Preserve SAFE', 'Keep the original GraphBPE boundaries, dots, connection labels and stereochemistry.'],
      ['02 · Define legal blocks', algorithm === 'A' ? 'Keep each top-level atom with its complete branches and ring labels. Bracket atoms and extended labels stay whole.' : 'Start with A, then keep in-fragment rings, amides and carboxylic-acid / ester groups together. Check graph connectivity.'],
      ['03 · Build the vocabulary', 'Learn legal spans on TRAIN. Reserve the same whole-fragment counts; use VALIDATION to select the ranking.'],
      ['04 · Encode & reconstruct', 'Whole entries take priority. A shortest-path search chooses spans or framed fallback. Decode exact SAFE before reconstructing the molecule.']
    ];
    const flow = el('sc-flow'); flow.replaceChildren();
    for (const [title, text] of steps) { const card = make('article'); card.append(make('h4',title),make('p',text)); flow.append(card); }
    el('sc-code').textContent = code[functionName](algorithm);
    root.querySelectorAll('[data-sc-function]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.scFunction === functionName)));
    el('sc-algorithm-title').textContent = names[algorithm];
  }
  function chip(id, model, scheme, inside = false) {
    const entry = model.dictionary[String(id)];
    const b = make('button', undefined, `sc-token sc-${entry.kind}${inside ? ' sc-primitive' : ''}`);
    b.type = 'button';
    b.append(make('span',entry.text),make('small',`ID ${id}`));
    b.title = `${names[scheme]} · ${budget} entries · ID ${id} · ${entry.text}`;
    b.addEventListener('click', () => {
      root.querySelectorAll('.sc-token.is-selected').forEach(n => n.classList.remove('is-selected'));
      b.classList.add('is-selected');
      el('sc-inspector').textContent = `${names[scheme]} · ${Number(budget).toLocaleString()} entries · ID ${id} → ${entry.text}. ${entry.kind === 'begin' || entry.kind === 'end' ? 'Logical-block marker: counted as one model token; contributes no SAFE characters.' : entry.kind === 'dot' ? 'Original SAFE fragment separator.' : inside ? 'Complete atomic or syntax lexeme inside a fallback block; not a standalone motif.' : 'Exact payload from this scheme’s dictionary; IDs are not comparable across dictionaries.'}`;
    });
    return b;
  }
  function drawTokens(parent, ids, model, scheme) {
    let target = parent, frame = null;
    for (const id of ids) {
      const entry = model.dictionary[String(id)];
      if (entry.kind === 'begin') {
        frame = make('div',undefined,'sc-frame'); frame.append(make('small','One logical block · multiple model IDs','sc-frame-label')); parent.append(frame); target = frame;
      }
      target.append(chip(id, model, scheme, Boolean(frame) && entry.kind === 'payload'));
      if (entry.kind === 'end') { target = parent; frame = null; }
    }
  }
  function renderExample() {
    const models = data.example.budgets[budget];
    el('sc-example-summary').textContent = `CHEMBL45183 · TEST record 39,224 · 8 original SAFE fragments · ${Number(budget).toLocaleString()} entries per dictionary`;
    el('sc-selected-safe').textContent = `${fragment + 1}. ${data.example.fragments[fragment]}`;
    const cards = el('sc-cards'); cards.replaceChildren();
    const sequences = el('sc-sequences'); sequences.replaceChildren();
    for (const scheme of ['hybrid','A','B']) {
      const model = models[scheme], f = model.fragments[fragment];
      const card = make('article',undefined,`sc-card sc-card-${scheme}`);
      const heading = make('div',undefined,'sc-card-heading'); heading.append(make('h4',names[scheme]),make('strong',`${model.length} total`)); card.append(heading);
      card.append(make('p',`Selected fragment: ${f.ids.length} model tokens`,'sc-subtle'));
      if (f.logical_units) {
        const details = make('details'); details.append(make('summary',`Allowed elementary blocks (${f.logical_units.length})`),make('code',f.logical_units.join(' | '),'sc-block-list')); card.append(details);
      } else card.append(make('p','Unrestricted text BPE; substrings need not preserve branches or motifs.','sc-subtle'));
      const tokenbox = make('div',undefined,'sc-token-box'); drawTokens(tokenbox,f.ids,model,scheme); card.append(tokenbox);
      if (f.fallback_blocks) card.append(make('p',`${f.fallback_blocks} unmatched protected block(s): the frame preserves the logical unit, but its contents occupy multiple model tokens.`,'sc-frame-note'));
      cards.append(card);
      const seq = make('article',undefined,'sc-sequence'); seq.append(make('h4',`${names[scheme]} · all ${model.length} model tokens`));
      drawTokens(seq,model.ids,model,scheme); sequences.append(seq);
    }
    const a = data.budgets[budget].A.all, b = data.budgets[budget].B.all;
    el('sc-example-note').textContent = `This molecule: ${models.hybrid.length} → ${models.A.length} → ${models.B.length} tokens (hybrid → A → B). Population means at this budget: A ${a.mean.toFixed(3)}, B ${b.mean.toFixed(3)}. The example is not a population-average compression claim.`;
    el('sc-allocation').textContent = budget === '16000' ? 'A/B 16k = 2,308 foundation entries + 6,932 whole fragments + 6,760 legal spans.' : 'A/B 32k = 2,308 foundation entries + 7,466 whole fragments + 22,226 legal spans.';
    el('sc-inspector').textContent = 'Select a token to inspect its actual dictionary ID and role. Each chip is one model token, including dots and block markers.';
    el('sc-decode-range').value = '0'; renderDecode();
  }
  function renderDecode() {
    const model = data.example.budgets[budget][decodeMethod];
    const slider = el('sc-decode-range'); slider.max = String(model.length);
    const n = Math.min(Number(slider.value),model.length); slider.value = String(n);
    const text = model.ids.slice(0,n).map(i => model.dictionary[String(i)].payload || '').join('');
    el('sc-decoded').textContent = text || '(No SAFE characters emitted yet.)';
    el('sc-decode-status').textContent = n === model.length ? `All ${n} IDs decoded: ${text === data.example.safe ? 'exact original SAFE verified in this browser. Molecular stereochemistry was verified offline.' : 'ERROR: SAFE mismatch.'}` : `${n} / ${model.length} IDs decoded. Partial text is not expected to be a valid molecule.`;
  }
  async function init() {
    try {
      const response = await fetch('structural-comparison.json?v=20261110-structural-v1');
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      data = await response.json();
      for (const models of Object.values(data.example.budgets)) for (const model of Object.values(models)) {
        if (model.ids.length !== model.length || model.ids.map(i => model.dictionary[String(i)].payload || '').join('') !== data.example.safe) throw new Error('Example dictionary audit failed');
      }
      tableRows(el('sc-stat-body'),'all'); tableRows(el('sc-test-body'),'test');
      el('sc-safe').textContent = data.example.safe;
      data.example.fragments.forEach((f,i) => { const option = make('option',`${i+1}. ${f}`); option.value = String(i); el('sc-fragment').append(option); });
      el('sc-budget').addEventListener('change', e => { budget = e.target.value; renderExample(); });
      el('sc-fragment').addEventListener('change', e => { fragment = Number(e.target.value); renderExample(); });
      el('sc-algorithm').addEventListener('change', e => { algorithm = e.target.value; renderAlgorithm(); });
      root.querySelectorAll('[data-sc-function]').forEach(b => b.addEventListener('click',() => { functionName = b.dataset.scFunction; renderAlgorithm(); }));
      el('sc-decode-method').addEventListener('change', e => { decodeMethod = e.target.value; el('sc-decode-range').value = '0'; renderDecode(); });
      el('sc-decode-range').addEventListener('input',renderDecode);
      el('sc-decode-next').addEventListener('click',() => { el('sc-decode-range').value = String(Number(el('sc-decode-range').value)+1); renderDecode(); });
      el('sc-decode-all').addEventListener('click',() => { el('sc-decode-range').value = String(data.example.budgets[budget][decodeMethod].length); renderDecode(); });
      renderAlgorithm(); renderExample();
      el('sc-loading').hidden = true; el('sc-content').hidden = false;
    } catch (error) { el('sc-loading').textContent = `Unable to load structural comparison: ${error.message}. Reload or open the verified JSON below.`; el('sc-loading').setAttribute('role','alert'); }
  }
  init();
})();
