function showOctalReader(state, packing, word) {
  const byId = id => document.getElementById('octal-' + id);
  const digits = ((BigInt(state.b) - BigInt(state.D)) / BigInt(state.S)).toString(8);
  const tape = byId('tape');
  tape.replaceChildren();
  byId('source').textContent = 'q = ' + digits + '₈. The decoder starts at the left. Use Next digit or select a digit to move the cursor.';
  const steps = [];
  let cursor = 0;
  function group(title) {
    const box = document.createElement('div');
    box.className = 'octal-group';
    const heading = document.createElement('strong');
    heading.textContent = title;
    box.append(heading);
    tape.append(box);
    return box;
  }
  function field(box, name, width, meaning, startNote = '') {
    const raw = digits.slice(cursor, cursor + width);
    if (raw.length !== width) throw new Error('Incomplete octal field.');
    const val = parseInt(raw, 8);
    const chunk = document.createElement('div');
    chunk.className = 'octal-field';
    const label = document.createElement('span');
    label.textContent = name;
    chunk.append(label);
    box.append(chunk);
    for (let j = 0; j < width; j++) {
      const pos = cursor++;
      const button = document.createElement('button');
      button.type = 'button';
      button.textContent = raw[j];
      button.setAttribute('aria-label', 'Digit ' + (pos + 1) + ', ' + box.firstChild.textContent + ', ' + name + ', ' + raw[j]);
      button.onclick = () => select(pos);
      chunk.append(button);
      const message = j + 1 < width
        ? 'Read ' + raw.slice(0, j + 1) + '. This field takes ' + width + ' digits; read the next digit to finish it.'
        : meaning(val, raw);
      steps.push({button, chunk, title: box.firstChild.textContent + ' · ' + name, message: (j === 0 ? startNote : '') + message});
    }
    return val;
  }
  if (digits[0] !== '1') throw new Error('Invalid octal prefix.');
  field(group('Prefix'), 'Marker', 1, () => 'Skip the leading 1. It preserves leading zeroes. Action 1 begins at the next digit.');
  let action = 0, previous = 0, recovered = '';
  while (cursor < digits.length) {
    action++;
    const box = group('Action ' + action);
    const offset = field(box, 'Offset', state.scaffold.length.toString(8).length,
      (v, raw) => 'Offset ' + raw + '₈ = ' + v + ' characters into the scaffold. Next read the action type.',
      action > 1 ? 'Action ' + action + ' starts here, immediately after the previous action. ' : '');
    const closing = field(box, 'Type', 1, v => v === 1
      ? '1 means CLOSE. After the bond, read ' + packing.rank_width + ' rank digits.'
      : '0 means OPEN. Only the bond field remains.');
    if (closing > 1) throw new Error('Invalid action type.');
    const bondCode = field(box, 'Bond', (packing.bonds.length - 1).toString(8).length, v =>
      'Bond ' + v + ' means ' + (packing.bonds[v] || 'implicit bond') + '.' +
      (closing ? ' The closing rank follows.' : ' An opening has no rank field.'));
    const rank = closing ? field(box, 'Rank', packing.rank_width,
      (v, raw) => 'Rank ' + raw + '₈ = ' + v + '. Close pending site ' + v + '.') : null;
    const marker = packing.bonds[bondCode] + (closing ? '{-' + rank + '}' : '{+}');
    recovered += state.scaffold.slice(previous,offset) + marker;
    previous = offset;
    steps[steps.length-1].message += ' Action ' + action + ' ends at digit ' + cursor + ': insert ' + marker + '.' +
      (cursor < digits.length ? ' Digit ' + (cursor+1) + ' starts the next offset.' : ' End of stream. Recovered fragment: ' + recovered + state.scaffold.slice(previous) + '.');
  }
  if (recovered + state.scaffold.slice(previous) !== word) throw new Error('Octal reader reconstruction mismatch.');
  let current = 0;
  function select(pos) {
    current = pos;
    for (const [i, step] of steps.entries()) {
      step.button.setAttribute('aria-pressed', String(i === pos));
      step.chunk.classList.toggle('current-field', step.chunk === steps[pos].chunk);
    }
    byId('position').textContent = 'Digit ' + (pos+1) + ' of ' + digits.length + ' · ' + steps[pos].title;
    byId('note').textContent = steps[pos].message;
    byId('prev').disabled = pos === 0;
    byId('next').disabled = pos === digits.length-1;
  }
  byId('prev').onclick = () => select(current-1);
  byId('next').onclick = () => select(current+1);
  select(0);
}
