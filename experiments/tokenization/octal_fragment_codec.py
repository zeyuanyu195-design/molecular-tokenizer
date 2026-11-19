"""Reversible integer coding of complete variants and fallback scaffold connections."""
def width(maximum):
    return len(format(maximum, 'o'))


def octal_field(value, size):
    if value < 0 or value >= 8**size:
        raise ValueError('field outside its octal width')
    return format(value, f'0{size}o')


def action_groups(scaffold, events, rules):
    groups = []
    for pos, operation, *rank in events:
        closing = operation.endswith('CLOSE')
        suffix = 'CLOSE' if closing else 'OPEN'
        if not operation.endswith(suffix) or len(rank) != int(closing):
            raise ValueError('invalid connection operation')
        bond = operation[:-len(suffix)]
        fields = [octal_field(pos, width(len(scaffold))), str(int(closing)),
                  octal_field(rules['bonds'].index(bond), width(len(rules['bonds'])-1))]
        if closing:
            if rank[0] < 1:
                raise ValueError('closing ranks start at one')
            fields.append(octal_field(rank[0], rules['rank_width']))
        groups.append(fields)
    return groups


def pack(scaffold_id, events, node, rules):
    groups = action_groups(node['scaffolds'][scaffold_id], events, rules)
    digits = '1' + ''.join(field for action in groups for field in action)
    return len(node['frequent']) + scaffold_id + len(node['scaffolds']) * int(digits, 8)


def unpack(b, node, rules):
    """Return the recovered relative-pairing word. No input molecule or per-row metadata."""
    direct = len(node['frequent'])
    if type(b) is not int or b < 0:
        raise ValueError('expected nonnegative integer b')
    if b < direct:
        return node['frequent'][b]
    q, s = divmod(b-direct, len(node['scaffolds']))
    scaffold = node['scaffolds'][s]
    digits = format(q, 'o')
    if not digits.startswith('1'):
        raise ValueError('missing leading digit')
    cursor, previous, chunks = 1, 0, []

    def read(size):
        nonlocal cursor
        field = digits[cursor:cursor+size]
        if len(field) != size:
            raise ValueError('incomplete action')
        cursor += size
        return int(field, 8)

    while cursor < len(digits):
        pos = read(width(len(scaffold)))
        closing = read(1)
        bond = rules['bonds'][read(width(len(rules['bonds'])-1))]
        if not previous <= pos <= len(scaffold) or closing not in (0, 1):
            raise ValueError('invalid position or action')
        chunks.extend([scaffold[previous:pos], bond])
        if closing:
            rank = read(rules['rank_width'])
            if not rank:
                raise ValueError('closing rank zero')
            chunks.append('{-' + str(rank) + '}')
        else:
            chunks.append('{+}')
        previous = pos
    return ''.join(chunks) + scaffold[previous:]
