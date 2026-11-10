import json
import unittest
from evaluate_structural_span_vocab import syntactic_units,partition_record,SpanCodec,lex,make_foundation
from types import SimpleNamespace


def codec(words):
    primitives=['C','N','O','(',')','=','1','.','[C@@H]','Cl']
    payloads=[p.encode().hex() for p in primitives]+[None,None]
    lookup={t:i for i,t in enumerate(primitives)}
    for word in words:
        if word not in lookup:
            lookup[word]=len(payloads)
            payloads.append(word.encode().hex())
    return SpanCodec(dict(K=len(payloads),content_ids=lookup,
        primitive_ids={t:i for i,t in enumerate(primitives)},begin_id=10,end_id=11,
        dot_id=7,payloads_hex=payloads))


class StructuralSpanTests(unittest.TestCase):
    def test_preexisting_rare_atom_needs_no_new_training_entry(self):
        base=SimpleNamespace(payloads=[b'C',b'[Mg+2]',b'.'],base_by_text={'C':0,'[Mg+2]':1,'.':2})
        payloads,lookup,primitives,begin,extras=make_foundation(base,{'C','.'})
        self.assertEqual(primitives['[Mg+2]'],1)
        self.assertEqual(extras,[])

    def test_branches_atoms_and_extended_labels_remain_whole(self):
        units,_=syntactic_units('C123CCC(C1)C45C26')
        self.assertEqual([u['text'] for u in units],['C123','C','C','C(C1)','C45','C26'])
        units,_=syntactic_units('[C@@H]%(100)(Cl)C=N')
        self.assertEqual([u['text'] for u in units],['[C@@H]%(100)(Cl)','C','=N'])
        self.assertEqual(lex('%(100)'),('%(100)',))

    def test_rings_amide_and_ester_protected(self):
        parts,stats=partition_record('CC(=O)NCCOC(=O)c1ccccc1')
        self.assertIn('C(=O)N',parts['B'][0])
        self.assertIn('OC(=O)',parts['B'][0])
        self.assertIn('c1ccccc1',parts['B'][0])
        self.assertEqual(stats['protected_ring'],1)
        self.assertEqual(stats['protected_amide'],1)
        self.assertEqual(stats['protected_carboxylic_acid_or_ester'],1)
        self.assertEqual(''.join(parts['A'][0]),''.join(parts['B'][0]))

    def test_cross_fragment_ring_does_not_change_original_boundaries(self):
        parts,stats=partition_record('C12CC.C1CC2')
        self.assertEqual(stats['cross_fragment_ring'],1)
        self.assertEqual(len(parts['B']),2)
        self.assertEqual([''.join(p) for p in parts['B']],['C12CC','C1CC2'])

    def test_dynamic_programming_beats_greedy_longest(self):
        c=codec(['CC','CNO'])
        # Greedy CC + N + O costs 3; global optimum C + CNO costs 2.
        ids=c.encode_fragment(['C','C','N','O'])
        self.assertEqual(len(ids),2)
        self.assertEqual(c.decode(ids),'CCNO')

    def test_fallback_frame_preserves_whole_branch_and_reloads(self):
        c=codec([])
        ids=c.encode_fragment(['C(N)','O'])
        self.assertEqual(ids.count(c.begin),1)
        self.assertEqual(len(ids),7)
        restored=SpanCodec(json.loads(json.dumps(c.config)))
        self.assertEqual(restored.decode(ids),'C(N)O')
        with self.assertRaises(ValueError):
            c.decode([c.begin,c.primitives['C'],c.primitives['('],c.end])
        with self.assertRaises(ValueError):
            c.decode([c.begin,c.primitives['C']])


if __name__=='__main__':
    unittest.main()
