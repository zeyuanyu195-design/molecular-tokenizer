import unittest
from evaluate_append_structural_vocab import append_codec


class AppendVocabularyTests(unittest.TestCase):
    def test_prefix_hits_budget_and_nested_vocabularies(self):
        frozen=[b'C',b'N',b'.',b'CN',None,b'\xff']
        lookup={'C':0,'N':1,'.':2,'CN':3}
        words=['CN','NN','CC','CNN','CCN','NNN']
        small=append_codec(frozen,lookup,{'C','N','.','%12'},words,5,'A','source')
        large=append_codec(frozen,lookup,{'C','N','.','%12'},words,7,'A','source')
        self.assertEqual(small.payloads[:6],frozen)
        self.assertEqual(small.config['extra_atomic_lexemes'],1)
        self.assertEqual(small.config['span_entries'],2)  # one lexeme + two markers
        self.assertEqual(large.payloads[:11],small.payloads)
        self.assertEqual(small.encode_fragment(['C','N']),[3])
        self.assertEqual(large.encode_fragment(['C','N']),[3])
        for model in (small,large):
            ids=model.encode([['C%12'],['N%12']])
            self.assertEqual(model.decode(ids),'C%12.N%12')


if __name__=='__main__':
    unittest.main()
