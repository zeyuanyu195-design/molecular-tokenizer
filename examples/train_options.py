"""Self-contained demonstration of all four options; creates no output files."""
from molecular_tokenizer import MolecularTokenizer

DATA=['CCO','CCN','CCC','CCCl','CCOC(=O)C','c1ccccc1','Cc1ccccc1','c1ccncc1']

def main():
    for fragmentation in ['brics','npe']:
        for representation in ['safe','demodiff']:
            model=MolecularTokenizer(fragmentation=fragmentation,representation=representation)
            options={}
            if fragmentation=='npe':
                options['ring_vocab_size']=2
                if representation=='safe': options['motif_vocab_size']=16
            report=model.train(DATA,128 if representation=='safe' else 20,**options)
            result=model.encode('CCO')
            print(fragmentation,representation,report)
            print(result.to_dict(),model.decode(result))

if __name__=='__main__': main()
