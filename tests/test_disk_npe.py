import pytest
from molecular_tokenizer import MolecularTokenizer, TokenizerError
from molecular_tokenizer.backends import NPEBackend
from molecular_tokenizer.disk_npe import DiskNPETrainer

DATA=['CCO','CCN','CCC','CCCC','CCCO','CCCN','CCOC(=O)C','CCN(CC)CC',
      'c1ccccc1','Cc1ccccc1','Oc1ccccc1','c1ccncc1','CCCl','CCF',
      'N[C@@H](C)C(=O)O','[13CH3]CO','CC(=O)[O-].[Na+]']*3


@pytest.mark.parametrize('batch_size',[1,7,64])
def test_disk_exactly_matches_in_memory_vocab_and_encoding(tmp_path,batch_size):
    memory=MolecularTokenizer(fragmentation='npe',representation='safe')
    disk=MolecularTokenizer(fragmentation='npe',representation='safe')
    opts=dict(motif_vocab_size=35,ring_vocab_size=2,num_workers=2)
    a=memory.train(DATA,256,**opts)
    b=disk.train(iter(DATA),256,**opts,npe_storage='sqlite',work_dir=tmp_path/'work',batch_size=batch_size)
    assert a==b
    assert memory.backend.npe.state()==disk.backend.npe.state()
    assert memory.backend.state()==disk.backend.state()
    assert memory.tokenizer_id==disk.tokenizer_id
    for s in DATA:
        assert memory.encode(s)==disk.encode(s)
    loaded=MolecularTokenizer.load(disk.save(tmp_path/'model'))
    assert loaded.encode('[13CH3]CO')==disk.encode('[13CH3]CO')


def test_resume_and_source_validation(tmp_path):
    class Interrupted(Exception): pass
    def interrupt(event):
        if event['event']=='merge': raise Interrupted()
    e=NPEBackend().engine
    with pytest.raises(Interrupted), DiskNPETrainer(tmp_path,35,2,workers=2,batch_size=5,progress=interrupt) as t:
        t.train(e,DATA)
    with DiskNPETrainer(tmp_path,35,2,workers=2,resume=True) as t:
        with pytest.raises(TokenizerError,match='source'):
            t.train(e,DATA+['CC'])
        t.train(e,DATA)
    reference=NPEBackend().engine
    reference.train_node(DATA,35,2,2)
    assert e.vocab_node_stats==reference.vocab_node_stats
    assert e.vocab_node==reference.vocab_node
    with pytest.raises(TokenizerError,match='configuration'):
        DiskNPETrainer(tmp_path,36,2,resume=True)


def test_disk_rejection_rows_are_original_indices(tmp_path):
    data=['CCO','not-smiles','CCN','*CC','[13CH3]CO']
    rejected=[]
    model=MolecularTokenizer(fragmentation='npe',representation='safe')
    report=model.train(data,128,motif_vocab_size=10,ring_vocab_size=0,npe_storage='sqlite',
                      work_dir=tmp_path,skip_invalid=True,num_workers=2,
                      on_reject_with_index=lambda s,r,i:rejected.append(i))
    assert rejected==[2,4]
    assert report.molecules==3 and report.rejected_molecules==2
    assert model.decode(model.encode('[13CH3]CO'))=='[13CH3]CO'


@pytest.mark.parametrize('phase',['remove','update'])
def test_resume_in_middle_of_merge_phase(tmp_path,phase):
    class Interrupted(Exception): pass
    data=DATA*2
    def interrupt(event):
        if event['event']=='checkpoint' and event.get('phase')==phase: raise Interrupted()
    engine=NPEBackend().engine
    with pytest.raises(Interrupted), DiskNPETrainer(tmp_path,20,2,workers=2,batch_size=1,progress=interrupt) as trainer:
        trainer.train(engine,data)
    with DiskNPETrainer(tmp_path,20,2,workers=2,batch_size=16,resume=True) as trainer:
        trainer.train(engine,data)
    reference=NPEBackend().engine
    reference.train_node(data,20,2,2)
    assert engine.vocab_node_stats==reference.vocab_node_stats
    assert engine.vocab_node==reference.vocab_node


def test_checkpoint_single_writer(tmp_path):
    with DiskNPETrainer(tmp_path,20,2):
        with pytest.raises(TokenizerError,match='Another process'):
            DiskNPETrainer(tmp_path,20,2,resume=True)
