import json
from molecular_tokenizer.__main__ import main


def test_cli_train_encode_evaluate(tmp_path,capsys):
    source=tmp_path/'input.smi'
    source.write_text('CCO\nCCN\nCCO\n',encoding='utf-8')
    model=tmp_path/'model'
    assert main(['train','--input',str(source),'--output',str(model),
                 '--fragmentation','npe','--vocab-size','128',
                 '--motif-vocab-size','6','--ring-vocab-size','0'])==0
    report=json.loads((model/'training_report.json').read_text())
    assert report['report']['fragmentation']=='npe'
    assert report['report']['representation']=='safe'
    assert report['report']['molecules']==3
    capsys.readouterr()
    assert main(['encode','--model',str(model),'--smiles','[13CH3]CO'])==0
    output=json.loads(capsys.readouterr().out)
    assert output['decoded']=='[13CH3]CO'
    result=tmp_path/'evaluation.json'
    source.write_text('CCO\ninvalid\n',encoding='utf-8')
    assert main(['evaluate','--model',str(model),'--input',str(source),'--output',str(result)])==0
    output=json.loads(result.read_text())
    assert output['molecules']==2 and output['strict_roundtrip_success']==1
    assert output['failures'][0]['source_row']==2
