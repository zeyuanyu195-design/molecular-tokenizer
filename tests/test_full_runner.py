"""Training-only runs must never enter audit or publication stages."""
import hashlib
import importlib
import json
from pathlib import Path
import sys

import pytest


def test_train_only_finishes_both_models_without_analysis(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / 'examples'))
    runner = importlib.import_module('run_full_experiment')
    source = tmp_path / 'source.smi'
    source.write_bytes(b'CCO\n')
    metadata = tmp_path / 'metadata.json'
    metadata.write_text(json.dumps({'exported_rows': 1, 'outputs': {
        'chembl_37_smiles.smi': {'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}}}))
    output = tmp_path / 'run'
    calls = []

    def supervise(command, directory, args):
        calls.append((command, directory.name))
        return {'exit_code': 0}

    def unexpected(*args, **kwargs):
        pytest.fail('Training-only mode attempted audit, rendering, or publication')

    monkeypatch.setattr(runner, 'supervise', supervise)
    monkeypatch.setattr(runner.subprocess, 'run', unexpected)
    monkeypatch.setattr(runner.subprocess, 'check_output', unexpected)
    monkeypatch.setattr(sys, 'argv', ['run_full_experiment.py', '--input', str(source),
        '--metadata', str(metadata), '--output', str(output), '--train-only'])
    runner.main()
    assert [name for _, name in calls] == ['npe_safe', 'brics_safe']
    assert all(command[1:4] == ['-m', 'molecular_tokenizer', 'train'] for command, _ in calls)
    assert all('--max-molecules' not in command for command, _ in calls)
    state = json.loads((output / 'run_state.json').read_text())
    assert state['scope'] == 'training_only' and state['status'] == 'complete'
    assert not (output / 'audit').exists()


def test_train_only_rejects_publish(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / 'examples'))
    runner = importlib.import_module('run_full_experiment')
    monkeypatch.setattr(sys, 'argv', ['run_full_experiment.py', '--input', 'unused',
        '--metadata', 'unused', '--output', str(tmp_path / 'unused'), '--train-only', '--publish'])
    with pytest.raises(SystemExit) as error:
        runner.main()
    assert error.value.code == 2
