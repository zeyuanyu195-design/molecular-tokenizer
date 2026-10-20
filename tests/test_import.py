import importlib
from pathlib import Path
import subprocess
import sys

import pytest


def test_import():
    assert importlib.import_module("safe") is not None
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import safe, sys; "
            "assert 'torch' not in sys.modules; "
            "assert 'transformers' not in sys.modules; "
            "assert safe.split('C%(100).N%99') == ['C', '%(100)', '.', 'N', '%99']; "
            "print(safe.__version__)",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip()

    subprocess.run(
        [
            sys.executable,
            "-c",
            "from safe import *; import sys; "
            "assert 'torch' not in sys.modules; "
            "assert callable(encode) and split('C.C') == ['C', '.', 'C']",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("module", ["batch_chembl", "batch_chembl_brics"])
@pytest.mark.parametrize("direct_script", [False, True])
def test_batch_cli_from_another_directory(module, direct_script, tmp_path):
    root = Path(__file__).resolve().parents[1]
    target = [str(root / "safe" / f"{module}.py")] if direct_script else ["-m", f"safe.{module}"]
    result = subprocess.run(
        [sys.executable, *target, "--help"],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "--input" in result.stdout
    assert "--output-dir" in result.stdout


def test_bpe_converter_imports_and_round_trips():
    import datamol as dm
    from safe import SAFEConverter

    prefix = (
        Path(__file__).resolve().parents[1] / "safe/data/tokenizer/vocab3000ring300/pretrain-token"
    )
    converter = SAFEConverter(slicer="bpe", bpe_vocab_path=prefix)
    smiles = "CCOC(=O)c1ccccc1"
    encoded = converter.encoder(smiles, allow_empty=True)
    assert dm.to_smiles(dm.to_mol(converter.decoder(encoded))) == smiles
