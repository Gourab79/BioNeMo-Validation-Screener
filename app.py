import os
import json
from flask import Flask, send_from_directory, request, jsonify

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKit_AVAILABLE = True
except Exception:
    RDKit_AVAILABLE = False


app = Flask(__name__, static_folder='.', static_url_path='')


@app.route('/')
def index():
    return send_from_directory('.', 'index.html')


def smiles_to_sdf_block(smiles: str) -> str | None:
    """Converts SMILES to a 3D mol block (SDF/MOL) using RDKit if available."""
    if not RDKit_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToMolBlock(mol)
    except Exception:
        return None


@app.route('/predict/<model>', methods=['POST'])
def predict(model: str):
    data = request.get_json(silent=True) or {}

    # Simple stubbed responses so the front-end can work immediately.
    if model in ('openfold3', 'alphafold2', 'diffdock'):
        pdb = (
            'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n'
            'ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n'
            'ATOM      3  C   ALA A   1       1.948   1.412   0.000  1.00  0.00           C\n'
        )
        return jsonify({'data': pdb, 'type': 'pdb'})

    if model == 'molmim':
        smiles = data.get('smiles') or data.get('smiles_seed') or 'Cc1ccccc1'
        molblock = smiles_to_sdf_block(smiles)
        if molblock:
            return jsonify({'data': molblock, 'type': 'sdf'})
        return jsonify({'data': f'SMILES:{smiles}', 'type': 'smi'})

    return jsonify({'error': 'unknown model'}), 400


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)