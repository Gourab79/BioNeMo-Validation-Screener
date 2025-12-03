import os
import json
import requests
from flask import Flask, request, jsonify, send_from_directory
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

app = Flask(__name__)

# --- CONFIGURATION ---
API_KEY = "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"
DATABASE_FILE = "drugs.json"

# Only the DiffDock model is required now
DIFFDOCK_URL = "https://health.api.nvidia.com/v1/biology/mit/diffdock/predict"

# --- RDKIT HELPERS ---

def smiles_to_sdf_block(smiles):
    """Converts a 1D chemical string into a 3D block (SDF) for NGL."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToMolBlock(mol)
    except:
        return None

# --- LOAD VECTOR DATABASE ---
drug_db = []
try:
    with open(DATABASE_FILE, 'r') as f:
        drug_db = json.load(f)
        for drug in drug_db:
            mol = Chem.MolFromSmiles(drug.get('smiles', ''))
            if mol:
                # Pre-calculate Morgan Fingerprints (Vectors)
                drug['fp'] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    print(f"✅ Loaded {len(drug_db)} drugs for screening.")
except Exception as e:
    print(f"⚠️ Error loading database: {e}")
    print("Please ensure drugs.json exists and is valid.")

# --- ROUTES ---

@app.route("/")
def index():
    return send_from_directory(".", "index.html")

@app.route("/screen", methods=["POST"])
def screen_database():
    """Endpoint for Molecular Similarity Search (RDKit)."""
    query_smiles = request.json.get("query_smiles", "")
    if not query_smiles: 
        return jsonify({"error": "Query SMILES required"}), 400

    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if not query_mol: return jsonify({"error": "Invalid SMILES"}), 400
        
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=1024)
        
        results = []
        for drug in drug_db:
            if 'fp' in drug:
                # Calculate Tanimoto Similarity (0.0 to 1.0)
                score = DataStructs.TanimotoSimilarity(query_fp, drug['fp'])
                if score > 0.2:
                    results.append({
                        "name": drug['name'],
                        "smiles": drug['smiles'],
                        "desc": drug['desc'],
                        "similarity": round(score, 4)
                    })
        
        results.sort(key=lambda x: x['similarity'], reverse=True)
        return jsonify({"matches": results[:7]}) # Return top 7

    except Exception as e:
        return jsonify({"error": str(e)}), 500
@app.route("/predict/diffdock", methods=["POST"])
def predict_diffdock():
    if not API_KEY:
        return jsonify({"error": "Server missing API key"}), 500

    data = request.json or {}
    pdb_text = data.get("pdb_string", "")
    sdf_or_smiles = data.get("smiles", "")

    if not pdb_text:
        return jsonify({"error": "Missing protein PDB"}), 400

    # Convert SMILES → SDF for DiffDock
    ligand_sdf = smiles_to_sdf_block(sdf_or_smiles)
    if not ligand_sdf:
        return jsonify({"error": "Invalid ligand SMILES"}), 400

    # ---- 1. UPLOAD ASSETS TO NVIDIA NIM ----
    def upload_asset(text_data):
        assets_url = "https://api.nvcf.nvidia.com/v2/nvcf/assets"

        meta_headers = {
            "Authorization": f"Bearer {API_KEY}",
            "Content-Type": "application/json",
            "accept": "application/json",
        }

        payload = {
            "contentType": "text/plain",
            "description": "diffdock-file"
        }

        # Create upload slot
        r1 = requests.post(assets_url, headers=meta_headers, json=payload, timeout=60)
        r1.raise_for_status()
        info = r1.json()

        upload_url = info["uploadUrl"]
        asset_id   = info["assetId"]

        # Upload actual file
        s3_headers = {
            "x-amz-meta-nvcf-asset-description": "diffdock-file",
            "content-type": "text/plain"
        }

        r2 = requests.put(upload_url, data=text_data, headers=s3_headers, timeout=300)
        r2.raise_for_status()

        return asset_id

    try:
        protein_id = upload_asset(pdb_text)
        ligand_id  = upload_asset(ligand_sdf)
    except Exception as e:
        return jsonify({"error": f"Upload failed: {e}"}), 500

    # ---- 2. CALL MIT DIFFDOCK ----
    url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"

    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {API_KEY}",
        "NVCF-INPUT-ASSET-REFERENCES": f"{protein_id},{ligand_id}",
    }

    payload = {
        "ligand": ligand_id,
        "ligand_file_type": "sdf",
        "protein": protein_id,
        "num_poses": 10,
        "time_divisions": 20,
        "steps": 18,
        "save_trajectory": False,
        "is_staged": True
    }

    try:
        r = requests.post(url, headers=headers, json=payload, timeout=600)
    except Exception as e:
        return jsonify({"error": f"Network error calling NVIDIA: {e}"}), 502

    if not r.ok:
        return jsonify({"error": r.text}), r.status_code

    return jsonify(r.json())

if __name__ == "__main__":
    app.run(debug=True, port=5000)