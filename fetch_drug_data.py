import requests
import json
import time

def fetch_1000_drugs():
    print("🚀 Connecting to ChEMBL Database (EBI)...")
    
    # ChEMBL API Endpoint for Approved Drugs (Phase 4)
    # We ask for 500 items, formatted as JSON
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule?max_phase=4&limit=10000&format=json"
    
    try:
        response = requests.get(url)
        data = response.json()
        
        molecules = data.get('molecules', [])
        clean_db = []
        
        print(f"✅ Downloaded raw data. Processing {len(molecules)} compounds...")

        for mol in molecules:
            # 1. Extract Name
            name = mol.get('pref_name')
            
            # 2. Extract SMILES (Skip if missing)
            structs = mol.get('molecule_structures')
            if not structs: continue
            smiles = structs.get('canonical_smiles')
            if not name or not smiles: continue
            
            # 3. Create a useful description
            # (ChEMBL IDs help if you want to look them up later)
            chembl_id = mol.get('molecule_chembl_id', 'Unknown')
            mol_type = mol.get('molecule_type', 'Small molecule')
            desc = f"FDA Approved {mol_type} (ID: {chembl_id})"
            
            entry = {
                "name": name,
                "smiles": smiles,
                "desc": desc
            }
            clean_db.append(entry)

        # Save to file
        output_file = "drugs.json"
        with open(output_file, 'w') as f:
            json.dump(clean_db, f, indent=2)
            
        print(f"🎉 Success! Saved {len(clean_db)} drugs to '{output_file}'.")
        print("You can now reload your BioNeMo Dashboard.")

    except Exception as e:
        print(f"❌ Error fetching data: {e}")

if __name__ == "__main__":
    fetch_1000_drugs()