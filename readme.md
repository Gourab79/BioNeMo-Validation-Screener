# 🔬 BioNeMo Validation Screener

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![NVIDIA NIM](https://img.shields.io/badge/NVIDIA-NIM-76b900?style=for-the-badge&logo=nvidia&logoColor=white)](https://build.nvidia.com/)
[![Flask](https://img.shields.io/badge/Flask-Web_App-000000?style=for-the-badge&logo=flask&logoColor=white)](https://flask.palletsprojects.com/)
[![RDKit](https://img.shields.io/badge/Cheminformatics-RDKit-orange?style=for-the-badge)](https://www.rdkit.org/)

> **BioNeMo Validation Screener** is a professional-grade Virtual Screening Studio that leverages AI to accelerate early-stage drug discovery. The application combines **molecular similarity search (RDKit)** with high-fidelity **AI molecular docking (NVIDIA DiffDock)** to simulate and validate how promising drug candidates fit into target proteins.

---

## 📽️ Live Demo & Preview

See the application in action. The video below demonstrates the full pipeline: from inputting a chemical query to visualizing the final 3D protein-ligand binding pose.

[![Watch the Demo](https://img.shields.io/badge/▶_Watch_Live_Demo-Click_Here-red?style=for-the-badge&logo=youtube)](https://github.com/Gourab79/BioNeMo-Validation-Screener/blob/main/Screen%20Recording%202025-12-02%20221304.mp4)

*(Click the badge above to play the video)*

---

## ✨ Key Features

* **🧪 AI Molecular Screening:** Rapid similarity search against a local `drugs.json` database using RDKit's chemical fingerprints (Tanimoto Similarity).
* **🔗 NVIDIA DiffDock Integration:** Validates top candidates by running complex 3D docking simulations via NVIDIA Inference Microservices (NIM).
* **🧬 Full Pipeline Simulation:** Seamless workflow from pasting a PDB target $\to$ molecular search $\to$ final 3D visualization.
* **📊 Professional UI/UX:** Clean, dark-mode interface with "Glassmorphism" design, real-time docking feedback, and high-fidelity molecular rendering (NGL.js).

---

## 🛠️ Project Setup & Installation

### Prerequisites
* **Python 3.8+**
* **Git**
* **NVIDIA API Key** (Required for the DiffDock AI model). [Get one here](https://build.nvidia.com/mit/diffdock).

### 1. Clone & Configure
```bash
# Clone the repository
git clone [https://github.com/Gourab79/BioNeMo-Validation-Screener.git](https://github.com/Gourab79/BioNeMo-Validation-Screener.git)
cd BioNeMo-Validation-Screener

# Create and activate virtual environment
python -m venv venv
source venv/bin/activate  # On Windows use: .\venv\Scripts\activate
```
🧬 Step-by-Step Guide: Finding a PDB Structure
You need the PDB ID or the protein's full name to search the Protein Data Bank (PDB), which is the global repository for 3D biological structures.

Step 1: Identify Your Target
First, determine the name of the protein you want to study.

Example (for your project): You need the structure for the enzyme Cyclooxygenase-2 (COX-2), ideally bound to an inhibitor.

Action: Choose a PDB ID, such as 5KIR.

Step 2: Search the PDB Database
You will use the official website to locate the structure file.

Go to RCSB PDB: Navigate to the main repository for PDB files (e.g., https://www.rcsb.org/).

Enter ID: Type the PDB ID (e.g., 5KIR) into the search bar and hit Enter.

Step 3: Download the File
Once you are on the structure's individual page, you can get the coordinate data.

Find the Download Option: Look for a button or link labeled "Download Files" (or similar).

Select Format: Choose the PDB Format (which will be a file named something like 5KIR.pdb).

Step 4: Obtain the PDB Text (Crucial for the Pipeline)
Your web application requires the raw text of the PDB file to be pasted into the input box, not the file itself.

Open the File: Open the downloaded .pdb file using any plain text editor (Notepad, VS Code, TextEdit, etc.). 

Copy All Content: Select and copy ALL the text in the file. The text will look like a long list of coordinates, starting with records like HEADER and having many ATOM records.

Paste into Dashboard: Paste the entire copied block into the "Target Protein Structure (PDB Text)" box in your BioAI dashboard.


The application now has the necessary 3D coordinates to begin the DiffDock simulation.
