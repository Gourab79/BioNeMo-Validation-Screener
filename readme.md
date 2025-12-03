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