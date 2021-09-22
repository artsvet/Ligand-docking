# Ligand-docking
Tools  to generate and convert ligand data, automation tools for protein-ligand binding screen.

Chemistry and biology involves data with many different representations, 
either variable contextual labels or formats with multiple levels of complexity.
The molecule classes in ligands.py use pre-defined patterns for lipid scaffolds and 
arbitrary branches to generate several types of chemical data. They provide an 
API to computational chemistry libraries for visualizing, formatting, or exporting 
for further use in binding screens. 

Proteins.py sets up Protein Data Bank structure files for calculating binding energies of ligands.
VinaDock.py runs Autodock VINA using DockParser.py to batch protein/ligand docking using a table I/O.
