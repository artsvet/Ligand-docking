# Ligand-docking
Tools  to generate and convert ligand data, automation tools for protein-ligand binding screen.

Chemistry and biology involves data with many different representations, 
either variable contextual labels or formats with multiple levels of complexity.
The molecule classes in lipidGen.py use pre-defined patterns for scaffolds and 
arbitrary branches to generate several types of chemical data. They provide an 
API to computational chemistry libraries for visualizing, formating, or exporting 
for further use in computations simulation.

Ligands.py and Proteins.py interface to common input formats and generate dependencies for docking.
VinaDock.py runs Autodock VINA using DockParser.py to batch docking using a table I/O.
