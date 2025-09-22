# Ligand-docking
Tools  to generate and convert ligand data, automation tools for protein-ligand binding screen.

-  Grab cleaned up targets using targets.Protein and varied input structures using ligands.Ligand.from_ objects
-  Make a Docker and run docks with VINA or Uni-Molv2, scoring with VINA
-  When the dock is finished, get output tables and files from _out properties for docks and ligands 

Chemistry and biology involves data with many different representations, 
either variable contextual labels or formats with multiple levels of complexity.
The Ligand docking library provides an API to computational chemistry libraries for 
visualizing, formatting, or exporting for further use in binding screens. 

The pipeline now implements both AutoDock VINA, the most widely used software for molecular docking simulation
using a genetic algorithm and united-atom scoring function, and the new state-of-the-art deep learning model 
Uni-Mol Docking v2!

Alcaide, E., Gao, Z., Ke, G., Li, Y., Zhang, L., Zheng, H., & Zhou, G. (2024). 
Uni-Mol Docking V2: Towards Realistic and Accurate Binding Pose Prediction. arXiv preprint arXiv:2405.11769.

J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli. (2021). 
AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.

O. Trott and A. J. Olson. (2010). 
AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-461.