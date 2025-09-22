import re
import os
from pathlib import Path
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdPartialCharges
from pymol import cmd

class sdfMol:
    """
    Base class for ligand molecules,
    takes SMILES formatted string representatons
    and implements interface to Rdkit and
    Openbabel to generate, convert,
    and display molecules.
    """
    def __init__(self, fPath):
        self.path = fPath
        self.name = 'CID' + fpath.split('_')[-1].split['.'][0]

    def show_mol(self):

        return Draw.ShowMol(Chem.MolFromSmiles(self.__str__()), size=(600, 300))

    def to_image(self):

        return Draw.MolToImage(self.__str__(), size=(600, 300))

    def to_pybel(self):

        return pybel.readstring('smi', self.__str__())

    def write_pdb(self, make_dir=False):

        if make_dir:
            target = os.path.join(
                os.getcwd(), self.type
            )
            os.makedirs(os.path.basename(target), exist_ok=True)
            os.chdir(target)
        else:
            target = os.getcwd()

        x = Chem.SDMolSupplier(self.path)
        x = Chem.AddHs(x)
        AllChem.EmbedMolecule(x, useRandomCoords=True, boxSizeMult=3.0,
                              maxAttempts=10000)
        try:
            AllChem.UFFOptimizeMolecule(x, 10000)
        except ValueError as err:
            print(err, self.name)

        AllChem.MolToPDBFile(x, self.name + '.pdb')
        self.pdb_path = Path(
            os.path.join(target, self.name + '.pdb')
            )

        return self.pdb_path

    def write_pdbqt(self):

        if self.pdb_path:
            pass
        else:
            self.write_pdb()

        self.pdbqt_path = Path(
            self.pdb_path.with_suffix('.pdbqt')
            )
        atoms = list(pybel.readfile('pdb', self.pdb_path.__str__()))
        writer = pybel.Outputfile(
            'pdbqt', self.pdbqt_path.__str__(),
            opt={'pdbqt': '-xh'}, overwrite=True
        )
        for atom in atoms:
            writer.write(atom)
            writer.close()

        '''pybel/pymol cache needs manual cleanup'''
        os.remove(self.pdb_path.__str__())
        self.pdb_path = None
        cmd.reinitialize()

        return self.pdbqt_path

    def construct(self):

        return self.smiles

    def compute_charges(self):

        chem = self.to_chem()
        rdPartialCharges.ComputeGasteigerCharges(chem)
        charges = [(atom.GetSymbol(), atom.GetDoubleProp('_GasteigerCharge'))
                   for atom in chem.GetAtoms()]
        return charges

    def __str__(self):

        return self.construct()

    def __repr__(self):

        return self.construct()