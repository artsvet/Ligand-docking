from __future__ import annotations
import re
import os
from pandas import Series
from pathlib import Path
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdPartialCharges
from pymol import cmd

'''
The most common representation for molecules is the SMILES
molecule format, a string of atom characters representing 
a depth-first tree traversal of the molecular graph.
The ligand classes are used to generate lipid molecule SMILES.
'''

PL_PATTERNS = {
    # template pattern for phospholipid molecule classes
    'pa': 'CC(COP(=O)(O)O)O',
    'pc': 'CC(COP(=O)([O-])OCC[N+](C)(C)C)O',
    'pe': 'CC(COP(=O)(O)OCCN)O',
    'pg': 'CC(COP(=O)(O)OCC(CO)O)O',
    'pi': 'CC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)O',
    'ps': 'CC(COP(=O)(O)OCC(C(=O)O)N)O',
    'cm': '(C(CO)N'
    }

class Ligand:
    """
    Base class for ligand molecules,
    takes SMILES formatted string representatons
    and implements interface to Rdkit and
    Openbabel to generate, convert,
    and display molecules.

    """
    def __init__(self, mol, name, path = '', cid = ''):

        self.path : Path = Path(path) if path else Path(
            os.getcwd() + '\\' + name 
        ) 
        self.pdb_path = self.path if self.path.suffix == '.pdb' else ''
        self.pdbqt_path = self.path if self.path.suffix == '.pdbqt' else ''
        self.sdf_path = self.path if self.path.suffix == '.sdf' else ''
        self.mol : Chem.rdchem.Mol = mol
        self.name : str = name
        self.cid : str = cid
        
    @classmethod
    def from_sdf(cls, sdf_path: str, name: str = '', cid: str = ''):

        path = Path(sdf_path)
        mol = Chem.SDMolSupplier(str(path), removeHs=False)[0]      
        if cid := re.search(r'CID_(\d*)', str(path)): cid = cid[0]
        if not name: name = path.stem

        return cls(mol, name, path=path)
    
    @classmethod
    def from_pdb(cls, pdb_path: str, name: str = '', cid: str = ''):

        path = Path(pdb_path)
        mol = Chem.rdmolfiles.MolFromPDBFile(str(path), removeHs=False)[0]      
        if cid := re.search(r'CID_(\d*)', str(path)): cid = cid[0]
        if not name: name = path.stem

        return cls(mol, name, path=path)

    @classmethod
    def from_smiles(cls, smiles: str, name: str = ''):

        mol = Chem.MolFromSmiles(smiles)
        if name: pass        
        else: name = smiles

        return cls(mol, name)

    def show_mol(self):

        return Draw.ShowMol(self.mol, size=(600, 300))

    def to_image(self):

        return Draw.MolToImage(self.mol, size=(600, 300))

    def to_pybel(self):

        return pybel.readstring('smi', self.__str__())

    def write_sdf(self):

        h = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(h, boxSizeMult=3.0, maxAttempts=10000)
        self.sdf_path = self.path.with_suffix('.sdf')
        Chem.rdmolfiles.SDWriter(str(self.sdf_path)).write(h)
        
        
        return self.sdf_path
    
    def write_pdb(self) -> Path:

        h = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(
            h, boxSizeMult=3.0, maxAttempts=10000
        )

        try:
            AllChem.UFFOptimizeMolecule(h, 10000)
        except ValueError as err:
            print(err, self.name)

        AllChem.MolToPDBFile(h, self.name + '.pdb')
        self.pdb_path = self.path.with_suffix('.pdb')

        return self.pdb_path

    def write_pdbqt(self) -> Path:

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
        cmd.reinitialize()
        
        return self.pdbqt_path
    
    def sdf_to_pdbqt(self) -> Path:

        self.pdb_path =  self.sdf_path.with_suffix('.pdb')
        mol = Chem.SDMolSupplier(str(self.sdf_path), removeHs=False)[0] 
        AllChem.MolToPDBFile(mol, str(self.pdb_path))
        self.pdbqt_path = self.pdb_path.with_suffix('.pdbqt')
        atoms = list(pybel.readfile('pdb', str(self.pdb_path)))
        writer = pybel.Outputfile(
            'pdbqt', str(self.pdbqt_path), overwrite=True
        )
        for atom in atoms:
            writer.write(atom)
            writer.close()
        cmd.reinitialize()
        
        return self.pdbqt_path

    def compute_charges(self) -> list:

        rdPartialCharges.ComputeGasteigerCharges(self.mol)
        charges = [(atom.GetSymbol(), atom.GetDoubleProp('_GasteigerCharge'))
                   for atom in self.mol.GetAtoms()]
        return charges

    def __str__(self):

        return Chem.MolToSmiles(self.mol)

    def __repr__(self):

        return self.name


class Lipid(Ligand):
    """
    Lipid class, generating a string
    representation of carbon chains with variable
    length, double-bonding, and acidity.
    ex - Lipid(8, E=[3,6]) : C/C=C/C/C=C/CC(=O)O
    """

    def __init__(self, len: int, E: list[int] = [], Z: list[int] = [], acid=True):

        self.len = len
        self.E, self.Z = E, Z
        self.acid = acid
        if self.acid:
            self.group = 'fatty acid'
        else:
            self.group = 'lipid'
        self.name = str(self.len)
        if self.E:
            self.name += 'E' + ''.join([str(i) for i in self.E])
        if self.Z:
            self.name += 'Z' + ''.join([str(i) for i in self.Z])
        self.mol = Chem.MolFromSmiles(self.construct())
        super().__init__(self.mol, self.name)


    @classmethod
    def from_series(cls, series: Series, acid=True):
        len = int(series['length'])
        E = [int(E)
             for E in series['E'].split(',')
             if E != ''
            ]
        Z = [int(Z)
             for Z in series['Z'].split(',')
             if Z != ''
            ]
        return cls(len, E=E, Z=Z, acid=acid)

    @classmethod
    def from_string(cls, string: str, acid=True): # string shorthand, i.e. 16E6Z12
        len = int(re.match(r'\d+', string)[0])
        E = []
        Z = []
        if 'E' in string:
            E = [int(E) for E in re.findall(r'\d+', re.split('E', string)[1])]
        if 'Z' in string:
            Z = [int(Z) for Z in re.findall(r'\d+', re.split('Z', string)[1])]

        return cls(len, E=E, Z=Z, acid=acid)

    def construct(self):

        chain = 'C' * self.len

        for i in sorted(self.E, reverse=True):
                chain = (chain[:i-1]
                         + '/C=C/'
                         + chain[i+1:]
                         )

        for i in sorted(self.Z, reverse=True):
                chain = (chain[:i-1]
                         + '/C=C\\'
                         + chain[i+1:]
                         )
        if self.acid:
            return chain[::-1] + '(=O)O'
        else:
            return chain[::-1]

    def __len__(self):
        return self.len


class Lea(Ligand):
    """
    Synthetic Lea drug class,
    generating a string representation using
    LEA template backbone and r1 Lipid.
    ex - Lea(Lipid(8, E=[3,6])) : C/C=C/C/C=C/CCOCC(CNC(C)C)O
    """

    def __init__(self, r1: Lipid):

        self.r1 = r1
        self.group = 'LEA'
        self.name = 'LEA' + self.r1.name
        self.mol = Chem.MolFromSmiles(self.construct())
        super().__init__(self.mol, self.name)

    @classmethod
    def from_series(cls, series: pd.Series):
            r1 = Lipid.from_series(series)
            return cls(r1)

    @classmethod
    def from_string(cls, string: str):
            r1 = Lipid.from_string(string[3:])
            return cls(r1)

    def construct(self):

        return self.r1.__str__().replace('(=O)O','') \
            + 'OCC(CNC(C)C)O'

    def __len__(self):
        return self.r1.__len__()


class Pl(Ligand):
    """
    Phospholipid compound class,
    generating a string representation using
    pattern backbone dict, r1, r2 lipid side chain
    """

    def __init__(self, pattern: str, r1: Lipid, r2: Lipid = ''):
        self.pattern = pattern
        self.r1 = r1
        self.r2 = r2
        if not r2:
            self.group = 'lysoPl'
        else:
            self.group = 'phospholipid'

        self.name = pattern \
                    + self.r1.name \
                    + self.r2.name
        self.mol = Chem.MolFromSmiles(self.construct())
        super().__init__(self.mol, self.name)


    @classmethod
    def from_string(cls, pattern: str, string: str):

        pattern = pattern
        split = re.split('  ', string)

        if len(split) > 1:
            r1 = Lipid.from_string(split[0])
            r2 = Lipid.from_string(split[1])

        else:
            r1 = Lipid.from_string(string)
            r2 = ''

        return cls(pattern, r1, r2)

    def construct(self):

        if self.r2:

            if self.pattern == 'sm':
                structure = self.r1.__str__()[:-1] \
                            + PL_PATTERNS[self.pattern] \
                            + self.r2.__str__()[::-1].replace(
                    'O)O=(C', 'C(=O)'
                )

            else:
                structure = self.r1.__str__() \
                    + PL_PATTERNS[self.pattern] \
                    + self.r2.__str__()[::-1].replace(
                    'O)O=(C', 'C(=O)'
                    )
                
        else: structure = self.r1.__str__() + PL_PATTERNS[self.pattern]

        return structure

    def __len__(self):

        if self.r2: return self.r1.__len__(), self.r2.__len__()
        else: return self.r1.__len__()

