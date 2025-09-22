import re
import os
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
    def __init__(self, mol, name=''):
        self.mol = mol
        self.name = name
        self.type = ''
        self.pdb_path = ''
        self.pdbqt_path = ''

    @classmethod
    def from_sdf(cls, sdf_path):
        path = Path(sdf_path)
        mol = Chem.SDMolSupplier(sdf_path)[0]
        cid = re.search(r'CID_(\d*)', sdf_path)
        if cid:
            name = cid
        else:
            name = path.stem

        return cls(mol, name)

    @classmethod
    def from_smiles(cls, smiles, name=''):
        mol = Chem.MolFromSmiles(smiles)

        return cls(mol, name)

    def show_mol(self):

        return Draw.ShowMol(self.mol, size=(600, 300))

    def to_image(self):

        return Draw.MolToImage(self.mol, size=(600, 300))

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

        h = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(h, useRandomCoords=True, boxSizeMult=3.0,
                              maxAttempts=10000)
        try:
            AllChem.UFFOptimizeMolecule(h, 10000)
        except ValueError as err:
            print(err, self.name)

        AllChem.MolToPDBFile(h, self.name + '.pdb')
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

    def compute_charges(self):

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

    def __init__(self, len, E=[], Z=[], acid=True):

        self.len = len
        self.E, self.Z = E, Z
        self.acid = acid
        if self.acid:
            self.type = 'fatty acid'
        else:
            self.type = 'lipid'
        self.name = str(self.len)
        if self.E:
            self.name += 'E' + ''.join([str(i) for i in self.E])
        if self.Z:
            self.name += 'Z' + ''.join([str(i) for i in self.Z])
        self.mol = Chem.MolFromSmiles(self.construct())

    @classmethod
    def from_series(cls, series, acid=True):
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
    def from_string(cls, string, acid=True):
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

    def __init__(self, r1):

        self.r1 = r1
        self.type = 'LEA'
        self.name = 'LEA' + self.r1.name
        self.mol = Chem.MolFromSmiles(self.construct())

    @classmethod
    def from_series(cls, series):
            r1 = Lipid.from_series(series)
            return cls(r1)

    @classmethod
    def from_string(cls, string):
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

    def __init__(self, pattern, r1='', r2=''):
        self.pattern = pattern
        self.r1 = r1
        self.r2 = r2
        if r2 == '':
            self.type = 'lysoPl'
        else:
            self.type = 'phospholipid'

        self.name = pattern \
                    + self.r1.name \
                    + self.r2.name
        self.mol = Chem.MolFromSmiles(self.construct())

    @classmethod
    def from_series(cls, pattern, series, lipid_table):
        r1 = ''
        r2 = ''
        if series['r1'] != '':
            r1 = Lipid.from_series(
                lipid_table.loc[
                    lipid_table['id'] == series['r1']
                    ].squeeze()
            )

        if series['r2'] != '':
            r2 = Lipid.from_series(
                lipid_table.loc[
                    lipid_table['id'] == series['r2']
                    ].squeeze()
            )

        return cls(pattern, r1, r2)

    @classmethod
    def from_string(cls, pattern, string):
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
        else:
            structure = self.r1.__str__() \
                    + PL_PATTERNS[self.pattern]

        return structure

    def __len__(self):
        if self.r2:
            return self.r1.__len__(), self.r2.__len__()
        else:
            return self.r1.__len__()

