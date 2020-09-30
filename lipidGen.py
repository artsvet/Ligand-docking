import os
import re
import pandas as pd
import pathlib
import threading
from pymol import cmd
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


class SafeThreads:
    """
    Performs parallel thread-safe execution
    of a function using a Semaphore
    on a list of arguments consumed to a Queue,
    for cleaning and conversion of Protein files.
    """

    def __init__(self, function, sema_buffer=4, q_size=0):
        self.function = function
        self.sema = threading.BoundedSemaphore(sema_buffer)
        self.q = Queue(maxsize=q_size)
        self.threads = []

    def wrap(self, arg):
        self.sema.acquire()
        print('In worker thread: ',
              "\n", arg,
              "\n", threading.current_thread()
        )
        self.function(arg)
        self.sema.release()
        print('done working: ', threading.current_thread())

    def __call__(self, *args, **kwargs):
        consume(self.q.put(arg) for arg in args)
        with self.sema:
            while self.q.qsize() > 0:
                t = threading.Thread(target=self.wrap, args=(self.q.get(),))
                t.start()
                self.threads = threading.enumerate()


PL_PATTERNS = {
    # template pattern for phospholipid molecule classes
    'pa':'CC(COP(=O)(O)O)O',
    'pc':'CC(COP(=O)([O-])OCC[N+](C)(C)C)O',
    'pe':'CC(COP(=O)(O)OCCN)O',
    'pg':'CC(COP(=O)(O)OCC(CO)O)O',
    'pi':'CC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)O',
    'ps':'CC(COP(=O)(O)OCC(C(=O)O)N)O',
    'sm':'(C(COP(=O)([O-])OCC[N+](C)(C)C)N',
    'cm':'(C(CO)N'
    }


class Protein(type(pathlib.Path())):
    """
    Protein class instantiated using path
    to .pdb, .pdbqt, or other structure file.
    Implements methods to clean and convert
    structure files, with thread-safe parallel
    execution if used with SafeThreads.
    """

    def __init__(self, path):

        super().__init__()
        
    def clean(self):

        if '.clean' in self.suffixes:
            return self
        else:
            self.path_clean = Protein(
                self.with_name(
                    self.stem + '.clean' + self.suffix
                )
            )
            with open(self, "r") as f:
                lines = f.readlines()
            with open(self.path_clean, "w") as f:
                for line in lines:
                    if line.strip("\n")[:4] == "ATOM":
                        f.write(line)
                f.close()

    def convert(self):

        if '.pdbqt' in self.suffixes:
            return self
        else:
            self.path_pdbqt = Protein(self.with_suffix('.pdbqt'))
            cmd.load(self.__str__())
            cmd.remove('resn HOH')
            cmd.h_add(selection='acceptors or donors')
            cmd.save(self.__str__())
            mols = list(pybel.readfile('pdb', self.path_clean.__str__()))
            writer = pybel.Outputfile(
                'pdbqt', self.path_pdbqt.__str__(), opt={'pdbqt': '-xh'}
            )
            for molecule in mols:
                writer.write(molecule)
                writer.close()

    def prepare(self):

        if hasattr(self, 'path_pdbqt'):
            self.path_pdbqt.clean()
        elif hasattr(self, 'path_clean'):
            self.path_clean.convert().clean()
        else:
            self.clean().convert().clean()

class Mol:
    """
    Base class for ligand molecules,
    takes SMILES formatted string representatons
    and implements interface to Rdkit and
    Openbabel to generate, convert,
    and display molecules.
    """

    def to_chem(self):

        return Chem.MolFromSmiles(self.__str__())

    def show_mol(self):

        return Draw.ShowMol(Chem.MolFromSmiles(self.__str__()), size=(600, 300))

    def to_image(self):

        return Draw.MolToImage(self.__str__(), size=(600, 300))

    def to_pybel(self):

        return pybel.readstring('smi', self.__str__())

    def write_pdb(self):

            target = os.path.join(
                os.getcwd(), self.type
            )
            os.makedirs(os.path.basename(target), exist_ok=True)
            os.chdir(target)
            x = self.to_chem()
            x = Chem.AddHs(x)
            AllChem.EmbedMolecule(x, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(x, 1000)
            AllChem.MolToPDBFile(x, self.name + '.pdb')

    def __str__(self):

        return self.construct()

    def __repr__(self):

        return self.construct()

class Lipid(Mol):
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
        self.name = str(self.len) \
                    + 'E' + str(self.E) \
                    + 'Z' + str(self.Z)

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
        if self.acid == True:
            return chain[::-1] + '(=O)O'
        else:
            return chain[::-1]

class Lea(Mol):
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

    def construct(self):

        return self.r1.__str__().replace('(=O)O','') \
            + 'OCC(CNC(C)C)O'

class Pl(Mol):
    """
    Phospholipid compound class,
    generating a string representation using
    pattern backbone dict, r1, r2 lipid side chain
    """

    def __init__(self, pattern, r1='', r2=''):
        self.pattern = pattern
        if r2:
            self.type = 'phospholipid'
        else:
            self.type = 'lysoPl'
        if self.pattern == 'sm' or self.pattern == 'cm':
            self.r1 = Lipid(16, E=[2], acid=False)
        else:
            self.r1 = r1
        self.r2 = r2
        self.name = pattern \
                    + self.r1.name \
                    + getattr(self.r2, 'name', default='')

    def construct(self):

        structure = self.r1.__str__() \
            + PL_PATTERNS[self.pattern] \
            + self.r2.__str__()[::-1].replace(
                'O)O=(C', 'C(=O)'
            )

        if pattern == 'sm' or self.pattern == 'cm':
            structure = structure + ')O'

        return structure

class LipidFromSeries(Lipid):
    """
    Parses series input to
    generate Lipid objects
    """

    def __init__(self, series, acid=True):
        self.len = series['length']
        self.E = [int(E)
                  for E in series['E'].split(',')
                  if E != ''
                  ]
        self.Z = [int(Z)
                  for Z in series['Z'].split(',')
                  if Z != ''
                  ]
        super().__init__(self.len, E=self.E, Z=self.Z, acid=acid)

class LeaFromSeries(Lea):
    """
    Parses series input to
    generate LEA objects
    """

    def __init__(self, series):
        self.name = series['id']
        self.r1 = LipidFromSeries(series)
        super().__init__(self.r1)

class PlFromSeries(Pl):
    """
    Parses series input to
    generate phospholipid objects
    """
    def __init__(self, pattern, series, lipidTable):
        self.pattern = pattern
        self.name = pattern + ' ' + series['id']
        self.r1 = LipidFromSeries(
            lipidTable.loc[
                lipidTable['id'] == series['r1']
                ].squeeze()
            )
        if series['r2']:
            self.r2 = LipidFromSeries(
                lipidTable.loc[
                    lipidTable['id'] == series['r2']
                    ].squeeze()
                ).__str__()[::-1].replace('O)O=(C', 'C(=O)')
        else:
            self.r2 = ''
        super().__init__(self.pattern, self.r1, self.r2)

class targetsParse:
    """
    Parses Series to set up dependencies for ligand
    docking. Contains Protein target attributes,
    specifies ligands that will be docked using
    reference lipid and pl tables. Will serve as
    primary interface for running dock and
    collecting outputs.
    """

    def __init__(self, targetsSeries, lipidTable, plTable):
        self.lipidTable = lipidTable
        self.plTable = plTable
        self.id = targetsSeries['id']
        self.name = targetsSeries['name']
        self.type = targetsSeries['type']
        self.species = targetsSeries['species']
        self.function = targetsSeries['function']
        self.ligands = re.split(':|,', targetsSeries['lipid'])
        self.ligandsIndexes = self.ligands[1:]
        self.box = targetsSeries['box']

    def build_lea_dict(self):
        leaDict = {}
        for index, series in self.lipidTable.iterrows():
            leaDict[series['id']] = LeaFromSeries(series)
        return leaDict

    def build_ligand_dict(self):
        ligandDict = {}

        if self.ligands[0] == 'all':
            ligandPatterns = [key for key in PL_PATTERNS]
        else:
            ligandPatterns = [self.ligands[0]]

        if self.ligandsIndexes:
            if ligandPatterns == 'lipid':
                lipidSeries = [self.lipidTable.iloc[i] for i in self.ligandsIndexes]
            else:
                plSeries = [self.plTable.iloc[i] for i in self.ligandsIndexes]

        else:
            plSeries = [series for index, series in self.plTable.iterrows()]

        for pattern in ligandPatterns:
            if pattern == 'lipid':
                for series in lipidSeries:
                    ligandDict[self.lipidTable['id']] = LipidFromSeries(series)
            else:
                for series in plSeries:
                    ligandDict[pattern + ' ' + series['id']] = PlFromSeries(
                    pattern, series, self.lipidTable
                )

        return ligandDict

    def write_ligands(self):
        ligands = {**self.build_lea_dict(), **self.build_ligand_dict()}
        root = os.getcwd()
        dir = os.path.join(
                root, self.name, self.id
            )
        os.makedirs(dir)
        os.chdir(dir)
        for ligand in ligands:
            ligands[ligand].write_pdb()
            os.chdir(dir)
        os.chdir(root)


targets = pd.read_csv(str(os.getcwd() + '\\targetsTable.csv')).fillna('')
lipidTable = pd.read_csv(str(os.getcwd() + '\\ligands\\lipid.csv')).fillna('')
plTable = pd.read_csv(str(os.getcwd() + '\\ligands\\phospholipid.csv')).fillna('')