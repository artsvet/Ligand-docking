import re
import os
import pandas as pd
import pathlib
import threading
import subprocess
from multiprocessing import Process, Queue, current_process
from more_itertools import consume
from pymol import cmd
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


class Multiprocess:
    """
    Performs parallel multiprocess execution
    of a function using a Semaphore
    on arguments consumed to a Queue,
    for cleaning and conversion of Protein files.
    """

    def __init__(self, function, sema_buffer=4, q_size=0):
        self.function = function
        self.sema = threading.BoundedSemaphore(sema_buffer)
        self.q = Queue(maxsize=q_size)

    def wrap(self, arg):
        self.sema.acquire()
        print('In worker process: ',
              "\n", arg,
              "\n", current_process()
        )
        self.function(arg)
        print('done working: ', current_process())
        self.sema.release()

    def __call__(self, *args, **kwargs):
        consume(self.q.put(arg) for arg in args)
        with self.sema:
            while self.q.qsize() > 0:
                a = self.q.get()
                p = Process(target=self.wrap, args=(a,))
                p.start()

class Docker:

    def __init__(self, receptor, ligand, log_path,
                 box=(0,0,0,30,30,30), exhaustiveness=10, run_count=1):
        self.receptor = receptor
        self.ligand = ligand
        self.box = box
        self.log = log_path
        self.exhaustiveness = exhaustiveness
        self.run_count = run_count

    def dock_args(self):

        return 'vina --receptor {0} --ligand {1} ' \
               '--center_x {2} --center_y {3} --center_z {4} ' \
               '--size_x {5} --size_y {6} --size_z {7}' \
               ' --log {8} --exhaustiveness {9} '.format(
                    self.receptor, self.ligand.pdbqt_path,
                    self.box[0], self.box[1], self.box[2],
                    self.box[3], self.box[4], self.box[5],
                    self.log, self.exhaustiveness,
                )

    def dock(self, run_number=1):

        self.ligand.write_pdbqt()
        r = subprocess.run(
            self.dock_args(), shell=True
        )
        list_out = []
        os.remove(self.ligand.pdbqt_path)
        delattr(self.ligand, 'pdbqt_path')
        with open(self.log, 'r') as log:
            for line in log:
                m = re.match(
                    r'(\d)\s*(-\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)',
                    line.lstrip()
                )
                if m:
                    list_out.append({'Receptor': self.receptor.stem,
                                     'Ligand': self.ligand,
                                     'Run number': run_number,
                                     'Rank': m.group(1),
                                     'Affinty': m.group(2),
                                     'Dist rmsd l.b.': m.group(3),
                                     'Dist rmsd u.b.': m.group(4)})
        return pd.DataFrame(list_out)

    def run(self):

        times_ran = 0
        run_outputs = []
        while times_ran < self.run_count:
            run_outputs.append(self.dock())
            times_ran += 1

        return pd.concat(run_outputs)

PL_PATTERNS = {
    # template pattern for phospholipid molecule classes
    'pa':'CC(COP(=O)(O)O)O',
    'pc':'CC(COP(=O)([O-])OCC[N+](C)(C)C)O',
    'pe':'CC(COP(=O)(O)OCCN)O',
    'pg':'CC(COP(=O)(O)OCC(CO)O)O',
    'pi':'CC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)O',
    'ps':'CC(COP(=O)(O)OCC(C(=O)O)N)O',
    'sm':'NC(COP(=O)(O)OCC[N+](C)(C)C)',
    'cm':'(C(CO)N'
    }


class Protein(type(pathlib.Path())):
    """
    Protein class instantiated using path
    to .pdb, .pdbqt, or other structure file.
    Implements methods to clean and convert
    structure files, with parallel
    execution when used with Multiprocess.
    """

    def __init__(self, path):

        super().__init__()

    def clean(self):

        if hasattr(self, 'path_clean'):
            return self.path_clean
        else:
            self.path_clean = Protein(
                self.with_name(
                    self.stem + '.clean' + self.suffix
                )
            )
            setattr(self.path_clean, 'path_clean', self.path_clean)
            with open(self, "r") as f:
                lines = f.readlines()
            with open(self.path_clean, "w") as f:
                for line in lines:
                    if line.strip("\n")[:4] == "ATOM":
                        f.write(line)
                f.close()

        return self.path_clean

    def convert(self):

        if hasattr(self, 'path_pdbqt'):
            return self.path_pdbqt
        else:
            self.path_pdbqt = Protein(self.with_suffix('.pdbqt'))
            setattr(self.path_pdbqt, 'path_pdbqt', self.path_pdbqt)
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
            cmd.reinitialize()

        return self.path_pdbqt

    def prepare(self):

        if hasattr(self, 'path_pdbqt'):
            self.path_pdbqt.clean()
        elif hasattr(self, 'path_clean'):
            self.path_clean.convert().clean()
        else:
            self.clean().convert().clean()

    def __repr__(self):

        return self.root

class Mol:
    """
    Base class for ligand molecules,
    takes SMILES formatted string representatons
    and implements interface to Rdkit and
    Openbabel to generate, convert,
    and display molecules.
    """
    def __init__(self, smiles):
        self.smiles = smiles

    def to_chem(self):

        return Chem.MolFromSmiles(self.__str__())

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

        x = self.to_chem()
        x = Chem.AddHs(x)
        AllChem.EmbedMolecule(x, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(x, 1000)
        AllChem.MolToPDBFile(x, self.name + '.pdb')
        self.pdb_path = pathlib.Path(os.path.join(target, self.name + '.pdb'))

        return self.pdb_path

    def write_pdbqt(self):

        if hasattr(self, 'pdb_path'):
            pass
        else:
            self.write_pdb()

        self.pdbqt_path = self.pdb_path.with_suffix('.pdbqt')
        mols = list(pybel.readfile('pdb', self.pdb_path.__str__()))
        writer = pybel.Outputfile(
            'pdbqt', self.pdbqt_path.__str__(),
            opt={'pdbqt': '-xh'}, overwrite=True
        )
        for molecule in mols:
            writer.write(molecule)
            writer.close()
        os.remove(self.pdb_path.__str__())
        delattr(self, 'pdb_path')
        cmd.reinitialize()

        return self.pdbqt_path

    def construct(self):

        return self.smiles

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
                    + 'E' + ''.join([str(i) for i in self.E]) \
                    + 'Z' + ''.join([str(i) for i in self.Z])

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

    def __len__(self):
        return self.len

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

    def __len__(self):
        return self.r1.__len__()

class Pl(Mol):
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
            self.type = 'phospholipid'
        else:
            self.type = 'lysoPl'
        if self.pattern == 'sm' or self.pattern == 'cm':
            if self.r1 == '':
                self.r1 = Lipid(16, E=[2], acid=False)
            self.r1 = r1
        self.name = pattern \
                    + self.r1.name \
                    + getattr(self.r2, 'name', '')

    def construct(self):

        structure = self.r1.__str__() \
            + PL_PATTERNS[self.pattern] \
            + self.r2.__str__()[::-1].replace(
                'O)O=(C', 'C(=O)'
            )

        if self.pattern == 'sm' or self.pattern == 'cm':
            structure = 'CCCCCCCCCCCCCC=CC(C(COP(=O)([O-])OCC[N+](C)(C)C)N' + ')O'

        return structure

    def __len__(self):
        if self.r2:
            return self.r1.__len__(), self.r2.__len__()
        else:
            return self.r1.__len__()


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
        if series['r2'] != '':
            self.r2 = LipidFromSeries(
                lipidTable.loc[
                    lipidTable['id'] == series['r2']
                    ].squeeze()
                )
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
        self.pdb = Protein(targetsSeries['pdb'])
        self.name = targetsSeries['name']
        self.type = targetsSeries['type']
        self.species = targetsSeries['species']
        self.function = targetsSeries['function']
        self.ligands = re.split(':|,', targetsSeries['lipid'])
        self.ligandsIndexes = [int(i) for i in self.ligands[1:]]
        self.selected = [int(i) for i in re.split(',', targetsSeries['selected']) if i != ',']
        self.box = eval(targetsSeries['box'])

    def build_lea_list(self):
        leaList = []
        for index, series in self.lipidTable.iterrows():
            leaList.append(LeaFromSeries(series))
        return leaList

    def build_ligand_list(self):
        ligandList = []

        if self.ligands[0] == 'all':
            ligandPatterns = [key for key in PL_PATTERNS]
        elif self.ligands[0] == '':
            return ligandList
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
                    ligandList.append(LipidFromSeries(series))
            else:
                for series in plSeries:
                    ligandList.append(PlFromSeries(
                    pattern, series, self.lipidTable
                ))

        return ligandList

    def build_select_list(self):
        ligandList = []
        for index in self.selected:
            ligand = Mol(selectedTable.iloc[index]['smiles'])
            setattr(ligand, 'name', selectedTable.iloc[index]['id'])
            ligandList.append(ligand)
        return ligandList

    def write_ligands(self):
        ligands = self.build_lea_list() + self.build_ligand_list()
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

    def dock(self, ligands, run_count=1, exhaustiveness=10):
        outputs = []
        for mol in ligands:
            outputs.append(
                Docker(receptor=self.pdb,
                ligand=mol,
                log_path=mol.name + '_log.txt',
                box=self.box, run_count=run_count,
                exhaustiveness=exhaustiveness).run())
        return outputs


targets = pd.read_csv(str(os.getcwd() + '\\targetsTable.csv')).fillna('')
lipidTable = pd.read_csv(str(os.getcwd() + '\\ligands\\lipid.csv')).fillna('')
plTable = pd.read_csv(str(os.getcwd() + '\\ligands\\phospholipid.csv')).fillna('')
selectedTable = pd.read_csv(str(os.getcwd() + '\\ligands\\selected.csv')).fillna('')

testTarget = targetsParse(targets.iloc[0], lipidTable, plTable)
print(testTarget.build_lea_list(), '\n', testTarget.build_ligand_list(), '\n', testTarget.build_select_list())
testLigands = testTarget.build_lea_list(), testTarget.build_ligand_list(), testTarget.build_select_list()
consume(print(mol.to_chem()) for mol in testLigands)


for index, row in targets.iloc[2].iterrows():

    if row['box'] != '':
        target = targetsParse(row, lipidTable, plTable)
        to_dock = target.build_ligand_list() + target.build_lea_list() + target.build_select_list()
        testTables = target.dock(to_dock, run_count=3, exhaustiveness=10)
        outputs_path = os.getcwd() + '\\' + row['name'] + '_' + row['id'] + '_outputs.csv'
        pd.concat.testTables.to_csv(outputs_path)

