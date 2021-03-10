import re
import os
import pandas as pd
import pathlib
import threading
import subprocess
import itertools
from datetime import datetime
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
                 box=(0,0,0,30,30,30), exhaustiveness=10, cpu=8, run_count=1):
        self.receptor = receptor
        self.ligand = ligand
        self.box = box
        self.log = log_path
        self.exhaustiveness = exhaustiveness
        self.cpu = cpu
        self.run_count = run_count
        self.times_ran = 0


    def dock_args(self):

        return 'vina --receptor {0} --ligand {1} ' \
               '--center_x {2} --center_y {3} --center_z {4} ' \
               '--size_x {5} --size_y {6} --size_z {7}' \
               ' --log {8} --cpu {8} --exhaustiveness {9} '.format(
                    self.receptor, self.ligand.pdbqt_path,
                    self.box[0], self.box[1], self.box[2],
                    self.box[3], self.box[4], self.box[5],
                    self.log, self.cpu, self.exhaustiveness
                )

    def dock(self):

        subprocess.run(
            self.dock_args(), shell=True
        )
        return self.scrape_log(self.log)

    def scrape_log(self, log):

        with open(log, 'r') as log:
            for line in log:
                m = re.match(
                    r'(\d)\s*(-\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)',
                    line.lstrip()
                )
                if m:
                    dict_out = {'Date_time': datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                                'Exhaustiveness': self.exhaustiveness,
                                'Run_number': self.times_ran,
                                'Receptor': self.receptor.stem,
                                'Ligand': self.ligand.name,
                                'Rank': m.group(1),
                                'Affinity': m.group(2),
                                'Dist_rmsd_l.b.': m.group(3),
                                'Dist_rmsd_u.b.': m.group(4)}
        return pd.DataFrame(dict_out)

    def run(self):

        run_outputs = []
        self.ligand.write_pdbqt()
        while self.times_ran < self.run_count:
            run_outputs.append(self.dock())
            self.times_ran += 1
        os.remove(self.ligand.pdbqt_path)
        delattr(self.ligand, 'pdbqt_path')
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
        os.remove(self)

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
        os.remove(self)

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
    def __init__(self, smiles, name=''):
        self.smiles = smiles
        self.name = name

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
        AllChem.EmbedMolecule(x, useRandomCoords=True, boxSizeMult=3.0,
                              maxAttempts=10000)
        try:
            AllChem.UFFOptimizeMolecule(x, 10000)
        except ValueError as err:
            print(err, self.name)


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
        self.name = str(self.len)
        if self.E:
            self.name += 'E' + ''.join([str(i) for i in self.E])
        if self.Z:
            self.name += 'Z' + ''.join([str(i) for i in self.Z])


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
            self.type = 'lysoPl'
        else:
            self.type = 'phospholipid'

        self.name = pattern \
                    + getattr(self.r1, 'name', '') \
                    + getattr(self.r2, 'name', '')

    def construct(self):
        if self.r2 != '':
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

class LipidFromSeries(Lipid):
    """
    Parses series input to
    generate Lipid objects
    """

    def __init__(self, series, acid=True):
        self.len = int(series['length'])
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
    def __init__(self, pattern, series, lipid_table):
        self.pattern = pattern
        self.name = pattern + ' ' + series['id']
        self.r1 = ''
        self.r2 = ''
        if series['r1'] != '':
            self.r1 = LipidFromSeries(
                lipid_table.loc[
                    lipid_table['id'] == series['r1']
                    ].squeeze()
                )

        if series['r2'] != '':
            self.r2 = LipidFromSeries(
                lipid_table.loc[
                    lipid_table['id'] == series['r2']
                    ].squeeze()
                )

        super().__init__(self.pattern, self.r1, self.r2)

class LipidFromString(Lipid):

    def __init__(self, string, acid=True):
        self.len = int(re.match('\d+', string)[0])
        self.E = []
        self.Z = []
        if 'E' in string:
            self.E = [int(E) for E in re.findall('\d+', re.split('E', string)[1])]
        if 'Z' in string:
            self.Z = [int(Z) for Z in re.findall('\d+', re.split('Z', string)[1])]

        super().__init__(self.len, E=self.E, Z=self.Z, acid=acid)

class PlFromString(Pl):

    def __init__(self, pattern, string):
        self.pattern = pattern
        if len(split := re.split('  ', string)) > 1:
            self.r1 = LipidFromString(split[0])
            self.r2 = LipidFromString(split[1])
        else:
            self.r1 = LipidFromString(string)
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

    def __init__(self, dock_series, lipid_table):
        self.lipid_table = lipid_table
        self.id = dock_series['id']
        self.pdb = Protein(dock_series['pdb'])
        self.name = dock_series['name']
        self.type = dock_series['type']
        self.species = dock_series['species']
        self.function = dock_series['function']
        if 'all' in dock_series['lipidPatterns']:
            self.lipid_patterns = PL_PATTERNS.keys()
        else:
            self.lipid_patterns = re.split(',', dock_series['lipidPatterns'])
        self.lipid_sidechains = [sidechains for sidechains
                                 in re.split('; ', dock_series['lipidSidechains'])
                                 ]
        self.selected_names = re.split('  ', dock_series['selectedNames'])
        self.selected_smiles = re.split('  ', dock_series['selectedSmiles'])
        assert len(self.selected_names) == len(self.selected_smiles)
        self.box = eval(dock_series['box'])

    def build_lea_list(self):
        lea_list = []
        for index, series in self.lipid_table.iterrows():
            lea_list.append(LeaFromSeries(series))
        return lea_list

    def build_lipid_list(self):
        lipid_list = []

        if self.lipid_patterns == ['']:
            return []

        for pattern, sidechain in itertools.product(
            self.lipid_patterns, self.lipid_sidechains
            ):

            if pattern == 'lipid':
                try:
                    assert ' ' not in sidechain
                except AssertionError:
                    continue
                lipid_list.append(LipidFromString(sidechain))
            else:
                lipid_list.append(PlFromString(pattern, sidechain))

        return lipid_list


    def build_select_list(self):
        select_list = []
        if self.selected_names == ['']:
            return []

        for i, name in enumerate(self.selected_names):
            ligand = Mol(self.selected_smiles[i])
            setattr(ligand, 'name', name)
            select_list.append(ligand)
        return select_list

    def write_ligands(self):
        ligands = self.build_lea_list() + self.build_ligand_list() + self.build_select_list()
        root = os.getcwd()
        dir = os.path.join(
            root, self.name, self.id
        )
        os.makedirs(dir)
        os.chdir(dir)
        for mol in ligands:
            mol.write_pdb()
            os.chdir(dir)
        os.chdir(root)

    def dock(self, ligands, run_count=1, exhaustiveness=10, write_dir=os.getcwd()):
        outputs = []
        start_time = datetime.datetime.now()
        for mol in ligands:
            outputs.append(
                Docker(receptor=self.pdb,
                       ligand=mol,
                       log_path=self.id + '_' + mol.name + '_log.txt',
                       box=self.box, run_count=run_count,
                       exhaustiveness=exhaustiveness).run())
            pd.concat(outputs).to_csv('{0}\\{1}_{2}_{3}.csv'.format(
                write_dir, self.name, self.id, str(start_time[:-10])))
        return outputs

targets = pd.read_csv(str(os.getcwd() + '\\targetsTable.csv')).fillna('')
lipidTable = pd.read_csv(str(os.getcwd() + '\\ligands\\lipid.csv')).fillna('')
