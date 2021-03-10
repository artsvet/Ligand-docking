import re
import os
import pandas as pd
import pathlib
import threading
import subprocess
import itertools
import time
import numpy as np
from datetime import datetime
from multiprocessing import Process, Queue, current_process
from more_itertools import consume
from pymol import cmd
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdPartialCharges


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
               ' --log {8} --cpu {9} --exhaustiveness {10} '.format(
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

        list_out = []
        with open(log, 'r') as log:
            for line in log:
                m = re.match(
                    r'(\d)\s*(-\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)',
                    line.lstrip()
                )
                if m:
                    list_out.append({'Date_time': datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                                     'Exhaustiveness': self.exhaustiveness,
                                     'Run_number': self.times_ran,
                                     'Receptor': self.receptor.stem,
                                     'Ligand': self.ligand.name,
                                     'Rank': m.group(1),
                                     'Affinity': m.group(2),
                                     'Dist_rmsd_l.b.': m.group(3),
                                     'Dist_rmsd_u.b.': m.group(4)})
        return list_out

    def run(self):

        run_outputs = []
        self.ligand.write_pdbqt()
        while self.times_ran < self.run_count:
            self.times_ran += 1
            run_outputs.extend(self.dock())
        delattr(self.ligand, 'pdbqt_path')

        return pd.DataFrame(run_outputs)


R = .0019872 #kcal mol^-1 K^-1
PL_PATTERNS = {
    # template pattern for phospholipid molecule classes
    'pa':'CC(COP(=O)(O)O)O',
    'pc':'CC(COP(=O)([O-])OCC[N+](C)(C)C)O',
    'pe':'CC(COP(=O)(O)OCCN)O',
    'pg':'CC(COP(=O)(O)OCC(CO)O)O',
    'pi':'CC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)O',
    'ps':'CC(COP(=O)(O)OCC(C(=O)O)N)O',
    #'sm':'NC(COP(=O)(O)OCC[N+](C)(C)C)',
    'cm':'(C(CO)N'
    }
DEFAULT_LEAS = ['LEA8', 'LEA16', 'LEA18Z9', 'LEA20Z581114']
DEFAULT_LIGANDS = {
    # default ligands for phospholipid molecule classes
    'pa':['PA'],
    'pc':['pc1818Z9'],
    'pe':['pe1820Z581114', 'pe1818Z9'],
    'pg':['pg1618Z9'],
    'pi':['pi1820Z581114'],
    'ps':['ps1818Z9'],
    'sm':['sm16', 'sm24']
    }

def KD(dG, T = 310):
    return 1 / np.e ** (-dG / (R * T))


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

    @classmethod
    def from_series(cls, series):
            r1 = Lipid.fromSeries(series)
            return cls(r1)

    @classmethod
    def from_string(cls, string):
            r1 = Lipid.fromString(string[3:])
            return cls(r1)

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

    @classmethod
    def from_series(cls, pattern, series, lipid_table):
        r1 = ''
        r2 = ''
        if series['r1'] != '':
            r1 = Lipid.fromSeries(
                lipid_table.loc[
                    lipid_table['id'] == series['r1']
                    ].squeeze()
            )

        if series['r2'] != '':
            r2 = Lipid.fromSeries(
                lipid_table.loc[
                    lipid_table['id'] == series['r2']
                    ].squeeze()
            )

        return cls(pattern, r1, r2)

    @classmethod
    def from_string(cls, pattern, string):
        pattern = pattern
        if len(split := re.split('  ', string)) > 1:
            r1 = Lipid.fromString(split[0])
            r2 = Lipid.fromString(split[1])
        else:
            r1 = Lipid.fromString(string)
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


class TargetsParse:
    """
    Parses Series to set up dependencies for ligand
    docking. Contains Protein target attributes,
    specifies ligands that will be docked using
    reference lipid and pl tables. Will serve as
    primary interface for running dock and
    collecting outputs.

    todo: abstract dock into class that takes targetsParse
    """

    def __init__(self, dock_series, lipid_table):
        self.lipid_table = lipid_table
        self.id = dock_series['id']
        self.pdb = Protein(dock_series['pdb'])
        self.name = dock_series['name']
        self.type = dock_series['type']
        self.species = dock_series['species']
        self.function = dock_series['function']
        self.lipid_patterns = re.split(', ', dock_series['lipidPatterns'])
        self.lipid_tails = [tails for tails
                                 in re.split('; ', dock_series['lipidtails'])
                                 ]
        self.selected_names = re.split('  ', dock_series['selectedNames'])
        self.selected_smiles = re.split('  ', dock_series['selectedSmiles'])
        assert len(self.selected_names) == len(self.selected_smiles)
        self.refLigands = [ligand for key in self.lipid_patterns
                           if key != ''
                           for ligand in DEFAULT_LIGANDS[key]
                           ]
        self.box = eval(dock_series['box'])
        self.lea = [LeaFromSeries(series)
                    for index, series in self.lipid_table.iterrows()]
        self.lipids = self.build_lipid_list()
        self.selects = self.build_select_list()

    def build_lipid_list(self):
        lipid_list = []

        if self.lipid_patterns == ['']:
            return []

        for pattern, sidechain in itertools.product(
            self.lipid_patterns, self.lipid_tails
            ):

            if pattern not in PL_PATTERNS:
                continue
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

    def annotate_long(self, longDf):

        speciesMask = 1 if self.species == 'H.sapiens' else 0
        leaMask = 1 if longDf.iloc[0]['Ligand'] in DEFAULT_LEAS else 0
        ligMask = 1 if longDf.iloc[0]['Ligand'][:2] in self.lipid_patterns else 0

        longDf['species'] = np.full(len(longDf.index), speciesMask)
        longDf['defaultLea'] = np.full(len(longDf.index), leaMask)
        longDf['defaultLigand'] = np.full(len(longDf.index), ligMask)

        return longDf

    def dock(self, ligands, run_count=1, cpu=8,
        exhaustiveness=10, write_dir=os.getcwd()):

        long = '{0}{1}_{2}_out.csv'.format(
            write_dir, self.name, self.id)
        trimmed = '{0}{1}_{2}_trimmed.csv'.format(
            write_dir, self.name, self.id)
        short = '{0}{1}_{2}_summary.csv'.format(
            write_dir, self.name, self.id)
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        os.chdir(write_dir)

        if os.path.exists(long):
            long_df = pd.read_csv(long)
        else:
            long_df = pd.DataFrame(
                columns=['Date_time', 'Exhaustiveness', 'Run_number',
                        'Receptor', 'Ligand', 'Rank', 'Affinity',
                        'Dist_rmsd_l.b.', 'Dist_rmsd_u.b.', 'species',
                        'defaultLea', 'defaultLigand'
                        ]
            )
            long_df.to_csv(long, index=False, mode='w')
            
        if os.path.exists(trimmed):
            trim_df = pd.read_csv(trimmed)
        else:
            trim_df = long_df
            trim_df.to_csv(trimmed, index=False, mode='w')
            
        if os.path.exists(short):
            short_df = pd.read_csv(short, header=None,
                                   names=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'])
            short_df.to_csv(short, header=False, index=False, mode='w')
        else:    
            short_df = pd.DataFrame(
                columns=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'],
                data=[['', self.id],
                      ['type', self.type],
                      ['species',  self.species],
                      ['lipids', ', '.join(self.lipid_patterns)],
                      [],
                      ['', 'dG mean', 'dG sd', 'KD mean', 'KD sd'],
                      ['ref', np.nan, np.nan, np.nan, np.nan]],
            )

        for mol in ligands:

            cpu_lock = time.localtime(time.time()).tm_hour
            if 21 > cpu_lock > 8:
                cpu = 6

            dock_output = self.annotate_long(
                Docker(receptor=self.pdb,
                       ligand=mol,
                       log_path='{0}\\{1}_{2}_log.txt'.format(
                            write_dir, self.id, mol.name),
                       box=self.box, run_count=run_count,
                       exhaustiveness=exhaustiveness,
                       cpu=cpu
                       ).run()
            )

            long_df = long_df.append(dock_output)
            trim = dock_output.loc[dock_output['Rank'] == '1']
            trim_df = trim_df.append(trim)

            dg = np.mean(trim_df.loc[
                             trim_df['Ligand'] == mol.name, 'Affinity'
                         ].astype(float)), \
                 np.std(trim_df.loc[
                            trim_df['Ligand'] == mol.name, 'Affinity'
                        ].astype(float))
            kd = KD(dg[0]), \
                 KD(dg[0]) * (-dg[1] / dg[0])

            short_df.loc[6, 'dGmean'] = (
                np.mean(trim_df.loc[
                            trim_df['defaultLigand'] == 1, 'Affinity'].astype(float))
                )
            short_df.loc[6, 'dGsd'] = (
                np.std(trim_df.loc[
                           trim_df['defaultLigand'] == 1, 'Affinity'].astype(float))
            )
            short_df.loc[6, 'KDmean'] = KD(short_df.loc[6, 'dGmean'])
            short_df.loc[6, 'KDsd'] = short_df.loc[6, 'KDmean'] * \
                -(short_df.loc[6, 'dGsd'] / short_df.loc[6, 'dGmean'])
            short_df = short_df.append(
                pd.DataFrame(
                    columns=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'],
                    data=[[mol.name, dg[0], dg[1], kd[0], kd[1]]]
                ), ignore_index=True
            )

            dock_output.to_csv(long, index=False, header=False, mode='a')
            trim.to_csv(trimmed, index=False, header=False, mode='a')
            short_df.to_csv(short, index=False, header=False, mode='w')

        return long_df, trim_df, short_df
