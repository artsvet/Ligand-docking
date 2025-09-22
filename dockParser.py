import itertools
import pathlib
import re
import os
import numpy as np
import pandas as pd
from proteins import Multiprocess, Protein
from ligands import PL_PATTERNS, Ligand, Lipid, Lea, Pl
from vinaDock import Docker

'''
Deprecated table IO for VINA docking.
'''

def KD(dG, T = 310):
    """converts enthalpy to coefficient of dissociation"""
    return 1 / np.e ** (-dG / (R * T))


R = .0019872 # kcal mol^-1 K^-1

DEFAULT_LIGANDS = {
    # default endogenous ligands for phospholipases
    'pc': ['pc1818Z9'],
    'pe': ['pe1820Z581114', 'pe1818Z9'],
    'pg': ['pg1618Z9'],
    'pi': ['pi1820Z581114'],
    'ps': ['ps1818Z9'],
    'sm': ['sm16', 'sm24']
    }

DEFAULT_LEAS = ['LEA8', 'LEA16', 'LEA18Z9', 'LEA20Z581114']  # default ligands

class DockerFactory:
    """
    Takes inputs, makes dockers
    """

    def __init__(self, target_id, pdb_path, pdb_name,
                 ligands: list, box: tuple, run_count=1, cpu=16,
                 exhaustiveness=10, write_dir=os.getcwd()):

        self.id = target_id
        if pathlib.Path(pdb_path).suffix == '.pdbqt':
            self.pdbqt = Protein(pdb_path)
        else:
            self.pdbqt = Protein(pdb_path).prepare()
        self.name = pdb_name
        self.ligands = ligands
        self.run_count = run_count
        self.cpu = cpu
        self.exhaustiveness = exhaustiveness
        self.box = box
        self.write_dir = write_dir
        self.didBreak = []
        self.dock_outputs = self.dock()
        self.sparsity = sparsity  # minimum affinity value to continue
        self.run_count = run_count  # iterations
        self.times_ran = 0  # iterations counter
        self.didBreak = 0



    def dock(self):
        write_dir = '{0}\\{1}'.format(
            self.write_dir, self.id)
        long = '{0}\\{1}_out.csv'.format(
            self.write_dir, self.id)
        trimmed = '{0}\\{1}_trimmed.csv'.format(
            self.write_dir, self.id)
        short = '{0}\\{1}_summary.csv'.format(
            self.write_dir, self.id)
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        os.chdir(write_dir)

        if os.path.exists(long):
            long_df = pd.read_csv(long)
        else:
            long_df = pd.DataFrame(
                columns=['Date_time', 'Exhaustiveness', 'Run_number',
                        'Receptor', 'Ligand', 'Rank', 'Affinity',
                        'Dist_rmsd_l.b.', 'Dist_rmsd_u.b.'
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
                      [self.name],
                      ['box', self.box],
                      [],
                      ['', 'dG mean', 'dG sd', 'KD mean', 'KD sd'],
                      ['ref', np.nan, np.nan, np.nan, np.nan]],
            )

        for mol in ligands:

            '''
            cpu_lock = time.localtime(time.time()).tm_hour
            if 21 > cpu_lock > 8:
                cpu = 6
            '''
            try:
                dock = Docker(receptor=self.pdb,
                    ligand=mol,
                    log_path='{0}\\{1}_{2}_log.txt'.format(
                    write_dir, self.id, mol.name),
                    box=self.box, run_count=run_count,
                    exhaustiveness=exhaustiveness,
                    cpu=cpu
                    ).run()
            except:
                logging.exception('exception!!!')
                continue


            dock_output, broke = dock[0], dock[1]
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

            short_df.loc[5, 'KDmean'] = KD(short_df.loc[5, 'dGmean'])
            short_df.loc[5, 'KDsd'] = short_df.loc[5, 'KDmean'] * \
                -(short_df.loc[5, 'dGsd'] / short_df.loc[5, 'dGmean'])
            short_df = short_df.append(
                pd.DataFrame(
                    columns=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'],
                    data=[[mol.name, dg[0], dg[1], kd[0], kd[1]]]
                ), ignore_index=True
            )

            dock_output.to_csv(long, index=False, header=False, mode='a')
            trim.to_csv(trimmed, index=False, header=False, mode='a')
            short_df.to_csv(short, index=False, header=False, mode='w')
            if broke:
                self.didBreak.append(mol.name)

        return long_df, trim_df, short_df, self.didBreak




class TargetsParse:
    """
    Parses Series to set up dependencies for batched ligand
    docking: combining and collecting various properties and formats
    of ligands as lists or reference tables, generating lists of ligands,
    retrieving Protein target attributes. Serves as primary
    interface for running dock and collecting outputs.

    todo: abstract dock into separate class for better flow control and interoperability
    """

    def __init__(self, dock_series):
        """Collects attributes for batched docking of ligands to a single Protein"""

        self.id = dock_series['id']
        self.pdb = Protein(dock_series['pdb'])
        self.name = dock_series['name']
        self.selected_names = re.split(' ', dock_series['selectedNames'])
        self.selected_smiles = re.split(' ', dock_series['selectedSmiles'])
        assert len(self.selected_names) == len(self.selected_smiles)
        self.box = eval(dock_series['box'])
        self.lea = [Lea.from_series(series)
                    for index, series in self.lipid_table.iterrows()]
        self.selects = self.build_select_list()

    def build_select_list(self):
        """Generates list of additional ligands as SMILES"""

        select_list = []
        if self.selected_names == ['']:
            return []

        for i, name in enumerate(self.selected_names):
            ligand = Ligand(self.selected_smiles[i])
            ligand.name = name
            select_list.append(ligand)
        return select_list

    def annotate_long(self, longDf):
        """annotates the long Dataframe with some extra variables for analysis"""

        speciesMask = 1 if self.species == 'H.sapiens' else 0
        leaMask = 1 if longDf.iloc[0]['Ligand'] in DEFAULT_LEAS else 0
        ligMask = 1 if longDf.iloc[0]['Ligand'][:2] in self.lipid_patterns else 0

        longDf['species'] = np.full(len(longDf.index), speciesMask)
        longDf['defaultLea'] = np.full(len(longDf.index), leaMask)
        longDf['defaultLigand'] = np.full(len(longDf.index), ligMask)

        return longDf

    def dock(self, ligands, run_count=1, cpu=8,
             exhaustiveness=10, write_dir=os.getcwd()):
        """Runs the dock, collects and writes output tables"""

        long = '{0}{1}_{2}_out.csv'.format(
            write_dir, self.name, self.id)
        trimmed = '{0}{1}_{2}_trimmed.csv'.format(
            write_dir, self.name, self.id)
        short = '{0}{1}_{2}_summary.csv'.format(
            write_dir, self.name, self.id)
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        os.chdir(write_dir)

        '''Long table format, all outputs'''
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

        '''trimmed table is highest ranking binding from each dock iteration only'''
        if os.path.exists(trimmed):
            trim_df = pd.read_csv(trimmed)
        else:
            trim_df = long_df
            trim_df.to_csv(trimmed, index=False, mode='w')

        '''short table format with summary statistics of all ligand docks'''
        if os.path.exists(short):
            short_df = pd.read_csv(short, header=None,
                                   names=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'])
            short_df.to_csv(short, header=False, index=False, mode='w')
        else:
            short_df = pd.DataFrame(
                columns=['', 'dGmean', 'dGsd', 'KDmean', 'KDsd'],
                data=[['', self.id],
                      ['type', self.type],
                      ['species', self.species],
                      ['lipids', ', '.join(self.lipid_patterns)],
                      [],
                      ['', 'dG mean', 'dG sd', 'KD mean', 'KD sd'],
                      ['ref', np.nan, np.nan, np.nan, np.nan]],
            )

        for mol in ligands:

            '''if you need cores for other things while docks are running'''
            # cpu_lock = time.localtime(time.time()).tm_hour
            # if 21 > cpu_lock > 8:
            #     cpu = 6

            '''vinaDock.py init and output long table'''
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

            '''collecting output and doing math to write to short table'''
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