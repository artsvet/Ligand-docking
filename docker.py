import re
import pandas as pd
import subprocess
import sys
import os
import json
from pathlib import Path
from targets import Protein
from ligands import Ligand
from datetime import datetime
from time import sleep
from typing import Dict, Tuple

UNIMOL_WEIGHTS = r""
UNIMOL_INTERFACE = r""


class Docker:
  
    def __init__(self, receptor, ligand, box = (0,0,0,20,20,20), log_path = ''): 

        self.receptor : Protein = receptor  # prepared protein pdbqt file
        self.ligand : Ligand = ligand  # ligand to dock
        self.box : Tuple = box  # search grid
        self.log_path : Path = log_path if log_path else \
            Path(f'{os.getcwd()}\\{self.receptor.name}\\{self.ligand.name}_log.txt')
        self.ligand.path = self.log_path.with_name(self.ligand.name)


    def grid_to_json(self):

        d = {"center_x": self.box[0], "center_y": self.box[1], "center_z": self.box[2],
             "size_x": self.box[3], "size_y": self.box[4], "size_z": self.box[5]
        }        
        with open('docking_grid.json', 'w') as f:
            json.dump(d, f)

        return os.getcwd() + '\\docking_grid.json'
    
    def vina_dock_args(self, cpu = 12, exhaustiveness = 10):

        return 'vina --receptor {0} --ligand {1} ' \
               '--center_x {2} --center_y {3} --center_z {4} ' \
               '--size_x {5} --size_y {6} --size_z {7}' \
               ' --log {8} --cpu {9} --exhaustiveness {10}'.format(
                    self.receptor.pdbqt_path, self.ligand.pdbqt_path,
                    self.box[0], self.box[1], self.box[2],
                    self.box[3], self.box[4], self.box[5],
                    self.log_path, cpu, exhaustiveness
                )
    
    def unimol_dock_args(self):

        return 'python {0} --mode single --conf-size 10 --cluster ' \
               '--input-protein {1} --input-ligand {2} --output-ligand-dir {3} ' \
               '--output-ligand-name {4} --model-dir {5} --input-docking-grid {6} ' \
               '--steric-clash-fix'.format(
                    UNIMOL_INTERFACE, self.receptor.pdb_path, self.ligand.sdf_path, 
                    self.ligand.sdf_path.parent, self.ligand.dock_out.stem, 
                    UNIMOL_WEIGHTS, self.grid_to_json()
               )

    def vina_score_args(self):
        
        return 'vina --receptor {0} --ligand {1} --log {2} --score_only'.format(
                    self.receptor.pdbqt_path, self.ligand.pdbqt_path, self.log_path
                )

    def unimol_dock(self):

        os.makedirs(self.log_path.parent, exist_ok=True)
        os.chdir(self.log_path.parent)
        sys.stdout.write('{0}  unimol_dock {1} on {2}: Box - {3} \n'.format(
            datetime.now(), self.ligand.name, self.receptor.name, self.box
            )
        )
        self.ligand.write_sdf()
        self.ligand.dock_out = self.ligand.sdf_path.with_name(
            self.ligand.name + '_unimol_dock_out.sdf'
        )
        subprocess.run(self.unimol_dock_args(), shell=True)
        self.ligand.sdf_path = self.ligand.dock_out
        self.ligand.pdbqt_path = self.ligand.sdf_to_pdbqt()
        os.remove(self.ligand.dock_out.with_suffix('.lmdb'))
        
        return self.vina_score()
    
    def vina_score(self):

        os.makedirs(self.log_path.parent, exist_ok=True)
        os.chdir(self.log_path.parent)
        sys.stdout.write('{0}  Scoring {1} on {2}: '.format(
            datetime.now(), self.ligand.name, self.receptor.name,
            )
        )
        subprocess.run(self.vina_score_args(), shell=True)
        self.score_out = self.scrape_score()
        sys.stdout.write(
            'Finished scoring {0} on {1}: \nScore -\n{2}\n\n'.format(
                self.ligand.name, self.receptor.name, str(self.score_out['Affinity'])
            )
        )

        return self.score_out
    
    def vina_dock(self, exhaustiveness = 10, cpu = 0, run_count = 3, 
        early_stop = False, stop_val: float = -7.0, clean_up = True) -> pd.DataFrame:
        
        assert ' ' not in str(self.receptor.path) + str(self.log_path), \
            'VINA does not accept whitespaces in docking args'
        os.makedirs(self.log_path.parent, exist_ok=True)
        os.chdir(self.log_path.parent)
        sys.stdout.write('{0}  Docking {1} on {2}: \nCPU - {3}, ' \
               'Exhaustiveness - {4}, Box - {5}, Run count - {6} \n'.format(
            datetime.now(), self.ligand.name, self.receptor.name,
            cpu, exhaustiveness,  self.box, run_count))
        
        if not hasattr(self, 'run_count'):
            self.run_count = 0
            self.dock_out = pd.DataFrame()
            self.did_break = False

        for run in range(run_count):
            self.run_count += 1
            self.ligand.write_pdbqt()
            subprocess.run(self.vina_dock_args(cpu, exhaustiveness), shell=True)
            self.ligand.dock_out = self.ligand.pdbqt_path.with_name(
                self.ligand.name + '_vina_out.pdbqt'
            )
            os.rename(
                f'{os.getcwd()}\\{self.receptor.name}\\{self.ligand.name}_out.pdbqt',
                str(self.ligand.dock_out)
            )
            log = self.scrape_log()
            log = log.assign(exhaustiveness = exhaustiveness)
            self.dock_out = pd.concat([self.dock_out, log])
            
            if early_stop and float(log['Affinity'][0]) > stop_val:
                self.did_break = True
                sys.stdout.write(
                    '{0} Docking {1} on {2} broke!!!\nTop affinity -\n{3}\n\n'.format(
                        datetime.now(), self.ligand.name, self.receptor.name, 
                        str(self.dock_out.loc[self.dock_out['Rank']==1, 'Affinity'])
                    )
                )
                return self.dock_out

        sys.stdout.write(
            '{0} Finished docking {1} on {2}: \nTop affinity -\n{3}\n\n'.format(
                datetime.now(), self.ligand.name, self.receptor.name, 
                str(self.dock_out.loc[self.dock_out['Rank']==1, 'Affinity'])
            )
        )
        
        if clean_up:
            os.remove(self.ligand.pdb_path)
            os.remove(self.ligand.pdbqt_path)

        self.ligand.pdbqt_path = self.ligand.dock_out

        return self.dock_out

    def scrape_score(self):

        with open(self.log_path, 'r') as log:
            log = ''.join([l for l in log])
            m = re.findall(r'\s*\w*\s*:\s*(-?\d+.\d*)', log)

        if m:
            score = pd.Series(
                data=[datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                      str(self.receptor.name), str(self.ligand.name)] + 
                      [float(match) for match in m],
                index=['Date_time', 'Receptor', 'Ligand',
                       'Affinity', 'Gauss 1', 'Gauss 2',
                       'Repulsion', 'Hydrophobic', 'Hydrogen']
            )
        else: return
                             
        return score

    def scrape_log(self) -> pd.DataFrame:
        
        df = pd.DataFrame()

        with open(self.log_path, 'r') as l:

            for line in l:
                m = re.match(
                    r'(\d)\s*(-?\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)',
                    line.lstrip()
                )

                if m:
                    df = df.append(pd.DataFrame(
                        data=[[datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                              [],
                              str(self.receptor.name),
                              str(self.ligand.name),
                              int(m.group(1)),
                              float(m.group(2)),
                              float(m.group(3)),
                              float(m.group(4))]],
                        columns=['Date_time',
                                 'Exhaustiveness',
                                 'Receptor',
                                 'Ligand',
                                 'Rank',
                                 'Affinity',
                                 'Dist_rmsd_l.b.',
                                 'Dist_rmsd_u.b.']), 
                                 ignore_index=True)
                    
        return df