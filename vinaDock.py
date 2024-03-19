import re
import pandas as pd
import subprocess
from pathlib import Path
from targets import Protein
from ligands import Ligand
from datetime import datetime
from time import sleep


class Docker:
    """
    Autodock VINA searches for lowest binding energy conformations for
    flexible ligand and rigid protein receptor.
    Outputs compiled vina rankings output as pandas dataframe.
    """
    def __init__(self, receptor: Protein, ligand: Ligand, log: Path,
                 box: tuple, exhaustiveness: int, cpu: int):
        
        self.receptor = receptor  # prepared protein pdbqt file
        self.ligand = ligand  # ligand to dock
        self.log = log  # vina subprocess output log.txt
        self.box = box  # search grid
        self.exhaustiveness = exhaustiveness  # vina search breadth
        self.cpu = cpu  # number of threads
        self.out = pd.DataFrame()

    def dock_args(self):

        return 'vina --receptor {0} --ligand {1} ' \
               '--center_x {2} --center_y {3} --center_z {4} ' \
               '--size_x {5} --size_y {6} --size_z {7}' \
               ' --log {8} --cpu {9} --exhaustiveness {10}'.format(
                    self.receptor.path, self.ligand.pdbqt_path,
                    self.box[0], self.box[1], self.box[2],
                    self.box[3], self.box[4], self.box[5],
                    self.log, self.cpu, self.exhaustiveness
                )

    def dock(self):

        self.ligand.write_pdbqt()
        subprocess.run(
            self.dock_args(), shell=True
        )
        self.out = self.scrape_log()

        return self.out

    def scrape_log(self) -> pd.DataFrame:

        list_out = []
        with open(self.log, 'r') as log:
            for line in log:
                m = re.match(
                    r'(\d)\s*(-?\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)',
                    line.lstrip()
                )
                if m:
                    list_out.append({'Date_time': datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                                     'Exhaustiveness': self.exhaustiveness,
                                     'Receptor': self.receptor.name,
                                     'Ligand': str(self.ligand.name),
                                     'Rank': m.group(1),
                                     'Affinity': m.group(2),
                                     'Dist_rmsd_l.b.': m.group(3),
                                     'Dist_rmsd_u.b.': m.group(4)})
        return pd.DataFrame(list_out)