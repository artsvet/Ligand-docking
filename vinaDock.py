import re
import pandas as pd
import subprocess
from datetime import datetime


class Docker:
    """
    Autodock VINA searches for lowest binding energy conformations for
    flexible ligand and rigid protein receptor.
    Outputs compiled vina rankings output as pandas dataframe.
    """
    def __init__(self, receptor, ligand, log_path,
                 box=(0, 0, 0, 30, 30, 30), exhaustiveness=10, cpu=8, run_count=1):
        self.receptor = receptor  # prepared protein pdbqt file
        self.ligand = ligand  # ligand to dock
        self.box = box  # search grid
        self.log = log_path  # vina subprocess output txt
        self.exhaustiveness = exhaustiveness
        self.cpu = cpu
        self.run_count = run_count  # iterations
        self.times_ran = 0  # iterations counter


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

