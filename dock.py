import subprocess
import re
import pandas as pd

R = 'C:\\Users\\starch\\PycharmProjects\\molDocking\\testDock1019\\5lia.pdbqt'
L = 'C:\\Users\\starch\\PycharmProjects\\molDocking\\testDock1019\\LPA.pdbqt'
W = 'C:\\Users\\starch\\PycharmProjects\\molDocking\\testDock1019\\LPA_log.txt'
B = (-42, 24, -14, 30, 30, 30)


class docker:

    def __init__(self, receptor, ligand, log,
                 box=(0,0,0,30,30,30), exhaustiveness=10, run_count=1):
        self.receptor = receptor
        self.ligand = ligand
        self.box = box
        self.log = log
        self.exhaustiveness = exhaustiveness
        self.run_count = run_count

    def dock_args(self):

        return 'vina --receptor {0} --ligand {1} ' \
               '--center_x {2} --center_y {3} --center_z {4} ' \
               '--size_x {5} --size_y {6} --size_z {7}' \
               ' --exhaustiveness {8} --log {9}'.format(
                    self.receptor.__str__(), self.ligand.__str__(),
                    self.box[0], self.box[1], self.box[2],
                    self.box[3], self.box[4], self.box[5],
                    self.exhaustiveness, self.log.__str__()
                )

    def dock(self, run_number):

        r = subprocess.run(
            self.dock_args(), shell=True
        )
        list_out = []
        with open(W, 'r') as log:
            for line in log:
                m = re.match(
                    r'(\d)\s*(-\d*.\d*)\s*(\d*.\d*)\s*(\d*.\d*)',
                    line.lstrip()
                )
                if m:
                    list_out.append({'Receptor': self.receptor.name,
                                     'Ligand': self.ligand.name,
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
            times_ran += 1
            run_outputs.append(self.dock(run_number=times_ran))


        return pd.concat(run_outputs)



print(docker(R, L, W, B, exhaustiveness=3, run_count=3).run())

