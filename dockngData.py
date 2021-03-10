import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from statistics import mean, stdev

R = .0019872 #kcal mol^-1 K^-1
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

def scrape_log(text_log):
    list_out = []
    with open(text_log, 'r') as log:
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
    return pd.DataFrame(list_out)


class SimpleBars:

    def __init__(self, title, x_label, y_label, labels, values):
        self.title, self.x_label, self.y_label = title, x_label, y_label
        self.labels, self.values = labels, values
        assert len(self.labels) == len(self.values)
        self.means = [np.mean(arr) for arr in self.values]
        self.stdev = [np.std(arr) for arr in self.values]
        self.bot_bar = tuple(0 for n in self.stdev)
        self.top_bar = tuple(n/2 for n in self.stdev)

    def plot(self):

        x_pos = np.arange(len(self.labels))
        fig, ax = plt.subplots()
        ax.bar(x_pos, self.means,
               yerr=[self.bot_bar, self.top_bar],
               align='center', alpha=0.5, ecolor='black', capsize=10)
        ax.set_ylabel(self.y_label)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(self.labels)
        ax.set_title(self.title)
        ax.yaxis.grid(True)

        # Save the figure and show
        plt.tight_layout()
        plt.savefig('bar_plot_with_error_bars.png')
        plt.show()

class DockInstance:

    def __init__(self, df):
        self.df = df.assign(Kd = KD(df['Affinity']))
        self.top_score = self.df.loc[df['Rank']==1]
        self.top_score_dG = self.top_score['Affinity']
        self.top_score_Kd = self.top_score['Kd']

class DockLigand:

    def __init__(self, df):
        self.df = df
        self.name = self.df['Ligand'].iloc[0]
        self.dock_instances = [DockInstance(i)
                               for t, i in self.df.groupby(['Date_time'])
                               ]
        self.top_scores = [i.top_score
                           for i in self.dock_instances
                           ]
        self.top_scores_dG = [i.top_score_dG.item()
                              for i in self.dock_instances
                              ]
        self.top_scores_Kd = [i.top_score_Kd.item()
                              for i in self.dock_instances
                              ]

class DockTarget:

    def __init__(self, df, targets_table):
        self.df = df
        self.targets_series = targets_table.loc[
            targets_table['id'] == self.df.Receptor[0]
            ].squeeze()
        self.name = self.targets_series['name']
        self.dock_ligands = {name: DockLigand(ligand)
                             for name, ligand in self.df.groupby(['Ligand'])
                             }

    def best_n(self, n):

        best_n = sorted([self.dock_ligands[key] for key in self.dock_ligands],
                        key=lambda x: mean(x.top_scores_dG)
                        )[:n]

        return best_n

    def get_default_ligands(self):
        try:
            ref_ligands = [ligand
                           for key in re.split(', ', self.targets_series['lipidPatterns'])
                           if key != 'lipid'
                           for ligand in DEFAULT_LIGANDS[key]
                           ]
            defaults_list = [self.dock_ligands[name]
                             for name in self.dock_ligands
                             if name in ref_ligands
                             ]
        except:
            return [self.dock_ligands[ligand]
                    for ligand in self.dock_ligands
                    if ligand in DEFAULT_LEAS]

        return [self.dock_ligands[ligand]
                for ligand in self.dock_ligands
                if ligand in defaults_list + DEFAULT_LEAS]

    def build_report(self, best_n=False):
        if not best_n:
            best_n = len(self.dock_ligands)

        defaults_list = self.get_default_ligands()

        defaults = pd.DataFrame(
            data=[(ligand.name, *ligand.top_scores_Kd,
                   mean(ligand.top_scores_Kd), stdev(ligand.top_scores_Kd)
                   ) for ligand in defaults_list
                  ],
            columns=['{0} LEA and reference Ligands'.format(self.name),
                     'Run 1', '2', '3', 'Mean Kd (moles/L)', 'StDev Kd']
        )

        best = pd.DataFrame(
            data=[(ligand.name, *ligand.top_scores_Kd,
                   mean(ligand.top_scores_Kd), stdev(ligand.top_scores_Kd)
                   ) for ligand in self.best_n(best_n)
                  ],
            columns=['{0} best {1} scores'.format(self.name, best_n),
                     'Run 1', '2', '3', 'Mean Kd (moles/L)', 'StDev Kd']
        )
        if defaults.empty:
            return [best]
        else:
            return [defaults, best]

    def report_bars(self, best_n=3):

        defaults_list = self.get_default_ligands()

        to_plot = defaults_list + [ligand for ligand in self.best_n(best_n)
                                   if ligand not in defaults_list]

        title = '{0} (PDB:{1}) ligand binding affinities'.format(
            self.targets_series['name'], self.targets_series['id'],)
        x_label = 'Ligand'
        y_label = 'Dissociation constant (Kd, mol/L)'
        labels = [ligand.name for ligand in to_plot]
        values = [ligand.top_scores_Kd for ligand in to_plot]

        return SimpleBars(title, x_label, y_label, labels, values)
    


targets = pd.read_csv(str(os.getcwd() + '\\targetsTableNew.csv')).fillna('')
sbir_reports = 'C:\\Users\\starch\\PycharmProjects\\molDocking\\finished docks\\'
report_df = [(stan_reports + file, pd.read_csv(stan_reports + file))
             for file in os.walk(stan_reports).__next__()[2]
             ]

reports = [(path, DockTarget(df, targets)) for path, df in report_df]
to_summary = [(path, target.build_report())
              for path, target in reports
              ]

export_csv = [df.to_csv(path[:-4] + '_summary.csv', mode='a')
              for path, dfs in to_summary
              for df in dfs
              ]

plots = []


