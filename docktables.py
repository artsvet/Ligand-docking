import pandas as pd
from vinaDock import Docker

class Dock_Report_Factory:
    def __init__(self, docks):
        self.docks = docks
        self.all = pd.concat([dock.out for dock.out in self.docks])

    def best(self):
        best_affinities = pd.concat([row for dock in self.docks 
                                     for i, row in dock.out.iterrows() 
                                     if row['Rank'] == 1])
        
        return best_affinities
        
    def summary(self):
        targets = self.all['Receptor'].unique()
        df = pd.DataFrame()
        for t in targets:
            df = pd.concat([pd.Series(df.groupby('Ligand')['Affinity'].mean()).rename('Mean_affinity_' + t), 
                            pd.Series(df.groupby('Ligand')['Affinity'].std()).rename('Standard_deviation_' + t)
                            ], axis = 1)
                  
        return data
