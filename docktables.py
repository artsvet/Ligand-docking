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
        data = pd.DataFrame()
        for t in targets:
            data = data.join(pd.concat([pd.DataFrame(data = [df.mean(), df.std()], columns = ['Mean_affinity_'+ t, 'Stdev_affinity_' + t])
                                        for df in self.all.groupby('Ligand')['Affinity']]
                             )
                   )
                  
        return data
