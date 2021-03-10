>>> class targetsParse:
...     """
...     Parses Series to set up dependencies for ligand
...     docking. Contains Protein target attributes,
...     specifies ligands that will be docked using
...     reference lipid and pl tables. Will serve as
...     primary interface for running dock and
...     collecting outputs.
...     """
... 
...     def __init__(self, targetsSeries, lipidTable, plTable, selectTable):
...         self.lipidTable = lipidTable
...         self.plTable = plTable
...         self.selectTable = selectTable
...         self.id = targetsSeries['id']
...         self.pdb = Protein(targetsSeries['pdb'])
...         self.name = targetsSeries['name']
...         self.type = targetsSeries['type']
...         self.species = targetsSeries['species']
...         self.function = targetsSeries['function']
...         self.ligands = re.split(':|,', targetsSeries['lipid'])
...         self.ligandsIndexes = [int(i) for i in self.ligands[1:]]
...         self.selected = [
...             int(i) for i in re.split(',', targetsSeries['selected'])
...             if i not in [',', '']
...         ]
...         self.box = eval(targetsSeries['box'])
... 
...     def build_lea_list(self):
...         lea_list = []
...         for index, series in self.lipidTable.iterrows():
...             lea_list.append(LeaFromSeries(series))
...         return lea_list
... 
...     def build_ligand_list(self):
...         ligandList = []
... 
...         if self.ligands[0] == 'all':
...             ligandPatterns = [key for key in PL_PATTERNS]
...         elif self.ligands[0] == '':
...             return ligandList
...         else:
...             ligandPatterns = [self.ligands[0]]
... 
...         if self.ligandsIndexes:
...             if ligandPatterns == 'lipid':
...                 lipidSeries = [
...                     self.lipidTable.iloc[i] for i in self.ligandsIndexes
...                 ]
...             else:
...                 plSeries = [self.plTable.iloc[i] for i in self.ligandsIndexes]
... 
...         else:
...             plSeries = [series for index, series in self.plTable.iterrows()]
... 
...         for pattern in ligandPatterns:
...             if pattern == 'lipid':
...                 for series in lipidSeries:
...                     ligandList.append(LipidFromSeries(series))
...             else:
...                 for series in plSeries:
...                     ligandList.append(PlFromSeries(
...                         pattern, series, self.lipidTable
...                     ))
... 
...         return ligandList
... 
...     def build_select_list(self):
...         ligandList = []
...         for index in self.selected:
...             ligand = Mol(self.selectTable.iloc[index]['smiles'])
...             setattr(ligand, 'name', selectedTable.iloc[index]['id'])
...             ligandList.append(ligand)
...         return ligandList
... 
...     def write_ligands(self):
...         ligands = self.build_lea_list() + self.build_ligand_list() + self.build_select_list()
...         root = os.getcwd()
...         dir = os.path.join(
...             root, self.name, self.id
...         )
...         os.makedirs(dir)
...         os.chdir(dir)
...         for mol in ligands:
...             mol.write_pdb()
...             os.chdir(dir)
...         os.chdir(root)
... 
...     def dock(self, ligands, run_count=1, exhaustiveness=10, write_dir=os.getcwd()):
...         outputs = []
...         for mol in ligands:
...             outputs.append(
...                 Docker(receptor=self.pdb,
...                        ligand=mol,
...                        log_path=self.id + '_' + mol.name + '_log.txt',
...                        box=self.box, run_count=run_count,
...                        exhaustiveness=exhaustiveness).run())
...             pd.concat(outputs).to_csv(write_dir + '\\' + self.name
...                                       + '_' + self.id + '_outputs.csv')
... 
...         return outputs