import os
import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from openbabel import pybel

plPatterns = {
    'pa':'CC(COP(=O)(O)O)O',
    'pc':'CC(COP(=O)([O-])OCC[N+](C)(C)C)O',
    'pe':'CC(COP(=O)(O)OCCN)O',
    'pg':'CC(COP(=O)(O)OCC(CO)O)O',
    'pi':'CC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)O',
    'ps':'CC(COP(=O)(O)OCC(C(=O)O)N)O',
    'sm':'(C(COP(=O)([O-])OCC[N+](C)(C)C)N',
    'cm':'(C(CO)N'
    }

def writeFile(mol_in, format, f_out):

    if type(mol_in) is not list:
    mol_in = [mol_in]
    writer = pybel.Outputfile(format, f_out)

    for molecule in mol_in:
        writer.write(molecule)
        writer.close()

class Mol():

    def toChem(self):
        return Chem.MolFromSmiles(self.__str__())

    def showMol(self):
        return Draw.ShowMol(Chem.MolFromSmiles(self.__str__()), size=(600, 300))

    def toImage(self):
        return Draw.MolToImage(self.__str__(), size=(600, 300))

    def toPybel(self):
        return pybel.readstring('smi', self.__str__())

class Lipid(Mol):

    def __init__(self, len, E=[], Z=[], acid=True):
        self.len = len
        self.E = E
        self.Z = Z
        self.acid = acid

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
        if self.acid == True:
            return chain[::-1] + '(=O)O'
        else:
            return chain[::-1]

    def __str__(self):
        return self.construct()

class LeaFactory(Mol):
    '''
    Parametrizes LEA structure generation
    '''
    def __init__(self, r1):
        self.r1 = r1.__str__().replace('(=O)O','')

    def __str__(self):
        return str(self.r1 + 'OCC(CNC(C)C)O')

class PlFactory(Mol):

    def __init__(self, name, r1='', r2=''):
        self.name = name
        self.r1 = r1.__str__()
        self.r2 = r2.__str__()[::-1].replace('O)O=(C', 'C(=O)')


    def construct(self):
        if self.name == 'sm' or self.name == 'cm':
            self.r1 = 'CCCCCCCCCCCCCC=CC'
            return str(self.r1
                       + plPatterns[self.name]
                       + self.r2
                       + ')O'
                       )

        else:
            return str(self.r1
                       + plPatterns[self.name]
                       + self.r2
                       )

    def __str__(self):
        return self.construct()

class LipidFromSeries(Lipid):

    def __init__(self, series):
        print(series)
        self.name = series['id']
        self.len = series['length']
        self.E = [int(E)
                  for E in series['E'].split(',')
                  if E != ''
                  ]
        self.Z = [int(Z) for Z in series['Z'].split(',') if Z != '']
        self.acid=True

class LeaFromSeries(LeaFactory):

    def __init__(self, series):
        self.name = series['id']
        self.r1 = LipidFromSeries(series).__str__()
        super().__init__(self.r1)

class PlFromSeries(PlFactory):
    def __init__(self, name, series, lipidTable):
        self.name = name
        self.r1 = LipidFromSeries(
            lipidTable.loc[
                lipidTable['id'] == series['r1']
                ].squeeze()
            ).__str__()
        if series['r2']:
            self.r2 = LipidFromSeries(
                lipidTable.loc[
                    lipidTable['id'] == series['r2']
                    ].squeeze()
                ).__str__()[::-1].replace('O)O=(C', 'C(=O)')
        else:
            self.r2 = ''
        super().__init__(self.name, self.r1, self.r2)

class targetsParse:
    def __init__(self, targetsSeries, lipidTable, plTable):
        self.lipidTable = lipidTable
        self.plTable = plTable
        self.id = targetsSeries['id']
        self.name = targetsSeries['name']
        self.type = targetsSeries['type']
        self.species = targetsSeries['species']
        self.function = targetsSeries['function']
        self.ligands = re.split(':|,',targetsSeries['lipid'])
        self.ligandsIndexes = self.ligands[1:]

    def buildLeaDict(self):
        leaDict = {}
        for index, series in self.lipidTable.iterrows():
            leaDict[series['id']] = LeaFromSeries(series)
        return leaDict

    def buildLigandDict(self):
        ligandDict = {}

        if self.ligands[0] == 'all':
            ligandNames = [key for key in plPatterns]
        else:
            ligandNames = [self.ligands[0]]

        if self.ligandsIndexes:
            if ligandNames == 'lipid':
                lipidSeries = [self.lipidTable.iloc[i] for i in self.ligandsIndexes]
            else:
                plSeries = [self.plTable.iloc[i] for i in self.ligandsIndexes]

        else:
            plSeries = [series for index, series in self.plTable.iterrows()]

        for name in ligandNames:
            if name == 'lipid':
                for series in lipidSeries:
                    ligandDict[self.lipidTable['id']] = LipidFromSeries(series)
            else:
                for series in plSeries:
                    ligandDict[name + ' ' + series['id']] = PlFromSeries(
                    name, series, self.lipidTable
                )

        return ligandDict

    def compileLigands(self):
        ligands = self.buildLeaDict(self.lipidTable) + self.buildLigandDict(
            self.lipidTable, self.plTable
        )
        return ligands

    def convertLigands(self):
        for key in (x := self.compileLigands()):
            x[key] = x[key].toPybel()


targets = pd.read_csv(str(os.getcwd() + '\\targetsTable.csv')).fillna('')
lipid = pd.read_csv(str(os.getcwd() + '\\ligands\\lipid.csv')).fillna('')
phospholipid = pd.read_csv(str(os.getcwd() + '\\ligands\\phospholipid.csv')).fillna('')
selected = pd.read_csv(str(os.getcwd() + '\\ligands\\selected.csv')).fillna('')