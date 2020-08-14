import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from openbabel import pybel


def writeFile(mol_in, format, f_out):

    if type(mol_in) is not list: mol_in = [mol_in]
    writer = pybel.Outputfile(format, f_out)

    for molecule in mol_in:
        writer.write(molecule)

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


class PlFactory(Mol):

    def __init__(self, name, r1, r2=''):
        self.name = name
        self.r1 = r1.__str__()
        self.r2 = r2.__str__()[::-1].replace('O)O=(C', 'C(=O)')
        self.head = {'pc':'CC(COP(=O)([O-])OCC[N+](C)(C)C)O',
                     'pe':'CC(COP(=O)(O)OCCN)O',
                     'pg':'CC(COP(=O)(O)OCC(CO)O)O',
                     'pi':'CC(COP(=O)(O)OC1C(C(C(C(C1O)O)O)O)O)O',
                     'sm':'(C(COP(=O)([O-])OCC[N+](C)(C)C)N'
                     }

    def construct(self):
        if self.name == 'sm':
            return str(self.r1[:-5]
                       + self.head[self.name]
                       + self.r2
                       + ')O'
                       )

        else:
            return str(self.r1
                       + self.head[self.name]
                       + self.r2
                       )

    def __str__(self):
        return self.construct()


class LeaFactory(Mol):
    '''
    Parametrizes LEA structure generation
    '''
    def __init__(self, r1):
        self.r1 = r1.__str__()
        self.head = 'OCC(CNC(C)C)O'

    def __str__(self):
        return str(self.r1 + self.head)

