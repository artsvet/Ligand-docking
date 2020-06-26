from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

class LeaFactory:
    '''
    Parametrizes LEA structure generation
    '''

    def __init__(self, r1, r2='OCC(CNC(C)C)',
                 r3='O', E = [], Z = [],
                ):
        self.r1 = 'C' * r1
        self.r2 = r2
        self.r3 = r3
        self.E = E
        self.Z = Z


    def to_mol(self):
        r4 = self.r1
        for i in sorted(self.E, reverse=True):
                r4 = r4[:i-1] + '\\C=C/' + r4[i+1:]
        for i in sorted(self.Z, reverse=True):
                r4 = r4[:i-1] + '\\C=C\\' + r4[i+1:]
        print(r4[::-1] + self.r2 + self.r3)
        return Chem.MolFromSmiles(r4[::-1] + self.r2 + self.r3)

