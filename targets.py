import os
import pathlib
from multiprocessing import Process, Queue, current_process
import threading
from more_itertools import consume
from openbabel import pybel
from pymol import cmd

aa3to1={
   'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
   'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
   'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
   'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
   'MSE':'M',
}

class Multiprocess:
    """
    to parallelize conversion of Protein files:
    Initialize with function and buffer size
    call with list of Protein arguments to clean/convert
    """

    def __init__(self, function, sema_buffer=4, q_size=0):
        self.function = function
        self.sema = threading.BoundedSemaphore(sema_buffer)
        self.q = Queue(maxsize=q_size)

    def __call__(self, *args, **kwargs):
        consume(self.q.put(arg) for arg in args)
        with self.sema:
            while self.q.qsize() > 0:
                a = self.q.get()
                p = Process(target=self.wrap, args=(a,))
                p.start()

    def wrap(self, arg):
        self.sema.acquire()
        print('In worker process: ',
              "\n", arg,
              "\n", current_process()
        )
        self.function(arg)
        print('done working: ', current_process())
        self.sema.release()


class Protein:
    """
    Protein class instantiated using path
    to .pdb, .pdbqt, or other structure file.
    Implements methods to clean and convert
    structure files, with parallel
    execution when used with Multiprocess.
    """

    def __init__(self, path, name = ''):

        self.path = pathlib.Path(path)
        self.pdb_path = self.path if self.path.suffix == '.pdb' else ''
        self.pdbqt_path = self.path if self.path.suffix == '.pdbqt' else ''
        self.name = name if name else self.path.stem

    @classmethod
    def to_fasta(cls, pdb_path: str) -> dict:

        ca_pattern=re.compile(
            '^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])'\
            '^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])')
        chain_dict=dict()
        chain_list=[]
        fp=open(pdb_path,'r')

        for line in fp.read().splitlines():
            if line.startswith("ENDMDL"):
                break
            match_list=ca_pattern.findall(line)
            if match_list:
                resn=match_list[0][0]
                chain=match_list[0][1]
                if chain in chain_dict:
                    chain_dict[chain]+=aa3to1[resn]
                else:
                    chain_dict[chain]=aa3to1[resn]
                    chain_list.append(chain)

        fp.close()

        return chain_dict

    def clean(self):
        
        if self.pdbqt_path:
            to_name = self.path.with_suffix('.pdbqt')
            with open(self.pdbqt_path, "r") as f:
                lines = f.readlines()
            with open(to_name, "w") as f:
                for line in lines:
                    if line.strip("\n")[:4] == "ATOM":
                        f.write(line)
                f.close()
            self.path = to_name
        else:           
            self.path_clean = self.path.with_suffix('.clean.pdb')
            with open(self.path, "r") as f:
                lines = f.readlines()
            with open(self.path_clean, "w") as f:
                for line in lines:
                    if line.strip("\n")[:4] == "ATOM":
                        f.write(line)
                f.close()
        
        return self

    def convert(self):

        if self.path.suffix == '.pdbqt':
            return self.path        
        else:
            self.pdbqt_path = self.path_clean.with_suffix('.pdbqt')
            cmd.load(self.path_clean)
            cmd.remove('resn HOH')
            cmd.h_add(selection='acceptors or donors')
            cmd.save(self.path_clean)
            mols = list(pybel.readfile('pdb', str(self.path_clean)))
            writer = pybel.Outputfile(
                'pdbqt', str(self.pdbqt_path), opt={'pdbqt': '-xh'}
            )

            for molecule in mols:
                writer.write(molecule)
                writer.close()
            cmd.reinitialize()
        
        return self

    def write_pdbqt(self):

        self.clean().convert().clean()
        os.remove(self.path_clean)
        os.remove(self.pdbqt_path)
        self.pdbqt_path = self.path

        return self

    def __repr__(self):

        return str(self.path)
