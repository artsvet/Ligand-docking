import os
import pathlib
from multiprocessing import Process, Queue, current_process
import threading
from more_itertools import consume
from openbabel import pybel
from pymol import cmd

"""Protein files are published as Protein Data Bank files (pdb)
and require additional editing to use with autodock vina"""

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


class Protein(type(pathlib.Path())):
    """
    Protein class reading path
    to .pdb, .pdbqt, or other structure file.
    Implements methods to clean and convert
    structure files, with parallel
    execution when used with Multiprocess.
    """

    def __init__(self, path):

        super().__init__()

    def clean(self):
        """Removes metadata/comments"""

        self.path_clean = Protein(
            self.with_name(
                self.stem + '.clean' + self.suffix
            )
        )
        with open(self, "r") as f:
            lines = f.readlines()
        with open(self.path_clean, "w") as f:
            for line in lines:
                if line.strip("\n")[:4] == "ATOM":
                    f.write(line)
            f.close()
        os.remove(self)

        return self.path_clean

    def convert(self):
        """adds Gasteiger partial charges and AutoDock atom-types with pybel"""

        self.path_pdbqt = Protein(self.path_clean.with_suffix('.pdbqt'))
        cmd.load(self.path_clean.__str__())
        cmd.remove('resn HOH')
        cmd.h_add(selection='acceptors or donors')
        cmd.save(self.path_clean.__str__())
        atoms = list(pybel.readfile('pdb', self.path_clean.__str__()))
        writer = pybel.Outputfile(
            'pdbqt', self.path_pdbqt.__str__(), opt={'pdbqt': '-xh'}
        )
        for atom in atoms:
            writer.write(atom)
            writer.close()
        cmd.reinitialize()
        os.remove(self)

        return self.path_pdbqt

    def prepare(self):
        """Gets .pdb file ready for docking"""

        return self.clean().convert().clean()

    def __repr__(self):

        return self.root
