import os
import pathlib
from multiprocessing import Process, Queue, current_process
import threading
from more_itertools import consume
from openbabel import pybel


class Multiprocess:
    """
    Initialize with function and buffer size
    call with list of arguments to clean/convert
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

        if hasattr(self, 'path_clean'):
            return self.path_clean
        else:
            self.path_clean = Protein(
                self.with_name(
                    self.stem + '.clean' + self.suffix
                )
            )
            setattr(self.path_clean, 'path_clean', self.path_clean)
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

        if hasattr(self, 'path_pdbqt'):
            return self.path_pdbqt
        else:
            self.path_pdbqt = Protein(self.with_suffix('.pdbqt'))
            setattr(self.path_pdbqt, 'path_pdbqt', self.path_pdbqt)
            cmd.load(self.__str__())
            cmd.remove('resn HOH')
            cmd.h_add(selection='acceptors or donors')
            cmd.save(self.__str__())
            mols = list(pybel.readfile('pdb', self.path_clean.__str__()))
            writer = pybel.Outputfile(
                'pdbqt', self.path_pdbqt.__str__(), opt={'pdbqt': '-xh'}
            )
            for molecule in mols:
                writer.write(molecule)
                writer.close()
            cmd.reinitialize()
        os.remove(self)

        return self.path_pdbqt

    def prepare(self):

        if hasattr(self, 'path_pdbqt'):
            self.path_pdbqt.clean()
        elif hasattr(self, 'path_clean'):
            self.path_clean.convert().clean()
        else:
            self.clean().convert().clean()

    def __repr__(self):

        return self.root
