import os
import pathlib
from multiprocessing import Process, Queue, current_process
import threading
from more_itertools import consume
from openbabel import pybel
from pymol import cmd


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

    def __init__(self, path, name=''):
        self.path = pathlib.Path(path)
        self.name = name

    def clean(self):

        if hasattr(self, 'path_clean'):
            return self.path_clean
        else:
            self.path_clean = self.path.with_name(
                    self.path.stem + '.clean' + self.path.suffix
                )
            with open(self.path, "r") as f:
                lines = f.readlines()
            with open(self.path_clean, "w") as f:
                for line in lines:
                    if line.strip("\n")[:4] == "ATOM":
                        f.write(line)
                f.close()

        os.remove(self.path)

        if self.path.suffix == '.pdbqt':
            os.rename(self.path_clean, self.path.stem + self.path.suffix)

        self.path = self.path_clean

        return self.path

    def convert(self):

        if self.path.suffix == '.pdbqt':
            return self.path
        
        else:
            self.path_pdbqt = self.path.with_suffix('.pdbqt')
            cmd.load(self.path_clean)
            cmd.remove('resn HOH')
            cmd.h_add(selection='acceptors or donors')
            cmd.save(self.path_clean)
            mols = list(pybel.readfile('pdb', self.path_clean))
            writer = pybel.Outputfile(
                'pdbqt', self.path_pdbqt, opt={'pdbqt': '-xh'}
            )
            for molecule in mols:
                writer.write(molecule)
                writer.close()
            cmd.reinitialize()

        os.remove(self.path)
        delattr(self, path_clean)
        self.path = self.path_pdbqt

        return self.path

    def prepare(self):

        if self.path.suffix == '.pdbqt':
            self.clean()
        else:
            self.clean().convert().clean()

        return self.path

    def __repr__(self):

        return self.path
