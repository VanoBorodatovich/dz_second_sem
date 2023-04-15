#! python

import sys
# from importlib import reload
# sys.path.append(r"C:/Users/vano/PycharmProjects/pythonProject/project1/py/great_work/test")
# import seq_class
# reload(seq_class)
from seq_class import *


class FastaFile:

    def __init__(self, file: str | None, name=None, arr=None):  # инициализирует только фасты с ДНА
        if file:
            self.file = file.rstrip()
            self.name = re.search(r"/?([^/]+)\.fa", file).groups()[0]
            self.arr = []

            with open(self.file) as file:
                file = file.read()
            seq_names = re.findall(r">(.+)\n", file)
            pre_seq = re.findall(rf">.+\n([^>]+)", file, re.MULTILINE)  # take all the characters
            if len(seq_names) != len(pre_seq):
                self.arr = None  # потом по этому нону можно делать проверки
                print('smt wrong with your fasta: ', self.name)
            else:
                seq = []
                for f in pre_seq:
                    seq.append(re.sub('\n', '', f))
                for f in range(len(seq_names)):
                    self.arr.append(Dnaseq(seq_names[f], seq[f]))
        else:  # if not file
            self.name = name
            if self.check_arr(arr):
                self.arr = arr
            else:
                self.arr = None

    @staticmethod
    def check_arr(arr: list):  # check are all the sequences the same type (DNA or Protein)
        b, sim = True, arr[0]
        for f in arr:
            if not isinstance(f, (Dnaseq, Protein)):
                print('your arr have not sequences object')
                return False
            if type(sim) != type(f):
                print('your arr has different object')
                return False
            sim = f
        return True

    def view(self):
        print(self.name)
        for f in self.arr:
            print(">" + f.name + "\n" + f.sequence + '\n')

    def seqs_len(self):
        for f in self.arr:
            print(">" + f.name + "\n" + str(len(f)))

    def write(self, new_name):
        with open(new_name, 'w') as file:
            for i in self.arr:
                file.write(">" + i.name + "\n")
                file.write(i.sequence + "\n")

    def get_prots(self):
        prots_arr = []
        for f in self.arr:
            prots_arr += f.get_prot()
        var = FastaFile(None, self.name + '_prots', prots_arr)
        var.write(var.name)
        return True
