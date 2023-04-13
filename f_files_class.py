#! python

import sys
from importlib import reload
sys.path.append(r"C:/Users/vano/PycharmProjects/pythonProject/project1/py/great_work/test")
import seq_class
reload(seq_class)
from seq_class import *


class FastaFile:

    def __init__(self, file: str | None, name=None, arr=None):
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

    def find_in(self, pattern: str):
        arr_pattern, other = [], []
        for f in self.arr:
            if re.search(rf"{pattern}", f.name):
                arr_pattern.append(f)
            else:
                other.append(f)
        return arr_pattern, other

    def write_on_top(self, pattern: str, add_to_name='_reor.fa'):
        arr_ref, other = self.find_in(pattern)
        new_name = self.name + add_to_name
        with open(new_name, 'w') as file:
            for f in range(len(arr_ref)):
                file.write(">" + arr_ref[f].name + "\n")
                file.write(arr_ref[f].sequence + "\n")
            for f in range(len(other)):
                file.write(">" + other[f].name + "\n")
                file.write(other[f].sequence + "\n")

    def split_by_pattern(self, pattern: str):
        arr_ref, other = self.find_in(pattern)
        if len(arr_ref) < 2:
            return "refs are not found"
        for f in range(len(arr_ref)):
            new_name = self.name + '_' + str(f) + '.fa'
            with open(new_name, 'w') as file:
                file.write(">" + arr_ref[f].name + "\n")
                file.write(arr_ref[f].sequence + "\n")
                for i in other:
                    file.write(">" + i.name + "\n")
                    file.write(i.sequence + "\n")

    def get_namespace(self):
        pattern = re.compile(r"^(?:_R_)?([A-Za-z\d]+)")
        dict_names = {}
        for f in self.arr:
            var = re.search(pattern, f.name)
            if var:
                if var.groups()[0] in dict_names:
                    dict_names[var.groups()[0]].append(f)
                else:
                    dict_names[var.groups()[0]] = [f]
            else:
                print('for ', f.name, ' pattern is not found')
        return dict_names

    def get_prots(self):
        prots_arr = []
        for f in self.arr:
            prots_arr += f.get_prot()
        var = FastaFile(None, self.name + '_prots', prots_arr)
        var.write(var.name)
        return True


class Aligning(FastaFile):

    def __init__(self, file: str | None, name=None, arr=None, ref=None):
        super().__init__(file, name, arr)
        print(self.check_length())
        # set the ref
        if not ref:
            self.ref = None
        else:
            set_ref, other = self.find_in(ref)
            if len(set_ref) != 1:
                print('__init__: too many/zero ref')
            else:
                self.ref = set_ref[0]
                self.arr = other

    def check_length(self):
        for f in self.arr:
            if len(self.arr[0]) != len(f):
                del self.arr
                return 'Your aligning has different length inside!!!'

    def view(self):
        print(self.name)
        if self.ref:
            print(">" + self.ref.name + "\n" + self.ref.sequence + '\n' + 'the ref is above\n')
        for f in self.arr:
            print(">" + f.name + "\n" + f.sequence + '\n')

    def __len__(self):  # just len of the alignment
        return len(self.arr[0])

    def non_mis_bases(self):
        count = 0
        for f in self.arr:
            count += f.non_mis_bases()
        self.score = count

    def completeness(self, per=0.5):
        bol = True
        for f in self.arr:
            if f.completeness < per:
                bol = False
            if not bol:
                return False

    @staticmethod
    def count_2seq_sim(first_seq: Dnaseq, second_seq: Dnaseq):
        if len(first_seq) != len(second_seq):  # это можно заменить на исключение
            print("sequences have different length")
            return
        else:
            count_len, n_dif = 1, 0
            for f in range(len(first_seq)):
                if first_seq.sequence[f] != '-' and second_seq.sequence[f] != '-':
                    count_len += 1
                    if first_seq.sequence[f] != second_seq.sequence[f]:
                        n_dif += 1
            return round(n_dif/count_len, 4)

    @staticmethod
    def merge_2seq(first_seq: str, second_seq: str):
        seq = ''
        for pos in range(len(first_seq)):
            if first_seq[pos] != '-' or second_seq[pos] != '-':
                if first_seq[pos] == second_seq[pos]: # если они равны, значит оба не равны "-"
                    seq += first_seq[pos]
                elif second_seq[pos] == '-':
                    seq += first_seq[pos]
                elif first_seq[pos] == '-':
                    seq += second_seq[pos]
                else:
                    seq += 'N'
            else:
                seq += '-'
        return seq

    def check_isoform(self, threshold: float):
        d, arr_problem, everyth_good = self.get_namespace(), [], True
        for isogr in d:
            if len(d[isogr]) == 1:  # if only one isoform
                continue
            for i in range(len(d[isogr])):  # compare each item in a list with the rest
                for j in range(i + 1, len(d[isogr])):
                    if Aligning.count_2seq_sim(d[isogr][i], d[isogr][j]) > threshold:
                        arr_problem.append(d[isogr][i].name + " vs " + d[isogr][j].name)
                        everyth_good = False
        return everyth_good, arr_problem

    def merge_isoforms(self):
        d, new_arr = self.get_namespace(), []
        for isogr in d:  # iterate over names of gr (keys)
            if len(d[isogr]) == 1:  # if only one isoform
                new_arr.append(d[isogr][0])
                continue
            seq = d[isogr][0].sequence
            for isof in range(1, len(d[isogr])):
                seq = Aligning.merge_2seq(seq, d[isogr][isof].sequence)
            new_arr.append(Dnaseq(isogr, seq))
        return new_arr

    def get_coverage(self):
        cov = []
        for pos in range(len(self.arr[0])):
            cov_pos = 0
            for seq in self.arr:
                if seq.sequence[pos] != '-':
                    cov_pos += 1
            cov.append(cov_pos)
        return tuple(cov)

    def crop_by_cov(self, cov_threshold: int):
        cov, new_arr = self.get_coverage(), []
        for seq in self.arr:
            new_seq = ''
            for pos in range(len(cov)):
                if cov[pos] >= cov_threshold:
                    new_seq += seq.sequence[pos]
            new_arr.append(Dnaseq(seq.name, new_seq))
        return new_arr

    def write_in_order(self, new_name):  # just write the existing object
        with open(new_name, 'w') as file:
            if self.ref:
                file.write(">" + self.ref.name + "\n")
                file.write(self.ref.sequence + "\n")
                for i in self.arr:
                    file.write(">" + i.name + "\n")
                    file.write(i.sequence + "\n")
            else:
                self.write(new_name)
        return True

    def trim_by_ref(self):  # trim all the sequences in the borders of nonmissing values of ref
        if not self.ref:
            return "ref is not set"
        new_seqs = []
        for query in self.arr:
            seq = ''
            for pos in range(len(self.ref.sequence)):
                if self.ref.sequence[pos] != '-':
                    seq += query.sequence[pos]
            new_seqs.append(Dnaseq(query.name, seq))
        new_ref = self.ref.remove_gaps()
        return new_seqs, new_ref

    def count_sim_to_ref(self, query: Dnaseq):
        if not self.ref:
            return "the ref is not set((( "
        else:
            count = 0
            for pos in range(len(self.ref.sequence)):
                if self.ref.sequence[pos] == query.sequence[pos]:
                    count += 1
            return count

    def choose_one(self):
        if not self.ref:
            return "the ref is not set((( "
        else:
            namespace = self.get_namespace()
            arr_elect = []
            for f in namespace:
                count, elect = 0, None
                for i in namespace[f]:  # for isoforms inside sample
                    current_count = self.count_sim_to_ref(i)
                    if count <= current_count:
                        count = current_count
                        elect = i
                arr_elect.append(elect)
            return arr_elect
        # then you just need reassign self.arr to the returned arr and use the write_on function


CHECK = 'f_file'