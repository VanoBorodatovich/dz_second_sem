#! python

import re


class Sequences:
    def __init__(self, name, sequence, missing_val='N', gap_val='-'):
        self.sequence = sequence
        self.name = name
        self.missing_val, self.gap_val = missing_val, gap_val
        if len(self) != 0:
            self.completeness = round(self.non_mis_bases() / len(self), 3)
        else:
            self.completeness = None

    def __str__(self):
        return '>' + self.name + '\n' + self.sequence

    def __len__(self):
        return len(self.sequence)

    def make_bigger(self, fromm, to):
        str = ''
        for f in self.sequence:
            if f in fromm:
                ind = fromm.index(f)
                str += to[ind]
            if f not in fromm:
                str += f
        return str

    def check_alphabet(self, alphabet):  # only for subclasses
        unexpected = []
        is_ok = True
        for f in self.sequence:
            if f not in alphabet:
                is_ok = False
                unexpected.append(f)
        if unexpected:
            print('unexpected symbols:' + ' '.join(unexpected) + ' your sequence will not be accepted')
        return is_ok

    def remove_gaps(self):
        new_seq, trimmed = '', None
        for pos in self.sequence:
            if pos != self.gap_val:
                new_seq += pos
        trimmed = Dnaseq(self.name, new_seq)
        return trimmed

    def non_mis_bases(self):
        count = 0
        for f in self.sequence:
            if f != self.missing_val and f != self.gap_val:
                count += 1
        return count


class Dnaseq(Sequences):
    main_forward = re.compile(r"[Aa][Tt][Gg](?:[AaTtGgCc]{3})+?[Tt](?:[Aa][Gg]|[Gg][Aa]|[Aa][Aa])")
    transl_dic = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'ATT': 'I',
                'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S',
                'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'GCT': 'A', 'GCC': 'A',
                'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '', 'TAG': '', 'TGA': '',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
                'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
                'GAG': 'E', 'TGT': 'C', 'TGC': 'C',
                'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'AGT': 'S',
                'AGC': 'S', 'GGT': 'G', 'GGC': 'G',
                'GGA': 'G', 'GGG': 'G'}
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'W': 'W', 'S': 'S', 'R': 'Y', 'Y': 'R', 'K': 'M',
                       'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', '-': '-'}
    deg_nucl_dict = {'S': 'CG', 'R': 'AG', 'Y': 'TC', 'K': 'GT', 'M': 'CA', 'W': 'AT', 'B': 'CGT', 'V': 'CGA', 'D': 'AGT',
                     'H': 'TCA', 'N': 'AGTC', 'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T'}

    def __init__(self, name, sequence, missing_val='N', gap_val='-'):
        super().__init__(name, sequence, missing_val, gap_val)  # нужно передавать все аргументы
        self.alphabet = ['A', 'T', 'C', 'G', 'W', 'S', 'R', 'Y', 'K',
                         'M', 'B', 'V', 'D', 'H', self.missing_val, self.gap_val]
        self.small_alphab = ['a', 't', 'c', 'g', 'w', 's', 'r', 'y', 'k', 'm', 'b', 'v', 'd', 'h']
        self.sequence = Sequences.make_bigger(self, fromm=self.small_alphab, to=self.alphabet)
        if not self.check_alphabet(self.alphabet):
            self.sequence = None

    def show_orfs(self):
        print(self)
        frames, end_before = '', False
        for f in re.finditer(Dnaseq.main_forward, self.sequence):
            if end_before:
                frames += (f.start() - end_before) * ' ' + '}' + (f.end() - f.start() - 2) * '+' + ']'
            else:
                frames += f.start() * ' ' + '}' + (f.end() - f.start() - 2) * '+' + ']'
            end_before = f.end()
        print(frames)

    def get_prot(self):  # ищет прямые не перекрывающиеся рамки
        arr_prots, count = [], 0
        for f in re.finditer(Dnaseq.main_forward, self.sequence):
            protein, i, j = '', 0, 3
            while j <= len(f.group()):
                if f.group()[i:j] not in Dnaseq.transl_dic:
                    protein += '-'
                else:
                    protein += Dnaseq.transl_dic[f.group()[i:j]]
                i += 3
                j += 3
            count += 1
            protein = Protein(self.name + '_v' + str(count), protein)
            arr_prots.append(protein)
        return arr_prots

    def get_complement(self):
        seq, new_name = '', self.name + '_comp'
        for f in self.sequence:
            seq += Dnaseq.complement_dict[f]
        return Dnaseq(new_name, seq)

    def get_reverse_comp(self):
        new_name = self.name + '_rc'
        f = self.get_complement().sequence[::-1]
        return Dnaseq(new_name, f)

    def get_regex(self):
        seq = ''
        for f in self.sequence:
            if f not in Dnaseq.deg_nucl_dict:
                return 'unusual character: ' + f
            elif f in ('A', 'G', 'T', 'C'):
                seq += f
            else:
                seq += '[' + Dnaseq.deg_nucl_dict[f] + ']'
        return f"{seq}"


class Protein(Sequences):
    def __init__(self, name, sequence, missing_val='?', gap_val='-'):
        super().__init__(name, sequence, missing_val, gap_val)
        self.alphabet = ['G', 'F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y', 'H', 'Q', 'N', 'K', 'D', 'C', 'W',
                    'R', 'S', 'E', self.missing_val, self.gap_val]
        if not self.check_alphabet(self.alphabet):
            self.sequence = None


if __name__ == '__main__':
    var = Dnaseq('test', 'ATGCTTAGATGACATGAATCTGAGGAGGATTTTCAGTAGATAAAGCCA')
    var.show_orfs()
    print(var.get_prot())
    for f in var.get_prot():
        print(f)




Check = 'seqs'