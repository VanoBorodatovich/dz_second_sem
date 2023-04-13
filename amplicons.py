#! python
import re

from f_files_class import *


# ref = FastaFile('ref.fa')
# if len(ref.arr) > 1:
#     print('your ref is hard')
forw = Dnaseq('forward', 'CTYATT')
rev = Dnaseq('reverse', 'ATCTYT')
seq = 'ACTCATTGGTCTGGTACTTATTAGTTTACTGAAGTACCAAAGAT'
# seq = ref.arr[0].sequence


def get_ampl_regex(forward: Dnaseq, reverse: Dnaseq):
    return re.compile(rf'.{forward.get_regex()}[AGTC]+?{Dnaseq.get_regex(reverse.get_reverse_comp())}')



def get_amplicons(amp_pattern: re.Pattern, ref):
    return re.findall(get_ampl_regex(forw, rev), ref)




with open('log_amp.txt', 'w') as log:
    log.write('start' + '\t' + 'end' + '\t' + 'length' + '\t' + 'orientation' + '\n')
print(get_ampl_regex(forw, rev))
print(get_amplicons(forw, rev, seq))

def rec_search(pattern, seq):
    for f in re.finditer(pattern, seq):
        print(f.groups()[0])
        rec_search(pattern, f.groups()[0])