#! python

from f_files_class import *
import click


@click.group()
def main():
    pass


@main.command()
@click.option('--filename', '-f', type=str, help='Path to file')
def get_sequences(filename):
    fasta = FastaFile(filename)
    fasta.view()


@main.command()
@click.option('--filename', '-f', type=str, help='Path to file')
def find_proteins(filename):
    fasta = FastaFile(filename)
    fasta.get_prots()


def get_ampl_coor(ref: Dnaseq, first_seq: Dnaseq, second_seq: Dnaseq):
    # create patterns
    first_p, second_p = re.compile(first_seq.get_regex()), re.compile(Dnaseq.get_regex(second_seq.get_reverse_comp()))
    # search amplicons
    pos_dict = {}
    for forw in re.finditer(first_p, ref.sequence):
        start_pos = forw.start()
        for rev in re.finditer(second_p, ref.sequence):
            end_pos = rev.end()
            if 0 < end_pos - start_pos < 24:
                print('your primers are overlapping')
                return
            elif 24 < end_pos - start_pos < 2000:
                pos_dict[tuple([forw.group(), rev.group()])] = tuple([start_pos, end_pos])
    return pos_dict


def write_amplicons(output: str, direction: str, dic: dict):
    with open(output, 'a') as file:
        file.write(f'{direction}\n')
        if len(dic) == 0:
            file.write('nothing_found\n')
        else:
            file.write('pos_start\tpos_end\tlen\tsequence\n')
            for f in dic:
                file.write(f'{dic[f][0]}\t{dic[f][1]}\t{dic[f][1]-dic[f][0]}\t{f}\n')


@main.command()
@click.option('--reference_file_name', '-ref', type=str, help='Path to reference file')
@click.option('--primers_file_name', '-pr', type=str, help='Path to primers file')
@click.option('--output', '-o', type=str, default='amplicons.txt', help='Path to output')
def get_amplicons(reference_file_name: str, primers_file_name: str, output: str):
    # read the files
    primers = FastaFile(primers_file_name)
    if len(primers.arr) != 2:
        print('must be two sequences in primer_file')
        return
    ref = FastaFile(reference_file_name)
    if len(ref.arr) != 1:
        print('must be exactly one ref')
        return
    # get the coordinates
    first_direction_amp = get_ampl_coor(ref.arr[0], primers.arr[0], primers.arr[1])
    second_direction_amp = get_ampl_coor(ref.arr[0], primers.arr[1], primers.arr[0])
    # write
    with open(output, 'w') as file:
        file.write('rabota\n')
    write_amplicons(output, 'first_direction', first_direction_amp)
    write_amplicons(output, 'second_direction', second_direction_amp)


if __name__ == '__main__':
    main()