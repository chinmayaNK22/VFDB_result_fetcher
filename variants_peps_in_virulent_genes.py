from itertools import islice
import read_fasta_file

variant_file = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_variants_identified_final_complete_info.txt"

virulent_fasta = "VFDB_Mavium_virulent_Proteins_nucleotide_sequences_Oct_18-9903385365.ffn"

virulent_file = "VFDB_Mavium_Virulence_Factor_Oct_16-79629633_Final.txt"


def get_header(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            return i.rstrip().split('\t')


dicts = {}
with open(variant_file) as file:
    for i in islice(file, 1, None):
        split_i = i.rstrip().split('\t')
        pep = split_i[0].split('.')[1].upper()
        snp_pos = split_i[10]
        if snp_pos not in dicts:
            dicts[snp_pos] = [split_i]
        else:
            dicts[snp_pos].append(split_i)

virulents = {}
with open(virulent_file) as file:
    for i in islice(file, 1, None):
        split_i = i.rstrip().split('\t')
        orfs = split_i[4].split(';')
        for orf in orfs:
            if orf.lstrip() not in virulents:
                virulents[orf.lstrip()] = [split_i]
            else:
                virulents[orf.lstrip()].append(split_i)

variant_vir = {}
for rows in read_fasta_file.read_fasta(virulent_fasta):
    header = rows[0]
    seq = rows[1]
    header_split = header.split(' ')
    orf_start = int(header_split[-4])
    orf_end = int(header_split[-3])
    if (orf_end - orf_start) > 0:
        for k, v in dicts.items():
            if int(k) >= orf_start and int(k) <= orf_end:
                for j in v:
                    if seq[int(k)-orf_start] == j[11]:
                        if header_split[0].split('\t')[0] not in variant_vir:
                            variant_vir[header_split[0].split('\t')[0]] = [j]
                        else:
                            variant_vir[header_split[0].split('\t')[0]].append(j)

    else:
        for k, v in dicts.items():
            if int(k) >= orf_end  and int(k) <= orf_start:
                for j in v:
                    if seq[int(k)-orf_end] == j[11]:
                        if header_split[0].split('\t')[0] not in variant_vir:
                            variant_vir[header_split[0].split('\t')[0]] = [j]
                        else:
                            variant_vir[header_split[0].split('\t')[0]].append(j)
                            
print (len(virulents))
print (len(variant_vir))
output = []
for k, v in variant_vir.items():
    #print (k)
    if k in virulents:
        for j in v:
            output.append(j + virulents[k][0])

header = get_header(variant_file) + get_header(virulent_file)
print (header)
outfile = "{0}_virulent_genes.txt".format(variant_file.rstrip('.txt'))
with open(outfile, 'w') as outf:
    outf.write('\t'.join(header) + '\n')
    outf.writelines('\t'.join(i) + '\n' for i in output)
