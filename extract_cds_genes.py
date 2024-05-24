from Bio import SeqIO
import argparse 
parser = argparse.ArgumentParser(

)
parser.add_argument('path', metavar='F', help='Path to a file in genbank format.')
parser.add_argument('name', metavar='N', help='Name of file to write.')
parser.add_argument('-D', '--DNA', help='Extract DNA sequences (stranded) from coding sequences.', action='store_true')
parser.add_argument('-P', '--protein', help='(Default) Extract the \'translation\' tag (protein sequence) from coding sequences.', action='store_false')

parser.print_help()

def main(args):
    gen_record = SeqIO.parse(args.path, 'genbank')
    fasta = ''
    for rec in gen_record:
        for feat in rec.features:
            if feat.type == 'CDS':
                if args.DNA:
                    if feat.strand == -1:
                        sequence = rec.seq[feat.location.start:feat.location.end].reverse_complement()
                    else:
                        sequence = rec.seq[feat.location.start:feat.location.end]
                    format = {
                        'sequence': str(sequence),
                        'chrom': rec.id,
                        'start': feat.location.start,
                        'end': feat.location.end,
                        'strand': feat.strand,
                        'locus_tag': 'No Locus Tag!'
                        }
                    if 'locus_tag' in feat.qualifiers:
                        format['locus_tag'] = feat.qualifiers['locus_tag'][0]

                    fasta += '>{locus_tag} locus:{chrom} {start}:{end} strand: {strand}\n{sequence}\n'.format(**format)
                elif args.protein:
                    if 'translation' not in feat.qualifiers:
                        print('CDS missing a translation, omitted:')
                        print(feat.qualifiers)
                    else:
                        format = {
                            'sequence': feat.qualifiers['translation'][0],
                            'locus_tag': feat.qualifiers['locus_tag'][0],
                            'transl_table': 'No Translation Table!',
                            'protein_id': '',
                            'product': 'No Product Defined!',
                        }
                        if 'protein_id' in feat.qualifiers:
                            format['protein_id'] = feat.qualifiers['protein_id'][0] + ' '
                        if 'transl_table' in feat.qualifiers:
                            format['transl_table'] = feat.qualifiers['transl_table'][0]
                        if 'product' in feat.qualifiers:
                            format['product'] = feat.qualifiers['product'][0]

                        fasta += '>{protein_id}locus:{locus_tag} trans:{transl_table} {product}\n{sequence}\n'.format(**format)

    f = open(args.name, 'a')
    f.write(fasta)
    f.close()
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
