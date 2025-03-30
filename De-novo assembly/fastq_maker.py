from glob import glob
import os, shutil

fastas = glob('**/*.fasta', recursive = True)
for fasta in fastas:
    gene = fasta.split('\\')[1]
    path = '/'.join(fasta.split('\\')[:-1]) + '/'
    cmd = r'.\ART\art_illumina -ss HS25 -i fasta -l 150 -f 30 -qs 36 -qs2 39 -o output'.replace('fasta', fasta)
    print(gene)
    os.system(cmd)
    os.rename('output.aln', gene + '.aln')
    
    aln_file = [i for i in open(f'{gene}.aln', 'r').read().split('\n')[4:] if i != '']
    read_ids = aln_file[::3]
    real_sequences = aln_file[1::3]
    simulated_sequences = aln_file[2::3]
    fastq_reads = []
    reads = []
    for n, read_id in enumerate(read_ids):
        if real_sequences[n] == simulated_sequences[n]:
            read_coordinate = read_id.split('\t')[1].strip()
            read_sequence = real_sequences[n].strip()
            phred_score = 'F'*len(read_sequence.strip())
            read = f'@{read_coordinate}\n{read_sequence}\n+\n{phred_score}'
            reads.append(read_sequence)
            fastq_reads.append(read)
    
    open(gene + '.sequences.txt', 'a').write('\n'.join(reads))
    open(gene + '.fastq', 'a').write('\n'.join(fastq_reads))
    os.remove('output.fq')
    shutil.move(gene + '.sequences.txt', path + '/' + gene + '.sequences.txt')
    for suffix in ('.aln', '.fastq'):
        shutil.move(gene + suffix, path + '/' + gene + suffix)

print('Done')
