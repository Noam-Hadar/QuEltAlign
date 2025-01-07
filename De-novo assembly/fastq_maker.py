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
    with open('output.fq', 'r') as f:
        n = 1
        for line in f:
            if n == 2:
                open(gene + '.sequences.txt', 'a').write(line)
                l = len(line.strip())
            if n == 4:
                line = 'F'*l + '\n'
            if n == 4:
                n = 0
            n += 1
            open(gene + '.fastq', 'a').write(line)
    os.remove('output.fq')
    shutil.move(gene + '.sequences.txt', path + '/' + gene + '.sequences.txt')
    for suffix in ('.aln', '.fastq'):
        shutil.move(gene + suffix, path + '/' + gene + suffix)

print('Done')
