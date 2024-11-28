# $1 - full fasta name
# $2 - K mer length

import pandas as pd
import sys
from os import system
import gzip
from itertools import product
from Levenshtein import distance

fasta = sys.argv[1]
L = int(sys.argv[2])
reference = fasta.split('.')[0]
sequenceFile = reference + '.sequence'

def complementary(sequence):
    #Returns the complementary DNA sequence
    return sequence.upper().replace('A', 't').replace('T', 'a').replace('C','g').replace('G', 'c').upper()
    
def lexicographically_smallest_form(sequence):
    #Because every sequence can be presented in 4 ways, e.g. ATCG = GCTA = CGAT = TAGC, this function return the lexicographically smallest form. 
    rev = sequence[::-1]
    comp = complementary(sequence)
    rev_comp = comp[::-1]
    result = sorted([sequence, rev, comp, rev_comp])[0]
    return result 


with open(fasta, 'r', encoding = 'utf-8-sig') as f:
    chromosomeName = None
    chromosomeSequence = ''
    for line in f:
        if line.startswith('>'):
            if chromosomeName != None:
                open(chromosomeName + '.' + reference, 'w', encoding = 'utf-8-sig').write(chromosomeSequence)
            chromosomeName = line.split('>')[1].split(' ')[0].strip()
            chromosomeSequence = ''
        else:
            if len(line.strip()) > 0:
                chromosomeSequence += line.strip().upper()

def generate_dna_sequences(k, mode):
    nucleotides = ['A', 'C', 'G', 'T']
    if mode == 'lexi':
        sequences = sorted(list(set([lexicographically_smallest_form(i) for i in [''.join(seq) for seq in product(nucleotides, repeat = k)]])))
    else:
        sequences = [''.join(seq) for seq in product(nucleotides, repeat = k)]
    return sequences

def CIGAR(str1, str2):
    m, n = len(str1), len(str2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize dynamic programming matrix
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    # Compute edit distance
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i - 1][j],      # Insertion
                                   dp[i][j - 1],      # Deletion
                                   dp[i - 1][j - 1])  # Substitution

    # Traceback for CIGAR string
    cigar = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and dp[i][j] == dp[i - 1][j] + 1:
            cigar.append('D')  # Deletion
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + 1:
            cigar.append('I')  # Insertion
            j -= 1
        else:
            cigar.append('M' if str1[i - 1] == str2[j - 1] else 'X')  # Match or Mismatch
            i -= 1
            j -= 1

    # Reverse and compress CIGAR string
    compressed_cigar = []
    last_char = cigar[0]
    count = 1
    for c in cigar[1:] + ['']:
        if c == last_char:
            count += 1
        else:
            compressed_cigar.append(f'{count}{last_char}')
            count = 1
            last_char = c

    return ''.join(compressed_cigar[::-1])

def scanGenome(genome, kmer):
    indices = []
    genomeSequence = open(sequenceFile, 'r', encoding = 'utf-8-sig').read()
    for i in range(len(genomeSequence) - len(kmer)):
        sequence = genomeSequence[i:i + len(kmer)]
        if sequence == kmer:
            indices.append(i)
    return indices



open(fasta.split('.')[0] + '.sequence', 'w', encoding = 'utf-8-sig').write(chromosomeSequence)

tempdir = 'temp' + '_' + reference + '_' + str(L)
system('mkdir ' + tempdir)
sequences = set()
chromosome = open(reference + '.sequence', 'r', encoding = 'utf-8-sig').read()
for n in range(len(chromosome) - L):
    sequence = chromosome[n : n + L]
    if (sequence in sequences) or ('N' in sequence) or len(sequence.strip()) != L:
        continue
    else:
        sequences.add(sequence)
sequences = '\n'.join(sorted(list(sequences))) + '\n'
open(tempdir + '/' + str(L) + '.' + reference, 'w').write(sequences)

system('cat ' + tempdir + '/' + str(L) + '.' + reference + ' | sort -u -T ' + tempdir + ' > Sequences/' + str(L) + 'L.' + reference + '.sequences')
system('rm -f ' + tempdir + '/' + str(L) + '.' + reference)
#system('gzip Sequences/' + str(L) + 'L.' + reference + '.sequences')

firstLine = True
with open('Sequences/' + str(L) + 'L.' + reference + '.sequences', 'r') as f:
    for line in f:
        if len(line) > 1:
            sequence = lexicographically_smallest_form(line.strip())
            if firstLine:
                open('Sequences/' + str(L) + 'L.' + reference + '.lexi_unsorted.sequences', 'w').write(sequence)
            else:
                open('Sequences/' + str(L) + 'L.' + reference + '.lexi_unsorted.sequences', 'a').write('\n' + sequence)
            firstLine = False

system('sort -u -T ' + tempdir + ' Sequences/' + str(L) + 'L.' + reference + '.lexi_unsorted.sequences > Sequences/' + str(L) + 'L.' + reference + '.lexi.sequences')
system('rm -f Sequences/' + str(L) + 'L.' + reference + '.lexi_unsorted.sequences')
system('rm -rf ' + tempdir)
#system('gzip Sequences/' + str(L) + 'L.' + reference + '.lexi.sequences')

df = pd.read_csv('Sequences/' + str(L) + 'L.' + reference  + '.sequences', header = None)
df.columns = ['Sequence']

df['Indices'] = df['Sequence'].apply(lambda x : scanGenome(sequenceFile, x))
df.set_index('Sequence', inplace = True)
df.to_pickle('kmerome/' + str(L) + '.indices.pkl')

df = pd.read_pickle('kmerome/' + str(L) + '.indices.pkl')
df['lexi'] = [lexicographically_smallest_form(i) for i in df.index]
df['lexi_indices'] = df.apply(lambda x : [item for sublist in df[df['lexi'] == x['lexi']]['Indices'].tolist() for item in sublist], axis = 1)
df = df.sort_values(by = ['lexi'])
df = df.drop_duplicates('lexi', keep = 'first')
del df['lexi']
del df['Indices']
df.to_pickle('kmerome_lexi/' + str(L) + '.lexi.indices.pkl')

for mode in ('lexi', 'normal'):
    indices = generate_dna_sequences(L, mode)
    if mode == 'lexi':
        kmers = pd.read_pickle('kmerome_lexi/' + str(L) + '.lexi.indices.pkl').index.tolist()
    else:
        kmers = pd.read_pickle('kmerome/' + str(L) + '.indices.pkl').index.tolist()
    #kmers = [kmer for kmer in kmers if kmer.startswith(prefix)]

    results = {}
    for i in indices:
        lowest_edit_distance = L//2
        closest_kmers = []
        for kmer in kmers:
            edit_distance = distance(i, kmer)
            if edit_distance < lowest_edit_distance:
                lowest_edit_distance = edit_distance
                closest_kmers = [kmer]
            elif edit_distance == lowest_edit_distance:
                closest_kmers.append(kmer)
            else:
                continue
        l = len(closest_kmers)
        if l > 0 and l < 4:
            results[i] = (lowest_edit_distance, tuple([(kmer, CIGAR(i, kmer)) for kmer in closest_kmers]))

    df = pd.DataFrame()
    df['Query'] = results.keys()
    df['Matches'] = [i[1] for i in results.values()]
    df['Edit distance'] = [i[0] for i in results.values()]
    df.set_index('Query', inplace = True)
    if df.empty == False:
        df.to_parquet('DBSM/' + '.'.join([str(L), 'DBSM', 'parquet']))
