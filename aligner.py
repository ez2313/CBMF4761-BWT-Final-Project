'''Alignment of given read set to given reference, using BWT to align in Protein space (with exact matching)
and then continue alignment in Nucleotide space.
Output is given as CSV file'''
import pandas as pd
import numpy as np
import sys
from csv import reader



'''Credit: Benjamin Langmead "FMIndex"
next cells are from open source notebook'''
def suffixArray(s):
    ''' Given T return suffix array SA(T).  Uses "sorted"
        function for simplicity, which is probably very slow. '''
    satups = sorted([(s[i:], i) for i in range(len(s))])
    return list(map(lambda x: x[1], satups)) # extract, return just offsets

def bwtFromSa(t, sa=None):
    ''' Given T, returns BWT(T) by way of the suffix array. '''
    bw = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(t)
    for si in sa:
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si-1])
    return ''.join(bw), dollarRow # return string-ized version of list bw

class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''

    def __init__(self, bw, cpIval=4):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {}        # checkpoints
        self.cpIval = cpIval # spacing between checkpoints
        tally = {}           # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Now build the checkpoints
        for i, c in enumerate(bw):
            tally[c] += 1 # up to *and including*
            if i % cpIval == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])

    def rank(self, bw, c, row):
        ''' Return # c's there are in bw up to and including row '''
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc

st = 'teststring'
#     0123456789
cps = FmCheckpoints(st)

# should get back list of integers, where elt i gives
# # times 't' appears up to and including offset i
[ cps.rank(st, 't', i) for i in range(10) ]

# likewise for 'g'
[ cps.rank(st, 'g', i) for i in range(10) ]

class FmIndex():
    ''' O(m) size FM Index, where checkpoints and suffix array samples are
        spaced O(1) elements apart.  Queries like count() and range() are
        O(n) where n is the length of the query.  Finding all k
        occurrences of a length-n query string takes O(n + k) time.

        Note: The spacings in the suffix array sample and checkpoints can
        be chosen differently to achieve different bounds. '''

    @staticmethod
    def downsampleSuffixArray(sa, n=4):
        ''' Take only the suffix-array entries for every nth suffix.  Keep
            suffixes at offsets 0, n, 2n, etc with respect to the text.
            Return map from the rows to their suffix-array values. '''
        ssa = {}
        for i, suf in enumerate(sa):
            # We could use i % n instead of sa[i] % n, but we lose the
            # constant-time guarantee for resolutions
            if suf % n == 0:
                ssa[i] = suf
        return ssa

    def __init__(self, t, cpIval=4, ssaIval=4):
        if t[-1] != '$':
            t += '$' # add dollar if not there already
        # Get BWT string and offset of $ within it
        sa = suffixArray(t)
        self.bwt, self.dollarRow = bwtFromSa(t, sa)
        # Get downsampled suffix array, taking every 1 out of 'ssaIval'
        # elements w/r/t T
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        self.slen = len(self.bwt)
        # Make rank checkpoints
        self.cps = FmCheckpoints(self.bwt, cpIval)
        # Calculate # occurrences of each character
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        # Calculate concise representation of first column
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count

    def count(self, c):
        ''' Return number of occurrences of characters < c '''
        if c not in self.first:
            # (Unusual) case where c does not occur in text
            for cc in sorted(self.first.keys()):
                if c < cc: return self.first[cc]
            return self.first[cc]
        else:
            return self.first[c]

    def range(self, p):
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.slen - 1 # closed (inclusive) interval
        for i in range(len(p)-1, -1, -1): # from right to left
            l = self.cps.rank(self.bwt, p[i], l-1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r)   + self.count(p[i]) - 1
            if r < l:
                break
        return l, r+1

    def resolve(self, row):
        ''' Given BWM row, return its offset w/r/t T '''
        def stepLeft(row):
            ''' Step left according to character in given BWT row '''
            c = self.bwt[row]
            return self.cps.rank(self.bwt, c, row-1) + self.count(c)
        nsteps = 0
        while row not in self.ssa:
            row = stepLeft(row)
            nsteps += 1
        return self.ssa[row] + nsteps

    def hasSubstring(self, p):
        ''' Return true if and only if p is substring of indexed text '''
        l, r = self.range(p)
        return r > l

    def hasSuffix(self, p):
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.resolve(l)
        return r > l and off + len(p) == self.slen-1

    def occurrences(self, p):
        ''' Return offsets for all occurrences of p, in no particular order '''
        l, r = self.range(p)
        return [ self.resolve(x) for x in range(l, r) ]
'''end credit'''

'''function to convert nucleotide sequence in amino acidsd,and to store information beyond codons (naively only accounts for offset of mod 3):'''
def convert2(nucleotide_sequence):
    #dictionary to convert nucleotide
    mapping={
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "CGT": "R",
    "CGC" : "R",
    "CGA" : "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "AAT": "N",
    "AAC": "N",
    "GAT": "D",
    "GAC": "D",
    "AAT": "B",
    "AAC": "B",
    "GAT": "B",
    "GAC": "B",
    "TGT": "C",
    "TGC": "C",
    "CAA":"Q",
    "CAG": "Q",
    "GAA": "E",
    "GAG": "E",
    "CAA": "Z",
    "CAG": "Z",
    "GAA": "Z",
    "GAG": "Z",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "CAT": "H",
    "CAC": "H",
    "ATG": "M",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
     "TTA": "L",
     "TTG" : "L",
     "AAA": "K",
    "AAG": "K",
    "TTT": "F",
    "TTC": "F",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG" : "T",
    "TGG": "W",
    "TAT": "Y",
    "TAC": "Y",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TAA": "U",
    "TGA": "U",
    "TAG": "U",#stop codons


    }
    protein=""
    counter=0
    while counter<=len(nucleotide_sequence)-3:
      protein=protein+mapping[nucleotide_sequence[counter:counter+3]]
      counter=counter+3;
    return protein,nucleotide_sequence[counter::]


'''function to convert given nucleotide sequence into sequence of amino acids and create ORFS'''
def translate_reads2(nucleotide_sequence):
  orfs=[]
  nucs=[]
  #given read, check 6 ORF and create protein sequence for each
  for i in range(3):
    protein,extra=convert2(nucleotide_sequence[i:])
    orfs.append(protein)
    nucs.append(extra)
  #complement direction
  reversed=nucleotide_sequence[::-1]
  for i in range(3):
    protein=convert2(reversed[i:])
    orfs.append(protein)
    nucs.append(extra)
  return(orfs,nucs)



'''alignment function
#FM index in protein space and nucleotide space
#attemps exact match of read in protein space, then performs leftover inexact matching in nucleotide space'''
def align_BWT_protein_nuc(nuc_read,nucleotide_sequence):
  #create 6 ORF of read
  reads,extras=translate_reads2(nuc_read)
  #convert ref sequence to protein space
  #create FM index
  protein,temp=convert2(nucleotide_sequence)
  ref_protein=FmIndex(protein)
  number_matches=0
  matched=[]
  for read in reads:
    matches = sorted(ref_protein.occurrences(read))
    #track all reads that matched
    if ref_protein.hasSubstring(read):
      matched.append(read)

    #choose read as best based on existence+number of matches
    #best match location is last
    if len(matches)!=0 and len(matches)>number_matches:
      number_matches=len(matches)
      best_read=read
      read_number=reads.index(best_read)
      best_matches=matches

      #if the best matched read has leftover nucleotide information, use to align further
      if extras[read_number]!="":
        newread=extras[read_number]
        newref=nucleotide_sequence[:best_matches[-1]*3]
        i=0
        j=0
        while i <len(newread) and j<len(newref):
            if newread[i]==newref[j]:
              number_matches=number_matches+1
              j=j+1
              i=i+1
            else:
              i=i+1


  return number_matches,best_matches


'''run alignment on given inputs, output to given location'''
def main():
    #check and set output file
    if len(sys.argv)>3:
        output=sys.argv[3]
    else:
        output='output.txt'
    #read in reference
    f=open(sys.argv[2],'r')
    reference=f.read()
    f.close()
    alignments=[]

    #input in read dataset
    with open(sys.argv[1],'r') as dataset:
        csv_reader=reader(dataset)
        for row in csv_reader:
            read=str(row)[2:-2]
            score, matches=align_BWT_protein_nuc(read,reference)
            alignments.append((score,matches))

    #output results to given  location
    results_df=pd.DataFrame(alignments)
    results_df.to_csv(output,header=None,index=None)




if __name__ == "__main__":

    main()
