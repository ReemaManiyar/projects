# Implemented the pigeonhole principle for finding the approximate
# match of given Alu sequences in the excerpt of human chromosome 1.
import bisect

def extract_genome(filename):
	genome=''
	with open (filename, 'r') as fh:
		for line in fh:
			if line[0]!= '>':
				genome=genome+line.rstrip()
	return genome
	
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def query_subseq_index(p,t,index):
	k=index.k
	offset=[]
	for j in index.query(p):
		offset.append(j)
	return offset
	

	#012012012012012012012012x
def approximate_match(p,t,n):
	#segment_length= int(round(len(p)/(n+1)))
	all_matches=set()
	index_hit=0
	i=0
	print(i)
	while i<3:
		sindex= SubseqIndex(t,8,3)
		matches= query_subseq_index(p[i:],t,sindex)
		index_hit=index_hit+len(matches)
		
		for m in matches:
			if m <i or m-i+len(p)>len(t):
				continue
				
			mismatches=0
			for j in range(len(p)):
				if not p[j]==t[m+j]:
					mismatches+=1
					if mismatches>n:
						break
						
			if mismatches <=n:
				all_matches.add(m)
		i=i+1
	return list(all_matches), index_hit
	
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
t = extract_genome('chr1.GRCh38.excerpt.fasta')
subseq_ind, index_hit = approximate_match(p,t,2)
print(subseq_ind)
print(index_hit)