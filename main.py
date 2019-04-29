import statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from matplotlib import pyplot as plt

def graphing(GC_cont, avg_qual, state):
    '''
    Input: 2 lists (GC content per read, and average quality score per read).
    Function: 2 histograms are created, one for average quality score per read and the other is GC content of the reads.
    Output: No return.
    '''
    plt.rcParams.update({'font.size': 24})
    plt.figure(1)
    plt.hist(avg_qual)
    plt.xlabel("Average Quality Score per Read")
    plt.ylabel("Number of Reads")
    plt.title(state + " Reads Averag Quality Score")

    plt.figure(2)
    plt.hist(GC_cont)
    plt.xlabel("Average GC content per Read")
    plt.ylabel("Number of Reads")
    plt.title(state+ " Reads Averag GC Content")
    plt.show()


def get_stats(reads, state):
    '''
    Input:
    Function:
    Output: No return.
    '''
    sizes = []
    GC_cont = []
    avg_qual = []
    for r in reads:
        sizes.append(len(r.seq))
        GC_cont.append(GC(r.seq))
        avg_qual.append(statistics.mean(r.letter_annotations["phred_quality"]))
    print("Total reads: %i" % len(sizes))
    print("Mean read length: %i" % statistics.mean(sizes))
    print("Max. read length: %i" % max(sizes))
    print("Min. read length: %i" % min(sizes))
    print("GC content: %i" % statistics.mean(GC_cont))
    print()
##    graphing(GC_cont, avg_qual, state)


def quality_filter(reads, qual=20):
    '''
    Input:
    Function:
    Output: 
    '''
    return (r for r in reads if int(statistics.mean(r.letter_annotations["phred_quality"])) > qual)

def len_filter(reads, min_len):
    '''
    Input:
    Function:
    Output: 
    '''
    ##init_len=int(r.description.split('length=')[1])
    return (r for r in reads if len(r)>= (min_len*int(r.description.split('length=')[1])))

def trim_adapter(reads, adapter=""):
    '''
    Input:
    Function:
    Output: 
    '''
    if adapter=="":
        return reads
    for read in reads:
        if read.seq.find(adapter)==0:
            yield read[len(adapter):]
        else:
            yield read

def leading(reads, min_score=15):
    '''
    Input:
    Function:
    Output: 
    '''
    for read in reads:
        qual = read.letter_annotations["phred_quality"]
        for i in range(len(qual)):
            if qual[i]>= min_score:
                break
        yield read[i:]

def trailing(reads):
    '''
    Input:
    Function:
    Output: 
    '''
    for read in reads:
        qual = read.letter_annotations["phred_quality"]
        for i in range(len(qual)-1,-1,-1):
            if qual[i]>= min_score:
                break
        yield read[:i]


def main(fq, adapt, avg_qual=20, min_len=0.65, min_score=15):
    '''
    Input:
    Function:
    Output: 
    '''
    # parse the fastq file
    reads = SeqIO.parse(fq, "fastq")
    # initial stats
    print("Statistics of the Raw reads:")
    get_stats(reads, "Raw")

    # parse the file again
    reads = SeqIO.parse(fq, "fastq")
    # adapter removal
    no_adapter = trim_adapter(reads, adapt)
    # leading and trailing
    trim_lead = leading(no_adapter, min_score)
    trim_trail = trailing(trim_lead, min_score)
    # average quality filtering
    avg_qual_filter = quality_filter(trim_trail, avg_qual)
    # min length filtering
    min_len_filter = len_filter(avg_qual_filter, min_len)

    # save to the output
    new_fq = "filtered_" + fq
    SeqIO.write(min_len_filter, new_fq , 'fastq')
    # final stats
    final_reads = SeqIO.parse(new_fq, "fastq")
    print("Statistics of the processed reads:")
    get_stats(final_reads, "Final")

### Testing
fq = 'SRR2079499_1.fastq'
adapter = "AGGCCTGTCTCCTCTGAGTGATTGAC"
##main(fq, adapt)
reads = SeqIO.parse(fq, "fastq")
get_stats(reads, "Raw")
reads = SeqIO.parse(fq, "fastq")
