import argparse
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from matplotlib import pyplot as plt

args = None

# ----------------------------------------------------------------------
def get_args():
	"""
	
	"""
    parser = argparse.ArgumentParser(
		description="RNA filtering tool that removes adapters, trims low quality bases from both ends, and filters low qiuality read and very short reads.",
		epilog="Example:\npython RNAfilter.py [-P <path> | -I <str:ID>] --read_qual <int (20)> --base_qual <int (15)> --read_len <float (0.65)> -A <str (ACGT)>"
    )

    # Arguments
    # input Run ID or FASTQ file local path
    parser.add_argument('-I', help='ID of the Run from NCBI', default='', metavar='<str:ID>')
    parser.add_argument('-P', help='Local path of the FASTQ file uncompressed', default='', metavar='<path>')

    # Quality measures needed
    parser.add_argument('--read_qual', help='the minimum average read quality score to filter out the read.', default=20, metavar='<int>')
    parser.add_argument('--base_qual', help='the minimum base quality score to trim the base from both ends.', default=15, metavar='<int>')
    parser.add_argument('--read_len', help='the minimum accepted fraction of the trimmed read to be considered in the output.', default=0.65, metavar='<float>')
    parser.add_argument('-A', help='The adapter sequence to be trimmed from the beginning of the read.', default='', metavar='<str>')
    
    arguments = vars(parser.parse_args())
    return arguments
# ----------------------------------------------------------------------


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


def main(adapt, avg_qual=20, min_len=0.65, min_score=15):
    '''
    Input:
    Function:
    Output: 
    '''
    fq = args['P']
    if fq == '':
        fq = download_sample(args['I'])
    
    
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

#------------------------------------------------------------------------------

if __name__ == '__main__':
	args = get_args()
	main()
