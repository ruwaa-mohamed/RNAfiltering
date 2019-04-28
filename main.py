import statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

def get_stats(reads):
    sizes = [len(r.seq) for r in reads]
    GC_cont = [GC(r.seq) for r in reads]
    print("Total reads: %i" % len(sizes))
    print("Mean read length: %i" % statistics.mean(sizes))
    print("Max. read length: %i" % max(sizes))
    print("Min. read length: %i" % min(sizes))
##    print(GC_cont)
##    print("GC content: %i" % statistics.mean(GC_cont) )
    print()

def graphing():
    pass

def quality_filter(reads, qual=20):
    return (r for r in reads if int(statistics.mean(r.letter_annotations["phred_quality"])) > qual)

def len_filter(reads):
    pass

def trim_adap(reads, adapter=""):
    for read in reads:
        if read.seq.find(adapter)==0:
            yield read[len(adapter):]
        else:
            yield read

def leading(reads):
    pass

def trailing(reads):
    pass


def main(fq, adapt, avg_qual=20, min_length=15):
    # parse the fastq file
    reads = SeqIO.parse(fq, "fastq")
    # initial stats
    get_stats(reads)
    # graphing?

    # parse the file again
    reads = SeqIO.parse(fq, "fastq")
    # adapter removal
    no_adapt = trim_adap(reads, adapt)
    # leading and trailing
    trim_lead = leading(no_adapt)
    trim_trail = trailing(trim_lead)
    # average quality filtering
    avg_qual_filter = quality_filter(trim_trail, avg_qual)
    # min length filtering
    min_len_filter = len_filter(avg_qual_filter, min_length)

    # save to the output
    new_fq = "filtered_" + fq
    SeqIO.write(min_len_filter, new_fq , 'fastq')
    # final stats
    final_reads = SeqIO.parse(new_fq, "fastq")
    get_stats(final_reads)
    # graphing?

###Testing
fq = 'SRR2079499_1.fastq'
adapt = "AGGCCTGTCTCCTCTGAGTGATTGAC"
main(fq, adapt)

