
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

def quality_filter(reads, qual=20):
##    filtered = [r for r in reads if statistics.mean(r.letter_annotations["phred_quality"]) > qual]
    return (r for r in reads if int(statistics.mean(r.letter_annotations["phred_quality"])) > qual)

def trim_polyA(records, numA, minLen):
    for record in records:
        if len(record) < minLen: continue
        record = record.seq.split("A"*numA, 1)[0]
        yield record

def _remove_adaptor(seq, region, right_side=True):
    if right_side:
        try:
            pos = seq.find(region)
        # handle Biopython SeqRecords
        except AttributeError:
            pos = seq.seq.find(region)
        return seq[:pos]
    else:
        try:
            pos = seq.rfind(region)
        # handle Biopython SeqRecords
        except AttributeError:
            pos = seq.seq.rfind(region)
        return seq[pos+len(region):]

def trim_adaptor(seq, adaptor, num_errors, right_side=True):
    gap_char = '-'
    exact_pos = str(seq).find(adaptor)
    if exact_pos >= 0:
        seq_region = str(seq[exact_pos:exact_pos+len(adaptor)])
        adapt_region = adaptor
    else:
        seq_a, adaptor_a, score, start, end = pairwise2.align.localms(str(seq),
                                                                      str(adaptor),
                                                                      5.0, -4.0, -9.0, -0.5,
                                                                      one_alignment_only=True,
                                                                      gap_char=gap_char)[0]
        adapt_region = adaptor_a[start:end]
        seq_region = seq_a[start:end]
    matches = sum((1 if s == adapt_region[i] else 0) for i, s in enumerate(seq_region))
    # too many errors -- no trimming
    if (len(adaptor) - matches) > num_errors:
        return seq
    # remove the adaptor sequence and return the result
    else:
        return _remove_adaptor(seq, seq_region.replace(gap_char, ""),
                right_side)

def process_reads(fq, qual, adapt, numA, minLen):
    with gzip.open(fq) as f:
        rawReads = SeqIO.parse(f, "fastq")
        # get_stats(rawReads) # When I run this, everything downstream fails..
        qualFil = quality_filter(rawReads, qual) # I think this work fine.
        trimmedPoly = trim_polyA(qualFil, numA, minLen)
        trimmedAdap = trim_adaptor(trimmedPoly, adapt, 2)

        # count = SeqIO.write(trimmedAdap, "good_quality.fastq", "fastq")
        # print(count)

# TEST PROCESSING
#fq = "test/TAAGGCGA_2.fq.gz"
#process_reads(fq, qual=50, adapt="AAGCAGTGGTATCAACGCAGAGTGAATGGG", numA=6, minLen=20)



# open the fastq file, read it, and do the pre- statistics
fq = 'SRR2079499_1.fastq'
reads = SeqIO.parse(fq, "fastq")
get_stats(reads)

# read the same file again, do the filtering and trimming, and write to a new file
reads = SeqIO.parse(fq, "fastq")
new_reads=quality_filter(reads)
#all other filters should be added here!
SeqIO.write(new_reads, 'new_reads.fastq', 'fastq')
# read that new output file to do the post- statistics
fq2 = 'new_reads.fastq'
reads_2 = SeqIO.parse(fq2, "fastq")
get_stats(reads_2)
