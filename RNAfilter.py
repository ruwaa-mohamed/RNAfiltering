import argparse
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from matplotlib import pyplot as plt
from fpdf import FPDF

args = None

# ----------------------------------------------------------------------
def get_args():
        """
	get_args() function defines the arguments the tha user will provide for the tool. Note that all arguments are optional but the tool will not work unless either -I or -P is used.
	used libraries: argparse
	"""
        parser = argparse.ArgumentParser(
                description="RNA filtering tool that removes adapters, trims low quality bases from both ends, and filters low qiuality read and very short reads.",
                epilog="Example:\npython RNAfilter.py [-P <path> | -I <str:ID>] --read_qual <int (20)> --base_qual <int (15)> --read_len <float (0.65)> -A <str (ACGT)>"
                )

        # Arguments
        # input Run ID or FASTQ file local path (one of them is required)
        parser.add_argument('-I', help='ID of the Run from NCBI', default='', metavar='<str:ID>')
        parser.add_argument('-P', help='Local path of the FASTQ file uncompressed', default='', metavar='<path>')

        # Quality measures needed (all optional)
        parser.add_argument('--read_qual', help='the minimum average read quality score to filter out the read.', default=20, metavar='<int>')
        parser.add_argument('--base_qual', help='the minimum base quality score to trim the base from both ends.', default=15, metavar='<int>')
        parser.add_argument('--read_len', help='the minimum accepted fraction of the trimmed read to be considered in the output.', default=0.65, metavar='<float>')
        parser.add_argument('-A', help='The adapter sequence to be trimmed from the beginning of the read.', default='', metavar='<str>')

        arguments = vars(parser.parse_args())
        return arguments
# ----------------------------------------------------------------------

def download_sample(ID):
        '''
        download_sample() function takes an input NCBI-SRA Run ID and download the sample in FASTQ format in the current working directory.
        No outpyt is expected from this function.
        used libraries: 
        '''
        pass

def graphing(GC_cont, avg_qual, state):
    '''
    Input: 2 lists (GC content per read, and average quality score per read).
    Function: 2 histograms are created, one for average quality score per read and the other is GC content of the reads.
    Output: No return. The 2 graphs are saved as PNG images in the current working directory.
    used libraries: matplotlib
    '''
    plt.rcParams.update({'font.size': 12})
    plt.figure()
    plt.hist(avg_qual)
    plt.xlabel("Average Quality Score per Read")
    plt.ylabel("Number of Reads")
    plt.title(state + " Reads Average Quality Score")
    plt.savefig(state+'_1.png')

    plt.figure()
    plt.hist(GC_cont)
    plt.xlabel("Average GC content per Read")
    plt.ylabel("Number of Reads")
    plt.title(state+ " Reads Average GC Content")
    plt.savefig(state+'_2.png')


def get_stats(reads, state):
    '''
    Input:  It takes a SeqRecord iterator object that contains all the reads as SeqRecords.
            the "state" that is given is just used for the titles of the graph that would be generated.
    Function: This function is intended to generate the statistics that will be outputed later in the report.
    Output: No return.
            this function calls the graphing() function to generate the required graphs from the generated statistics.
    used libraries: Bio.SeqUtils, statistics
    '''
    sizes = []
    GC_cont = []
    avg_qual = []
    for r in reads:
        sizes.append(len(r.seq))
        GC_cont.append(GC(r.seq))
        avg_qual.append(statistics.mean(r.letter_annotations["phred_quality"]))
    d = {}
    d['Total reads'] = len(sizes)
    d['Mean read length'] = round(statistics.mean(sizes), 2)
    d['Max. read length'] = max(sizes)
    d['Min. read length'] = min(sizes)
    d['GC content'] = round(statistics.mean(GC_cont), 2)
    graphing(GC_cont, avg_qual, state)
    return d


def quality_filter(reads, qual=20):
    '''
    Input:  It takes a SeqRecord iterator object that contains all the reads as SeqRecords.
            the qual score is an optional argument that defines the minimum accepted average read quality score of the read.
    Function: for every read (loop), if the average quality score of the read is below this cutoff point, the whole read is discarded from the output
    Output: a SeqRecord iterator object containing all the input SeqRecord objects except those that disconfirm the criteria of selection.
    used libraries: statistics
    '''
    return (r for r in reads if int(statistics.mean(r.letter_annotations["phred_quality"])) > qual)

def len_filter(reads, min_len):
    '''
    Input:  It takes a SeqRecord iterator object that contains all the reads as SeqRecords.
            the min_len argument is an optional argument that defines the minimum accepted trimming percentage of any read. It's a float number.
    Function: for every read (loop), if the read sequence length is less than the percentage specified of the initial read sequences legth, the whole read is discarded from the output
    Output: a SeqRecord iterator object containing all the input SeqRecord objects except those that disconfirm the criteria of selection.
    used libraries: --
    '''
    ##init_len=int(r.description.split('length=')[1])
    return (r for r in reads if len(r)>= (min_len*int(r.description.split('length=')[1])))

def trim_adapter(reads, adapter=""):
    '''
    Input:  It takes a SeqRecord iterator object that contains all the reads as SeqRecords.
            the adapter is an optional argument that is set to empty string as default not trim anything if no adapter sequence is specified.
    Function: cuts the adapter sequence from the beginning of the read only if the whole adapter sequence is present.
    Output: the same input SeqRecord iterator object but whenever a read has the adapter specified, the adapter is removed.
    used libraries: --
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
    Input:  a SeqRecord iterator object that contains all the reads as SeqRecords.
            the min_score argument specifies the minimum accepted score for the base at the left end.
    Function: for every SeqRecord object, cut the bases from the beginning of the reads if their quality score is below the min_score given.
    Output: the same input SeqRecord iterator object but the reads the have low quality bases at the beginning are trimmed.
    used libraries: --
    '''
    for read in reads:
        qual = read.letter_annotations["phred_quality"]
        for i in range(len(qual)):
            if qual[i]>= min_score:
                break
        yield read[i:]

def trailing(reads, min_score=15):
    '''
    Input:  a SeqRecord iterator object that contains all the reads as SeqRecords.
            the min_score argument specifies the minimum accepted score for the base at the right end.
    Function: for every SeqRecord object, cut the bases from the end of the reads if their quality score is below the min_score given.
    Output: the same input SeqRecord iterator object but the reads the have low quality bases at the end are trimmed.
    used libraries: --
    '''
    for read in reads:
        qual = read.letter_annotations["phred_quality"]
        for i in range(len(qual)-1,-1,-1):
            if qual[i]>= min_score:
                break
        yield read[:i]

def report(d1, d2):
        '''
        creates the quality control report from the statistics generated un the get_stats() and graphing() functions. no return but a pdf is created in the work directory.
        '''
        # create pdf file object
        pdf = FPDF()
        pdf = FPDF('P', 'mm', 'A4')
        pdf.add_page()
        pdf.set_font('Times', 'B', 16)
        pdf.cell(60)
        pdf.cell(75, 10, "Report of project 8", 0, 1, 'C')
        pdf.cell(90, 10, " ", 0, 2, 'C')

        # specify the column width of the table
        epw = pdf.w -2*pdf.l_margin
        col_width =epw/3

        # the header of the table
        pdf.set_font('Times', 'B', 12)
        pdf.cell(col_width, 10, 'Parameters', 1, 0, 'C')
        pdf.cell(col_width, 10, 'Raw Data', 1, 0, 'C')
        pdf.cell(col_width, 10, 'Processed Data', 1, 1, 'C')

        # the rows of the table from the two dictionaries previously created
        pdf.set_font('Times', '', 12)
        for key in d1.keys():
                pdf.cell(col_width, 10, str(key), 1, 0, 'C')
                pdf.cell(col_width, 10, str(d1[key]), 1, 0, 'C')
                pdf.cell(col_width, 10, str(d2[key]), 1, 1, 'C')
                pdf.cell(90, 10, " ", 0, 2, 'C')

        # adding the four histograms
        pdf.image('Raw_1.png', x=None, y=100, w=100, h=0, type='', link='')
        pdf.image('Final_1.png', x=110, y=100, w=100, h=0, type='', link='')
        pdf.image('Raw_2.png', x=None, y=180, w=100, h=0, type='', link='')
        pdf.image('Final_2.png', x=110, y=180, w=100, h=0, type='', link='')
        pdf.output('Quality control report.pdf', 'F')


def main():
        '''
        this is the main function that orchestrates the whole work pipeline
        '''

        # get the values from the argparse function
        fq = args['P']
        if fq == '':
                fq = download_sample(args['I'])
        adapt = args['A']
        avg_qual = args['read_qual']
        min_len = args['read_len']
        min_score = args['base_qual']

        # parse the fastq file
        reads = SeqIO.parse(fq, "fastq")
        # initial stats
        print("Statistics of the Raw reads:")
        d1 = get_stats(reads, "Raw")
        print(d1)

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
        d2 = get_stats(final_reads, "Final")
        print(d2)

        report(d1, d2)

### Testing
##fq = 'SRR2079499_1.fastq'
##adapter = "AGGCCTGTCTCCTCTGAGTGATTGAC"

#------------------------------------------------------------------------------

if __name__ == '__main__':
	args = get_args()
	main()
