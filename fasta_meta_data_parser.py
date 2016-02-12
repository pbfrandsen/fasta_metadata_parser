import sys
# from xml.dom.minidom import getDOMImplementation
import skbio
from skbio.sequence import Sequence

# Reads a genome fasta file into a Sequence object
def read_genome(filename):
    genome = []
    for seq in skbio.io.read(filename, format='fasta'):
        genome.append(seq)
    return genome

# Generates the total number of base pairs and esimates the contig n50 for a
# genome
def get_stats(genome):
    total_bps = 0
    cs = 0
    gs = 0
    seq_lens = []
    # Find out how many total base pairs and make a list of the contig sizes
    for seq in genome:
        this_len = len(seq)
        cs += seq.count("G")
        gs += seq.count("C")
        total_bps += this_len
        seq_lens.append(this_len)

    # Calculate GC content
    gc = float(cs + gs)
    gc_cont = (gc/total_bps) * 100
    # Create a sorted list of the individual contig sizes
    seq_lens = sorted(seq_lens, reverse = True)
    counter = 0
    half_bps = 0
    # Look for the contig n50 (the contig length for which half of the total
    # bases are in a contig at equal to or greater than the current contig)
    while half_bps <= (total_bps / 2):
        half_bps += seq_lens[counter]
        counter += 1

    n50 = seq_lens[counter - 1]
    length = sum(seq_lens[0:counter])
    print "Here are the contigs before the contig n50: " + str(seq_lens[0:(counter + 1)])
    # print "There are " + str(len(seq_lens)) + " contigs!"
    # print "we are looking at contig number " + str(counter-1) + "!"
    # print "The length of the contigs before this one is: " + str(length)
    # Return the total number of base pairs and the contig n50 in a tuple
    return total_bps, n50, gc_cont

def write_stats_to_xml(outfilename, stats_list):
    # this_xml = getDOMImplementation()
    # stats_doc = this_xml.createDocument(None, "genome_stats", None)
    outfile = open(outfilename, "w")
    outfile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"+\
                  "<genome_stats>\n"+\
                  "  <stat>\n"+\
                  "    <name>Total base pairs</name>\n"+\
                  "    <description>The number of base pairs in the assembled genome</description>\n"+\
                  "    <value>"+str(stats_list[0])+"</value>\n"+\
                  "  </stat>\n"+\
                  "  <stat>\n"+\
                  "    <name>Contig n50</name>\n"+\
                  "    <description>The contig N50 statistic of the assembled genome</description>\n"+\
                  "    <value>"+str(stats_list[1])+"</value>\n"+\
                  "  </stat>\n"+\
                  "  <stat>\n"+\
                  "    <name>GC content</name>\n"+\
                  "    <description>The percentage of GC in the assembled genome</description>\n"+\
                  "    <value>"+str(stats_list[2])+"</value>\n"+\
                  "  </stat>\n"+\
                  "</genome_stats>")
    outfile.close()



if __name__ == "__main__":
    infilename = sys.argv[1]
    outfile_xml = sys.argv[2]
    genome = read_genome(infilename)
    # Grab the stats in this order: total_bps, n50, & gc_content
    stats_list = get_stats(genome)
    write_stats_to_xml(outfile_xml, stats_list)
    print "The total number of base pairs is: " + str(stats_list[0])
    print "Your contig N50 is: " + str(stats_list[1])
    print "Your GC content is: " + str(float("{0:.2f}".format(stats_list[2]))) + "%"


