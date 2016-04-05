import sys
import math
import numpy as np
import pylab as pl
# from xml.dom.minidom import getDOMImplementation
import skbio
from skbio.sequence import Sequence
from bokeh.plotting import figure, show, output_file, vplot

class ContigStats(object):

    def __init__(self, genome):
        self.genome = genome

    def change_seq_lens(self, seq_lens):
        self.seq_lens = seq_lens

    # Generates the total number of base pairs and the lengths of each contig
    def get_lens(self):
        total_bps = 0
        cs = 0
        gs = 0
        seq_lens = []
        # Find out how many total base pairs and make a list of the contig sizes
        for seq in self.genome:
            this_len = len(seq)
            cs += seq.count("G")
            gs += seq.count("C")
            total_bps += this_len
            seq_lens.append(this_len)
        # Calculate GC content
        gc = float(cs + gs)
        self.total_bps = total_bps
        self.gc_cont = (gc/total_bps) * 100
        # Create a sorted list of the individual contig sizes
        self.seq_lens = sorted(seq_lens, reverse = True)
        # print("these are your seq lens " + str(seq_lens))

    # Calculate N and L statistics and mean/median contig length
    def get_stats(self):
        total_bps = self.total_bps
        seq_lens = self.seq_lens
        counter = 0
        half_bps = 0
        # Initialize a dictionary to store the stats
        stats_dict = {}
        # Look for the contig n* (the contig length for which *% of the total
        # bases are in a contig at equal to or greater than the current contig)
        while half_bps <= (total_bps / 2):
            if half_bps <= (total_bps * 0.1):
                stats_dict["n90"] = seq_lens[counter]
                stats_dict["l90"] = counter
            if half_bps <= (total_bps * 0.2):
                stats_dict["n80"] = seq_lens[counter]
                stats_dict["l80"] = counter
            if half_bps <= (total_bps * 0.3):
                stats_dict["n70"] = seq_lens[counter]
                stats_dict["l70"] = counter
            if half_bps <= (total_bps * 0.4):
                stats_dict["n60"] = seq_lens[counter]
                stats_dict["l60"] = counter
            stats_dict["n50"] = seq_lens[counter]
            stats_dict["l50"] = counter
            half_bps += seq_lens[counter]
            counter += 1

        stats_dict["median_contig"] = np.median(seq_lens)
        stats_dict["mean_contig"] = np.mean(seq_lens)
        stats_dict["total_bps"] = total_bps
        stats_dict["gc_cont"] = self.gc_cont
        stats_dict["num_contigs"] = len(seq_lens)
        self.stats_dict = stats_dict
    # Write some of the stats to xml for ingestion into SIdora (work in
    # progress)
    def write_stats_to_xml(outfilename, stats_list):
        # this_xml = getDOMImplementation()
        # stats_doc = this_xml.createDocument(None, "genome_stats", None)
        outfile = open(outfilename, "w")
        outfile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"+\
                      "<genome_stats>\n"+\
                      "  <stat>\n"+\
                      "    <name>Total base pairs</name>\n"+\
                      "    <description>The number of base pairs in the assembled genome</description>\n"+\
                      "    <value>"+str(self.stat_dict["total_bps"])+"</value>\n"+\
                      "  </stat>\n"+\
                      "  <stat>\n"+\
                      "    <name>Contig n50</name>\n"+\
                      "    <description>The contig N50 statistic of the assembled genome</description>\n"+\
                      "    <value>"+str(self.stats_dict["n50"])+"</value>\n"+\
                      "  </stat>\n"+\
                      "  <stat>\n"+\
                      "    <name>GC content</name>\n"+\
                      "    <description>The percentage of GC in the assembled genome</description>\n"+\
                      "    <value>"+str(self.stats_dict["gc_cont"])+"</value>\n"+\
                      "  </stat>\n"+\
                      "</genome_stats>")
        outfile.close()
    # Print the stats to the screen
    def print_stats(self):
        stats_dict = self.stats_dict
        print("Total number of base pairs: " + str(stats_dict["total_bps"]))
        print("Total number of contigs: " + str(stats_dict["num_contigs"]))
        these_stats = ["90", "80", "70", "60", "50"]
        for i in these_stats:
            this_stat = "n" + i
            print("N" + i + ": " + str(stats_dict[this_stat]))
        for i in these_stats:
            this_stat = "l" + i
            print("L" + i + ": " + str(stats_dict[this_stat]))
        print("GC content: " + str(float("{0:.2f}".format(stats_dict["gc_cont"]))) + "%")
        print("Median contig size: " + str(stats_dict["median_contig"]))
        print("Mean contig size: " + str(float("{0:.2f}".format(stats_dict["mean_contig"]))))

    def write_stats(self, filename):
        stats_dict = self.stats_dict
        outfile = open(filename, "w")
        outfile.write("Total number of base pairs: " +
          str(stats_dict["total_bps"]) + "\n")
        outfile.write("Total number of contigs: " +
          str(stats_dict["num_contigs"]) + "\n")
        these_stats = ["90", "80", "70", "60", "50"]
        for i in these_stats:
            this_stat = "n" + i
            outfile.write("N" + i + ": " + str(stats_dict[this_stat]) + "\n")
        for i in these_stats:
            this_stat = "l" + i
            outfile.write("L" + i + ": " + str(stats_dict[this_stat]) + "\n")
        outfile.write("GC content: " +
          str(float("{0:.2f}".format(stats_dict["gc_cont"]))) + "%" + "\n")
        outfile.write("Median contig size: " + str(stats_dict["median_contig"])
          + "\n")
        outfile.write("Mean contig size: " +
          str(float("{0:.2f}".format(stats_dict["mean_contig"]))) + "\n")

    def create_histogram(self, num_bins=50):
        # This is a histogram in bokeh, for some reason, I can't sort out the
        # log scale, so I'm working on it

        # output_file("hist.html")
        # h1 = figure(title="Contig Length Histogram", tools="save",
        #   y_axis_type="log")
        # hist, edges = np.histogram(self.seq_lens, density=True, bins=num_bins)
        # # Take the log of the tops, if you wish...
        # # other_hist = [math.log(x) for x in hist if x > 0]
        # h1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
        #   fill_color="#036564", line_color="#033649")
        # show(vplot(h1))

        # Simple pyplot histogram with log yscale
        pl.hist(self.seq_lens, bins=num_bins)
        pl.gca().set_yscale("log")
        pl.savefig("hist.png", bbox_inches="tight")
        pl.show()



# Reads a genome fasta file into a Sequence object
def read_genome(filename):
    genome = []
    for seq in skbio.io.read(filename, format='fasta'):
        genome.append(seq)
    return genome





if __name__ == "__main__":
    infilename = sys.argv[1]
    # outfile_xml = sys.argv[2]
    genome = read_genome(infilename)
    # Create a ContigStats object
    stats = ContigStats(genome)
    # Get the lengths and stats
    stats.get_lens()
    stats.get_stats()
    # Print out the stats
    stats.print_stats()
    stats.write_stats("genome_stats.txt")
    stats.create_histogram()
