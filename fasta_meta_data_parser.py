#! /usr/bin/env python
import sys
import math
import numpy as np
import pylab as pl
# from xml.dom.minidom import getDOMImplementation
import skbio
from skbio.sequence import Sequence
# from bokeh.plotting import figure, show, output_file, vplot

class FastaStats(object):

    def __init__(self, genome):
        self.genome = genome

    def change_scaf_eq_lens(self, seq_lens):
        self.seq_lens_scaf = seq_lens

    def get_contig_lens(self):
        total_bps = 0
        cs = 0
        gs = 0
        seq_lens = []
        # Find out how many total base pairs and make a list of the contig 
        # sizes
        for seq in self.genome:
            seq = str(seq)
            if "NN" in seq:
                contig_list = seq.split("NN")
                for contig in contig_list:
                    this_len = len(contig)
                    cs += contig.count("G")
                    gs += contig.count("C")
                    if this_len != 0:
                        total_bps += this_len
                        seq_lens.append(this_len)
            else:
                this_len = len(seq)
                cs += seq.count("G")
                gs += seq.count("C")
                total_bps += this_len
                seq_lens.append(this_len)

        # Calculate GC content
        gc = float(cs + gs)
        self.total_bps_contig = total_bps
        self.gc_cont_contig = (gc/total_bps) * 100
        # Create a sorted list of the individual contig sizes
        self.seq_lens_contig = sorted(seq_lens, reverse = True)

    # Generates the total number of base pairs and the lengths of each scaffold
    def get_scaf_lens(self):
        total_bps = 0
        cs = 0
        gs = 0
        seq_lens = []
        # Find out how many total base pairs and make a list of the scaffold 
        # sizes
        for seq in self.genome:
            seq = str(seq)
            # seq = seq.lower()
            this_len = len(seq)
            cs += seq.count("G")
            gs += seq.count("C")
            total_bps += this_len
            seq_lens.append(this_len)
        # Calculate GC content
        gc = float(cs + gs)
        self.total_bps_scaf = total_bps
        self.gc_cont_scaf = (gc/total_bps) * 100
        # Create a sorted list of the individual contig sizes
        self.seq_lens_scaf = sorted(seq_lens, reverse = True)
        # print("these are your seq lens " + str(seq_lens))

    # Calculate N and L statistics and mean/median contig length
    def get_stats_contig(self):
        total_bps = self.total_bps_contig
        seq_lens = self.seq_lens_contig
        counter = 0
        half_bps = 0
        # Initialize a dictionary to store the stats
        stats_dict = {}
        # Look for the contig n* (the contig length for which *% of the total
        # bases are in a contig at equal to or greater than the current contig)
        while half_bps <= (total_bps / 2):
            if half_bps <= (total_bps * 0.1):
                stats_dict["n10"] = seq_lens[counter]
                stats_dict["l10"] = counter
            if half_bps <= (total_bps * 0.2):
                stats_dict["n20"] = seq_lens[counter]
                stats_dict["l20"] = counter
            if half_bps <= (total_bps * 0.3):
                stats_dict["n30"] = seq_lens[counter]
                stats_dict["l30"] = counter
            if half_bps <= (total_bps * 0.4):
                stats_dict["n40"] = seq_lens[counter]
                stats_dict["l40"] = counter
            stats_dict["n50"] = seq_lens[counter]
            stats_dict["l50"] = counter
            half_bps += seq_lens[counter]
            counter += 1

        stats_dict["median_contig"] = np.median(seq_lens)
        stats_dict["mean_contig"] = np.mean(seq_lens)
        stats_dict["total_bps"] = total_bps
        stats_dict["gc_cont"] = self.gc_cont_contig
        stats_dict["num_contigs"] = len(seq_lens)
        stats_dict["largest_contig"] = seq_lens[0]
        stats_dict["shortest_contig"] = seq_lens[-1]
        self.contig_stats_dict = stats_dict
    # Write some of the stats to xml for ingestion into SIdora (work in
    # progress)

    def get_stats_scaf(self):
        total_bps = self.total_bps_scaf
        seq_lens = self.seq_lens_scaf
        counter = 0
        half_bps = 0
        # Initialize a dictionary to store the stats
        stats_dict = {}
        # Look for the contig n* (the contig length for which *% of the total
        # bases are in a contig at equal to or greater than the current contig)
        while half_bps <= (total_bps / 2):
            if half_bps <= (total_bps * 0.1):
                stats_dict["n10"] = seq_lens[counter]
                stats_dict["l10"] = counter
            if half_bps <= (total_bps * 0.2):
                stats_dict["n20"] = seq_lens[counter]
                stats_dict["l20"] = counter
            if half_bps <= (total_bps * 0.3):
                stats_dict["n30"] = seq_lens[counter]
                stats_dict["l30"] = counter
            if half_bps <= (total_bps * 0.4):
                stats_dict["n40"] = seq_lens[counter]
                stats_dict["l40"] = counter
            stats_dict["n50"] = seq_lens[counter]
            stats_dict["l50"] = counter
            half_bps += seq_lens[counter]
            counter += 1

        stats_dict["median_scaf"] = np.median(seq_lens)
        stats_dict["mean_scaf"] = np.mean(seq_lens)
        stats_dict["total_bps"] = total_bps
        stats_dict["gc_cont"] = self.gc_cont_scaf
        stats_dict["num_scaf"] = len(seq_lens)
        stats_dict["largest_scaf"] = seq_lens[0]
        stats_dict["shortest_scaf"] = seq_lens[-1]
        self.scaf_stats_dict = stats_dict

    # def write_stats_to_xml(self, outfilename):
    #     # this_xml = getDOMImplementation()
    #     # stats_doc = this_xml.createDocument(None, "genome_stats", None)
    #     stats_dict = self.stats_dict
    #     outfile = open(outfilename, "w")
    #     outfile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    #                   "<document>\n" +
    #                   "    <TotalNumberOfBasePairs>"+str(stats_dict["total_bps"])+"</TotalNumberOfBasePairs>\n" +
    #                   "    <TotalNumberOfContigs>"+str(stats_dict["num_contigs"])+"</TotalNumberOfContigs>\n" +
    #                   "    <N10>"+str(stats_dict["n10"])+"</N10>\n" +
    #                   "    <N20>"+str(stats_dict["n20"])+"</N20>\n" +
    #                   "    <N30>"+str(stats_dict["n30"])+"</N30>\n" +
    #                   "    <N40>"+str(stats_dict["n40"])+"</N40>\n" +
    #                   "    <N50>"+str(stats_dict["n50"])+"</N50>\n" +
    #                   "    <L10>"+str(stats_dict["l10"])+"</L10>\n" +
    #                   "    <L20>"+str(stats_dict["l20"])+"</L20>\n" +
    #                   "    <L30>"+str(stats_dict["l30"])+"</L30>\n" +
    #                   "    <L40>"+str(stats_dict["l40"])+"</L40>\n" +
    #                   "    <L50>"+str(stats_dict["l50"])+"</L50>\n" +
    #                   "    <GCcontent>"+str(float("{0:.2f}".format(stats_dict["gc_cont"]))) + "%"+"</GCcontent>\n" +
    #                   "    <MedianContigSize>"+str(stats_dict["median_contig"])+"</MedianContigSize>\n" +
    #                   "    <MeanContigSize>"+str(float("{0:.2f}".format(stats_dict["mean_contig"])))+"</MeanContigSize>\n" +
    #                   "    <LongestContigIs>"+str(float("{0:.2f}".format(stats_dict["largest_contig"])))+"</LongestContigIs>\n" +
    #                   "    <ShortestContigIs>"+str(float("{0:.2f}".format(stats_dict["shortest_contig"])))+"</ShortestContigIs>\n" +
    #                   "</document>\n")
    #     outfile.close()

    # Print the stats to the screen
    def print_stats_contig(self):
        stats_dict = self.contig_stats_dict
        print("Contig statistics:")
        print("Total number of base pairs: " + str(stats_dict["total_bps"]))
        print("Total number of contigs: " + str(stats_dict["num_contigs"]))
        these_stats = ["10", "20", "30", "40", "50"]
        for i in these_stats:
            this_stat = "n" + i
            print("N" + i + ": " + str(stats_dict[this_stat]))
        for i in these_stats:
            this_stat = "l" + i
            print("L" + i + ": " + str(stats_dict[this_stat]))
        print("GC content: " + str(float("{0:.2f}".format(stats_dict["gc_cont"]))) + "%")
        print("Median contig size: " + str(stats_dict["median_contig"]))
        print("Mean contig size: " + str(float("{0:.2f}".format(stats_dict["mean_contig"]))))
        print("Longest contig is: " + str(float("{0:.2f}".format(stats_dict["largest_contig"]))))
        print("Shortest contig is: " + str(float("{0:.2f}".format(stats_dict["shortest_contig"]))))

    def print_stats_scaf(self):
        stats_dict = self.scaf_stats_dict
        print("\nScaffold statistics:")
        print("Total number of base pairs: " + str(stats_dict["total_bps"]))
        print("Total number of scaffolds: " + str(stats_dict["num_scaf"]))
        these_stats = ["10", "20", "30", "40", "50"]
        for i in these_stats:
            this_stat = "n" + i
            print("N" + i + ": " + str(stats_dict[this_stat]))
        for i in these_stats:
            this_stat = "l" + i
            print("L" + i + ": " + str(stats_dict[this_stat]))
        print("GC content: " + str(float("{0:.2f}".format(stats_dict["gc_cont"]))) + "%")
        print("Median scaffold size: " + str(stats_dict["median_scaf"]))
        print("Mean scaffold size: " + str(float("{0:.2f}".format(stats_dict["mean_scaf"]))))
        print("Longest scaffold is: " + str(float("{0:.2f}".format(stats_dict["largest_scaf"]))))
        print("Shortest scaffold is: " + str(float("{0:.2f}".format(stats_dict["shortest_scaf"]))))

    def write_contig_stats(self, filename):
        stats_dict = self.contig_stats_dict
        outfile = open(filename, "w")
        outfile.write("Total number of base pairs: " +
          str(stats_dict["total_bps"]) + "\n")
        outfile.write("Total number of contigs: " +
          str(stats_dict["num_contigs"]) + "\n")
        these_stats = ["10", "20", "30", "40", "50"]
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
        outfile.write("Longest contig is: " +
          str(float("{0:.2f}".format(stats_dict["largest_contig"]))) + "\n")
        outfile.write("Shortest contig is: " +
          str(float("{0:.2f}".format(stats_dict["shortest_contig"]))) + "\n")

    def write_scaf_stats(self, filename):
        stats_dict = self.scaf_stats_dict
        outfile = open(filename, "w")
        outfile.write("Total number of base pairs: " +
          str(stats_dict["total_bps"]) + "\n")
        outfile.write("Total number of scaffolds: " +
          str(stats_dict["num_scaf"]) + "\n")
        these_stats = ["10", "20", "30", "40", "50"]
        for i in these_stats:
            this_stat = "n" + i
            outfile.write("N" + i + ": " + str(stats_dict[this_stat]) + "\n")
        for i in these_stats:
            this_stat = "l" + i
            outfile.write("L" + i + ": " + str(stats_dict[this_stat]) + "\n")
        outfile.write("GC content: " +
          str(float("{0:.2f}".format(stats_dict["gc_cont"]))) + "%" + "\n")
        outfile.write("Median scaffold size: " + str(stats_dict["median_scaf"])
          + "\n")
        outfile.write("Mean scaffold size: " +
          str(float("{0:.2f}".format(stats_dict["mean_scaf"]))) + "\n")
        outfile.write("Longest scaffold is: " +
          str(float("{0:.2f}".format(stats_dict["largest_scaf"]))) + "\n")
        outfile.write("Shortest scaffold is: " +
          str(float("{0:.2f}".format(stats_dict["shortest_scaf"]))) + "\n")

    def create_contig_histogram(self, num_bins=50):
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
        pl.hist(self.seq_lens_contig, bins=num_bins)
        pl.gca().set_yscale("log")
        pl.savefig("hist.png", bbox_inches="tight")
        pl.show()

    def create_scaf_histogram(self, num_bins=50):
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
        pl.hist(self.seq_lens_scaf, bins=num_bins)
        pl.gca().set_yscale("log")
        pl.savefig("hist.png", bbox_inches="tight")
        pl.show()


# Reads a genome fasta file into a Sequence object
def read_genome(filename):
    genome = []
    for seq in skbio.io.read(filename, format='fasta', lowercase=True):
        genome.append(seq)
    return genome





if __name__ == "__main__":
    infilename = sys.argv[1]
    # outfile_xml = sys.argv[2]
    genome = read_genome(infilename)
    # Create a ContigStats object
    stats = FastaStats(genome)
    # Get the lengths and stats
    stats.get_scaf_lens()
    stats.get_contig_lens()
    stats.get_stats_scaf()
    stats.get_stats_contig()
    # Print out the stats
    stats.print_stats_contig()
    stats.print_stats_scaf()
    # stats.write_stats_to_xml("genome_stats.xml")
    # stats.write_stats("genome_stats.txt")
    # stats.create_histogram()
