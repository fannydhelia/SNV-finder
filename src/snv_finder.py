__author__ = 'Fanny-Dhelia'

import argparse
import time
from subprocess import call, Popen, PIPE
from datetime import datetime
from math import sqrt

import scipy.stats as st



# Estimates the variance for k-mer frequencies
def get_variance(kmers, l, readlen, k):
    disp = 0
    for i in range(0, k - 1):
        disp += sum(kmers[:(i + 1)]) ** 2
        disp += sum(kmers[(i + 1):]) ** 2
    disp += (readlen - 2 * k + 2) * sum(kmers)
    disp *= l
    return disp


# Finds p-values for z-tests for three genotypes
def pvalues(frequencies, kmers, l, readlen, k):
    # Test statistic
    t = frequencies[0] - frequencies[1]
    myy = l * (readlen - k + 1) * sum(kmers)
    disp = get_variance(kmers, l, readlen, k)
    # Z values
    z1 = -abs((t - myy) / sqrt(disp))
    z2 = -abs(t / sqrt(disp))
    z3 = -abs((t + myy) / sqrt(disp))
    # P-values for the z-tests
    p1 = 2 * st.norm.cdf(z1)
    p2 = 2 * st.norm.cdf(z2)
    p3 = 2 * st.norm.cdf(z3)
    p = [p1, p2, p3]
    return p


# Get a k-mer frequency from query file
def get_kmer_frequency(query):
    line2 = query.readline().strip()
    results = line2.split("\t")
    freq = results[1]
    return freq


# Get the k-mer frequencies for a SNV
def get_frequencies(n, query):
    s2 = ""
    frequencies = [0, 0]
    # For every unique k-mer
    for i in range(n):
        count = []
        # For both (two) allele variants
        for j in range(2):
            freq = get_kmer_frequency(query)
            frequencies[j] += int(freq)
            count.append(freq)
        s2 += "\t" + "/".join(count)
    return frequencies, s2


# Identify the genotype of the SNV
def get_genotype(k, alleles, loc, frequencies, l, readlen, a):
    binary = bin(loc)[2:]
    present = [0] * (k - len(binary)) + [int(x) for x in list(binary)]
    p = pvalues(frequencies, present, l, readlen, k)
    pp = [str(round(x, 4)) for x in p]
    res = [p[0] > a, p[1] > a, p[2] > a]
    decision = ""
    if sum(res) == 1:
        if res[0]:
            decision = 2 * alleles[0] + "\t" + pp[0]
        elif res[1]:
            decision = alleles[0] + alleles[1] + "\t" + pp[1]
        else:
            decision = 2 * alleles[1] + "\t" + pp[2]
    return decision


# Identifies the SNV genotypes
def find_snvs(k, a, l, readlen, snv_file, out_file, query_file):
    with open(snv_file, "r") as snvs, open(query_file, "r") as query, open(out_file, "w") as w:
        for line in snvs:
            parts = line.split("\t")
            # rs_number chr pos ref_allele/alt_allele
            s = "\t".join(parts[:4])
            n = len(parts) - 5
            frequencies, s2 = get_frequencies(n)
            s += s2
            # If none of the k-mers were found, the genotypes cannot be identified
            if sum(frequencies) > 0:
                alleles = parts[3].split("/")
                loc = int(parts[n + 4])
                decision = get_genotype(k, alleles, loc, frequencies, l, readlen, a)
                if decision:
                    # Only writes to file if the genotype is identified
                    s += "\t" + decision + "\n"
                    w.write(s)


# Get command line arguments and create help
def parse_arguments():
    parser = argparse.ArgumentParser(description='SNP finder')
    parser.add_argument('-s', '--snp', help='the file containing SNP information', required=True)
    parser.add_argument('-q', '--sequences', help='the file containing SNP kmer sequences', required=True)
    parser.add_argument('-sl', '--snplist', help='the file of kmers of snps from glistmaker', required=True)
    parser.add_argument('-l', '--list', help='list file of kmers of reads from glistmaker')
    parser.add_argument('-f', '--fastq', help='the fastq file(s) containing reads', nargs="+", required=True)
    parser.add_argument('-a', '--alpha', help='the alpha value (significance level), default=0.1',
                        type=float, default=0.1)
    parser.add_argument('-if', '--infofile', help='the file for information (default stdout)')
    parser.add_argument('-o', '--output', help='the output file for the solution')
    parser.add_argument('-i', '--info', help='show information about the process', action='store_true')
    args = parser.parse_args()
    return args


def main():
    time_start = time.time()
    k = 25

    # Time string added to the output files
    now = datetime.now()
    time_string = now.strftime('%Y%m%d_%H%M%S')

    # Parse the command line arguments
    args = parse_arguments()

    info = args.info
    if args.infofile:
        infofile = open(args.infofile, "w")

    if args.alpha:
        alpha = args.alpha

    if args.output:
        out_file = args.output
    else:
        out_file = "found_snps_" + time_string + ".txt"

    # Get the read length and count for estimating the mean and variance of k-mer frequencies
    fastq = args.fastq
    p = Popen((['wc', '-l'] + fastq), stdout=PIPE, stderr=PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    result = result.decode()
    lines = result.split("\n")
    lastline = lines[len(lines) - 2]
    n = round(int(lastline.strip().split()[0]) / 4)
    # n = 952864636

    l = n / 3000000000

    file1 = open(fastq[0], "r")
    file1.readline()
    readlen = len(file1.readline().strip())

    if info:
        s = "lambda " + str(l) + "\nreadlen " + str(readlen) + "\nnr of reads " + str(n) + "\n"
        if args.infofile:
            infofile.write(s)
        else:
            print(s)

    # Create the binary list files for sequencing data with glistmaker
    time_list_1 = time.time()
    if args.list:
        list_small = args.list
    else:
        if info:
            s = "Creating the k-mer list from fastq file(s)\n"
            if args.infofile:
                infofile.write(s)
            else:
                print(s)
            call("/mambakodu/fanny/maka/glistmaker " + " ".join(fastq) + " -w " + str(k) + " -o " + out_file[:-4],
                 shell=True)
            list_file = out_file[:-4] + "_" + str(k) + ".list"
            s = "Creating a smaller list file using the snp kmers\n"
            if args.infofile:
                infofile.write(s)
            else:
                print(s)
            list_small = out_file[:-4] + "_small_" + str(k) + "_intrsec.list"
            call("/mambakodu/fanny/maka/glistcompare " + list_file + " " + args.snplist + " -i -r first -o " + out_file[
                                                                                                               :-4] + "_small",
                 shell=True)
    time_list_2 = time.time()

    # Find the frequencies for the unique k-mers from the list files
    time_query_1 = time.time()
    if info:
        s = "Running glistquery to find the kmer frequencies from the reads\n"
        if args.infofile:
            infofile.write(s)
        else:
            print(s)
    query_result = "query_results_" + time_string + ".txt"
    call("/mambakodu/fanny/maka/glistquery " + list_small + " -f " + args.sequences + " > " +
         query_result, shell=True)
    time_query_2 = time.time()

    # Find the SNV genotypes
    time_snps_1 = time.time()
    if info:
        s = "Finding the SNPs\n"
        if args.infofile:
            infofile.write(s)
        else:
            print(s)
    snp_file = args.snp
    find_snvs(k, alpha, l, readlen, snp_file, out_file, query_result)
    time_snps_2 = time.time()

    # Output program time information
    time_finish = time.time()
    time_total = time_finish - time_start
    time_list = time_list_2 - time_list_1
    time_query = time_query_2 - time_query_1
    time_snps = time_snps_2 - time_snps_1
    if info:
        s = "Glistmaker+glistcompare time: " + str(round(time_list / 60)) + " minutes\n"
        s += "Glistquery time: " + str(round(time_query / 60)) + "minutes\n"
        s += "SNV genotyping time: " + str(round(time_snps / 60)) + "minutes\n"
        s += "Total time: " + str(round(time_total / 60)) + "minutes\n"
        if args.infofile:
            infofile.write(s)
            infofile.close()
        else:
            print(s)


if __name__ == "__main__":
    main()
