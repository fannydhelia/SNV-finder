__author__ = 'Fanny-Dhelia'

import argparse
import time
from subprocess import call, Popen, PIPE
from datetime import datetime
import scipy.stats as st
from math import sqrt


def pvalues(frequencies, kmers, l, readlen, k):
    f1 = frequencies[0]
    f2 = frequencies[1]
    t = f1 - f2
    myy = l * (readlen - k + 1) * sum(kmers)
    disp = 0
    for i in range(0, k-1):
        disp += sum(kmers[:(i+1)])**2
        disp += sum(kmers[(i+1):])**2
    disp += (readlen - 2*k + 2) * sum(kmers)
    z1 = -abs((t-myy)/sqrt(disp))
    z2 = -abs(t/sqrt(disp))
    z3 = -abs((t+myy)/sqrt(disp))
    p1 = 2*st.norm.cdf(z1)
    p2 = 2*st.norm.cdf(z2)
    p3 = 2*st.norm.cdf(z3)
    p = [p1, p2, p3]
    return p


# Find the SNVs
def find_snps(k, a, l, readlen, snp_file, out_file, query_file):
    snps = open(snp_file, "r")
    query = open(query_file, "r")
    w = open(out_file, "w")
    w2 = open(out_file[:-4] + "_2" + out_file[-4:], "w")
    line = snps.readline()
    while line:
        parts = line.split("\t")
        # rs_number chr pos ref_allele/alt_allele1/alt_allelele2
        s = "\t".join(parts[:4])
        any_kmer = False
        s2 = s
        alleles = parts[3].split("/")
        for i in range(4, len(parts)-1):
            count = []
            frequencies = [0, 0]
            for j in range(len(alleles)):
                line2 = query.readline().strip()
                results = line2.split("\t")
                kmer = results[0]
                freq = results[1]
                frequencies[j] += int(freq)
                count.append(freq)
                s3 = "\t" + kmer + "\t" + freq
                w2.write(s2 + s3 + "\n")
                if int(freq) > 0:
                    any_kmer = True
            s += "\t" + "/".join(count)
        if any_kmer:
            number = int(parts[len(parts)-1].strip())
            binary = bin(number)[2:]
            present = [0]*(k-len(binary)) + [int(x) for x in list(binary)]
            p = pvalues(frequencies, present, l, readlen, k)
            pp = [str(round(x,4)) for x in p]
            res = [p[0]>a, p[1]>a, p[2]>a]
            if sum(res) == 1:
                if res[0]:
                    dec = 2 * alleles[0] + "\t" + pp[0]
                elif res[1]:
                    dec = alleles[0] + alleles[1] + "\t" + pp[1]
                else:
                    dec = 2 * alleles[1] + "\t" + pp[2]
            elif sum(res) == 2:
                if not res[0]:
                    dec = 2 * alleles[1] + " " + alleles[0] + alleles[1] + "\t" + pp[2] + " " + pp[1]
                elif not res[1]:
                    dec = 2 * alleles[0] + " " + 2 * alleles[2] + "\t" + pp[0] + " " + pp[2]
                else:
                    dec = 2 * alleles[0] + " " + alleles[0] + alleles[1] + "\t" + pp[0] + " " + pp[1]
            else:
                dec = "N/A\t" + " ".join(pp)
            w.write(s + "\t" + dec + "\n")
        line = snps.readline()
    snps.close()
    query.close()
    w.close()
    w2.close()


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

    # Parse the command line arguments
    args = parse_arguments()

    info = args.info
    if args.infofile:
        infofile = open(args.infofile, "w")

    now = datetime.now()
    time_string = now.strftime('%Y%m%d_%H%M%S')

    if args.alpha:
        alpha = args.alpha

    if args.output:
        out_file = args.output
    else:
        out_file = "found_snps_" + time_string + ".txt"

    # Get the read length and number for statistics
    fastq = args.fastq
    p = Popen((['wc', '-l'] + fastq), stdout=PIPE, stderr=PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    result = result.decode()
    lines = result.split("\n")
    lastline = lines[len(lines)-2]
    n = round(int(lastline.strip().split()[0])/4)
    #n = 952864636
    l = n/3000000000

    file1 = open(fastq[0], "r")
    file1.readline()
    readlen = len(file1.readline().strip())

    if info:
        s = "lambda " + str(l) + "\nreadlen " + str(readlen) + "\nnr of reads " + str(n) + "\n"
        if infofile:
            infofile.write(s)
        else:
            print(s)

    time_list_1 = time.time()
    # Find the kmer frequences from the fastq file reads
    if args.list:
        list_small = args.list
    else:
        if info:
            s = "Creating the k-mer list from fastq file(s)\n"
            if infofile:
                infofile.write(s)
            else:
                print(s)
            call("/mambakodu/fanny/maka/glistmaker " + " ".join(fastq) + " -w " + str(k) + " -o " + out_file[:-4], shell=True)
            list_file = out_file[:-4] + "_" + str(k) + ".list"
            s = "Creating a smaller list file using the snp kmers\n"
            if infofile:
                infofile.write(s)
            else:
                print(s)
            list_small = out_file[:-4] + "_small_" + str(k) + "_intrsec.list"
            call("/mambakodu/fanny/maka/glistcompare " + list_file + " " + args.snplist + " -i -r first -o " + out_file[:-4] + "_small", shell=True)
    time_list_2 = time.time()

    time_query_1 = time.time()
    if info:
        s = "Running glistquery to find the kmer frequencies from the reads\n"
        if infofile:
            infofile.write(s)
        else:
            print(s)
    query_result = "query_results_" + time_string + ".txt"
    call("/mambakodu/fanny/maka/glistquery " + list_small + " -f " + args.sequences + " > " +
         query_result, shell=True)
    time_query_2 = time.time()

    time_snps_1 = time.time()
    # Find the SNPs
    if info:
        s = "Finding the SNPs\n"
        if infofile:
            infofile.write(s)
        else:
            print(s)
    snp_file = args.snp
    find_snps(k, alpha, l, readlen, snp_file, out_file, query_result)
    time_snps_2 = time.time()

    time_finish = time.time()
    time_total = time_finish - time_start
    time_list = time_list_2 - time_list_1
    time_query = time_query_2 - time_query_1
    time_snps = time_snps_2 - time_snps_1
    if info:
        s = "Glistmaker+glistcompare time: " + str(round(time_list/60)) + " minutes\n"
        s += "Glistquery time: " + str(round(time_query/60)) + "minutes\n"
        s += "SNV genotyping time: " + str(round(time_snps/60)) + "minutes\n"
        s += "Total time: " + str(round(time_total/60)) + "minutes\n"
        if infofile:
            infofile.write(s)
        else:
            print(s)
    infofile.close()


if __name__ == "__main__":
    main()