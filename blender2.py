#!/usr/bin/env python
import pysam
import multiprocessing as mp
import sys
import argparse
import warnings

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="autoBLENDER",
        description = "find Cas on- and off-targets using AutoDisco")
    parser.add_argument('-f', '--file', required=True, help='experimental BAM (required)')
    parser.add_argument('-c', '--control', help='control BAM (optional, but highly recommended)')
    parser.add_argument('-g', '--guide', required=True, help="on-target guide RNA sequence, provided 5'-3' without the PAM sequence")
    parser.add_argument('-t', '--threshold', type=int, default=3, help='Number of reads to consider for threshold (dfault 3)')
    parser.add_argument('-p', '--pams', metavar='PAM', default=["GG", "AG"], nargs='+', help='PAMs to look for (default NGG NAG)')
    parser.add_argument('-r', '--reference', required=True, help='Indexed (faidx) reference genome (fasta format). Index should be called <reference>.fai')
    parser.add_argument('-m', '--max_mismatches', type=int, default=8, help="Maximum number of mismatches to allow to the guide sequence (default 8)")
    parser.add_argument('-s', '--score_min', type=int, default=3, help="Minimum score to consider a hit (default 3)")
    parser.add_argument('-b', '--blacklist', help='Blacklist to use for filtering hits, e.g. from ENCODE (BED3 format)')
    parser.add_argument('--verbose', action='store_true', default=False, help="verbose output")
    parser.add_argument('--debug', action='store_true', default=False, help='debug output')
    args = parser.parse_args()
    return args

def check_read(read, min_MQ=25):
    if read.is_unmapped:
        return False
    if read.mapping_quality <= min_MQ:
        return False
    return True

def read_in_blacklist( read, blacklist ): # backlist format = {str(chr): (start,end)}
    try:
        b_locations = blacklist[read.reference_name]
    except KeyError:
        return False
    else:
        for location in b_locations:
            (b_start, b_end) = location
            if (read.reference_start+1 >= b_start) and (read.reference_start+1 <= b_end): # add 1 to match python 0-indexing to bedfile 1-indexing
                return True
        return False

def location_in_blacklist( chromosome, start, blacklist ): # backlist format = {str(chr): (start,end)}
    try:
        b_locations = blacklist[chromosome]
    except KeyError:
        return False
    else:
        for location in b_locations:
            (b_start, b_end) = location
            if (start+1 >= b_start) and (start+1 <= b_end): # add 1 to match python 0-indexing to bedfile 1-indexing
                return True
        return False

def combine_starts(for_starts, rev_starts, threshold):
    both_starts = {}
    for start in sorted(for_starts.keys()):
        if start == 0:
            continue
        both = for_starts.get(start, 0) + rev_starts.get(start-1, 0)
        if both >= threshold:
            both_starts[start-1] = both
            both_starts[start] = both
    return both_starts

def get_pam(chromosome, location, direction, fastaref):
    ref_pam = ""
    s = None
    e = None
    if direction == "left":
        s = location-5
    elif direction == "right":
        s = location+5
    e = s+2 # python closed end notation
    if s < 0:
        return
    ref_pam = fastaref.fetch(reference=chromosome, start=s, end=e).upper()
    if direction == "left":
        ref_pam = revcomp(ref_pam)
    return ref_pam

def get_sequence(chromosome, start, end, fastaref):
    return fastaref.fetch(reference=chromosome, start=start, end=end+1).upper()

def sum_window(for_starts, rev_starts, start, window_size):
    x = 0
    for i in range( (start-window_size-1), (start+1+1) ):
        x += rev_starts.get(i, 0)
    for i in range( start, (start + window_size+1+1)):
        x += for_starts.get(i, 0)
    return x

def n_mm(seq1, seq2): # number of mismatches between two sequences
    if len(seq1) != len(seq2):
        if verbose:
            warnings.warn("Length of sequences are not equal: " + seq1 + " " + str(len(seq1)) + " / " + seq2 + " " + str(len(seq2)), RuntimeWarning)
    mm = sum(c1!=c2 for c1,c2 in zip(seq1,seq2))
    return mm

def revcomp(seq):
    rc = seq[::-1]; # reverse slicing
    table = str.maketrans("ABCDGHMNRSTUVWXYabcdghmnrstuvwxy", "TVGHCDKNYSAABWXRtvghcdknysaabwxr")
    rc = rc.translate(table)
    return rc

def get_blacklist(blacklist_fname):
    blacklist = {}
    if blacklist_fname == None:
        return {}
    f = open(blacklist_fname)
    for line in f:
        (chromosome, start, end) = line.rstrip().split()[:3]
        if chromosome not in blacklist:
            blacklist[chromosome] = [(int(start),int(end))]
        else:
            blacklist[chromosome].append((int(start),int(end)))
    if debug:
        print(blacklist)
    f.close()
    return blacklist

if __name__ == '__main__':
    args = parse_arguments()
    edited_fname = args.file
    control_fname = args.control
    input_guide = args.guide
    pams = args.pams
    reference_fname = args.reference
    threshold = args.threshold
    verbose = args.verbose
    debug = args.debug
    max_mismatches = args.max_mismatches
    score_min = args.score_min
    blacklist_fname = args.blacklist

    if (verbose):
        print (args)
    
    blacklist = get_blacklist(blacklist_fname)

    # Note that all output needs to be 1-offset to change from python 0-indexing to 1-indexing!    
    print("Chr:Start-End\tCutsite\tDiscoscore\tCutsite Ends\tStrand\tPAM\tGuide sequence\tMismatches")

    edited_bamfile = pysam.AlignmentFile(edited_fname, "rb")
    for chromosome in edited_bamfile.references:
        if chromosome[:4] == "chrUn" or chromosome[:4] == "chrM":
            continue
        if debug:
            print("Working on " + chromosome)
        for_starts = {}
        rev_starts = {}
        both_starts = {}
        for read in edited_bamfile.fetch(contig=chromosome, multiple_iterators=True):
            goodread = check_read(read)
            if goodread:
                if read.template_length == 0: # mate is unmapped
                    if not (read.flag and read.mate_is_unmapped): # this combo only happens if the read is not the 2nd in the pair
                        for_starts[read.reference_start] = for_starts.get(read.reference_start, 0) + 1
                    else:
                        rev_starts[read.reference_end-1] = rev_starts.get(read.reference_end-1, 0) + 1
                elif read.template_length > 0: # first in pair
                    for_starts[read.reference_start] = for_starts.get(read.reference_start, 0) + 1
                elif read.template_length < 0: # second in pair
                    rev_starts[read.reference_end-1] = rev_starts.get(read.reference_end-1, 0) + 1
        both_starts = combine_starts(for_starts, rev_starts, threshold)
        if debug:
            print(for_starts)
            print(rev_starts)
            print(both_starts)
        if debug:
            print(chromosome + " testing " + str(len(both_starts.keys())) + " starts")

        edited_count = {}
        if control_fname:
            ctrl_count = {}
            control_bamfile = pysam.AlignmentFile(control_fname, "rb")
        for start in both_starts.keys():
            edited_count[start] = edited_bamfile.count(contig=chromosome, start=start, stop=start)
            if control_fname:
                ctrl_count[start] = control_bamfile.count(contig=chromosome, start=start, stop=start)
        if control_fname:
            control_bamfile.close()

        reference_fasta = pysam.Fastafile(filename=reference_fname)
        
        output = {}
        for start in both_starts.keys():
            if blacklist != {}:
                if location_in_blacklist(chromosome, start, blacklist ):
                    if verbose: 
                        print (read.reference_name + ":" + str(read.reference_start) + "-" + str(read.reference_end-1) + "\t read FILTERED:blacklisted")
                    continue
            if control_fname:
                if ctrl_count[start] > 10:
                    if verbose:
                        print("CONTROL skipping " + chromosome + ":" + str(start+1) + " " + str(ctrl_count[start]))
                    continue

            if edited_count[start] > 0:
                if(both_starts[start]/edited_count[start] < 0.25):
                    if verbose: 
                        print(chromosome + ":x-x\t" +
                                str(start+1) + "\t" +
                                str(score) + "\t" + 
                                str(both_starts[start]) + "\t" +
                                "\tFILTERED: deep area")
                    continue

            score = sum_window(for_starts, rev_starts, start, window_size=5)
            if score < score_min: # doesn't pass discover-score cutoff
                if verbose: 
                        print(chromosome + ":" + "x" + "-" + "x" + "\t" +
                            str(start+1) + "\t" +
                            str(score) + "\t" + 
                            str(both_starts[start]) + "\t" + 
                            "\tFILTERED:fails disco score " + str(score_min))
                continue

            pamleft = get_pam(chromosome, start, "left", reference_fasta)
            pamright = get_pam(chromosome, start, "right", reference_fasta)
            if debug:
                print(chromosome, start, both_starts[start], pamleft, pamright, pams)

            if pamleft in pams:
                s = start - 2
                e = start + 17
                if e < 0 or s < 0: # too close to and end to be out sequence
                    continue
                guide = get_sequence(chromosome, s, e, reference_fasta)
                guide = revcomp(guide)
                mm = n_mm(input_guide, guide)
                if debug:
                    print(input_guide, guide, mm)
                if mm > max_mismatches:
                    if verbose: 
                        print(chromosome + ":" + str(s+1) + "-" + str(e+1) + "\t" +
                            str(start+1) + "\t" +
                            str(score) + "\t" + 
                            str(both_starts[start]) + "\t" + 
                            "antisense\t" +
                            "N"+pamleft+"\t" +
                            guide + 
                            " FILTERED: " + str(mm) +  " mismatches")
                else:
                    outstr = chromosome + ":" + str(s+1) + "-" + str(e+1) + "\t" + str(start+1) + "\t" + str(score) + "\t" +  str(both_starts[start]) + "\t" + "antisense\t" + "N"+pamleft+"\t" + guide + "\t" + str(mm)
                    output[chromosome+str(start)+guide] = outstr

            if pamright in pams:
                e = start + 3
                s = start - 16
                if e < 0 or s < 0: # cannot be our sequence
                    continue
                guide = get_sequence(chromosome, s, e, reference_fasta)
                mm = n_mm(input_guide, guide)
                if debug:
                    print(input_guide, guide, mm)
                if mm > max_mismatches:
                    if verbose: 
                        print(chromosome + ":" + str(s+1) + "-" + str(e+1) + "\t" +
                        str(start+1) + "\t" +
                        str(score) + "\t" + 
                        str(both_starts[start]) + "\t" + 
                        "sense\t" +
                        "N"+pamright+"\t" +
                        guide + 
                        " FILTERED: " + str(mm) +  " mismatches")
                else:
                    outstr = chromosome + ":" + str(s+1) + "-" + str(e+1) + "\t" + str(start+1) + "\t" + str(score) + "\t" + str(both_starts[start]) + "\t" + "sense\t" + "N"+pamright+"\t" + guide + "\t" + str(mm)
                    output[chromosome+str(start)+guide] = outstr
        for site in sorted(output.keys()):
            print(output[site])