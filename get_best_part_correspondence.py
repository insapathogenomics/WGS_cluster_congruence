#!/usr/bin/env	python3

"""
This script was designed to determine the partition correspondence between two pipelines.
It takes as input the indication of a directory with all the *final_score.tsv files (obtained with comparing_partitions_v2.py). 

By Veronica Mixao
@INSA
"""

import sys
import argparse
import textwrap
import glob
import pandas

version = "1.0.1"
last_updated = "2025-01-09"

def get_correspondence_all(input_dir, score):
    final = {"comparison": [], "method1": [], "method2": []}
    print("\n" + "Checking:")

    for filename in glob.glob(input_dir + "/*final_score.tsv"):
        with open(filename) as infile:
            counter = 0
            f_lines = infile.readlines()
            if "/" in filename:
                fname = filename.split("/")[-1]
            else:
                fname = filename
            comparison = fname.split("_final_score.tsv")[0]
            print("\t" + comparison)
            for line in f_lines:
                lin = line.split("\n")
                l = lin[0].split("\t")
                if counter > 0:
                    max_value = max(l[1:])
                    pos = [i for i, j in enumerate(l[1:]) if j == max_value]
                    if len(pos) == 1:
                        method1_pos = counter - 1
                        method2_pos = pos[0]
                    else:
                        method1_pos = counter - 1
                        dists = []
                        for val in pos:
                            dist = abs(counter - int(val))
                            dists.append(dist)
                        min_dist = min(dists)
                        pos2 = [i for i, j in enumerate(dists) if j == min_dist]
                        method2_pos = pos[pos2[0]]
                    if float(max_value) >= score:
                        final["comparison"].append(comparison)
                        final["method1"].append(method1_pos)
                        final["method2"].append(method2_pos)
                    else:
                        final["comparison"].append(comparison)
                        final["method1"].append(method1_pos)
                        final["method2"].append("-")
                counter += 1    
        mx = pandas.read_table(filename)
        counter2 = 0
        if "/" in filename:
            fname = filename.split("/")[-1]
        else:
            fname = filename
        comparison = fname.split("_final_score.tsv")[0] + "_rev"
        for info in mx.columns[1:]:
            max_value = max(mx[info])
            pos = [i for i, j in enumerate(mx[info]) if j == max_value]
            if len(pos) == 1:
                method1_pos = counter2
                method2_pos	= pos[0]
            else:
                method1_pos	= counter2
                dists = []
                for val in pos:
                    dist = abs(counter2 - int(val))
                    dists.append(dist)
                min_dist = min(dists)
                pos2 = [i for i, j in enumerate(dists) if j == min_dist]
                method2_pos = pos[pos2[0]]
            if float(max_value) >=	score:
                final["comparison"].append(comparison)
                final["method1"].append(method1_pos)
                final["method2"].append(method2_pos) 
            else:
                final["comparison"].append(comparison) 
                final["method1"].append(method1_pos) 
                final["method2"].append("-")
            counter2 += 1
            
    final_results = pandas.DataFrame(data = final)
    final_results.to_csv(input_dir + "/ALL_CORRESPONDENCE.tsv", index = False, header=True, sep = "\t")

def main():
	# argument options
    
    parser = argparse.ArgumentParser(prog="get_best_part_correspondence_threshold.py ", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                  get_best_part_correspondence_threshold.py                  #
									#                                                                             #
									############################################################################### 
									
									This script was designed to determine the partition correspondence between two 
									pipelines. It takes as input the indication of a directory with all the 
									*final_score.tsv files (obtained with comparing_partitions_vs.py) for which 
									the partition correspondence must be determined. In the end, it generates the
									ALL_CORRESPONDENCE.tsv file with the best partition correspondence for each
									pairwise pipeline comparison. A threshold of 2.85 is the default minimum score
									to determine this correspondence. It can be modified with "--score".
                                                                                                                                                   
									-------------------------------------------------------------------------------"""))
										
    parser.add_argument("-i", "--inputdir", dest="inputdir", required=True, type=str, help="[MANDATORY] Input directory with all the *final_score.tsv files from the all-to-all partition comparison between two pipelines (obtained with comparing_partitions_v2.py).")
    parser.add_argument("-s", "--score", dest="score", required=False, default=2.85, type=float, help="[OPTIONAL] Minimum score to consider two partitions as a correspondence. Default: 2.85.")

    args = parser.parse_args()
	
    print("\n******************************\n")
    print("Running get_best_part_correspondence_threshold.py")
    print(version, "last updated on", last_updated)

    get_correspondence_all(args.inputdir, args.score)

if __name__ == "__main__":
    main()
