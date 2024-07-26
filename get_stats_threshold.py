#!/usr/bin/env  python3

"""
This script sumarizes the distribution of a variable within the clusters at a given threshold

By Veronica Mixao
@INSA
"""

import argparse
import textwrap
import pandas

version = "1.0.0"
last_updated = "2024-07-26"

def get_info_col_summary(fname, mx, col, thr):
    """ Get stats per column and cluster """
    
    v_lst = mx[col].values.tolist()
    v = set(v_lst)
    N = 0
    N_single_cluster = 0
    N_single_cluster_100pct = 0
    N_split_2 = 0
    N_split_3_more = 0
    clusters2check = set()
    for val in v:
        if v_lst.count(val) > 1:
            N += 1
            flt_mx = mx[mx[col] == val]
            clusters_lst = flt_mx[thr].values.tolist()
            clusters = set(clusters_lst)
            if len(clusters) == 1: # all grouped into 1 cluster
                N_single_cluster += 1
                clusters2check.add(clusters_lst[0])
            elif len(clusters) == 2:
                N_split_2 += 1
            else:
                N_split_3_more += 1
    
    for cluster in clusters2check:
        flt_mx = mx[mx[thr] == cluster]
        sts_lst = flt_mx[col].values.tolist()
        sts = set(sts_lst)
        if len(sts) == 1:
            N_single_cluster_100pct += 1
    
    print(fname + "\t" + str(col) + "\t" + str(N) + "\t" + str(N_single_cluster) + "\t" + str(N_single_cluster_100pct) + "\t" + str(N_split_2) + "\t" + str(N_split_3_more))

def get_info_col(fname, mx, col, thr):
    with open(fname + "_" + col + "_" + thr + ".tsv", "w+") as out:
        print("fname\tcluster\tcluster_size\tmost_represented\tN_most_represented_dataset\tcounter_most_represented_cluster\tn_clusters_most_represented\tpct_most_represented_cluster\tpct_most_represented", file = out)
        v_lst = mx[col].values.tolist()
        v = set(v_lst)
        N = 0
        N_single_cluster = 0
        N_multiple_clusters = 0
        st_cluster_counter = {}
        for val in v:
            if v_lst.count(val) > 0:
                N += 1
                flt_mx = mx[mx[col] == val]
                clusters_lst = flt_mx[thr].values.tolist()
                clusters = set(clusters_lst)
                if len(clusters) == 1: # all grouped into 1 cluster
                    N_single_cluster += 1
                else:
                    N_multiple_clusters += 1
                st_cluster_counter[val] = len(clusters)
        
        all_clusters_lst = mx[thr].values.tolist()
        all_clusters = set(all_clusters_lst)
        N_clusters = 0
        N_single_info = 0
        N_multiple_info = 0
        for cluster in all_clusters:
            cluster_size = all_clusters_lst.count(cluster)
            if cluster_size > 0:
                N_clusters += 1
                flt_mx = mx[mx[thr] == cluster]
                info_lst = flt_mx[col].values.tolist()
                info = set(info_lst)
                if len(info) == 1:
                    N_single_info += 1
                    most_represented = info_lst[0]
                    counter = cluster_size
                else:
                    N_multiple_info += 1
                    counter = 0
                    most_represented = ""
                    for info_val in info:
                        counter_new = info_lst.count(info_val)
                        if counter_new > counter:
                            counter = counter_new
                            most_represented = info_val
                N_most_represented = len(mx[mx[col] == most_represented].values.tolist())
                if N_most_represented > 0:
                    pct_most_represented_cluster = int(counter)/float(cluster_size)
                    pct_most_represented = int(counter)/float(N_most_represented)
                    print(str(fname) + "\t" + str(cluster) + "\t" + str(cluster_size) + "\t" + str(most_represented) + "\t" + str(N_most_represented) + "\t" + str(counter) + "\t" + str(st_cluster_counter[most_represented]) + "\t" + str(pct_most_represented_cluster) + "\t" + str(pct_most_represented), file = out)

def main():
		# argument options
    
    parser = argparse.ArgumentParser(prog="get_stats_threshold.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                             get_stats_threshold.py                          #
									#                                                                             #
									############################################################################### 
									
									This script sumarizes the cluster distribution of a variable (e.g. MLST) within 
									the clustering results obtained with ReporTree at a given threshold.
                                                                                                                                                   
									-------------------------------------------------------------------------------"""))
										
    parser.add_argument("-m", "--metadata", dest="metadata", required=True, type=str, help="[MANDATORY] Metadata with partitions file (can be obtained with ReporTree).")
    parser.add_argument("-t", "--threshold", dest="threshold", required=True, help="Column name with the threshold used for reporting (e.g., single-150x1.0).")
    parser.add_argument("-col", "--column", dest="column", required=False, type=str, help="Column name with the variable to summarize (e.g. MLST_ST).")
    parser.add_argument("-a", "--analysis", dest="analysis", required=False, type=str, help="Type of analysis (summary or detail). Default: summary")

    args = parser.parse_args()
	
    print("\n******************************\n")
    print("Running get_stats_threshold.py")
    print(version, "last updated on", last_updated)

    mx = pandas.read_table(args.metadata)
    if "/" in args.metadata:
        name_info = args.metadata.split("/")[-1]
    else:
        name_info = args.metadata
    name = name_info.split(".")[0]
    col = args.column
    
    if args.analysis == "summary":
        get_info_col_summary(name, mx, col, args.threshold)
    elif args.analysis == "detail":
        get_info_col(name, mx, col, args.threshold)

if __name__ == "__main__":
    main()