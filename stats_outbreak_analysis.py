#!/usr/bin/env  python3

"""
This script was designed to assess the percentage of clusters detected by one pipeline that are also detected, 
with the exact same composition, by another pipeline.

By Veronica Mixao
@INSA
"""

import sys
import argparse
import textwrap
import pandas

version = "1.0.0"
last_updated = "2024-07-26"

def main():
    
	# argument options	----------
    
    parser = argparse.ArgumentParser(prog="stats_outbreak_analysis.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                          stats_outbreak_analysis.py                         #
									#                                                                             #
									############################################################################### 
									                            
									This script was designed to assess the percentage of clusters detected by one 
									pipeline that are also detected, with the exact same composition, by another 
									pipeline.
									                  
									-------------------------------------------------------------------------------"""))
	
    ## parameters
	
    parser.add_argument("-i", "--input", dest="infile", required=True, type=str, help="Input *tsv file containing the path to each *clusterComposition.tsv file that has to be compared (these files \
                        can be obtained with ReporTree)")
    parser.add_argument("-o", "--out", dest="output", required=True, type=str, help="Tag for output file")
    parser.add_argument("-c", "--comparison", dest="comparison", required=True, type=str, help="Type of comparison ('equal' for the assessment of whether a cluster is detected at a given threshold by \
                        another pipeline; 'lower_equal' for the assessment of wether a cluster is detected at up to a given threshold by another pipeline)")
    parser.add_argument("-t1", "--threshold1", dest="threshold1", required=True, type=int, help="Threshold at which the genetic clusters must be identified for the pipeline of interest")
    parser.add_argument("-t2", "--threshold2", dest="threshold2", required=True, type=int, help="Threshold at which the genetic clusters must be searched in the other pipelines")
	
    args = parser.parse_args()
	
    print("\n******************************\n")
    print("Running stats_outbreak_analysis.py")
    print(version, "last updated on", last_updated)

    clusters = {}
    clusters7 = {}
    samples = {}

    files = []
    thr_min = args.threshold1
    thr_max = args.threshold2
    comp_type = args.comparison # options: lower_equal or equal
    input_file = args.infile
    tag_output = args.output

    with open(input_file) as infile:
        lines = infile.readlines()
        for line in lines:
            l = line.split("\n")[0]
            files.append(l)

    for filename in files:
        if "/" in filename:
            fname = filename.split("/")[-1]
            pipeline_name = fname.split("_clusterComposition")[0]
        else:
            pipeline_name = filename.split("_clusterComposition")[0]
        
        info = {}
        all_samples = set()
        with open(filename) as cluster_composition_file:
            lines = cluster_composition_file.readlines()

            # determine all samples used in the analysis
            line_interest = lines[-1]
            lin = line_interest.split("\n")[0]
            l = lin.split("\t")
            cluster = l[3]
            samples_pipeline = set()
            for sample in cluster.split(","):
                samples_pipeline.add(sample)
            samples[pipeline_name] = samples_pipeline

            # get clusters at 7
            for line in lines[1:]:
                lin = line.split("\n")[0]
                l = lin.split("\t")
                partition_name = l[0]
                cluster = l[3]
                cluster_len = int(l[2])
                if "single" in partition_name:
                    partition = partition_name.split("single-")[1]
                    partition = partition.split("x")[0]
                elif "MST" in partition_name:
                    partition = partition_name.split("MST-")[1]
                    partition = partition.split("x")[0]
                else:
                    print("Unknown method in partition name of pipeline " + str(pipeline_name))
                    sys.exit(1)

                if comp_type == "lower_equal":    
                    if int(partition) <= thr_max and cluster_len >= 2:
                        if pipeline_name not in clusters.keys():
                            clusters[pipeline_name] = set()
                        clusters[pipeline_name].add(cluster)
                    if int(partition) == thr_min and cluster_len >= 2:
                        if pipeline_name not in clusters7.keys():
                            clusters7[pipeline_name] = set()
                        clusters7[pipeline_name].add(cluster)
                elif comp_type == "equal":
                    if int(partition) == thr_max and cluster_len >= 2:
                        if pipeline_name not in clusters.keys():
                            clusters[pipeline_name] = set()
                        clusters[pipeline_name].add(cluster)
                    if int(partition) == thr_min and cluster_len >= 2:
                        if pipeline_name not in clusters7.keys():
                            clusters7[pipeline_name] = set()
                        clusters7[pipeline_name].add(cluster)
                else:
                    print("I do not know what to do... no valid comparison type provided!")
                    sys.exit()

    # pairwise comparisons
                            
    pairwise_comparison = {"pipeline": []}
    pairwise_comparison_pct = {"pipeline": []}
    n_pipelines = len(clusters.keys())
    pipeline_info = {}
    col_order = ["pipeline"]
    cluster_info = {}

    with open(tag_output + "_missing_clusters_" + str(thr_min) + "_" + str(comp_type) + "_" + str(thr_max) + ".tsv", "w+") as out_fail:
        for pipeline1 in clusters7.keys():
            cluster_info[pipeline1] = {}
            pairwise_comparison["pipeline"].append(pipeline1)
            pairwise_comparison_pct["pipeline"].append(pipeline1)
            col_order.append(pipeline1)
            clusters_pipeline1 = clusters7[pipeline1]
            total_clusters = len(clusters_pipeline1)
            shared_at_least_one = set()

            for pipeline2 in clusters.keys():
                shared_by_two = 0
                if pipeline2 not in pairwise_comparison.keys():
                    pairwise_comparison[pipeline2] = []
                    pairwise_comparison_pct[pipeline2] = []
                clusters_pipeline2 = clusters[pipeline2]

                for cluster1 in clusters_pipeline1:
                    found = False
                    if cluster1 not in cluster_info[pipeline1].keys():
                        cluster_info[pipeline1][cluster1] = 0
                    cluster1_composition = set(cluster1.split(",")) & samples[pipeline2]
                    for cluster2 in clusters_pipeline2:
                        cluster2_composition = set(cluster2.split(",")) & samples[pipeline1]
                        if len(cluster1_composition - cluster2_composition) == 0 and len(cluster2_composition - cluster1_composition) == 0: # same cluster
                            if len(cluster1_composition) > 1:
                                shared_by_two += 1
                                cluster_info[pipeline1][cluster1] += 1
                                found = True
                                continue
                    if not found:
                        print(pipeline1, pipeline2, cluster1, file = out_fail)
                            
                pairwise_comparison[pipeline2].append(shared_by_two)
                shared_by_two_pct = float(int(shared_by_two) / int(total_clusters))
                pairwise_comparison_pct[pipeline2].append(shared_by_two_pct)
            pipeline_info[pipeline1] = total_clusters

    pairwise_df = pandas.DataFrame(data = pairwise_comparison, columns = col_order)
    pairwise_pct_df = pandas.DataFrame(data = pairwise_comparison_pct, columns = col_order)
    pairwise_df.to_csv(tag_output + "_pairwise_comparison_" + str(thr_min) + "_" + str(comp_type) + "_" + str(thr_max) + ".tsv", index = False, header=True, sep ="\t")
    pairwise_pct_df.to_csv(tag_output + "_pairwise_comparison_" + str(thr_min) + "_" + str(comp_type) + "_" + str(thr_max) + "_pct.tsv", index = False, header=True, sep ="\t")

    with open(tag_output + "_summary_" + str(thr_min) + "_" + str(comp_type) + "_" + str(thr_max) + ".tsv", "w+") as out:
        header = ["pipeline\ttotal_clusters\tn\texact_composition_all\tshared\texclusive"]
        print("\t".join(header), file = out)
        for pipeline in pipeline_info.keys():
            total_clusters = pipeline_info[pipeline]
            n = 0
            exact = 0
            shared = 0
            exclusive = 0
            for cluster in cluster_info[pipeline].keys():
                n += 1
                n_pipelines_cluster = cluster_info[pipeline][cluster]
                if n_pipelines_cluster == len(pipeline_info.keys()): # all
                    exact += 1
                elif n_pipelines_cluster < len(pipeline_info.keys()) and n_pipelines_cluster > 1: # shared by at least one
                    shared += 1
                else:
                    exclusive += 1
            info = [pipeline, str(total_clusters), str(n), str(exact), str(shared), str(exclusive)]
            print("\t".join(info), file = out)

if __name__ == "__main__":
    main()