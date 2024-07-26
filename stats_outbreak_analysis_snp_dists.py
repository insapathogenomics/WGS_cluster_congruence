#!/usr/bin/env  python3

"""
This script was designed to assess the allele and SNP distances within outbreak clusters.

By Veronica Mixao
@INSA
"""

import sys
import argparse
import textwrap
import pandas
import glob
import statistics

version = "1.0.0"
last_updated = "2024-07-26"

def main():
    
	# argument options	----------
    
    parser = argparse.ArgumentParser(prog="stats_outbreak_analysis_snp_dists.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                    stats_outbreak_analysis_snp_dists.py                     #
									#                                                                             #
									############################################################################### 
									                            
									This script was designed to assess the allele and SNP distances within
									outbreak clusters.
									                  
									-------------------------------------------------------------------------------"""))
	
    ## parameters
	
    parser.add_argument("-c", "--cluster-composition", dest="cluster_composition", required=True, type=str, help="Directory with all the *clusterComposition.tsv files obtained with ReporTree and from which \
                        the clusters to analyze will be identified.")
    parser.add_argument("-snp", "--snp-dist", dest="snp_dist", required=True, type=str, help="Directory with all the SNP distance matrices to analyze (distance files must end with *_dist.tsv).")
    parser.add_argument("-allele", "--allele-dist", dest="allele_dist", required=True, type=str, help="Directory with all the allele distance matrices to analyze (distance files must end with *_dist.tsv).")
    parser.add_argument("-t", "--threshold", dest="threshold", required=True, type=str, help="Maximum threshold at which the clusters will be identified.")
    parser.add_argument("-o", "--out", dest="output", required=True, type=str, help="Tag for output files")
	
    args = parser.parse_args()

    print("\n******************************\n")
    print("Running stats_outbreak_analysis_snp_dists.py")
    print(version, "last updated on", last_updated)

    clusters = {}
    threshold = int(args.threshold)

    pipeline_allele_counter = 0
    for filename in glob.glob(args.cluster_composition + "/*_clusterComposition.tsv"):
        if "Serotype" not in filename or "ST" not in filename or "serotype" not in filename:
            if "/" in filename:
                fname = filename.split("/")[-1]
                pipeline_name = fname.split("_clusterComposition")[0]
            else:
                pipeline_name = filename.split("_clusterComposition")[0]
            
            pipeline_allele_counter += 1
            info = {}
            all_samples = set()
            with open(filename) as cluster_composition_file:
                lines = cluster_composition_file.readlines()

                # get clusters
                for line in lines[1:]:
                    lin = line.split("\n")[0]
                    l = lin.split("\t")
                    partition_name = l[0]
                    cluster = l[3]
                    cluster_len = int(l[2])
                    cluster_name = l[1]
                    if "single" in partition_name:
                        partition = partition_name.split("single-")[1]
                        partition = partition.split("x")[0]
                    elif "MST" in partition_name:
                        partition = partition_name.split("MST-")[1]
                        partition = partition.split("x")[0]
                    else:
                        print("Unknown method in partition name of pipeline " + str(pipeline_name))
                        sys.exit(1)
                        
                    if int(partition) <= threshold and cluster_len >= 2:
                        if partition not in clusters.keys():
                            clusters[partition] ={}
                        if pipeline_name not in clusters[partition].keys():
                            clusters[partition][pipeline_name] = []
                        info = cluster_name,cluster
                        clusters[partition][pipeline_name].append(info)

    with open(args.output + "_dists.tsv", "w+") as out3:
        print("partition\tpipeline_detected\tcluster_name\tcluster_len\tpipeline_dist\ttype_pipeline\tmin\tmax\tmed", file = out3)
        for partition in clusters.keys():
            for pipeline in clusters[partition].keys():
                for cluster in clusters[partition][pipeline]:
                    snp_pip = False
                    cluster_name,cluster_composition = cluster
                    samples = cluster_composition.split(",")
                    cluster_length = len(samples)
                    for filename in glob.glob(args.snp_dist + "*_dist.tsv"):
                        pipeline_snp = filename.split("/")[-1]
                        pipeline_snp_name = pipeline_snp.split("_ST")[0]
                        mx = pandas.read_table(filename)
                        dist_samples = mx[mx.columns[0]].values.tolist()
                        overlap = set(samples) & set(dist_samples)
                        if len(overlap) >= 1:
                            snp_pip = True
                            dists = []
                            for i in range(0,len(dist_samples)):
                                if dist_samples[i] in samples:
                                    flt_mx = mx[mx[mx.columns[0]] == dist_samples[i]]
                                    for s in mx.columns[1:]:
                                        if s in samples:
                                            d = flt_mx[s].values.tolist()[0]
                                            dists.append(d)
                            minimum = min(dists)
                            maximum = max(dists)
                            median = statistics.median(dists)
                            type_pipeline = "snp"
                            print(str(partition) + "\t" + pipeline + "\t" + cluster_name + "\t" + str(cluster_length) + "\t" + str(pipeline_snp_name) + "\t" + type_pipeline + "\t" + str(minimum) + "\t" + str(maximum) + "\t" + str(median))
                            print(str(partition) + "\t" + pipeline + "\t" + cluster_name + "\t" + str(cluster_length) + "\t" + str(pipeline_snp_name) + "\t" + type_pipeline + "\t" + str(minimum) + "\t" + str(maximum) + "\t" + str(median), file = out3)
                    for filename in glob.glob(args.allele_dist + "*_dist.tsv"):
                        pipeline_allele = filename.split("/")[-1]
                        pipeline_allele_name = pipeline_allele.split("_dist.tsv")[0]
                        if pipeline_allele_name != "allele3":
                            #if snp_pip:
                                mx = pandas.read_table(filename)
                                dist_samples = mx[mx.columns[0]].values.tolist()
                                overlap = set(samples) & set(dist_samples)
                                if len(overlap) >= 1:
                                    dists = []
                                    for i in range(0,len(dist_samples)):
                                        if dist_samples[i] in samples:
                                            flt_mx = mx[mx[mx.columns[0]] == dist_samples[i]]
                                            for s in mx.columns[1:]:
                                                if s in samples:
                                                    d = flt_mx[s].values.tolist()[0]
                                                    dists.append(d)
                                    minimum = min(dists)
                                    maximum = max(dists)
                                    median = statistics.median(dists)
                                    type_pipeline = "allele"
                                    print(str(partition) + "\t" + pipeline + "\t" + cluster_name + "\t" + str(cluster_length) + "\t" + str(pipeline_allele_name) + "\t" + type_pipeline + "\t" + str(minimum) + "\t" + str(maximum) + "\t" + str(median))
                                    print(str(partition) + "\t" + pipeline + "\t" + cluster_name + "\t" + str(cluster_length) + "\t" + str(pipeline_allele_name) + "\t" + type_pipeline + "\t" + str(minimum) + "\t" + str(maximum) + "\t" + str(median), file = out3)
if __name__ == "__main__":
    main()