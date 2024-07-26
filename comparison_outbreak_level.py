#!/usr/bin/env	python3

"""
This script was designed to determine pipeline comparability at outbreak-level. 
It takes as input the cluster composition file of each pipeline and reports, 
for each cluster at a given level, what is the level at which the same cluster 
is detected in another pipeline.

By Veronica Mixao
@INSA
"""

import sys
import argparse
import textwrap
import pandas

version = "1.0.0"
last_updated = "2024-07-26"

def get_pipelines(infile):
	""" Get the list of pipelines to compare
	input: infile with the path to the cluster composition files
	output: dictionary """

	pipelines = {}
	types = {}

	with open(infile) as path:
		lines = path.readlines()
		for line in lines:
			f = line.split("\n")[0]
			f = f.split("\t")
			if "/" in f[0]:
				fname = f[0].split("/")[-1]
				pipeline_name = fname.split("_clusterComposition")[0]
			else:
				pipeline_name = f[0].split("_clusterComposition")[0]
			print("\t" + pipeline_name)
			pipelines[pipeline_name] = f[0]
			types[pipeline_name] = f[1]
	
	return pipelines, types

def get_samples_exclude(pipelines, types):
	""" Determine samples to exclude from the analysis, i.e. absent from at least one pipeline
	input: dictionary
	output: set """

	info = {}
	all_samples = set()
	exclude_samples = set()

	for pipeline in pipelines.keys():
		if types[pipeline] == "allele":
			with open(pipelines[pipeline]) as cluster_composition_file:
				lines = cluster_composition_file.readlines()
				line_interest = lines[-1]
				lin = line_interest.split("\n")[0]
				l = lin.split("\t")
				cluster = l[3]
				info[pipeline] = set()
				for sample in cluster.split(","):
					info[pipeline].add(sample)
					all_samples.add(sample)

	for pipeline in info.keys():
		pipeline_samples = info[pipeline]
		for sample in all_samples:
			if sample not in pipeline_samples:
				exclude_samples.add(sample)

	return exclude_samples

def get_clusters(pipelines,exclude_samples):
	""" Identify the clusters which will be used for comparison
	input: dictionary with the files to retrieve clusters from
	output: dictionary """

	clusters = {} # clusters[pipeline][clusters] = partition
	partitions = {} # partitions[pipeline][partition] = [clusters]

	for pipeline in pipelines.keys():
		clusters[pipeline] = {}
		partitions[pipeline] = {}
		with open(pipelines[pipeline]) as cluster_composition_file:
			lines = cluster_composition_file.readlines()
			for line in lines[1:]:
				lin = line.split("\n")[0]
				l = lin.split("\t")
				partition_name = l[0]
				cluster = l[3]
				if "single" in partition_name:
					partition = partition_name.split("single-")[1]
					partition = partition.split("x")[0]
				elif "MST" in partition_name:
					partition = partition_name.split("MST-")[1]
					partition = partition.split("x")[0]
				else:
					print("Unknown method in partition name of pipeline " + str(pipeline))
					sys.exit(1)
				if int(partition) not in partitions[pipeline].keys():
					partitions[pipeline][int(partition)] = []
				
				samples = cluster.split(",")
				cluster = list(set(samples) - exclude_samples)
				cluster.sort()
				
				if len(cluster) > 1:
					cluster = ",".join(cluster)
					partitions[pipeline][int(partition)].append(cluster)
					if cluster not in clusters[pipeline].keys():
						clusters[pipeline][cluster] = partition

	return clusters, partitions

def get_clusters_of_interest(partitions,threshold,types):
	""" Get the set of clusters of interest
	input: partitions dictionary and the threshold where to search
	output: set """

	clusters_of_interest = set()
	for pipeline in partitions.keys():
		if types[pipeline] == "allele":
			for cluster in partitions[pipeline][int(threshold)]:
				clusters_of_interest.add(cluster)
	
	return clusters_of_interest

def get_matrix(clusters_of_interest,clusters,types):
	""" generate the final matrix
	input: clusters of interest and the clusters dictionary
	output: matrix """
	
	mx = {}
	mx["cluster"] = []
	mx["cluster_length"] = []
	mx["n_allele_pipelines"] = []
	pipeline_lst = []

	for cluster in clusters_of_interest:
		mx["cluster"].append(cluster)
		cluster_counter = 0
		for pipeline in clusters.keys():
			if pipeline not in mx.keys():
				mx[pipeline] = []
				pipeline_lst.append(pipeline)
			if cluster in clusters[pipeline]:
				mx[pipeline].append(clusters[pipeline][cluster])
				if types[pipeline] == "allele":
					cluster_counter += 1
			else:
				mx[pipeline].append("-")
		mx["cluster_length"].append(len(cluster.split(",")))
		mx["n_allele_pipelines"].append(str(cluster_counter))
	
	return mx, pipeline_lst

def get_pivot(mx, pipeline_lst):
	""" generate the pivot of the final matrix
	input: dictionary
	output: matrix 
	"""

	pivot_mx = {"n_allele_pipelines": [], "cluster_length": [], "pipeline": [], "threshold": []}

	for i in range(0,len(mx["cluster_length"])):
		for pipeline in pipeline_lst:
			if mx[pipeline][i] != "-":
				pivot_mx["n_allele_pipelines"].append(mx["n_allele_pipelines"][i])
				pivot_mx["cluster_length"].append(mx["cluster_length"][i])
				pivot_mx["pipeline"].append(pipeline)
				pivot_mx["threshold"].append(mx[pipeline][i])
	
	pivot_mx_df = pandas.DataFrame.from_dict(pivot_mx)

	return pivot_mx_df

def main():
    
	# argument options	----------
    
	parser = argparse.ArgumentParser(prog="comparison_outbreak_level.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                        comparison_outbreak_level.py                         #
									#                                                                             #
									############################################################################### 
									                            
									This script was designed to determine pipeline comparability at outbreak-level. 
									It takes as input the cluster composition file of each pipeline and reports, 
									for each cluster at a given level, what is the level at which the same cluster 
									is detected in another pipeline.

									Briefly:
									1. Create a list of genetic clusters
									1.1 Pipeline 1: identify all clusters at threshold X
									1.2 Pipeline 2: identify all clusters at threshold X
									1.3 ...

									2. For each genetic cluster, determine the threshold at which it is identified
									by each pipeline

									3. Present the results in a tsv matrix

									                  
									-------------------------------------------------------------------------------"""))
	
	## parameters
	
	parser.add_argument("-i", "--input", dest="infile", required=True, type=str, help="Input *tsv file containing the path to each *clusterComposition.tsv file that has to be compared and the type of pipeline, i.e. allele or snp (cluster composition files can be obtained with ReporTree)")
	parser.add_argument("-o", "--out", dest="output", required=True, type=str, help="Tag for output file")
	parser.add_argument("-t", "--threshold", dest="threshold", required=True, type=int, help="Threshold at which the genetic clusters must be identified")
	
	args = parser.parse_args()
	
	print("\n******************************\n")
	print("Running comparison_outbreak_level.py")
	print(version, "last updated on", last_updated)

	## pipeline
	
	print("\nDetecting pipelines...")
	pipelines, types = get_pipelines(args.infile)

	print("Identifying samples that were absent from at least one allele pipeline...")
	exclude_samples = get_samples_exclude(pipelines, types)

	print("Getting cluster information...")
	clusters, partitions = get_clusters(pipelines,exclude_samples)

	print("Determining clusters of interest...")
	clusters_of_interest = get_clusters_of_interest(partitions,args.threshold,types)

	print("Pipeline comparison...")
	mx, pipeline_lst = get_matrix(clusters_of_interest,clusters,types)
	mx_df = pandas.DataFrame.from_dict(mx)

	print("Generating pivot table...")
	mx_pivot = get_pivot(mx, pipeline_lst)

	mx_df.to_csv(args.output + ".tsv", index = False, header=True, sep ="\t")
	mx_pivot.to_csv(args.output + "_pivot.tsv", index = False, header=True, sep ="\t")

if __name__ == "__main__":
    main()