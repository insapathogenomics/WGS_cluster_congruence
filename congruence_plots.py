#!/usr/bin/env  python3

"""
This script generates all the plots used for the BeONE cluster congruence analysis

By Veronica Mixao
@INSA
"""

import argparse
import textwrap
import glob
import pandas
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

version = "1.0.0"
last_updated = "2024-07-26"

def get_n_clusters(directory):
    """ This function determines the number of clusters found in each partition for all partitions file in a directory 
    input: directory
    output: dataframe
    """
    
    info_partitions = {"pipeline": [], "threshold": [], "partitions": []}
    order_col = ["pipeline", "threshold", "partitions"]

    for filename in glob.glob(directory + "/*_partitions.tsv"):
        fname = filename.split("/")[-1]
        pipeline = fname.split("_partitions.tsv")[0]
        partitions = pandas.read_table(filename)
        for i in range(1,len(partitions.columns)):
            clusters = pandas.unique(partitions[partitions.columns[i]])
            info_partitions["pipeline"].append(pipeline)
            info_partitions["threshold"].append(i)
            info_partitions["partitions"].append(len(clusters))
    
    matrix = pandas.DataFrame(data = info_partitions, columns = order_col)

    return matrix

def get_log_stability(directory):
	""" This function transforms the stableRegion files into a single matrix with the region coordenated represented in log2 
	input: directory
	output: dataframe
	"""

	info_blocks = {"pipeline": [], "block": [], "pipeline_code": [], "partition": []}
	order_col = ["pipeline", "block", "pipeline_code", "partition"]
	order_blocks = []
	for filename in glob.glob(directory + "/*_stableRegions.tsv"):
		fname = filename.split("/")[-1]
		pipeline = fname.split("_stableRegions.tsv")[0]
		block_code = 1
		with open(filename) as infile:
			lines = infile.readlines()
			for line in lines:
				if "#" not in line:
					lin = line.split("\n")[0]
					l = lin.split("\t")
					block = l[0]
					pipeline_code = pipeline + "_" + str(block_code)
					start_info = l[1]
					start_info = start_info.split("-")[1]
					start_info = start_info.split("x")[0]
					end_info = l[2]
					end_info = end_info.split("-")[1]
					end_info = end_info.split("x")[0]

					if int(l[3]) >= 10:
						start_log2 = math.log2(int(start_info))
						end_log2 = math.log2(int(end_info))
						#start_log2 = int(start_info)
						#end_log2 = int(end_info)

						# line with start
						info_blocks["pipeline"].append(pipeline)
						info_blocks["block"].append(block)
						info_blocks["pipeline_code"].append(pipeline_code)
						info_blocks["partition"].append(start_log2)

						# line with "median"
						info_blocks["pipeline"].append(pipeline)
						info_blocks["block"].append(block)
						info_blocks["pipeline_code"].append(pipeline_code)
						info_blocks["partition"].append(end_log2)		

						# line with "median"
						info_blocks["pipeline"].append(pipeline)
						info_blocks["block"].append(block)
						info_blocks["pipeline_code"].append(pipeline_code)
						info_blocks["partition"].append(end_log2)
						block_code += 1
						order_blocks.append(pipeline_code)
	
	matrix = pandas.DataFrame(data = info_blocks, columns = order_col)
	
	return matrix, order_blocks

def get_ST_info(filename, col, col_n):
	""" This function counts the number of samples per ST and generates a matrix with this info
	input: directory
	output: dataframe
	"""

	counts = {col: [], col_n: []}

	info = pandas.read_table(filename, dtype={col:"object"})
	info[col].replace('', np.nan, inplace = True)
	info.dropna(subset=[col], inplace = True)
	ST_lst = info[col].values.tolist()
	ST = set(ST_lst)
	for st in ST:
		c = ST_lst.count(st)
		if c >= 50:
			counts[col].append(st)
			counts[col_n].append(c)
	
	matrix = pandas.DataFrame(data = counts)
	matrix = matrix.sort_values(by=[col_n], ascending=False)

	return matrix

def get_polityping_info(filename, min_length):
	""" This function transforms the polityping output
	input: directory
	output: dataframe
	"""

	info = {}

	i = 0
	with open(filename) as infile:
		lines = infile.readlines()
		for line in lines:
			lin = line.split("\n")[0]
			l = lin.split("\t")

			if i == 0:
				main_col = l[0]
				size_col = l[1]
				info[main_col] = []
				info[size_col] = []
				info["pipeline"] = []
				info["threshold"] = []
				pipelines = l[2:]
			else:
				if int(l[1]) >= min_length:
					for k in range(0,len(pipelines)):
						if l[k+2] != "-":
							info[main_col].append(l[0])
							info[size_col].append(int(l[1]))
							info["pipeline"].append(pipelines[k])
							info["threshold"].append(int(l[k+2]))
			i += 1

	matrix = pandas.DataFrame(data = info)
	matrix = matrix.sort_values(by=[size_col], ascending=False)

	return matrix, main_col, pipelines


def main():
		# argument options
    
	parser = argparse.ArgumentParser(prog="congruence_plots.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                             congruence_plots.py                             #
									#                                                                             #
									############################################################################### 
									This script generates useful plots for the WGS cluster congruence analysis.
																																				
									Required input files:
									- Plot Clusters_per_partition: directory with *_partitions.tsv from ReporTree
									- Plot Stability: directory with *_stableRegions.tsv from ReporTree
									- ST_barplot: *metadata_w_partitions.tsv file in "-f" and column name in "-col"
									- ST_swarmplot: *_poli.tsv file from poli_typing.py
									- ST_boxplot: *_poli.tsv from poli_typing.py
									-------------------------------------------------------------------------------"""))
										
	parser.add_argument("-d", "--directory", dest="directory", required=True, type=str, help="[MANDATORY] Directory with the input files.")
	parser.add_argument("-f", "--file", dest="filename", required=False, type=str, help="[ONLY REQUIRED FOR ST_barplot or ST_swarmplot or ST_boxplot] Filename.")
	parser.add_argument("-c", "--color", dest="color", required=False, default="Blues", type=str, help="Seaborn color map. Default: Blues")
	parser.add_argument("-t", "--tag", dest="tag", required=False, default="Congruence", type=str, help="Tag for output")
	parser.add_argument("-p", "--plot", dest="plot", required=False, default="Clusters_per_partition", type=str, help="Type of plot. Options: Clusters_per_partition, Stability, ST_barplot, ST_swarmplot, ST_boxplot")   
	parser.add_argument("-col", "--column", dest="column", required=False, type=str, help="[ONLY REQUIRED FOR ST_barplot] Column name (e.g. MLST_ST).")   

	args = parser.parse_args()

	print("\n******************************\n")
	print("Running beone_congruence_plots.py")
	print(version, "last updated on", last_updated)

	### get line plot with the number of clusters per partition
	if args.plot == "Clusters_per_partition":
		matrix = get_n_clusters(args.directory)
		matrix.to_csv(args.tag + "_clusters_partition.tsv", index = False, header=True, sep ="\t")
		lineplot = sns.lineplot(data=matrix,  x="threshold", y="partitions", hue="pipeline", palette="deep")

		fig = lineplot.get_figure()
		fig.savefig(args.tag + "_lineplot.svg", format = "svg")

	### get boxplot with the stability blocks
	if args.plot == "Stability":
		matrix, order_blocks = get_log_stability(args.directory)
		matrix = matrix.sort_values(by=["pipeline_code"])
		box_plot = sns.boxplot(data=matrix, x="partition", y="pipeline_code", order=order_blocks)
		fig = box_plot.get_figure()
		fig.savefig(args.tag + "_stability.svg", format = "svg") 
    
	### get ST barplot
	elif args.plot == "ST_barplot":
		col = args.column
		col_n = args.column + "_n"
		matrix = get_ST_info(args.directory + "/" + args.filename, col, col_n)
		order_ST = matrix[col].values.tolist()
		bar_plot = sns.barplot(data=matrix, x=col, y=col_n, order=order_ST)
        
		fig = bar_plot.get_figure()
		fig.savefig(args.tag + "_barplot.svg", format = "svg")
	
	### get ST swarmplot
	elif args.plot == "ST_swarmplot":
		matrix, main_col, order_pipelines = get_polityping_info(args.directory + "/" + args.filename, 50)
		matrix = matrix.sort_values(by=matrix.columns[1], ascending=False)
		#matrix = matrix.sort_values(by=["pipeline"], key=lambda col: col.str.lower())
		swarmplot = sns.swarmplot(data=matrix, x=main_col, y="threshold", hue="pipeline", hue_order = order_pipelines, palette="deep")
		fig = swarmplot.get_figure()
		fig.savefig(args.tag + "_swarmplot.svg", format = "svg")

    ### get ST boxplot
	elif args.plot == "ST_boxplot":
		matrix, main_col, order_pipelines = get_polityping_info(args.directory + "/" + args.filename, 2)
		box_plot = sns.boxplot(data=matrix, x="threshold", y="pipeline", order=order_pipelines)
        
		fig = box_plot.get_figure()
		fig.savefig(args.tag + "_boxplot.svg", format = "svg")

if __name__ == "__main__":
    main()