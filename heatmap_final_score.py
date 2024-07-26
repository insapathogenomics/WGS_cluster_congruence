#!/usr/bin/env  python3

"""
This script can be used to generate the heatmap of the congruence score obtained 
for the comparison of two pipelines at all threshold levels, as well as 
determine the tendency of the respective corresponding points.
It takes as input the ALL_CORRESPONDENCE.tsv file (obtained with 
get_best_part_correspondence.py) and the *final_score.tsv file obtained with 
comparing_partitions_v2.py. 

By Veronica Mixao
@INSA
"""

import sys
import argparse
import textwrap
import pandas
import seaborn as sn
import matplotlib.pyplot as plt
from scipy import stats

version = "1.0.0"
last_updated = "2024-07-26"

# functions	----------

def get_tendency(correspondence, pipeline1, pipeline2, threshold):
	""" This function generates the heatmap 
	input: filename
	output: seaborn heatmap
	"""

	# TENDENCY

	possible_comparison_names = [pipeline1 + "_vs_" + pipeline2, pipeline2 + "_vs_" + pipeline1]

	mx = pandas.read_table(correspondence)
	all_comparisons = pandas.unique(mx[mx.columns[0]])
	
	comparison = ""
	for comp1 in all_comparisons:
		if "_rev" not in comp1:
			if comp1 in possible_comparison_names:
				if threshold != "":
					extension = comp1 + "_thr" + str(threshold) + "_slope.tsv"
				else:
					extension = comp1 + "_slope.tsv"
				with open(extension, "w+") as out:
					print("#comp1\tcomp2\tslope\tintercept\tr_value\tp_value\tstd_err", file = out)
					for comp2 in all_comparisons:
						if "_rev" in comp2:
							if comp2.split("_rev")[0] == comp1:
								comps = [comp1, comp2]
								flt_mx = mx.loc[mx[mx.columns[0]].isin(comps)]
								if threshold != "":
									flt_mx = flt_mx[(flt_mx["method1"] <= int(threshold)) & (flt_mx["method2"] <= int(threshold))]
								#sn.lmplot(x="method1", y="method2", hue=mx.columns[0], data=flt_mx)
								#plt.savefig(comp1 + ".png") 
								if len(flt_mx[flt_mx.columns[0]].values.tolist()) == 0:
									print("No trend line will be provided as no congruence point was found!!")
								else:
									if len(flt_mx[flt_mx[flt_mx.columns[0]] == comp1]["method1"].values.tolist()) != 0:
										slope, intercept, r_value, p_value, std_err = stats.linregress(flt_mx[flt_mx[flt_mx.columns[0]] == comp1]["method1"],flt_mx[flt_mx[flt_mx.columns[0]] == comp1]["method2"])
										print(comp1 + "\t" + comp2 + "\t" + str(slope) + "\t" + str(intercept) + "\t" + str(r_value) + "\t" + str(p_value) + "\t" + str(std_err), file = out)
									if len(flt_mx[flt_mx[flt_mx.columns[0]] == comp2]["method1"].values.tolist()) != 0:
										slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(flt_mx[flt_mx[flt_mx.columns[0]] == comp2]["method1"],flt_mx[flt_mx[flt_mx.columns[0]] == comp2]["method2"])
										print(comp2 + "\t" + comp1 + "\t" + str(slope2) + "\t" + str(intercept2) + "\t" + str(r_value2) + "\t" + str(p_value2) + "\t" + str(std_err2), file = out)
									sn.lmplot(x="method1", y="method2", hue=mx.columns[0], data=flt_mx)
									if threshold != "":
										extension = "_thr" + str(threshold) + "_tendency.png"
									else:
										extension = "_tendency.png"
									plt.savefig(comp1 + extension)
									comparison = comp1
									plt.close()
	return comparison

def get_heatmap(filename,color,p1,p2,threshold):
	# HEATMAP
	mx = pandas.read_table(filename)
	mx = mx.drop(mx.columns[0], axis=1)
	comparison = str(p1) + "_vs_" + str(p2)
	if threshold != "":
		mx = mx.iloc[:, : int(threshold) + 1]
		mx = mx.head(int(threshold) + 1)
	mx = mx.set_axis(range(0,len(mx.columns)), axis=1)
	hm = sn.heatmap(data = mx, cmap = color)
	plt.xlabel(p2) 
	plt.ylabel(p1)
	hm.set(title="Congruence Score (" + str(p1) + " vs " + str(p2) + ")")
	#plt.savefig(comparison + "_final_score.png")
	#plt.show()
	fig = hm.get_figure()
	fig.set_size_inches(8, 6, forward=True)
	if threshold != "":
		fig.savefig(comparison + "_thr" + str(threshold) + "_final_score.png", format = "png", dpi=500)
	else:
		fig.savefig(comparison + "_final_score.png", format = "png", dpi=500)

def main():
	
	# argument options
    
	parser = argparse.ArgumentParser(prog="heatmap_final_score.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									################################################################################             
									#                                                                              #
									#                            heatmap_final_score.py                            #
									#                                                                              #
									################################################################################ 
									This script can be used to generate the heatmap of the congruence score obtained 
									for the comparison of two pipelines at all threshold levels, as well as 
									determine the tendency of the respective corresponding points.
									It takes as input the ALL_CORRESPONDENCE.tsv file (obtained with 
									get_best_part_correspondence.py) and the *final_score.tsv file obtained with 
									comparing_partitions_v2.py. 
									
									-------------------------------------------------------------------------------"""))
										
	parser.add_argument("--final-score", dest="final_score", required=False, type=str, help="[MANDATORY] *final_score.tsv output of comparing_partitions_v2.py with the CS matrix between two WGS pipelines.")
	parser.add_argument("--correspondence", dest="correspondence", required=True, type=str, help="[MANDATORY] Pivot table with the threshold correspondence (e.g. ALL_CORRESPONDENCE.tsv obtained with get_best_part_correspondence.py).")
	parser.add_argument("-t", "--threhold", dest="threshold", required=False, default="", help="[OPTIONAL] Maximum threshold to be used for plot")

	parser.add_argument("-c", "--color", dest="color", required=False, default="Blues", type=str, help="Seaborn color map. Default: Blues")
				
	args = parser.parse_args()
    
	filename = str(args.final_score).split("/")[-1]
	comparison_name = filename.split("_final_score.tsv")[0]
	pipeline1 = comparison_name.split("_vs_")[0]
	pipeline2 = comparison_name.split("_vs_")[1]

	get_tendency(args.correspondence, pipeline1, pipeline2, args.threshold)
	get_heatmap(filename, args.color, pipeline1, pipeline2, args.threshold)

if __name__ == "__main__":
    main()
