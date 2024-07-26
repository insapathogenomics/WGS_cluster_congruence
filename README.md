# WGS cluster congruence
This repository contains the scripts used for the [comprehensive cluster congruence analysis of multiple bioinformatics pipelines](https://www.medrxiv.org/content/10.1101/2024.07.24.24310933v1) used in WGS-based surveillance of important foodborne bacterial species.

## Input
These analyses rely on the usage of two other important tools, which are the basis of all work:
- [ReporTree](https://github.com/insapathogenomics/ReporTree) - a flexible tool for the identification of genetic clusters at all possible threshold levels and the integration of this information with any metadata of interest. 
- [ComparingPartitions](https://github.com/insapathogenomics/ComparingPartitions) - a tool originally developed by [João Carriço](https://github.com/jacarrico) and that we adapted in order to respond to the needs of this and other extensive works. The script [comparing_partitions_v2.py](https://github.com/insapathogenomics/ComparingPartitions/blob/master/comparing_partitions_v2.py) was used to obtain the Congruence Score (CS) score between any two typing methods.

 The scripts present in this repository take many of the outputs of these tools as input.

 ## Scripts description

### Cluster congruence between WGS pipelines and traditional typing
- _poli_typing.py_ - This script can be used to identify the minimum partition required to cluster all samples within the same typing category (e.g. ST, serotype...).
- _get_stats_threshold.py_ - This script sumarizes the distribution of a variable within the clusters at a given threshold.

### Cluster congruence between WGS pipelines
- _get_best_part_correspondence.py_ - This script can be used for each pairwise pipeline comparison in order to assess what is the threshold that provides the most similar clustering results in the other pipeline (i.e., the best “correspondence point”), as assessed by CS scores.
- _heatmap_final_score.py_ - This script can be used to generate the heatmap of the congruence score obtained for the comparison of two pipelines at all threshold levels, as well as determine the tendency of the respective corresponding points.

### Concordance at high resolution level 
- _comparison_outbreak_level.py_ - This script can be used to assess for each cluster at a given level, what is the level at which the same cluster is detected in another pipeline.
- _stats_outbreak_analysis.py_ - This script can be used to determine the number of clusters identified by a given pipeline at a given threshold that are also detected with the exact same composition by another pipeline at a (or up to a) similar or even higher threshold.
- _stats_outbreak_analysis_snp_dists.py_ - This script can be used to determine, for each cluster identified by a given pipeline at a given threshold, the maximum allelic or SNP distance within the same cluster that the same or an alternative allele- or SNP-based determined.
- _wgmlst_exercise.py_ - This script was used to determine the gain in resolution provided by ReporTree dynamic zoom-in analysis in potential outbreak clusters.

### General plots
- _congruence_plots.py_ - This script can be used to generate some of the plots presented in the work

 
## Citation
If you use any of the scripts available in this repository, please cite our work:

[Mixão, V. et al. (2024) Multi-country and intersectoral assessment of cluster congruence between different bioinformatics pipelines for genomics surveillance of foodborne bacterial pathogens. medRxiv 2024.07.24.24310933. doi: https://doi.org/10.1101/2024.07.24.24310933](https://www.medrxiv.org/content/10.1101/2024.07.24.24310933v1)
