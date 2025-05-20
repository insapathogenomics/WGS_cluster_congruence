#!/usr/bin/env	python3

"""
Created on Mon Sep  9 19:01:22 2024

@author: joana.gomes

The remove-hifen.py script is part of the WGS_cluster_congruence repository and takes as input ALL_CORRESPONDENCE.tsv file, obtained from get_best_part_correspondence.py.
The goal of this script is to remove the hyphens from the "method1" and "method2" columns in the ALL_CORRESPONDENCE.tsv file.
"""
#V1 JOANA
import argparse

def remove_hifen (input_file, output_file):

    try:
        with open(input_file, "r" ) as infile:
            lines=infile.readlines() #string list
            
            filtred_lines = []
            for line in lines:
                new_list = line.strip().split('\t') 
                if '-' not in new_list:
                    filtred_lines.append(line)
                   
            with open(output_file, "w") as outfile:
                outfile.writelines(filtred_lines)
       
    except FileNotFoundError:
        print("File was not found.")
        
def main ():

    parser = argparse.ArgumentParser(description="Remove hifen from ALL_CORRESPONDENCE.tsv file.")

    #input arguments
    parser.add_argument('-i','--input', action='store', required=True, help='[MANDATORY] Path for input ALL_CORRESPONDENCE.tsv file.')
    parser.add_argument('-o', '--output', action='store', required=True, help='[MANDATORY] Path and tag for output file name.')

    arguments = parser.parse_args()

    #function_call
    remove_hifen(arguments.input, arguments.output)

if __name__ == "__main__":
    main()