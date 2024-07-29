#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import math
from Genetic_coevolution_function import *


    
def Genetic_coevolution_program(ArgsVal):
    global args
    description = """ Study of phylogenetic coevolution of specific genes in brassicaceaes """
    parser = argparse.ArgumentParser(prog="Coevolution between genes",description=description)
    ## name of softwares used
    parser.add_argument("-m", "--muscle", default="muscle")
    parser.add_argument("-fi", "--phyml", default="phyml")
    parser.add_argument("-bedt", "--bedtools", default="bedtools")
    parser.add_argument("-bns", "--blastn", default="blastn")
    ## other parameters
    parser.add_argument("-D", "--directory", default="/auto/plzen1/home/levevea/pipeline_genetic_phylogeny_coevolution/")
    parser.add_argument("-P", "--parameter_file", default="Parameter_file.txt")
    parser.add_argument("-B", "--bed", default="candidate.bed")
    parser.add_argument("-R", "--ref", default="A.lyrata")
    parser.add_argument("-g", "--gff", default="Alyrata_384_v2.1.gene_exons.gff3")
    parser.add_argument("-S", "--sister", default="A.halleri")
    parser.add_argument("-e", "--eval_blastn", default="1e-5")
    parser.add_argument("-c", "--coverage", default="70")
    
    ## name of steps on pipeline
    parser.add_argument("-fcc", "--filter_candidate_coordonate", default="T")
    parser.add_argument("-ac", "--annotation_candidate", default="T")
    parser.add_argument("-bs", "--basic_stat", default="T")
    parser.add_argument("-plb", "--phyml_length_branch", default="T")
    parser.add_argument("-pp", "--phyml_prog", default="T")
    parser.add_argument("-pa", "--phyml_all", default="T")
    parser.add_argument("-pprep", "--phyml_preparation", default="T")
    parser.add_argument("-pe", "--presence_estimator", default="T")
    parser.add_argument("-ce", "--conservation_estimator", default="T")
    parser.add_argument("-mfs", "--muscle_fasta_species", default="T")
    parser.add_argument("-cfs", "--concatenate_fasta_species", default="T")
    parser.add_argument("-ebp", "--extraction_bedtool_species", default="T")
    parser.add_argument("-bbp", "--blastn_bed_species", default="T")
    parser.add_argument("-bc", "--blastn_candidate", default="T")
    parser.add_argument("-ecb", "--extraction_candidate_bedtool", default="T")

    args = parser.parse_args(ArgsVal)
    
    bed=str(args.bed)
    gff=str(args.gff)
    directory=str(args.directory)
    phyml_prog=str(args.phyml_prog)
    bedtools=str(args.bedtools)
    parameter=str(args.parameter_file)
    muscle=str(args.muscle)
    blastn=str(args.blastn)
    ref=str(args.ref)
    sister=str(args.sister)
    eval_blastn=str(args.eval_blastn)
    coverage=float(args.coverage)

    filter_candidate_coordonate_function=str(args.filter_candidate_coordonate)
    basic_stat_function=str(args.basic_stat)
    annotation_candidate_function=str(args.annotation_candidate)
    phyml_length_branch_function=str(args.phyml_length_branch)
    phyml_function=str(args.phyml_prog)
    phyml_all_function=str(args.phyml_all)
    phyml_preparation_function=str(args.phyml_preparation)
    presence_estimator_function=str(args.presence_estimator)
    conservation_estimator_function=str(args.conservation_estimator)
    muscle_fasta_species_function=str(args.muscle_fasta_species)
    concatenate_fasta_species_function=str(args.concatenate_fasta_species)
    extraction_bedtool_species_function=str(args.extraction_bedtool_species)
    blastn_bed_species_function=str(args.blastn_bed_species)
    blastn_candidate_function=str(args.blastn_candidate)
    extraction_candidate_bedtool_function=str(args.extraction_candidate_bedtool)
    

    Parameter=open(str(parameter),"r")
    Parameter=Parameter.read()
    Parameter=Parameter.split('#')
    REF=Parameter[1]
    REF=REF.split('\n')
    REF=REF[1:]
    
    Bed=open(str(directory)+"data/"+str(bed),"r")
    Bed=Bed.read()
    Bed=Bed.split('\n')
    Bed=Bed[1:]
    
    filter_candidate_coordonate(REF,ref,directory)
    if annotation_candidate_function=="T":annotation_candidate(REF,ref,directory)
    if extraction_candidate_bedtool_function=="T":extraction_candidate_bedtool(REF,ref,directory,bedtools)
    if blastn_candidate_function=="T":blastn_candidate(REF,ref,directory,blastn,eval_blastn)
    if basic_stat_function=="T":basic_stat(ref,gff,bedtools,REF,directory,Bed)
    if blastn_bed_species_function=="T":blastn_bed_species(REF,ref,directory,coverage,Bed)
    if extraction_bedtool_species_function=="T":extraction_bedtool_species(REF,ref,directory,bedtools)
    if concatenate_fasta_species_function=="T":concatenate_fasta_species(REF,ref,directory,Bed)
    if muscle_fasta_species_function=="T":muscle_fasta_species(muscle,Bed)
    if conservation_estimator_function=="T" :conservation_estimator(ref,Bed,sister)
    if presence_estimator_function=="T" :presence_estimator(ref,Bed,REF)
    if phyml_preparation_function=="T":phyml_preparation(ref,Bed,REF)
    if phyml_all_function=="T":phyml_all(ref,Bed,phyml)
    if phyml_function=="T":phyml(ref,Bed,phyml)
    if phyml_length_branch_function=="T" and phyml_all_function=="T":phyml_length_branch(REF,Bed,phyml)
    
if __name__=='__main__':
    
    ArgsVal = sys.argv[1:]
    Genetic_coevolution_program(ArgsVal)
