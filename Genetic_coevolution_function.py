#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import os
import sys
import argparse
import math
import random

# functions used for coevolution of genes analysis

##

## filter coordonate of candidate

def filter_candidate_coordonate(REF,ref,directory):
    for species in REF:
      ref2=species.split("\t")
      if str(ref2[0])==ref :
        x2=str(ref2[0])+"_filtrered.gff"
        fichier=open(str(x2), "w")
        x="Chrom\tstart\tend"
        fichier.write(str(x))  
        gff=str(directory)+"Database/"+str(ref2[2])
        gff=open(str(gff), "r")
        gff=gff.read()
        gff=gff.split('\n')
        for i in gff:
          if i.count("#")==0 and i!="":
            seq2=i.split("\t")
            if len(seq2)>3 :
              if str(seq2[2])=="gene"  :
                x=str(seq2[0])+"\t"+str(seq2[3])+"\t"+str(seq2[4])
                fichier.write("\n"+str(x) )   
        fichier.close() 
      
    print("filter_candidate_coordonate ok")

## annotation of sites in candidate genes in reference

def annotation_candidate(REF,ref,directory):
    
    for species in REF:
        ref2=species.split("\t")
        if str(ref2[0])== ref:
            cmd="python NewAnnotateRef2.py "+str(directory)+"Database/"+str(ref2[1])+" "+str(directory)+"Database/"+str(ref2[2])+" -o > "+str(ref)+"_annotation.csv"
            os.system(cmd)
            DP = open(str(ref)+"_annotation.csv", "r")
            DP=DP.read()
            DP = DP.split("\n")
            fichier=open(str(ref)+"_annotation_filtered.csv", "w")
            fichier.write("Chrom pos ref gene sens annot")
            for j in range(1,len(DP)):
                annot=DP[j]
                annot=annot.split(" ")
                if len(annot)>=5:
                    annot2=annot[5]
                    annot2=annot2.split(",")
                if annot2.count("3")!=0:fichier.write("\n"+str(DP[j]))
                if annot2.count("4")!=0 and annot2.count("3")==0:fichier.write("\n"+str(DP[j]))
            fichier.close()

## extraction of sequences of candidate genes with bedtools

def extraction_candidate_bedtool(REF,ref,directory,bedtools):
    for species in REF:
        Species=species.split("\t")
        if Species[0]==ref:
          cmd=str(bedtools)+" getfasta -fi "+str(directory)+"Database/"+str(Species[1])+" -bed "+str(directory)+"/data/"+str(bed)+" > "+str(ref)+"_study_filtrered.fa"
          os.system(cmd)

## blastn of each candidate genes in the 25 others species
          
def blastn_candidate(REF,ref,directory,blastn,eval_blastn):
    for species in REF:
        Species=species.split("\t")
        if Species[0]!=ref:
            cmd=str(blastn)+" -subject "+str(directory)+"Database/"+str(Species[1])+" -query "+str(directory)+"/"+str(ref)+"_study_filtrered.fa"+" -out "+str(ref)+"_"+str(Species[0])+"_study_filtrered.txt"+" -evalue "+str(eval_blastn)+"  -outfmt 6"
            os.system(cmd)



## generate bed files with better candidate in the 25 others species
            
def blastn_bed_species(REF,ref,directory,coverage,Bed):
    for species in REF:
        ref2=species.split("\t")
        if str(ref2[0])!= ref:
            x2=str(ref2[0])+".bed"
            fichier=open(str(x2), "w")
            fichier.write("#chrom\tstart\tend\tsens\trefname")
    
            blastn_file=str(str(directory)+str(ref)+"_"+str(ref2[0])+"_study_filtrered.txt")
            blastn_file=open(str(blastn_file), "r")
            blastn_file=blastn_file.read()
            blastn_file = blastn_file.split("\n")
            candidate={}
            for i in blastn_file:
                if i !="":
                  gene=i.split("\t")
                  ID=gene[0]
                  Gene=gene[0]
                  Gene=Gene.split(":")
                  Gene=Gene[1]
                  Gene=Gene.split("-")
                  lenght=(int(Gene[1])-int(Gene[0]))
                  cover=int(gene[7])-int(gene[6])-int(gene[5])
                  cover=cover*100/lenght
                  if cover >= coverage:
                    x=[i]
                    if ID in candidate.keys():
                      x=candidate[ID]
                      x.append(i)
                    candidate[ID]=x
            
            for i in Bed:
                seq=i.split("\t")
                ID=str(seq[0])+":"+str(seq[1])+"-"+str(seq[2])
                if ID in candidate.keys() :
                    candidate2=candidate[ID]
                    if len(candidate2)>1:
                        Eval=[]
                        Length_align=[]
                        IDD=[]
                        for c in candidate2:
                          gene=c.split("\t")
                          Eval.append(float(gene[10]))
                          Length_align.append(int(int(gene[7])-int(gene[6])-int(gene[5])))
                          IDD.append(float(gene[2]))
                        candidate3=[]
                        for c in candidate2:
                          gene=c.split("\t")
                          if float(gene[10])<= min(Eval):candidate3.append(c)
                        candidate2=candidate3[:]
                        if len(candidate2)>1:
                          candidate3=[]
                          for c in candidate2:
                            gene=c.split("\t")
                            if int(int(gene[7])-int(gene[6])-int(gene[5]))>= max(Length_align):candidate3.append(c)
                            IDD.append(float(gene[2]))
                          candidate2=candidate3[:]
                          if len(candidate2)>1:
                            candidate3=[]
                            for c in candidate2:
                              gene=c.split("\t")
                              if float(gene[2])>= min(IDD):candidate3.append(c)
                            candidate2=candidate3[:]
                            if len(candidate2)>1:
                              random.shuffle(candidate2)
                    if len(candidate2)>0:         
                      candidate2=candidate2[0]    
                      gene=candidate2.split("\t")            
                      deb=int(gene[8])
                      fin=int(gene[9])
                      sens="+"
                      if deb>fin:
                          fin=int(gene[8])
                          deb=int(gene[9])
                          sens="-"
                      x=str(gene[1])+"\t"+str(deb)+"\t"+str(fin)+"\t"+str(sens)+"\t"+str(gene[0])
                      fichier.write("\n"+str(x))
            fichier.close()           
          
## extraction of sequences of candidate genes in the 25 others species with bedtools
            
def extraction_bedtool_species(REF,ref,directory,bedtools):
    for species in REF:
        ref2=species.split("\t")
        if str(ref2[0])!= ref:
            x2=str(ref2[0])+".bed"
            cmd=str(bedtools)+" getfasta -fi "+str(directory)+"Database/"+str(ref2[1])+" -bed "+str(directory)+str(x2)+" > "+str(directory)+str(ref2[0])+"_study_filtrered.fa"
            os.system(cmd)
    
## creation of fasta files with all seq of candidate genes found in the 26 species

def concatenate_fasta_species(REF,ref,directory,Bed):
    
    ref_fasta_file=str(str(ref)+"_study_filtrered.fa")
    ref_fasta_file=open(str(ref_fasta_file), "r")
    ref_fasta_file=ref_fasta_file.read()
    ref_fasta_file = ref_fasta_file.split(">")
    dict_annotation={}
    for line in ref_fasta_file[1:]:
      if line!="":
        line_splitted = line.split("\n") 
        test2=str(line_splitted[0])
        dict_annotation[test2]=str(line_splitted[1])

    for i in Bed:
      seq=i.split("\t")
      x2=str(seq[3])+".fa"
      fichier=open(str(x2), "w")
      seq=str(seq[0])+":"+str(seq[1])+"-"+str(seq[2])
      seq2=dict_annotation[seq]
      fichier.write(">"+str(ref))
      fichier.write("\n"+str(seq2))
      for species in REF:
        ref2=species.split("\t")
        if str(ref2[0])!= ref:
            bed_file=str(str(ref2[0])+".bed")
            bed_file=open(str(bed_file), "r")
            bed_file=bed_file.read()
            if bed_file.count(str(seq))!=0:
              bed_file = bed_file.split("\n")
              i2=1
              seq3=bed_file[i2]
              while seq3.count((str(seq)))==0 and i2<(len(bed_file)-1):
                i2=i2+1
                seq3=bed_file[i2]
              seq3=seq3.split("\t")
              sens=seq3[3]
              sens=sens.split(";")
              sens=str(sens[0])
              seq3=str(seq3[0])+":"+str(seq3[1])+"-"+str(seq3[2])
              ref_fasta_file=str(str(directory)+str(ref2[0])+"_study_filtrered.fa")
              ref_fasta_file=open(str(ref_fasta_file), "r")
              ref_fasta_file=ref_fasta_file.read()
              ref_fasta_file = ref_fasta_file.split(">")
              i2=1
              seq4=ref_fasta_file[i2]
              while seq4.count((str(seq3)))==0 and i2<(len(ref_fasta_file)-1):
                i2=i2+1
                seq4=ref_fasta_file[i2]
              fichier.write("\n"+">"+str(ref2[0]))
              seq4=seq4.split("\n")
              seq4=seq4[1]
              seq5=str(seq4)
              if sens=="-":
                if str(seq4[0])=="A" or str(seq4[0])=="a" :seq5="T"
                elif str(seq4[0])=="T" or str(seq4[0])=="t" :seq5="A"
                elif str(seq4[0])=="C" or str(seq4[0])=="c" :seq5="G"
                elif str(seq4[0])=="G" or str(seq4[0])=="g" :seq5="C"
                else:seq5="-"
                for i2 in seq4[1:]:
                  if str(i2)=="A" or str(i2)=="a" :seq5="T"+str(seq5)
                  elif str(i2)=="T" or str(i2)=="t" :seq5="A"+str(seq5)
                  elif str(i2)=="C" or str(i2)=="c" :seq5="G"+str(seq5)
                  elif str(i2)=="G" or str(i2)=="g" :seq5="C"+str(seq5)
                  else:seq5="-"+str(seq5)
              fichier.write("\n"+str(seq5))
      fichier.close()
    
## alignement muscle of each fasta of candidates genes 

def muscle_fasta_species(muscle,Bed):
    j=0
    for i in Bed:
      seq=i.split("\t")
      x2=str(seq[3])+".fa"
      cmd=str(muscle)+" -in "+str(x2)+" -out "+str(seq[3])+"_alig.fa"
      os.system(cmd)
      j=j+1

## conservation estimator
      
def conservation_estimator(ref,Bed,sister):
    
    file1=str(ref)+"_annotation_filtered.csv"
    file1 = open(str(file1), "r")
    file1=file1.read()
    file1 = file1.split("\n")
    dict_annotation2={}

    for line in file1[1:]:
      if line!="":
        line_splitted = line.split(" ")
        test2=str(str(line_splitted[0])+"_"+str(line_splitted[1]))
        exec("dict_annotation2['"+str(test2)+"']=str(line_splitted[5])")

    
    name_ref=ref
    name_sister=sister
    
    fichier=open(str(ref)+"_conservation_by_pos_all.csv", "w")
    fichier.write("Gene;Chrom;Position;ref;Annotation;alt_sister;Divergence_D sister;Number_of_seq_considered;Divergence_D all")
    
    for i in Bed:
      seq=i.split("\t")
      Gene_name=str(seq[3])
      Chrom=str(seq[0])
      Start=int(seq[1])
      End=int(seq[2])
      x2=str(seq[3])+"_alig.fa"
      if os.path.exists(str(x2)):
        fasta_alig_file = open(str(x2), "r")
        fasta_alig_file=fasta_alig_file.read()
        fasta_alig_file = fasta_alig_file.split(">")
        seq_ref=""
        seq_other=[]
        seq_sister=""
        for k in fasta_alig_file[1:]:
          seq=k
          seq2=seq.split("\n")
          name=seq2[0]
          seq3=seq.replace(name,"")
          seq3=seq3.replace("\n","")
          if seq2[0]==ref:seq_ref=seq3
          elif seq2[0]==name_sister:
            seq_sister=seq3
            seq_other.append(seq3)
          else: seq_other.append(seq3)
          
        Position=Start
        for j in range(0,len(seq_ref)):
          SEQ=str(seq_ref[j])
          if SEQ!="-":
            annot="Other"
            search_annotation=str(Chrom)+"_"+str(Position)
            if search_annotation in dict_annotation2:
              x=dict_annotation2[search_annotation]
              if x.count("3")!=0:annot="0fold"
              elif x.count("3")==0 and x.count("4")!=0:annot="4fold"
            fichier.write("\n"+str(Gene_name)+";"+str(Chrom)+";"+str(Position)+";"+str(SEQ)+";"+str(annot))
          
            if seq_sister!="":
              if str(seq_sister[j])!=SEQ:
                if  str(seq_sister[j])!="-": fichier.write(";"+str(seq_sister[j])+";1")
                elif  str(seq_sister[j])=="-": fichier.write(";"+str(seq_sister[j])+";NA")
              if str(seq_sister[j])==SEQ:fichier.write(";"+str(seq_sister[j])+";0")
            if seq_sister=="":fichier.write(";NA;NA")
            
            D_pos=[]
            for other in seq_other:
              if other[j]!=SEQ and str(other[j])!="-":D_pos.append(1)
              elif other[j]!=SEQ:D_pos.append(0)
            if len(D_pos)==0: D_pos="NA"
            elif len(D_pos)!=0: D_pos=sum(D_pos)/len(D_pos)
            fichier.write(";"+str(len(seq_other))+";"+str(D_pos))
            Position=Position+1
      print(Gene_name)
    fichier.close()
    
    file1=str(ref)+"_conservation_by_pos_all.csv"
    file1 = open(str(file1), "r")
    file1=file1.read()
    file1 = file1.split("\n")
    dict_annotation2={}
    dict_annotation_all={}
    dict_annotation_NS={}
    dict_annotation_NS_all={}
    dict_annotation_S={}
    dict_annotation_S_all={}
    dict_annotation_gene={}

    for line in file1[1:]:
      if line!="":
        line_splitted = line.split(";")
        test2=str(line_splitted[0])
        annot=str(line_splitted[4])
        sister=str(line_splitted[6])
        all=str(line_splitted[8])
        nb=int(line_splitted[7])
        
        if test2 not in dict_annotation_gene.keys():
          if sister != "NA" and nb>0  :dict_annotation_gene[test2]=[1,nb]
          elif sister == "NA" and nb>0  :dict_annotation_gene[test2]=[0,nb]
          elif sister == "NA" and nb==0  :dict_annotation_gene[test2]=[0,0]
        
        if sister != "NA":
          x=[float(line_splitted[6])]
          if test2 in dict_annotation2.keys():
            x=dict_annotation2[test2]
            x.append(float(line_splitted[6]))
          dict_annotation2[test2]=x
          if annot=="0fold":
            x=[float(line_splitted[6])]
            if test2 in dict_annotation_NS.keys():
              x=dict_annotation_NS[test2]
              x.append(float(line_splitted[6]))
            dict_annotation_NS[test2]=x
          if annot=="4fold":
            x=[float(line_splitted[6])]
            if test2 in dict_annotation_S.keys():
              x=dict_annotation_S[test2]
              x.append(float(line_splitted[6]))
            dict_annotation_S[test2]=x
            
        if all != "NA":
          x=[float(line_splitted[8])]
          if test2 in dict_annotation_all.keys():
            x=dict_annotation_all[test2]
            x.append(float(line_splitted[8]))
          dict_annotation_all[test2]=x
          if annot=="0fold":
            x=[float(line_splitted[8])]
            if test2 in dict_annotation_NS_all.keys():
              x=dict_annotation_NS_all[test2]
              x.append(float(line_splitted[8]))
            dict_annotation_NS_all[test2]=x
          if annot=="4fold":
            x=[float(line_splitted[8])]
            if test2 in dict_annotation_S_all.keys():
              x=dict_annotation_S_all[test2]
              x.append(float(line_splitted[8]))
            dict_annotation_S_all[test2]=x
    
    print("annotok")
    
    fichier=open(str(ref)+"_conservation_by_gene_all.csv", "w")
    fichier.write("Gene;Chrom;Start;End;Divergence_D sister;Divergence_DNS sister;Divergence_DS sister;Number_of_seq_considered;Divergence_D all;Divergence_DNS all;Divergence_DS all")
    
    for i in Bed:
      seq=i.split("\t")
      Gene_name=str(seq[3])
      Chrom=str(seq[0])
      Start=int(seq[1])
      End=int(seq[2])
      fichier.write("\n"+str(Gene_name)+";"+str(Chrom)+";"+str(Start)+";"+str(End))
      
      #sister
      D="NA"
      DNS="NA"
      DS="NA"
      if Gene_name in dict_annotation2.keys() :
        x=dict_annotation2[Gene_name]
        if len(x)>0:
          D=sum(x)/len(x)
      if Gene_name in dict_annotation_NS.keys() :
        x=dict_annotation_NS[Gene_name]
        if len(x)>0:
          DNS=sum(x)/len(x)
      if Gene_name in dict_annotation_S.keys() :
        x=dict_annotation_S[Gene_name]
        if len(x)>0:
          DS=sum(x)/len(x)
      fichier.write(";"+str(D)+";"+str(DNS)+";"+str(DS))
      
      # all
      nb=0
      if Gene_name in dict_annotation_gene.keys() :
        x=dict_annotation_gene[Gene_name]
        if x[1]>0:
          nb=x[1]
      fichier.write(";"+str(nb))
      D="NA"
      DNS="NA"
      DS="NA"
      if Gene_name in dict_annotation_all.keys() :
        x=dict_annotation_all[Gene_name]
        #print(x)
        if len(x)>0:
          D=sum(x)/len(x)
      if Gene_name in dict_annotation_NS_all.keys() :
        x=dict_annotation_NS_all[Gene_name]
        #print(x)
        if len(x)>0:
          DNS=sum(x)/len(x)
      if Gene_name in dict_annotation_S_all.keys() :
        x=dict_annotation_S_all[Gene_name]
        #print(x)
        if len(x)>0:
          DS=sum(x)/len(x)
      fichier.write(";"+str(D)+";"+str(DNS)+";"+str(DS))

      print(Gene_name)
    fichier.close()

## presence
      
def presence_estimator(ref,Bed,REF):

    fichier2=open(str(ref)+"_presence_by_gene.csv", "w")
    fichier2.write("Gene;Chrom;Start;End")
    for species in REF:
        ref2=species.split("\t")
        if str(ref2[0])!= ref:fichier2.write(";"+str(ref2[0]))
    
    for i in Bed:
      seq=i.split("\t")
      Gene_name=str(seq[3])
      Chrom=str(seq[0])
      Start=int(seq[1])
      End=int(seq[2])
      name_ref=str(Chrom)+":"+str(Start)+"-"+str(End)
      fichier2.write("\n"+str(Gene_name)+";"+str(Chrom)+";"+str(Start)+";"+str(End))
      for species in REF:
        presence="0"
        ref2=species.split("\t")
        if str(ref2[0])!= ref:
          file_temp=open(str(ref2[0])+".bed","r")
          file_temp=file_temp.read()
          if file_temp.count(name_ref)!=0:presence="1"
          fichier2.write(";"+str(presence))
      
    fichier2.close()    
    
## create phyml file

def phyml_preparation(ref,Bed,REF):
    
    for candidate in Bed:
      seq=candidate.split("\t")
      x2=str(seq[3])+"_alig.fa"
      if os.path.exists(str(x2)):
        fasta_alig_file = open(str(x2), "r")
        fasta_alig_file=fasta_alig_file.read()
        Gene_name=str(seq[3])
        fasta_alig_file = fasta_alig_file.split(">")
        fichier=open(str(seq[3])+"_alig_ordered_phyml_ready.fas", "w") # a supprimer a la fin
        nb_seq=len(fasta_alig_file)-1
        seq2=fasta_alig_file[1]
        seq3=seq2.split("\n")
        seq3=seq3[0] #name
        seq2=seq2.replace("\n","")
        seq2=seq2.replace(seq3,"") #seq real
        fichier.write(str(nb_seq)+" "+str(len(seq2)))

        for seq in fasta_alig_file[1:]:
          seq2=seq.split("\n")
          seq3=seq2[0] #name
          seq2=str(seq)
          seq2=seq2.replace("\n","")
          seq2=seq2.replace(seq3,"") #seq real
          fichier.write("\n"+str(seq3)+" "+str(seq2))

        fichier.close()

    dict_species={}
    nb_gene=0
    for Gene in Bed:
      gene=Gene.split("\t")
      #if gene[4]=="Leaf" and gene[4]=="Root":
      if gene[4]!="Control" :
        x2=str(gene[3])+"_alig_ordered_phyml_ready.fas"
        if os.path.exists(str(x2)):
          fasta_alig_file = open(str(x2), "r")
          fasta_alig_file=fasta_alig_file.read()
          fasta_alig_file2 = fasta_alig_file.split("\n")
          nb_seq=fasta_alig_file2[0]
          nb_seq=nb_seq.split(" ")
          false_seq="-"*int(nb_seq[1])
          nb_seq=int(nb_seq[0])

          if nb_seq>=3:
            nb_gene=nb_gene+1
            i=1
            for species in REF:
              name=species.split("\t")
              name_species=str(name[0])
              seq=""
              if fasta_alig_file.count(name_species)==0:
                if name_species in dict_species.keys():seq=str(dict_species[name_species])
                seq=str(seq)+str(false_seq)
                dict_species[name_species]=seq
            for species in fasta_alig_file2[1:]:
              name=species.split(" ")
              name_species=str(name[0])
              seq=""
              if name_species in dict_species.keys():seq=str(dict_species[name_species])
              seq2=name[1]
              seq=str(seq)+str(seq2)
              dict_species[name_species]=seq
    print(nb_gene)
    
    fichier=open("all_candidate_alig_ordered_phyml_ready.fas", "w")
    nb_seq=len(REF)
    species=REF[0]
    name=species.split("\t")
    name_species=str(name[0])
    seq2=str(dict_species[name_species])
    fichier.write(str(nb_seq)+" "+str(len(seq2)))
    for species in REF:
      name=species.split("\t")
      name_species=str(name[0])
      seq2=str(dict_species[name_species])
      fichier.write("\n"+str(name_species)+" "+str(seq2))

    fichier.close()


    
## phylogenetic three of align fasta of genes in all species
def phyml_all(ref,Bed,phyml):
    
    x2="all_candidate_alig_ordered_phyml_ready.fas"
    cmd="phyml -i "+str(x2)+" -b 1000 -o n > output_"+str(x2)
    os.system(cmd)
    
    cmd="rm "+str(x2)+"_phyml_boot_trees.txt"
    os.system(cmd) 
    cmd="rm "+str(x2)+"_phyml_boot_stats.txt"
    os.system(cmd) 
    cmd="rm output_"+str(x2)
    os.system(cmd) 
    
    x2="all_candidate_alig_ordered_phyml_ready.fas_phyml_tree.txt"
    fasta_alig_file = open(str(x2), "r")
    fasta_alig_file=fasta_alig_file.read()
    fichier=open("all_candidate_alig_ordered_phyml_topology.nwk", "w")
    fasta_alig_file=str(fasta_alig_file)
    j=0
    for i in fasta_alig_file:
      if str(i)!="0" and str(i)!="1" and str(i)!="2" and str(i)!="3" and str(i)!="4" and str(i)!="5" and str(i)!="6" and str(i)!="7" and str(i)!="9" and str(i)!="8" and str(i)!=":" and str(i)!="." : fichier.write(str(i))
      elif str(i)=="." and str(fasta_alig_file[j-1])!="0"  : fichier.write(str(i))
      j=j+1
    fichier.close()


## phylogenetic tree of each corrected align fasta
    
def phyml(ref,Bed,phyml):    
    for i in Bed:
      seq=i.split("\t")
      x2=str(seq[3])+"_alig_ordered_phyml_ready.fas"
      if os.path.exists(str(x2)):
        cmd="phyml -i "+str(x2)+" -b 1000 -o n> output_"+str(x2)
        os.system(cmd) 
    
    for i in Bed:
      seq=i.split("\t")
      x2=str(seq[3])+"_alig_ordered_phyml_ready.fas"    
      cmd="rm "+str(x2)+"_phyml_boot_trees.txt"
      os.system(cmd) 
      cmd="rm "+str(x2)+"_phyml_boot_stats.txt"
      os.system(cmd) 
 
      cmd="cp "+str(x2)+"_phyml_tree.txt "+str(x2)+"_save_alig_ordered_phyml_ready.nwk"
      os.system(cmd) 
      cmd="rm output_"+str(x2)
      os.system(cmd) 
    print("ok phyml")




## length of branch

def phyml_length_branch(REF,Bed,phyml):    
    for candidate in Bed:
      seq=candidate.split("\t")
      x2=str(seq[3])+"_alig_ordered_phyml_ready.fas"
      if os.path.exists(str(x2)):
        fasta_alig_file = open(str(x2), "r")
        fasta_alig_file=fasta_alig_file.read()
        fichier=open(str(seq[3])+"_alig_ordered_phyml_ready_all_species.fas", "w")
        fasta_alig_file2 = fasta_alig_file.split("\n")
        seq2=fasta_alig_file2[1]
        seq3=seq2.split(" ")
        seq2=seq3[1] #name
        len_seq=len(seq2)
        false_seq="-"*len(seq2)
        fichier.write(str(len(REF))+"\t"+str(len_seq))
        for seq in fasta_alig_file2[1:]:fichier.write("\n"+str(seq))

        for species in REF:
          nb_gene=0
          name=species.split("\t")
          name_species=str(name[0])
          if fasta_alig_file.count(str(name_species))==0: 
            fichier.write("\n"+str(name_species)+"\t"+str(false_seq))
        
        fichier.close()

    for i in Bed:
      seq=i.split("\t")
      x2=str(seq[3])+"_alig_ordered_phyml_ready_all_species.fas" 
      cmd="phyml -i "+str(seq[3])+"_alig_ordered_phyml_ready_all_species.fas -u all_candidate_alig_ordered_phyml_topology.nwk -b 0 -o n> "+str(seq[3])+"_phylo"
      os.system(cmd)
      cmd="rm "+str(seq[3])+"_phylo"
      cmd="cp "+str(x2)+"_phyml_tree.txt "+str(x2)+"_save_all_phylo_phyml_tree_all_species.nwk"
      os.system(cmd)  
    print("ok phyml_length_branch")
          


## Basic stat

def basic_stat(ref,gff,bedtools,REF,directory,Bed):
    
    cds=str(ref)+"_cds.bed"
    fichier=open(str(cds), "w")
    x="#Chrom\tstart\tend"
    fichier.write(str(x))
    GFF=open(str(directory)+"Database/"+str(gff), "r")
    GFF=GFF.read()
    GFF=GFF.split('\n')
    for i in GFF:
      if i.count("#")==0 and i!="":
        seq2=i.split("\t")
        if len(seq2)>3 :
          if str(seq2[2])=="CDS"  :
            x=str(seq2[0])+"\t"+str(seq2[3])+"\t"+str(seq2[4])
            fichier.write("\n"+str(x) )   
    fichier.close()

    x2=str(ref)+"_basic_stat.csv"
    fichier=open(str(x2), "w")
    fichier.write("gene;chrom;start;end;group;gene_size;cds_percentage;GC_percentage;GC_tot in CDS; size CDS; GC_percent_cds")
    
    x2=str(directory)+str(ref)+"_study_filtrered.fa"
    Fasta=open(str(x2),"r")
    Fasta=Fasta.read()
    Fasta=Fasta.split('>')
    Fasta=Fasta[1:]
    tot_gene_fasta={}
    for seq in Fasta:
      seq2=seq.split("\n")
      name=seq2[0]
      seq2=seq2[1]
      tot_gene_fasta[name]=seq2

    for species in REF:
        ref2=species.split("\t")
        if str(ref2[0])== ref:
            cmd=str(bedtools)+" getfasta -fi "+str(directory)+"Database/"+str(ref2[1])+" -bed "+str(cds)+" > "+str(ref)+"_CDS_basic_analysis.fa"
            os.system(cmd)

    x2=str(ref)+"_CDS_basic_analysis.fa"
    Fasta_cds=open(str(x2),"r")
    Fasta_cds=Fasta_cds.read()
    Fasta_cds=Fasta_cds.split('>')
    Fasta_cds=Fasta_cds[1:]
    CDS_fasta={}
    for seq in Fasta_cds:
      seq2=seq.split("\n")
      name=seq2[0]
      seq2=seq2[1]
      CDS_fasta[name]=seq2


    GFF=open(str(directory)+"Database/"+str(gff),"r")
    GFF=GFF.read()
    GFF=GFF.split('\n')
    GFF=GFF[3:]

    DIC_gene={}
    DIC_cds={}
    for gene in GFF:
      Gene=gene.split("\t")
      if gene!="":
        name=str(Gene[8])
        name=name.split(".")
        care=name[1]
        name=str(name[0])
        name=str(name[3:])
        if str(Gene[2])=="gene" :
          size=int(Gene[4])-int(Gene[3])
          DIC_gene[str(name)]=[size,0]

        if str(Gene[2])=="CDS" and care=="t1":
          if name in DIC_gene.keys():
            x=DIC_gene[str(name)]
            size=int(x[0])
            CDS=int(Gene[4])-int(Gene[3])
            CDS=int(x[1])+int(CDS)
            DIC_gene[str(name)]=[size,CDS] # taille cds cummulée
            start=int(Gene[3])
            end=int(Gene[4])
            chrom=str(Gene[0])
            CDS_name=str(chrom)+":"+str(start)+"-"+str(end)
            if name not in DIC_cds.keys():DIC_cds[str(name)]=[CDS_name]
            elif name in DIC_cds.keys():
              x=DIC_cds[str(name)]
              x.append(CDS_name)
              DIC_cds[str(name)]=x


    i=0
    for gene in Bed:
      Gene=gene.split("\t")
      start=int(Gene[1])
      end=int(Gene[2])
      chrom=str(Gene[0])
      ID=str(chrom)+":"+str(start)+"-"+str(end)
      group=str(Gene[4])
      name=str(Gene[3])
      fasta2=str(tot_gene_fasta[str(ID)])
      GC=float((fasta2.count("G")+fasta2.count("C"))/len(fasta2))
      GC=GC*100
      if name in DIC_gene.keys():
        x=DIC_gene[str(name)]
        CDS=int(x[1])*100
        size=int(x[0])
      if name not in DIC_gene.keys():
        size=len(fasta2)
        CDS="NA"
      
      
      CDS_seq=DIC_cds[str(name)]
      GC2=0
      size2=0
      for j in CDS_seq:
        fasta2=str(CDS_fasta[j])
        GC2=GC2+int((fasta2.count("G")+fasta2.count("C")))
        size2=size2+len(fasta2)

      following_line="\n"+str(name)+";"+str(chrom)+";"+str(start)+";"+str(end)+";"+str(group)+";"+str(size)+";"+str(CDS/size)+";"+str(GC)+";"+str(GC2)+";"+str(size2)+";"+str(GC2*100.0/size2)
      fichier.write(str(following_line))

    fichier.close()



