import csv
import sys
from operator import itemgetter
##########################################################################################
#######################      Translation  Dictionary     #################################
##########################################################################################

translate = {'ATG' : 'M', 'ATA' : 'I', 'ATC' : 'I', 'ATT': 'I', 'ACG' : 'T', 'ACA' : 'T', 'ACC' : 'T',
 'ACT' : 'T', 'AAG' : 'K', 'AAA' : 'K', 'AAC' : 'N', 'AAT' : 'N', 'AGG' : 'R', 'AGA' : 'R', 'AGC' : 'S',
 'AGT' : 'S', 'GTG' : 'V', 'GTA' : 'V', 'GTC' : 'V', 'GTT' : 'V', 'GCG' : 'A', 'GCA' : 'A', 'GCC' : 'A',
 'GCT' : 'A', 'GAG' : 'E', 'GAA' : 'E', 'GAC' : 'D', 'GAT' : 'D', 'GGG' : 'G', 'GGA' : 'G', 'GGC' : 'G',
 'GGT' : 'G', 'TTG' : 'L', 'TTA' : 'L', 'TTC' : 'F', 'TTT' : 'F', 'TCG' : 'S', 'TCA' : 'S', 'TCC' : 'S',
 'TCT' : 'S', 'TAG' : 'x', 'TAA' : 'x', 'TAC' : 'Y', 'TAT' : 'Y', 'TGG' : 'W', 'TGA' : 'x', 'TGC' : 'C',
 'TGT' : 'C', 'CTG' : 'L', 'CTA' : 'L', 'CTC' : 'L', 'CTT' : 'L', 'CCG' : 'P', 'CCA' : 'P', 'CCC' : 'P',
 'CCT' : 'P', 'CAG' : 'Q', 'CAA' : 'Q', 'CAC' : 'H', 'CAT' : 'H', 'CGG' : 'R', 'CGA' : 'R', 'CGC' : 'R',
 'CGT' : 'R'}
 
 
##########################################################################################
#######################       Codon count Dictionary     #################################
##########################################################################################

codon_count = {'ATG' : 0, 'ATA' : 0, 'ATC' : 0, 'ATT': 0, 'ACG' : 0, 'ACA' : 0, 'ACC' : 0,
 'ACT' : 0, 'AAG' : 0, 'AAA' : 0, 'AAC' : 0, 'AAT' : 0, 'AGG' : 0, 'AGA' : 0, 'AGC' : 0,
 'AGT' : 0, 'GTG' : 0, 'GTA' : 0, 'GTC' : 0, 'GTT' : 0, 'GCG' : 0, 'GCA' : 0, 'GCC' : 0,
 'GCT' : 0, 'GAG' : 0, 'GAA' : 0, 'GAC' : 0, 'GAT' : 0, 'GGG' : 0, 'GGA' : 0, 'GGC' : 0,
 'GGT' : 0, 'TTG' : 0, 'TTA' : 0, 'TTC' : 0, 'TTT' : 0, 'TCG' : 0, 'TCA' : 0, 'TCC' : 0,
 'TCT' : 0, 'TAG' : 0, 'TAA' : 0, 'TAC' : 0, 'TAT' : 0, 'TGG' : 0, 'TGA' : 0, 'TGC' : 0,
 'TGT' : 0, 'CTG' : 0, 'CTA' : 0, 'CTC' : 0, 'CTT' : 0, 'CCG' : 0, 'CCA' : 0, 'CCC' : 0,
 'CCT' : 0, 'CAG' : 0, 'CAA' : 0, 'CAC' : 0, 'CAT' : 0, 'CGG' : 0, 'CGA' : 0, 'CGC' : 0,
 'CGT' : 0}
 
 
##########################################################################################
#######################             File Names            ################################
##########################################################################################
 
fwd_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/hht_try1_plus_edit.wig"
rev_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/hht_try1_minus_edit.wig"
annotation_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/GCF_000025685.1_ASM2568v1_genomic_update_080119_main_chrom_only_for_script.txt"
print annotation_file
seq_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/GCF_000025685.1_ASM2568v1_genomic_update_080119_main_chrom_only_single_line_edit_for_riboRET_script.txt"
#f = open (input_file, 'r')


##########################################################################################
#####################       Import the wig to memory        ##############################
##########################################################################################

Gene_starts = []
Gene_ends = []
strands = []
genes = []
gene_density={}
gene_counter=0



e=open(seq_file,'r')
with open(annotation_file,'rU') as d:
	genes=zip(*[line.split('\t') for line in d])[0]

with open(annotation_file,'rU') as d:
    Gene_starts=(zip(*[line.split('\t') for line in d])[1])
    
with open(annotation_file,'rU') as d:
	Gene_ends=(zip(*[line.split('\t') for line in d])[2])

with open(annotation_file,'rU') as d:
	strands=(zip(*[line.split('\t') for line in d])[3])

reader=csv.reader(open(annotation_file), delimiter="\t")

g=open('Harringtonine_Hfx_all_gene_density.txt','w')
##########################################################################################
#Open the forward file with the header previously removed. Then for each line in the file#
#check to see if the peak falls within an annotated gene of the correct direction. If it #
#does then add that density to that gene. Repeat this for the forward and reverse 		 #
#directions. Once this is done write the newly created dictionary to a file.             #
##########################################################################################
def density_countfwd():
	with open(fwd_file,'rU') as f:
		start=0
		end=0
		gene_name=''
		frame=4
		start_position=0
		direction="+"
		count=1
		global gene_counter
		for line in f:
			print line
			count+=1
			fields=line.split()
			position=int(fields[0])
			rpm=float(fields[1])
			if position>=start and position<=end:
				gene_density[gene_name]+=float(rpm);
				print position
			else:
				for i in range(len(genes)):
					if position >= int(Gene_starts[i]) and position <= int(Gene_ends[i]) and strands[i]=="+":
						start=int(Gene_starts[i])
						end=int(Gene_ends[i])
						gene_name=genes[i]
						gene_density[genes[i]]=float(rpm)
						gene_counter+=1;
					else:
						pass;
						
def density_countrev():
	with open(rev_file,'rU') as f:
		start=0
		end=0
		gene_name=''
		frame=4
		start_position=0
		direction="-"
		global gene_counter
		for line in f:
			print line
			fields=line.split()
			position=int(fields[0])
			rpm=float(fields[1])
			if position<=start and position>=end:
				gene_density[gene_name]+=float(rpm);
				print position
			else:
				for i in range(len(genes)):
					if position <= int(Gene_starts[i]) and position >= int(Gene_ends[i]) and strands[i]=="-":
						start=int(Gene_starts[i])
						end=int(Gene_ends[i])
						gene_name=genes[i]
						gene_density[genes[i]]=float(rpm)
						gene_counter+=1;
					else:
						pass;
density_countfwd()
density_countrev()
#print gene_density

for i in gene_density:
	g.write(i+'\t'+str(gene_density[i])+'\n')
#print gene_density