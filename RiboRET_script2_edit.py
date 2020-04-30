import csv
import sys
from operator import itemgetter
from Bio.Seq import Seq
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
 
reverse_complement = {'CAT' : 'M', 'TAT' : 'I', 'GAT' : 'I', 'AAT': 'I', 'CGT' : 'T', 'TGT' : 'T', 'GGT' : 'T',
 'AGT' : 'T', 'CTT' : 'K', 'TTT' : 'K', 'GTT' : 'N', 'ATT' : 'N', 'CCT' : 'R', 'TCT' : 'R', 'GCT' : 'S',
 'ACT' : 'S', 'CAC' : 'V', 'TAC' : 'V', 'GAC' : 'V', 'AAC' : 'V', 'CGC' : 'A', 'TGC' : 'A', 'GGC' : 'A',
 'AGC' : 'A', 'CTC' : 'E', 'TTC' : 'E', 'GTC' : 'D', 'ATC' : 'D', 'CCC' : 'G', 'TCC' : 'G', 'GCC' : 'G',
 'ACC' : 'G', 'CAA' : 'L', 'TAA' : 'L', 'GAA' : 'F', 'AAA' : 'F', 'CGA' : 'S', 'TGA' : 'S', 'GGA' : 'S',
 'AGA' : 'S', 'CTA' : 'x', 'TTA' : 'x', 'GTA' : 'Y', 'ATA' : 'Y', 'CCA' : 'W', 'TCA' : 'x', 'GCA' : 'C',
 'ACA' : 'C', 'CAG' : 'L', 'TAG' : 'L', 'GAG' : 'L', 'AAG' : 'L', 'CGG' : 'P', 'TGG' : 'P', 'GGG' : 'P',
 'AGG' : 'P', 'CTG' : 'Q', 'TTG' : 'Q', 'GTG' : 'H', 'ATG' : 'H', 'CCG' : 'R', 'TCG' : 'R', 'GCG' : 'R',
 'ACG' : 'R'}
 
 
##########################################################################################
#######################             File Names            ################################
##########################################################################################
 
fwd_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/hht_try1_plus_edit.wig"
rev_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/hht_try1_minus_edit.wig"
annotation_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/GCF_000025685.1_ASM2568v1_genomic_update_080119_main_chrom_only_for_script.txt"
seq_file="/Users/DRG/Downloads/2019-smORF-master/all_TSS/GCF_000025685.1_ASM2568v1_genomic_update_080119_main_chrom_only_single_line_edit_for_riboRET_script.txt"
#f = open (input_file, 'r')


##########################################################################################
##################       Import the annotation to memory        ##########################
##########################################################################################

Gene_starts = []
Gene_ends = []
strands = []
genes = []


e=open(seq_file,'r')
gene_density={}
with open(annotation_file,'rU') as d:
	genes=zip(*[line.split('\t') for line in d])[0]

with open(annotation_file,'rU') as d:
    Gene_starts=(zip(*[line.split('\t') for line in d])[1])
    
with open(annotation_file,'rU') as d:
	Gene_ends=(zip(*[line.split('\t') for line in d])[2])

with open(annotation_file,'rU') as d:
	strands=(zip(*[line.split('\t') for line in d])[3])


with open('Harringtonine_Hfx_all_gene_density.txt', 'rU') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        gene_density[row[0]] = row[1]
        #print row[0],
        #print row[1]

reader=csv.reader(open(annotation_file), delimiter="\t")

##########################################################################################
##################       Write header of the output file        ##########################
##########################################################################################

g=open('Harringtonine_Hfx_annotation_output_all.txt','w')
g.write('type')
g.write("\t")
g.write('start_position')
g.write("\t")
g.write('stop_position')
g.write('\t')
g.write('strand')
g.write('\t')
g.write('gene')
g.write('\t')
g.write('aa_length')
g.write("\t")
g.write('peak_height')
g.write("\t")
g.write('codon')
g.write("\t")
g.write('frame')
g.write("\t")
g.write('aa_seq')
g.write("\t")
g.write('nt_seq')
g.write("\t")
g.write('relative_density')
g.write("\t")
g.write('5\'_Distance')
g.write("\t")
g.write('3\'_Distance')
g.write("\n")
#print gene_density

##########################################################################################
####################       Counter for naming new genes        ###########################
##########################################################################################
new_gene_counter=1


##########################################################################################
#################       Function for reading forward file        #########################
##########################################################################################               
def readfileforward():
	new_gene_counter=1
	
##########################################################################################
#Open the sequence file, 'n' checks for the first start site, 'b' is a switch that 		 #
#controls checking for rRNA genes. The seq_chunk is the sequence around a peak in which  #
#we look for a start codon. Next we open the wig file as .txt. The wig header has already#
#been removed. For each line we first blank the gene, frame, and start position. We read #
#the line and separate the position and rpm.											 #
##########################################################################################	
	
	with open(seq_file,'rU') as d:
		n=1
		b=0
		last_rpm=0
		seq_chunk=''
		f=open(fwd_file,'rU')
		for line in f:
			gene=''
			frame=4
			start_position=0
			direction="+"
			fields=line.split()
			position=int(fields[0])
			rpm=float(fields[1])
			
##########################################################################################
#If the rpm is greater than 5 we check if there is a start codon near that position. We  #
#grab a chunk from -3 to +3 around that peak and look for a start codon. The start 		 #
#position in that chunk is adjusted to map to the genome start position. The codon is	 #
#designated																				 #
##########################################################################################
			if rpm>5:
				#print "Forward check"
				codon='DDD'
				e.seek(position-4)
				seq_chunk=e.read(7)
				temp_position= -1
				if	seq_chunk.find('ATG')!=-1:
					temp_position=seq_chunk.find('ATG')
					codon='ATG'
					start_position=position+temp_position-3;
				elif	seq_chunk.find('GTG')!=-1:
					temp_position=seq_chunk.find('GTG')
					codon='GTG'
					start_position=position+temp_position-3;
				elif	seq_chunk.find('TTG')!=-1:
					temp_position=seq_chunk.find('TTG')
					codon='TTG'
					start_position=position+temp_position-3;
				else:
					start_position=position
					codon='DDD';	
##########################################################################################
#Here we ignore rrl, rrs, or rra genes. We also check for the first peak in the genome.  #
#If the last start position is the same as the current start position we sum the peaks.  #
#If they are different we prepare to write them to the file								 #
##########################################################################################								
				
				if b==1:
					b=0;
				elif 'rrl' in gene or 'rrs' in gene or 'rra' in gene:
					b=0;	
				elif n==1:
					last_position=position
					last_rpm=rpm
					last_start_position=start_position
					n+=1;		
				elif last_start_position==start_position:
					last_rpm+=rpm;
##########################################################################################
#Read the first codon, and translate it. Add that codon/aa to the mRNA and aa sequences. #
#First check if the start position matches any annotated sites or if it falls within 10  #
#nt of an annotated site. If the gene shares the same stop codon of an annotated gene 	 #
#check if it is either an internal initiation site or an N-terminal extension. If none of#
#these tests are positive then check if it initiation is within the region of another 	 #
#gene but out of frame. If the initiation site is external to annotated genes and doesn't#
#share the same stop codon then it is an unannotated gene.								 #
##########################################################################################
				else:
					d.seek(last_start_position-1)
					nt=d.read(3)
					aa=translate[nt]
					start_codon=nt
					if nt=='ATG' or nt=='GTG' or nt=='TTG': 
						aa='M'
					nt_seq=nt
					aa_seq=aa
					while aa != 'x':
						nt=d.read(3)
						aa=translate[nt]
						nt_seq+=nt
						aa_seq+=aa
					stop_position=d.tell()
					frame=0	
					type='Unannotated'
					for i in range(len(genes)):
						if last_start_position==int(Gene_starts[i]) and strands[i]=="+":
							type='Annotated'
							frame=0
							gene=genes[i]
							break
						elif abs(last_start_position-int(Gene_starts[i]))<10 and abs(int(Gene_starts[i])-last_start_position)<10 and strands[i]=='+':
							type='Near_Annotated'
							gene=genes[i]
							if 0== (int(Gene_starts[i])-last_start_position)%3:
								frame=0;
							elif 1== (int(Gene_starts[i])-last_start_position)%3:
								frame=2;
							else:
								frame=1;
							break
						elif stop_position==int(Gene_ends[i]) and strands[i]=="+":
							if last_start_position > int(Gene_starts[i]):
								type='Internal_Inframe'
								frame=0
								gene=genes[i]
								break
							else:
								type='N-terminal_extension'
								gene=genes[i]
								frame=0
								break
						else:
							if last_start_position >= int(Gene_starts[i]) and last_start_position <= int(Gene_ends[i]) and strands[i]=='+':
								type='Internal_OutofFrame'
								gene=genes[i]
								if 0== (int(Gene_starts[i])-last_position)%3:
									frame=0;
								elif 1== (int(Gene_starts[i])-last_position)%3:
									frame=2;
								else:
									frame=1;
								break;

					if type=='Unannotated':
						gene="gene_F_"+str(new_gene_counter)
						new_gene_counter+=1
##########################################################################################
#Write the information. If it initiation is annotated or internal then calculate the 	 #
#relative densities and distance from the 3' and 5' end. If it is an N-terminal extension#
#then only count the distance but not relative density.									 #
##########################################################################################
					g.write(type)
					g.write("\t")
					g.write(str(last_start_position))
					g.write("\t")
					g.write(str(stop_position))
					g.write("\t")
					g.write('+')
					g.write("\t")
					g.write(gene)
					g.write("\t")
					g.write(str((stop_position-last_start_position-2)/3))
					g.write("\t")
					g.write(str(last_rpm))
					g.write("\t")
					g.write(start_codon)
					g.write('\t')
					g.write(str(frame))
					g.write('\t')
					g.write(aa_seq)
					g.write("\t")
					g.write(nt_seq)
					g.write("\t")
					if type=='Annotated' or type=="Internal_OutofFrame":
						g.write(str(last_rpm/float(gene_density[gene])))
						g.write("\t")
						g.write(str(last_start_position-int(Gene_starts[genes.index(gene)])))
						g.write("\t")
						g.write(str(int(Gene_ends[genes.index(gene)])-last_start_position+1))
					elif type=='Internal_Inframe':
						g.write(str(last_rpm/float(gene_density[gene])))
						g.write("\t")
						g.write(str(last_start_position-int(Gene_starts[genes.index(gene)])))
						g.write("\t")
						g.write(str(int(Gene_ends[genes.index(gene)])-last_start_position+1))
					elif type=='N-terminal_extension':
						g.write('N/A')
						g.write("\t")
						g.write(str(last_start_position-int(Gene_starts[genes.index(gene)])))
						g.write("\t")
						g.write(str(int(Gene_ends[genes.index(gene)])-last_start_position+1))						
					else:
						g.write("N/A")
						g.write("\t")
						g.write("N/A")
						g.write("\t")
						g.write("N/A")
					g.write("\n")
					last_start_position=start_position
					last_rpm=rpm;

#######################################################################################
#######################################################################################
#######################################################################################

									
def readfilereverse():
	new_gene_counter=1
	with open(seq_file,'rU') as d:
		n=1
		b=0
		seq_chunk=''
		l=open(rev_file,'rU')
		for line in l:
			gene=''
			frame=4
			start_position=0
			direction="-"
			fields=line.split()
			position=int(fields[0])
			rpm=float(fields[1])
			if rpm>5:
				codon='DDD'
				e.seek(position-4)
				seq_chunk=e.read(7)
				temp_position= -1
				if	seq_chunk.find('CAT')!=-1:
					temp_position=seq_chunk.find('CAT')-1
					codon='ATG'
					start_position=position+temp_position;
				elif	seq_chunk.find('CAC')!=-1:
					temp_position=seq_chunk.find('CAC')-1
					codon='GTG'
					start_position=position+temp_position;
				elif	seq_chunk.find('CAA')!=-1:
					temp_position=seq_chunk.find('CAA')-1
					codon='TTG'
					start_position=position+temp_position;
				elif	seq_chunk.find('CAG')!=-1:
					temp_position=seq_chunk.find('CAG')-1
					codon='CTG'
					start_position=position+temp_position;
				elif	seq_chunk.find('AAT')!=-1:
					temp_position=seq_chunk.find('AAT')-1
					codon='ATT'
					start_position=position+temp_position;
				elif	seq_chunk.find('GAT')!=-1:
					temp_position=seq_chunk.find('GAT')-1
					codon='ATC'
					start_position=position+temp_position;

				else:
					start_position=position
					codon='DDD';	
				
				
				if b==1:
					b=0;
				elif 'rrl' in gene or 'rrs' in gene or 'rra' in gene:
					b=0;	
				elif n==1:
					last_position=position
					last_rpm=rpm
					last_start_position=start_position
					n+=1;		
				elif last_start_position==start_position:
					last_rpm+=rpm;
					#print rpm
				else:
					d.seek(last_start_position-3)
					nt=str(Seq(d.read(3)).reverse_complement())
					aa=translate[nt]
					start_codon=nt
					if nt=='ATG' or nt=='GTG' or nt=='TTG': 
						aa='M'
					nt_seq=nt
					aa_seq=aa
					while aa != 'x':
						d.seek(-6,1)
						nt=str(Seq(d.read(3)).reverse_complement())
						aa=translate[nt]
						nt_seq+=nt
						aa_seq+=aa
					if last_start_position==4389135:
						print codon	
					stop_position=d.tell()-2
					frame=0	
					#nt_seq=Seq(nt_seq).reverse_complement()
					type='Unannotated'
					for i in range(len(genes)):
						if last_start_position==int(Gene_starts[i]) and strands[i]=="-":
							type='Annotated'
							frame=0
							gene=genes[i]
							break
						elif abs(last_start_position-int(Gene_starts[i]))<10 and abs(int(Gene_starts[i])-last_start_position)<10 and strands[i]=='-':
							type='Near_Annotated'
							gene=genes[i]
							if 0== (int(Gene_starts[i])-last_start_position)%3:
								frame=0;
							elif 1== (int(Gene_starts[i])-last_start_position)%3:
								frame=2;
							else:
								frame=1;
							break
						elif stop_position==int(Gene_ends[i]) and strands[i]=="-":
							if last_start_position < int(Gene_starts[i]):
								type='Internal_Inframe'
								frame=0
								gene=genes[i]
								break
							else:
								type='N-terminal_extension'
								gene=genes[i]
								frame=0
								break
						else:
							if last_start_position <= int(Gene_starts[i]) and last_start_position >= int(Gene_ends[i]) and strands[i]=='-':
								type='Internal_OutofFrame'
								gene=genes[i]
								if 0== (int(Gene_starts[i])-last_position)%3:
									frame=0;
								elif 1== (int(Gene_starts[i])-last_position)%3:
									frame=1;
								else:
									frame=2;
								break;

					if type=='Unannotated':
						gene="gene_R_"+str(new_gene_counter)
						new_gene_counter+=1
					g.write(type)
					g.write("\t")
					g.write(str(last_start_position))
					g.write("\t")
					g.write(str(stop_position))
					g.write("\t")
					g.write('-')
					g.write("\t")
					g.write(gene)
					g.write("\t")
					g.write(str((last_start_position-stop_position-2)/3))
					g.write("\t")
					g.write(str(last_rpm))
					g.write("\t")
					g.write(start_codon)
					g.write('\t')
					g.write(str(frame))
					g.write('\t')
					g.write(aa_seq)
					g.write("\t")
					g.write(str(nt_seq))
					g.write("\t")
					if type=='Annotated' or type=="Internal_OutofFrame":
						g.write(str(last_rpm/float(gene_density[gene])))
						g.write("\t")
						g.write(str(int(Gene_starts[genes.index(gene)])-last_start_position))
						g.write("\t")
						g.write(str(last_start_position-int(Gene_ends[genes.index(gene)])+1))
					elif type=='Internal_Inframe':
						g.write(str(last_rpm/float(gene_density[gene])))
						g.write("\t")
						g.write(str(int(Gene_starts[genes.index(gene)])-last_start_position))
						g.write("\t")
						g.write(str(last_start_position-int(Gene_ends[genes.index(gene)])+1))
					elif type=="N-terminal_extension":
						g.write('N/A')
						g.write("\t")
						g.write(str(int(Gene_starts[genes.index(gene)])-last_start_position))
						g.write("\t")
						g.write(str(last_start_position-int(Gene_ends[genes.index(gene)])+1))
					else:
						g.write("N/A")
						g.write("\t")
						g.write("N/A")
						g.write("\t")
						g.write("N/A")
					g.write("\n")
					last_start_position=start_position
					last_rpm=rpm
					
    		
readfileforward()
readfilereverse()

