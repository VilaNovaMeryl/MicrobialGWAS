####################################################
#####  Converting a VCF file to a binary file  #####
####################################################

def construct_file(fileNameIn, fileNameOut): 

	with open (fileNameIn,"r") as inputfile:  # open the input file in read mode

		with open (fileNameOut,"w") as outputfile: # opening the output file in writing mode

			total = inputfile.readlines() # read and save file lines 

			for line in total:

				if line.startswith('##'): # description line of a VCF file
					pass

				elif line.startswith('#CHROM'): # line with the names of the different genomes analyzes
					firstline = line
					firstlines = firstline.replace('Resultat/BAM/','') # remove the path from the genome names
					firstliness = firstlines.split() # retrieve the names of the different genomes in a list
					
					comment1 = 'Variant' + '\t' + '\t'.join(firstliness[9:]) + '\n' # create the first line of the binary file containing the names of the different genomes at the head of the column
					
					outputfile.write(comment1)

				else :
					linesplit = line.split("\t") # lines with the genotype reference / alternative information of the different genomes
					chromosome = linesplit[0] # name of the chromosome on which the variant is located
					position = linesplit[1] # position of the variant on the reference genome
					ref = linesplit[3] # reference allele
					alt = linesplit[4] # alternate allele

					if ',' in alt : # if ',' in the alternative allele -> bi, sort, ... allelic -> here only variants with only one allele are taken into account
						pass

					else :
						listeinfo=[] # list with genotype 0 (allele ref) or 1 (allele alt) for each genome of VCF

						for genome in range (9, len(firstliness)): # the genomes concern columns 9 until the end of the file
							info = linesplit[genome].split(":")[0] # the genotype is in the first part of the genome column
							listeinfo.append(info) # we add each genotype of each genome in a list for one and the same variant

						comment2 = ''.join(chromosome) + '_' + ''.join(position) + '_' + ''.join(ref) + '/' + ''.join(alt) + '\t' + '\t'.join([unicode(i) for i in listeinfo]) + '\n' # create each line of the file with form specific name of the variant (chromosome_position_reference / alternative) and genotype associated with each genome
						
						outputfile.write(comment2)

fileNameIn = '/Users/Meryl/Desktop/Analyse_mini_dataset_inra/iVARCall2/Resultat/VCF/Resultat_SNP_INDEL_filtered_ann.vcf' # path of the input file
fileNameOut = '/Users/Meryl/Desktop/Resultat_SNP_INDEL_filtered_ann_binaire.txt' # output file path


#####  MAIN  #####

construct_file(fileNameIn,fileNameOut)
