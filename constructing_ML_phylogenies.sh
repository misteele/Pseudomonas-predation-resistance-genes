# Steele et al., 2024 BioRxiv: https://doi.org/10.1101/2024.08.09.607352
# Instructions for constructing Maximum Likelihood phylogenies

### Step 1: Make a BLAST database containing all proteins from all of the strains that will be included in the phylogeny

# Download all protein sequences for each genome from NCBI

# Rename sequences in .faa files so there won't be duplicate IDs when they are concatenated later
for i in *.faa; do
	id=$(echo $i | sed 's/_GCF_.*//g' | sed 's/.*\///g');
	name=$(echo $i | sed 's/.faa//g' | sed 's/.*\///g');
	echo "id is $id";
	echo "name is $name";
	more $i | sed "s/>/>${id}_protein_/g" > ${name}.2.faa;
done

# Concatenate sequences
cat *.2.faa > Pseudomonas_faa_concat_for_ML_tree.faa

# Make the BLAST database
makeblastdb -in Pseudomonas_faa_concat_for_ML_tree.2.faa -out Pseudomonas_20240617_db -dbtype prot -parse_seqids


### Step 2: Prepare a fasta file and list of proteins of interest
# I constructed my list of reference proteins by searching NCBI for Pseudomonas fluorescens proteins from the UBCG gene list (http://leb.snu.ac.kr/ubcg2/genes) and copy-pasting the amino acid sequences into a single text document

# Remove linebreaks, spaces, brackets in the reference FASTA files
more Pseudomonas_reference_protein_sequences_for_phylogeny.fa | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' | sed 's/\ \[.*\]//g' | sed 's/\ /_/g' > Pseudomonas_reference_protein_sequences_for_phylogeny.2.fa
# Make list of proteins to loop through
more Pseudomonas_reference_protein_sequences_for_phylogeny.2.fa | grep '>' | sed 's/>//g' > Pseudomonas_reference_protein_list.txt

# Remove linebreaks, spaces, brackets in the reference FASTA files
more all_predation_resistance_proteins.fa | awk '/^>/ {if(N>0) printf("\n"); printf("%s\n",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' | sed 's/\ \[.*\]//g' | sed 's/\ /_/g' | sed 's/\//-/g'> all_predation_resistance_proteins.2.fa
# Make list of proteins to loop through
more all_predation_resistance_proteins.2.fa | grep '>' | sed 's/>//g' > Pseudomonas_predation_resistance_protein_list.txt


### Step 3: Run script to build alignments
./make_alignments_for_ML_phylogenies.sh

# Move aligned sequences into new directory
mv *.aligned.renamed.fa Protein_Alignments  
mv *_seqs.fa *_seqs.aligned.fa *_blastp.txt *_blastp.hits.txt Other_Alignment_Files

# Import alignments into Geneious and concatenate, then export concatenated alignments in phylip format

### Step 4: Build phylogenies with RAxML or IQ-TREE
# Use RAxML
raxmlHPC-AVX2 -f a -s Pseudomonas_conserved_protein_concatenated_alignment.phy -m PROTGAMMAAUTO -p 423643 -x 87295 -N 100 -n Pseudomonas_conserved.raxml_tree -T 8

# Use IQ-TREE
iqtree -s Pseudomonas_78core.phy -m MFP -bb 1000 -alrt 1000 -nt AUTO

# Loop through a list of phylip alignments and use IQ-TREE to construct phylogenies for each
LIST1="7A10 7G12 17D1 17F4"

for i in $LIST1; do
	echo 'Running analysis for ${i}'
	iqtree -s ${i}.phy -m MFP -bb 1000 -alrt 1000 -nt AUTO 			# This generates prot_tree for tanglegrams (see Steele_et_al_2024.R)
	iqtree -s ${i}_core.phy -m MFP -bb 1000 -alrt 1000 -nt AUTO		# This generates spec_tree for tanglegrams (see Steele_et_al_2024.R)
done


