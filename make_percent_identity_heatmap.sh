
# Steele et al., 2024 BioRxiv: https://doi.org/10.1101/2024.08.09.607352

# How to construct a heatmap showing the distribution of predation-resistance genes among Pseudomonas strains
# In bash:
# 1. Make a BLAST database
# 2. Generate a list/fasta file of sequences disrupted by transposon insertions in edible mutants
# 3. Loop through fasta headers:
# 	BLAST each protein against a database of all Pseudomonas proteins
# 	Filter out hits with >70% coverage and >40% identity

# In R (See: Steele_et_al_2024.R): 
# 1. Import csv file with % identity for top BLAST hits
# 2. Convert data to matrix
# 3. Melt data
# 4. Make heatmap


### Step 1: Prepare a BLAST database containing all protein sequences from all genomes of interest
# Download all protein sequences for each genome from NCBI
# Make a CSV file containing all of the species/strain names (Pseudomonas_strain_list.csv)

# Add species/strain names to the beginning of protein IDs to prevent duplicate IDs in BLAST database 
for i in *.faa; do
	id=$(echo $i | sed 's/_GCF_.*//g' | sed 's/.*\///g');
	name=$(echo $i | sed 's/.faa//g' | sed 's/.*\///g');
	echo "id is $id";
	echo "name is $name";
	more $i | sed "s/>/>${id}_protein_/g" > ${name}.2.faa;
done

# Concatenate sequences
cat *.2.faa > all_pseudomonas_proteins.faa

# Make the BLAST database
makeblastdb -in all_pseudomonas_proteins.faa -out all_pseudomonas_protein_db -dbtype prot -parse_seqids


### Step 2: Prepare list of predation resistance genes

more predation_resistance_protein_sequences.fa | grep '>' > predation_resistance_gene_protein_sequence_list.txt


### Step 3: Loop through predation resistance genes, performing a blastp search for each against the database of all protein sequences from reference genomes

while read line; do id=$(echo $line | sed 's/>//g' | sed 's/ .*//g'); grep -A1 "$line" predation_resistance_protein_sequences.fa > ${id}.fa; echo ${id}.fa; blastp -db all_pseudomonas_protein_db -query ${id}.fa -outfmt "6 sseqid qseqid qcovs pident bitscore" -out ${id}_hits.txt -num_threads 8; done<predation_resistance_gene_protein_sequence_list.txt

# For each predation resistance gene, filter out the best hit from each genome with > 70% coverage

for i in *_hits.txt; do id=$(echo $i | sed 's/_hits.txt//g'); more $i| awk '$3>=70{print $1,$4}' | sed 's/_protein_/ /g' | sort -k1,1 -k3,3nr | sort -u -k1,1 | awk '{print $1","$3}' > ${id}_best_hits.csv; done

# Join percent identity for best blast hits to list of reference genome species/strain names creating a CSV file
# Note: Strain names in key need to be sorted in the same order as strain names in the best blast hits file

for i in *_best_hits.csv; do id=$(echo $i | sed 's/_best_hits.csv//g'); join -t, -1 1 -2 1 -a 2 -e 0 -o 2.1,1.2 $i Pseudomonas_strain_list.csv  > ${id}_best_hits.key.csv; done

mv *_best_hits.key.csv Pseudomonas_predation_resistance_homologs

# Combine results into a single CSV file 

### Visualize heatmap in R using ggplot2
library(RColorBrewer)
library(reshape2)
library(ggplot2)

pident_table <- read.csv("Pseudomonas_predation_resistance_homologs.csv", row.names=1, na.strings="0")
data <- as.matrix(pident_table)
data[is.na(data)] <- 0
data_melt <- melt(data)

ggp <- ggplot(data_melt, aes(Var2, Var1, fill = value)) +                         
  geom_tile(color="black") + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 60)) +
  scale_fill_gradient(low = "white", high = "black") + # #330066
  xlab("Genes important for predation resistance from Tn screen") + ylab("Strain") + labs(fill = "% identity")
ggp

