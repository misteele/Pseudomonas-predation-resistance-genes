while read line; do		
    echo "Reading line $line"

 # Use grep to extract FASTA sequence for protein in wkB2
    grep -A1 "$line" all_predation_resistance_proteins.2.fa > ${line}.fa
    #grep -A1 "$line" Pseudomonas_reference_protein_sequences_for_phylogeny.2.fa > ${line}.fa
    echo "Query file: ${line}.fa Sequence:" >&1
    head ${line}.fa >&1

    blastp -db Pseudomonas_20240617_db -query ${line}.fa -outfmt "6 sseqid qseqid qcovs pident bitscore" -out ${line}_blastp.txt -num_threads 8 -max_target_seqs 500
    echo "blastp output:" >&1
    head ${line}_blastp.txt >&1

### Filter out any BLAST hits with less than 70% coverage
### Put a comma and tab between the genome assembly accession number and protein accession number to make it possible to filter out the best matching protein from each genome

    echo "Filtering blastp output:" >&1
    echo "Split genome accession and protein accession:" >&1 
    awk '$3>=70{print $1","$5}' ${line}_blastp.txt | sed 's/_protein/,/g' | sort -t, -k1,1 -k3,3nr | head >&1 # print to output to make sure sorting worked as expected
    echo " " >&1
                
### Next, filter BLAST hits and split qseqids as above, then filter out the best hit per genome (based on bitscore)
### The protein accession number is then rejoined to the genome assembly accession number
### The sseqid identifiers have to be returned to their original state to use blastdbcmd in the next step
                
    echo "Fuse genome accession and protein accession after filtering best match for genome:" >&1
    awk '$3>=70{print $1","$5}' ${line}_blastp.txt | sed 's/_protein/,/g' | sort -t, -k1,1 -k3,3nr | sort -t, -u -k1,1 | sed 's/,_/_protein_/g' | sed 's/,\./_protein./g' | awk -F, '{print $1}' > ${line}_blastp.hits.txt

    echo "Filtered blastp hits:" >&1
    head ${line}_blastp.hits.txt >&1 # checking to make sure that the filtered hits file looks okay

### Use blastdbcmd to extract the protein sequences for the top blast hits
	blastdbcmd -db Pseudomonas_20240617_db -entry_batch ${line}_blastp.hits.txt -out ${line}_seqs.fa
    echo "blastdbcmd output:" >&1
    head ${line}_seqs.fa >&1
 				
### Use muscle to align proteins, then remove the protein accession numbers from the alignment (keeping only the genome accession number, which will be used to concatenate the alignments of all proteins encoded by genes within the reference T6SS locus (i.e., a single T6SS)
    muscle -in ${line}_seqs.fa -out ${line}_seqs.aligned.fa
    sed 's/_protein.*$//g' ${line}_seqs.aligned.fa > ${line}_seqs.aligned.renamed.fa
                
done<Pseudomonas_predation_resistance_protein_list.txt #Pseudomonas_reference_protein_list.txt

