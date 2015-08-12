# non-canonical miRNA target predictor - ncMirTar

All objective c++ files and perl files can be executed under any Linux system, and windows system with Perl installed.



## Part I: searching targets from the database

### Inputs

1. miRNA id from miRBase: e.g. hsa-miR-103a-3p
2. mRNA accession id from RefSeq: e.g. NM_000546 or mRNA symbol from Genbank: e.g. TP53 
3. directory of the ncMirTar database (A compressed file including all data could be downloaded at 
[http://biomine-ws.ece.ualberta.ca/ncMirTar/downloads/ncMirTar.tar.gz](http://biomine-ws.ece.ualberta.ca/ncMirTar/downloads/ncMirTar.tar.gz) 
please note that the file size is over 2GB)
4. directory of the result

### Execution and Outputs

1. Only Input 1 is given, return all target genes in the given species genome
2. Only Input 2 is given, return all miRNAs that target the given gene or transcript
3. Both Input 1 and Input 2 are given, return their prediction if there is, otherwise say 'no target is found'.

### Examples

**Command line:**

	perl search_ncMirTar.pl hsa-miR-103a-3p NM_000546 ../ncMirTar/ ../result/

**Outputs:**

	@hsa-miR-103a-3p::TP53,NM_000546,6-2,78.96  
	>1187,1.67551  
	AGUAUCGGGACAUGUUACGACGA  
	|: ||   :   ||  ||||||   
	UUUUACAAUAAAACUUUGCUGCC  

**Output format (output_search.txt):**

	@miRNA_id::gene_name,transcript,seed_type,propensity_gene  
	>position_binding,proensity_duplex  
	Target sequence from 5' end to 3' end  
	interaction (| for Watson-Crick base pair, : for GU wobble, and space for mismatch)  
	miRNA sequence from 3' end to 5' end  



## Part II: predicting targets using ncMirTar

### Inputs

1. species (drop down list: human and mouse)
2. miRNA sequence from 5' end to 3' end
3. directory of the ncMirTar database (A compressed file including all data could be downloaded at 
[http://biomine-ws.ece.ualberta.ca/ncMirTar/downloads/ncMirTar.tar.gz](http://biomine-ws.ece.ualberta.ca/ncMirTar/downloads/ncMirTar.tar.gz) 
please note that the file size is over 2GB)
4. directory of the result
5. procedure mode (0: sequence similarity alignment for user's first request; 1: if the user still wants the prediction for the new miRNA after alignment)

### Execution and Outputs

1. If the sequence matches one miRNA in our database, return the prediction of this miRNA
2. If the sequence is similar with one miRNA in our database (E-value<0.001), ask the user if s/he wants the prediction for this similar miRNA or wait for a few hours for prediction for the new miRNA
3. If the sequence is not similar with any miRNA in our database (E-value>=0.001), predict the targets for this new miRNA

### Examples

#### Example 1

**Command line:**

	perl predict_ncMirTar.pl human CAAAGUGCUUACAGUGCAGGUAG ../ncMirTar/ ../result/ 0

**Outputs:**

	>hsa-let-7a-5p, Length=22, 42.1, 9e-08  
	Query  1   TGAGGTAGTAGGTTGTATAGT  21  
	           |||||||||||||||||||||      
	Sbjct  1   TGAGGTAGTAGGTTGTATAGT  21  

**Output format (output_alignment before prediction.txt):**

	>most_similar_miRNA_id, length_miRNA, alignment_score, alignment_evalue  
	Query starting_position_alignment Query_sequence end_position_alignment  
	alignment  
	Sbjct starting_position_alignment most_similar_sequence end_position_alignment  

#### Example 2

**Command line:**

	perl predict_ncMirTar.pl human CAAAGUGCUUACAGUGCAGGUAG ../ncMirTar/ ../result/ 1

**Outputs:**

	@hsa-miR-new::A2ML1,NM_001282424,6-3,39.11  
	>172,0.687165  
	UAUAAUACUUUCUACUACCUUU  
	||: ||||  :|||||||||:   
	AUGAUAUGUUGGAUGAUGGAGU  

**Output format (output_alignment before prediction.txt):**

	@miRNA_id::gene_name,transcript,seed_type,propensity_gene  
	>position_binding,proensity_duplex  
	Target sequence from 5' end to 3' end  
	interaction (| for Watson-Crick base pair, : for GU wobble, and space for mismatch)  
	miRNA sequence from 3' end to 5' end  