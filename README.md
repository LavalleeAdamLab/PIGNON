# PIGNON
PIGNON is a protein-protein interaction (PPI)-guided functional enrichment analysis for quantitative proteomics. 
The bases of the algorithm measures the clustering of proteins with a shared Gene Ontology (GO) annotation within the provided PPI network weighted with quantitative proteomics data. The significance of this clustering measure is then estimated from a normal distribution approximated from a Monte Carlo Sampling Distribution. To correct for multiple hypothesis testing, we assess the false discovery rate at various thresholds against a null model. We tested PIGNON using a breast cancer dataset generated by Tyanova et al. 

PIGNON is a Java application that can be run from the command line. You will need to download the **PIGNON.jar** file. 

## Dependencies
* Java Version 8+
* Required library: The Apache Commons Mathematics Library ([commons-math3-3.6.1.jar](http://commons.apache.org/proper/commons-math/download_math.cgi))

## File Descriptions
### Required input files
Examples files can be found under: *input_files*
1. **Protein-protein interaction network (either [BioGRID](https://downloads.thebiogrid.org/BioGRID) or [STRING](https://string-db.org/cgi/download.pl) PPI network format unzipped)**
2. **String ID to Entrez ID mapping file (required to run STRING network))**

   This tab-delimited text file can be generated using [BioMart](https://bioconductor.org/packages/release/bioc/tml/biomaRt.html) and formatted as the following:
   | ensembl_peptide_id | entrezgene_id | hgnc_symbol |
   | -------------------|---------------|-------------|
   |ENSP00000338785     |90627          |STARD13      |
  
3. **Propagated [Gene Ontology terms](https://git.dhimmel.com/gene-ontology/)**
4. **Quantitative proteomics dataset (ncomms10259-BreastCancerProteinExpression.txt)**

   This tab delimited text file should contain the HGNC_symbol for genes in the 1st column. Columns 2 to n will be the quantified values in the 2 studied conditions. The order of these columns are not important. More than 2 conditions can be included in the same file. It is important that the name of the conditions correspond to those identified in the properties file (eg. if in the properties file condition1 = Her2, in this file all columns for condition1 should be labelled Her2.n) Any missing values should be represented by **NA**.

   |HGNC_symbol | Condition1.1 | Condition1.2 | Condition2.1 | Condition2.2 | ConditionX.n | 
   |------------|--------------|--------------|--------------|--------------|--------------|
   |STARD13     |2.5           |1.8           |0.7           |NA            | ...          |

5. **Properties file (params.txt)** 


### Intermediate files PIGNON will generate
* Initial distance matrix
* Final distance matrix of fully connected component
* Monte Carlo distribution
* Normal Distribution parameters calculated from the Monte Carlo Distribution File
* shuffle Gene Ontology set

### Files PIGNON will output 
* false discovery rates at significant thresholds mapping (.tsv)
* Stats summary of tested GO terms (.tsv)
* Detailed results for every GO annotation (.tsv)

## To run PIGNON
1. Download the repository.
2. You will need to update the working directory in the **Main.java** class.
3. Compile the algorithm: navigate to the folder and input command `javac Main.java graph/*.java utils/*.java opt/*.java`
4. Run the program: `java -Xmx8g Main 3 10 10000`
 
