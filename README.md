# PIGNON
PIGNON is a protein-protein interaction (PPI)-guided functional enrichment analysis for quantitative proteomics. 
The bases of the algorithm measures the clustering of proteins with a shared Gene Ontology (GO) annotation within the provided PPI network weighted with quantitative proteomics data. The significance of this clustering measure is then estimated from a normal distribution approximated from a Monte Carlo Sampling Distribution. To correct for multiple hypothesis testing, we assess the false discovery rate at various thresholds against a null model. We tested PIGNON using a breast cancer dataset generated by Tyanova et al. 

PIGNON is a Java application that can be run from the command line. You will need to download the **PIGNON.jar** file. 

Note: We recommend running a first instance of PIGNON on your chosen PPI network with your quantitative data and running a second instance without the quantitative data in order to eliminate results that are significant only due to the innate network topology. 

## Dependencies
* Java Version 8+
* Required library: The Apache Commons Mathematics Library ([commons-math3-3.6.1.jar](http://commons.apache.org/proper/commons-math/download_math.cgi))

## File Descriptions
### Required input files
Examples files can be found under: *input_files*
1. **Protein-protein interaction network (either [BioGRID](https://downloads.thebiogrid.org/BioGRID) or [STRING](https://string-db.org/cgi/download.pl) PPI network format unzipped)**
   
   PIGNON is currently set up to run on either the BioGRID or STRING networks. In the params file: you will need to specify the network type either BioGRID (0) or STRING (1) and the taxonomy ID of the species eg. human (9606).
   
   For an alternative PPI network, you can format your network as tab delimited file where each row is an interaction formatted as specified below. In the params file: you will need to specify the network type either BioGRID (0), you should leave the taxonomy ID blank. Required information (ie. the other columns can be blank columns): 
   * column 2: Entrez ID for interactor 1
   * column 3: Entrez ID for interactor 2
   * column 8: HGNC symbol for interactor 1
   * column 9: HGNC symbol for interactor 2
   * column 16: Species identifier interactor 1
   * column 17: Species identifier interactor 2
    
   |   | EntrezID 1 | EntrezID 2 |   |   |   |   | HGNC symbol 1 | HGNC symbol 2 |   |   |   |   |   |   | SpeciesID 1 | SpeciesID 2 | 
   |---|------------|------------|---|---|---|---|---------------|---------------|---|---|---|---|---|---|-------------|-------------|
   |   |6416	     |2318        |   |   |   |   | MAP2K4        |FLNC           |   |   |   |   |   |   |9606         |	9606

   
2. **String ID to Entrez ID mapping file (required to run STRING network)**

   This tab-delimited text file can be generated using [BioMart](https://bioconductor.org/packages/release/bioc/tml/biomaRt.html) and formatted as the following:
   | ensembl_peptide_id | entrezgene_id | hgnc_symbol |
   | -------------------|---------------|-------------|
   |ENSP00000338785     |90627          |STARD13      |
  
3. **Propagated [Gene Ontology terms](https://git.dhimmel.com/gene-ontology/)**

   PIGNON is currently set up to run using only this type of annotation file. 
   
   Alternatively you can format your annotations as a tab delimited file where every row is a new annotation. Required information: 
   * column 1: Annotation ID
   * column 2: Annotation Name (can be blank)
   * column 3: Annotation descriptor (can be left blank)
   * column 7: List of EntrezGene IDs, where elements are separated by a pipe (|)
   * column 8: List of HGNC symbols, where elements are separated  by a pipe (|) 
   
   | AnnotationID | Annotation Name | Annotation descriptor |   |   |   | EntrezGene IDs | hgnc_symbols |
   | ------------ | --------------- | --------------------- |---|---|---|----------------|--------------|
   |GO:0000015 |phosphopyruvate hydratase complex| cellular_component |   |   |   | 2023\|2026\|2027\|387712 |	ENO1\|ENO2\|ENO3\|ENO4
   
4. **Quantitative proteomics dataset**

   This is a tab delimited text file where each row represents the quantitative information for a given gene/protein in 2 or more conditions. Required information:
   * column 1: HGNC_symbol
   * columns 2-n: quantified values in the 2 studied conditions. It is important that the column labels for each condition corresponds to the labels identified in the params file (ie. if in the params file **condition1 = Her2**, in this file all columns for condition1 should be labelled **Her2.n**). The order of the columns is not important. 
   
   Any missing values should be represented by **NA**.

   |HGNC_symbol | Condition1.1 | Condition1.2 | Condition2.1 | Condition2.2 | ConditionX.n | 
   |------------|--------------|--------------|--------------|--------------|--------------|
   |STARD13     |2.5           |1.8           |0.7           |NA            | ...          |

5. **Params file** 
   
   A template of this text file must be usedto run PIGNON. It is passed to the program as command line argument. 
   
   It is important to specify the working directory, file paths and the proper parameters. A detailed explanation of these parameters can be found [here](file/path/to/doc). 

### Intermediate files PIGNON will generate (fyi)
Note: These files will be generated in a sub-directory of your working directory labelled `IO_files` which will be automatically generated by the program
* Initial distance matrix
* Final distance matrix of fully connected component
* Monte Carlo distribution
* Normal Distribution parameters calculated from the Monte Carlo Distribution File
* shuffle Gene Ontology set

### Files PIGNON will output (fyi)
Note: These files will be generated in a sub-directory of your working directory labelled `output_files` which will be automatically generated by the program
* false discovery rates at significant thresholds mapping (.tsv)
* Stats summary of tested GO terms (.tsv)
* Detailed results for every GO annotation (.tsv)

## To run PIGNON
Note: we recommend running PIGNON on a computer with a minimum of 8GB of RAM. The program can run for up to 24hrs. 

1. Download the `PIGNON.jar` file.
2. Prepare your input files as specified above. 
4. Open a terminal (mac/linux) or command prompt (windows) and navigate to where `PIGNON.jar` is stored. 
5. Enter command: 

   `java -Xmx8g -jar PIGNON.jar file/path/params.txt` 
