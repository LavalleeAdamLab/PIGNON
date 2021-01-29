import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Properties;

import graph.Annotation;
import graph.FalseDiscoveryRate;
import graph.Interaction;
import graph.Protein;
import utils.*;


public class Main {

	public static final boolean removeOverlyConnectedProteins = true; // default = TRUE
	public static final int numOfExcessInteractions = 1000;	// default = 1000

	/*** Constants method variables ***/
	public static final boolean toMonteCarlo = false;

	public static final boolean testGoTerms = false;
	public static final boolean cleanClusters = false;

	/**** Constant output variables ***/
	public static final boolean printNumberOfInteractionsPerProtein = false;
	public static final boolean printPPInetwork = false;

	/**** Constant  distribution variables ****/ 
	public static final boolean computeDistributionParams = true;
	public static final boolean useNormalDistribution = true; //if false uses monte carlo values: default = True

	public static void main(String[] args) throws FileNotFoundException, IOException {

		/**
		 * PIGNON stands for a Protein-protein Interaction-Guided fuNctiOnal eNrichment analysis for quantitative proteomics. 
		 * This algorithm assesses the clustering of proteins associated to Gene Ontology terms within the BioGRID or String PPI network. 
		 * This algorithm was implemented on a PPI network weighted with proteomics expression data of various breast cancer
		 * subtype (study by Tyanova et al. 2015). 
		 * 
		 * Input
		 *  - BioGRID (Homo sapiens). version 3.4.158 release March 1st 2018.
		 *  - String (Homo sapiens) 
		 *  - Protein expression file (Tyanova et al. study)
		 *  - GO annotations (Homo sapiens). Gene Ontology Consortium obtained June 2018.
		 *  
		 *  Output
		 *  - FDR to p-value mapping file
		 *  - Summary of GO annotations and clustering statistics
		 *  
		 *  TO UPDATE: 
		 *  - Format input file
		 ***/
	
		/* Load Java properties file which contains parameters for running of the tool */
		System.out.println("Loading parameters file \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));

		String condition1 = params.getProperty("condition1").replace("\\s+", "");
		String condition2 = params.getProperty("condition2").replace("\\s+", "");

		int nProtToSampleLowerBound = Integer.parseInt(params.getProperty("lowerBound").replaceAll("\\s+",""));				// Lower limit of number of proteins sampled
		int nProtToSampleUpperBound = Integer.parseInt(params.getProperty("upperBound").replaceAll("\\s+",""));				// Upper limit of number of proteins sampled
		int numOfTimesNetworkIsSampled = Integer.parseInt(params.getProperty("numOfIterations").replaceAll("\\s+",""));		// Number of times the network is sampled during Monte Carlo Method

		int networkWeights = Integer.parseInt(params.getProperty("networkWeights").replaceAll("\\s+",""));		// weither to use protein expression weights or not
		int networkType = Integer.parseInt(params.getProperty("networkType").replaceAll("\\s+",""));			// specify use of BioGRID or String network
		String taxonomyID = params.getProperty("taxonomyID").replaceAll("\\s+", "");
		
		int samplingType = Integer.parseInt(params.getProperty("samplingType").replace("\\s+", ""));
		
		String projectName = params.getProperty("project_name");	// Specify project name
		
		/* UPDATE working directory */ 
		String working_directory = params.getProperty("working_directory");
		
		/* Create directory for IO files ands output files */ 
		File directory = new File(working_directory + "/IO_files");
	    if (! directory.exists()){
	    	System.out.println("creating IO_files directory");
	        directory.mkdir();
	    }
	    
		File directory2 = new File(working_directory + "/output_files");
	    if (! directory2.exists()){
	    	System.out.println("creating output_files directory\n");
	        directory2.mkdir();
	    }

		/* input files */
		String interactionNetwork_inputFile = working_directory+ params.getProperty("protein_interaction_repository").replaceAll("\\s+", "");
		String ensembleToEntrezIdMap_inputFile = working_directory + params.getProperty("ensembleIdToEntrezIdFile").replaceAll("\\s+", "");		// required for String network
		String protein_expression_file = working_directory+ params.getProperty("protein_expression_data").replaceAll("\\s+", "");
		String annotationGO_inputFile = working_directory + params.getProperty("gene_ontology_file").replaceAll("\\s+", "");
		String shuffledGOFile = working_directory +"IO_files/" +"shuffled_gene_ontologies.txt";
		
		/* Specify network type name for output file writing */ 
		String ppiNetwork = "";
		switch(networkType) {
		case 0:
			ppiNetwork = "BioGRID";
			break;
		case 1: 
			ppiNetwork = "String";
			break;
		}

		String samplingName = "";
		switch(samplingType) {
		case 0:
			samplingName = "unweightedSampling";
			break;
		case 1: 
			samplingName = "weightedSampling";
			break;
		}
		
		String networkWeighted = "";
		switch(networkWeights) {
		case 0 :
			networkWeighted = "not weighted";
			break;
		case 1: 
			networkWeighted = "weighted";
			break;
		}
		
		/* intermediate files */
		String distanceMatrixFile = working_directory + "IO_files/" +projectName + "_" + ppiNetwork + "_DistanceMatrix.txt";
		String distanceMatrix2File = working_directory + "IO_files/" + projectName + "_" + ppiNetwork + "_DistanceMatrix_FullConnected.txt";

		String distributionFile = working_directory + "IO_files/" + projectName + "_" + ppiNetwork + "_" + samplingName + "_s" + numOfTimesNetworkIsSampled + "_" +
				nProtToSampleLowerBound + "_" + nProtToSampleUpperBound + ".txt";

		String normalDistributionParametersFile = working_directory + "IO_files/" + projectName + "_" + ppiNetwork + "_" + samplingName + "_s" + numOfTimesNetworkIsSampled+"_normalDistributionParams.txt";

		/* output files  */
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss.");	
		Timestamp timestamp = new Timestamp(System.currentTimeMillis());		// get current time stamp in above specified format

		String fdrExportFile = working_directory + "output_files/" + sdf.format(timestamp) + projectName + "_" + ppiNetwork + "_" + samplingName + "_MonoTransf_FDRvPval_s" + numOfTimesNetworkIsSampled + ".txt";
		String goExportFile = working_directory + "output_files/" + sdf.format(timestamp) + projectName + "_" + ppiNetwork +  "_" + samplingName + "_SummaryGoTerms_s" + numOfTimesNetworkIsSampled + ".txt";
		String goDetailsFile = working_directory + "output_files/" + sdf.format(timestamp) + projectName + "_" + ppiNetwork + "_" + samplingName + "_DetailedGoTerms_s" + numOfTimesNetworkIsSampled + ".txt";;
		
		System.out.println("Running PIGNON job: " + projectName + "\n" + "PPI network: " + ppiNetwork + " is " + networkWeighted + " with quantification data");
		if(networkWeights == 1) {
			System.out.println("condition1 = " + condition1 + "| condition 2 = " + condition2);
		}
		System.out.println("Null model using " + samplingName + " for annotations of size " + nProtToSampleLowerBound + " - " + nProtToSampleUpperBound + " with " + numOfTimesNetworkIsSampled + " sampling\n");
		//String mclGraphFile = working_directory + "output_files/unweightedMCLgraph.txt";
		//String proteinListFile = working_directory + "output_files/2020.10.15.proteinsInNetworkList_wOverConnectedProteins.txt";
		//		String nInteractionsPerProtFile = "/Users/Rachel/eclipse-files/network2.0/output_files/proteins.txt";


		/* Read PPI repository. Extract all possible interactions between proteins and store in Interaction object */
		System.out.println("Loading interaction repository");
		ArrayList<Interaction> networkInteractionsList = NetworkInteractionsLoader.importInteractionNetwork(interactionNetwork_inputFile, networkType, ensembleToEntrezIdMap_inputFile, removeOverlyConnectedProteins, numOfExcessInteractions, taxonomyID);
		System.out.println("number of interactions: " + networkInteractionsList.size());

		/* Extract all proteins in the network (ie. all proteins within the interaction file) */
		ArrayList<Protein> networkProteinList = NetworkProteins.getProteinsInNetwork(networkInteractionsList);
		System.out.println("Number of proteins: " + networkProteinList.size() + "\n");

		if(networkWeights == 1) {
			System.out.println("Loading protein quantification data");
			/* Load protein expression data */
			Loader.loadProteinFoldChangeData(protein_expression_file, networkProteinList, condition1, condition2);

			/* Modify interaction weight based on protein expression */
			Modifier.modifyInteractionWeight(networkProteinList, networkInteractionsList);
		}

		File f = new File(distanceMatrixFile);
		if(!f.exists() && !f.isDirectory()) {
			/* Compute distance matrix using the list of proteins in the network and the list of known interactions in the network. */
			System.out.println("Computing distance matrix");
			DistanceMatrix.computeDistanceMatrix(networkInteractionsList, networkProteinList, distanceMatrixFile);
		}

		/* Load distance matrix (Note: it contains disconnected components) */ 
		System.out.println("Loading distance matrix");
		double[][] distanceMatrix = Loader.loadDistanceMatrix(distanceMatrixFile, networkProteinList);

		/* Determine which proteins are disconnected*/ 
		boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);

		/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
		ArrayList<Protein> networkProteinsList3 = NetworkProteins.modifyNetworkProteinsList(networkProteinList, proteinsToKeep);

		if(networkProteinList.size() != networkProteinsList3.size()) {
			File f1 = new File(distanceMatrix2File);
			if(!f1.exists() && !f1.isDirectory()) {
				/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
				System.out.println("Updating distance matrix");
				DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, distanceMatrix2File);
			}
			/* Load distance matrix representing fully connected component */
			System.out.println("Loading updated distance matrix");
			distanceMatrix = Loader.loadDistanceMatrix(distanceMatrix2File, networkProteinsList3);
		} 

		/* Load annotations */ 
		System.out.println("Loading annotations file");
		ArrayList<Annotation> annotationGoList = Loader.importAnnotationGo(annotationGO_inputFile, networkProteinsList3);
		System.out.println("Loaded annotations: " + annotationGoList.size() + "\n");

		File f4 = new File(distributionFile);
		if(!f4.exists() && !f4.isDirectory()) {
			/* Measure distributions for given number of proteins within the bounds for a given number of sampling */
			System.out.println("Performing Monte Carlo Sampling procedure\n");
			Sampling sampling = new Sampling(annotationGoList, distanceMatrix, samplingType);
			sampling.computeMultipleDistributions(nProtToSampleLowerBound, nProtToSampleUpperBound, numOfTimesNetworkIsSampled, distributionFile);
		}

		ArrayList<Annotation> shuffled_goAnnotationList = new ArrayList<Annotation>();
		File f3 = new File(shuffledGOFile);

		if(!f3.exists() && !f3.isDirectory()) {
			/* Shuffle proteins associated to annotations */
			System.out.println("Shuffling annotations");
			shuffled_goAnnotationList = Calculator.shuffleGoProtAssociations2(annotationGoList);

			/* Export shuffled file to not repeat */ 
			System.out.println("Printing shuffled annotations\n");
			Exporter.printShuffledNetwork(shuffled_goAnnotationList, shuffledGOFile);

		} else {
			/* Load shuffled annotations */ 
			Loader.loadShuffledGoAnnotations(shuffled_goAnnotationList, shuffledGOFile);
		}

		/* Initialize FdrCalculator object; compute clustering of proteins for all GO annotation */ 
		FdrCalculator fdrCalculator = new FdrCalculator(annotationGoList, shuffled_goAnnotationList);
		fdrCalculator.modifyGoAnnotationsWithTPD(distanceMatrix);

		/* Compute mean and standard deviation from Monte Carlo Distribution file */ 
		File f2 = new File(normalDistributionParametersFile);
		if(!f2.exists() && !f2.isDirectory()) {
			System.out.println("Loading Monte Carlo distribution\n");
			fdrCalculator.computeNormalDistributionParameters(distributionFile, nProtToSampleLowerBound, nProtToSampleUpperBound, normalDistributionParametersFile);
		}

		/* Assess the significance of the clustering measure */
		double minimum_pval = 1; 
		System.out.println("Computing significance scores");
		if(useNormalDistribution) {
			minimum_pval = fdrCalculator.modifyGoAnnotationsWithPvalueFromNormalApproximation(normalDistributionParametersFile, numOfTimesNetworkIsSampled);
		} else {
			minimum_pval = Loader.setAnnotationPvaluesFromMonteCarloDistribution(distributionFile, nProtToSampleLowerBound, nProtToSampleUpperBound, annotationGoList, shuffled_goAnnotationList);
		}

		/* Estimate the false discovery rate at various significant thresholds and apply monotonic transformations */ 
		System.out.println("Performing FDR estimation procedure");
		ArrayList<FalseDiscoveryRate> fdrs = fdrCalculator.computeFdr(minimum_pval);
		fdrs = fdrCalculator.monotonicTransformationForFdr(fdrs);

		Modifier.modifyClustersFdr(annotationGoList, fdrs);

		/* Print stats summary for each tested GO annotation */
		System.out.println("Print GO results"); 
		Exporter.printGO_results(annotationGoList, shuffled_goAnnotationList, goExportFile);

		FormatOutput.printAnnotationDetails(annotationGoList, networkProteinsList3, fdrs, annotationGO_inputFile, goDetailsFile);


		/* Print FDR calculated at various p-value thresholds mapping */ 

		System.out.println("Print FDR"); 
		Exporter.testGoFDR(fdrs, fdrExportFile);



	} // close main

} // end Network
