import java.util.ArrayList;

import graph.Annotation;
import graph.FalseDiscoveryRate;
import graph.Interaction;
import graph.Protein;
import utils.*;


public class Main {

	/*** Constant network variables ***/
	public static final boolean testBreastCancerWeightedNetwork = true;
	public static final boolean testHer2vsTN = true;
	public static final boolean testHer2vsHR = false;
	public static final boolean testHRvsTN = false;
	
	public static final boolean removeOverlyConnectedProteins = false; 
	public static final int numOfExcessInteractions = 1000;

	/*** Constants method variables ***/

	public static final boolean toComputeDistanceMatrix = false;
	public static final boolean toUpdateDistanceMatrix = false;

	public static final boolean toMonteCarlo = false;
	public static final boolean weightedSampling = true;
	public static final boolean testNormalApprox  = true;

	public static final boolean testGoTerms = false;
	public static final boolean printFDRresults = true;
	public static final boolean printSignificantGOterms = true;
	public static final boolean cleanClusters = false;

	/*** Constant running variables ***/
	public static final boolean runComputeCanada = false;

	/**** Constant output variables ***/
	public static final boolean printNumberOfInteractionsPerProtein = false;
	public static final boolean printPPInetwork = false;

	public static final boolean toShuffle = false;
	public static final boolean printShuffledGOList = true;

	/**** Constant  distribution variables ****/ 
	public static final boolean computeDistributionParams = true; //if true don't recalculate parameters
	public static final boolean useNormalDistribution = true; //if false uses monte carlo values

	public static void main(String[] args) {

		/**
		 * PIGNON stands for a Protein-protein Interaction-Guided fuNctiOnal eNrichment analysis for quantitative proteomics. 
		 * This algorithm assesses the clustering of proteins associated to Gene Ontology terms within the BioGRID PPI network. 
		 * This algorithm was implemented on a PPI network weighted with proteomics expression data of various breast cancer
		 * subtype (study by Tyanova et al. 2015). 
		 * 
		 * Input
		 *  - BioGRID (Homo sapiens). version 3.4.158 release March 1st 2018.
		 *  - Protein expression file (Tyanova et al. study)
		 *  - GO annotations (Homo sapiens). Gene Ontology Consortium obtained June 2018.
		 *  
		 *  Output
		 *  - FDR to p-value mapping file
		 *  - Summary of GO annotations and clustering statistics
		 ***/

		int nProtToSampleLowerBound = Integer.parseInt(args[0]);		// Lower limit of number of proteins sampled
		int nProtToSampleUpperBound = Integer.parseInt(args[1]);		// Upper limit of number of proteins sampled
		int numOfTimesNetworkIsSampled = Integer.parseInt(args[2]);		// Number of times the network is sampled during Monte Carlo Method
		
		/* UPDATE working directory */ 
		String working_directory = "/Users/Rachel/eclipse-files/network2020/";
		
		/* File paths when running on remote cluster */
		String bioGRID_inputFile = "BIOGRID-ORGANISM-Homo_sapiens-3.4.161.tab2.txt";
		String protein_expression_file = "ncomms10259-BreastCancerProteinExpression.txt";
		String annotationGO_inputFile = "GO_annotations-9606-inferred-allev-2.tsv";

		String distanceMatrixFile = "unweighted_DistanceMatrix.txt";
		String distanceMatrix2File = "unweighted_DistanceMatrix_FullConnected.txt";

		String distributionFile = "unweightedNetwork2.0_s10^7_n1_1000.txt";

		if (testBreastCancerWeightedNetwork) {
			if (testHer2vsTN) {
				distanceMatrixFile = "her2vTNegWeighted_DistanceMatrix.txt";
				distanceMatrix2File = "her2vTNegWeighted_DistanceMatrix2_FullConnected.txt";
				
				if(!removeOverlyConnectedProteins) {
					distanceMatrixFile = "her2vTNegWeighted_wOverConnectedProteins_DistanceMatrix.txt";
					distanceMatrix2File = "her2vTNegWeighted_wOverConnectedProteins_DistanceMatrix2_FullConnected.txt";
				}
				
				if(weightedSampling) {
					distributionFile = "weightedNetwork2.0_s10^7_n1_1000.txt";
				} else {
					distributionFile = "unweightedNetwork2.0_s10^7_n1_1000.txt";
				}

			} else if(testHer2vsHR) {
				distanceMatrixFile = "her2vHRWeighted_DistanceMatrix.txt";
				distanceMatrix2File = "her2vHRWeighted_DistanceMatrix_FullConnected.txt";
			} else if(testHRvsTN) {
				distanceMatrixFile = "HRvTNegWeighted_DistanceMatrix.txt";
				distanceMatrix2File = "HRvTNegWeighted_DistanceMatrix_FullConnected.txt";
			}
		}
		
		/* UPDATE File paths when running on local computer */
		if (!runComputeCanada) {
			bioGRID_inputFile = working_directory + "input_files/BIOGRID-ORGANISM-Homo_sapiens-3.4.161.tab2.txt";
			protein_expression_file = working_directory + "input_files/ncomms10259-BreastCancerProteinExpression.txt";
			annotationGO_inputFile = working_directory + "input_files/GO_annotations-9606-inferred-allev-2.tsv";

			distanceMatrixFile = working_directory + "IO_files/unweighted_DistanceMatrix.txt";
			distanceMatrix2File = working_directory + "IO_files/unweighted_DistanceMatrix_FullConnected.txt";

			if(weightedSampling) {
				distributionFile = working_directory + "IO_files/unweightNetwork_weightedSampling_s10X7_n3_1000.txt";
			} else {
				distributionFile = working_directory + "IO_files/unweightedNetwork_unweightedSampling_s10X7_n3_1000.txt";
			}

			if (testBreastCancerWeightedNetwork) {
				if (testHer2vsTN) {
					distanceMatrixFile = working_directory + "IO_files/her2vTNegWeighted_DistanceMatrix.txt";
					distanceMatrix2File = working_directory + "IO_files/her2vTNegWeighted_DistanceMatrix2_FullConnected.txt";
					
					if(!removeOverlyConnectedProteins) {
						distanceMatrixFile = working_directory + "IO_files/her2vTNegWeighted_wOverConnectedProteins_DistanceMatrix.txt";
						distanceMatrix2File = working_directory + "IO_files/her2vTNegWeighted_wOverConnectedProteins_DistanceMatrix2_FullConnected.txt";
					}

					if(weightedSampling) {
						distributionFile = working_directory + "IO_files/her2vTNeg_weightedSampling_updatedGO_s10X7_n3_1000.txt";
						
						if(!removeOverlyConnectedProteins) {
							distributionFile = working_directory + "IO_files/her2vTNegWeightNetworkFullConnected_weightedSampling_s10X7_n3_1000.txt";
						}
					} else {
						distributionFile = working_directory + "IO_files/her2vTNegWeightedNetwork_unweightedSampling_s10X7_n3_1000.txt";
					}

				}  else if(testHer2vsHR) {
					distanceMatrixFile = working_directory + "IO_files/her2vHRWeighted_DistanceMatrix.txt";
					distanceMatrix2File = working_directory + "IO_files/her2vHRWeighted_DistanceMatrix_FullConnected.txt";

					distributionFile = working_directory + "IO_files/her2vHRWeightNetwork_weightedSampling_s10X7_n3_1000.txt";
				} else if(testHRvsTN) {
					distanceMatrixFile = working_directory + "IO_files/HRvTNegWeighted_DistanceMatrix.txt";
					distanceMatrix2File = working_directory + "IO_files/HRvTNegWeighted_DistanceMatrix_FullConnected.txt";
					
					distributionFile = working_directory + "IO_files/HRvTNegWeightNetwork_weightedSampling_s10X7_n3_1000.txt";
				}
			}
		}

		/* Output files */
		String fdrExportFile = working_directory + "output_files/2020.02.19.MonoTransf_FDRvPval_unweightedNetwork_weightedSampling10X7.txt";
		String goExportFile = working_directory + "output_files/2020.02.19.MonoTransf_GoTerms_unweightedNetwork_weightedSampling10X7_newGOFile.txt";

		String shuffledGOFile = working_directory + "IO_files/20202.02.19_geneOntologiesShuffled2.txt";
		String normalDistributionParametersFile = working_directory + "IO_files/unweightedNetwork_weightedSampling_normalDistributionParams.txt";

		if(testHer2vsTN) {
			if(weightedSampling) {
				fdrExportFile = working_directory + "output_files/2020.02.19.MonoTransf_FDRvPval_her2vTNweightedSampling10X7_newGOFile.txt";
				goExportFile = working_directory + "output_files/2020.02.19.MonoTransf_GoTerms_her2vTNweightedSampling10X7_newGOFile.txt";
				normalDistributionParametersFile = working_directory + "IO_files/her2vTNeg_weightedSampling_normalDistributionParams.txt";
				if(!removeOverlyConnectedProteins) {
					fdrExportFile = working_directory + "output_files/2020.08.07.MonoTransf_FDRvPval_her2vTNweightedSampling_wOverConnectedProteins10X7.txt";
					goExportFile = working_directory + "output_files/2020.08.07.MonoTransf_GoTerms_her2vTNweightedSampling_wOverConnectedProteins10X7.txt";
					normalDistributionParametersFile = working_directory + "IO_files/her2vTNeg_weightedSampling__wOverConnectedProteins_normalDistributionParams.txt";
					shuffledGOFile =  working_directory + "IO_files/overconnectedNetwork_shuffledGOannotations.tsv";
				}
			} else {
				fdrExportFile = working_directory + "output_files/2020.03.31.MonoTransf_FDRvPval_her2vTN_unweightedSampling10X7.txt";
				goExportFile = working_directory + "output_files/2020.03.31.MonoTransf_GoTerms_her2vTN_unweightedSampling10X7.txt";
				normalDistributionParametersFile = working_directory + "IO_files/her2vTNeg_unweightedSampling_normalDistributionParams.txt";
			}
		} else if(testHer2vsHR) {
			fdrExportFile = working_directory + "output_files/2020.07.24.MonoTransf_FDRvPval_her2vHRweightedSampling10X7.txt";
			goExportFile = working_directory + "output_files/2020.07.24.MonoTransf_GoTerms_her2vHRweightedSampling10X7.txt";
			normalDistributionParametersFile = working_directory + "IO_files/her2vHR_weightedSampling_normalDistributionParams.txt";
		} else if(testHRvsTN) {
			fdrExportFile = working_directory + "output_files/2020.06.23.MonoTransf_FDRvPval_HRvTNweightedSampling10X7.txt";
			goExportFile = working_directory + "output_files/2020.06.23.MonoTransf_GoTerms_HRvTNweightedSampling10X7.txt";
			normalDistributionParametersFile = working_directory + "IO_files/HRvTNeg_weightedSampling_normalDistributionParams.txt";
		} else {
			if(weightedSampling) {
				fdrExportFile = working_directory + "output_files/2020.02.19.MonoTransf_FDRvPval_unweightedNetwork_weightedSampling10X7.txt";
				goExportFile = working_directory + "output_files/2020.02.19.MonoTransf_GoTerms_unweightedNetwork_weightedSampling10X7_newGOFile.txt";
				normalDistributionParametersFile = working_directory + "IO_files/unweightedNetwork_weightedSampling_normalDistributionParams.txt";
			} else { 
				fdrExportFile = working_directory + "output_files/2020.04.16.MonoTransf_FDRvPval_unweightedNetwork_unweightedSampling10X7.txt";
				goExportFile = working_directory + "output_files/2020.04.16.MonoTransf_GoTerms_unweightedNetwork_unweightedSampling10X7.txt";
				normalDistributionParametersFile = working_directory + "IO_files/unweightedNetwork_unweightedSampling10X7_normalDistributionParams.txt";
			}
		}
		//String mclGraphFile = working_directory + "output_files/unweightedMCLgraph.txt";
		String proteinListFile = working_directory + "output_files/2020.10.15.proteinsInNetworkList_wOverConnectedProteins.txt";
		//		String nInteractionsPerProtFile = "/Users/Rachel/eclipse-files/network2.0/output_files/proteins.txt";


		/* Read BioGRID file. Extract all possible interactions between proteins and store in Interaction object */
		ArrayList<Interaction> networkInteractionsList = NetworkInteractionsLoader.importInteractionNetwork(bioGRID_inputFile, runComputeCanada, removeOverlyConnectedProteins, numOfExcessInteractions);

		/* Extract all proteins in the network (ie. all proteins within the interaction file) */
		ArrayList<Protein> networkProteinList = NetworkProteins.getProteinsInNetwork(networkInteractionsList);

		if(testBreastCancerWeightedNetwork) {
			/* Load protein expression data */
			Loader.loadProteinFoldChangeData(protein_expression_file, networkProteinList, testHer2vsTN, testHer2vsHR, testHRvsTN, runComputeCanada);

			/* Modify interaction weight based on protein expression */
			Modifier.modifyInteractionWeight(networkProteinList, networkInteractionsList);
		}

		if(toComputeDistanceMatrix) {
			/* Compute distance matrix using the list of proteins in the network and the list of known interactions in the network. */
			DistanceMatrix.computeDistanceMatrix(networkInteractionsList, networkProteinList, distanceMatrixFile);
		}

		/* Load distance matrix (Note: it contains disconnected components) */ 
		double[][] distanceMatrix = Loader.loadDistanceMatrix(distanceMatrixFile, networkProteinList, runComputeCanada);

		/* Determine which proteins are disconnected*/ 
		boolean[] proteinsToKeep = Calculator.determineConnectedProteins(distanceMatrix);

		/* Update proteins to keep in the network (ie. those that contribute to the most connected component) */
		ArrayList<Protein> networkProteinsList3 = NetworkProteins.modifyNetworkProteinsList(networkProteinList, proteinsToKeep);

		if(toUpdateDistanceMatrix) {
			/* Update distance matrix to only include proteins that contribute to the most connected component of the graph */
			DistanceMatrix.updateDistanceMatrix(proteinsToKeep, distanceMatrix, distanceMatrix2File);
		}

		/* Load distance matrix representing fully connected component */
		double[][] distanceMatrix2 = Loader.loadDistanceMatrix(distanceMatrix2File, networkProteinsList3, runComputeCanada);

		/* Load annotations */ 
		ArrayList<Annotation> annotationGoList = Loader.importAnnotationGo(annotationGO_inputFile, networkProteinsList3, runComputeCanada);

		if(toMonteCarlo) {
			/* Measure distributions for given number of proteins within the bounds for a given number of sampling */
			Sampling sampling = new Sampling(annotationGoList, distanceMatrix2, weightedSampling);
			sampling.computeMultipleDistributions(nProtToSampleLowerBound, nProtToSampleUpperBound, numOfTimesNetworkIsSampled);
		}

		if(testGoTerms) {		
			ArrayList<Annotation> shuffled_goAnnotationList = new ArrayList<Annotation>();
			if(toShuffle) {	
				/* Shuffle proteins associated to annotations */
				shuffled_goAnnotationList = Calculator.shuffleGoProtAssociations2(annotationGoList);
				if (printShuffledGOList) {
					
					/* Export shuffled file to not repeat */ 
					Exporter.printShuffledNetwork(shuffled_goAnnotationList, shuffledGOFile);
				}
			} else {
				/* Load shuffled annotations */ 
				Loader.loadShuffledGoAnnotations(shuffled_goAnnotationList, shuffledGOFile);
			}

			/* Initialize FdrCalculator object; compute clustering of proteins for all GO annotation */ 
			FdrCalculator fdrCalculator = new FdrCalculator(annotationGoList, shuffled_goAnnotationList);
			fdrCalculator.modifyGoAnnotationsWithTPD(distanceMatrix2);
			
			/* Compute mean and standard deviation from Monte Carlo Distribution file */ 
			if(computeDistributionParams) {
				fdrCalculator.computeNormalDistributionParameters(distributionFile, nProtToSampleLowerBound, nProtToSampleUpperBound, normalDistributionParametersFile);
			}
			
			/* Assess the significance of the clustering measure */
			double minimum_pval = 1; 
			if(useNormalDistribution) {
				minimum_pval = fdrCalculator.modifyGoAnnotationsWithPvalueFromNormalApproximation(normalDistributionParametersFile);
			} else {
				minimum_pval = Loader.setAnnotationPvaluesFromMonteCarloDistribution(distributionFile, nProtToSampleLowerBound, nProtToSampleUpperBound, annotationGoList, shuffled_goAnnotationList);
			}

			/* Estimate the false discovery rate at various significant thresholds and apply monotonic transformations */ 
			ArrayList<FalseDiscoveryRate> fdrs = fdrCalculator.computeFdr(minimum_pval);
			fdrs = fdrCalculator.monotonicTransformationForFdr(fdrs);

			Modifier.modifyClustersFdr(annotationGoList, fdrs);

			/* Print stats summary for each tested GO annotation */
			if (printSignificantGOterms) {
				System.out.println("Print GO results"); 
				Exporter.printGO_results(annotationGoList, shuffled_goAnnotationList, goExportFile);
			}
			
			/* Print FDR calculated at various p-value thresholds mapping */ 
			if (printFDRresults) {
				System.out.println("Print FDR"); 
				Exporter.testGoFDR(fdrs, fdrExportFile);
			}

		}

	} // close main

} // end Network
