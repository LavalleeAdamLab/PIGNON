package utils;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import graph.Annotation;
import graph.Interaction;
import graph.Protein;

public class Loader {
	/***
	 * @deprecated
	 * Read BioGRID repository (from inputFile). Extract name of all possible interactions 
	 * as well as their number of recurrence. Ensures that there are no repeat interactions.
	 * Ensures that all proteins involved in interactions are human (homo sapiens; 9606) 
	 ***/
	public static ArrayList<Interaction> loadInteractionRepository(String inputRepositoryFile, boolean runOnComputeCanada) {


		// To contain interaction name and ID, and number of occurence
		HashMap<String, Integer> InteractionToNumber = new HashMap<String, Integer>();

		try {
			/* Read BioGrid ; get all possible human protein-protein interactions */

			BufferedReader input = null;

			if (runOnComputeCanada) {
				InputStream in = Loader.class.getClassLoader().getResourceAsStream(inputRepositoryFile);
				input = new BufferedReader(new InputStreamReader(in));
			} else {
				input = new BufferedReader(new FileReader(new File(inputRepositoryFile)));
			}

			String line = input.readLine(); // read first line
			line = input.readLine(); // read second line

			while (line != null) { // stops when there are no more lines in file (in)

				String[] col = line.split("\t"); // split line by tabs; obtain individual columns

				/* Do not include protein interactions that involve non human proteins. 
				 * The code for homo sapiens is 9606 */

				int species1 = Integer.parseInt(col[15]); // get species of prot1 
				int species2 = Integer.parseInt(col[16]); // get species of prot2

				if (species1 == 9606 && species2 == 9606) {

					String interactor1 = col[7]; // get name of prot1
					String interactor2 = col[8]; // get name of prot2

					String interactorID1 = col[1]; // get ID of prot1
					String interactorID2 = col[2]; // get ID of prot2

					// establish the two possible interaction combinations 
					String interaction1 = interactor1 + "\t" + interactor2 + "\t" + interactorID1 + "\t"
							+ interactorID2;
					String interaction2 = interactor2 + "\t" + interactor1 + "\t" + interactorID2 + "\t"
							+ interactorID1;

					/* Look through existing list of interactions to see if interaction1 or interaction2 exists
					 * if it the interaction is already present we increment the numberOfOccurence, otherwise the 
					 * interaction is initialized */
					if (InteractionToNumber.containsKey(interaction1)) {
						int numberOfOccurences = InteractionToNumber.get(interaction1);
						InteractionToNumber.put(interaction1, numberOfOccurences + 1);
					} else if (InteractionToNumber.containsKey(interaction2)) {
						int numberOfOccurences = InteractionToNumber.get(interaction2);
						InteractionToNumber.put(interaction2, numberOfOccurences + 1);
					} else {
						InteractionToNumber.put(interaction1, 1); // initialize occurrence to 1
					}
				}
				line = input.readLine(); // read next line
			}
			input.close(); // close BufferedReader
		} catch (Exception e) {
			e.printStackTrace();
		}

		/* Store interactions as objects of type Interaction */
		ArrayList<Interaction> interaction_list = new ArrayList<Interaction>(); // new ArrayList of class Interaction

		for (String inter : InteractionToNumber.keySet()) {

			String[] interaction_identifier = inter.split("\t"); // split interaction identifier by tab

			Interaction inter1 = new Interaction(interaction_identifier[0], interaction_identifier[1], Integer.parseInt(interaction_identifier[2]),
					Integer.parseInt(interaction_identifier[3])); // initialize new interaction object

			interaction_list.add(inter1); // store interaction object in arrayList
		}
		return interaction_list;
	}

	/**
	 * @param foldChangeDataFile file with breast cancer dataset
	 * @param proteinsInNetworkList list of proteins in network
	 * @param testHer2vsTN
	 * @param runOnComputeCanada
	 */
	public static void loadProteinFoldChangeData(String foldChangeDataFile, ArrayList<Protein> proteinsInNetworkList, boolean testHer2vsTN, boolean testHer2vsHR, boolean testHRvsTN, boolean runOnComputeCanada) {
		/* Compute Fold change between Her2 and Triple Negative patients from the index
		 * list of each group and the file path. File contains protein expression.*/

		// InteractionToNumber will contain the number of each possible interaction
		HashMap<String, Double> foldChangeMap = new HashMap<String, Double>();

		try {
			BufferedReader input = null;

			if (runOnComputeCanada) {
				InputStream in = Loader.class.getClassLoader().getResourceAsStream(foldChangeDataFile);
				input = new BufferedReader(new InputStreamReader(in));
			} else {
				input = new BufferedReader(new FileReader(new File(foldChangeDataFile)));
			}

			String line = input.readLine();	// read first line
			line = input.readLine();
			line = input.readLine();
			String[] col = line.split("\t");

			ArrayList<Integer> her2Indexes = new ArrayList<Integer>();
			ArrayList<Integer> tnIndexes = new ArrayList<Integer>();
			ArrayList<Integer> erprIndexes = new ArrayList<Integer>();

			for(int i=0; i<col.length; i++) {
				if(col[i].matches("Her2(.*)")) {
					her2Indexes.add(i);
				} else if(col[i].matches("TN(.*)")) {
					tnIndexes.add(i);
				} else if(col[i].matches("ERPR(.*)")) {
					erprIndexes.add(i);
				}
			}

			// read next line
			line = input.readLine();
			// TODO for testing
			//			int missed_rows = 0;
			//			int total_rows = 0;
			//TODO for testing
			while (line != null) {

				// split new line in string array
				//				total_rows++;
				col = line.split("\t");
				try {
					double averageCondition1 = 1;
					double averageCondition2 = 1;
					
					if(testHer2vsTN) {
						averageCondition1 = Calculator.computeAvg(col, her2Indexes);
						averageCondition2 = Calculator.computeAvg(col, tnIndexes);
					} else if(testHer2vsHR) {
						averageCondition1 = Calculator.computeAvg(col, her2Indexes);
						averageCondition2 = Calculator.computeAvg(col, erprIndexes);
					} else if(testHRvsTN) {
						averageCondition1 = Calculator.computeAvg(col, erprIndexes);
						averageCondition2 = Calculator.computeAvg(col, tnIndexes);
					}
					
					// expression = fold change for Her2
					double foldChange = averageCondition1 / averageCondition2;

					// store protein name(col[0]) and fold change (expression ratio) in HashMap
					// (FoldChange)

					if(!Double.isNaN(foldChange)){
						if(col.length > 132 && col[132] != null) {
							foldChangeMap.put(col[132], foldChange);
						}
					}


					// next line
				} catch (IllegalStateException ex){
					//					missed_rows++;
				}
				line = input.readLine();
			}
			input.close();
			//System.out.println(String.format("Total rows: %s, missed rows: %s", total_rows, missed_rows));

		} catch (Exception e) {
			e.printStackTrace();
		}

		//int countprot = 0; // initialize counter to check for proteins found in the network from the project data

		/* Extract info from fold change map to be stored in Protein object */
		for (int i = 0; i < proteinsInNetworkList.size(); i++) {	

			String protein_name = proteinsInNetworkList.get(i).getProteinName(); // get protein name

			// if protein is found in the fold change map, store fold change in Protein object
			if (foldChangeMap.containsKey(protein_name)) {
				double fold_change = foldChangeMap.get(protein_name);
				proteinsInNetworkList.get(i).setFoldChange(fold_change);
				//countprot++;
			}
		}
		/*	System.out.println("proteins from project found in network " +countprot +"; total proteins in project data " +foldChangeMap.size()
		+"; represents " +(countprot/(double)foldChangeMap.size())); */
	}

	public static double[][] loadDistanceMatrix(String distance_matrixFile, ArrayList<Protein> ProteinList, boolean runComputeCanada) {
		/* Import distance matrix from text file, it's dimensions are based on the size of the 
		 * proteinNetwork List initially used to build the distance Matrix file */

		double[][] distanceMatrix = new double[ProteinList.size()][ProteinList.size()]; // Initialize distance matrix

		try {

			BufferedReader input;

			if(runComputeCanada) {
				InputStream in = Loader.class.getClassLoader().getResourceAsStream(distance_matrixFile);
				input = new BufferedReader(new InputStreamReader(in));
			} else {
				input = new BufferedReader(new FileReader(new File(distance_matrixFile)));
			}

			String line = input.readLine(); // read first line

			int x = 0; // initialize row counter

			while (line != null) {

				String[] col = line.split("\t"); // split columns
				int y = 0; // initialize column counter (resets at the end of every row)

				for (String str : col) { // str is the index to go through all element of col

					distanceMatrix[x][y] = Double.parseDouble(str);;// set the value (distance) at the appropriate coordinates
					y++; // adds one to the value of y (change column)
				}
				x++; // adds one to the vale of x (change row)
				line = input.readLine(); // read next line

			}
			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return distanceMatrix;
	} // end import distance matrix

	public static ArrayList<Annotation> importAnnotationGo(String inputFile, ArrayList<Protein> networkProteinList, boolean runComputeCanada) {
		/* Maps protein IDs (values) in the form of a string array with it's
		 * corresponding GO ID (key) from given annotation GO file. Only keeps proteins found in the network
		 */

		// GoTerms will contain all GO terms with a list of associated protein IDs
		ArrayList<Annotation> goTermsList = new ArrayList<Annotation>();

		try {
			BufferedReader input;

			if(runComputeCanada) {
				InputStream in = Loader.class.getClassLoader().getResourceAsStream(inputFile);
				input = new BufferedReader(new InputStreamReader(in));
			} else {
				input = new BufferedReader(new FileReader(new File(inputFile)));
			}

			String line = input.readLine(); // Skip header
			line = input.readLine(); // read first line
			// int goCount = 0;
			while (line != null) {

				// System.out.println("go term = " +goCount);
				String[] col = line.split("\t");

				String term = col[0]; // first column holds GO terms (index = 0)
				String protein_ids = col[6]; // 7th column holds gene IDs (proteins IDs) (index = 6)
				String protein_symbols = col[7]; // 8th column holds gene symbols (protein symbols) (index = 7)

				if(!term.startsWith("GO")) {
					System.out.println("hi");
				}

				String[] proteinIDList = protein_ids.split("\\|"); // split for names of individual proteins
				String[] proteinSymbolsList = protein_symbols.split("\\|");

				/* Check if proteins associated to the goTerm are found in the network. Only
				 * proteins found in the network will be saved. Also store index of proteins
				 * within the network. */
				ArrayList<String> goProteinsInNetworkList = new ArrayList<String>();
				ArrayList<String> goSymbolsInNetworkList = new ArrayList<String>();
				ArrayList<Integer> idxProteinsInNetworkList = new ArrayList<Integer>();

				for (int i = 0; i < proteinIDList.length; i++) { // search for all IDs associated to the GO term

					boolean foundProt = false; // initialize that protein is not found in the network
					int protCounter = 0; // indicates what Protein Object we are looking at

					while (!foundProt && (protCounter < networkProteinList.size())) {

						if (networkProteinList.get(protCounter).getProteinId() == Integer.parseInt(proteinIDList[i])) {
							goProteinsInNetworkList.add(proteinIDList[i]); // save protein ID number
							goSymbolsInNetworkList.add(proteinSymbolsList[i]);
							idxProteinsInNetworkList.add(protCounter); // save index of protein in network
							foundProt = true; // mark as protein found
						} else {
							protCounter++; // increase counter
						}
					}
				}

				/* Check prevalence of proteins from annotation GO in network */
				Double goPrevalence = goProteinsInNetworkList.size() / (double) (proteinIDList.length);

				if (goPrevalence > 0.5 && goProteinsInNetworkList.size() > 2 && goProteinsInNetworkList.size() < 1000) {

					Annotation goAnnotation = new Annotation(term, goProteinsInNetworkList, goSymbolsInNetworkList, idxProteinsInNetworkList);

					goTermsList.add(goAnnotation);

				}
				// goCount++;
				line = input.readLine(); // next line
			}

			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return goTermsList;
	}
	/**************************************************************************************************
	 * Load a Monte Carlos distribution file for all possible TPD generated on compute canada. 
	 * File format : [header] TPD (n=#) "\t" Freq/NotFound
	 * 
	 * Use data to compute p-values.
	 * 
	 * 
	 **************************************************************************************************/
	public static void loadDistributions(String distributionInputFile, int go_start, int go_stop, String normalDistributionParametersFile) {


		try {
			BufferedReader input = new BufferedReader(new FileReader(new File(distributionInputFile)));

			String line = input.readLine();
			HashMap<Double, Double> distributionMap = new HashMap<Double, Double>(); /* Store distribution information for nProt in a hash map */

			double[][] normalDistributionParameters = new double[go_stop][2];

			String[] col = line.split("\\s+");
			String[] splitTPD = col[3].split("\\)");


			int nProt = Integer.parseInt(splitTPD[0]);
			distributionMap = new HashMap<Double, Double>();

			line = input.readLine();

			while(line!=null) {
				if (!line.equals("")) {

					col = line.split("\\s+");

					if (col[0].equals("TPD") && col[4].equals("Frequency")) {

						if (nProt != 1) {
							normalDistributionParameters[nProt][0] = NormalApproximation.computeDistributionMean(distributionMap);

							normalDistributionParameters[nProt][1] = NormalApproximation.computeDistributionStandardDeviation(distributionMap, normalDistributionParameters[nProt][0]);

						}
						/* Reset hashMap */
						distributionMap = new HashMap<Double, Double>();

						/* Reset protein number */
						splitTPD = col[3].split("\\)");
						nProt = Integer.parseInt(splitTPD[0]);
						if(nProt==701) {
							System.out.println("nProt=701");
						}

						//System.out.println("nProt = " + nProt);


					} else if (col[0].equals("TPD") && col[4].equals("Not")) {
						line = input.readLine();
						continue;

					} else {
						/* Import Frequency map */

						double tpd = Double.parseDouble(col[0]);
						double freq = Double.parseDouble(col[1]);


						distributionMap.put(tpd, freq);
					}
				}
				line = input.readLine();
			}	

			NormalApproximation.exportNormalDistributionParameters(normalDistributionParameters, normalDistributionParametersFile, go_start);

			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	/**
	 * Load a Monte Carlos distribution file for all possible TPD generated on compute canada. 
	 * File format : [header] TPD (n=#) "\t" Freq/NotFound
	 * The mean and standard deviation are computed from the distributions and outputed to a file. 
	 * 
	 * => output of file should be separate from this method
	 * 
	 * @param distributionInputFile				File path for distributions 
	 * @param nProtToSampleLowerBound			Lower limit of proteins that were sampled
	 * @param nProtToSampleUpperBound			Upper limit of proteins that were sampled
	 * @param normalDistributionParametersFile	
	 */
	public static double[][] loadMonteCarloDistributions(String distributionInputFile, int nProtToSampleLowerBound, int nProtToSampleUpperBound, String normalDistributionParametersFile) {

		/* Initialize array to contain normal distribution parameters {mean, standard deviation} */
		double[][] normalDistributionParameters = new double[nProtToSampleUpperBound-nProtToSampleLowerBound][2];
		int countLine =1;
		try {
			BufferedReader input = new BufferedReader(new FileReader(new File(distributionInputFile)));
			HashMap<Double, Double> distributionMap = new HashMap<Double, Double>();; /* Store distribution information for nProt in a hash map */
			String line = input.readLine();
			int nProt = 0;
			while(line!=null) {


				String[] col = line.split("\\s+");



				if(line.contains("TPD")) {
					String[] splitTPD = col[3].split("\\)");
					nProt = Integer.parseInt(splitTPD[0]);
					distributionMap = new HashMap<Double, Double>();

					// NEED to take into account that some distributions aren't computed and will have a frequency of 0.
				}else if(line.isEmpty()) {
					/* if line is empty; confirm that distribution totals 1; compute mean and standard deviation of distribution= */
					boolean frequenciesEqualOne = Loader.verifySumOfFrequencies(distributionMap);

					if(frequenciesEqualOne) {
						normalDistributionParameters[nProt-nProtToSampleLowerBound][0] = NormalApproximation.computeDistributionMean(distributionMap);
						normalDistributionParameters[nProt-nProtToSampleLowerBound][1] = NormalApproximation.computeDistributionStandardDeviation(distributionMap, normalDistributionParameters[nProt-nProtToSampleLowerBound][0]);
					}
				} else {
					/* Import Frequency map */
					
					double tpd = Double.parseDouble(col[0]);
					double freq = Double.parseDouble(col[1]);


					distributionMap.put(tpd, freq);
				}
				line = input.readLine();
				countLine++;
	
				
				
			}	

			NormalApproximation.exportNormalDistributionParameters(normalDistributionParameters, normalDistributionParametersFile, nProtToSampleLowerBound);

			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

		return normalDistributionParameters;
	}

	public static double setAnnotationPvaluesFromMonteCarloDistribution(String distributionInputFile, int nProtToSampleLowerBound, int nProtToSampleUpperBound, ArrayList<Annotation> annotationList, ArrayList<Annotation> shuffledList) {
		double minimum_pValue = 1;	
		try {
			BufferedReader input = new BufferedReader(new FileReader(new File(distributionInputFile)));
			HashMap<Double, Double> distributionMap = new HashMap<Double, Double>();; /* Store distribution information for nProt in a hash map */
			String line = input.readLine();
			int nProt = 0;


			while(line!=null) {


				String[] col = line.split("\\s+");



				if(line.contains("TPD")) {
					String[] splitTPD = col[3].split("\\)");
					nProt = Integer.parseInt(splitTPD[0]);
					distributionMap = new HashMap<Double, Double>();

					// NEED to take into account that some distributions aren't computed and will have a frequency of 0.
				}else if(line.isEmpty()) {
					/* if line is empty; confirm that distribution totals 1; compute mean and standard deviation of distribution= */
					boolean frequenciesEqualOne = Loader.verifySumOfFrequencies(distributionMap);

					if(frequenciesEqualOne) {

						for(int i = 0; i<annotationList.size(); i++) {
							//			System.out.println("hello3");
							if(annotationList.get(i).getNumberOfProteins() == nProt) {
								//										System.out.println("hello4");
								double tpd = annotationList.get(i).getTPD();
								double pval = Calculator.compute_pValue(distributionMap, tpd);

								annotationList.get(i).set_pValue(pval);
								//											System.out.println("moyenne = " + nd.getMean() + "	std = " + nd.getStandardDeviation());
								//											System.out.println("pvalue = " + p_val + "	TPD = " + tpd);

								if(pval < minimum_pValue && pval != 0 ) {
									minimum_pValue = pval;
								} 
							}
							if(shuffledList.get(i).getNumberOfProteins() == nProt) {
								//										System.out.println("hello4");
								double tpd2 = shuffledList.get(i).getTPD();
								double pval2 = Calculator.compute_pValue(distributionMap, tpd2);

								shuffledList.get(i).set_pValue(pval2);
								//											System.out.println("moyenne = " + nd.getMean() + "	std = " + nd.getStandardDeviation());
								//											System.out.println("pvalue = " + p_val + "	TPD = " + tpd);

								if(pval2 < minimum_pValue && pval2 != 0 ) {
									minimum_pValue = pval2;
								}

							} 
						}
					} else {
						System.out.println("Nprot " + nProt + " does NOT exist" );
					}

				} else {
					/* Import Frequency map */

					double tpd = Double.parseDouble(col[0]);
					double freq = Double.parseDouble(col[1]);


					distributionMap.put(tpd, freq);
				}
				line = input.readLine();
			}	


			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return minimum_pValue;
	}

	/**
	 * Compute the sum of frequencies for a given distribution, if it isn't equal to 1, make error. 
	 * 
	 * @param distributionMap	
	 * @return frequenciesVerified	 boolean; if true, sum of frequencies is 1
	 */
	public static boolean verifySumOfFrequencies(HashMap<Double, Double> distributionMap) { 

		double countFrequencies = 0;
		boolean frequenciesVerified = true;

		/* compute the sum of frequencies in the distribution */
		for(Double tpd:distributionMap.keySet()) {
			countFrequencies += distributionMap.get(tpd);
		}

		/* check if sum of frequencies equals 1 */
		if(Math.rint(countFrequencies)!=1) {
			//System.out.println("ERROR loading distribution; frequency = " + countFrequencies);
			frequenciesVerified = false;
		}
		return frequenciesVerified;
	}

	public static void loadShuffledGoAnnotations(ArrayList<Annotation> shuffledAnnotationGoList, String shuffled_annotationGo_inputFile) {

		try {
			//InputStream in = Loader.class.getClassLoader().getResourceAsStream(shuffled_annotationGo_inputFile);
			//BufferedReader input = new BufferedReader(new InputStreamReader(in));

			BufferedReader input = new BufferedReader(new FileReader(new File(shuffled_annotationGo_inputFile)));

			String line;

			line = input.readLine(); // read header
			line = input.readLine();

			while(line != null) {

				String[] col = line.split("\t");	
				String[] proteinIDs = col[1].split("\\|");
				String[] proteinIdxs = col[2].split("\\|");

				//convert proteinIDs String[] into an arrayList
				ArrayList<String> proteinIDsList = new ArrayList<String>(Arrays.asList(proteinIDs)); 

				//convert proteinIDxs String[] into ArrayList<Integer>
				ArrayList<Integer> proteinIdxList = General.convertStringArrayToIntArrayList(proteinIdxs);

				//Add go annotation to Cluster List
				Annotation shuffledGo = new Annotation(col[0], proteinIDsList, proteinIdxList);
				shuffledAnnotationGoList.add(shuffledGo);

				line = input.readLine();
			}


			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	/**
	 * Loads frequency distribution for a certain sampling number.
	 *
	 * @param fileName
	 * @return hash map, where key is tpd and value is frequency
	 */
	public static HashMap<Double, Double> loadDistributionByProteinsNumber(String fileName) {
		HashMap<Double, Double> pairs = new HashMap<>();
		try {
			BufferedReader input = new BufferedReader(new FileReader(fileName));
			input.readLine();
			String line = input.readLine();
			while (line != null) {
				line = line.trim();
				double tpd = Double.parseDouble(line.split("\t")[0]);
				double frequency = Double.parseDouble(line.split("\t")[1]);
				pairs.put(tpd, frequency);
				line = input.readLine();
			}

			input.close();
		} catch (IOException ex) {
			System.out.println("The file with distribution was not found.");
		}
		return pairs;
	}


	public static void loadDistributions2(String file, int lowerBound, int upperBound, String distParamsFile){
		ArrayList<String> lines = new ArrayList<String>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = reader.readLine();
			while (line != null){
				lines.add(line);
				line = reader.readLine();
			}
			lines.add(line);
			reader.close();
		} catch (IOException ex){
			System.out.println("File was not found");
		}
		HashMap<Double, Double> distributionMap = new HashMap<>();
		double[][] parameters = new double[upperBound-lowerBound][2];
		int currentProteinsAmount = lowerBound;
		for (String line:lines){
			if (line.contains("=")) {
				parameters[currentProteinsAmount-lowerBound][0] = NormalApproximation.computeDistributionMean(distributionMap);
				parameters[currentProteinsAmount-lowerBound][1] = NormalApproximation.computeDistributionStandardDeviation(distributionMap, parameters[currentProteinsAmount-lowerBound][0]);
				distributionMap = new HashMap<>();
				currentProteinsAmount = Integer.parseInt(line.substring(line.indexOf('='), line.indexOf(')')));
			} else if (line.contains("\t") && line.split("\t").length == 2){
				Double key = Double.parseDouble(line.split("\t")[0]);
				Double value = Double.parseDouble(line.split("\t")[1]);
				distributionMap.put(key, value);
			}
		}
		NormalApproximation.exportNormalDistributionParameters(parameters, distParamsFile, lowerBound);
	}

}
