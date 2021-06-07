package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import graph.Interaction;

public class NetworkInteractionsLoader {

	/**
	 * Loads the PPI repository as a HashMap<String, Integer>, to identify unique PPIs. 
	 * - BioGRID : keeps only human to human interactions
	 * - String : first loads ensemble to entrez Id map
	 * Removes proteins that contribute to > 1000 interactions
	 * Returns Interactions in network as an ArrayList of Interaction Objects
	 * 
	 * @param inputRepositoryFile				String - file path for PPI repository 
	 * @param networkType						int - specify network type: 0 = BioGRID, 1 = STRING
	 * @param ensembleToEntrezIdMapFile			String - file path for ensemble to entrez Id map *required for STRING*
	 * @param removeOverConnectedProteins		boolean - specify to remove overly connected proteins
	 * @param numberOfExcessiveInteractions		int - number of maximum interactions a protein can take part in
	 * 
	 * @return interactionList 		ArrayList<Interaction> - formated interaction list
	 */
	public static ArrayList<Interaction> importInteractionNetwork(String inputRepositoryFile, int networkType, String ensembleToEntrezIdMapFile, boolean removeOverConnectedProteins, int numberOfExcessiveInteractions, String taxonomyID, int combinedScore){
		
		HashMap<String, Integer> interactionToNumber = new HashMap<>(); 

		switch(networkType) {
		case 0 : 	// BioGRID
			boolean checkTaxID = false;
			if(!taxonomyID.isEmpty()) {
				checkTaxID = true;
			}
			if(checkTaxID) {
				interactionToNumber = loadBioGRIDinteractionsWithTaxId(inputRepositoryFile, taxonomyID);		
			} else {
				interactionToNumber = loadBioGRIDinteractions(inputRepositoryFile);
			}
			break;
		case 1 : 	// STRING links-detailed
			HashMap<String, String> ensembleToEntrezIdMap = loadStringMap(ensembleToEntrezIdMapFile); 
			interactionToNumber = loadDetailedStringNetwork(inputRepositoryFile, ensembleToEntrezIdMap, combinedScore);
			break;
		case 2 : // STRING links / physical
			HashMap<String, String> ensembleToEntrezIdMap2 = loadStringMap(ensembleToEntrezIdMapFile);
			interactionToNumber = loadStringNetwork(inputRepositoryFile, ensembleToEntrezIdMap2, combinedScore);
		}
		
		if(removeOverConnectedProteins) {
			interactionToNumber = removeOverConnectedProteins(interactionToNumber, numberOfExcessiveInteractions);
		}
		ArrayList<Interaction> interactionList = storeInteractionList(interactionToNumber);
		return interactionList;  
	}

	/**
	 * Read BioGRID repository (from inputFile) to extract ALL names of possible interactions and their number of re-occurrence.
	 * Ensures that there are no repeat interactions and that all proteins involved in interactions are human (homo sapiens; 9606) 
	 */
	public static HashMap<String, Integer> loadBioGRIDinteractions(String inputRepositoryFile) {

		HashMap<String, Integer> interactionToNumber = new HashMap<String, Integer>();// To contain interaction name and ID, and number of occurence
		try {
			/* Read BioGrid ; get all possible human protein-protein interactions */

			InputStream in = new FileInputStream(new File(inputRepositoryFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // read first line
			line = input.readLine(); // read second line

			while (line != null) { // stops when there are no more lines in file (in)

				String[] col = line.split("\t"); // split line by tabs; obtain individual columns

				/* Do not include protein interactions that involve non human proteins. 
				 * The code for homo sapiens is 9606 */

					String interactor1 = col[7]; // get name of prot1
					String interactor2 = col[8]; // get name of prot2

					String interactorID1 = col[1]; // get ID of prot1
					String interactorID2 = col[2]; // get ID of prot2

					/* establish the two possible interaction combinations */
					String interaction1 = interactor1 + "\t" + interactor2 + "\t" + interactorID1 + "\t"
							+ interactorID2;
					String interaction2 = interactor2 + "\t" + interactor1 + "\t" + interactorID2 + "\t"
							+ interactorID1;

					/* Look through existing list of interactions to see if interaction1 or interaction2 exists
					 * if it the interaction is already present we increment the numberOfOccurence, otherwise the 
					 * interaction is initialized */
					if (interactionToNumber.containsKey(interaction1)) {
						int numberOfOccurences = interactionToNumber.get(interaction1);
						interactionToNumber.put(interaction1, numberOfOccurences + 1);
					} else if (interactionToNumber.containsKey(interaction2)) {
						int numberOfOccurences = interactionToNumber.get(interaction2);
						interactionToNumber.put(interaction2, numberOfOccurences + 1);
					} else {
						interactionToNumber.put(interaction1, 1); // initialize occurrence to 1
					}
				}
				line = input.readLine(); // read next line
			
			input.close(); // close BufferedReader
		} catch (Exception e) {
			e.printStackTrace();
		}
		return interactionToNumber;
	}

	/**
	 * Read BioGRID repository (from inputFile) to extract ALL names of possible interactions and their number of re-occurrence.
	 * Ensures that there are no repeat interactions and that all proteins involved in interactions are human (homo sapiens; 9606) 
	 */
	public static HashMap<String, Integer> loadBioGRIDinteractionsWithTaxId(String inputRepositoryFile, String taxonomyID) {

		HashMap<String, Integer> interactionToNumber = new HashMap<String, Integer>();// To contain interaction name and ID, and number of occurence
		try {
			/* Read BioGrid ; get all possible human protein-protein interactions */

			InputStream in = new FileInputStream(new File(inputRepositoryFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // read first line
			line = input.readLine(); // read second line
			
			int taxId = Integer.parseInt(taxonomyID);
			while (line != null) { // stops when there are no more lines in file (in)

				String[] col = line.split("\t"); // split line by tabs; obtain individual columns

				/* Do not include protein interactions that involve non human proteins. 
				 * The code for homo sapiens is 9606 */

				int species1 = Integer.parseInt(col[15]); // get species of prot1 
				int species2 = Integer.parseInt(col[16]); // get species of prot2

				if (species1 == taxId && species2 == taxId) {

					String interactor1 = col[7]; // get name of prot1
					String interactor2 = col[8]; // get name of prot2

					String interactorID1 = col[1]; // get ID of prot1
					String interactorID2 = col[2]; // get ID of prot2

					/* establish the two possible interaction combinations */
					String interaction1 = interactor1 + "\t" + interactor2 + "\t" + interactorID1 + "\t"
							+ interactorID2;
					String interaction2 = interactor2 + "\t" + interactor1 + "\t" + interactorID2 + "\t"
							+ interactorID1;

					/* Look through existing list of interactions to see if interaction1 or interaction2 exists
					 * if it the interaction is already present we increment the numberOfOccurence, otherwise the 
					 * interaction is initialized */
					if (interactionToNumber.containsKey(interaction1)) {
						int numberOfOccurences = interactionToNumber.get(interaction1);
						interactionToNumber.put(interaction1, numberOfOccurences + 1);
					} else if (interactionToNumber.containsKey(interaction2)) {
						int numberOfOccurences = interactionToNumber.get(interaction2);
						interactionToNumber.put(interaction2, numberOfOccurences + 1);
					} else {
						interactionToNumber.put(interaction1, 1); // initialize occurrence to 1
					}
				}
				line = input.readLine(); // read next line
			}
			input.close(); // close BufferedReader
		} catch (Exception e) {
			e.printStackTrace();
		}
		return interactionToNumber;
	}

	
	private static HashMap<String, Integer> loadDetailedStringNetwork(String inputRepositoryFile, HashMap<String, String> ensembleToEntrezIdMap, int combinedScore){
		HashMap<String, Integer> interactionToNumber = new HashMap<String, Integer>();// To contain interaction name and ID, and number of occurence
		try {
			/* Read String ; get all possible human protein-protein interactions */

			InputStream in = new FileInputStream(new File(inputRepositoryFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // read first line
			line = input.readLine(); // read second line

			//int countLine = 2;
			while (line != null) { // stops when there are no more lines in file (in)
				//System.out.println("line : " + countLine);
				String[] col = line.split("\\s+"); // split line by tabs; obtain individual columns

				/* Check interaction conditions: needs experimental or database validation and confidence score >= 0.4 */

				/*int exp = Integer.parseInt(col[6]);
				int dat = Integer.parseInt(col[7]);
				int score = Integer.parseInt(col[9]);*/

				//if((exp >= 1 || dat >=1) && score >= 400) {
				if (Integer.parseInt(col[9]) >= combinedScore) {
					String interactor1 = col[0]; // get ensemble ID of prot1
					String interactor2 = col[1]; // get ensemble ID of prot2

					String[] interactorID1 = null, interactorID2 = null;
					String entrezID1, entrezID2 = null;
					String geneName1, geneName2 = null;
					
					if(ensembleToEntrezIdMap.containsKey(interactor1)) {
						interactorID1 = ensembleToEntrezIdMap.get(interactor1).split("\t"); // get ID of prot1
						geneName1 = interactorID1[0];
						entrezID1 = interactorID1[1];
					} else {
						geneName1 = interactor1;
						entrezID1 = String.valueOf(Integer.MAX_VALUE);
					}
					
					if(ensembleToEntrezIdMap.containsKey(interactor2)) {
						interactorID2 = ensembleToEntrezIdMap.get(interactor2).split("\t"); // get ID of prot2
						geneName2 = interactorID2[0];
						entrezID2 = interactorID2[1];
					} else {
						geneName2 = interactor2;
						entrezID2 = String.valueOf(Integer.MAX_VALUE);
					}

					/* establish the two possible interaction combinations */
					String interaction1 = geneName1 + "\t" + geneName2 + "\t" + entrezID1 + "\t"
							+ entrezID2;
					String interaction2 = geneName2 + "\t" + geneName1 + "\t" + entrezID2 + "\t"
							+ entrezID1;
					/* Look through existing list of interactions to see if interaction1 or interaction2 exists
					 * if it the interaction is already present we increment the numberOfOccurence, otherwise the 
					 * interaction is initialized */
					if (interactionToNumber.containsKey(interaction1)) {
						int numberOfOccurences = interactionToNumber.get(interaction1);
						interactionToNumber.put(interaction1, numberOfOccurences + 1);
					} else if (interactionToNumber.containsKey(interaction2)) {
						int numberOfOccurences = interactionToNumber.get(interaction2);
						interactionToNumber.put(interaction2, numberOfOccurences + 1);
					} else {
						interactionToNumber.put(interaction1, 1); // initialize occurrence to 1
					}
				}

				//}
				//countLine++;
				line = input.readLine(); // read next line
			}
			input.close(); // close BufferedReader

		} catch (Exception e) {
			e.printStackTrace();
		}


		return interactionToNumber;	
	}

	private static HashMap<String, Integer> loadStringNetwork(String inputRepositoryFile, HashMap<String, String> ensembleToEntrezIdMap, int combinedScore){
		HashMap<String, Integer> interactionToNumber = new HashMap<String, Integer>();// To contain interaction name and ID, and number of occurence
		try {
			/* Read String ; get all possible human protein-protein interactions */

			InputStream in = new FileInputStream(new File(inputRepositoryFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // read first line
			line = input.readLine(); // read second line

			//int countLine = 2;
			while (line != null) { // stops when there are no more lines in file (in)
				//System.out.println("line : " + countLine);
				String[] col = line.split("\\s+"); // split line by tabs; obtain individual columns

				/* Check interaction conditions: needs experimental or database validation and confidence score >= 0.4 */

				/*int exp = Integer.parseInt(col[6]);
				int dat = Integer.parseInt(col[7]);
				int score = Integer.parseInt(col[9]);*/

				//if((exp >= 1 || dat >=1) && score >= 400) {
				if (Integer.parseInt(col[2]) >= combinedScore) {
					String interactor1 = col[0]; // get string ID of prot1
					String interactor2 = col[1]; // get string ID of prot2

					String[] interactorID1 = null, interactorID2 = null;
					String entrezID1, entrezID2 = null;
					String geneName1, geneName2 = null;
					
					if(ensembleToEntrezIdMap.containsKey(interactor1)) {
						interactorID1 = ensembleToEntrezIdMap.get(interactor1).split("\t"); // get ID of prot1
						geneName1 = interactorID1[0];
						entrezID1 = interactorID1[1];
					} else {
						geneName1 = interactor1;
						entrezID1 = String.valueOf(Integer.MAX_VALUE);
					}
					
					if(ensembleToEntrezIdMap.containsKey(interactor2)) {
						interactorID2 = ensembleToEntrezIdMap.get(interactor2).split("\t"); // get ID of prot2
						geneName2 = interactorID2[0];
						entrezID2 = interactorID2[1];
					} else {
						geneName2 = interactor2;
						entrezID2 = String.valueOf(Integer.MAX_VALUE);
					}

					/* establish the two possible interaction combinations */
					String interaction1 = geneName1 + "\t" + geneName2 + "\t" + entrezID1 + "\t"
							+ entrezID2;
					String interaction2 = geneName2 + "\t" + geneName1 + "\t" + entrezID2 + "\t"
							+ entrezID1;

					/* Look through existing list of interactions to see if interaction1 or interaction2 exists
					 * if it the interaction is already present we increment the numberOfOccurence, otherwise the 
					 * interaction is initialized */
					if (interactionToNumber.containsKey(interaction1)) {
						int numberOfOccurences = interactionToNumber.get(interaction1);
						interactionToNumber.put(interaction1, numberOfOccurences + 1);
					} else if (interactionToNumber.containsKey(interaction2)) {
						int numberOfOccurences = interactionToNumber.get(interaction2);
						interactionToNumber.put(interaction2, numberOfOccurences + 1);
					} else {
						interactionToNumber.put(interaction1, 1); // initialize occurrence to 1
					}
				}

				//}
				//countLine++;
				line = input.readLine(); // read next line
			}
			input.close(); // close BufferedReader

		} catch (Exception e) {
			e.printStackTrace();
		}


		return interactionToNumber;	
	}
	
	public static HashMap<String, String> loadStringMap(String inputFile){
		HashMap<String, String> ensembleToEntrezIdMap = new HashMap<>();

		try {

			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); //header
			line = input.readLine();

			while(line!= null) {
				String[] col = line.split("\t");
				String stringID, geneSymbol, entrezId = null;
				if(col[1] == null) {
					stringID = "NA";
				}else {
					stringID = col[1];
				}
				if(col[0] == null) {
					geneSymbol = "NA";
				}else {
					geneSymbol = col[0];
				}
				if(col.length<3) {
					entrezId = String.valueOf(Integer.MAX_VALUE);
				} else { 
					if(col[2].equals("NA")) {
						entrezId = String.valueOf(Integer.MAX_VALUE);
					}else {
						entrezId = col[2];
					}
					
				}
				ensembleToEntrezIdMap.put(stringID, geneSymbol + "\t" + entrezId); // Ensemble ID; Entrez ID \t GeneName
				line = input.readLine();
			}

			input.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return ensembleToEntrezIdMap;
	}

	/**
	 * @param interactionToNumber
	 * @param numberOfExcessiveInteractions
	 * @return newInteractionToNumberMap
	 * 
	 * Takes initial interaction to number map and removes interactions that contain proteins that are involved in excess interactions (>1000)
	 */
	public static HashMap<String, Integer> removeOverConnectedProteins(HashMap<String, Integer> interactionToNumber, int numberOfExcessiveInteractions){

		HashMap<String, Integer> interactionToNumber2 = new HashMap<String, Integer>(); // new interaction list without proteins that are overly connected
		HashMap<String, Integer> proteinToNumber = new HashMap<String, Integer>(); // stores number of interactions proteins are involved in

		HashSet<String> proteinToRemove = new HashSet<>();
		/* Count number of interactions each protein is part of */
		//int countProtsToRemove = 0;
		for(String inter : interactionToNumber.keySet()) {

			String[] interaction_identifier = inter.split("\t"); // split interaction identifier by tab
			String prot1 = interaction_identifier[0];
			String prot2 = interaction_identifier[1];

			/* Check interactions for prot1 */
			
			if(proteinToNumber.containsKey(prot1)) {
				int nInteractionsProt1 = proteinToNumber.get(prot1);
				proteinToNumber.put(prot1, nInteractionsProt1+1);
				if(nInteractionsProt1 == 999) {
					proteinToRemove.add(prot1);
				}
			} else {
				proteinToNumber.put(prot1, 1);
			}

			/* Check interactions for prot2 */
			if(proteinToNumber.containsKey(prot2)) {
				int nInteractionsProt2 = proteinToNumber.get(prot2);
				proteinToNumber.put(prot2, nInteractionsProt2+1);
				if(nInteractionsProt2 == 999) {
					proteinToRemove.add(prot2);
				}
			} else {
				proteinToNumber.put(prot2, 1);
			}
		}
		
		/* Remove interactions which contain overly connected proteins */ 
		int interactionsToRemove = 0;
		for(String inter : interactionToNumber.keySet()) {

			String[] interaction_identifier = inter.split("\t"); // split interaction identifier by tab
			String prot1 = interaction_identifier[0];
			String prot2 = interaction_identifier[1];

			/* Both interacting proteins must be present in less than 1000 interactions */
			if(proteinToNumber.get(prot1) < numberOfExcessiveInteractions && proteinToNumber.get(prot2) < numberOfExcessiveInteractions) {
				interactionToNumber2.put(inter, interactionToNumber.get(inter));
			} else {
				interactionsToRemove++;
			}
		}
		System.out.println("Removing overly connected proteins: " + proteinToRemove.size() + "; removed interactions: " + interactionsToRemove); 
		
//		try {
//			BufferedWriter out = new BufferedWriter(new FileWriter(new File ("C://Users//Rachel//Documents//PIGNON//proteinsInStringNetworkRemoved.tsv")));
//
//			for(String id: proteinToRemove) {
//				out.write(id + "\n");
//				out.flush();
//			}
//
//			out.close();		
//		}catch (Exception e) {
//			e.printStackTrace();
//		} 

		return interactionToNumber2;
	}

	/**
	 * @param interactionToNumber
	 * @return interaction_list 
	 *
	 * Takes interaction to number map and stores elements as Interactions Objects in an ArrayList
	 */
	public static ArrayList<Interaction> storeInteractionList(HashMap<String, Integer> interactionToNumber){

		/* Store interactions as objects of type Interaction */
		ArrayList<Interaction> interaction_list = new ArrayList<Interaction>(); // new ArrayList of class Interaction

		for (String inter : interactionToNumber.keySet()) {

			String[] interaction_identifier = inter.split("\t"); // split interaction identifier by tab


			Interaction inter1 = new Interaction(interaction_identifier[0], interaction_identifier[1], Integer.parseInt(interaction_identifier[2]),
					Integer.parseInt(interaction_identifier[3])); // initialize new interaction object

			interaction_list.add(inter1); // store interaction object in arrayList
		}

		return interaction_list;
	}


}
