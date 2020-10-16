package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

import graph.Interaction;

public class NetworkInteractionsLoader {

	/**
	 * Inputs the BioGRID repository for human interactions, removes interactions that contain proteins that are overly connected (>1000 interactions)
	 * Returns Interactions in network as an ArrayList of Interaction Objects
	 * 
	 * @param inputRepositoryFile
	 * @param runOnComputeCanda
	 * @param numberOfExcessiveInteractions
	 * @return interactionList 
	 */
	public static ArrayList<Interaction> importInteractionNetwork(String inputRepositoryFile, boolean runOnComputeCanda, boolean removeOverConnectedProteins, int numberOfExcessiveInteractions){
		
		HashMap<String, Integer> interactionToNumber = loadInteractionRepository(inputRepositoryFile, runOnComputeCanda);
		if(removeOverConnectedProteins) {
			interactionToNumber = removeOverConnectedProteins(interactionToNumber, numberOfExcessiveInteractions);
		}
		ArrayList<Interaction> interactionList = storeInteractionList(interactionToNumber);
		
		return interactionList;  
	}
	
	/**
	 * Read BioGRID repository (from inputFile) to extract ALL names of possible interactions and their number of re-occurrence.
	 * Ensures that there are no repeat interactions and that all proteins involved in interactions are human (homo sapiens; 9606) 
	 **/
	public static HashMap<String, Integer> loadInteractionRepository(String inputRepositoryFile, boolean runOnComputeCanada) {
		
		HashMap<String, Integer> interactionToNumber = new HashMap<String, Integer>();// To contain interaction name and ID, and number of occurence
		
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
		
		/* Count number of interactions each protein is part of */
		for(String inter : interactionToNumber.keySet()) {
			
			String[] interaction_identifier = inter.split("\t"); // split interaction identifier by tab
			String prot1 = interaction_identifier[0];
			String prot2 = interaction_identifier[1];
			
			/* Check interactions for prot1 */
			if(proteinToNumber.containsKey(prot1)) {
				int nInteractionsProt1 = proteinToNumber.get(prot1);
				proteinToNumber.put(prot1, nInteractionsProt1+1);
			} else {
				proteinToNumber.put(prot1, 1);
			}
			
			/* Check interactions for prot2 */
			if(proteinToNumber.containsKey(prot2)) {
				int nInteractionsProt2 = proteinToNumber.get(prot2);
				proteinToNumber.put(prot2, nInteractionsProt2+1);
			} else {
				proteinToNumber.put(prot2, 1);
			}
		}
		
		/* Remove interactions which contain overly connected proteins */ 
		for(String inter : interactionToNumber.keySet()) {

			String[] interaction_identifier = inter.split("\t"); // split interaction identifier by tab
			String prot1 = interaction_identifier[0];
			String prot2 = interaction_identifier[1];

			/* Both interacting proteins must be present in less than 1000 interactions */
			if(proteinToNumber.get(prot1) < numberOfExcessiveInteractions && proteinToNumber.get(prot2) < numberOfExcessiveInteractions) {
				interactionToNumber2.put(inter, interactionToNumber.get(inter));
			}
		}
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
