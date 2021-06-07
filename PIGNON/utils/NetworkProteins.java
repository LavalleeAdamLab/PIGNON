package utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import graph.Interaction;
import graph.Protein;

public class NetworkProteins {
	
    public static ArrayList<Protein> getProteinsInNetwork(ArrayList<Interaction> InteractionList) {
        /* Find all proteins in the network form the Interaction list and store them as protein
         * objects in an ArrayList */
        ArrayList<Protein> proteinList = new ArrayList<Protein>();

        /* Find all proteins involved in protein network and store in a HashSet.
         * HashSets can only contain unique values, we obtain a list of unique protein names */
        HashSet<String> proteinNetworkSet = new HashSet<String>();

        for (int i = 0; i < InteractionList.size(); i++) {
            Interaction inter = InteractionList.get(i);

            proteinNetworkSet.add(inter.getProtein1());
            proteinNetworkSet.add(inter.getProtein2());
        }

        /*  Store protein names and IDs from HashSet in the ArrayList */
        Iterator<String> iterator = proteinNetworkSet.iterator(); //create an iterator for hash set

        while (iterator.hasNext()) {
            String proteinName = iterator.next();

            int protId = getProteinId(proteinName, InteractionList);
            Protein protein1 = new Protein(proteinName, protId); // call Protein Class constructor, store protein name (iterator.next())

            proteinList.add(protein1); // add protein object to protein list
        }
       
		
//        try {
//			BufferedWriter out = new BufferedWriter(new FileWriter(new File ("C://Users//Rachel//Documents//PIGNON//proteinsInStringNetwork2.tsv")));
//
//			for(String id: proteinNetworkSet ){
//				out.write(id + "\n");
//				out.flush();
//			}
//
//			out.close();		
//		}catch (Exception e) {
//			e.printStackTrace();
//		} 

        return proteinList;
    } // close getProteinsInNetwork
    
    private static int getProteinId(String protName, ArrayList<Interaction> interactionList) {

        int protId = 0;
        int j = 0; // Initialize counter for interactions
        boolean foundID = false; // condition to run through interactions of the InteractionList until the ID
        // corresponding to the protein name is found

        while (foundID == false) {
            Interaction inter = interactionList.get(j); // get first interaction object of the list

            // compare proteins of the Interaction to the query protein. If they match,
            // reset counter & make while loop condition true. Otherwise, go to the next
            // interaction.
            if (inter.getProtein1() == protName) {
                protId = inter.getID1();
                j = 0;
                foundID = true;
            } else if (inter.getProtein2() == protName) {
                protId = inter.getID2();
                j = 0;
                foundID = true;
            } else {
                j++;
            }
        }
        return protId;
    }
    
	/***
	 * Modify the initial list of proteins in the network based on the number the MaxValue counts in the distance matrix. 
	 * Proteins that are disconnected from the network will be removed; =proteins are disconnected when half their distance values are MAX_Values
	 *
     * @param networkProteinsList
     * @param proteinsToKeep
     * @return
     */
	public static ArrayList<Protein> modifyNetworkProteinsList(ArrayList<Protein> networkProteinsList, boolean[] proteinsToKeep) {
		

		ArrayList<Protein> networkProteinsListUpdate = new ArrayList<Protein>();

		// Run through the rows of the matrix
		for (int i = 0; i < proteinsToKeep.length; i++) {

			// if the row has less than 1/2 it's values at max value it will remain in the
			// distance matrix
			if (proteinsToKeep[i]) {
				// create new array list that will be stored in the main array
				networkProteinsListUpdate.add(networkProteinsList.get(i));
			}
		}
		return networkProteinsListUpdate;
	}
	
	/**
	 * Generate easy look up for the index of protein in the existing network. 
	 * 
	 * @param proteinsInNetworkList 
	 * @return proteinsIdxInNetworkMap
	 */
	public static HashMap<String, Integer> getProteinIdexInNetwork(ArrayList<Protein> proteinsInNetworkList){
		
		HashMap<String, Integer> proteinsIdxInNetworkMap = new HashMap<String, Integer>();
		
		for(int i=0; i<proteinsInNetworkList.size(); i++) {
			/* Map protein idx (value) to protein name (key) */
			proteinsIdxInNetworkMap.put(proteinsInNetworkList.get(i).getProteinName(), i); 
		}
		
		return proteinsIdxInNetworkMap;
	}

}
