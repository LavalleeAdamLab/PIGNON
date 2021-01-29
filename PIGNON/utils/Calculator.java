package utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import graph.Annotation;
import graph.Interaction;
import graph.Protein;
import graph.Association;

public class Calculator {

    public static boolean[] determineConnectedProteins(double[][] distance_matrix) {

        // Initialize counter
        int[] count_MaxValue = new int[distance_matrix.length];
        for (int i = 0; i < distance_matrix.length; i++) {

            // run through row i and count the number of Max Values found
            for (int j = 0; j < distance_matrix.length; j++) {
                if (distance_matrix[i][j] == Double.MAX_VALUE) {
                    count_MaxValue[i]++;
                }
            }
        }

        boolean[] proteinsToKeep = new boolean[count_MaxValue.length];

        for (int i = 0; i < count_MaxValue.length; i++) {

            if (count_MaxValue[i] < (count_MaxValue.length / 2)) {
                proteinsToKeep[i] = true;
            } else {
                proteinsToKeep[i] = false;
            }
        }
        return proteinsToKeep;
    }

    public static double computeAvg(String[] dataLine, ArrayList<Integer> indices) {
        /*
         * Calculate sum of elements contained in string [] at given indexes contained
         * in an arrayList
         */
        // array to store actual values (not null)
        ArrayList<Double> values = new ArrayList<>();
        for (int idx: indices){
            // if value of the cell is not null, we add it to the values table
            if (idx<dataLine.length && !dataLine[idx].equals("") && !dataLine[idx].equals("NA")) {
                values.add(Double.parseDouble(dataLine[idx]));
            }
        }
        // if we find more than 75% of possible values, we can compute the average.
        if ( values.size() < indices.size() * 0.75){
            throw new IllegalStateException();
        }
        double sum = 0.0;
        for (double value: values) {
            sum += value;
        }
        double average = sum / values.size();

        return average;
    }

    public static int[] findInteractingProteinsIdxsInProteinNetworkList(Interaction interactionOfInterest, ArrayList<Protein> proteinsInNetworkList) {
        /*******************************************************************************
         * Search proteinInNetworkList to identify the indexes corresponding to the
         * protein interactors in the given interaction.
         *******************************************************************************/

        int[] proteinIdxs = new int[2];

        String interactor1 = interactionOfInterest.getProtein1();
        String interactor2 = interactionOfInterest.getProtein2();

        /* Find the fold change associated to each protein in the interaction */
        boolean prot1found = false;
        boolean prot2found = false;
        int protCount = 0;

        /* Search through the list of network proteins, for the two proteins in this interaction
         * proteinCount < size of protein list exist since disconnect proteins have been removed from the
         * network  */

        while ((prot1found == false || prot2found == false) && protCount < proteinsInNetworkList.size()) {

            String currentProtein = proteinsInNetworkList.get(protCount).getProteinName();

            if (currentProtein.equals(interactor1)) {
                prot1found = true;
                proteinIdxs[0] = protCount;
            } else if (currentProtein.equals(interactor2)) {
                prot2found = true;
                proteinIdxs[1] = protCount;
            }

            protCount++;
        }

        return proteinIdxs;
    }

    public static int countNumberOfInteractionsPerProtein(Protein proteinOfInterest, ArrayList<Interaction> networkPPIList) {
        /***************************************************************
         * Method computes the number of interactions for a protein (prot) in the network
         * Input : Protein and PPIs
         ***************************************************************/
        int interactionCounter = 0;

        for (int i = 0; i < networkPPIList.size(); i++) {

            int proteinInteractor1 = networkPPIList.get(i).getID1();
            int proteinInteractor2 = networkPPIList.get(i).getID2();

            if (proteinOfInterest.getProteinId() == proteinInteractor1 || proteinOfInterest.getProteinId() == proteinInteractor2) {
                interactionCounter++;
            }
        }

        return interactionCounter;
    }
    
    /**
     * Computes the total pairwise distance from a list of protein indexes of
     * interest corresponding to proteins in the network, using the distances found
     * in the distance matrix
     * 
     * @param distance_matrix
     * @param proteinIdxList : list of protein indexes in distance matrix
     * @return total pairwise distance
     **/
    public static double computeTPD(double[][] distance_matrix, ArrayList<Integer> proteinIdxList) {
        

        double distance = 0; // initialize distance

        // read distance matrix to get sum of distances between significant proteins
        for (int i = 0; i < proteinIdxList.size(); i++) {
            for (int j = i + 1; j < proteinIdxList.size(); j++) {
            	
                if (i > distance_matrix.length || j > distance_matrix.length) {
                    System.out.println("error");
                }
                
                if (distance_matrix[proteinIdxList.get(i)][proteinIdxList.get(j)] == Double.MAX_VALUE || distance == Double.MAX_VALUE) {
                    distance = Double.MAX_VALUE;
                } else {
                    // indexes of proteins of interests are found in the array list indexProt
                    distance += distance_matrix[proteinIdxList.get(i)][proteinIdxList.get(j)];
                }
                
            }
        }
        if(distance == 0.0) {
        	System.out.println("tpd = 0");
        }
        return (double) Math.round(distance * 100d) / 100d;

    }

    /**
     * Compute the p-value of a cluster from the previously generated probOfTPDmap
     */
    public static double compute_pValue(HashMap<Double, Double> probOfTPDmap, double TPD) {
        
        double p_value1 = 0; // initialize p-value

        for (double distance : probOfTPDmap.keySet()) {
            if (distance <= TPD) { // frequencies of distances smaller then the obtained TPD contribute to the p-value
                p_value1 += probOfTPDmap.get(distance);
            }

        }
        return p_value1;
    } // end compute p_value

    public static int selectRandomIndexFromList(int listSize) {
        /********************************************
         * Select a random index, from an arrayList
         ********************************************/

        Random r = new Random(); // generate random element
        int index = r.nextInt(listSize); // select a random index of an association

        return index;
    }


    public static void shuffleGoProtAssociations(ArrayList<Annotation> annotationGoList) {
        /******************************************************************************************************
         * Select 2 proteins to swap in list of goTerms.
         * Check that the chosen goTerms and proteins aren't the same.
         * Check that proteins to swap don't already exist in the next goTerm.
         * Swap proteins: previous association will no longer exist in the list and the new association will be added.
         *
         * Note : association = {goTerm		protein_idxOfProt}
         *******************************************************************************************************/

        /* Initialize list of goterm-protein association list, Association {GoTerm	ProteinName	ProteinIdx}
         * Initialize shuffled go-term list, to contain only the goTerm name */
        ArrayList<Association> goTerm_ProtAssociationList = new ArrayList<Association>();
        ArrayList<Annotation> shuffledAnnotationGoList = new ArrayList<Annotation>();
        HashMap<String, int[]> goTerm_associationIntervalMap = new HashMap<String, int[]>();

        int indexCount = 0;

        for (int i = 0; i < annotationGoList.size(); i++) {

            Annotation annotation1 = annotationGoList.get(i);

            Annotation shuffleAnnotation = new Annotation(annotation1.getName(), annotation1.getNumberOfProteins()); // new Cluster containing only goTerm
            shuffledAnnotationGoList.add(shuffleAnnotation);

            int[] associationInterval = new int[2];
            associationInterval[0] = indexCount; // initial association bound
            associationInterval[1] = indexCount + annotation1.getNumberOfProteins(); // final association bound

            indexCount += annotation1.getNumberOfProteins(); // set indexCount for next association

            goTerm_associationIntervalMap.put(annotation1.getName(), associationInterval);

            for (int j = 0; j < annotation1.getProteinIds().size(); j++) {

                Association association1 = new Association(annotation1.getName(), annotation1.getProteinIds().get(j), annotation1.getIdxProteinsList().get(j));
                goTerm_ProtAssociationList.add(association1);
            }

        }

        /*****************
         * Time concern!
         *****************/

        int numSwaps = 1000 * goTerm_ProtAssociationList.size();

        for (int i = 0; i < numSwaps; i++) {

			/*if (i%(goTerm_ProtAssociationList.size()*10) == 0) {
				System.out.println("Swap = " + i/(goTerm_ProtAssociationList.size()*10) +"%");
			}*/

            int indexAssociation1, indexAssociation2;
            Association association1, association2;
            String go1, go2;
            String prot1 = "";
            String prot2 = "";
            boolean checkGo = true;
            boolean checkProt = true;
            boolean checkProt1inAssociation2 = true;
            boolean checkProt2inAssociation1 = true;

            do {
                /* Select two pseudo-random association */
                indexAssociation1 = selectRandomIndexFromList(goTerm_ProtAssociationList.size());
                indexAssociation2 = selectRandomIndexFromList(goTerm_ProtAssociationList.size());

                association1 = goTerm_ProtAssociationList.get(indexAssociation1);
                association2 = goTerm_ProtAssociationList.get(indexAssociation2);

                go1 = association1.getGoTerm();
                go2 = association2.getGoTerm();

                checkGo = go1.equals(go2); // check selected go terms aren't the same

                if (!checkGo) {

                    prot1 = association1.getGeneID();
                    prot2 = association2.getGeneID();

                    checkProt = prot1.equals(prot2); // check selected proteins aren't the same

                    if (!checkProt) {

                        /* Check go Id's don't contain the others selected protein : obtain a list of all proteins associated to selected go term */
                        ArrayList<String> currentProteinsAssociatedGo1 = new ArrayList<String>();
                        ArrayList<String> currentProteinsAssociatedGo2 = new ArrayList<String>();

                        int[] intervalAssociation1 = goTerm_associationIntervalMap.get(go1);
                        int[] intervalAssociation2 = goTerm_associationIntervalMap.get(go2);

                        for (int a1 = intervalAssociation1[0]; a1 < intervalAssociation1[1]; a1++) {
                            currentProteinsAssociatedGo1.add(goTerm_ProtAssociationList.get(a1).getGeneID());
                        }

                        for (int a2 = intervalAssociation2[0]; a2 < intervalAssociation2[1]; a2++) {
                            currentProteinsAssociatedGo2.add(goTerm_ProtAssociationList.get(a2).getGeneID());
                        }

                        checkProt2inAssociation1 = currentProteinsAssociatedGo1.contains(prot2);
                        checkProt1inAssociation2 = currentProteinsAssociatedGo2.contains(prot1);
                    }
                }

            } while (checkGo || checkProt || checkProt1inAssociation2 || checkProt2inAssociation1); // if any of the checks are true, re-do selection

            /* Define new associations */
            association1.setGeneID(prot2);
            association1.setProteinIdx(indexAssociation2);

            association2.setGeneID(prot1);
            association2.setProteinIdx(indexAssociation1);

        }

        /* Re-structure association list into a cluster list. Searches existing shuffled goTerm list to
         * add individual element of new associations */

        System.out.println("GoTerm" + "\t" + "ProteinNames" + "\t" + "ProteinIdxes");

        for (int j = 0; j < shuffledAnnotationGoList.size(); j++) {

            Annotation shuffledGoAnnotation = shuffledAnnotationGoList.get(j);

            String goTerm = shuffledGoAnnotation.getName();
            System.out.print(goTerm + "\t");

            int[] intervalAssociation1 = goTerm_associationIntervalMap.get(goTerm);

            for (int i = intervalAssociation1[0]; i < intervalAssociation1[1]; i++) {

                Association currentAssociation = goTerm_ProtAssociationList.get(i);
                shuffledGoAnnotation.addGeneSymbol(currentAssociation.getGeneID());
                shuffledGoAnnotation.addProteinIdx(currentAssociation.getProteinIdx());
            }

            ArrayList<String> geneSymbolsList = shuffledGoAnnotation.getProteinIds();
            ArrayList<Integer> proteinIdxList = shuffledGoAnnotation.getIdxProteinsList();

            for (int h = 0; h < geneSymbolsList.size(); h++) {
                System.out.print(geneSymbolsList.get(h) + "|");
            }

            System.out.print("\t");

            for (int k = 0; k < proteinIdxList.size(); k++) {
                System.out.print(proteinIdxList.get(k) + "|");
            }

            System.out.print("\n");
        }



		/*for(int i=0; i<goTerm_ProtAssociationList.size(); i++) {

			Association currentAssociation = goTerm_ProtAssociationList.get(i);
			String currentGoTerm = currentAssociation.getGoTerm();
			String currentGoProt = currentAssociation.getGeneID();
			int currentProtIdx = currentAssociation.getProteinIdx();

			boolean goTermFound = false;
			int goTermCount = 0;

				while(!goTermFound) {

					Cluster shuffledGoAnnotation = shuffledAnnotationGoList.get(goTermCount);
					String shuffledGoTerm = shuffledGoAnnotation.getName();

					if(shuffledGoTerm == currentGoTerm) {

						shuffledGoAnnotation.addGeneSymbol(currentGoProt);
						shuffledGoAnnotation.addProteinIdx(currentProtIdx);

						goTermFound = true;
					}
					goTermCount++;
				}
		}  */


//		return shuffledAnnotationGoList;
    }

    public static ArrayList<Annotation> shuffleGoProtAssociations2(ArrayList<Annotation> annotationGoList) {
        /******************************************************************************************************
         * Select 2 proteins to swap in list of goTerms.
         * Check that the chosen goTerms and proteins aren't the same.
         * Check that proteins to swap don't already exist in the next goTerm.
         * Swap proteins: previous association will no longer exist in the list and the new association will be added.
         *
         * Note : association = {goTerm		protein_idxOfProt}
         *******************************************************************************************************/

        /* Initialize list of goterm-protein association list, Association {GoTerm	ProteinName	ProteinIdx}
         * Initialize shuffled go-term list, to contain only the goTerm name */
        ArrayList<Association> goTerm_ProtAssociationList = new ArrayList<Association>();
        ArrayList<Annotation> shuffledAnnotationGoList = new ArrayList<Annotation>();
        HashMap<String, HashSet<String>> goTerm_associationMap = new HashMap<String, HashSet<String>>();

        for (int i = 0; i < annotationGoList.size(); i++) {

            Annotation annotation1 = annotationGoList.get(i);

            Annotation shuffleAnnotation = new Annotation(annotation1.getName(), annotation1.getNumberOfProteins()); // new Cluster containing only goTerm
            shuffledAnnotationGoList.add(shuffleAnnotation);

            HashSet<String> goProteinsSet = new HashSet<String>();

            for (int j = 0; j < annotation1.getProteinIds().size(); j++) {

                // Association(String _goTerm, String _geneID, int _networkProteinIdx)
                Association association1 = new Association(annotation1.getName(), annotation1.getProteinIds().get(j), annotation1.getIdxProteinsList().get(j));
                goTerm_ProtAssociationList.add(association1);

                goProteinsSet.add(annotation1.getProteinIds().get(j));
            }

            goTerm_associationMap.put(annotation1.getName(), goProteinsSet);
        }

        /*****************
         * Time concern!
         *****************/

        int numSwaps = 1000 * goTerm_ProtAssociationList.size();

        for (int i = 0; i < numSwaps; i++) {


			/*if (i%(goTerm_ProtAssociationList.size()*10) == 0) {
				System.out.println("Swap = " + i/(goTerm_ProtAssociationList.size()*10) +"%");
			}  */

            int indexAssociation1, indexAssociation2;
            Association association1, association2;
            String go1, go2;
            String prot1 = "";
            String prot2 = "";
            int prot1idx = 0;
            int prot2idx = 0;
            boolean checkGo = true;
            boolean checkProt = true;
            boolean checkProt1inAssociation2 = true;
            boolean checkProt2inAssociation1 = true;

            HashSet<String> proteinsInAssociation1;
            HashSet<String> proteinsInAssociation2;

            do {
                /* Select two pseudo-random association */
                indexAssociation1 = selectRandomIndexFromList(goTerm_ProtAssociationList.size());
                indexAssociation2 = selectRandomIndexFromList(goTerm_ProtAssociationList.size());

                association1 = goTerm_ProtAssociationList.get(indexAssociation1);
                association2 = goTerm_ProtAssociationList.get(indexAssociation2);

                go1 = association1.getGoTerm();
                go2 = association2.getGoTerm();

                checkGo = go1.equals(go2); // check selected go terms aren't the same

                if (!checkGo) {

                    prot1 = association1.getGeneID();
                    prot2 = association2.getGeneID();

                    prot1idx = association1.getProteinIdx();
                    prot2idx = association2.getProteinIdx();

                    checkProt = prot1.equals(prot2); // check selected proteins aren't the same

                    if (!checkProt) {
                        checkProt2inAssociation1 = goTerm_associationMap.get(go1).contains(prot2);
                        checkProt1inAssociation2 = goTerm_associationMap.get(go2).contains(prot1);
                    }
                }

            } while (checkGo || checkProt || checkProt1inAssociation2 || checkProt2inAssociation1); // if any of the checks are true, re-do selection

            /* Define new associations */
            association1.setGeneID(prot2);
            association1.setProteinIdx(prot2idx);

            association2.setGeneID(prot1);
            association2.setProteinIdx(prot1idx);

            /* Update protein Set */
            proteinsInAssociation1 = goTerm_associationMap.get(go1);
            proteinsInAssociation2 = goTerm_associationMap.get(go2);

            proteinsInAssociation1.remove(prot1);
            proteinsInAssociation1.add(prot2);

            proteinsInAssociation2.remove(prot2);
            proteinsInAssociation2.add(prot1);

        }

		/* Re-structure association list into a cluster list. Searches existing shuffled goTerm list to
		 * add individual element of new associations

		System.out.println("GoTerm" +"\t"+ "ProteinNames" + "\t" + "ProteinIdxes");

		for(int j=0; j<shuffledAnnotationGoList.size(); j++) {

			Cluster shuffledGoAnnotation = shuffledAnnotationGoList.get(j);

			String goTerm = shuffledGoAnnotation.getName();
			System.out.print(goTerm + "\t");

			int[] intervalAssociation1 = goTerm_associationMap.get(goTerm);

			for(int i=intervalAssociation1[0]; i<intervalAssociation1[1]; i++) {

				Association currentAssociation = goTerm_ProtAssociationList.get(i);
				shuffledGoAnnotation.addGeneSymbol(currentAssociation.getGeneID());
				shuffledGoAnnotation.addProteinIdx(currentAssociation.getProteinIdx());
			}

			ArrayList<String> geneSymbolsList = shuffledGoAnnotation.getGeneSymbols();
			ArrayList<Integer> proteinIdxList = shuffledGoAnnotation.getIdxProteinsList();

			for(int h=0; h<geneSymbolsList.size(); h++) {
				System.out.print(geneSymbolsList.get(h) + "|");
			}

			System.out.print("\t");

			for(int k=0; k<proteinIdxList.size(); k++) {
				System.out.print(proteinIdxList.get(k) + "|");
			}

			System.out.print("\n");
		} */


        int associationCount = 0;

        for (int i = 0; i < shuffledAnnotationGoList.size(); i++) {

            Annotation shuffledGoAnnotation = shuffledAnnotationGoList.get(i);
            String shuffledGoTerm = shuffledGoAnnotation.getName();

            boolean newGoTerm = false;
            int protCount = 0;

            do {
                Association currentAssociation = goTerm_ProtAssociationList.get(associationCount);
                String currentGoTerm = currentAssociation.getGoTerm();

                if (shuffledGoTerm.equals(currentGoTerm)) {

                    shuffledGoAnnotation.addGeneSymbol(currentAssociation.getGeneID());
                    shuffledGoAnnotation.addProteinIdx(currentAssociation.getProteinIdx());

                    associationCount++;
                    protCount++;
                } else {
                    newGoTerm = true;
                }

            } while (!newGoTerm && associationCount < goTerm_ProtAssociationList.size());

            if (protCount != shuffledAnnotationGoList.get(i).getNumberOfProteins()) {
                System.out.println("error in protCount");
            }
        }

        return shuffledAnnotationGoList;
    }

    public static double computeFDR(ArrayList<Annotation> goAnnotationList, ArrayList<Annotation> goAnnotationShuffleList, double pThreshold) {
        /***************************************************************************************************
         * Compute the false discovery rate for goTerms surpassing a certain p-value threshold (pThreshold)
         *
         * FDR = nOfShuffledGo_pval that passed / nGo_pval that passed
         ***************************************************************************************************/

        double fdr;

        int goCount = 0;
        int shuffledCount = 0;
        
        for (int i = 0; i < goAnnotationList.size(); i++) {

            /* Count goTerms that have a p-value smaller than threshold */
            Annotation annotation = goAnnotationList.get(i);
            double pval = annotation.getPvalue();

            if (pval <= pThreshold) {
                goCount++;
            }


            /* Count shuffled goTerms that have a p-value smaller than the threshold*/
            Annotation shuffledAnnotation = goAnnotationShuffleList.get(i);
            double shuffled_pval = shuffledAnnotation.getPvalue();

            if (shuffled_pval <= pThreshold) {
                shuffledCount++;
            }

        }

        fdr = shuffledCount / (double) goCount;

        return fdr;
    }

    public static HashSet<Integer> determineUniProtIDsInNetwork(ArrayList<Protein> proteinsInNetworkList) {

        HashSet<Integer> uniProtIDsInNetwork = new HashSet<Integer>();

        for (int i = 0; i < proteinsInNetworkList.size(); i++) {
            uniProtIDsInNetwork.add(proteinsInNetworkList.get(i).getProteinId());
        }

        return uniProtIDsInNetwork;
    }

    public static int computeGoPassFDRthreshold(double p_val, ArrayList<Annotation> goAnnotationList) {

        int goCount = 0;

        for (int i = 0; i < goAnnotationList.size(); i++) {
            Annotation goAnnotation = goAnnotationList.get(i);
            //System.out.println(goAnnotation.getPvalue());
            if (goAnnotation.getPvalue() <= p_val) {
                goCount++;
            }
        }


        return goCount;
    }

}