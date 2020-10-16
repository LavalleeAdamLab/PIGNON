package opt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import graph.Interaction;
import graph.Protein;
import utils.Calculator;

public class NetworkInformation {

	/***
	 * Update the number of protein interactions a given protein in the network is involved in.
	 * Export the name of proteins and their number of interactions in a text file. 
	 ***/
	public static void updateNumberOfProteinInteractions(ArrayList<Protein> networkProteins, ArrayList<Interaction> networkInteractions) {

		for(int i=0; i<networkProteins.size(); i++) {

			Protein proteinOfInterest = networkProteins.get(i);
			proteinOfInterest.setNumberOfInteractions(Calculator.countNumberOfInteractionsPerProtein(proteinOfInterest, networkInteractions));
		}
		
	}
	
    /**********************************************************************************
     * Method outputs .txt file which contains the number of interaction a protein has
     * and weither or not it was quantified in the protein expression file.
     **********************************************************************************/
    public static void printNumberOfInteractionsPerProtein(ArrayList<Protein> proteinsInNetworkList, String exportFile) {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(exportFile)));

            out.write("Protein" + "\t" + "NumberOfInteractions" + "\n");
            //out.write("Protein" + "\t" + "NumberOfInteractions" + "\t" + "Quantified" + "\n");

            for (int i = 0; i < proteinsInNetworkList.size(); i++) {
                out.write(proteinsInNetworkList.get(i).getProteinName() + "\t" + proteinsInNetworkList.get(i).getNumberOfInteractions());
                //out.write(proteinsInNetworkList.get(i).getProteinName() + "\t" + proteinsInNetworkList.get(i).getNumberOfInteractions() + "\t");

				/*if(proteinsInNetworkList.get(i).getFoldChange() == 1.0) {
					out.write("No");
				} else {
					out.write("Yes");
				}*/

                out.write("\n");
                out.flush();
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    
    public static void exportNetwork(ArrayList<Protein> proteinsInNetworkList, ArrayList<Interaction> networkInteractionsList, String networkFile) {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(networkFile));

            out.write("Protein1" + "\t" + "Protein2" + "\t" + "FC1" + "\t" + "FC2" + "\n");

            for (int i = 0; i < 50000; i++) {

                Interaction networkInteraction1 = networkInteractionsList.get(i);
                int prot1ID = networkInteraction1.getID1();
                int prot2ID = networkInteraction1.getID2();

                boolean foundProt1 = false;
                boolean foundProt2 = false;

                int protCount = 0;

                double protFC1 = 1;
                double protFC2 = 1;

                while (!foundProt1 && !foundProt2 && protCount < proteinsInNetworkList.size()) {

                    Protein protX = proteinsInNetworkList.get(protCount);

                    if (protX.getProteinId() == prot1ID) {
                        protFC1 = protX.getFoldChange();
                        foundProt1 = true;
                    } else if (protX.getProteinId() == prot2ID) {
                        protFC2 = protX.getFoldChange();
                        foundProt2 = true;
                    } else {
                        protCount++;
                    }
                }

                out.write(prot1ID + "\t" + prot2ID + "\t" + protFC1 + "\t" + protFC2 + "\n");
                out.flush();
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
    
    public static void outputGraphForMCL(ArrayList<Interaction> interactionList, ArrayList<Protein> proteinList, String outputFile) {
    
    	HashSet<String> keptProteins = new HashSet<String>();
    	for(int i=0; i<proteinList.size(); i++) {
    		keptProteins.add(proteinList.get(i).getProteinName());
    	}
    	
    	try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int i=0; i<interactionList.size(); i++) {
				Interaction inter = interactionList.get(i);
				if(keptProteins.contains(inter.getProtein1()) && keptProteins.contains(inter.getProtein2())) {
					out.write(inter.getProtein1() + "\t" + inter.getProtein2() +"\n");
					out.flush();
				}
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
    	
    	
    }
    
    public static void outputProteinList(ArrayList<Protein> proteinList, String outputFile) {
    	
    	BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int i=0; i<proteinList.size(); i++) {
    			out.write(proteinList.get(i).getProteinId() + "\t" + proteinList.get(i).getProteinName() + "\n");
    			out.flush();
    		}
    		out.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
    		
    		
    }
	
}
