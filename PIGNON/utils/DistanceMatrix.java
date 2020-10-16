package utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import graph.Interaction;
import graph.Protein;

public class DistanceMatrix {

	public static void computeDistanceMatrix(ArrayList<Interaction> networkInteractionsList, ArrayList<Protein> networkProteinsList, String outputFile) {
		/* Compute distance matrix using the list of proteins in the network and the list of known interactions in the network.
		 * This distance matrix assumes a distance of 1 if proteins interact. Outputs file : ./IO_files/DistanceMatrix.txt */
		/* Generate a weighted distance matrix based on the fold change of proteins in the network  */
		
		double[][] distance_matrix = new double[networkProteinsList.size()][networkProteinsList.size()];
		
		for (int i=0; i<networkProteinsList.size(); i++) {
			
			HashMap<String, Double> interactingProteins = new HashMap<String, Double>(); // initialize Map to contain names of interacting proteins and their weight
			Protein prot1 = networkProteinsList.get(i); // get protein object
			
			
			for (int h = 0; h < networkInteractionsList.size(); h++) {
			/* find all proteins that protein(i) interacts with by going through all 
			 * possible interactions stored in the interaction list and store in map */
			 
				Interaction inter = networkInteractionsList.get(h); // get interaction
				
				if (inter.getProtein1().equals(prot1.getProteinName())) {
					interactingProteins.put(inter.getProtein2(), inter.getWeight());
				} else if (inter.getProtein2().equals(prot1.getProteinName())) {
					interactingProteins.put(inter.getProtein1(), inter.getWeight());
				}
			}

			for (int j = 0; j < distance_matrix.length; j++) {
			/* Initialize distance matrix */
				if (i == j) { // if it's the same proteins
					distance_matrix[i][j] = 0;
					continue;
				}

				if (interactingProteins.containsKey(networkProteinsList.get(j).getProteinName())) { // if protein i and j are connected
					distance_matrix[i][j] = interactingProteins.get(networkProteinsList.get(j).getProteinName()); // set strength
				} else { // otherwise set to max_value
					distance_matrix[i][j] = Double.MAX_VALUE;
				}
			}
		}
		
		/* Implement the floyd-warshall algorithm */
		for (int k = 0; k < distance_matrix.length; k++) {
			
			if(k%100 == 0) {
				System.out.println("k = " + k);
			}
			
			for (int i = 0; i < distance_matrix.length; i++) {
				for (int j = 0; j < distance_matrix.length; j++) {
					if (distance_matrix[i][j] > (distance_matrix[i][k] + distance_matrix[k][j])) {
						distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j];
					}
				}
			}
		}

		//Print distance matrix
		try {
			BufferedWriter o = new BufferedWriter(new FileWriter(outputFile));

			for (int i = 0; i < distance_matrix.length; i++) {
				for (int j = 0; j < distance_matrix.length; j++) {
					o.write(distance_matrix[i][j] + "\t");
				}
				o.write("\n");
				o.flush();
			}
			o.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void updateDistanceMatrix(boolean[] proteinsToKeep, double[][] distance_matrix, String filePath) {

		int protCount = 0;

		// Determine number of proteins that are connected to the network
		for(int i=0; i<proteinsToKeep.length; i++) {
			if(proteinsToKeep[i]) {
				protCount++;
			} 
		}

		System.out.println("Number proteins to keep : " + protCount);

		double[][] updatedDistanceMatrix = new double[protCount][protCount];

		int row = 0;
		for(int i=0; i<distance_matrix.length; i++) {
			int col = 0;

			if(proteinsToKeep[i]) {
				for (int j=0; j<distance_matrix.length; j++){
					if(proteinsToKeep[i] == true && proteinsToKeep[j] == true) {
						//System.out.println("row=" +row +"; col=" +col +"; i=" +i +"; j=" +j);
						updatedDistanceMatrix[row][col] = distance_matrix[i][j];
						col++;
					}
				}
				row++;
			}

		}


		try {
			// Define new text file

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(filePath)));

			for (int i = 0; i < updatedDistanceMatrix.length; i++) {
				for (int j = 0; j < updatedDistanceMatrix.length; j++) {

					out.write(updatedDistanceMatrix[i][j] + "\t");
				}
				out.write("\n");
				out.flush();
			}

			out.close();

		} catch (Exception e) {
			e.printStackTrace();
		}	
	}
}
