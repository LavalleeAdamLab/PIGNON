package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

import graph.Protein;

public class WeightInteractions {

	/**
	 * Compute Fold change between Her2 and Triple Negative patients from the index list of each group and the file path. 
	 * File contains protein expression.
	 * 
	 * @param foldChangeDataFile file with breast cancer dataset
	 * @param proteinsInNetworkList list of proteins in network
	 * @param testHer2vsTN
	 * @param runOnComputeCanada
	 */
	public static void loadProteinExpressionDataAccrossSamples(String foldChangeDataFile, ArrayList<Protein> proteinsInNetworkList, String condition1, String condition2) {


		// InteractionToNumber will contain the number of each possible interaction
		HashMap<String, Double> foldChangeMap = new HashMap<String, Double>();

		try {

			InputStream in = new FileInputStream(new File(foldChangeDataFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();	// header which contains Sample names

			String[] col = line.split("\t");

			ArrayList<Integer> condition1idx = new ArrayList<Integer>();
			ArrayList<Integer> condition2idx = new ArrayList<Integer>();

			for(int i=1; i<col.length; i++) {
				if(col[i].matches(condition1 + "(.*)")) {
					condition1idx.add(i);
				} else if(col[i].matches(condition2 +"(.*)")) {
					condition2idx.add(i);
				} 
			}

			// read next line
			line = input.readLine();

			while (line != null) {

				// split new line in string array
				//				total_rows++;
				col = line.split("\t");
				try {
					double averageCondition1 = 1;
					double averageCondition2 = 1;

					averageCondition1 = Calculator.computeAvg(col, condition1idx);
					averageCondition2 = Calculator.computeAvg(col, condition2idx);


					// expression = fold change for Her2
					double foldChange = averageCondition1 / averageCondition2;

					// store protein name(col[0]) and fold change (expression ratio) in HashMap
					// (FoldChange)

					if(!Double.isNaN(foldChange)){
						if(col[0] != null) {
							foldChangeMap.put(col[0], foldChange);
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

		int countprot = 0; // initialize counter to check for proteins found in the network from the project data

		/* Extract info from fold change map to be stored in Protein object */
		for (int i = 0; i < proteinsInNetworkList.size(); i++) {	

			String protein_name = proteinsInNetworkList.get(i).getProteinName(); // get protein name

			// if protein is found in the fold change map, store fold change in Protein object
			if (foldChangeMap.containsKey(protein_name)) {
				double fold_change = foldChangeMap.get(protein_name);
				proteinsInNetworkList.get(i).setFoldChange(fold_change);
				countprot++;
			}
		}
		System.out.println("proteins from project found in network " +countprot +"; total proteins in project data " +foldChangeMap.size()
		+"; represents " +(countprot/(double)foldChangeMap.size())); 
	}

	public static void loadFoldChangeData(String foldChangeDataFile, ArrayList<Protein> proteinsInNetworkList) {

		HashMap<String, Integer> proteinIdxMap = getProteinIndexes(proteinsInNetworkList);
		int countprot = 0; 
		int totalprot = 0;
		try {
			InputStream in = new FileInputStream(new File(foldChangeDataFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			System.out.println("Proteins not found:");
			while(line !=null) {
				
				String geneName = line.split("\t")[0];
				double fc = Double.parseDouble(line.split("\t")[1]);
				if(!Double.isNaN(fc)){
					if(proteinIdxMap.containsKey(geneName)) {
						countprot++;
						proteinsInNetworkList.get(proteinIdxMap.get(geneName)).setFoldChange(fc);
					} else {
						System.out.print(geneName + "|");
					}
				}
				line = input.readLine();
				totalprot++;
			}
			input.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("\nproteins from project found in network " +countprot +"; total proteins in project data " + totalprot
		+"; represents " +(countprot/(double)totalprot));
	}

	private static HashMap<String, Integer> getProteinIndexes(ArrayList<Protein> proteinsInNetworkList){

		HashMap<String, Integer> proteinIdxMap = new HashMap<>();

		for(int i=0; i<proteinsInNetworkList.size(); i++) {
			proteinIdxMap.put(proteinsInNetworkList.get(i).getProteinName(), i);
		}

		return proteinIdxMap;
	}

}
