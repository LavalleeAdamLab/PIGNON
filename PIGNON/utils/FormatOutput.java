package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import graph.Protein;
import graph.Annotation;
import graph.FalseDiscoveryRate;

public class FormatOutput {

	/**
	 * Output details of significant GO annotations deemed by statistical approach:
	 * Load initially used annotations; output information only for GO annotations deemed significant
	 * Consideration: insure only protein information that was studied in our approach
	 * 
	 * @param annotationFile		String - file path for all annotations tested
	 * @param outputFile			String - file path for output
	 * @param significantGOMap		HashMap<String, Double> - maps {significant GO-term: p-value}
	 * @param proteinsInNetworkMap	HashMap<String, String> - maps {accession: protein}
	 */
	public static void printAnnotationDetails(ArrayList<Annotation> annotationList, ArrayList<Protein> networkProteinList, ArrayList<FalseDiscoveryRate> fdrs, String annotationFile, String outputFile) {

		HashMap<String, Integer> annotationIdxMap = indexAnnotationList(annotationList);
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(annotationFile)));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			String line = in.readLine(); //HEADER
			line = in.readLine();

			out.write("GO-term\tName\tCategory\tNumberOfProteinsInNetwork\tP-value\tFDR\tProteinAccessions\tProteinNames\n"); // header output file

			while(line!=null) {	
				String[] col = line.split("\t");

				if(annotationIdxMap.containsKey(col[0])) {
					/* get annotation */
					int annotationIdx = annotationIdxMap.get(col[0]);
					Annotation go1 = annotationList.get(annotationIdx); 
					
					double fdr = 1;
					String fdrToPrint = "";
					boolean zeroFDR = false;
					/* iterate through fdr - pval associations */
					for(int i=0; i<fdrs.size(); i++) {
						FalseDiscoveryRate thresholds = fdrs.get(i);
						
						/* update current fdr while threshold is greater than annotations pval */
						if(thresholds.getPvalue() >= go1.getPvalue()) {
							fdr = thresholds.getFalseDiscoveryRate();
							if (fdr == 0) {
								zeroFDR = true;
							}
							if (fdr > 0 && !zeroFDR) {
								fdrToPrint = Double.toString(fdr);
								break;	
							} else if(fdr > 0 && zeroFDR) {
								fdrToPrint = "<" + Double.toString(fdr);
								break;
							}
						}
					}

					out.write(col[0] + "\t" + col[1] + "\t" + col[2]+ "\t" + go1.getNumberOfProteins() + "\t" + 
							+ go1.getPvalue() + "\t" + fdrToPrint + "\t");

					for(int i=0; i<go1.getProteinIds().size(); i++) {
						out.write(go1.getProteinIds().get(i) + "|");
					}
					out.write("\t");

					for(int i=0; i<go1.getProtein_symbols().size(); i++) {
						out.write(go1.getProtein_symbols().get(i) + "|");
					}
					out.write("\n");
					out.flush();
				}
				line = in.readLine();
			}
			in.close();
			out.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

	public static HashMap<String, Integer> indexAnnotationList(ArrayList<Annotation> annotationList){
		HashMap<String, Integer> idxAnnotationsMap = new HashMap<>();

		for(int i=0; i<annotationList.size(); i++) {
			idxAnnotationsMap.put(annotationList.get(i).getName(), i);
		}
		return idxAnnotationsMap;
	}

}
