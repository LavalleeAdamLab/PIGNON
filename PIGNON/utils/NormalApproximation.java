package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.distribution.NormalDistribution;

import graph.Annotation;

public class NormalApproximation {

	public static NormalDistribution computeProbabilityDensityFunction(HashMap<Double, Double> distributionMap) {


		double mean = computeDistributionMean(distributionMap); 
		double sd = computeDistributionStandardDeviation(distributionMap, mean); 

		NormalDistribution distribution = new NormalDistribution(mean, sd);
		//printNormalDistribution(distribution, nProts, distributionMap);

		return distribution;
	}

	public static double computeDistributionMean(HashMap<Double, Double> distributionMap) { 

		double mean = 0;

		for (double distance : distributionMap.keySet()) {		
			mean += (distance*distributionMap.get(distance));
		}
		//mean = mean/(double) nSampling;

		return mean;
	}

	public static double computeDistributionStandardDeviation(HashMap<Double, Double> distributionMap, double mean) {

		double variance = 0;

		for (double distance : distributionMap.keySet()) {		
			variance += (Math.pow((distance-mean), 2) * distributionMap.get(distance)) ;
		}

		double sd = Math.sqrt(variance); 
		return sd;
	}

	public static void printNormalDistribution(NormalDistribution nd, int nProteins, HashMap<Double, Double> distributionMap) {

		ArrayList<Integer> nProtsToTest = new ArrayList<Integer>(
				Arrays.asList(5, 10, 50, 100, 250));

		if(nProtsToTest.contains(nProteins)) {

			String outputFileName = "/Users/Rachel/eclipse-files/network2.0/output_files/normalDistribution_n" + nProteins + ".txt";
			File outputFile = new File(outputFileName);

			if(!outputFile.exists()) {

				/* Find upperbound of Monte Carlo Distribution */
				double upperbound = 0.0;

				for(double tpd : distributionMap.keySet()) {
					if(tpd>upperbound) {
						upperbound = tpd;
					}
				}

				//				double[] sampleDist = new double[sampleSize];
				//				
				//				for(int i=0; i<sampleSize; i++) {
				//					sampleDist[i] = (double) Math.round(nd.sample() * 100d) / 100d;
				//				}
				//				
				//				HashMap<Double, Double> sampleFreqMap = computeFreqOfTPD(sampleDist, sampleSize);

				try {
					BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

					out.write("TPD" + "\t" + "Freq" + "\n");

					//					for(double tpd : sampleFreqMap.keySet()) {
					//						out.write(tpd + "\t" + sampleFreqMap.get(tpd) + "\n");
					//						out.flush();
					//					}

					for(double i=0; i<upperbound; i = i + 0.01) {
						out.write(i + "\t" + nd.probability(i - 0.005, i+0.005) + "\n");
					}

					out.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

		}
	}

	public static HashMap<Double, Double> computeFreqOfTPD(double[] sampleDist, int sampleSize){

		HashMap<Double, Double> probabilityOfTPDMap = new HashMap<Double, Double>(); // map to contain frequency of TPD 

		/*********************************************************
		 * Map the occurrence of total pairwise distance observed during the random sampling 
		 *		 Key = total pairwise distance (double)
		 *		 Value = occurrence (double) 
		 **********************************************************/
		for (int i = 0; i < sampleDist.length; i++) {
			if (probabilityOfTPDMap.containsKey(sampleDist[i])) {
				double freq = probabilityOfTPDMap.get(sampleDist[i]);
				probabilityOfTPDMap.put(sampleDist[i], freq + 1);
			} else {
				probabilityOfTPDMap.put(sampleDist[i], 1.0);
			}
		}

		/* Replace the occurrence of total pairwise distance with it's frequency */
		for (double TPD : probabilityOfTPDMap.keySet()) {

			double occurence = probabilityOfTPDMap.get(TPD);
			double frequency = occurence / (double) (sampleSize);

			probabilityOfTPDMap.replace(TPD, frequency); // replace occurrence by frequency
		}

		/* Check that the sum of events (total pairwise distance) is equal to 1. Otherwise throw error */	
		float probability = 0; // initialize probability

		for (double distance : probabilityOfTPDMap.keySet()) {
			probability += probabilityOfTPDMap.get(distance); // add frequency of TPD to probability
		}

		if (Math.rint(probability) != 1) {
			System.out.println("Cluster error, probability = " + probability);
		}
		return probabilityOfTPDMap;
	}

	/* Make function to export mean and variance in file 
	 * Some of the array will not have values, print NA? */

	public static void exportNormalDistributionParameters(double[][] ndParams, String outputFile, int lowerBound) {

		try {
			File of = new File(outputFile);

			System.out.println("printing params");
			BufferedWriter out = new BufferedWriter(new FileWriter(of));

			out.write("nProt" + "\t" + "mean" + "\t" + "stdev" + "\n");


			for(int i=0; i<ndParams.length; i++) {

				out.write((i+lowerBound) + "\t" + ndParams[i][0] + "\t" + ndParams[i][1] + "\n");
				out.flush();
			}

			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static double importNormalDistributionParameters(ArrayList<Annotation> annotationGoList, String inputFile, int numOfSampling) {

		double minimum_pValue = 1.0;
		ArrayList<Integer> idxOfAnnotationsWithZeroPvalues = new ArrayList<>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));

			String line = in.readLine();
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				int nProt = Integer.parseInt(col[0]);
				double mean = Double.parseDouble(col[1]);
				double stdev = Double.parseDouble(col[2]);
				//	System.out.println("hello1");
				if(stdev != 0) { 
					//		System.out.println("hello2");
					NormalDistribution nd = new NormalDistribution(mean, stdev);

					for(int i = 0; i<annotationGoList.size(); i++) {
						//			System.out.println("hello3");

						if(annotationGoList.get(i).getNumberOfProteins() == nProt) {

							//						System.out.println("hello4");
							double tpd = annotationGoList.get(i).getTPD();
							double p_val = nd.probability(0, tpd);

							if(p_val == 0) {
								idxOfAnnotationsWithZeroPvalues.add(i);
							}
							annotationGoList.get(i).set_pValue(p_val);
							//							System.out.println("moyenne = " + nd.getMean() + "	std = " + nd.getStandardDeviation());
							//							System.out.println("pvalue = " + p_val + "	TPD = " + tpd);

							if(p_val < minimum_pValue && p_val !=0) {
								minimum_pValue = p_val;
							}

						}
					}

				}
				line = in.readLine();

			}
			in.close();	
		} catch (IOException e) {
			e.printStackTrace();
		}
		for(int i=0; i<idxOfAnnotationsWithZeroPvalues.size(); i++) {
			annotationGoList.get(idxOfAnnotationsWithZeroPvalues.get(i)).set_pValue(minimum_pValue);
		}
		
		return minimum_pValue;
	}

	public static HashMap<Integer, NormalDistribution> importNormalDistributionParameters(String inputFile) {
		HashMap<Integer, NormalDistribution> distributions = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(inputFile)));

			String line = in.readLine();
			line = in.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				int nProt = Integer.parseInt(col[0]);
				double mean = Double.parseDouble(col[1]);
				double stdev = Double.parseDouble(col[2]);

				if (!(mean == 0 && stdev == 0)) {
					NormalDistribution nd = new NormalDistribution(mean, stdev);
					distributions.put(nProt, nd);
				}


				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e){
			e.printStackTrace();
		}
		return distributions;
	}

	public static void exportDistributionParameters(HashMap<Integer, HashMap<Double, Double>> distributionsList, String fileName){
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			writer.write("Number of proteins\tMean\tSd\n");
			for (int key:distributionsList.keySet()){
				double mean = computeDistributionMean(distributionsList.get(key));
				double sd = computeDistributionStandardDeviation(distributionsList.get(key), mean);
				writer.write(String.format("%s\t%s\t%s\n", key, mean, sd));
			}
			writer.flush();
			writer.close();
		} catch (IOException ex){
			System.out.println("Error while exporting distributions parameters: export file was not found");
		}
	}
}
