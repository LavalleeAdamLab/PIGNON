package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class ApproximateNormalDistribuiton {

	/**
	 * Load Monte Carlo Distribution for all n (proteins) sampled. Compute mean and variance of the distribution.
	 * Output computed parameters to text file.
	 * Each distribution is in it's own file. 
	 * 
	 * @param mcDistPrefix			String - file path prefix to Monte Carlo distributions
	 * @param lowerBound			int - lower bound of proteins sampled 
	 * @param upperBound			int - upper bound of proteins sampled
	 * @param numSampling			int - number of times the network was sampled
	 * @param normalDistParamsFile	String - file path to output normal distribution parameters
	 */
	public static void getNormalDistributionParams(String mcDistPrefix, int lowerBound, int upperBound, int numSampling, String normalDistParamsFile) {

		HashMap<Integer, Double[]> normalDistributionParamsMap = new HashMap<>();
		List<Integer> missingDistributions = new ArrayList<>();

		System.out.println("Loading mc distributions: ");
		for(int i=lowerBound; i<=upperBound; i++) {

			String mcDistFile = mcDistPrefix + i;
			File f = new File(mcDistFile);

			/* Calculate params if distribution file exists */
			if(f.exists()) {			

				/* Load each distribution individually && check that freq = 1 */

				HashMap<Double, Double> distributionMap = loadMonteCarloDistributions(mcDistFile);
				boolean frequenciesEqualOne = Loader.verifySumOfFrequencies(distributionMap);

				if(frequenciesEqualOne) {
					/* Calculate normal distribution parameters */
					double mean = computeDistributionMean(distributionMap);
					double stdv = computeDistributionStandardDeviation(distributionMap, mean);

					/* Store parameters */
					normalDistributionParamsMap.put(i, new Double[]{mean, stdv});
				}
			} else {
				missingDistributions.add(i);	
			}
		}
		System.out.print("Done\n");
		/* Output normal distribution parameters */
		exportDistributionParameters(normalDistributionParamsMap, lowerBound, upperBound, normalDistParamsFile);

		/* Output missing distributions to console */
		if(!missingDistributions.isEmpty()) {
			System.out.println("Missing file for distribution, where n = " + Arrays.toString(missingDistributions.toArray()));
		}
	}

	/**
	 * Load a Monte Carlo distribution from file. Return the distribution in HashMap
	 * 
	 * @param mcDistributionFile	String - file path for distributions 
	 * @param distributionMap		HashMap<Double, Double> - map of tpd = freq
	 */
	private static HashMap<Double, Double> loadMonteCarloDistributions(String mcDistributionFile) {

		/* Store distribution information for nProt in a hash map */
		HashMap<Double, Double> distributionMap = new HashMap<Double, Double>(); 

		try {

			InputStream in = new FileInputStream(new File(mcDistributionFile));	
			BufferedReader input = new BufferedReader(new InputStreamReader(in));


			String line = input.readLine(); // header
			String nprot = line.split("\\s+")[3].split("\\)")[0];
			System.out.print(nprot + ".");

			if(Integer.parseInt(nprot)%50 ==0) {
				System.out.println();
			}

			line = input.readLine();

			while(line!=null && !line.isEmpty()) {

				String[] col = line.split("\t");

				double tpd = Double.parseDouble(col[0]);
				double freq = Double.parseDouble(col[1]);

				distributionMap.put(tpd, freq);
				line = input.readLine();
			}

			input.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

		return  distributionMap;
	}

	/**
	 * Compute the TPD mean for a given Monte Carlo distribution 
	 * 
	 * @param distributionMap	HashMap<Double, Double> - Distribution map where tpd = frequency
	 * @return mean				Double - calculated mean 
	 */
	private static double computeDistributionMean(HashMap<Double, Double> distributionMap) { 

		double mean = 0;

		for (double distance : distributionMap.keySet()) {		
			mean += (distance*distributionMap.get(distance));
		}

		return mean;
	}

	/**
	 * Compute the variance for a given Monte Carlo Distribution. 
	 * 
	 * @param distributionMap	HashMap<Double, Double> - map of distribution tdp = freq
	 * @param mean				Double - calculated mean
	 * @return variance			Double - calculate variance
	 */
	private static double computeDistributionStandardDeviation(HashMap<Double, Double> distributionMap, double mean) {

		double variance = 0;

		for (double distance : distributionMap.keySet()) {		
			variance += (Math.pow((distance-mean), 2) * distributionMap.get(distance)) ;
		}

		double sd = Math.sqrt(variance); 
		return sd;
	}

	/**
	 * Export the normal distribution parameters for each n (number of protein) sampled.
	 * 
	 * @param distributionParamsMap		HashMap<Integer, Double[]> - map of normal distribution params for given n (n = [mean, variance])
	 * @param lowerBound				Integer - lower bound of proteins sampled
	 * @param upperBound				Integer - upper bound of proteins sampled
	 * @param fileName					String - file path to contain normal distribution parameters
	 */
	private static void exportDistributionParameters(HashMap<Integer, Double[]> distributionParamsMap, int lowerBound, int upperBound, String fileName){
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(fileName)));
			writer.write("nProt\tMean\tSd\n");
			for (int nProt=lowerBound; nProt<=upperBound; nProt++){
				Double[] params = new Double[] {0.0, 0.0};

				if(distributionParamsMap.containsKey(nProt)) {
					params = distributionParamsMap.get(nProt);
				} 

				writer.write(nProt + "\t" + params[0] + "\t" + params[1] + "\n");
				writer.flush();
			}

			writer.close();
		} catch (IOException ex){
			System.out.println("Error while exporting distributions parameters: export file was not found");
		}
	}

}
