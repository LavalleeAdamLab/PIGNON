package utils;

import graph.Annotation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Sampling {

	/* goTerms */
	private ArrayList<Annotation> annotation_goList;
	/* matrix of TPD between proteins */
	double[][] distanceMatrix;
	private ArrayList<Integer> unweightedList;
	/* weights for every protein calculated as a number of occurrences in goTerms
	 * Key: index of the protein, value: weight */
	private HashMap<Integer, Integer> weights;
	/* represents the distribution of proteins' occurrences in goTerms.
	 * If protein was found in 3 goTherms, it will be added to this list 3 times */
	private ArrayList<Integer> weightedList;
	/* flag to set the type of random protein selection */
	private int samplingType;

	/**
	 * Construct a new Sampling object
	 *
	 * @param goAnnotations   array of goTerms
	 * @param distance_matrix distance matrix between proteins to calculate the TPD
	 * @param _weighted       flag to set the type of random protein selection
	 */
	public Sampling(ArrayList<Annotation> goAnnotations, double[][] distance_matrix, int _samplingType) {
		annotation_goList = goAnnotations;
		distanceMatrix = distance_matrix;
		samplingType = _samplingType;
		switch(samplingType) {
		case 1:
			weights = computeProteinFrequenciesListBySystemId();
			weightedList = getWeightedProteinsList();
		case 0:
			unweightedList = getAnnotatedProteinsList();
		}
	}

	/**
	 * Computes distributions for multiple number of proteins in range from go_start to go_stop. For each amount of
	 * proteins checks if exists a goTerm with this amount of associated proteins. If exists, then computes the
	 * distribution
	 *
	 * @param nProtToSampleLowerBound		beginning of the range of proteins that will be sampled
	 * @param nProtToSampleUpperBound		end of the range of proteins that will be sampled
	 * @param numOfTimesNetworkIsSampled 	number of times to calculate the distribution for each amount of proteins
	 */
	public void computeMultipleDistributions(int nProtToSampleLowerBound, int nProtToSampleUpperBound, int numOfTimesNetworkIsSampled, String mcDistFile) {
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(mcDistFile)));

			for (int n = nProtToSampleLowerBound; n <= nProtToSampleUpperBound; n++) { // range of proteins to sample

				/* Check if a GO term is associated to the number of proteins to be sampled in the network
				 * Prevents the need of computing a distribution that will not be used. */ 
				boolean numberOfAssociatedProteinsExists = false;
				int goTermIndex = 0;

				/* Iterate over annotation list until a GO term is found with the given lenght of proteins
				 * or the list of GO terms is exhausted  */
				while (numberOfAssociatedProteinsExists == false && goTermIndex < annotation_goList.size()) {
					if (annotation_goList.get(goTermIndex).getNumberOfProteins() == n) {
						numberOfAssociatedProteinsExists = true;
					} else {
						goTermIndex++;
					}
				}

				/* Given a GO Term is associated with n number of proteins, we perform Monte Carlo Sampling
				 * and output the resulting distribution to the console */


				if (numberOfAssociatedProteinsExists) {
					HashMap<Double, Double> distribution = compute_distributionTPD(n, numOfTimesNetworkIsSampled);

					out.write("TPD (n = " + n + ")" + "\t" + "Frequency" + "\n");
					for (double dist : distribution.keySet()) {
						out.write(dist + "\t" + distribution.get(dist) + "\n");
						out.flush();
					}
				} else {
					out.write("TPD (n = " + n + ") Not found\n");
					out.flush();
				}
				out.write("\n");
				out.flush();
			}
			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * Computes the distribution of TPD for certain amount of randomly selected proteins.
	 *
	 * @param distanceMatrix 			distance matrix of proteins shortest pairwise distances
	 * @param timesToSampleNetwork   	number of iterations to sample Network
	 * @return the distribution of TPD for certain amount of randomly selected proteins
	 */
	private HashMap<Double, Double> compute_distributionTPD(int numProteinsToSample, int timesToSampleNetwork) {

		HashMap<Double, Double> distribution;

		/* Sample network for x amounts of proteins, y amount of times */
		double[] tpdSampleList = sampleNetwork(numProteinsToSample, timesToSampleNetwork);

		/* Measure the frequencies of sampled total pairwise distances from x amount of proteins */
		distribution = computeFrequenciesOfSampledTPDs(tpdSampleList, timesToSampleNetwork);

		/* Sanity check to ensure the sum of TPD frequencies equals 1 */
		checkFrequencyTotal(distribution, numProteinsToSample);

		return distribution;
	}

	/**
	 * Sample network X amount of times for Y number of proteins. 
	 * 
	 * @param numProteinsToSample		Proteins to sample in network
	 * @param timesToSampleNetwork		Number of iterations to sample network
	 * @return	list of TPDs measured in network
	 */
	private double[] sampleNetwork(int numProteinsToSample, int timesToSampleNetwork) {
		double[] tpdSampleList = new double[timesToSampleNetwork]; // array to store the results of sampling {list of TPDs}

		/* Randomly select proteins from the network and compute their total pairwise distance
		 * as many times as the network is to be sampled (e.g. 10X7) */
		for (int i = 0; i < timesToSampleNetwork; i++) {
			ArrayList<Integer> randomProteins = new ArrayList<>();

			switch(samplingType) {
			case 1:
				/* select proteins from the weighted list (ie. proteins are proportional to their occurrence in annotation list */
				randomProteins = getRandomWeightedProteins(numProteinsToSample);
			case 0:
				/* select proteins from the network (ie. all proteins are represented only once) */
				randomProteins = getRandomAnnotatedProteins(numProteinsToSample);
			}

			/* compute the total pairwise distance from the proteins selected above */
			tpdSampleList[i] = Calculator.computeTPD(distanceMatrix, randomProteins);
		}

		return tpdSampleList;
	}

	/**
	 * Compute the Frequencies of TPD that were sampled from network proteins
	 * 
	 * @param sampledTPDList			List of total pairwise distances measured from sampled proteins in network
	 * @param timesNetworkWasSampled	Number of times the network was sampled (eg. 10X7)
	 * @return Monte Carlo Distribution as a map of {TPD: frequency} 
	 */
	private HashMap<Double, Double> computeFrequenciesOfSampledTPDs(double[] sampledTPDList, int timesNetworkWasSampled){
		HashMap<Double, Double> probabilityOfTPDMap = new HashMap<Double, Double>(); //map to store {TPD: frequency}

		/* calculate the occurrence of obtained total pairwise distances */
		for (int i = 0; i < sampledTPDList.length; i++) {
			if (probabilityOfTPDMap.containsKey(sampledTPDList[i])) {
				double occurrence = probabilityOfTPDMap.get(sampledTPDList[i]);
				probabilityOfTPDMap.put(sampledTPDList[i], occurrence + 1);
			} else {
				probabilityOfTPDMap.put(sampledTPDList[i], 1.0);
			}
		}

		/* calculate frequencies of obtained total pairwise distance */
		for (double TPD : probabilityOfTPDMap.keySet()) {
			double occurrence = probabilityOfTPDMap.get(TPD);
			double frequency = occurrence / (double) (timesNetworkWasSampled);
			probabilityOfTPDMap.replace(TPD, frequency); // replace occurrence by frequency
		}
		return probabilityOfTPDMap;
	}

	/**
	 * Sanity check: ensure that the sum of TPD frequencies equals 1
	 * 
	 * @param probabilityMapOfTPD
	 */
	private void checkFrequencyTotal(HashMap<Double, Double> probabilityOfTPDMap, int numProteinsToSample) {

		float probability = 0; // initialize probability

		for (double distance : probabilityOfTPDMap.keySet()) {
			probability += probabilityOfTPDMap.get(distance); // add frequency of TPD to probability
		}

		if (Math.rint(probability) != 1) { // allow slight variation in the number given the nature of small numbers
			throw new IllegalArgumentException("Distribution probability isn't 1. Number of sampled Proteins: " + numProteinsToSample + ", probability = " + probability);
		}
	}


	@SuppressWarnings("unused")
	private double calculateClosestDistancesSum(int[] distances, int l) {
		Arrays.sort(distances);
		double sum = 0;
		for (int i = 0; i < l; i++) {
			sum += distances[i];
		}
		return sum;
	}

	/**
	 * Randomly selects proteins from weighted list
	 * Update 14.01.2020: Fixed bug, when checking if protein was already selected to sample. 
	 *
	 * @param numProteinsToSample number of proteins to select from network
	 * @return array of selected proteins
	 */
	private ArrayList<Integer> getRandomWeightedProteins(int numProteinsToSample) {
		Random r = new Random(); 
		ArrayList<Integer> randomProteins = new ArrayList<>();

		/* Selection process occurs until the number of selected proteins equals number of proteins to sample from */ 
		while (randomProteins.size() < numProteinsToSample) {
			int idx = r.nextInt(weightedList.size()); // select a random protein from the weighted list in 

			/* if this protein is not already in the list, add it's index to selected */
			if (!randomProteins.contains(weightedList.get(idx))) { 
				randomProteins.add(weightedList.get(idx));
			}
		}
		return randomProteins;
	}

	/**
	 * Selects proteins randomly from network proteins
	 *
	 * @param numProteinsToSample number of proteins to select
	 * @param bound bound of random number (based on the number of proteins in the network)
	 * @return array of selected proteins
	 */
	private ArrayList<Integer> getRandomAnnotatedProteins(int numProteinsToSample) {
		Random r = new Random();
		ArrayList<Integer> randomProteins = new ArrayList<Integer>();

		/* Selection process occurs until the number of selected proteins equals number of proteins to sample from */
		while (randomProteins.size() < numProteinsToSample) {
			Integer proteinIndex = r.nextInt(unweightedList.size()); // protein as an index in the distance Matrix

			/* if this protein is not already in the list, add it's index to selected */
			if (!randomProteins.contains(unweightedList.get(proteinIndex))) {
				randomProteins.add(unweightedList.get(proteinIndex));
			}
		}
		return randomProteins;
	}

	/***
	 * Generates a list of proteins annotated by a GO term. A protein is only present ONCE. 
	 * 
	 * @return annotatedProteinList 	ArrayList<Integer> list of protein indexes in network
	 */
	private ArrayList<Integer> getAnnotatedProteinsList(){
		HashSet<Integer> annotatedProteinSet = new HashSet<>();
		/* Store index of protein in network annotated by a GO term in HashSet */
		for(Annotation go: annotation_goList) {
			for(Integer proteinIdx: go.getIdxProteinsList()) {
				annotatedProteinSet.add(proteinIdx);
			}
		}
		/* Convert HashSet to ArrayList */
		ArrayList<Integer> annotatedProteinList = new ArrayList<Integer>(annotatedProteinSet);
		return annotatedProteinList;
	}


	/**
	 * Computes the weighted list (If protein was found in 3 goTherms, it will be added to
	 * this list 3 times)
	 *
	 * @return weighted list
	 */
	private ArrayList<Integer> getWeightedProteinsList() {
		ArrayList<Integer> proteins = new ArrayList<>();
		for (int key : weights.keySet()) {
			// add every protein to the weighted list n times, where n is equal to the weight of the protein
			for (int i = 0; i < weights.get(key); i++) {
				proteins.add(key);
			}
		}
		return proteins;
	}

	/**
	 * Computes the weight for every protein based on the number of occurrences in goTerms.
	 *
	 * @return dictionary, where key is an official id of the protein in goTerms and
	 * value is the weight of this protein.
	 */
	@SuppressWarnings("unused")
	private HashMap<String, Integer> computeProteinFrequenciesListByGoTermsId() {
		HashMap<String, Integer> counts = new HashMap<>();
		for (Annotation cluster : annotation_goList) {
			for (int i = 0; i < cluster.getProteinIds().size(); i++) {
				String id = cluster.getProteinIds().get(i);
				if (!counts.containsKey(id)) {
					counts.put(id, 1);
				} else {
					int oldValue = counts.get(id);
					counts.put(id, oldValue + 1);
				}
			}
		}
		return counts;
	}

	/**
	 * Computes the weight for every protein based on the number of occurrences in goTerms.
	 *
	 * @return dictionary, where key is an id of the protein in the program and
	 * value is the weight of this protein.
	 */
	private HashMap<Integer, Integer> computeProteinFrequenciesListBySystemId() {
		HashMap<Integer, Integer> counts = new HashMap<>();
		for (Annotation cluster : annotation_goList) {
			for (int i = 0; i < cluster.getIdxProteinsList().size(); i++) {
				int id = cluster.getIdxProteinsList().get(i);
				if (!counts.containsKey(id)) {
					counts.put(id, 1);
				} else {
					int oldValue = counts.get(id);
					counts.put(id, oldValue + 1);
				}
			}
		}
		return counts;
	}

}
