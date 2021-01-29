package utils;

import graph.Annotation;
import graph.FalseDiscoveryRate;

import java.util.ArrayList;

public class FdrCalculator {

    /* original goTerms list */
    private ArrayList<Annotation> goAnnotations;
    /* shuffled goTerms list */
    private ArrayList<Annotation> shuffledGoAnnotations;

    /**
     * Constructor of the class.
     *
     * @param goAnnotations list of goTerms
     * @param shuffledGoAnnotations list of shuffled goTerms
     */
    public FdrCalculator(ArrayList<Annotation> goAnnotations, ArrayList<Annotation> shuffledGoAnnotations) {
        this.goAnnotations = goAnnotations;
        this.shuffledGoAnnotations = shuffledGoAnnotations;
    }

    /**
     * Modifies every goTerm in both shuffled and original lists with TPD.
     *
     * @param distanceMatrix distance matrix which is used to calculate TPD.
     */
    public void modifyGoAnnotationsWithTPD(double[][] distanceMatrix) {
    	System.out.println("computing annotation TPDs:");
        Modifier.setClusterTPD(goAnnotations, distanceMatrix); // Annotation Go List
        System.out.println("");
        System.out.println("computing shuffled annotation TPDs:");
        Modifier.setClusterTPD(shuffledGoAnnotations, distanceMatrix); // Shuffled Annotations
        System.out.println("");
    }

    /**
     * Modifies every goTerm in both shuffled and original lists with a p-value. P-value is calculated using normal
     * approximation (which is built according to previously computed parameters), so it is obligatory to have a file
     * with normal distribution parameters before calling this function. If you don't have this file, you can call
     * computeNormalDistributionParameters and then pass the name of the new file to this function.
     *
     * @param distributionParametersFilePath path to the file with distribution parameter
     */
    public double modifyGoAnnotationsWithPvalueFromNormalApproximation(String distributionParametersFilePath, int numOfSampling) {
    	
    	double[] minimum_pvals = new double[2];
        minimum_pvals[0] = NormalApproximation.importNormalDistributionParameters(goAnnotations, distributionParametersFilePath, numOfSampling);
        minimum_pvals[1] = NormalApproximation.importNormalDistributionParameters(shuffledGoAnnotations, distributionParametersFilePath, numOfSampling);
        
        double min_pval = Math.min(minimum_pvals[0], minimum_pvals[1]);
        
        return min_pval;
    }

    /**
     * Computes all normal distributions parameters and stores them in the file.
     *
     * @param distributionFilePath 				path to the file with distributions
     * @param nProtToSampleUpperBound 			largest amount of protein for which we will compute distribution parameters
     * @param distributionParametersFilePath 	path to the file with output (computed parameters)
     */
    public void computeNormalDistributionParameters(String distributionFilePath, int nProtToSampleLowerBound, int nProtToSampleUpperBound, String distributionParametersFilePath) {
        Loader.loadMonteCarloDistributions(distributionFilePath, nProtToSampleLowerBound, nProtToSampleUpperBound, distributionParametersFilePath);
    	//Loader.loadDistributions2(distributionFilePath, nProtToSampleLowerBound, nProtToSampleUpperBound, distributionParametersFilePath);
    }
    
    
    /**
     * Computes the false discovery rate table.
     *
     * @return array of FalseDiscoveryRates for p-values from 1e-180 to 1.
     */
    public ArrayList<FalseDiscoveryRate> computeFdr(double min_pval) {
        ArrayList<FalseDiscoveryRate> fdrs = new ArrayList<>();
        for (double pVal = min_pval; pVal < 0.1; pVal = pVal * 5) {
            fdrs.add(computeFdrForPvalue(pVal));
        }
        for (double pVal = 0.1; pVal <= 1.01; pVal += 0.0005) {
            fdrs.add(computeFdrForPvalue(pVal));
        }
        return fdrs;
    }

    /**
     * Computes the false discovery rate for the certain p-value.
     * @param pvalue - p-value, for which we need to compute the FDR
     * @return false discovery rate
     */
    public FalseDiscoveryRate computeFdrForPvalue(double pvalue){
        double fdr = Calculator.computeFDR(goAnnotations, shuffledGoAnnotations, pvalue);
        int nGoThatPassFDR = Calculator.computeGoPassFDRthreshold(pvalue, goAnnotations);
        return new FalseDiscoveryRate(fdr, pvalue, nGoThatPassFDR);
    }

    /**
     * Modifies the FDR list so that FDR increases monotonously.
     *
     * @param fdrs list of already computed FDR values.
     * @return list of FDR values after transformation.
     */
    public ArrayList<FalseDiscoveryRate> monotonicTransformationForFdr(ArrayList<FalseDiscoveryRate> fdrs){
        for(int i = fdrs.size()-1; i > 0; i--){
            if (fdrs.get(i-1).getFalseDiscoveryRate() > fdrs.get(i).getFalseDiscoveryRate()){
                fdrs.get(i-1).setFalseDiscoveryRate(fdrs.get(i).getFalseDiscoveryRate());
            }
        }
        return fdrs;
    }
}
