package utils;

import graph.Annotation;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Cleaner {

    /* distance matrix for tpd computing */
    private double[][] distanceMatrix;
    /* goTerm */
    private Annotation cluster;
    /* path of the distributions files */
    private String distributionFilePath;
    /* criterion for removing remote proteins */
    private double minRemovalPValue;
    private HashMap<Integer, NormalDistribution> distributions;
    private FdrCalculator fdrCalculator;

    /**
     * Creates new Cleaner object.
     *
     * @param cluster                        goTerm which is going to be checked and cleaned from remote proteins
     * @param distanceMatrix                 distance matrix for tpd computing
     * @param distributionParametersFilePath location of distribution files, needed for distributions loading
     * @param removalPValueCondition         criterion for removing remote proteins
     */
    public Cleaner(Annotation cluster, double[][] distanceMatrix, String distributionParametersFilePath, double removalPValueCondition, FdrCalculator fdrCalculator) {

        this.cluster = cluster;
        this.distanceMatrix = distanceMatrix;
        this.distributionFilePath = distributionParametersFilePath;
        this.minRemovalPValue = removalPValueCondition;
        this.distributions = NormalApproximation.importNormalDistributionParameters(distributionParametersFilePath);
        this.fdrCalculator = fdrCalculator;
    }


    /**
     * Step by step removes all remote proteins from a temporary list closeProteins. In the beginning, all proteins in
     * the goTerm considered to be close, so the closeProteins list is initialized with protein_idxList of the cluster.
     * While the size of the new list is greater than 2 (it's nothing to remove from such list) and p-value of tpd of
     * the closeProteins list is greater than removal criterion, we are selecting the most remote protein and removing
     * it from the list. After all remote proteins will be removed, we store the list of closest proteins in the cluster,
     * field closeProteinIdxList.
     */
    public void cleanGoTerm() {
        ArrayList<Integer> closeProteins = new ArrayList<>(cluster.getIdxProteinsList());
        if (cluster.getPvalue() <= 0.5) {
            ArrayList<Integer> tempProteins = null;
            System.out.println(cluster.getName());
            while (closeProteins.size() > 2) {
                try {
                    if (distributions.keySet().contains(closeProteins.size())) {
                        if (tempProteins != null && getPValue(tempProteins) < getPValue(closeProteins)) {
                            if (getPValue(tempProteins) < getPValue(closeProteins)) {
                                closeProteins = new ArrayList<>(tempProteins);
                                tempProteins = null;
                            } else
                            {
                                break;
                            }
                        }
                        if (getPValue(closeProteins) < minRemovalPValue) {
                            break;
                        }
                        Integer remote = findRemoteProtein(closeProteins);
                        closeProteins.remove(remote);
                    } else {
                        if (tempProteins == null) {
                            tempProteins = closeProteins;
                        }
                        Integer remote = findRemoteProtein(closeProteins);
                        closeProteins.remove(remote);
                    }

                } catch (Exception e) {
                    System.out.println("Problems with cluster " + cluster.getName());
                    e.printStackTrace();
                }
            }
        }
        cluster.setCloseProteinsIdxList(closeProteins);
        cluster.setInitialPValue(getPValue(cluster.getIdxProteinsList()));
        cluster.setClosestPValue(getPValue(cluster.getCloseProteinsIdxList()));
        cluster.setInitial_false_discovery_rate(fdrCalculator.computeFdrForPvalue(cluster.getInitialPValue()).getFalseDiscoveryRate());
        cluster.setClean_false_discovery_rate(fdrCalculator.computeFdrForPvalue(cluster.getClosestPValue()).getFalseDiscoveryRate());
    }

    /**
     * Gets p-value of the given tpd for certain amount of proteins.
     *
     * @param proteinIdxs proteins list, which tpd we need to compute
     * @return p-value of getting this tpd for given amount of proteins
     */
    private double getPValue(ArrayList<Integer> proteinIdxs) {
        if (distributions.keySet().contains(proteinIdxs.size())) {
            double tpd = Calculator.computeTPD(distanceMatrix, proteinIdxs);
            NormalDistribution distribution = distributions.get(proteinIdxs.size());
            return distribution.probability(0, tpd);
        } else {
            return -1;
        }
    }

    /**
     * Finds the most remote protein from the given list of proteins. Distance between a target protein and a group of
     * other proteins (called protein-group distance in other comments) is defined as: (tpd between all proteins in the
     * list) - (tpd between proteins in the list without the target protein).
     *
     * @param proteinIdxs list of proteins, for which we need to find the most remote protein.
     * @return index of the most remote protein.
     */
    private int findRemoteProtein(ArrayList<Integer> proteinIdxs) {
        // maxDistance will store the maximum found protein-group distance
        double maxDistance = Double.MIN_VALUE;
        // proteinIdx will store the index of the potential most remote protein
        int proteinIdx = -1;
        // used to remove and add target proteins
        ArrayList<Integer> tempProteinIdxs = new ArrayList<>(proteinIdxs);
        // tpd between all given proteins
        double totalDistance = Calculator.computeTPD(distanceMatrix, proteinIdxs);
        for (Integer idx : proteinIdxs) {
            tempProteinIdxs.remove(idx);
            // difference is protein-group distance
            double difference = totalDistance - Calculator.computeTPD(distanceMatrix, tempProteinIdxs);
            // if computed protein-group distance (difference) is greater than current protein-group distance than
            // overwrite maxDistance with difference and current proteinIdx with idx
            if (difference > maxDistance) {
                maxDistance = difference;
                proteinIdx = idx;
            }
            tempProteinIdxs.add(idx);
        }
        return proteinIdx;
    }
}
