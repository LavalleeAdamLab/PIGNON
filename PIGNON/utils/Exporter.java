package utils;

import java.io.*;
import java.util.ArrayList;

import graph.Annotation;
import graph.FalseDiscoveryRate;
import graph.Interaction;

public class Exporter {

    /*********************************************************************************
     * Method prints out information of interest regarding goTerms
     * Current information : goTerm name, calculated p-value and false discovery rate
     *********************************************************************************/
    public static void printGoTermFDR(ArrayList<Annotation> goAnnotationList, String exportFileName) {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(exportFileName)));

            out.write("GoTerm" + "\t" + "p-value" + "\t" + "FDR" + "\n"); // header

            for (int i = 0; i < goAnnotationList.size(); i++) {

                Annotation annotation1 = goAnnotationList.get(i);
                out.write(annotation1.getName() + "\t" + annotation1.getPvalue() + "\t" + annotation1.getFDR() + "\n");
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void testGoFDR(ArrayList<FalseDiscoveryRate> falseDiscoveryRates, String fdrExportFile) {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fdrExportFile));
            out.write("FDR" + "\t" + "Pval" + "\t" + "#Go" + "\n");

            for (FalseDiscoveryRate fdr : falseDiscoveryRates) {
                out.write(fdr.getFalseDiscoveryRate() + "\t" + fdr.getPvalue() + "\t" + fdr.getPassingGoTerms() + "\n");
                out.flush();
            }

            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void printGO_results(ArrayList<Annotation> annotationGoList, ArrayList<Annotation> shuffledAnnotationList, String goExportFile) {

//		ArrayList<Double> pvalToTest = new ArrayList<Double>(5);
//		pvalToTest.addAll(Arrays.asList(0.1, 0.05, 0.01, 0.005, 0.001));

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(goExportFile));
            out.write("GO Term\tSize\tRealTPD\tRealPvalue\tShuffledTPD\tShuffledPvalue\n");

            for (int i = 0; i < annotationGoList.size(); i++) {
            	Annotation go1 = annotationGoList.get(i);

            	for(int j=0; j<shuffledAnnotationList.size(); j++) {
            		if(shuffledAnnotationList.get(j).getName().equals(go1.getName())){
            			out.write(go1.getName() + "\t" + go1.getNumberOfProteins() + "\t" + go1.getTPD() + "\t" + go1.getPvalue() + "\t" +
            					shuffledAnnotationList.get(j).getTPD() + "\t" + shuffledAnnotationList.get(j).getPvalue() + "\n");
            			out.flush();
            		}
            	}
            	
                
                
                
            }

            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public static void printShuffledNetwork(ArrayList<Annotation> shuffledAnnotationGoList, String shuffledGOExportFile) {

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(shuffledGOExportFile));
            out.write("Shuffled_Annotation" + "\t" + "ProteinNames" + "\t" + "ProteinIdxList" + "\n");

            for (int i = 0; i < shuffledAnnotationGoList.size(); i++) {

                Annotation shuffledGo = shuffledAnnotationGoList.get(i);
                out.write(shuffledGo.getName() + "\t");

                for (int prot = 0; prot < shuffledGo.getNumberOfProteins(); prot++) {
                    out.write(shuffledGo.getProteinIds().get(prot) + "|");
                }
                out.write("\t");

                for (int idx = 0; idx < shuffledGo.getNumberOfProteins(); idx++) {
                    out.write(shuffledGo.getIdxProteinsList().get(idx) + "|");
                }
                out.write("\n");
                out.flush();
            }

            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }



   

    public static void exportCleanGoTerms(ArrayList<Annotation> clusters, String outputPath) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));
            writer.write("goTerm id\tinitial proteins count\tclosest proteins count\tclosest proteins ratio\t" +
                    "initial p-value\tclosest p-value\tinitial fdr\tclosest fdr\tinitial_proteins\tclosest_proteins\n");

            for (Annotation cluster : clusters) {
                ArrayList<String> initialProteins = getOfficialProteinNames(cluster, cluster.getIdxProteinsList());
                ArrayList<String> closestProteins = getOfficialProteinNames(cluster, cluster.getCloseProteinsIdxList());
                double ratio = cluster.getCloseProteinsIdxList().size() / (double) cluster.getNumberOfProteins();
                writer.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", cluster.getName(), cluster.getNumberOfProteins(),
                        cluster.getCloseProteinsIdxList().size(), ratio, cluster.getInitialPValue(), cluster.getClosestPValue(),
                        cluster.getInitial_false_discovery_rate(), cluster.getClean_false_discovery_rate(), String.join(",", initialProteins),
                        String.join(",", closestProteins)));
            }
            writer.flush();
            writer.close();

        } catch (IOException exc) {
            exc.printStackTrace();
        }
    }

    private static ArrayList<String> getOfficialProteinNames(Annotation cluster, ArrayList<Integer> idxs) {
        ArrayList<String> names = new ArrayList<>();
        for (Integer idx : idxs) {
            names.add(cluster.getOfficialIdBySystemId(idx));
        }
        return names;
    }

    public static void exportNetworkInteractions(ArrayList<Interaction> interactions, String outputPath) {
        try {
            BufferedWriter reader = new BufferedWriter(new FileWriter(outputPath));
            for (Interaction interaction : interactions) {
                reader.write(String.format("%s\t%s\t%s\n", interaction.getProtein1(), interaction.getProtein2(), interaction.getWeight()));
            }
            reader.flush();
            reader.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }



}