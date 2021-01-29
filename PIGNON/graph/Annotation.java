package graph;

import java.util.ArrayList;

public class Annotation {

    // instance variables
    private String name;
    private ArrayList<String> protein_names; // protein ids(!!)
    private ArrayList<String> protein_symbols;
    private ArrayList<Integer> protein_idxList;
    private ArrayList<Integer> closeProteinsIdxList;
    private int numberOfproteins;
    private double total_pairwise_distance;
    private double top_percent_pairwise_distance;
    private double p_value;
    private double adj_pvalue;
    private double false_discovery_rate;
    private double initial_false_discovery_rate;
    private double clean_false_discovery_rate;

    private double initialPValue;
    private double closestPValue;


    // constructor
    public Annotation(String _name, int _protNum) {
        this.name = _name;
        this.protein_names = new ArrayList<String>();
        this.protein_idxList = new ArrayList<Integer>();
        this.closeProteinsIdxList = new ArrayList<>();
        this.numberOfproteins = _protNum;
    }

    public Annotation(String _name, ArrayList<String> _gene_ids, ArrayList<String> _gene_symbols, ArrayList<Integer> _proteinIdxList) {
        this.name = _name;
        this.protein_names = _gene_ids;
        this.protein_symbols = _gene_symbols;
        this.protein_idxList = _proteinIdxList;
        this.numberOfproteins = _gene_ids.size();
        this.closeProteinsIdxList = new ArrayList<>();
    }

    public Annotation(String _name, ArrayList<String> _gene_ids, ArrayList<Integer> _proteinIdxList) {
        this.name = _name;
        this.protein_names = _gene_ids;
        this.protein_idxList = _proteinIdxList;
        this.numberOfproteins = _gene_ids.size();
        this.closeProteinsIdxList = new ArrayList<>();
    }

    public Annotation(String _name, ArrayList<String> _gene_symbols) {
        this.name = _name;
        this.protein_symbols = _gene_symbols;
        this.numberOfproteins = _gene_symbols.size();
    }

    public String getName() {
        return this.name;
    }

    public ArrayList<String> getProteinIds() {
        return this.protein_names;
    }

    public ArrayList<Integer> getIdxProteinsList() {
        return this.protein_idxList;
    }

    public int getNumberOfProteins() {
        return this.numberOfproteins;
    }

    public double getTPD() {
        return this.total_pairwise_distance;
    }

    public double getPvalue() {
        return this.p_value;
    }

    public double getAdjPvalue() {
        return this.adj_pvalue;
    }

    public double getFDR() {
        return this.false_discovery_rate;
    }

    public void setGeneSymbols(ArrayList<String> geneList) {
        protein_names = geneList;
    }

    public void setNumberOfProteins(int nProteins) {
        numberOfproteins = nProteins;
    }

    public void setTPD(double TPD) {
        total_pairwise_distance = TPD;
    }

    public void set_pValue(double p_value1) {
        this.p_value = p_value1;
    }

    public void set_adj_pvalue(double adj_pvalue1) {
        this.adj_pvalue = adj_pvalue1;
    }

    public void setFDR(double fdr) {
        this.false_discovery_rate = fdr;
    }

    public void addGeneSymbol(String _geneSymbol) {
        this.protein_names.add(_geneSymbol);
    }

    public void addProteinIdx(int protIdx) {
        this.protein_idxList.add(protIdx);
    }

    public void setCloseProteinsIdxList(ArrayList<Integer> closeProteinsIdxList) {
        this.closeProteinsIdxList = closeProteinsIdxList;
    }

    public ArrayList<Integer> getCloseProteinsIdxList() {
        return closeProteinsIdxList;
    }

    public String getOfficialIdBySystemId(Integer id){
        int idx = protein_idxList.indexOf(id);
        return protein_symbols.get(idx);
    }

    public double getInitialPValue() {
        return initialPValue;
    }

    public void setInitialPValue(double initialPValue) {
        this.initialPValue = initialPValue;
    }

    public double getClosestPValue() {
        return closestPValue;
    }

    public void setClosestPValue(double closestPValue) {
        this.closestPValue = closestPValue;
    }

    public double getTTPD(){
        return top_percent_pairwise_distance;
    }

    public void setTTPD(double ttpd){
        top_percent_pairwise_distance = ttpd;
    }

    public ArrayList<String> getProtein_symbols() {
        return protein_symbols;
    }

    public double getClean_false_discovery_rate() {
        return clean_false_discovery_rate;
    }

    public void setClean_false_discovery_rate(double clean_false_discovery_rate) {
        this.clean_false_discovery_rate = clean_false_discovery_rate;
    }

    public double getInitial_false_discovery_rate() {
        return initial_false_discovery_rate;
    }

    public void setInitial_false_discovery_rate(double initial_false_discovery_rate) {
        this.initial_false_discovery_rate = initial_false_discovery_rate;
    }



    @Override
    public String toString() {
        return "Cluster{cluster=" + name + "; totalPairwiseDistance=" + total_pairwise_distance
                + "; numberOfProteins=" + numberOfproteins + "; p_value=" + p_value + "}";
    }
}
