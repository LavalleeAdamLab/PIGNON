package graph;

public class FalseDiscoveryRate {
    private double FalseDiscoveryRate;
    private double Pvalue;
    private int PassingGoTerms;

    public FalseDiscoveryRate(double falseDiscoveryRate, double pvalue, int passingGoTerms) {
        FalseDiscoveryRate = falseDiscoveryRate;
        Pvalue = pvalue;
        this.PassingGoTerms= passingGoTerms;
    }

    public double getFalseDiscoveryRate() {
        return FalseDiscoveryRate;
    }

    public void setFalseDiscoveryRate(double falseDiscoveryRate) {
        FalseDiscoveryRate = falseDiscoveryRate;
    }

    public double getPvalue() {
        return Pvalue;
    }

    public void setPvalue(double pvalue) {
        Pvalue = pvalue;
    }

    public int getPassingGoTerms() {
        return PassingGoTerms;
    }

    public void setPassingGoTerms(int passingGoTerms) {
        this.PassingGoTerms = passingGoTerms;
    }
}
