package graph;

public class Association {

	// fields
	private String GoTerm;
	private String GeneID;
	private int NetworkProteinIdx;

			// constructors
			public Association(String _goTerm, String _geneID, int _networkProteinIdx) {
				GoTerm = _goTerm;
				GeneID = _geneID;
				NetworkProteinIdx = _networkProteinIdx;
			}
	
	public String getGoTerm() {
		return this.GoTerm;
	}
	
	public String getGeneID() {
		return this.GeneID;
	}
	
	public int getProteinIdx() {
		return this.NetworkProteinIdx;
	}

	public void setGeneID(String _geneID) {
		this.GeneID = _geneID;
	}
	
	public void setProteinIdx(int _proteinIdx) {
		this.NetworkProteinIdx = _proteinIdx;
	}
	
}
