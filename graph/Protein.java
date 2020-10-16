package graph;
import java.util.ArrayList;

public class Protein {
		
	private final static double DEFAULT_FOLDCHANGE = 1.0;
	
		// fields
		private String proteinName;
		private int proteinId;
		private int numberOfInteractions;
		private double fold_change;

		// constructor
		public Protein(String _name, int _id) {
			proteinName = _name;
			proteinId = _id;
			fold_change = DEFAULT_FOLDCHANGE;
		}
		
		
		// methods 
		
		public String getProteinName() {
			return this.proteinName;
		}
		
		public int getProteinId() {
			return this.proteinId;
		}
		
		public int getNumberOfInteractions() {
			return this.numberOfInteractions;
		}
		
		public double getFoldChange() {
			return this.fold_change;
		}
		
		public void setName(String protein_name) {
			proteinName = protein_name;
		}
		
		public void setId(ArrayList<Interaction> InteractionList) {
			
			int j = 0; // Initialize counter for interactions
			boolean foundID = false; // condition to run through interactions of the InteractionList until the ID
										// corresponding to the protein name is found

			while (foundID == false) {
				Interaction inter = InteractionList.get(j); // get first interaction object of the list

				// compare proteins of the Interaction to the query protein. If they match,
				// reset counter & make while loop condition true. Otherwise, go to the next
				// interaction.
				if (inter.getProtein1() == this.proteinName) {
					proteinId = inter.getID1();
					j = 0;
					foundID = true;
				} else if (inter.getProtein2() == this.proteinName) {
					proteinId = inter.getID2();
					j = 0;
					foundID = true;
				} else {
					j++;
				}
			}

		}
		
		public void setNumberOfInteractions(int nInteractions) {
			numberOfInteractions = nInteractions;
		}
		
		public void setFoldChange(double proteinFC) {
			fold_change = proteinFC;
		}
	
}
