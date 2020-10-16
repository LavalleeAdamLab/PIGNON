package utils;

import java.util.ArrayList;

public class General {
	
	public static ArrayList<Integer> convertStringArrayToIntArrayList(String[] stringArrayToConvert){
	
		ArrayList<Integer> convertedList = new ArrayList<Integer>();
		
		for(int i=0; i<stringArrayToConvert.length; i++) {
			convertedList.add(Integer.parseInt(stringArrayToConvert[i]));
		}
		
		return convertedList;
	}
	
}
