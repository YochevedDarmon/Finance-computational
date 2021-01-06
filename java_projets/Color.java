package devoir4;

import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import java.util.*; 

public class Color {
	int[] ordre= {1,2,3,4,5};
	int[] initial= {2,1,1,1,3};
	int[][] voisins={{0,1,0,1,1},{1,0,0,0,0},
			{0,0,0,0,1},{1,0,0,0,1},{1,0,1,1,0}};
	int color_number=4;
	int sommets_number=initial.length;
	int n=216000;
	//int n=100;
	RandomStream stream= new MRG32k3a();
	public Color() {
	}
		
	public int[] voisins_color(int[]echantillon,int i){
		int[]voisinsColor=new int[sommets_number];
			for (int j=0; j<sommets_number;j++) {
				if (voisins[i][j]>0) {voisinsColor[j]=echantillon[j];}
		}
		return voisinsColor;
	}
	
	public int[] admissible_color(int[] voisinsColor) {
	    int[] admcolor0= {1,2,3,4};
	    List<Integer> voisinsColorList = new ArrayList<Integer>(voisinsColor.length);
        for (int i : voisinsColor)voisinsColorList.add(i);
        List<Integer> admcolorList = new ArrayList<Integer>(admcolor0.length);
        for (int i : admcolor0)admcolorList.add(i);
        admcolorList.removeAll(voisinsColorList);
        int[] admcolor = new int[admcolorList.size()];
        int i= 0;
        for (int color : admcolorList) admcolor[i++] = color;
		return admcolor;
	}
	
	public int [] simulate(int[] echantillon){
		for (int j=0; j<5;j++) {
			int [] colorvoisin=voisins_color(echantillon,j);
			int[] possibility_number=admissible_color(colorvoisin);
			UniformIntGen gen=new UniformIntGen(stream,0,
					possibility_number.length-1);
			echantillon[j]=possibility_number[gen.nextInt()];
		}
		return echantillon;
	}
	
	public int[] run() {
		Map<String, Integer> possible = new HashMap<String, Integer>();
		String str=String.valueOf(initial[0]);
		for (int j=1; j<5;j++) {str+=String.valueOf(initial[j]);}
		possible.put(str, 1);
		int[] echantillon=initial;
		for (int i=1; i<n;i++) { 
			//int[]simulation=simulate(echantillon);
			echantillon=simulate(echantillon);
			str=String.valueOf(echantillon[0]);
			for (int j=1; j<5;j++) {str+=String.valueOf(echantillon[j]);}
			//possible.put(str, possible.get(str) + 1);
			if (possible.containsKey(str)) {
				possible.put(str, possible.get(str) + 1);}
			else {possible.put(str,1);}
		}
		int[]echantillon_number= new int[possible.size()];
		int i=0;
		for (String key:possible.keySet()) { 
			echantillon_number[i]=possible.get(key);
			i++;
		}
		return echantillon_number;
	}
	
	public String[] keyy() {
		Map<String, Integer> possible = new HashMap<String, Integer>();
		String str=String.valueOf(initial[0]);
		for (int j=1; j<5;j++) {str+=String.valueOf(initial[j]);}
		possible.put(str, 1);
		int[] echantillon=initial;
		for (int i=1; i<n;i++) { 
			//int[]simulation=simulate(echantillon);
			echantillon=simulate(echantillon);
			str=String.valueOf(echantillon[0]);
			for (int j=1; j<5;j++) {str+=String.valueOf(echantillon[j]);}
			//possible.put(str, possible.get(str) + 1);
			if (possible.containsKey(str)) {
				possible.put(str, possible.get(str) + 1);}
			else {possible.put(str,1);}
		}
		String[]echantillon_key= new String[possible.size()];
		int i=0;
		for (String key:possible.keySet()) { 
			echantillon_key[i]=key;
			i++;
		}
		return echantillon_key;
	}
	

	
	public static void main (String[] args) {
		System.out.println("Voici les résultats de l'exercice 4");
		Color color=new Color();
		ChiSquareDist Y= new ChiSquareDist(215);
		for (int j=0; j<100;j++) color.stream.nextDouble();
		int[]echantillon_number=color.run();
		System.out.println(echantillon_number.length);
		double x=0.0;
		for (int j=0; j<echantillon_number.length;j++) {
			System.out.println(echantillon_number[j]);
			x+=Math.pow(echantillon_number[j]-1000.0, 2)/1000.0;
		}
		System.out.println(Y.barF(x));
	}

}

