package varianceGamma;
// EuropeanOptionVGIS {
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.*;
import umontreal.ssj.functions.*;

import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.util.Chrono;

public class EuropeanOptionVGIS {
	int d = 1;
	double T = 1.0;
	double mu = -0.1436; 			
	double sigma = 0.12136; 		
	double r = 0.1; 				
	double nu = 0.2;	
	double K1 = 130.0;
	double K2 = 160.0;
	double s0 = 100.0;
	double omega = Math.log (1 - mu*nu - sigma*sigma*nu / 2.0) / nu;
	int n = 16 * 1024;  // 2^14  for Monte Carlo
	
	double theta=9.75;
	double theta1=11.9005607518060;
	double theta2=23.5189440032670;
	
	double muplus=(Math.sqrt(Math.pow(mu, 2.0)+2*Math.pow(sigma, 2.0)/nu)+mu)/2;
	double mumoins=(Math.sqrt(Math.pow(mu, 2.0)+2*Math.pow(sigma, 2.0)/nu)-mu)/2;
	double nuplus=Math.pow(muplus, 2.0)*nu;
	double numoins=Math.pow(mumoins, 2.0)*nu;
	double alphaplus=1/nu;
	double alphamoins=1/nu;
	double lambdaplus=muplus/nuplus;
	double lambdamoins=mumoins/numoins;
	
	double lambdatp1=lambdaplus-theta1;
	double lambdatm1=lambdamoins+theta1;
	double lambdatp2=lambdaplus-theta2;
	double lambdatm2=lambdamoins+theta2;
	
	double muptheta1=alphaplus/lambdatp1;
	double muptheta2=alphaplus/lambdatp2;
	double nuptheta1=alphaplus/Math.pow(lambdatp1, 2);
	double nuptheta2=alphaplus/Math.pow(lambdatp2, 2);
	double mumtheta1=alphamoins/lambdatm1;
	double mumtheta2=alphamoins/(lambdamoins+theta2);
	double numtheta1=alphamoins/Math.pow(lambdatm1,2);
	double numtheta2=alphamoins/Math.pow( lambdatm2, 2);

	double L = (lambdaplus*lambdamoins)/((lambdaplus-theta1)*(lambdamoins+theta1));
	
	RandomStream noise = new MRG32k3a();
	Chrono timer = new Chrono();
	Tally statValue = new Tally("Stats on value of Asian option for MC");
	Tally statMc = new Tally ("Mc averages for Asian option");
	Tally statIs1 = new Tally ("Is averages for Asian option");
	Tally statIs2 = new Tally ("Is averages for Asian option");
	Tally statIs_test = new Tally ("Is averages for Asian option");
	
	public EuropeanOptionVGIS () {	
	}
	
	//////////////////////    Simulation Monte-Carlo ///////////////////////////
	public double simul_mc(double K) {   
		RandomVariateGen genYp = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdaplus)); 
		RandomVariateGen genYm = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdamoins));
		double S=s0*Math.exp(r+omega+genYp.nextDouble()-genYm.nextDouble());
		double y= Math.exp(-r)*Math.max(S-K, 0.0);
		return y;
	}
//////////////////////Simulation Monte-Carlo Importance sampling ///////////////////////////
	public double simul_is(double K, double theta) {
		RandomVariateGen genYp = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdaplus-theta)); 
		RandomVariateGen genYm = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdamoins+theta));
		double G =genYp.nextDouble()-genYm.nextDouble();
		double L=Math.pow((lambdaplus*lambdamoins/((lambdaplus-theta)*(lambdamoins+theta))),1/nu)*Math.exp(-theta*G);
		double S=s0*Math.exp(r+omega+G);
		double y= L*Math.exp(-r)*Math.max(S-K, 0.0);
		return y;
	}
	
	public static void main (String[] args) {
		/////////////////////          Results                ////////////////////
		/*Comparison between Monte-carlo simulation
		  and Monte carlo simulation with Importance sampling method 
		  for the folowing parameters: mean, variance and confidence interval
		  
		  Goal to quantify the benefit of using Monte carlo with is(importance sampling)
		  to reduce the variance */
		
		EuropeanOptionVGIS a=new EuropeanOptionVGIS ();
	     int n= a.n;
	     
	     Tally statMc= a.statMc;
	     for (int i=0; i < n; i++) {
	    	 statMc.add (a.simul_mc(a.K2));
	        }
	     //System.out.println (statMc.report (0.99, 1));
	     System.out.println (statMc.formatCIStudent (0.99));
	     System.out.println (" Option value with Mc = " +
	        PrintfFormat.format (1, 2, 2, statMc.average()));
	     System.out.println (" Variance with Mc = " + statMc.variance());
	     double p = statMc.average();
	     System.out.println ();
	     	     
	     Tally statIs_test= a.statIs_test;
	     for (int i=0; i < n; i++) {
	    	 statIs_test.add (a.simul_is(a.K2, a.theta2));
	        }
	     //System.out.println (statIs_test.report (0.99, 1));
	     System.out.println ();
	     System.out.println (statIs_test.formatCIStudent (0.99));
	     System.out.println (" Option value with IS = " +
	        PrintfFormat.format (1, 2, 2, statIs_test.average()));
	     System.out.println (" Variance with IS = " + statIs_test.variance());
	     p = statIs_test.average();
	     System.out.println ();
		
	}
}