package varianceGamma;

//public class EuropeanOptionVG_ISandRQMC extends EuropeanOptionVG {
import java.io.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.DigitalNet;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloExperiment;
import umontreal.ssj.mcqmctools.RQMCExperiment;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.TruncatedDist;
import umontreal.ssj.randvar.GammaAcceptanceRejectionGen;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.stat.*;

public class EuropeanOptionVG_ISandRQMC extends EuropeanOptionVG {

	// Constructor.
	public EuropeanOptionVG_ISandRQMC(double s0, double K, double mu, double sigma, double r, double nu, double theta) {
		super(s0, K, mu, sigma, r, nu,theta); 
	}

	// Generates and returns X, with IS.
	public void simulate(RandomStream stream) {
		double gm= genYm.nextDouble();
		double b=Math.log(K/s0)-r-omega;
		double F2 = distgp.cdf(b+gm);  // U2 must be larger than this.
		double u2 = F2 + (1.0 - F2) * stream.nextDouble();
		double gpc=distgp.inverseF(u2);
		double S=s0*Math.exp(r+omega+gpc-gm);
		double Lthetap=Math.pow(lambdamoins /(lambdamoins + theta),1/nu)*Math.exp(theta*gm);
		payoff = Math.exp(-r)*Math.max(0.0, (S - K)) * (1.0 - F2)*Lthetap;
		
	}

	// Descriptor of model
	public String toString() {
		return "Simplified financial option with barriers, with IS";
	}
	
	public static void main(String[] args) throws IOException {
		/////////////////////          Results                ////////////////////
		/*Comparison between Monte-carlo simulation,
		  Monte carlo simulation with is(Importance sampling) method,
		  and Monte carlo simulation with Importance sampling method and RQMC(randomized quasi-MC) method
		  for the folowing parameters: mean, variance and confidence interval
		  
		  Goal to quantify the benefit of using Monte carlo with is and RQMC
		  by comparison using Monte carlo with is only to reduce the variance */

		EuropeanOptionVG pb = new EuropeanOptionVG( 0,  0,  0,  0,  0,  0, 0);
		double mu = -0.1436; 			
		double sigma = 0.12136; 		
		double r = 0.1; 				
		double nu = 0.2;	
		double K1 = 130.0;
		double theta = pb.theta1;
		double K2 = 160.0;
		//double theta= pb.theta2;
		double s0 = 100.0;
		double omega = Math.log (1 - mu*nu - sigma*sigma*nu / 2.0) / nu;
		int n = 16 * 1024;  // 2^14  for Monte Carlo
		int m = 32;
		
		double varMC=1;
		double varIs=1;
		double ratio=1;
		
		RandomStream stream = new LFSR113();
		Tally statX = new TallyStore("Option payoffs");  // To store the n observations of X.
		
		
		
		/////////////////   first results     /////////////////
		
		RandomStream streamRQMC = new LFSR113();
		DigitalNet p = new SobolSequence(16, 31, 2); // n = 2^{16} points in s dim.
		PointSetRandomization rand = new LMScrambleShift(streamRQMC);
		

		//theta=1.0E7;
		theta=pb.theta1;
		
		pb = new EuropeanOptionVG( s0,  K1,  mu,  sigma,  r,  nu, theta);
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
		varMC=statX.variance();
		pb = new EuropeanOptionVG_ISandRQMC(s0,  K1,  mu,  sigma,  r,  nu, theta);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varIs=statX.variance();
		ratio=varMC/varIs;		
		System.out.println ("ratio MC vs IS :"+ratio+"varmc /varis "+varMC+" "+ varIs);
		System.out.println ();
		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
			    (new EuropeanOptionVG_ISandRQMC( s0,  K1,  mu,  sigma,  r,  nu, theta), streamRQMC, p, rand, n, m));
		
		
		theta= pb.theta2;
		
		pb = new EuropeanOptionVG(s0,  K2,  mu,  sigma,  r,  nu, theta);
		System.out.println (MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX));
		varMC=statX.variance();
		pb = new EuropeanOptionVG_ISandRQMC(s0,  K2,  mu,  sigma,  r,  nu, theta);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varIs=statX.variance();
		ratio=varMC/varIs;
		System.out.println ("ratio MC vs IS :"+ratio);
		System.out.println ();
		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
			    (new EuropeanOptionVG_ISandRQMC( s0,  K2,  mu,  sigma,  r,  nu, theta), streamRQMC, p, rand, n, m));
		
		
		
		/////////////////  optimisation of theta star value /////////////////
		double[] thetas=new double [100];
		thetas[0]=1;
		pb = new EuropeanOptionVG( s0,  K1,  mu,  sigma,  r,  nu, thetas[0]);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varMC=statX.variance();
		pb = new EuropeanOptionVG_ISandRQMC(s0,  K1,  mu,  sigma,  r,  nu, thetas[0]);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varIs=statX.variance();
		ratio=varMC/varIs;		
		System.out.println ("ratio MC vs IS :"+ratio+" avec theta ="+thetas[0]);
		for (int j = 1; j < 100; j++) {
			thetas[j]=thetas[j-1]+1;
			pb = new EuropeanOptionVG( s0,  K1,  mu,  sigma,  r,  nu, thetas[j]);
			MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
			varMC=statX.variance();
			pb = new EuropeanOptionVG_ISandRQMC(s0,  K1,  mu,  sigma,  r,  nu, thetas[j]);
			MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
			varIs=statX.variance();
			ratio=varMC/varIs;		
			System.out.println ("ratio MC vs IS :"+ratio+" avec theta ="+thetas[j]);
		}
		
		System.out.println ("");
		System.out.println ("");
		thetas=new double [100];
		thetas[0]=1;
		pb = new EuropeanOptionVG( s0,  K2,  mu,  sigma,  r,  nu, thetas[0]);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varMC=statX.variance();
		pb = new EuropeanOptionVG_ISandRQMC(s0,  K2,  mu,  sigma,  r,  nu, thetas[0]);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varIs=statX.variance();
		ratio=varMC/varIs;		
		System.out.println ("ratio MC vs IS :"+ratio+" avec theta ="+thetas[0]);
		for (int j = 1; j < 100; j++) {
			thetas[j]=thetas[j-1]+1;
			pb = new EuropeanOptionVG( s0,  K2,  mu,  sigma,  r,  nu, thetas[j]);
			MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
			varMC=statX.variance();
			pb = new EuropeanOptionVG_ISandRQMC(s0,  K2,  mu,  sigma,  r,  nu, thetas[j]);
			MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
			varIs=statX.variance();
			ratio=varMC/varIs;		
			System.out.println ("ratio MC vs IS :"+ratio+" avec theta ="+thetas[j]);
		}
				
				
    	///////////////// final comparison with theta star /////////////////

		theta=20;
		pb = new EuropeanOptionVG_ISandRQMC(s0,  K1,  mu,  sigma,  r,  nu, theta);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varIs=statX.variance();
		ratio=varMC/varIs;		
		System.out.println ("ratio MC vs IS :"+ratio+"varmc /varis "+varMC+" "+ varIs);
		System.out.println ();
		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
			    (new EuropeanOptionVG_ISandRQMC( s0,  K1,  mu,  sigma,  r,  nu, theta), streamRQMC, p, rand, n, m));
		
		theta=24;
		pb = new EuropeanOptionVG_ISandRQMC(s0,  K2,  mu,  sigma,  r,  nu, theta);
		MonteCarloExperiment.simulateRunsDefaultReportStudent(pb, n, stream, statX);
		varIs=statX.variance();
		ratio=varMC/varIs;
		System.out.println ("ratio MC vs IS :"+ratio);
		System.out.println ();
		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
			    (new EuropeanOptionVG_ISandRQMC( s0,  K2,  mu,  sigma,  r,  nu, theta), streamRQMC, p, rand, n, m));
		
  }
}