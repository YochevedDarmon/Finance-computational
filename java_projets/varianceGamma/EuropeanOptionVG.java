package varianceGamma;
// EuropeanOptionVG {
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.GammaAcceptanceRejectionGen;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.mcqmctools.*;

public class EuropeanOptionVG implements MonteCarloModelDouble {

	double mu = -0.1436; 			
	double sigma = 0.12136; 		
	double r = 0.1; 				
	double nu = 0.2;	
	double K1 = 130.0;
	double K2 = 160.0;
	double K;
	double s0 = 100.0;
	double omega = Math.log (1 - mu*nu - sigma*sigma*nu / 2.0) / nu;
	int n = 16 * 1024;  // 2^14  for Monte Carlo
	
	double theta1=11.9005607518060;
	double theta2=23.5189440032670;
	double theta;
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
	double payoff;   // Value of X to return.
	RandomVariateGen genYp;
	RandomVariateGen genYm;
	RandomVariateGen genYm_mc;
	RandomStream noise = new MRG32k3a();
	ContinuousDistribution distgp;
	ContinuousDistribution distu;
	ContinuousDistribution distpt;

	// Constructor.
	public EuropeanOptionVG(double s0, double K, double mu, double sigma, double r, double nu, double theta ) {
		this.s0 = s0;
		this.theta=theta;
		this.K=K;
		genYp = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdaplus)); 
		genYm = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdamoins+theta));
		genYm_mc = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdamoins));
		
		
		distgp = new GammaDist(alphaplus, lambdaplus);
	}

	// Generates payoff X, without IS.
	public void simulate(RandomStream stream) {
		/*payoff = 0.0;
		double gm= genYm_mc.nextDouble();
		double S=s0*Math.exp(r+omega+genYp.nextDouble()-gm);
		if (S-K> 0)
			payoff = Math.exp(-r)*Math.max(0, (S - K));*/
		RandomVariateGen genYp = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdaplus)); 
		RandomVariateGen genYm = new GammaAcceptanceRejectionGen
		          (noise, new GammaDist (alphaplus, lambdamoins));
		double S=s0*Math.exp(r+omega+genYp.nextDouble()-genYm.nextDouble());
		double y= Math.exp(-r)*Math.max(S-K, 0.0);
		payoff= y;
			
	}
	

	// Returns payoff X
	public double getPerformance() {
		return payoff;
	}

	// Descriptor of model
	public String toString() {
		return " financial option with barriers, MC";
	}
	
}