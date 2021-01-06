// TypesVG
package varianceGamma;

import java.io.IOException;

import umontreal.ssj.probdist.BetaDist;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class TypesVG {
	double S_bgss;
	double S_bgbs;
	double S_dgbs;
	double S_bgss_moy;
	double S_bgbs_moy;
	double S_dgbs_moy;
	double mu,mu_plus, mu_minus;
	double nu,nu_plus, nu_minus;
	double lambda, lambda_plus, lambda_minus;
	double alpha,alpha_plus, alpha_minus;
	double thetha;
	double sigma;
	double[] time;
	int n;
	int T;
	RandomStream stream1= new MRG32k3a();
	RandomStream stream2= new MRG32k3a();
	RandomStream stream3= new MRG32k3a();
	RandomStream stream4= new MRG32k3a();
	RandomStream stream5= new MRG32k3a();
	RandomStream stream6= new MRG32k3a();
	
	public TypesVG() {
		mu=1.0;
		nu=0.3;
		alpha=Math.pow(mu, 2)/nu;
		lambda=mu/nu;
		thetha=-0.1436;
		sigma=0.12136;
		
		mu_plus=(Math.sqrt(Math.pow(thetha, 2)+2.0*Math.pow(sigma, 2)/nu)+thetha)/2.0;
		nu_plus=Math.pow(mu_plus, 2)*nu;
		alpha_plus=Math.pow(mu_plus, 2)/nu_plus;
		lambda_plus=mu_plus/nu_plus;
		
		mu_minus=(Math.sqrt(Math.pow(thetha, 2)+2.0*Math.pow(sigma, 2)/nu)-thetha)/2.0;
		nu_minus=Math.pow(mu_minus, 2)*nu;
		alpha_minus=Math.pow(mu_minus, 2)/nu_minus;
		lambda_minus=mu_minus/nu_minus;
		
		T=5;
		n=1000;
		time=new double[T];
		time[0]=0.0;
		time[1]=0.6;
		for (int i=2; i<T; i++) {
			time[i]=time[i-1]+0.2;
		}
		
	}
	public double simulate_bridgeGamma(
			double t,double t1, double t2,
			double alpha,RandomStream stream) {
		double simulation=new BetaDist(alpha*(t-t1),
				alpha*(t2-t)).inverseF(
						stream3.nextDouble());;
		return simulation;
	}
	public double simulate_sequentielGamma(
			double t,double t0, double alpha,
			double lambda,RandomStream stream) {
		double simulation=new GammaDist(alpha*(t-t0),
				lambda).inverseF(stream.nextDouble());;
		return simulation;
	}
	
	public double BGSS() {
		double []tau = new double[T];
		double []Y= new double [T];
		double S_bgss = 0.0;
		Y[0]=0;
		tau[0]=0;
		for (int i=1; i<T; i++) {
			tau[i]=tau[i-1]+new GammaDist(alpha*(time[i]-time[i-1]),
					lambda).inverseF(stream1.nextDouble());
			Y[i]= Y[i-1]+NormalDist.inverseF(thetha*(tau[i]-tau[i-1]),
					sigma * Math.sqrt(tau[i]-tau[i-1]), stream2.nextDouble());
			S_bgss+=Y[i];
		}
		return S_bgss;
	}
	public double simulateRunsBGSS(int n,  Tally statBGSS) {
		statBGSS.init();
		S_bgss_moy=0.0;
		for (int i=0; i<n; i++) {
			S_bgss=BGSS();
			statBGSS.add(S_bgss);
			S_bgss_moy+=S_bgss;
		}
		return 	S_bgss_moy;	
    }
	public double BGBS() {
		double []tau = new double[T];
		double []Y= new double [T];
		double S_bgss = 0.0;
		Y[0]=0;
		tau[0]=0;
		//en t4
		tau[4]=tau[0]+new GammaDist(alpha*(time[4]-time[0]),
				lambda).inverseF(stream3.nextDouble());
		Y[4]= Y[0]+NormalDist.inverseF(thetha*(tau[4]-tau[0]),
				sigma * Math.sqrt(tau[4]-tau[0]), stream4.nextDouble());
		//en t2
		tau[2]=tau[0]+(tau[4]-tau[0])*
				new BetaDist(alpha*(time[2]-time[0]),
						alpha*(time[4]-time[2])).inverseF(
								stream3.nextDouble());
		Y[2]= NormalDist.inverseF(Y[0]+(tau[2]-tau[0])*(Y[4]-Y[0])/(tau[4]-tau[0]),
				Math.sqrt((tau[2]-tau[0])*(tau[4]-tau[2])/(tau[4]-tau[0]))*sigma,
				stream4.nextDouble());
		//en t1
		tau[1]=tau[0]+(tau[2]-tau[0])*
				new BetaDist(alpha*(time[1]-time[0]),
						alpha*(time[2]-time[1])).inverseF(
								stream3.nextDouble());
		Y[1]=NormalDist.inverseF(Y[0]+(tau[1]-tau[0])*(Y[2]-Y[0])/(tau[2]-tau[0]),
				Math.sqrt((tau[1]-tau[0])*(tau[2]-tau[1])/(tau[2]-tau[0]))*sigma,
				stream4.nextDouble());
		//en t3
		tau[3]=tau[2]+(tau[4]-tau[2])*
				new BetaDist(alpha*(time[3]-time[2]),
						alpha*(time[4]-time[3])).inverseF(
								stream3.nextDouble());
		Y[3]=NormalDist.inverseF(Y[2]+(tau[3]-tau[2])*(Y[4]-Y[2])/(tau[4]-tau[2]),
				Math.sqrt((tau[3]-tau[2])*(tau[4]-tau[3])/(tau[4]-tau[2]))*sigma,
				stream4.nextDouble());
		S_bgbs=Y[1]+Y[2]+Y[3]+Y[4];
		return S_bgbs;
	}
	public double simulateRunsBGBS(int n,  Tally statBGSS) {
		statBGSS.init();
		S_bgbs_moy=0.0;
		for (int i=0; i<n; i++) {
			S_bgss=BGBS();
			statBGSS.add(S_bgbs);
			S_bgbs_moy+=S_bgbs;
		}
		return 	S_bgbs_moy;	
    }
	
	public double DGBS() {
		double []tau_plus = new double[T];
		double []tau_minus = new double[T];
		double []Y= new double [T];
		double S_bgss = 0.0;
		tau_plus[0]=0.0;
		tau_minus[0]=0.0;
		Y[0]=0.0;
		//en t4
		tau_plus[4]=tau_plus[0]+new GammaDist(alpha_plus*(time[4]-time[0]),
				lambda_plus).inverseF(stream5.nextDouble());
		tau_minus[4]=tau_minus[0]+new GammaDist(alpha_minus*(time[4]-time[0]),
				lambda_minus).inverseF(stream6.nextDouble());
		Y[4]=tau_plus[4]-tau_minus[4];
		
		
		//en t2
		
		tau_plus[2]=tau_plus[0]+(tau_plus[4]-tau_plus[0])*
				new BetaDist(alpha_plus*(time[2]-time[0]),
						alpha_plus*(time[4]-time[2])).inverseF(
								stream5.nextDouble());
		tau_minus[2]=tau_minus[0]+(tau_minus[4]-tau_minus[0])*
				new BetaDist(alpha_minus*(time[2]-time[0]),
						alpha_minus*(time[4]-time[2])).inverseF(
								stream6.nextDouble());
		Y[2]=tau_plus[2]-tau_minus[2];
		//en t1
		tau_plus[1]=tau_plus[0]+(tau_plus[2]-tau_plus[0])*
				new BetaDist(alpha_plus*(time[1]-time[0]),
						alpha_plus*(time[2]-time[1])).inverseF(
								stream5.nextDouble());
		tau_minus[1]=tau_minus[0]+(tau_minus[2]-tau_plus[0])*
				new BetaDist(alpha_minus*(time[1]-time[0]),
						alpha_minus*(time[2]-time[1])).inverseF(
								stream6.nextDouble());
		Y[1]=tau_plus[1]-tau_minus[1];
		//en t3
		tau_plus[3]=tau_plus[2]+(tau_plus[4]-tau_plus[2])*
				new BetaDist(alpha_plus*(time[3]-time[2]),
						alpha_plus*(time[4]-time[3])).inverseF(
								stream5.nextDouble());
		tau_minus[3]=tau_minus[2]+(tau_minus[4]-tau_minus[2])*
				new BetaDist(alpha_minus*(time[3]-time[2]),
						alpha_minus*(time[4]-time[3])).inverseF(
								stream6.nextDouble());
		Y[3]=tau_plus[3]-tau_minus[3];
		S_dgbs=Y[1]+Y[2]+Y[3]+Y[4];
		return S_dgbs;
	}
	public double simulateRunsDGBS(int n,  Tally statDGBS) {
		statDGBS.init();
		S_dgbs_moy=0.0;
		for (int i=0; i<n; i++) {
			S_dgbs=DGBS();
			statDGBS.add(S_dgbs);
			S_dgbs_moy+=S_dgbs;
		}
		return 	S_bgbs_moy;	
    }
	public static void main (String[] args)throws IOException {
		TypesVG vg =new TypesVG();
		
		Tally statBGSS = new Tally("BGSS");
		vg.simulateRunsBGSS(vg.n, statBGSS);
		statBGSS.setConfidenceIntervalStudent();
	    System.out.println(statBGSS.report(0.95, 3));
	    
	    Tally statBGBS = new Tally("BGBS");
		vg.simulateRunsBGBS(vg.n, statBGBS);
		statBGBS.setConfidenceIntervalStudent();
	    System.out.println(statBGBS.report(0.95, 3));
	    
	    Tally statDGBS = new Tally("DGBS");
		vg.simulateRunsDGBS(vg.n, statDGBS);
		statDGBS.setConfidenceIntervalStudent();
	    System.out.println(statDGBS.report(0.95, 3));
	}
	
	

}
