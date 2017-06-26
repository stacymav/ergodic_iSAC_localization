#ifndef ERG_SETTINGS_HPP
#define ERG_SETTINGS_HPP

namespace sac {

	typedef std::vector<int> Vi;
	typedef std::vector<Vi> Vvi;
  
  
	
	/* Terrain dimensions set as [0, Lmax]. 
		Dynamics and measurements will have to be scaled down to this terrain for real systems.
		 Maximum range should be 0-1 for best results (even for rectangular real terrains) !!!*/
	const double ranges[][2] = {{0, 1}, {0, 1}};

	/*Indices of system states that correspond to the range above. Number of indices also defines dimensionality of problem.*/
	const int indices[] = {0, 1};
	const int erg_dim(sizeof(indices)/sizeof(*indices));//don't change

	/*max number of harmonics for Fourier coefficients*/
	const int K = 15;
	
	/*Ergodic cost weight*/
	const double Q_erg = 55;
	
	/*Used to place larger weight on lower frequency information - Lambda_k = 1/(1+|k|^coef_weight_decay)^((erg_dim+1)/2.0)
		The smaller this value, the more the higher frequencies contribute to the cost*/
	const double coef_weight_decay = 1.0/2.0;
	
	
	/********************************************
	Ergodic memory: Used to create memory of visited places (e.g. for time-varying EIDs)*/
	
	/*Set following to *ZERO* if you want the system to remember *EVERYTHING* from initial time.
		A value of e.g. 10 means that the system has a (rolling) memory window of 10*ts which will be equal
		to the difference between upper and lower integral limit in ck*/
	const int ergodic_memory_steps = 0;

	/*Expected Information Density 
		Bimodal gaussian for this example*/	
	const double EID_max = 35;
		
		
	Eigen::Matrix< double, erg_dim, erg_dim > Sigma;//covariance matrix
	Eigen::Matrix< double, erg_dim, erg_dim > Sigma_inv;//inverse covariance matrix
	Eigen::Matrix< double, erg_dim, 1 > mu;
	
	Eigen::Matrix< double, erg_dim, erg_dim > Sigma2;//covariance matrix
	Eigen::Matrix< double, erg_dim, erg_dim > Sigma_inv2;//inverse covariance matrix
	Eigen::Matrix< double, erg_dim, 1 > mu2;
	
	
	void initialize_cov_means() {
		Sigma << pow(0.07,2), 0,        
				 0, pow(0.07,2); 
		Sigma_inv = Sigma.inverse();
	}

	inline void EID( const double t, Eigen::Matrix< double, erg_dim, 1 > & x, double & eid_t ) {
		mu(0,0) = 0.2*cos(t/7.5) + 0.5;
		mu(1,0) = 0.2*sin(t/7.5) + 0.5;
		eid_t = (1.0/sqrt(pow(2.0*PI, erg_dim)*Sigma.determinant()))*exp(-0.5*(x-mu).transpose()*Sigma_inv*(x-mu));
	}
	 
	

}

#endif  // ERG_SETTINGS_HPP
