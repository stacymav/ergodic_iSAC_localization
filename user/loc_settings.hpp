#ifndef LOC_SETTINGS_HPP
#define LOC_SETTINGS_HPP

namespace sac {

	/*Number of parameters that are being measured (equal to length of measurement vector)*/
	const int meas_vars = 2; 
	
	/*Number of parameters to be estimated for each target (at least equal to erg_dim - could be larger if you want to use the filter to estimate else too)*/
	const int local_num = erg_dim; 
		
	/*Number of targets to be localized*/
	const int no_targets = 2;
	
   /* Update EID every "update_steps" time steps */
	const int update_steps = 1;

   /* Take measurements every 'meas_freq' seconds*/
	const double meas_freq = ts;
	
	/* How far the quadrotor can take measurements */
	const double range_radius = 0.2;	
	
	/* This parameter can be used to slow down convergence of the EKF in case you are running the original dynamics of the system without scaling them down
		to the ergodic terrain. By doing that you can "simulate" a more realistic task. Set equal to 1 for no scaling.*/
	const double Kscale = 2.5;
	
	/* Indexes of measured vars to be angle wrapped - leave brace blank if no measurements are to be wrapped */
	const int meas_wrap_loc[] = {0, 1};
	const int meas_wrap_size(sizeof(meas_wrap_loc)/sizeof(meas_wrap_loc[0]));//do not change
	
	/*Measurement model*/
	inline void Y( Eigen::Matrix< double, local_num, 1 > & a, state_type & x, Eigen::Matrix< double, meas_vars, 1 > &result ) {
		result << atan2 (x[0]-a(0,0), x[1]-a(1,0)),        
				atan2(x[2], sqrt((x[1]-a(1,0))*(x[1]-a(1,0))+(x[0]-a(0,0))*(x[0]-a(0,0)))); 
	}
	
	/*Derivative of Measurement model (used in calculation of fisher information)*/
	inline void dY_da( Eigen::Matrix< double, local_num, 1 > & a, state_type & x, Eigen::Matrix< double, meas_vars, local_num > &dyda ) {
		if ((std::abs(x[0]-a(0,0))<epsilon) && (std::abs(x[1]-a(1,0))<epsilon)) {
			dyda = Eigen::Matrix< double, meas_vars, local_num>::Zero(); 
		} 
		else {
			double temp1 = ((a(0,0) - x[0])*(a(0,0) - x[0]) + (a(1,0) - x[1])*(a(1,0) - x[1]));
			dyda(0,0) =   (a(1,0) - x[1])/temp1;
			dyda(0,1) =  -(a(0,0) - x[0])/temp1;
			temp1 = (2*sqrt(pow(a(0,0) - x[0], 2.0) + pow(a(1,0) - x[1], 2.0))*(x[2]*x[2] + pow(a(0,0) - x[0], 2.0) + pow(a(1,0) - x[1], 2.0)));
			dyda(1,0) = -(x[2]*(2*a(0,0) - 2*x[0]))/temp1;
			dyda(1,1) = -(x[2]*(2*a(1,0) - 2*x[1]))/temp1;
			}
	
	}	
	
	
	/*Function that calculates the motion (if any) of targets (unknown to the ergodic algorithm)*/
	inline void a_now_i(const double t, const int i, state_type & anow ) {		
		//Pick the appropriate function depending on value of i
		switch ( i ) {
			case 0: //XYZ coords of target  
				anow[0] = 0.2*cos(2*PI/45.0*t)+0.5;
				anow[1] = 0.5+0.2*sin(2.0*PI/45.0*t);
				anow[2] = 0.2;
				break;
			case 1:            
				anow[0] = 0.2*cos(2.0*PI/45.0*(t-45.0/2.0))+0.5;
				anow[1] = 0.5+0.2*sin(2.0*PI/45.0*(t-45.0/2.0));
				anow[2] = 0.2;
				break;
			default:            
				std::cout<<"Error, bad input, quitting\n";
			break;
		}		
	}
	

	Eigen::Matrix< double, local_num, local_num > meas_Sigma;//measurement noise	
	Eigen::Matrix< double, local_num, local_num > meas_Sigma_inv;//inverse measurement noise	
	inline void initialize_meas_noise( ) {		
		meas_Sigma << 0.01, 0.0,        
					0.0, 0.01; 	

		meas_Sigma_inv = meas_Sigma.inverse();			
	}
	
	
	/*Function called when targets are created. Initializes the transition model and transition noise of a target etc*/
	inline void target_params(const int i, 	Eigen::Matrix< double, local_num, 1 > &mu, Eigen::Matrix< double, local_num, local_num > &Sigma, 
							Eigen::Matrix< double, local_num, local_num > & A_f, Eigen::Matrix< double, local_num, local_num >  &trans_sigma ) {		
		//Pick the appropriate function depending on value of i
		switch ( i ) {
			case 0: 
				trans_sigma << 0.001, 0.0,        
								0.0, 0.001;
				//model moving targets by a diffusion process
				A_f << 1.0, 0.0,        
						0.0, 1.0; 
				//Initial values for belief params
				Sigma << pow(10.0,2.0), 0.0,        
						0.0, pow(10.0,2.0); 
				mu << 0.5, 0.5;
				break;
			case 1:            
				trans_sigma << 0.001, 0.0,        
								0.0, 0.001;
				//model moving targets by a diffusion process
				A_f << 1.0, 0.0,        
						0.0, 1.0; 
				//Initial values for belief params
				Sigma << pow(10.0,2.0), 0.0,        
						0.0, pow(10.0,2.0); 
				mu << 0.5, 0.5;
				break;
			default:            
				std::cout<<"Error, bad input, quitting\n";
			break;
		}		
	}

}

#endif  // LOC_SETTINGS_HPP
