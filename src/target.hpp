#ifndef TARGET_HPP
#define TARGET_HPP

namespace sac {

  //[ The class is instantiated for each target to be localized. It must be modified by user 
  class target {
		//State transition model
		Eigen::Matrix< double, local_num, local_num > A_f_;
		Eigen::Matrix< double, local_num, local_num > trans_sigma_;//transition noise (different for each target generally, here it's the same)
		//EKF prediction vars
		Eigen::Matrix< double, local_num, local_num > Sigma_pred;//covariance
		Eigen::Matrix< double, local_num, 1 > mu_pred;
		
		Eigen::Matrix< double, meas_vars, local_num > dyda_;
		Eigen::Matrix< double, meas_vars, 1 > y_, innovation;
		Eigen::Matrix< double, local_num, meas_vars > K_;//EKF gain
		size_t indx_;

		//protected:
  public:	
		//Gaussian belief params
		Eigen::Matrix< double, local_num, local_num > Sigma_;//covariance
		Eigen::Matrix< double, local_num, 1 > mu_;
		state_type EID_values_t;//current EID values for this target
		double EID_t_max;
    
    target( int i ) { target_params(i, mu_, Sigma_, A_f_, trans_sigma_ ); }//initialize matrices for current target
			
			
	/*Calculates and returns belief for a single target (Gaussian)*/
	double p_of_a( Eigen::Matrix< double, erg_dim, 1 > & x, Eigen::Matrix< double, local_num, local_num > & Sigma_temp,
					Eigen::Matrix< double, local_num, 1 > & mu_temp) {
		return (1.0/sqrt(pow(2.0*PI, erg_dim)*Sigma_temp.determinant()))*exp(-0.5*(x-mu_temp).transpose()*Sigma_temp.inverse()*(x-mu_temp));
	}
	
	
	/*EKF filter - upfates mu_ and Sigma_*/
	void EKF( state_type & x, Eigen::Matrix< double, meas_vars, 1 > & z, bool measurement ) {
		//Prediction update
		mu_pred = A_f_*mu_;
		Sigma_pred = A_f_*Sigma_*A_f_.transpose() + trans_sigma_;
		
		//Measurement update
		if(measurement) {
			//call measurement model and its derivative
			dY_da( mu_pred, x, dyda_ );
			Y( mu_pred, x, y_ );
			K_ = Sigma_pred*dyda_.transpose()*(dyda_*Sigma_pred*dyda_.transpose() + meas_Sigma).inverse()/Kscale;
			innovation = z - y_;
			//Anglewrap innovation if applicable
			for (indx_ = 0; indx_ < meas_wrap_size; ++indx_ ) {AngleWrap( innovation(meas_wrap_loc[indx_],0) );} // Angle wrapping (if any)
			//
			mu_ = mu_pred + K_*innovation;
			Sigma_ = (Eigen::Matrix< double, local_num, local_num >::Identity(local_num, local_num) - K_*dyda_) * Sigma_pred;			
		}
		else 
		{
			mu_ = mu_pred;
			Sigma_ = Sigma_pred;
		}
	}	
  
  };
  
}

#endif  // TARGET_HPP
