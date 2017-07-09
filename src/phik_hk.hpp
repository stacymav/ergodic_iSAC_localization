/*
    ergodic_iSAC_localization: Real-time receding horizon ergodic localization of multiple targets. 
    Copyright (C) 2017 Anastasia Mavrommati
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PHIK_HK_HPP
#define PHIK_HK_HPP

namespace sac {

  //[ The class that calculates phik and ck
  class phik_hk {
	std::vector<target> &r_mytargets_;
	size_t i_, j_, o_, l_, m_, n_;
  	const size_t grid;	
	Vi outputTemp;
	std::vector<state_type> integrands;//used in trapezoid		
	Vi counters; //set counters to zero
	Eigen::Matrix< double, erg_dim, 1 > x_vec_erg;
	double EID_t_max;//max value of current EID (for normalization)
	Eigen::Matrix< double, local_num, 1 > mu_temp;
	Eigen::Matrix< double, 3, 1 > x_temp;
	state_type x_temp_st;
	double a1, a2, dt_a1, dt_a2;
	state_type phi11a1, phi12a1, phi21a1, phi22a1;
	state_type phi11a2, phi12a2, phi21a2, phi22a2;
	state_type a1_v, a2_v;
	Eigen::Matrix< double, meas_vars, local_num > dyda_;
	Eigen::Matrix< double, local_num, local_num > PHI_temp;
	Eigen::Matrix< double, local_num, 1 > a_vector;
	double pa;
	Eigen::Matrix< double, 2,2 > FI_expected;
	state_type det_FI_all_targets;
	double temp;
	
  //protected:
  public:	
	const int no_k_comb, no_terr_comb;
	state_type EID_values_t;//current EID values
	std::vector<state_type> space;//meshgrid of terrain for calculation of hk and phik
	Vvi terr_indices;//indices of terrain (0-30)
	Vvi K_vector;//indices of terrain (0-30)
	Vvi k, all_terr_indices;//vectors that store all combinations of k and terrain indices
	state_type hk, phik;//vector of normalizing values
	
    
    phik_hk( std::vector<target> &r_mytargets ) : r_mytargets_(r_mytargets),phi11a1(6), phi12a1(6), phi21a1(6), phi22a1(6),
					phi11a2(6), phi12a2(6), phi21a2(6), phi22a2(6),x_temp_st(3), a1_v(6), a2_v(6), det_FI_all_targets(no_targets),
					 grid(31), no_terr_comb(pow(grid, erg_dim)), space(erg_dim, state_type(grid)), 
					 terr_indices(erg_dim, Vi(grid)), hk(no_k_comb), phik(no_k_comb), EID_values_t(no_terr_comb),
					 K_vector(erg_dim, Vi(K+1)), no_k_comb(pow(K+1, erg_dim)),
					 integrands(erg_dim+1, state_type(grid)), counters(erg_dim+1) { 
						initialize();			
						//Find all combinations
						cart_product(k, outputTemp, K_vector.begin(), K_vector.end());
						outputTemp.clear();
						cart_product(all_terr_indices, outputTemp, terr_indices.begin(), terr_indices.end());	
						//calculate normalizing factor for fourier coefficients
						calc_hk( );		
						//Calculate Fourier coefficients
						update_phik( t_init );
						}
			

	inline void initialize( );		
		
	inline void calc_hk( );
	
	inline void calc_norm_EID( const double t, state_type & EID_values_curr );
	
	inline void update_phik( const double t );
	
	inline void hk_integrand( double & integrand, Vi & space_indices, Vi & k_curr );
	
	inline void phik_integrand( double & integrand, const size_t i, const size_t o );
	
	inline void EID_FI( Eigen::Matrix< double, erg_dim, 1 >  x, state_type & result );

  };
  
  
  		/*Function that finds all combinations of harmonics, and indices of terrain.
			Also calculates the meshgrid values for the terrain (equivalent to linspace)*/
	inline void phik_hk::initialize( ) { 			
		for(int j=0; j < erg_dim; ++j) {				
			for (int i = 0; i < grid; ++i) {			
				//Meshgrid values
				space[j][i] = ranges[j][0] + i*(ranges[j][1]-ranges[j][0])/(grid-1.0);
				//Terrain indices
				terr_indices[j][i] = i;					
			}
			for (int i = 0; i <= K; ++i) { 
				K_vector[j][i] = i;
			}
		}
	}
		
	/*Calculate normalizing factors for Fourier coefficients*/
	inline void phik_hk::calc_hk( ) { 
		for (o_ = 0; o_ < no_k_comb; ++o_) {
			for (i_ = 0; i_ <= erg_dim; ++i_) {counters[i_] = 0;}//reinitialize counters
			for (i_ = 0; i_ < no_terr_comb; ++i_) {
				hk_integrand( integrands[0][counters[0]], all_terr_indices[i_], k[o_] );
				//std::cout << integrands[0][counters[0]] <<"\n";
				counters[0]++;
				for (j_ = 0; j_ < erg_dim; ++j_) {				    
					if (fmod(i_+1.0,pow(grid, j_+1)) < epsilon) {
						 trapezoid( integrands[j_+1][counters[j_+1]], integrands[j_], space[j_] );
						counters[j_] = 0; counters[j_+1]++;
						}
				}
			}			
			hk[o_] = pow(integrands[erg_dim][0], 1.0/erg_dim);//final result for each k vector
			//std::cout << hk[o_]  <<"\n";
		}
			
	}

	/*Calculates EID for all x in terrain and normalizes its values*/
	inline void phik_hk::calc_norm_EID( const double t, state_type & EID_values_curr ) { 
		for (l_ = 0; l_ < no_targets; ++l_) {
			r_mytargets_[l_].EID_values_t.resize(no_terr_comb);
			r_mytargets_[l_].EID_t_max = 0;//initialize
		}	
		
		for (i_ = 0; i_ < no_terr_comb; ++i_) {
			for (j_ = 0; j_ < erg_dim; ++j_) { x_vec_erg(j_, 0) = space[j_][all_terr_indices[i_][j_]]; }		
			EID_FI( x_vec_erg, det_FI_all_targets );
			for (l_ = 0; l_ < no_targets; ++l_) {
				r_mytargets_[l_].EID_values_t[i_] = det_FI_all_targets[l_];
				//Find max of current target EID
				if(r_mytargets_[l_].EID_values_t[i_] > r_mytargets_[l_].EID_t_max) { r_mytargets_[l_].EID_t_max = r_mytargets_[l_].EID_values_t[i_]; }
			}
		}
		/*Normalize EID*/
		for (i_ = 0; i_ < no_terr_comb; ++i_) { 
			EID_values_curr[i_] = 0;
			for (l_ = 0; l_ < no_targets; ++l_) {
				r_mytargets_[l_].EID_values_t[i_] = r_mytargets_[l_].EID_values_t[i_]*EID_max/r_mytargets_[l_].EID_t_max; 
				if(r_mytargets_[l_].EID_values_t[i_] > EID_values_curr[i_]) { EID_values_curr[i_] = r_mytargets_[l_].EID_values_t[i_]; }
			}			
		}
		//std::cout << "The largest element is "  << EID_t_max << '\n';
	}
	
	/*Calculate Fourier coefficients*/
	inline void phik_hk::update_phik( const double t ) {
		calc_norm_EID( t, EID_values_t );//calculate EID	
		for (o_ = 0; o_ < no_k_comb; ++o_) {
			for (i_ = 0; i_ <= erg_dim; ++i_) {counters[i_] = 0;}//reinitialize counters
			for (i_ = 0; i_ < no_terr_comb; ++i_) {
				phik_integrand( integrands[0][counters[0]], i_, o_ );
				//std::cout << integrands[0][counters[0]] <<"\n";
				counters[0]++;
				for (j_ = 0; j_ < erg_dim; ++j_) {				    
					if (fmod(i_+1.0,pow(grid, j_+1)) < epsilon) {
						 trapezoid( integrands[j_+1][counters[j_+1]], integrands[j_], space[j_] );
						counters[j_] = 0; counters[j_+1]++;
						}
				}
			}			
			phik[o_] = integrands[erg_dim][0];//final result for each k vector
			//std::cout << phik[o_]  <<"\n";
		}
	}
			
	/*Calculates integrand value for hk specified by space_indices and k and stores it in integrand.*/
	inline void phik_hk::hk_integrand( double & integrand, Vi & space_indices, Vi & k_curr ) { 
		integrand = 1.0;
		for (j_ = 0; j_ < erg_dim; ++j_) {
			integrand = integrand * pow(cos(k_curr[j_]*PI*space[j_][space_indices[j_]]/ranges[j_][1]), erg_dim);}			
	}
	
	/*Calculates integrand value for phik specified by space_indices and k and stores it in integrand.*/
	inline void phik_hk::phik_integrand( double & integrand, const size_t i, const size_t o ) { 
		integrand = 1.0/hk[o];
		for (j_ = 0; j_ < erg_dim; ++j_) {
			integrand = integrand * cos(k[o][j_]*PI*space[j_][all_terr_indices[i][j_]]/ranges[j_][1]);}			
		integrand = integrand*EID_values_t[i];
	}
	
	/*Calculates and returns EID value using Fisher Information at given point for each target.*/
	inline void phik_hk::EID_FI( Eigen::Matrix< double, erg_dim, 1 >  x, state_type & result ) { 
		
		for(m_ = 0; m_ < erg_dim; ++m_) {
			x_temp(m_,0) = x(m_,0) + 2.0*range_radius;
			x(m_,0) = x(m_,0) + 2.0*range_radius;
			x_temp_st[m_] = x_temp(m_,0);
			}
		x_temp(2,0) = des_height;
		x_temp_st[2] = x_temp(2,0);
		for (l_ = 0; l_ < no_targets; ++l_) {//loop through all targets
			for(m_ = 0; m_ < local_num; ++m_) {mu_temp(m_,0) = r_mytargets_[l_].mu_(m_,0) + 2.0*range_radius;}
			dt_a1 = 2.0*range_radius/5.0; a1 = x_temp(0,0) - range_radius;
			for(m_ = 0; m_ < 6; ++m_) { //for a1
				temp = range_radius*range_radius - pow((a1-x_temp(0,0)),2.0);
				temp = temp > 0 ? sqrt(temp) : 0;				
				dt_a2 = (2*temp + 2*x_temp(1,0))/5.0;
				

				
				if(dt_a2 > epsilon) {
					a2 = -temp + x_temp(1,0);
					for(n_ = 0; n_ < 6; ++n_) { //for a2
						a_vector(0,0) = a1; a_vector(1,0) = a2;
						dY_da( a_vector, x_temp_st, dyda_ );
						pa = r_mytargets_[l_].p_of_a( x, r_mytargets_[l_].Sigma_, mu_temp);//pa value for this a vector					
											
											
						PHI_temp = dyda_.transpose()*meas_Sigma_inv*dyda_* pa;
						phi11a2[n_] = PHI_temp(0,0);
						phi12a2[n_] = PHI_temp(0,1);
						phi21a2[n_] = PHI_temp(1,0);
						phi22a2[n_] = PHI_temp(1,1);
						/*phi11a2[n_] = (dyda_.col(0)).transpose()*meas_Sigma_inv*dyda_.col(0)* pa;
						phi12a2[n_] = (dyda_.col(0)).transpose()*meas_Sigma_inv*dyda_.col(1)* pa;
						phi21a2[n_] = (dyda_.col(1)).transpose()*meas_Sigma_inv*dyda_.col(0)* pa;
						phi22a2[n_] = (dyda_.col(1)).transpose()*meas_Sigma_inv*dyda_.col(1)* pa;*/
						a2_v[n_] = a2;
						a2 = a2 + dt_a2;
					}				
					trapezoid( phi11a1[m_], phi11a2, a2_v );
					trapezoid( phi12a1[m_], phi12a2, a2_v );
					trapezoid( phi22a1[m_], phi22a2, a2_v );
					trapezoid( phi21a1[m_], phi21a2, a2_v );
				}
				else
				{
					phi11a1[m_] = 0;
					phi12a1[m_] = 0;
					phi22a1[m_] = 0;
					phi21a1[m_] = 0;
				}
				a1_v[m_] = a1;
				a1 = a1 + dt_a1;
				
			}			
			trapezoid( FI_expected(0,0), phi11a1, a1_v );
			trapezoid( FI_expected(0,1), phi12a1, a1_v );
			trapezoid( FI_expected(1,1), phi22a1, a1_v );
			trapezoid( FI_expected(1,0), phi21a1, a1_v );	

			result[l_] = FI_expected.determinant();					
		}
	}
	
}


#endif  // PHIK_HK_HPP
