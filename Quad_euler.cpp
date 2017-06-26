#include <master.hpp>               // Master include file

#include <ctime> //for timing

#include <boost/random/normal_distribution.hpp> //for timing
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

using namespace sac;

/* iSACstep() function operator
  input:  initial state and time 
  return: Does not explicitly return anything but the following fields of class "sac_step" can be accessed
  
  iSACstep.xnext - integrated state at time t0+ts

  iSACstep.u_switch - vector of SAC action values applied from [t_i, t_f] which is a subset of [t0, t0+ts].
          If [t_i, t_f] is not equal to [t0, t0+ts] then the default control is applied over the remaining interval. 
  t_i - initial time for application of the control.  t0 <= t_i <= t0+ts
  t_f - final time for control application.  t0 <= t_f <= t0+ts

  WARNING: iSACstep.u_switch is only applied when t_f-t_i > 0, otherwise u_default is applied.
  WARNING: If [t_i, t_f] is not equal to [t0, t0+ts] then u_default is applied 
           over the remaining interval.
  NOTE: for speed return and input types should be changed and passed as
        references / pointers
*/


int main(int /* argc */ , char** /* argv */ )
{
	using namespace std;
	/*********************************************/
	/* Vars etc*/
	std::vector<target> mytargets; //empty vector that will hold target instances
	for (int i = 0; i < no_targets; ++i) {
		mytargets.push_back(target( i )); }// ith element is a copy of this
	isac_step iSACstep(mytargets);//instance
	Eigen::Matrix< double, meas_vars, 1 > z;//measurement vector
	Eigen::Matrix< double, local_num, 1 > anow;
	state_type x0(xlen);
	state_type x_temp(erg_dim), a_now(3);
	Eigen::Matrix< double, xlen, 1 > xnext;//for prints
	int i, indx_; double temp;
	ofstream myfile;
  	myfile.open ("./data/states.csv");//open file to save stuff
	
	// This is the underlying integer random number generator
	boost::mt19937 igen;
	// The second template parameter is the actual floating point
	// distribution that the user wants
	boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
	gen(igen,
	boost::normal_distribution<>(0,sqrt(0.01)));

	/*********************************************/
	/* Initializations*/
	//state	
	for (i = 0; i < xlen; ++i) { x0[i] = x_init[i]; }
		

	clock_t begin = clock();//for timing

	/*********************************************/
	/* Receding horizon loop*/
	for (double t0 = t_init; t0 < t_final; t0 = t0 + ts)
	{
		/* Perform SAC iteration - updates: J0, Jn, u, x_intp */
		iSACstep( t0, x0 );

		//update state
		for (i=0; i < xlen; ++i) { x0[i] = iSACstep.xnext[i]; }	
		
		
		
		/* Take measurements */
		for (i = 0; i < no_targets; ++i) {	
			z(0,0) = 0; z(1,0) = 0;//initialize
			a_now_i(t0+ts, i, a_now );//calculate target position (unknown to the algorithm)
			anow(0,0) = a_now[0]; anow(1,0) = a_now[1];			
			x_temp[0] = a_now[0] - x0[0]; x_temp[1] = a_now[1] - x0[1];			
			calc_norm( temp, x_temp );
			if(temp <= range_radius) {				
				Y( anow, x0, z );//call measurement model	
				z(0,0) = z(0,0) + gen(); z(1,0) = z(1,0) + gen();
				//Anglewrap measurement if applicable
				for (indx_ = 0; indx_ < meas_wrap_size; ++indx_ ) {AngleWrap( z(meas_wrap_loc[indx_],0) );} // Angle wrapping (if any)
				//		
				mytargets[i].EKF( x0, z, true ); //run EKF for this target			
			}
			else
			{
				mytargets[i].EKF( x0, z, false ); //run EKF for this target
			}
			
			//std::cout << i << " " << anow.transpose() << "\n";
			//std::cout << i << " " << z.transpose() << "\n";
		}		
		
		
		//Prints
		State2Mat( iSACstep.xnext, xnext ); // convert state to matrix form to be able to print state directly
		myfile << t0 << " " << mytargets[0].mu_.transpose()<<" "<< mytargets[0].Sigma_(0,0) << " " << mytargets[1].mu_.transpose()<<" "<< mytargets[1].Sigma_(0,0) << " " <<
		iSACstep.xnext[0] << " " << iSACstep.xnext[1] <<"\n";//write to file
		//cout << t0 << ", " << iSACstep.xnext[0] << ", " << iSACstep.xnext[2] << "\n";
		cout << t0 << "\n";
		
		for(i=0; i<erg_dim; ++i) {
			if((iSACstep.xnext[indices[i]] > ranges[i][1])||(iSACstep.xnext[indices[i]] < ranges[i][0])) {
				std::cout << "Out of bounds" <<"\n";
				system("pause");}
		}

		
		//system("pause");
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Elapsed time (s): " << elapsed_secs << "\n";//print elapsed time


	myfile.close();//close file



	return 0;
}
