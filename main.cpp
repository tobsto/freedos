// ****************************************************************************
// Tobias Stollenwerk, stollenwerk@th.physik.uni-bonn.de
// ****************************************************************************

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<fstream>
#include<vector>
#include<string>
#include<complex>
#include<limits>
#include"grid.h"
#include"multigrid.h"
#include"mesh.h"
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

// redefine exponetial because of limited interval (double precision)
double  expo(double x)
{
	double limit=700;
	if (x>700)  
	{
		return exp(700);
	}  
	else if (x<-700) 
	{
		return 0;
	}
	else
	{
		return exp(x);
	} 
}
// Fermi function
double  fermi(double z, double beta)
{
        double x=beta*z;
        if (x>700)
        {
                return 0;
        }
        else if (x<-700)
        {
                return 1;
        }
        else 
        {
                return 1/(exp(x)+1);
        }
}
// bare conduction electron density of states
double dos2d(double omega, double center, double width)
{	
	// box without exponential tails
	if (fabs(omega-center)>width)
	{
		return 0;
	}
	else if (fabs(omega-center)<=width)
	{
		return 1.0/(2.0*width);
	}
}
double dos3d(double omega, double center, double width)
{	
	// half circle
	if (1-pow(((omega-center)/width),2)<=0)
	{
		return 0;
	}
	else if (1-pow(((omega-center)/width),2)>0)
	{
		return 2.0/(M_PI*width) * sqrt(1-pow(((omega-center)/width),2));
	}
}

double dn(double & x, double & beta, double & n_cc, double & D0)
{
	// Integration grid
	double domega_max=0.005;
	double omegalog1=0.1;
	mesh egrid;
	egrid.add_gr_equi(x-D0-2*omegalog1, x+D0+2*omegalog1, domega_max);
	// resolve band edges
	egrid.add_sgr_log(x-D0, 0.1, 1E-7, domega_max);
	egrid.add_sgr_log(x+D0, 0.1, 1E-7, domega_max);
	// position of the maximal point density in the Fermi function (Maximum of the 2nd derivative)
	double omegab=1.31696/beta;
	// maximal and minimal resolution for the grid which resolves the Fermi function (fermi grid)
	double domegab_min=4.0/(400*beta);
	// resolve fermi edge
	egrid.add_sgr_log(-omegab, omegalog1, domegab_min, domega_max, "fermi_l");
	egrid.add_sgr_log( omegab, omegalog1, domegab_min, domega_max, "fermi_r");
	egrid.add_lendpoint(x-D0);
	egrid.add_rendpoint(x+D0);
	egrid.create();

	double result=0.0;
	for (int n=0; n<=egrid.M; n++)
	{
		result+=fermi(egrid.omega[n], beta)*dos3d(egrid.omega[n], x, D0)*egrid.domega[n];
	}
	result-=n_cc;
	return result;
}

void plotFunctions(double & x, double & beta, double & n_cc, double & D0)
{

	// Integration grid
	double domega_max=0.005;
	double omegalog1=0.1;
	mesh egrid;
	egrid.add_gr_equi(x-D0-2*omegalog1, x+D0+2*omegalog1, domega_max);
	// resolve band edges
	egrid.add_sgr_log(x-D0, 0.1, 1E-7, domega_max);
	egrid.add_sgr_log(x+D0, 0.1, 1E-7, domega_max);
	// position of the maximal point density in the Fermi function (Maximum of the 2nd derivative)
	double omegab=1.31696/beta;
	// maximal and minimal resolution for the grid which resolves the Fermi function (fermi grid)
	double domegab_min=4.0/(400*beta);
	// resolve fermi edge
	egrid.add_sgr_log(-omegab, omegalog1, domegab_min, domega_max, "fermi_l");
	egrid.add_sgr_log( omegab, omegalog1, domegab_min, domega_max, "fermi_r");
	egrid.add_lendpoint(x-D0);
	egrid.add_rendpoint(x+D0);
	egrid.create();

	ofstream out_fermi("fermi.dat");
	ofstream out_dos("dos.dat");
	double result=0.0;
	for (int n=0; n<=egrid.M; n++)
	{
		out_fermi << scientific << setprecision(15) << egrid.omega[n] << "\t" << fermi(egrid.omega[n], beta) << endl;
		out_dos   << scientific << setprecision(15) << egrid.omega[n] << "\t" << dos3d(egrid.omega[n], x, D0) << endl;
	}
	out_fermi.close();
	out_dos.close();
}

void plotDn(double & beta, double & n_cc, double & D0)
{
	ofstream out ("dn.dat");
	for (double d=-1.5; d<=1.5; d+=0.001)
	{
		out << scientific << setprecision(15) << d << "\t" << dn(d, beta, n_cc, D0) << endl;
	}
	out.close();
}
class xTolerance
{
	public:
	double tol;
	xTolerance(double t)
	{
		this->tol=t;
	}
};

void getDelta(double & x, double & beta, double & n_cc, double & D0)
{
	double ftol=1E-8;
	double xa, xb, xc, fa, fb, fc;
	double dx0=0.1;
	double dx=dx0;
	double dx1=0.1;
	// precision
	double xtol=numeric_limits<double>::epsilon()*10*D0;
	int max1=100;
	int max2=100;
	double lowerBound=-10*D0+1.0;
	double upperBound= 10*D0-1.0;

	xa=x;
	fa=dn(xa, beta, n_cc, D0);
	//ofstream out ("delta_search.dat");
	//cout << "\t* mu:\t" << xa << "\t" << "delta_n:\t" << fa << endl;
	//out << scientific << setprecision(15) << xa << "\t" << fa << endl;

	if (fabs(fa)>ftol)
	{
		xb=xa - fa/abs(fa) * dx;
		fb=dn(xb, beta, n_cc, D0);
		//cout << "\t* mu:\t" << xb << "\t" << "delta_n:\t" << fb << endl;
		//out << scientific << setprecision(15) << xb << "\t" << fb << endl;
		if (fabs(fb)<=ftol)
		{
			x=xb;
			return;
		}	

		int counter=0;
		bool constFuncFlag=false;
		double factor=1.5;
		while (fa*fb>0.0 && counter<max2 && fabs(fb)>ftol)
		{
			if (fabs(fb-fa)>xtol)
			{
				xc=xb - fb * (xb-xa)/(fb-fa);
				constFuncFlag=false;
			}
			else
			{
				constFuncFlag=true;
			}
			if (constFuncFlag || xc>=upperBound || xc<=lowerBound)
			{
				xc=xb - fb/abs(fb) * dx1;
				while (xc>=upperBound)
				{
					xc-=dx1;
				}
				while (xc<=lowerBound)
				{
					xc+=dx1;
				}
			}
			
			fc=dn(xc, beta, n_cc, D0);
			//cout << "\t* mu:\t" << xc << "\t" << "delta_n:\t" << fc << endl;
			//out << scientific << setprecision(15) << xc << "\t" << fc << endl;
	
			xa=xb;
			fa=fb;
			xb=xc;
			fb=fc;
			counter++;
		}
		if (counter==max2)
		{
			// Did not found two points with different signs
			if (fa<0.0)
			{
				x=1.0;
				return;
			}
			else
			{
				x=-1.0;
				return;
			}
		}
		if (fabs(fb)<=ftol)
		{
			x=xb;
			return;
		}
			
		// inital values for bounded secant (bisection) method
		double x1, x2, f1, f2, f;
		if (fa<fb)
		{
			x1=xa;
			f1=fa;
			x2=xb;
			f2=fb;
		}
		else
		{
			x1=xb;
			f1=fb;
			x2=xa;
			f2=fa;
		}

		// container for difference in function values above and below the x-axis 
		vector<double> df;
		// container storing if element of df is below or above x-axis
		vector<bool> below;
		double avdfp;
		double avdfm;
		double fp_before=f2;
		double fm_before=f1;
		bool contFlag=true;
		int N_av=10;
		double df_tol=1E-6;

		// bisection method to find root
		counter=0;
		do
		{
			if (fabs(f2-f1)>xtol)
			{
				// secant method
				x=x1-f1*(x2-x1)/(f2-f1);
			}
			// if there is a division by 0 due to equal function values f1 and f2, use bisection
			else
			{
				x=0.5*(x1+x2);
			}
			// get function value
			f=dn(x, beta, n_cc, D0);
			//cout << "\t* mu:\t" << x << "\t" << "delta_n:\t" << f << endl;
			//out << scientific << setprecision(15) << x << "\t" << f << endl;


			// Stop iteration for the case of switching between two function values above and below the x-axis
			// this happens, since the numerically evaluated function can not be calculated with arbitrary 
			// accuracy. The algorithm stops if the difference of the function values below and above the x-axis
			// have only tiny variations (<1E-6) over some iterations (10).
			// Therefore one has to keep track of the difference in the function values below (df[0]) and above (df[1])
			// the x-axis and calculate the average over the last iterations (avdfm and avdfp).
			if (f<=0)
			{
				below.push_back(true);
				df.push_back(fabs(f-fm_before));
				fm_before=f;
			}
			else
			{
				df.push_back(fabs(f-fp_before));
				below.push_back(false);
				fp_before=f;
			}
	

			if (counter>=N_av)
			{
				// calculate average of difference in function values over the
				// N_av iterations. Above and below the x-axis. If both values
				// are below df_tol, stop iteration and set root to the last 
				// x value of that branch (above or below) which occured more
				// often in the last N_av iterations. 
				avdfm=0.0;
				avdfp=0.0;
				int N_below=0;
				for (int i=df.size()-N_av; i<df.size(); i++)
				{
					if (below[i])
					{
						avdfm+=df[i];
						N_below++;
					}
					else
					{
						avdfp+=df[i];
					}
				}
				avdfm/=max(N_below,1);
				avdfp/=max(N_av-N_below,1);

				if (avdfm<=df_tol && avdfp<=df_tol)
				{
					contFlag=false;	
				}
			}
			if (contFlag==false)
			{
				//cerr << endl;
				//cerr << "Warning: Error in fixMu routine. Search for root was not successful." << endl;
				//cerr << "Difference in x points fell under machine precision." << endl;
				throw xTolerance(f);
			}
			if (f*f1<0.0)
			{
				x2=x;
				f2=f;
			}
			else
			{
				x1=x;
				f1=f;
			}
			counter++;	
		}		
		while (fabs(f)>ftol && counter<max2 && fabs(x2-x1)>=xtol);

		if (counter==max2)
		{
			cerr << endl;
			cerr << "Warning: Error in fixMu routine. Search for root was not successful."<< endl;
			cerr << "Not successful after " << counter << " iterations. Tolerance " << ftol << " to small?" << endl;
			throw xTolerance(f);
		}
		if (fabs(x2-x1)<xtol)
		{
			//cerr << endl;
			//cerr << "Warning: Error in fixMu routine. Search for root was not successful." << endl;
			//cerr << "Difference in x points fell under machine precision." << endl;
			throw xTolerance(f);
		}
	}
}

int main(int argc, char * argv[])
{
	// Boltzmann constant in eV/Kelvin
	double k_b=8.617343E-5;
	// parameter
	double T;
	double n_cc;
	double D0;
	// read in parameters
	try 
	{
		// Declare a group of options that will be 
		// allowed only on command line
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h", "produce help message")
		;
		
		// Declare a group of options that will be 
		// allowed both on command line and in
		// config file
		po::options_description options("Standard options");
		options.add_options()
			("temperature,t"     , po::value<double>(       &T)->default_value(40  )       , "Temperature in Kelvin")
			("occupation,n"      , po::value<double>(    &n_cc)->default_value(0.01)       , "Occupation")
			("half_band_width,D" , po::value<double>(      &D0)->default_value(1.0)        , "Half band width")
		;

		po::options_description all("Allowed options");
		all.add(generic).add(options);
		
		po::options_description cmdline_options;
		cmdline_options.add(generic).add(options);
		
		po::options_description visible("Allowed options");
		visible.add(generic).add(options);
		
		po::variables_map vm;
		store(parse_command_line(argc, argv, all), vm);
		notify(vm);
	}
	catch(exception& e)
	{
		cerr << "Something went wrong during the parameter input. Exception:" << endl;
		cerr << e.what() << endl;
		cerr << "Break." << endl;
		exit(1);
	}
	// Inverse temperature (dimensionless)
	double beta=D0/(k_b*T);
	double Delta=0.0;
	//plotDn(beta, n_cc, D0);
	//plotFunctions(Delta, beta, n_cc, D0);
	getDelta(Delta, beta, n_cc, D0);
	cout << "Delta is: " << Delta << endl;

	ofstream out ("delta_vs_ncc_left.dat");
	for (double n=0.00; n<=0.0005; n+=5E-7)
	{
		Delta=0.0;
		try
		{
			getDelta(Delta, beta, n, D0);
		}
		catch(xTolerance & xtol)
		{
			cerr << "Warning: getDelta: Occupation number accuracy was " << fabs(xtol.tol) << endl;
		}
		catch(...)
		{
			cerr << "Error: getDelta_r: Something went wrong. Break." << endl;
			exit(1);
		}
		out << scientific << setprecision(15) << n << "\t" << Delta << endl;
		cout << scientific << setprecision(15) << n << "\t" << Delta << endl;
	}
	out.close();
}
