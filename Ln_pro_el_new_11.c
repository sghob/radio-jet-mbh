/* Include the appropriate libraries */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>


int n_phase = 4;
/* r in kpc */
/*double a1=0.1, a2=0.05, a3=0.01, a4=0.001;*/
double a1=0.001, a2=0.001, a3=0.001, a4=0.001;
double rstart = 1.0;
double e = 0.5;
double redshift[] = {0, 1, 2, 3, 4, 5};
/*double z, q, fg, Mtot;*/
double PhaseTimeArray(double z, double mass_q, double fg, double M_tot);
/*double total_mass[] = {1.0E+7, 1.0E+8, 1.0E+9}; *//* in solar mass */ 
/*double discgas_fraction[] = {0.8, 0.85, 0.9};*//* gas_fraction = M_gas disc/(M_gas disc + M_stellar disc) within 1kpc*/
/*double v_gas[] = {0.2, 0.5, 0.8};*/
double v_gas[] = {0.2, 0.5, 0.8};
double discgas_fraction[] = {0.3, 0.5, 0.9}; /* rc in  kpc */

double total_mass[] = {3.0E+6, 1.0E+7, 1.0E+8};
/*double total_mass[] = {1.0E+6, 3.0E+6, 1.0E+7};*/
/*double total_mass[] = {3.0E+7, 1.0E+8, 3.0E+8};*/
/*double total_mass[] = {2.0E+5, 3.0E+5, 5.0E+5};*/


double rho_gas[] = {1.0E+2, 2.0E+2, 3.0E+2}; /* central gas particle /cm^3 at rc */

double mass_ratio[] = {1.0/4.0, 1.0/8.0, 1.0/10.0};
/*double mass_ratio[] = {1.0/4.0, 1.0/5.0, 1.0/6.0};*/
/*double mass_ratio[] = {1.0/7.0, 1.0/8.0, 1.0/9.0};*/
/*double mass_ratio[] = {1.0/10.0, 1.0/11.0, 1.0/12.0};*/
/*double mass_ratio[] = {1.0/3.0, 1.0/3.5, 1.0/2.5};*/



double type[] = {1.0, 0.0}; /* 1.0 for early type galaxy, 0.0 for late type galaxy */
/* low mass BHBs m≤ 10^7 solar mass */

double G = 2.80210835E+11; /* in N kpc^2/ solar mass ^2 */
double PI = 3.14159265358979323846;

/*  Functions used in Main  */
	
/* Functions for ddr */

	
/* S`T`E`L`L`A`R`  `B`U`L`G`E  */

/* S`T`E`L`L`A`R`  `B`U`L`G`E   D`E`N`S`I`T`Y*/

/* define density profile of the stellar bulge */
double gama = 1.8;/*>0.6 */
double gama_e = 1.8;
double alpha = 4.0;
double alpha_b = 1.8;
double q_b = 0.6;


/* define the enclosed mass of the stellar bulge */
double mass_integral_f1(double s, double r_b, double a_b) 
{
	/*double f = 4.0*PI*s*s*pow((s/r0), (-1*gama)) *pow((1 + pow(s/r0, alpha)),((gama-gama_e)/alpha));*/
	double f = 4.0*PI*s*s*(pow((s/a_b),-1.0*alpha_b)* exp(-1.0*pow((s/r_b),2)));
	return f;
}
	
	
	
double sumintegral_1(double radius, double r_b, double a_b)
{
	int n = 20;
	double lowbound = 0.0;
	double upbound = radius;
	double dx = (double) (upbound - lowbound)/ n;
	double cumsum = 0;
	for (int i = 1; i<n+1; i++)
	{
		double xi = lowbound + i*dx;
		double function_value = mass_integral_f1 (xi, r_b, a_b);
		double rectangle_area = function_value*dx;
		cumsum += rectangle_area;
	}
	return cumsum;
}
	

	
/* This needs to be larger than the maximum radius considered in the system */



/* define the enclosed potential due to the stellar bulge */
double starbulge_potential_f(double s, double r_b, double a_b) 
{
	/*double f = 4.0*PI*s*s*pow((s/r0), (-1*gama)) *pow((1 + pow(s/r0, alpha)),((gama-gama_e)/alpha));*/
	double f = 4.0*PI*s*(pow((s/a_b),-1*alpha_b)* exp(-1*pow((s/r_b),2)));
	return f;
}
	
	
	
double starbulge_potential_integral(double radius, double r_b, double a_b, double R_bulge)
{
	int n = 20;
	double lowbound = radius;
	double upbound = R_bulge;
	double dx = (double) (upbound - lowbound)/ n;
	double cumsum = 0;
	for (int i = 1; i<n+1; i++)
	{
		double xi = lowbound + i*dx;
		double function_value = starbulge_potential_f(xi, r_b, a_b);
		double rectangle_area = function_value*dx;
		cumsum += rectangle_area;
	}
	return cumsum;
}



/* D`F  F`O`R`C`E  */

/* Define dynamical friction by collisionless particles */
/* using GM_BH/velocity dispersion^2 as the maximum impact parameter of the galaxy and using*/
 /* M_BH - velocity dispersion relation for the velocity dispersion */

/* ln (lambda) ~ 10 */

	\
/* S`L`O`W` `M`O`V`I`N`G  S`T`A`R  D`F`  `F`O`R`C`E */
	
/* set up integration for f_DF_integral1 */
double f_df_f1(double s, double velocity_c, double velocity, double p_max, double M2, double vg) 
{
	if (abs(velocity) < pow(2.0, 0.5)* velocity_c){
		/* g is the field star velocity distribution*/
		double g = (tgamma(gama+1)/ tgamma(gama-0.5)) * pow((2*pow(velocity_c,2) - pow(s,2)), (gama-1.5))/ (pow(2,gama) * pow(PI, 1.5) * pow(velocity_c,(2*gama)));
		/*double g = exp(-1*pow(s, 2) / pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),2)) / ( pow(PI, 1.5)* pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),3) ) ;*/
		double f = 4*PI* g * pow(s, 2) * log((p_max/((6.67E-11)* M2*(2E+30))) * (pow(velocity, 2)-pow(s, 2)));
		return f;
		
	}
	else{
		double f = 0.0;
		return f;
	}
}
	
	
	
double sumintegral_10(double velocity_c, double velocity, double p_max,double M2, double vg)
{
	int n = 20;
	double lowbound = 0.0;
	double upbound = abs(velocity);
	double dx = (double) (upbound - lowbound)/ n;
	double cumsum = 0;
	for (int i = 1; i<n; i++)
	{
		double xi = lowbound + i*dx;
		double function_value = f_df_f1 (xi, velocity_c, abs(velocity), p_max, M2, vg);
		double rectangle_area = function_value*dx;
		cumsum += rectangle_area;
	}
	return cumsum;
}

	/*AM2012 DF equations with Maxwellian velocity distribution */
/* set up integration for f_DF_integral1 */
double f_df_f1_maxw(double s, double velocity_c, double velocity, double p_max, double M2, double vg,double sigma_star) 
{
	if (abs(velocity) < pow(2.0, 0.5)* velocity_c){
		/* g is the field star velocity distribution*/
		double g = exp(-1.0*pow(s,2.0)/(2.0*pow(sigma_star,2.0)))/pow(2.0*PI*pow(sigma_star,2.0),1.5);
		/*double g = exp(-1*pow(s, 2) / pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),2)) / ( pow(PI, 1.5)* pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),3) ) ;*/
		double f = 4*PI* g * pow(s, 2) * log((p_max/((6.67E-11)* M2*(2E+30))) * (pow(velocity, 2)-pow(s, 2)));
		return f;
		
	}
	else{
		double f = 0.0;
		return f;
	}
}

double sumintegral_10_maxw(double velocity_c, double velocity, double p_max,double M2, double vg, double sigma_star)
{
	int n = 20;
	double lowbound = 0.0;
	double upbound = abs(velocity);
	double dx = (double) (upbound - lowbound)/ n;
	double cumsum = 0;
	for (int i = 1; i<n; i++)
	{
		double xi = lowbound + i*dx;
		double function_value = f_df_f1_maxw(xi, velocity_c, abs(velocity), p_max, M2, vg, sigma_star);
		double rectangle_area = function_value*dx;
		cumsum += rectangle_area;
	}
	return cumsum;
}
	
/* F`A`S`T  `M`O`V`I`N`G  S`T`A`R   D`F` `F`O`R`C`E */

/* set up integration for f_DF_integral2 */
double f_df_f2(double s, double velocity_c, double velocity, double vg) 
{
	double g = (tgamma(gama+1)/ tgamma(gama-0.5)) * pow((2*pow(velocity_c,2) - pow(s,2)), (gama-1.5)) / (pow(2,gama) * pow(PI, 1.5) * pow(velocity_c,(2*gama)));
	/*double g = exp(-1*pow(s, 2) / pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),2)) / ( pow(PI, 1.5)* pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),3) ) ;*/
	double f = 4*PI* g * pow(s, 2) * (log((s + abs(velocity))/(s - abs(velocity))) - 2* abs(velocity)/s);
	return f;
}
	
	
	
double sumintegral_11(double velocity_c, double velocity, double vg)
{
	int n = 20;
	double lowbound = abs(velocity);
	double upbound = pow(2,0.5)*velocity_c;
	if (lowbound >= upbound) {
		double cumsum = 0;
		return cumsum;
	}
	else {
		double dx = (double) (upbound - lowbound)/ n;
		double cumsum = 0;
		for (int i = 1; i<n; i++)
		{
			double xi = lowbound + i*dx;
			double function_value = f_df_f2 (xi, velocity_c, abs(velocity), vg);
			double rectangle_area = function_value*dx;
			cumsum += rectangle_area;
		}
		return cumsum;
	}
	
}
	
	/* AM2012 DF with Maxwellian velocity distribution : FAST MOVING STAR  */
/* set up integration for f_DF_integral2 */
double f_df_f2_maxw(double s, double velocity_c, double velocity, double vg, double sigma_star) 
{
	double g = exp(-1.0*pow(s,2.0)/(2.0*pow(sigma_star,2.0)))/pow(2.0*PI*pow(sigma_star,2.0),1.5);
	/*double g = exp(-1*pow(s, 2) / pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),2)) / ( pow(PI, 1.5)* pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),3) ) ;*/
	double f = 4*PI* g * pow(s, 2) * (log((s + abs(velocity))/(s - abs(velocity))) - 2* abs(velocity)/s);
	return f;
}
	
	
	
double sumintegral_11_maxw(double velocity_c, double velocity, double vg, double sigma_star)
{
	int n = 20;
	double lowbound = abs(velocity);
	double upbound = pow(2,0.5)*velocity_c;
	if (lowbound >= upbound) {
		double cumsum = 0;
		return cumsum;
	}
	else {
		double dx = (double) (upbound - lowbound)/ n;
		double cumsum = 0;
		for (int i = 1; i<n; i++)
		{
			double xi = lowbound + i*dx;
			double function_value = f_df_f2_maxw(xi, velocity_c, abs(velocity), vg, sigma_star);
			double rectangle_area = function_value*dx;
			cumsum += rectangle_area;
		}
		return cumsum;
	}
	
}


/* BULGE DF */
/* Chandrasekhar DF with Maxwellian velocity distribution */


double f_df_maxw(double s,double sigma_star, double velocity, double M2, double p_max) 
{
	double f = (exp(-1.0*pow(s,2.0)/(2.0*pow(sigma_star,2.0)))/pow(2.0*PI*pow(sigma_star,2.0),1.5))*log((p_max/((6.67E-11)* M2*(2E+30))) * (pow(velocity, 2)-pow(s, 2)));
	return f;
}
	
	
	
double sumintegral_ch_maxw(double velocity_c, double velocity,double sigma_star, double M2, double p_max)
{
	int n = 20;
	double lowbound = 0.0;
	double upbound = abs(velocity);
	if (lowbound >= upbound) {
		double cumsum = 0;
		return cumsum;
	}
	else {
		double dx = (double) (upbound - lowbound)/ n;
		double cumsum = 0;
		for (int i = 1; i<n; i++)
		{
			double xi = lowbound + i*dx;
			double function_value = f_df_maxw(xi, sigma_star, velocity, M2, p_max)*pow(xi,2.0);
			double rectangle_area = function_value*dx;
			cumsum += rectangle_area;
		}
		return cumsum;
	}
	
}
	
/* ODE solver description:                                                    */
/*     The 11 stage 8th order Runge-Kutta method for solving a differential   */
/*    equation y'(x) = f(x,y) with initial condition y = c when x = x0       */
/*    evaluates f(x,y) eleven times per step. For step i+1,                  */
/*     y[i+1] = y[i] + ( 9*k1 + 49*k8 + 64*k9 + 49*k10 + 9*k11 ) / 180        */
/*     where, after setting s21 = sqrt(21),                                   */
double *ddr(double M1, double M2, double rho, double fg, double vg, double t,double r_vector_0, double r_vector_1, double dr_vector_0, double dr_vector_1){
	
	/* t in s */
	/* dr in m / s */
	double r_vector[2] = {r_vector_0, r_vector_1};
	double dr_vector[2] = {dr_vector_0, dr_vector_1};
	
	/*printf("r: %E, phi: %E\n", r_vector_0, r_vector_1);*/
	/*printf("v_r: %E, v_phi: %E\n", dr_vector_0, dr_vector_1);*/
	double theta;
	double nomerge;
	/*theta = r_vector[2];*/
	theta = PI/2.0; /* in radian */
	double r = r_vector[0];
	double phi = r_vector[1];
	/* dr in m/s */
	double dr = pow((pow(dr_vector[0], 2) + pow((dr_vector[1]), 2)), 0.5);
	/* r in (r, phi, theta), phi is the polar angle in the disc plan */
	/*printf ("at radius r (kpc, radian): (%E, %E)\n", r_vector[0], r_vector[1]);*/
	/*printf ("dr (m/s, m/s): (%E, %E)\n", dr_vector[0], dr_vector[1]);*/
	double p_max = (6.67E-11)*M1*(2E+30) / pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),2); /* maximum impact parameter in m*/
	double sigma_star= 2E+5*pow((M1/(1.9E+8)),(1.0/5.1)); /* bulge stellar dispersion in m/s */
	
	/*double vc_10kpc = 3.0E+5*(M1/1.0E+6);*/ /* circular velocity at 10 kpc is 250 km/s */
	/*double vc_10kpc = 2.0E+5;*/
	
	/* build the galaxy */
	
	/*defining the size parameters of the galaxy componants */
	/*double R_bulge=4.0;*/ /* kpc */ 
	double R_bulge=1.0 * log10(M1/1.0E+5); /* kpc */ 
	double R_stellardisc = 10.0* log10(M1/1.0E+5);
	double R_gasdisc = 15.0* log10(M1/1.0E+5);
	double r_b, a_b;
	double Rd, z0, z1;
	double sigma_d;
	double alpha0 = 0.5;
	double alpha1 = 0.5; /* alpha0 + alpha1 = 1 */
	double sigma_g;
	double R0, Rm, zg, Rg;
	
	Rd = 1.0 * log10(M1/1.0E+5); /* in kpc*/
	Rg = 2.0* Rd; /* in kpc*/
	z0 = (0.3/3.0)* log10(M1/1.0E+5); /* kpc */
	z1 = (1.0/3.0)* log10(M1/1.0E+5); /* kpc */
	r_b = (1.9/4.0)*R_bulge;
	a_b = (1.0/4.0)*R_bulge;
	R0 = 8.0* log10(M1/1.0E+5);
	Rm = 0.00001* log10(M1/1.0E+5);
	zg = 0.08;
	
	
	
	
	
	/* G`A`S`  D`I`S`C */
	/*double rc= 1.0* log10(M1/1.0E+5); *//* central: where the particle number density is equaled to rho_gas parameter */
	double rc= 0.0; /* central: where the particle number density is equaled to rho_gas parameter */
	/*double bar_ratio = 0.98;*/
	/*double bar_ratio = 1.0;*/
	/*sigma_g= (rho*(2.0 * zg*3.086E+19)/(1.0E-6/(1.67E-27)))/ (2E+30* exp(-1*(rc * sin(theta)/ Rg) - (Rm / (rc * sin (theta))) - (rc * cos(theta)/ zg) ) );*/ /* M_solar/m^2 */
	/*sigma_g= ((rho/exp(-1.0*rc/Rd))* (1.67E-27*5.027E-31/(1.0E-6))) * 2.0 * zg*(3.086E+19);*/ /* M_solar/m^2 */
	sigma_g= ((rho/exp(-1.0*rc/Rd))* (1.67E-27*5.027E-31/(1.0E-6))) * 2.0* zg*3.086E+19; /* M_solar/m^2 */
	double sigma_g_n = (rho/exp(-1.0*rc/Rd))* 2.0* zg*3.086E+21; /* num/cm^2*/
	double rho_gas_n = (sigma_g_n/(2.0*zg*3.086E+21))*exp(-1.0*r/Rd);
	
	/*double gas_fraction= (M_gas_disc_1kpc * (1.0-bar_ratio))/((bar_ratio)*  M_DM_1kpc);*//* gas fraction within the central 1kpc */
	/*printf("gas fraction is: %E\n", gas_fraction);*/
	/*printf("M_gas_disc_1kpc:%E, sigma_g: %E\n", M_gas_disc_1kpc ,sigma_g);*/
	/*printf("M_DM_1kpc:%E\n", M_DM_1kpc);*/
	double M_gasdisc;
	/* define the enclosed gas disc mass */
	if ((r <= R_gasdisc ) && (r >= 0.0))
	{
		M_gasdisc = 2.0*PI*sigma_g *pow(Rg*3.086E+19, 2) * (1.0- exp(-1.0*r/Rg) * (1.0 + r/Rg)) ;
		/*printf ("mass of the gas disc at r is: %E\n", M_gasdisc);*/
	}
	else if ((r > R_gasdisc ) && (r >= 0.0))
	{
		M_gasdisc = 2.0*PI*sigma_g *pow(Rg*3.086E+19, 2) * (1.0 - exp(-1.0*R_gasdisc/Rg) * (1.0 + R_gasdisc/Rg)) ;
		/*printf ("mass of the gas disc at r is: %E\n", M_gasdisc);*/
	}
	else 
	{	
		M_gasdisc = 0.0;
	}
	
	double M_gas_disc_1kpc;
	/* define the enclosed gas disc mass */
	if ((1.0 <= R_gasdisc ) && (1.0 >= 0.0))
	{
		 M_gas_disc_1kpc = 2.0*PI*sigma_g *pow(Rg*3.086E+19, 2) * (1.0- exp(-1.0*1.0/Rg) * (1.0 + 1.0/Rg)) ;
		/*printf ("mass of the gas disc at r is: %E\n", M_gasdisc);*/
	}
	else if ((1.0 > R_gasdisc ) && (1.0 >= 0.0))
	{
		 M_gas_disc_1kpc = 2.0*PI*sigma_g *pow(Rg*3.086E+19, 2) * (1.0 - exp(-1.0*R_gasdisc/Rg) * (1.0 + R_gasdisc/Rg)) ;
		/*printf ("mass of the gas disc at r is: %E\n", M_gasdisc);*/
	}
	else 
	{	
		 M_gas_disc_1kpc = 0.0;
	}
	
	
	
	/* G`A`S` D`I`S`C` `S`P`S`E`F`I`C `P`O`T`E`N`T`I`A`L*/
	double tiny_r = 1.0E-12*r;
	/*double next_r = r - tiny_r;*/
	double next_r =abs((1.0-1.0E-12)*r);
	double y = abs(r/(2.0*Rg));
	double next_y = abs(next_r/(2.0*Rg));
	/* in spcific potential in J/kg */
	/*printf("y: %E\n", y);*/
	/*printf("next y: %E\n", next_y);*/
	if (y>= 100.0)
	{
		printf("ERROR: not possile merger: (vg: %E, M_total:%E, disc gas fraction: %E, rho gas: %E, mass ratio: %E\n)", vg, M1+M2, fg, rho, M2/(M1+M2));
		nomerge = 1.0;
	}
	else
	{
		nomerge = 0.0;
	}
	double sp_gas_disc= -1.0*PI*6.67E-11*sigma_g*2.0E+30*r*3.086E+19*(boost::math::cyl_bessel_i(0.0, y)*boost::math::cyl_bessel_k(1.0, y)- boost::math::cyl_bessel_i(1.0, y)*boost::math::cyl_bessel_k(0.0, y));
	double next_sp_gas_disc= -1.0*PI*6.67E-11*sigma_g*2.0E+30*next_r*3.086E+19*(boost::math::cyl_bessel_i(0.0, next_y)*boost::math::cyl_bessel_k(1.0, next_y)- boost::math::cyl_bessel_i(1.0, next_y)*boost::math::cyl_bessel_k(0.0, next_y));
	
	double gravity_gas_disc = -1.0*M2*2.0E+30*abs(sp_gas_disc - next_sp_gas_disc)/(tiny_r*3.086E+19);
		
	/* S`T`E`L`L`A`R`  D`I`S`C */
	double M_star_disc_1kpc = M_gas_disc_1kpc * ((1.0-fg)/fg);
	double M_star_disc_1kpc_rest = 2*PI *pow(Rd, 2) * (1.0 - exp(-1.0/Rd) * (1.0 + 1.0/Rd)) ;/* stellar disc mass enclosed within the evaluating r */
	
	/*sigma_d = (((1.0 - gas_fraction) * (bar_ratio)* M_DM_1kpc/ (1.0-bar_ratio)) - M_bulge_1kpc)/ M_star_disc_1kpc;*/ /* solar mass/ kpc^2 */
	
	sigma_d = M_star_disc_1kpc/ M_star_disc_1kpc_rest;
	/*sigma_d = sigma_g*(1.0/pow(3.24E-20, 2)) *((1.0-fg)/fg);*/ /* solar mass / kpc^2 */
	/*printf("check: %E\n", (((1.0 - gas) * (bar_ratio)* M_DM_1kpc/ (1.0-bar_ratio)) - M_bulge_1kpc));*/
	/*printf("M_bulge_1kpc:%E\n", M_bulge_1kpc);*/
	double go;
	if (sigma_d <= 0.0)
	{
		printf("ERROR: not possile galaxy model: (M_total:%E, disc gas fraction: %E, rho gas: %E, mass ratio: %E\n)", M1+M2, fg, rho, M2/(M1+M2));
		go = 0.0;
	}
	else
	{
		go = 1.0;
	}
	
	
	/*printf("M_star_disc_1kpc:%E\n", M_star_disc_1kpc);*/
	
	/* define the enclosed mass of the stellar disc */
	
	double M_stellardisc ;
	if ((r <= R_stellardisc ) && (r >= 0.0))
	{
		M_stellardisc  = 2.0*PI*sigma_d *pow(Rd, 2) * (1.0 - exp(-1.0*r/Rd) * (1.0 + r/Rd)) ;/* stellar disc mass enclosed within the evaluating r */
		/*printf ("mass of the stellar disc at r is: %E\n", M_stellardisc);*/
	}
	else if ((r > R_stellardisc ) && (r >= 0.0))
	{
		M_stellardisc  = 2.0*PI*sigma_d *pow(Rd, 2) * (1.0 - exp(-1.0*R_stellardisc/Rd) * (1.0 + R_stellardisc/Rd)) ;
		/*printf ("mass of the stellar disc at r is: %E\n", M_stellardisc);*/
	}
	else 
	{
		M_stellardisc = 0.0;
	}
	
	
		
	/* S`T`E`L`L`A`R  `D`I`S`C   D`E`N`S`I`T`Y*/
	
	/* define the density profile of stellar disc */
	
	double rho_star_disc = sigma_d * exp(-1.0*r*sin(theta)/Rd) * ((alpha0/(2.0*z0))*exp(-1.0*r*cos(theta)/z0) + (alpha1/(2.0*z1))*exp(-1.0*r*cos(theta)/z1)); /* density of stellar disc at radius r, polar angle theta */

	
	/* S`T`E`L`L`A`R` D`I`S`C` `S`P`S`E`F`I`C `P`O`T`E`N`T`I`A`L*/
	double ys = abs(r/(2.0*Rd));
	double next_ys = abs(next_r/(2.0*Rd));
	/* in spcific potential in J/kg */
	double sp_star_disc= -1.0*PI*6.67E-11*sigma_d*(2.0E+30/pow(3.086E+19,2.0))*r*3.086E+19*(boost::math::cyl_bessel_i(0.0, ys)*boost::math::cyl_bessel_k(1.0, ys)- boost::math::cyl_bessel_i(1.0, ys)*boost::math::cyl_bessel_k(0.0, ys));
	double next_sp_star_disc= -1.0*PI*6.67E-11*sigma_d*(2.0E+30/pow(3.086E+19,2))*next_r*3.086E+19*(boost::math::cyl_bessel_i(0.0, next_ys)*boost::math::cyl_bessel_k(1.0, next_ys)- boost::math::cyl_bessel_i(1.0, next_ys)*boost::math::cyl_bessel_k(0.0, next_ys));
	
	double gravity_star_disc = -1.0*M2*2.0E+30*abs(sp_star_disc - next_sp_star_disc)/(tiny_r*3.086E+19);
	/* S`T`E`L`L`A`R  `B`U`L`G`E */
	
	
	/* Normalization of the stellar bulge */
	double M_bulge_total = sumintegral_1(R_bulge, r_b, a_b);
	double rho0;
	rho0 = 1000.0*(M1) / M_bulge_total;
	
	double M_bulge ;
	if ((r <= R_bulge)&& (r >=0.0))
	{
		M_bulge = rho0* sumintegral_1(r, r_b, a_b); /* stellar bulge mass enclosed within r */
	}
	else if ((r > R_bulge)&& (r >=0.0))
	{
		M_bulge = rho0* sumintegral_1(R_bulge, r_b, a_b);
	}
	
	else 
	{
		M_bulge = 0.0;
	}
	
	double M_bulge_1kpc ;
	if ((1.0 <= R_bulge)&& (1.0 >=0.0))
	{
		M_bulge_1kpc = rho0* sumintegral_1(1.0, r_b, a_b); /* stellar bulge mass enclosed within r */
	}
	else if ((1.0 > R_bulge)&& (1.0 >=0.0))
	{
		M_bulge_1kpc = rho0* sumintegral_1(R_bulge, r_b, a_b);
	}
	
	else 
	{
		M_bulge_1kpc = 0.0;
	}
	
	double fg_1= M_gas_disc_1kpc/(M_gas_disc_1kpc+ M_star_disc_1kpc+ M_bulge_1kpc);
	
	double rho_bulge = rho0 * (pow((r/a_b),-1*alpha_b)* exp(-1*pow((r/r_b),2))); /* stellar bulge density at radius r */
	
	/* B`U`L`G`E `S`P`E`S`I`F`I`C `P`O`T`E`N`T`I`A`L */
	
	/*spesific potential of star bulge in J/kg*/
	double sp_star_bulge_in = -1.0*6.67E-11* (M_bulge*2.0E+30/(r*3.086E+19)); 
	/*double next_sp_star_bulge= -1.0*6.67E-11* (rho0*starbulge_potential_integral(next_r, r_b, a_b))*(2.0E+30/(3.086E+19)); */
	
	/*double gravity_star_bulge = -1.0*M2*2.0E+30*abs(sp_star_bulge - next_sp_star_bulge)/(tiny_r*3.086E+19);*/
	
	double sp_star_bulge_out = -1.0*6.67E-11* (rho0*starbulge_potential_integral(r, r_b, a_b, R_bulge))*(2.0E+30/(3.086E+19));
	
	double sp_star_bulge = sp_star_bulge_in + sp_star_bulge_out; 
	
	/* T`O`T`A`L  D`E`N`S`I`T`Y`  &  M`A`S`S  */
	
	/* total density profile of collisionless particles in the galaxy at the evaluating r and theta */
	/* in solar mass / kpc^3 */
	double rho_collisionless =  rho_bulge  + rho_star_disc;
		
	/* total enclosed mass with in r*/
	/* in solar mass */
	double M_total = M1 + M_stellardisc + M_gasdisc + M_bulge;
	
	
	

	
	/* P`O`T`E`N`T`I`A`L */
	
	/* define potential of the galaxy */


	/* S`M`B`H`  P`O`T`E`N`T`I`A`L */
	
	/* potential of the SMBH at the evaluating r*/
	double potential_SMBH = -1.0*(6.67E-11)*M1*2.0E+30/(r* 3.086E+19);

	/* T`O`T`A`L  `P`O`T`E`N`T`I`A`L */
	/* total potential of the galaxy */
	double p_total = potential_SMBH + sp_star_bulge+ sp_gas_disc+ sp_star_disc; /* at r sepesific potential in J/kg */
	double potential_total = p_total * (M2*2.0E+30); /* total potential in J*/
	
	/* T`O`T`A`L  `E`N`E`R`G`Y */
	
	/*double vc = pow((6.67E-11*M_total*2E+30/(r*3.086E+19)),0.5);*/
	double vc = sqrt(((6.67E-11)*M1*2.0E+30/(r* 3.086E+19)) + ((6.67E-11)*M_bulge*2.0E+30/(r* 3.086E+19)) + (r*3.086E+19)*(abs(sp_star_disc - next_sp_star_disc)/(tiny_r*3.086E+19)) + (r*3.086E+19)*(abs(sp_gas_disc - next_sp_gas_disc)/(tiny_r*3.086E+19)));
	double vd= vg*vc;
	double gravity = -1.0*(M2*2.0E+30 *pow(vc, 2.0)/(r*3.086E+19));
	/*printf("vc: %E\n", vc);*/
	double s_energy = 0.5*pow(dr,2.0) + p_total; /* in J/kg */
	double energy = s_energy*(M2*2.0E+30);
	double keplerian_energy = 0.5*pow(dr,2.0) - ((6.67E-11)*M_total*2.0E+30/(r* 3.086E+19));
	/* T`O`T`A`L  `A`N`G`U`L`A`R   M`O`E`N`T`U`M` */
	
	double s_l = (r*3.086E+19)*dr_vector[1]; /* in m^2/s*kg */
	double angular_momentum = s_l *(M2*2.0E+30);
	
	double parameter_alpha = -1.0*(M2*2.0E+30)*pow(vc, 2.0)*r*3.086E+19;
	double reduced_mass = (M_total*M2)*2.0E+30/(M_total+M2);
	double miu = 6.67E-11*(M_total*2.0E+30);
	double enc = pow((1.0 + 2.0*keplerian_energy*pow(s_l, 2.0)/pow(miu,2.0)),0.5);
	/*double enc = pow((1.0 + 2.0*energy*pow(angular_momentum, 2.0)/(pow(parameter_alpha,2.0)*reduced_mass)),0.5);*/
	double major_a= abs(-1.0*miu/(2.0*keplerian_energy))/3.086E+19; /* semi-major axis in unit of kpc */
	
	
	/*printf("check thinks inside en: %E\n", check);*/
	
	
	
	/* T`O`T`A`L  C`O`L`L`I`S`I`O`N-L`E`S`S` P`A`R`T`I`C`L`E  `D`F  `F`O`R`C`E */
	/*double f_DF_collisionless  = (-4*PI*pow(6.67E-11,2) * pow(M2*(2E+30),2) * rho_collisionless*(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) * (f_DF_integral1 + f_DF_integral2); */
	/* at r, at theta, at dr)*/
	/*printf("rho_collisionless:%E\n", rho_collisionless);*/
	/*printf(" (f_DF_integral1 + f_DF_integral2):%E\n",  (f_DF_integral1 + f_DF_integral2));*/
	/*printf("dr:%E\n", dr);*/
	
	
	/*double facss1 = (4*PI*pow(6.67E-11,2) * pow(M2*(2E+30),2) * rho_collisionless*(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) * (f_DF_integral1 + f_DF_integral2); */
	/*printf("fac1: %E\n", fac1);*/
	/*double f_DF_collisionless;*/
	/*if (dr_vector[1] >= vd){	*/
		/*f_DF_collisionless = -1*abs(facss1);*/
	/*}*/
	/*else{*/
		/*f_DF_collisionless = abs(facss1);*/
	/*}*/
	double xi;
	double vr=dr_vector[0];
	double vp=dr_vector[1];
	double avr=abs(dr_vector[0]);
	double avp=abs(dr_vector[1]);
	double angle = atan(avr/avp);
	if ((vp > 0.0) && (vr>= 0.0))
	{
		xi= 1.5*PI - angle ; 
	}
	else if ((vp > 0.0) && (vr < 0.0))
	{
		xi= 1.5*PI + angle ; 
	}
	
	else if ((vp < 0.0) && (vr >= 0.0))
	{
		xi = 0.5*PI + angle ;
	}
	else if ((vp < 0.0) && (vr < 0.0))
	{
		xi= 0.5*PI - angle ;
		/*printf("phi2: %E\n", phi2);*/
	}
	else if ((vp == 0.0)&& (vr >= 0.0))
	{
		xi = 0.0;
	}
	else if ((vp == 0.0) && (vr < 0.0))
	{
		xi = PI;
	}
	else
	{
		printf("xi is not calculated correctly");
		printf("vp: %E, vr: %E, xi: %E\n", vp,vr,xi);
	}
	
	
	/* D`F`  F`O`R`C`E */
	
	/* D`F` F`O`R`C`E  O`F`  D`M` `B`U`L`G`E */
	
	
	/* S`L`O`W` `M`O`V`I`N`G  S`T`A`R  D`F`  `F`O`R`C`E */
	/* AM2012 DF */
	double f_DF_integral1 = sumintegral_10(vc, dr, p_max, M2, vg); /* at r in meter at theta */
	/* AM2012 DF with Mawellian velocity distribution */
	double f_DF_integral1_maxw= sumintegral_10_maxw(vc, dr, p_max, M2, vg, sigma_star);
	
	/* F`A`S`T  `M`O`V`I`N`G  S`T`A`R   D`F` `F`O`R`C`E */
	/* AM2012 DF */
	double f_DF_integral2 = sumintegral_11(vc, dr, vg); /* at r at theta */
	/* AM2012 DF with Mawellian velocity distribution */
	double f_DF_integral2_maxw= sumintegral_11_maxw(vc, dr, vg, sigma_star);
	
	/* Chandrasekhar DF with Maxwellian velocity distribution */
	double f_DF_integral_chan_maxw = sumintegral_ch_maxw(vc, dr, sigma_star, M2, p_max); /* at r at theta */
	
	
	/* AM2012 DF */
	
	double facss2 = (4*PI*pow(6.67E-11,2) * pow(M2*(2E+30), 2) * rho_bulge *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) * (f_DF_integral1 + f_DF_integral2);
	/*printf("fac1: %E\n", fac1);*/
	double f_DF_star_bulge;
	if (dr_vector[1] > 0.0){	
		f_DF_star_bulge = abs(facss2);
	}
	else if (dr_vector[1] <= 0.0){
		f_DF_star_bulge = abs(facss2);
	}
	else{
		f_DF_star_bulge = 0.0;
	}
	
	double f_DF_star_bulge_r, f_DF_star_bulge_phi;
	
	f_DF_star_bulge_r = f_DF_star_bulge *cos(xi);
	f_DF_star_bulge_phi = f_DF_star_bulge *sin(xi);
	
	/*printf("facss2_1: %E, facss2_2: %E, facss2_3: %E\n", (4*PI*pow(6.67E-11,2) * pow(M2*(2E+30), 2) * rho_bulge *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) , f_DF_integral1 , f_DF_integral2);*/

	/* AM2012 DF with Maxwellian velocity distribution */
	
	double facss_am_maxw = (4*PI*pow(6.67E-11,2) * pow(M2*(2E+30), 2) * rho_bulge *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) * (f_DF_integral1_maxw + f_DF_integral2_maxw);
	/*printf("fac1: %E\n", fac1);*/
	double f_DF_star_bulge_maxw;
	if (dr_vector[1] > 0.0){	
		f_DF_star_bulge_maxw = abs(facss_am_maxw);
	}
	else if (dr_vector[1] <= 0.0){
		f_DF_star_bulge_maxw = abs(facss_am_maxw);
	}
	else{
		f_DF_star_bulge = 0.0;
	}
	
	double f_DF_star_bulge_maxw_r, f_DF_star_bulge_maxw_phi;
	
	f_DF_star_bulge_maxw_r = f_DF_star_bulge_maxw *cos(xi);
	f_DF_star_bulge_maxw_phi = f_DF_star_bulge_maxw *sin(xi);
	
	
	
	/* CHandrasekhar D`F` F`O`R`C`E  O`F`  D`M` `B`U`L`G`E with Maxwellian velocity distribution */
	/*double f_DF_star_bulge = (-4*PI*pow(6.67E-11,2) * pow(M2*(2E+30), 2) * rho_bulge *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) * (f_DF_integral1 + f_DF_integral2);*/
	double facss_ch_maxw;
	if (abs(dr) < pow(2,0.5)*vc){	
		facss_ch_maxw = (16.0*pow(PI,2.0)*pow(6.67E-11,2.0) *pow(M2*(2E+30),2.0) * rho_bulge *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) *f_DF_integral_chan_maxw;
	}
	else{
		facss_ch_maxw = 0.0;
	}
	
	/*printf("fac1: %E\n", fac1);*/
	double f_DF_star_bulge_ch_maxw;
	if (dr_vector[1] > 0.0){	
		f_DF_star_bulge_ch_maxw = abs(facss_ch_maxw);
	}
	else if (dr_vector[1] <= 0.0){
		f_DF_star_bulge_ch_maxw = abs(facss_ch_maxw);
	}
	else{
		f_DF_star_bulge_ch_maxw= 0.0;
	}
	
	double f_DF_star_bulge_ch_maxw_r, f_DF_star_bulge_ch_maxw_phi;
	
	f_DF_star_bulge_ch_maxw_r = f_DF_star_bulge_ch_maxw *cos(xi);
	f_DF_star_bulge_ch_maxw_phi = f_DF_star_bulge_ch_maxw *sin(xi);
	
	
	/* D`F` F`O`R`C`E  O`F` S`T`A`R` D`I`S`K` */
	/*double f_DF_star_disc = (-4*PI*pow(6.67E-11,2) * pow(M2*(2E+30), 2) * rho_star_disc *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)) * (f_DF_integral1 + f_DF_integral2); */
	
	/* AM2012*/
	double f_DF_integral3;
	double eps = 0.01;
	double art_dv = eps*vc;
	if ((dr_vector[1]>= 0.0)&&((abs(dr)-vd )>= art_dv)){	
		f_DF_integral3 = abs(abs(log((p_max/((6.67E-11)* M2*(2E+30))) * (pow(abs(dr), 2)-pow(vd, 2)))));
	}
	else if ((dr_vector[1]>= 0.0)&&((abs(dr) - vd)< -1.0*art_dv)){
		f_DF_integral3 = abs(abs(log((vd + abs(dr))/(vd - abs(dr))) - 2* abs(dr)/vd));
	}
	
	else if ((dr_vector[1]>= 0.0)&&((abs(dr) - vd)>= -1.0*art_dv)&&((abs(dr) - vd)< art_dv)){
		f_DF_integral3 = 0.0;
	}
	
	
	else if ((dr_vector[1]< 0.0)&&((abs(dr) - vd )> art_dv )){
		f_DF_integral3 = abs(abs(log((p_max/((6.67E-11)* M2*(2E+30))) * (pow(abs(dr), 2)-pow(vd, 2)))));
	}
	else if ((dr_vector[1]< 0.0)&&((abs(dr) - vd )< -1.0*art_dv)){
		f_DF_integral3 = abs(abs(log((vd + abs(dr))/(vd - abs(dr))) - 2* abs(dr)/vd));
	}
	else if ((dr_vector[1]< 0.0)&&((abs(dr) - vd) <= art_dv)&&((abs(dr) - vd )>= -1.0*art_dv)){
		f_DF_integral3 = 0.0;
	}
	else{
		f_DF_integral3 =0.0;
	}
	double f_DF_star_disc =(4.0*pow(PI,1.0)*pow(6.67E-11,2) * pow(M2*(2E+30), 2) * rho_star_disc *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2.0)) * (f_DF_integral3);
	
	double f_DF_star_disc_r, f_DF_star_disc_phi;
	
	f_DF_star_disc_r = f_DF_star_disc* cos(xi);
	f_DF_star_disc_phi = f_DF_star_disc*sin(xi);
	


	
	
	/* Chandrasekhar D`F` F`O`R`C`E  O`F` S`T`A`R` D`I`S`K` */
	double f_DF_integral4;
	if (abs(dr)>vd){	
		f_DF_integral4 = abs(abs(log((p_max/((6.67E-11)* M2*(2E+30))) * abs(pow(abs(dr), 2)-pow(vd, 2)))));
	}
	else if (abs(dr)<vd){
		f_DF_integral4 = 0.0;
	}
	
	else {
		f_DF_integral4 =0.0;
	}
	
	double facss_ch_disk = (4.0*pow(PI,1.0)*pow(6.67E-11,2) *pow(M2*(2E+30),2.0) * rho_star_disc *(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2))* (f_DF_integral4);
	/*printf("fac1: %E\n", fac1);*/
	double f_DF_star_disk_ch;
	if (dr_vector[1] > 0.0){	
		f_DF_star_disk_ch = abs(facss_ch_disk);
	}
	else if (dr_vector[1] <= 0.0){
		f_DF_star_disk_ch = abs(facss_ch_disk);
	}
	else{
		f_DF_star_disk_ch= 0.0;
	}
	
	double f_DF_star_disk_ch_r, f_DF_star_disk_ch_phi;
	
	f_DF_star_disk_ch_r = f_DF_star_disk_ch *cos(xi);
	f_DF_star_disk_ch_phi = f_DF_star_disk_ch *sin(xi);
	
	
		
	double mache;
	double mache_new;
	/*mache = pow(dr, 2) / ((PI/2.0)*r * 3.086E+19 * 6.67E-11 * ((sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg)- (Rm / (r * sin(theta)))))));*/
	double cs;
	cs = r * 3.086E+19 *(PI/2.0)* 6.67E-11 * ((sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg))))/(vd);
	double T0 = (pow(cs,2.0)* (1.0*1.67E-27)/ ((5.0/3.0)*1.38E-23));
	/*printf("T0: %E\n", T0);*/
	double T = 1.0E+4+ (pow(cs,2.0)* (1.0*1.67E-27)/ ((5.0/3.0)*1.38E-23));
	/*double mache;*/
	/*double T = 1.0E+4;*/ /*in K */
	cs = pow(((5.0/3.0)*1.38E-23*T/(1.0*1.67E-27)),0.5); /* in m/s*/
	
	
	double cs_c;
	/*double r_c = 0.0001;*/ /* in kpc */
	double r_c = 1.0E-6; /* in kpc */
	cs_c = r_c * 3.086E+19 *(PI/2.0)* 6.67E-11 * ((sigma_g*2E+30 * exp(-1*(r_c * sin(theta)/ Rg))))/(vd);
	double T0_c = (pow(cs_c,2.0)* (1.0*1.67E-27)/ ((5.0/3.0)*1.38E-23));
	/*printf("T0: %E\n", T0);*/
	double T_c = 1.0E+4+ (pow(cs_c,2.0)* (1.0*1.67E-27)/ ((5.0/3.0)*1.38E-23));
	/*double mache;*/
	/*double T = 1.0E+4;*/ /*in K */
	cs_c = pow(((5.0/3.0)*1.38E-23*T_c/(1.0*1.67E-27)),0.5); /* in m/s*/
	double rho_rc = rho*exp(-1*(r_c * sin(theta)/ Rg)); /* in cm^-3 */
	double rho_r = rho*exp(-1*(r * sin(theta)/ Rg)); /* in cm^-3 */
	double TQ =(cs*(vd/(r* 3.086E+19)))/( 6.67E-11* sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg)));
	
	mache = pow((pow((dr_vector[1]-vd),2.0)+pow(dr_vector[0],2.0)),0.5)/cs;
	mache_new = pow((pow((abs(dr_vector[1])+vd),2.0)+pow(dr_vector[0],2.0)),0.5)/cs;
	
	/* Eddington limit of M1 */
	double L1edd = 3.2E+4*M1; /* in L_sun */
	double M1_edd = 2.2E-8*M1; /* SOLAR MASS /YR*/
	double M1_bondi = 8.7E-3*(rho_rc)*pow(T_c/1.0E+4, -1.5)*pow(M1/1.0E+6, 2.0);/* M_sun/yr */
	/* Bondi accretion luminosity of M1 */
	double L_sun = 3.826E+26; /*J/s */
	double L1 = (0.1*M1_bondi*(2.0E+30/(365*24*3600.0))*pow(3.0E+8, 2.0))/L_sun ; /* in L_sun */
	double M1_rate, M2_rate;
	if (L1 < 0.1*L1edd)
	{
		L1=L1;
		M1_rate=M1_bondi;
	}
	else if (L1>=0.1*L1edd)
	{
		L1=0.1*L1edd;
		M1_rate=0.1*M1_edd;
	}
	
	
	
	/* Eddington limit for M2 */
	double L2edd = 3.2E+4*M2; /* in L_sun */
	double M2_edd= 2.2E-8*M2; /* SOLAR MASS /YR*/
	double M2_bondi = 8.7E-3*(rho_r)*pow(T/1.0E+4, -1.5)*pow(M2/1.0E+6, 2.0);/* M_sun/yr */
	double M2_BHL = M2_bondi/pow((1.0+ (pow((dr_vector[1]-vd),2.0)+pow(dr_vector[0],2.0))/pow(cs, 2.0)), 1.5); /* in M_sun/yr */
	double L2 = (0.1*M2_BHL*(2.0E+30/(365*24*3600.0))*pow(3.0E+8, 2.0))/L_sun ; /* in L_sun */
	if (L2 < L2edd)
	{
		L2=L2;
		M2_rate=M2_BHL;
	}
	else if (L2>=L2edd)
	{
		L2=L2edd;
		M2_rate=M2_edd;
	}
	
	
	double I_phi;
	
	double I_phi_new;
		
	
	if ((mache < 1.0) && (mache >= 0))
	{
		I_phi = abs(0.7706 * log((1.0 + mache)/ (1.0004 - 0.9185 * mache)) - 1.4703 * mache) ; 
	}
	else if ((mache < 4.4) && (mache >= 1.0))
	{
		I_phi = log(330.0 * 10.0 * pow((mache- 0.71),5.72) * pow(mache, -9.58));
	}
	else if ((mache >= 4.4) && (mache < 75.0))
	{
		I_phi = abs(log(10.0 / (0.11 * mache + 1.65)));
	}
	else
	{
		I_phi = 0.0;
	}
	
	if ((mache_new< 1.0) && (mache_new >= 0))
	{
		I_phi_new = abs(0.7706 * log((1.0 + mache_new)/ (1.0004 - 0.9185 * mache_new)) - 1.4703 * mache_new ); 
	}
	else if ((mache_new < 4.4) && (mache_new >= 1.0))
	{
		I_phi_new = log(330.0 * 10.0 * pow((mache_new- 0.71),5.72) * pow(mache_new, -9.58));
	}
	else if ((mache_new >= 4.4) && (mache_new < 75.0))
	{
		I_phi_new = abs(log(10.0 / (0.11 * mache_new + 1.65)));
	}
	else
	{
		I_phi_new = 0.0;
	}
	
	/* I_r */
	double I_r;
	
	double I_r_new;
	
	if ((mache< 1.1) && (mache >= 0))
	{
		I_r = pow(mache, 2) * pow(10.0, (3.51 * mache - 4.22)) ; 
	}
	else if ((mache < 4.4) && (mache >= 1.1))
	{
		I_r = 0.5 * log(9.33 * pow(mache, 2) * (pow(mache, 2) - 0.95));
	}
	else
	{
		I_r = 0.3 * pow(mache, 2);
	}
	
	if ((mache_new< 1.1) && (mache_new >= 0))
	{
		I_r_new = pow(mache_new, 2) * pow(10.0, (3.51 * mache_new - 4.22)) ; 
	}
	else if ((mache_new < 4.4) && (mache_new >= 1.1))
	{
		I_r_new = 0.5 * log(9.33 * pow(mache_new, 2) * (pow(mache_new, 2) - 0.95));
	}
	else
	{
		I_r_new = 0.3 * pow(mache_new, 2);
	}
	
	/*printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~I_phi: %E\n", I_phi);*/
	/*printf("0.5vc: %E\n", 0.5*vc);*/
	/*printf("sigma_g: %E\n", sigma_g);*/
	/*printf("wtf: %E\n", (-4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg) - (Rm / (r * sin (theta))) - (r * cos(theta)/ zg) ) / (2.0 * zg*3.086E+19)) * pow(6.67E-11*M2*(2E+30),2) / pow((dr_vector[1]-0.5*vc), 2)));*/
	/*double f_DF_gas_phi = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg)) / (1.6*3.086E+17)) * pow(6.67E-11*M2*(2E+30),2) / pow(dr, 2))* I_phi;*/
	/*double f_DF_gas_r = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg)) / (1.6*3.086E+17)) * pow(6.67E-11*M2*(2E+30),2) / pow(dr, 2))* I_r;*/
	/*printf("mache_r: %E\n", mache_r);*/
	/*printf("mache_phi: %E\n", mache_phi);*/
	/*printf("I_r: %E\n", I_r);*/
	/*printf("I_phi: %E\n", I_phi);*/
	double f_DF_gas_phi;
	double f_DF_gas_r;
	
	double fac1 = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg)) / (2.0 * zg*3.086E+19)) * pow(6.67E-11*M2*(2E+30),2) / (pow((dr_vector[1]-vd),2.0)+pow(dr_vector[0],2.0)))* I_phi;
	/*printf("fac1: %E\n", fac1);*/
	double fac2;
	if ((dr_vector[1]>= 0.0)&&(dr_vector[1] >= vd)){	
		f_DF_gas_phi = -1*fac1;
	}
	else if ((dr_vector[1]>= 0.0)&&(dr_vector[1] < vd)){
		f_DF_gas_phi = fac1;
	}
	else if (dr_vector[1]< 0.0){
		f_DF_gas_phi = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg) ) / (2.0 * zg*3.086E+19)) * pow(6.67E-11*M2*(2E+30),2) / (pow((abs(dr_vector[1])+vd),2.0)+pow(dr_vector[0],2.0)))* I_phi_new;
	}
	else {
		f_DF_gas_phi =  0.0;
	}
	
	if ((dr_vector[0]>= 0.0)&&(dr_vector[1] >= 0.0)){	
		fac2 = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg) ) / (2.0 * zg*3.086E+19)) * pow(6.67E-11*M2*(2E+30),2) / (pow((dr_vector[1]-vd),2.0)+pow(dr_vector[0],2.0)))* I_r;
		/*printf("fac2: %E\n", fac2);*/
		f_DF_gas_r = -1*abs(fac2);
	}
	else if ((dr_vector[0] < 0.0)&&(dr_vector[1] >= 0.0)){
		fac2 = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg) ) / (2.0 * zg*3.086E+19)) * pow(6.67E-11*M2*(2E+30),2) /(pow((dr_vector[1]-vd),2.0)+pow(dr_vector[0],2.0)))* I_r;
		/*printf("fac2: %E\n", fac2);*/
		f_DF_gas_r = -1*abs(fac2);
	}
	
	else if ((dr_vector[0]>= 0.0)&&(dr_vector[1] < 0.0)){	
		fac2 = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg) ) / (2.0 * zg*3.086E+19)) * pow(6.67E-11*M2*(2E+30),2) / (pow((abs(dr_vector[1])+vd),2.0)+pow(dr_vector[0],2.0)))* I_r_new;
		/*printf("fac2: %E\n", fac2);*/
		f_DF_gas_r = -1*abs(fac2);
	}
	else if ((dr_vector[0] < 0.0)&&(dr_vector[1] < 0.0)){
		fac2 = (4 * PI *(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg) ) / (2.0 * zg*3.086E+19)) * pow(6.67E-11*M2*(2E+30),2) /(pow((abs(dr_vector[1])+vd),2.0)+pow(dr_vector[0],2.0)))* I_r_new;
		/*printf("fac2: %E\n", fac2);*/
		f_DF_gas_r = -1*abs(fac2);
	}

	else {
		f_DF_gas_r = 0.0;
	}
	
	
	/* L` C` S`C`A`T`T`E`R`I`N`G */
	double f_LC;
	
	double f_DF_collisionless = f_DF_star_bulge+f_DF_star_disc;
	double f_DF_collisionless_maxw = f_DF_star_bulge_maxw+f_DF_star_disc;
	double f_DF_collisionless_ch_maxw = f_DF_star_bulge_ch_maxw+f_DF_star_disk_ch;
	
	/* define total force magnitude in each direction on the perturber at r at dr */
	double f_phi  = f_DF_gas_phi+f_DF_star_bulge_phi+ f_DF_star_disc_phi;
	double f_phi_maxw  = f_DF_gas_phi+f_DF_star_bulge_maxw_phi+ f_DF_star_disc_phi;
	double f_phi_ch_maxw  = f_DF_gas_phi+f_DF_star_bulge_ch_maxw_phi+ f_DF_star_disk_ch_phi;
	
	
	/*double f_r = -1*(6.67E-11 * M_total*(2E+30)*M2*(2E+30)/ pow(r*3.086E+19,2)) + f_DF_gas_r+ f_DF_star_bulge_r+ f_DF_star_disc_r ;*/

	/*double f_r = -1*(6.67E-11 * M_total*(2E+30)*M2*(2E+30)/ pow(r*3.086E+19,2));*/
	/*double f_r = gravity_star_bulge + gravity_star_disc + gravity_gas_disc;*/
	double f_r = gravity+ f_DF_gas_r+ f_DF_star_bulge_r+ f_DF_star_disc_r;
	double f_r_maxw  = gravity+f_DF_gas_r+f_DF_star_bulge_maxw_r+ f_DF_star_disc_r;
	double f_r_ch_maxw  = gravity+f_DF_gas_r+f_DF_star_bulge_ch_maxw_r+ f_DF_star_disk_ch_r;
	
	/*double f_r = -1*(6.67E-11 * M1*(2E+30)*M2*(2E+30)/ pow(r*3.086E+19,2)) + f_DF_gas_r+ f_DF_star_bulge_r+ f_DF_star_disc_r ;*/
	/*printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~r: %E\n", r);*/
	/*printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~f_phi: %E\n", f_phi);*/
	/*printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~f_r: %E\n", f_r);*/
	/*printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~f_DF_gas_phi: %E\n", f_DF_gas_phi);*/
	/*printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~f_DF_gas_r: %E\n", f_DF_gas_r);*/
	/*printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~f_DF_collisionless : %E\n", f_DF_collisionless );*/
	/*printf("(phi: f_DF_collisionless:%E, f_DF_gas_phi: %E)\n", f_DF_collisionless, f_DF_gas_phi);*/
	/*printf("(r: Gravity: %E, f_DF_gas_r: %E)\n", -1*(6.67E-11 * M_total*(2E+30)*M2*(2E+30)/ pow(r*3.086E+19,2)), f_DF_gas_r);*/

	double a_phi = f_phi / (M2*(2E+30));
	double a_r = f_r / (M2*(2E+30));
	double a_x = a_r*cos(phi) - a_phi * sin(phi);
	double a_y = a_r* sin(phi) + a_phi * cos(phi);
	
	double a_phi_maxw = f_phi_maxw / (M2*(2E+30));
	double a_r_maxw = f_r_maxw / (M2*(2E+30));
	double a_x_maxw = a_r_maxw*cos(phi) - a_phi_maxw * sin(phi);
	double a_y_maxw = a_r_maxw* sin(phi) + a_phi_maxw * cos(phi);
	
	double a_phi_ch_maxw = f_phi_ch_maxw / (M2*(2E+30));
	double a_r_ch_maxw = f_r_ch_maxw / (M2*(2E+30));
	double a_x_ch_maxw = a_r_ch_maxw*cos(phi) - a_phi_ch_maxw * sin(phi);
	double a_y_ch_maxw = a_r_ch_maxw* sin(phi) + a_phi_ch_maxw * cos(phi);
	
	/*printf ("[ax: %E, ay: %E]\n", a_x, a_y);*/

	/*double fly_go;*/
	/*static double fly[13];*/
	/*if (r > 20.0)*/
	/*{*/
		/*printf("ERROR: not possile merger: (M_total:%E, rho_gas: %E, disc gas fraction: %E,relative gas disc velocity: %E, mass ratio: %E)\n", M1+M2, rho, fg, vg, M2/(M1+M2));*/
		/*printf("(phi: f_DF_collisionless:%E, f_DF_gas_phi: %E)\n", f_DF_collisionless, f_DF_gas_phi);*/
		/*printf("(r: Gravity: %E, f_DF_gas_r: %E)\n", -1*(6.67E-11 * M_total*(2E+30)*M2*(2E+30)/ pow(r*3.086E+19,2)), f_DF_gas_r);*/
		/*printf("(f_DF_collisionless: f_DF_collisionless_1: %E, f_DF_collisionless_2: %E)\n", (-4*PI*pow(6.67E-11,2) * pow(M2*(2E+30),2) * rho_collisionless*(2E+30/(pow(3.086E+19, 3))) / pow(dr, 2)),  (f_DF_integral1 + f_DF_integral2));*/
		/*printf("（ dr: %E, vc: %E, p_max: %E, f_DF_integral1: %E, f_DF_integral2: %E)\n", dr, vc, p_max,  f_DF_integral1,  f_DF_integral2);*/
		/*fly[0] = 1.0;*/
		/*fly[1] = 1.0;*/
		/*fly[2] = 1.0;*/
		/*fly[3] = 1.0;*/
		/*fly[4] = 1.0;*/
		/*fly[5] = 1.0;*/
		/*fly[6] = 1.0;*/
		/*fly[7] = 1.0;*/
		/*fly[8] = 1.0;*/
		/*fly[9] = 1.0;*/
		/*fly[10] = 1.0;*/
		/*fly[11] = 1.0;*/
		/*fly[12] = 1.0;*/
		/*return (fly);*/
		/*exit(0);*/
		/*fly_go = 1.0;*/
		
	/*}*/
	/*else*/
	/*{*/
		/*fly_go = 0.0;*/
	/*}*/
	/*printf("ax: %E, ay: %E\n", a_x, a_y);*/
	/*printf("L1: %E, L2: %E\n", L1, L2);*/
	static double output_a[47];
	output_a[0] = a_x_maxw;
	output_a[1] = a_y_maxw;
	output_a[2] = go;
	output_a[3] =  cs;
	output_a[4] =  vc;
	output_a[5] = mache;
	output_a[6] = I_r/pow(mache, 2);
	output_a[7]  = I_phi/pow(mache, 2);
	output_a[8]  =(sigma_g*2E+30 * exp(-1*(r * sin(theta)/ Rg) ) / (2.0 * zg*3.086E+19))*(1.0E-6/(1.67E-27)); /* number density in particle/cm^3*/
	output_a[9]  = 2.0*PI*pow(r*3.086E+19, 1.5)/pow(6.67E-11*M_total*2E+30,0.5); /* in s */
	output_a[10] = TQ;
	output_a[11] = M_gasdisc/(M_stellardisc + M_gasdisc + M_bulge);
	output_a[12] = nomerge;
	output_a[13] = pow((pow((dr_vector[1]-vd),2.0)+pow(dr_vector[0],2.0)),0.5);
	output_a[14] = 1.0;
	output_a[15] = f_DF_gas_r;
	output_a[16] = -1*(6.67E-11 * M_total*(2E+30)*M2*(2E+30)/ pow(r*3.086E+19,2));
	output_a[17] = (dr_vector[1]-vd);
	output_a[18] = f_DF_star_bulge_r;
	output_a[19] = f_DF_star_disc_r;
	output_a[20] = f_DF_gas_phi;
	output_a[21] = 0.0;
	output_a[22] = rho_bulge;
	output_a[23] = rho_star_disc;
	output_a[24] =M_total;
	output_a[25] =T;
	output_a[26] =M_bulge;
	output_a[27] =M_stellardisc;
	output_a[28] =M_gasdisc;
	output_a[29] =f_DF_star_bulge_phi;
	output_a[30] =f_DF_star_disc_phi;
	output_a[31] =s_energy;
	output_a[32] =s_l;
	output_a[33] =enc;
	output_a[34] =major_a;
	/*output_a[34] =0.0;*/
	output_a[35] =miu;
	output_a[36] =fg_1;
	output_a[37] =f_DF_star_bulge;
	output_a[38] =f_DF_star_disc;
	output_a[39] =pow(pow(f_DF_gas_phi,2.0)+pow(f_DF_gas_r,2.0), 0.5);
	output_a[40] =f_DF_gas_r;
	output_a[41] =L1;
	output_a[42] =L2;
	output_a[43] =M2_rate; /* solar mass/yr*/
	output_a[44] =M2_edd; /* solar mass/yr*/
	output_a[45] =M1_rate; /* solar mass/yr*/
	output_a[46] =M1_edd; /* solar mass/yr*/
	return output_a;
	

}





double sqrt21 = 4.58257569495584000680;
static const double c1 = 1.0 / 2.0;    
static const double c2 = (7.0 + sqrt21 ) / 14.0;
static const double c3 = (7.0 - sqrt21 ) / 14.0;

static const double a21 =  1.0 / 2.0;
static const double a31 =  1.0 / 4.0;
static const double a32 =  1.0 / 4.0;
static const double a41 =  1.0 / 7.0;
static const double a42 = -(7.0 + 3.0 * sqrt21) / 98.0;
static const double a43 =  (21.0 + 5.0 * sqrt21) / 49.0;
static const double a51 =  (11.0 + sqrt21) / 84.0;
static const double a53 =  (18.0 + 4.0 * sqrt21) / 63.0;
static const double a54 =  (21.0 - sqrt21) / 252.0;
static const double a61 =  (5.0 + sqrt21) / 48.0;
static const double a63 =  (9.0 + sqrt21) / 36.0;
static const double a64 =  (-231.0 + 14.0 * sqrt21) / 360.0;
static const double a65 =  (63.0 - 7.0 * sqrt21) / 80.0;
static const double a71 =  (10.0 - sqrt21) / 42.0;
static const double a73 =  (-432.0 + 92.0 * sqrt21) / 315.0;
static const double a74 =  (633.0 - 145.0 * sqrt21) / 90.0;
static const double a75 =  (-504.0 + 115.0 * sqrt21) / 70.0;
static const double a76 =  (63.0 - 13.0 * sqrt21) / 35.0;
static const double a81 =  1.0 / 14.0;
static const double a85 =  (14.0 - 3.0 * sqrt21) / 126.0;
static const double a86 =  (13.0 - 3.0 * sqrt21) / 63.0;
static const double a87 =  1.0 / 9.0;
static const double a91 =  1.0 / 32.0;
static const double a95 =  (91.0 - 21.0 * sqrt21) / 576.0;
static const double a96 =  11.0 / 72.0;
static const double a97 = -(385.0 + 75.0 * sqrt21) / 1152.0;
static const double a98 =  (63.0 + 13.0 * sqrt21) / 128.0;
static const double a10_1 =  1.0 / 14.0;
static const double a10_5 =  1.0 / 9.0;
static const double a10_6 = -(733.0 + 147.0 * sqrt21) / 2205.0;
static const double a10_7 =  (515.0 + 111.0 * sqrt21) / 504.0;
static const double a10_8 = -(51.0 + 11.0 * sqrt21) / 56.0;
static const double a10_9 =  (132.0 + 28.0 * sqrt21) / 245.0;
static const double a11_5 = (-42.0 + 7.0 * sqrt21) / 18.0;
static const double a11_6 = (-18.0 + 28.0 * sqrt21) / 45.0;
static const double a11_7 = -(273.0 + 53.0 * sqrt21) / 72.0;
static const double a11_8 =  (301.0 + 53.0 * sqrt21) / 72.0;
static const double a11_9 =  (28.0 - 28.0 * sqrt21) / 45.0;
static const double a11_10 = (49.0 - 7.0 * sqrt21) / 18.0;

static const double  b1  = 9.0 / 180.0;
static const double  b8  = 49.0 / 180.0;
static const double  b9  = 64.0 / 180.0;
double tcos = 4.41504E+17; /*hubble time 14 billion yrs in s */
/* use Runge_Kutta_Verner method to solve the equation of motion of the decay */
/* r, dr in (r, phi)*/
double *solve_eom(double M1, double M2, double gas, double fg, double vg, double t0, double r0_r, double r0_phi, double dr0_r, double dr0_phi, double rf_r, double rf_phi, int knd, int jnd, int nnd, int ind, int vnd, double st, double af, double ai,double lf) 
{
	double e_mean=1000.0;
	 /* r0, dr0 , rf are all 3D vector (r, phi, theta) */
	 /* the multiplication between h and ddr output will be dot product not plain multiply */
	double k1_x, k2_x, k3_x, k4_x, k5_x, k6_x, k7_x, k8_x, k9_x, k10_x, k11_x;
	double l1_x, l2_x, l3_x, l4_x, l5_x, l6_x, l7_x, l8_x, l9_x, l10_x, l11_x;
	double k1_y, k2_y, k3_y, k4_y, k5_y, k6_y, k7_y, k8_y, k9_y, k10_y, k11_y;
	double l1_y, l2_y, l3_y, l4_y, l5_y, l6_y, l7_y, l8_y, l9_y, l10_y, l11_y;
	/*printf("v0_r:%E, v0_phi:%E\n", dr0_r, dr0_phi);*/
	double r0[2] = {r0_r, r0_phi};
	double rf[2] = {rf_r, rf_phi};
	double dr0[2] = {dr0_r, dr0_phi};
	double a_mean[2] = {100.0, 100.0 };
	double en;
	double sE;
	double ma;
	double sL = r0[0]*3.086E+19*dr0[1];
	double h;
	double tau;
	double fg_1;
	double major_a = ai;
	double c1h, c2h, c3h;
	double sign_keeper=0.0;
	double sign_change;
	double en_kep;
	double e0;
	double a0;
	st = st*(3600.0*24*365);
	double M2_new = M2;
	double d11=0.9; /*kpc*/
	double d10=0.8; /*kpc*/
	double d9=0.7;
	double d8=0.6;
	double d7=0.5;
	double d6=0.4;
	double d5=0.3;
	double d4=0.2;
	double d3=0.1;
	double d2=0.05;
	double d1=0.01;
	double tau_d1=0.0;
	double tau_d2=0.0;
	double tau_d3=0.0;
	double tau_d4=0.0;
	double tau_d5=0.0;
	double tau_d6=0.0;
	double tau_d7=0.0;
	double tau_d8=0.0;
	double tau_d9=0.0;
	double tau_d10=0.0;
	
	double M2_ratio1= 0.1;
	double M2_ratio2= 0.2;
	double M2_ratio3= 0.3;
	double M2_ratio4= 0.4;
	double M2_ratio5= 0.5;
	double M2_ratio6= 0.6;
	double M2_ratio7= 0.7;
	double M2_ratio8= 0.8;
	double M2_ratio9= 0.9;
	double tau_M2_dot_ratio0=0.0;
	double tau_M2_dot_ratio1=0.0;
	double tau_M2_dot_ratio2=0.0;
	double tau_M2_dot_ratio3=0.0;
	double tau_M2_dot_ratio4=0.0;
	double tau_M2_dot_ratio5=0.0;
	double tau_M2_dot_ratio6=0.0;
	double tau_M2_dot_ratio7=0.0;
	double tau_M2_dot_ratio8=0.0;
	double tau_M2_dot_ratio9=0.0;
	
	double L_ratio1=-3.0;
	double L_ratio2= -2.0;
	double L_ratio3= -1.5;
	double L_ratio4= -1.0;
	double L_ratio5= -0.5;
	double L_ratio6= 0.0;
	double L_ratio7= 0.5;
	double L_ratio8= 1.0;
	double L_ratio9= 2.0;
	double tau_L0=0.0;
	double tau_L1=0.0;
	double tau_L2=0.0;
	double tau_L3=0.0;
	double tau_L4=0.0;
	double tau_L5=0.0;
	double tau_L6=0.0;
	double tau_L7=0.0;
	double tau_L8=0.0;
	double tau_L9=0.0;
	
	double ecc1=0.1;
	double ecc2=0.2;
	double ecc3=0.3;
	double ecc4=0.4;
	double ecc5=0.5;
	double ecc6=0.6;
	double ecc7=0.7;
	double ecc8=0.8;
	double ecc9=0.9;
	
	double tau_e0=0.0;
	double tau_e1=0.0;
	double tau_e2=0.0;
	double tau_e3=0.0;
	double tau_e4=0.0;
	double tau_e5=0.0;
	double tau_e6=0.0;
	double tau_e7=0.0;
	double tau_e8=0.0;
	double tau_e9=0.0;
	
	double prob[10][10][10]={0.0};
	double prob_dl[10][10]={0.0};
	double prob_de[10][10]={0.0};
	double prob_le[10][10]={0.0};
	
	double d123[11]={0.01,0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9};
	/*double l123[11]={-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0};*/
	double l123[11]={-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0};
	double e123[11]={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	
	double M2_dot_ratio_ave=0.0;
	double M1_edd, M1_bondi;
	/*char filename[20];*/
	/*sprintf(filename, "lnr0.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *f6 = fopen(filename, "a");*/
	/*char filename3[16];*/
	/*sprintf(filename3, "Ind0.%d.%d.%d.%d.%d.txt", vnd, nnd,knd, jnd, ind);*/
	/*FILE *f7 = fopen(filename3, "a");*/
	/*char filename4[22];*/
	/*sprintf(filename4, "rho_gasnd0.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *f8 = fopen(filename4, "a");*/
	/*char filename10[15];*/
	/*sprintf(filename10, "tnd0.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *f10 = fopen(filename10, "a");*/
	/*char filename12[20];*/
	/*sprintf(filename12, "el2.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *f12 = fopen(filename12, "a");*/
	/*char filenamema[20];*/
	/*sprintf(filenamema, "mnd0.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *fma = fopen(filenamema, "a");*/
	/*char filenamea[20];*/
	/*sprintf(filenamea, "and0.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *fa = fopen(filenamea, "a");*/
	double fill = 0.0;
	double orb = 0.0;
	double pre_period= 3.15E+21;/* 1E+14 yr in s */
	double current_period=0.0;
	double tw_DF_star_bulge=0.0;
	double tw_DF_star_disc=0.0;
	double tw_DF_gas_disc=0.0;
	double pw_DF_star_bulge=0.0;
	double pw_DF_star_disc=0.0;
	double pw_DF_gas_disc=0.0;
	
	double tw_DF_star_bulge_phi=0.0;
	double tw_DF_star_disc_phi=0.0;
	double tw_DF_gas_disc_phi=0.0;
	double pw_DF_star_bulge_phi=0.0;
	double pw_DF_star_disc_phi=0.0;
	double pw_DF_gas_disc_phi=0.0;
	
	double tw_DF_star_bulge_r=0.0;
	double tw_DF_star_disc_r=0.0;
	double tw_DF_gas_disc_r=0.0;
	double pw_DF_star_bulge_r=0.0;
	double pw_DF_star_disc_r=0.0;
	double pw_DF_gas_disc_r=0.0;
	double flyout_error = 0.0;
	double L1;
	double L2;
	double M2_dot;
	double M2_edd;
	double E1_ave= 0.0;
	double E2_ave= 0.0;
	double M2_delta = 0.0;
	double M1_delta = 0.0;
	double detect = 0.0;
	while ((r0[0]>= af) && (st < tcos)) {
	/*while (((abs(a_mean[0]+a_mean[1])/2.0)>= af) && (st < tcos)) {*/
	
	/*while ((sL>= lf) && (st < tcos)) {*/
		double phi = r0[1];
		double v_circ = dr0[1];
		
		double *pointer;
		pointer= ddr(M1, M2, gas, fg, vg, st, r0[0], r0[1], dr0[0], dr0[1]);
		if (*(pointer+12) == 1.0){
			printf ("ERROR: This combination is not going to merger: (total mass: %E, gas density: %E, disc gas fraction: %E, gas velocity: %E, mass ration :%E\n)", M1+M2, gas, fg, vg, M2/(M1+M2));
			flyout_error = 1.0;
			break;
		}
		h=0.005*(*(pointer+9));
		
		/*period=*(pointer+9);*//* in s */
		c1h = c1 * h;
		c2h = c2 * h;
		c3h = c3 * h;
		k1_x = h * (dr0[0]* cos(phi) -dr0[1] * sin(phi));
		/* ddr() output output_a [] = {a_r, a_phi} in N/kg */
		l1_x = h * (*pointer);
		
		k1_y = h * (dr0[0] * sin(phi) +dr0[1]*cos(phi));
		/* ddr() output output_a [] = {a_r, a_phi} in N/kg */
		l1_y = h * (*(pointer+1));
		pointer = NULL;
		
		k2_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + a21 * l1_x);
		double x2 = (r0[0]*3.086E+19* cos(phi) + a21 * k1_x);
		double y2 = (r0[0]*3.086E+19* sin(phi) + a21 * k1_y);
		double vx2 = ((dr0[0]* cos(phi)-dr0[1] * sin(phi)) + a21 * l1_x);
		double vy2 = ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + a21 * l1_y);
		double r0_2r [2];
		double dr0_2r [2];
		double phi2;
	
		if ((x2 > 0.0) && (y2 >= 0.0))
		{
			phi2 = atan(y2/x2) ; 
		}
		else if ((x2 < 0.0) && (y2 >= 0.0))
		{
			phi2 = PI - atan(y2/abs(x2));
		}
		else if ((x2 < 0.0) && (y2 < 0.0))
		{
			phi2 = PI + atan(abs(y2)/ abs(x2));
		}
		else if ((x2 > 0.0) && (y2 < 0.0))
		{
			phi2 = 2.0*PI - atan(abs(y2)/ abs(x2));
			/*printf("phi2: %E\n", phi2);*/
		}
		else if ((x2 == 0.0) && (y2 > 0.0))
		{
			phi2 = PI/2.0;
		}
		else if ((x2 == 0.0) && (y2 < 0.0))
		{
			phi2 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x2: %E, y2: %E, phi2: %E\n", x2,y2, phi2);
		}
		/*printf("x2: %E, y2: %E, phi2: %E\n", x2,y2, phi2);*/
		r0_2r[0] = pow((pow(x2,2)+pow(y2,2)),0.5)/3.086E+19;
		r0_2r[1] = phi2;
		dr0_2r[0] = vx2 * cos(phi2) + vy2 * sin(phi2);
		dr0_2r[1] = -1*vx2 *sin(phi2) + vy2 * cos(phi2);
		/*printf("r: %E, phi: %E, v_r: %E, v_phi:%E\n",r0_2r[0], r0_2r[1], dr0_2r[0], dr0_2r[1]);*/
		double *pointer2;
		pointer2= ddr(M1, M2, gas, fg, vg, st + c1h, r0_2r[0],r0_2r[1], dr0_2r[0], dr0_2r[1]);
		l2_x = h * (*pointer2);
		k2_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi))  + a21 * l1_y);
		l2_y = h * (*(pointer2+1));
		pointer2 = NULL;
		/*printf("l2_x: %E, l2_y: %E, k2_y: %E\n", l2_x, l2_y, k2_y);*/
		k3_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a31 * l1_x + a32 * l2_x) );
		double x3 = (r0[0]*3.086E+19* cos(phi)) + ( a31 * k1_x + a32 * k2_x );
		double y3 = (r0[0]*3.086E+19* sin(phi)) + ( a31 * k1_y + a32 * k2_y);
		double vx3 = (dr0[0]* cos(phi)-dr0[1] * sin(phi))+ ( a31 * l1_x + a32 * l2_x); 
		double vy3 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a31 * l1_y + a32 * l2_y);
		double r0_3r [2];
		double dr0_3r [2];
		double phi3;
		if ((x3 > 0.0) && (y3 >= 0.0))
		{
			phi3 = atan(y3/x3) ; 
		}
		else if ((x3 < 0.0) && (y3 >= 0.0))
		{
			phi3 = PI - atan(y3/abs(x3));
		}
		else if ((x3 < 0.0) && (y3 < 0.0))
		{
			phi3 = PI + atan(abs(y3)/ abs(x3));
		}
		else if ((x3 > 0.0) && (y3 < 0.0))
		{
			phi3 = 2*PI - atan(abs(y3)/ x3);
		}
		else if ((x3 == 0.0) && (y3 > 0.0))
		{
			phi3 = PI/2.0;
		}
		else if ((x3 == 0.0) && (y3 < 0.0))
		{
			phi3 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x3: %E, y3: %E\n", x3,y3);
		}
		
		r0_3r[0] = pow((pow(x3,2)+pow(y3,2)),0.5)/3.086E+19;
		r0_3r[1] = phi3;
		dr0_3r[0] = vx3 * cos(phi3) + vy3 * sin(phi3);
		dr0_3r[1] = -1*vx3 *sin(phi3) + vy3 * cos(phi3);

		double *pointer3;
		pointer3 = ddr(M1, M2, gas, fg, vg, st + c1h,  r0_3r[0], r0_3r[1],dr0_3r[0],dr0_3r[1]);
		l3_x = h * (*pointer3);
		k3_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a31 * l1_y + a32 * l2_y) );
		l3_y = h * (*(pointer3+1));
		pointer3 = NULL;
	
		k4_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a41 * l1_x + a42 * l2_x + a43 * l3_x) );
		double x4 =(r0[0]*3.086E+19* cos(phi)) +  ( a41 * k1_x + a42 * k2_x + a43 * k3_x); 
		double y4 = (r0[0]*3.086E+19* sin(phi)) +  ( a41 * k1_y + a42 * k2_y + a43 * k3_y);
		double vx4 = (dr0[0]* cos(phi)-dr0[1] * sin(phi)) +  ( a41 * l1_x + a42 * l2_x + a43 * l3_x);
		double vy4 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) +  ( a41 * l1_y + a42 * l2_y + a43 * l3_y);
		double r0_4r [2];
		double dr0_4r [2];
		double phi4;
		if ((x4 > 0.0) && (y4 >= 0.0))
		{
			phi4 = atan(y4/x4) ; 
		}
		else if ((x4 < 0.0) && (y4 >= 0.0))
		{
			phi4 = PI - atan(y4/abs(x4));
		}
		else if ((x4 < 0.0) && (y4 < 0.0))
		{
			phi4 = PI + atan(abs(y4)/ abs(x4));
		}
		else if ((x4 > 0.0) && (y4 < 0.0))
		{
			phi4 = 2*PI - atan(abs(y4)/ x4);
		}
		else if ((x4 == 0.0) && (y4 > 0.0))
		{
			phi4 = PI/2.0;
		}
		else if ((x4 == 0.0) && (y4 < 0.0))
		{
			phi4 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x4: %E, y4: %E\n", x4,y4);
		}
		
		r0_4r[0] = pow((pow(x4,2)+pow(y4,2)),0.5)/3.086E+19;
		r0_4r[1] = phi4;
		dr0_4r[0] = vx4 * cos(phi4) + vy4 * sin(phi4);
		dr0_4r[1] = -1*vx4 *sin(phi4) + vy4 * cos(phi4);
		double *pointer4;
		pointer4 = ddr(M1, M2, gas, fg,vg, st + c2h, r0_4r[0],r0_4r[1], dr0_4r[0],dr0_4r[1]);
		l4_x = h * (*pointer4);
		k4_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a41 * l1_y + a42 * l2_y + a43 * l3_y) );
		l4_y = h * (*(pointer4+1));
		pointer4 = NULL;
		
		k5_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a51 * l1_x + a53 * l3_x + a54 * l4_x) );
		double x5 =(r0[0]*3.086E+19* cos(phi)) + ( a51 * k1_x + a53 * k3_x + a54 * k4_x ); 
		double y5 = (r0[0]*3.086E+19* sin(phi)) + ( a51 * k1_y + a53 * k3_y + a54 * k4_y );
		double vx5 = (dr0[0]* cos(phi)-dr0[1] * sin(phi)) + ( a51 * l1_x + a53 * l3_x + a54 * l4_x ); 
		double vy5 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a51 * l1_y + a53 * l3_y + a54 * l4_y );
		double r0_5r [2];
		double dr0_5r [2];
		double phi5;
		if ((x5 > 0.0) && (y5 >= 0.0))
		{
			phi5 = atan(y5/x5) ; 
		}
		else if ((x5 < 0.0) && (y5 >= 0.0))
		{
			phi5 = PI - atan(y5/abs(x5));
		}
		else if ((x5 < 0.0) && (y5 < 0.0))
		{
			phi5 = PI + atan(abs(y5)/ abs(x5));
		}
		else if ((x5 > 0.0) && (y5 < 0.0))
		{
			phi5 = 2*PI - atan(abs(y5)/ x5);
		}
		else if ((x5 == 0.0) && (y5 > 0.0))
		{
			phi5 = PI/2.0;
		}
		else if ((x5 == 0.0) && (y5 < 0.0))
		{
			phi5 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x5: %E, y5: %E\n", x5,y5);
		}
		
		r0_5r[0] = pow((pow(x5,2)+pow(y5,2)),0.5)/3.086E+19;
		r0_5r[1] = phi5;
		dr0_5r[0] = vx5 * cos(phi5) + vy5* sin(phi5);
		dr0_5r[1] = -1*vx5 *sin(phi5) + vy5 * cos(phi5);
		double *pointer5;
		pointer5 = ddr(M1, M2, gas, fg,vg, st + c2h, r0_5r[0],r0_5r[1], dr0_5r[0], dr0_5r[1]);
		l5_x = h * (*pointer5);
		k5_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a51 * l1_y + a53 * l3_y + a54 * l4_y) );
		l5_y = h * (*(pointer5+1));
		pointer5 = NULL;
		
		k6_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a61 * l1_x + a63 * l3_x + a64 * l4_x + a65 * l5_x ) );
		double x6 =(r0[0]*3.086E+19* cos(phi)) + ( a61 * k1_x + a63 * k3_x + a64 * k4_x + a65 * k5_x ); 
		double y6 = (r0[0]*3.086E+19* sin(phi)) + ( a61 * k1_y + a63 * k3_y + a64 * k4_y + a65 * k5_y );
		double vx6 = (dr0[0]* cos(phi)-dr0[1] * sin(phi)) + ( a61 * l1_x + a63 * l3_x + a64 * l4_x + a65 * l5_x ); 
		double vy6 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a61 * l1_y + a63 * l3_y + a64 * l4_y + a65 * l5_y );
		double r0_6r [2];
		double dr0_6r [2];
		double phi6;
		if ((x6 > 0.0) && (y6 >= 0.0))
		{
			phi6 = atan(y6/x6) ; 
		}
		else if ((x6 < 0.0) && (y6 >= 0.0))
		{
			phi6 = PI - atan(y6/abs(x6));
		}
		else if ((x6 < 0.0) && (y6 < 0.0))
		{
			phi6 = PI + atan(abs(y6)/ abs(x6));
		}
		else if ((x6 > 0.0) && (y6 < 0.0))
		{
			phi6 = 2*PI - atan(abs(y6)/ x6);
		}
		else if ((x6 == 0.0) && (y6 > 0.0))
		{
			phi6 = PI/2.0;
		}	
		else if ((x6 == 0.0) && (y6 < 0.0))
		{
			phi6 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x6: %E, y6: %E\n", x6,y6);
		}	
		
		r0_6r[0] = pow((pow(x6,2)+pow(y6,2)),0.5)/3.086E+19;
		r0_6r[1] = phi6;
		dr0_6r[0] = vx6 * cos(phi6) + vy6* sin(phi6);
		dr0_6r[1] = -1*vx6 *sin(phi6) + vy6 * cos(phi6);
		double *pointer6;
		pointer6 = ddr(M1, M2, gas, fg, vg, st + c1h, r0_6r[0],r0_6r[1] ,dr0_6r[0], dr0_6r[1]);
		l6_x = h * (*pointer6);
		k6_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a61 * l1_y + a63 * l3_y + a64 * l4_y + a65 * l5_y) );
		l6_y = h * (*(pointer6+1));
		pointer6 = NULL;
		
		k7_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a71 * l1_x + a73 * l3_x + a74 * l4_x + a75 * l5_x + a76 * l6_x) );
		double x7 =(r0[0]*3.086E+19* cos(phi)) + ( a71 * k1_x + a73 * k3_x + a74 * k4_x + a75 * k5_x + a76 * k6_x ); 
		double y7 = (r0[0]*3.086E+19* sin(phi))+ ( a71 * k1_y + a73 * k3_y + a74 * k4_y + a75 * k5_y + a76 * k6_y ) ;
		double vx7 = (dr0[0]* cos(phi)-dr0[1] * sin(phi)) + ( a71 * l1_x + a73 * l3_x + a74 * l4_x + a75 * l5_x + a76 * l6_x ); 
		double vy7 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a71 * l1_y + a73 * l3_y + a74 * l4_y + a75 * l5_y + a76 * l6_y);
		double r0_7r [2];
		double dr0_7r [2];
		double phi7;
		if ((x7 > 0.0) && (y7 >= 0.0))
		{
			phi7 = atan(y7/x7) ; 
		}
		else if ((x7 < 0.0) && (y7 >= 0.0))
		{
			phi7 = PI - atan(y7/abs(x7));
		}
		else if ((x7 < 0.0) && (y7 < 0.0))
		{
			phi7 = PI + atan(abs(y7)/ abs(x7));
		}
		else if ((x7 > 0.0) && (y7 < 0.0))
		{
			phi7 = 2*PI - atan(abs(y7)/ x7);
		}
		else if ((x7 == 0.0) && (y7 > 0.0))
		{
			phi7 = PI/2.0;
		}
		else if ((x7 == 0.0) && (y7 < 0.0))
		{
			phi7 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x7: %E, y7: %E\n", x7,y7);
		}
		
		r0_7r[0] = pow((pow(x7,2)+pow(y7,2)),0.5)/3.086E+19;
		r0_7r[1] = phi7;
		dr0_7r[0] = vx7 * cos(phi7) + vy7* sin(phi7);
		dr0_7r[1] = -1*vx7 *sin(phi7) + vy7 * cos(phi7);
		double *pointer7;
		pointer7 = ddr(M1, M2, gas, fg, vg, st + c3h, r0_7r[0] ,r0_7r[1] , dr0_7r[0], dr0_7r[1]);
		l7_x = h * (*pointer7);
		k7_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a71 * l1_y + a73 * l3_y + a74 * l4_y + a75 * l5_y + a76 * l6_y) );
		l7_y = h * (*(pointer7+1));
		pointer7 = NULL;
		
		k8_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a81 * l1_x + a85 * l5_x + a86 * l6_x + a87 * l7_x ) );
		double x8 =(r0[0]*3.086E+19* cos(phi)) + ( a81 * k1_x + a85 * k5_x + a86 * k6_x + a87 * k7_x );
		double y8 = (r0[0]*3.086E+19* sin(phi))+ ( a81 * k1_y + a85 * k5_y + a86 * k6_y + a87 * k7_y );
		double vx8 = (dr0[0]* cos(phi)-dr0[1] * sin(phi)) + ( a81 * l1_x + a85 * l5_x + a86 * l6_x + a87 * l7_x ); 
		double vy8 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a81 * l1_y + a85 * l5_y + a86 * l6_y + a87 * l7_y ) ;
		double r0_8r [2];
		double dr0_8r [2];
		double phi8;
		if ((x8 > 0.0) && (y8 >= 0.0))
		{
			phi8 = atan(y8/x8) ; 
		}
		else if ((x8 < 0.0) && (y8 >= 0.0))
		{
			phi8 = PI - atan(y8/abs(x8));
		}
		else if ((x8 < 0.0) && (y8 < 0.0))
		{
			phi8 = PI + atan(abs(y8)/ abs(x8));
		}
		else if ((x8 > 0.0) && (y8 < 0.0))
		{
			phi8 = 2*PI - atan(abs(y8)/ x8);
		}
		else if ((x8 == 0.0) && (y8 > 0.0))
		{
			phi8 = PI/2.0;
		}
		else if ((x8 == 0.0) && (y8 < 0.0))
		{
			phi8 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x8: %E, y8: %E\n", x8,y8);
		}
		
		r0_8r[0] = pow((pow(x8,2)+pow(y8,2)),0.5)/3.086E+19;
		r0_8r[1] = phi8;
		dr0_8r[0] = vx8 * cos(phi8) + vy8* sin(phi8);
		dr0_8r[1] = -1*vx8 *sin(phi8) + vy8 * cos(phi8);
		double *pointer8;
		pointer8 = ddr(M1, M2, gas, fg, vg, st + c3h, r0_8r[0],r0_8r[1], dr0_8r[0], dr0_8r[1]);
		l8_x = h * (*pointer8);
		k8_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi))+ ( a81 * l1_y + a85 * l5_y + a86 * l6_y + a87 * l7_y ) );
		l8_y = h * (*(pointer8+1));
		pointer8 = NULL;
		
		k9_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a91 * l1_x + a95 * l5_x + a96 * l6_x + a97 * l7_x + a98 * l8_x ) );
		double x9 = (r0[0]*3.086E+19* cos(phi)) + ( a91 * k1_x + a95 * k5_x + a96 * k6_x + a97 * k7_x + a98 * k8_x ); 
		double y9 = (r0[0]*3.086E+19* sin(phi)) + ( a91 * k1_y + a95 * k5_y + a96 * k6_y + a97 * k7_y + a98 * k8_y );
		double vx9 = (dr0[0]* cos(phi)-dr0[1] * sin(phi))+ ( a91 * l1_x + a95 * l5_x + a96 * l6_x + a97 * l7_x + a98 * l8_x ); 
		double vy9 = (dr0[0] * sin(phi) +dr0[1]*cos(phi))+ ( a91 * l1_y + a95 * l5_y + a96 * l6_y + a97 * l7_y + a98 * l8_y );
		double r0_9r [2];
		double dr0_9r [2];
		double phi9;
		if ((x9 > 0.0) && (y9 >= 0.0))
		{
			phi9 = atan(y9/x9) ; 
		}
		else if ((x9 < 0.0) && (y9 >= 0.0))
		{
			phi9 = PI - atan(y9/abs(x9));
		}
		else if ((x9 < 0.0) && (y9 < 0.0))
		{
			phi9 = PI + atan(abs(y9)/ abs(x9));
		}
		else if ((x9 > 0.0) && (y9 < 0.0))
		{
			phi9 = 2*PI - atan(abs(y9)/ x9);
		}
		else if ((x9 == 0.0) && (y9 > 0.0))
		{
			phi9 = PI/2.0;
		}
		else if ((x9 == 0.0) && (y9 < 0.0))
		{
			phi9 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x9: %E, y9: %E\n", x9,y9);
		}
		
		r0_9r[0] = pow((pow(x9,2)+pow(y9,2)),0.5)/3.086E+19;
		r0_9r[1] = phi9;
		dr0_9r[0] = vx9 * cos(phi9) + vy9* sin(phi9);
		dr0_9r[1] = -1*vx9 *sin(phi9) + vy9 * cos(phi9);
		double *pointer9;
		pointer9 = ddr(M1, M2, gas, fg, vg, st + c1h, r0_9r[0],r0_9r[1], dr0_9r[0], dr0_9r[1]);
		l9_x = h * (*pointer9);
		k9_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a91 * l1_y + a95 * l5_y + a96 * l6_y + a97 * l7_y + a98 * l8_y ) );
		l9_y = h * (*(pointer9+1));
		pointer9 = NULL;
		
		k10_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a10_1 * l1_x + a10_5 * l5_x + a10_6 * l6_x + a10_7 * l7_x + a10_8 * l8_x + a10_9 * l9_x ) );
		double x10 = (r0[0]*3.086E+19* cos(phi))+ ( a10_1 * k1_x + a10_5 * k5_x + a10_6 * k6_x + a10_7 * k7_x + a10_8 * k8_x + a10_9 * k9_x );
		double y10 =  (r0[0]*3.086E+19* sin(phi)) + ( a10_1 * k1_y + a10_5 * k5_y + a10_6 * k6_y + a10_7 * k7_y + a10_8 * k8_y + a10_9 * k9_y );
		double vx10 = (dr0[0]* cos(phi)-dr0[1] * sin(phi)) + ( a10_1 * l1_x + a10_5 * l5_x + a10_6 * l6_x + a10_7 * l7_x + a10_8 * l8_x + a10_9 * l9_x ); 
		double vy10 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a10_1 * l1_y + a10_5 * l5_y+ a10_6 * l6_y + a10_7 * l7_y + a10_8 * l8_y + a10_9 * l9_y );
		double r0_10r [2];
		double dr0_10r [2];
		double phi10;
		if ((x10 > 0.0) && (y10 >= 0.0))
		{
			phi10 = atan(y10/x10) ; 
		}
		else if ((x10 < 0.0) && (y10 >= 0.0))
		{
			phi10 = PI - atan(y10/abs(x10));
		}
		else if ((x10 < 0.0) && (y10 < 0.0))
		{
			phi10 = PI + atan(abs(y10)/ abs(x10));
		}
		else if ((x10 > 0.0) && (y10 < 0.0))
		{
			phi10 = 2*PI - atan(abs(y10)/ x10);
		}
		else if ((x10 == 0.0) && (y10 > 0.0))
		{
			phi10 = PI/2.0;
		}
		else if ((x10 == 0.0) && (y10 < 0.0))
		{
			phi10 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x10: %E, y10: %E\n", x10,y10);
		}
		
		r0_10r[0] = pow((pow(x10,2)+pow(y10,2)),0.5)/3.086E+19;
		r0_10r[1] = phi10;
		dr0_10r[0] = vx10 * cos(phi10) + vy10* sin(phi10);
		dr0_10r[1] = -1*vx10 *sin(phi10) + vy10 * cos(phi10);
		double *pointer10;
		pointer10 = ddr(M1, M2, gas, fg, vg, st + c2h, r0_10r[0], r0_10r[1], dr0_10r[0],dr0_10r[1] );
		l10_x = h * (*pointer10);
		k10_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi))  + ( a10_1 * l1_y + a10_5 * l5_y + a10_6 * l6_y + a10_7 * l7_y + a10_8 * l8_y + a10_9 * l9_y ) );
		l10_y = h * (*(pointer10+1));
		pointer10 = NULL;
		
		k11_x = h * ((dr0[0]* cos(phi) -dr0[1] * sin(phi)) + ( a11_5 * l5_x + a11_6 * l6_x + a11_7 * l7_x + a11_8 * l8_x + a11_9 * l9_x + a11_10 * l10_x ));
		double x11 = (r0[0]*3.086E+19* cos(phi)) + ( a11_5 * k5_x + a11_6 * k6_x + a11_7 * k7_x + a11_8 * k8_x + a11_9 * k9_x + a11_10 * k10_x); 
		double y11 = (r0[0]*3.086E+19* sin(phi)) + ( a11_5 * k5_y + a11_6 * k6_y + a11_7 * k7_y + a11_8 * k8_y + a11_9 * k9_y + a11_10 * k10_y );
		double vx11 = (dr0[0]* cos(phi)-dr0[1] * sin(phi)) + ( a11_5 * l5_x + a11_6 * l6_x + a11_7 * l7_x + a11_8 * l8_x + a11_9 * l9_x + a11_10 * l10_x ); 
		double vy11 = (dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a11_5 * l5_y + a11_6 * l6_y + a11_7 * l7_y + a11_8 * l8_y + a11_9 * l9_y + a11_10 * l10_y );
		double r0_11r [2];
		double dr0_11r [2];
		double phi11;
		if ((x11 > 0.0) && (y11 >= 0.0))
		{
			phi11 = atan(y11/x11) ; 
		}
		else if ((x11 < 0.0) && (y11 >= 0.0))
		{
			phi11 = PI - atan(y11/abs(x11));
		}
		else if ((x11 < 0.0) && (y11 < 0.0))
		{
			phi11 = PI + atan(abs(y11)/ abs(x11));
		}	
		else if ((x11 > 0.0) && (y11 < 0.0))
		{
			phi11 = 2*PI - atan(abs(y11)/ x11);
		}	
		else if ((x11 == 0.0) && (y11 > 0.0))
		{
			phi11 = PI/2.0;
		}
		else if ((x11 == 0.0) && (y11 < 0.0))
		{
			phi11 = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
		}
		
		r0_11r[0] = pow((pow(x11,2)+pow(y11,2)),0.5)/3.086E+19;
		r0_11r[1] = phi11;
		dr0_11r[0] = vx11 * cos(phi11) + vy11* sin(phi11);
		dr0_11r[1] = -1*vx11 *sin(phi11) + vy11 * cos(phi11);
		double *pointer11;
		pointer11 = ddr(M1, M2, gas, fg, vg, st, r0_11r[0],r0_11r[1], dr0_11r[0], dr0_11r[1]);
		l11_x = h * (*pointer11);
		k11_y = h * ((dr0[0] * sin(phi) +dr0[1]*cos(phi)) + ( a11_5 * l5_y + a11_6 * l6_y + a11_7 * l7_y + a11_8 * l8_y + a11_9 * l9_y + a11_10 * l10_y ));
		l11_y = h * (*(pointer11+1));
		pointer11 = NULL;
		st = st + h;
		t0 = t0 + h;/* in second*/
		/*printf ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!t is %E yrs\n", t0/(3600.0*24*365));*/
		double x, y, vx, vy;
		x = (r0[0]*3.086E+19* cos(phi)) + (b1 * k1_x + b8 * k8_x + b9 * k9_x + b8 * k10_x + b1 * k11_x);
		vx = (dr0[0]* cos(phi)-dr0[1] * sin(phi))+ (b1 * l1_x + b8 * l8_x + b9 * l9_x + b8 * l10_x + b1 * l11_x);
		y = (r0[0]*3.086E+19* sin(phi) ) + (b1 * k1_y + b8 * k8_y + b9 * k9_y + b8 * k10_y + b1 * k11_y);
		vy = (dr0[0] * sin(phi) +dr0[1]*cos(phi))  + (b1 * l1_y + b8 * l8_y + b9 * l9_y + b8 * l10_y + b1 * l11_y);
		/*convert (x,y) into (r, phi) */
		
		/* first get phi */
		if ((x > 0.0) && (y >= 0.0))
		{
			phi = atan(y/x) ; 
		}
		else if ((x < 0.0) && (y >= 0.0))
		{
			phi = PI - atan(y/abs(x));
		}
		else if ((x < 0.0) && (y < 0.0))
		{
			phi = PI + atan(abs(y)/ abs(x));
		}
		else if ((x > 0.0) && (y < 0.0))
		{
			phi = 2*PI - atan(abs(y)/ x);
		}
		else if ((x == 0.0) && (y > 0.0))
		{
			phi = PI/2.0;
		}
		else if ((x == 0.0) && (y < 0.0))
		{
			phi = 1.5*PI;
		}
		else
		{
			printf("phi is not calculated correctly");
			printf("x: %E, y: %E\n", x,y);
		}
		/*printf("phi is: %E:\n", phi);*/
		r0[0] = pow((pow(x,2)+ pow(y,2)),0.5)/3.086E+19;
		r0[1] = phi;
		dr0[0] = vx * cos(phi) + vy * sin(phi);
		dr0[1] = -1*vx *sin(phi) + vy * cos(phi);
		sign_change = dr0[0]*sign_keeper;
		sign_keeper = dr0[0];
		orb = orb + 1.0; 
		/* this is angular momentum */
		double angular_momentum;
		angular_momentum = r0[0]*3.086E+19 * M2 *2E+30 *dr0[1];
		
		double *pointer101;
		double dr = pow((pow(dr0[0], 2) + pow((dr0[1]), 2)), 0.5);
		
		pointer101 = ddr(M1, M2, gas, fg, vg, st ,r0[0], r0[1], dr0[0], dr0[1]);
		
		/* ENCENTRICITY */
	
		/* specific angular momentum */
		sL = *(pointer101+32);
		/* specific energy */
		sE = *(pointer101+31);
		
		en = *(pointer101+33);
		
		major_a= *(pointer101+34); /* semi-major axis in unit of kpc */
		
		tau = 0.005*(*(pointer101+9));
		
		fg_1 = *(pointer101+36);
		
		ma = *(pointer101+5);
		
		L1 = *(pointer101+41);
		L2= *(pointer101+42);
		M2_dot = *(pointer101+43); /* in solar mass/yr*/
		M2_edd = *(pointer101+44); /* in solar mass/yr*/
		M1_bondi = *(pointer101+45); /* in solar mass/yr*/
		M1_edd = *(pointer101+46); /* in solar mass/yr*/
		
		E1_ave = E1_ave + L1*tau; /* L*s*/
		E2_ave = E2_ave + L2*tau; 
		
		/*M2 = M2_new;*/
		
		M2_dot_ratio_ave=M2_dot_ratio_ave+(M2_dot/(M2_edd))*(tau/(3600*24*365.0)); 
		double ratio_now=(M2_dot/(M2_edd));
		
		if ((ratio_now<M2_ratio1)&&((ratio_now>=0.0))){
			tau_M2_dot_ratio0=tau_M2_dot_ratio0+tau;
		}
		else if ((ratio_now<M2_ratio2)&&((ratio_now>=M2_ratio1))){
			tau_M2_dot_ratio1=tau_M2_dot_ratio1+tau;
		}

		else if ((ratio_now<M2_ratio3)&&((ratio_now>=M2_ratio2))){
			tau_M2_dot_ratio2=tau_M2_dot_ratio2+tau;
		}
		
		else if ((ratio_now<M2_ratio4)&&((ratio_now>=M2_ratio3))){
			tau_M2_dot_ratio3=tau_M2_dot_ratio3+tau;
		}
		
		else if ((ratio_now<M2_ratio5)&&((ratio_now>=M2_ratio4))){
			tau_M2_dot_ratio4=tau_M2_dot_ratio4+tau;
		}
		
		else if ((ratio_now<M2_ratio6)&&((ratio_now>=M2_ratio5))){
			tau_M2_dot_ratio5=tau_M2_dot_ratio5+tau;
		}
		
		else if ((ratio_now<M2_ratio7)&&((ratio_now>=M2_ratio6))){
			tau_M2_dot_ratio6=tau_M2_dot_ratio6+tau;
		}
		
		else if ((ratio_now<M2_ratio8)&&((ratio_now>=M2_ratio7))){
			tau_M2_dot_ratio7=tau_M2_dot_ratio7+tau;
		}
		else if ((ratio_now<M2_ratio9)&&((ratio_now>=M2_ratio8))){
			tau_M2_dot_ratio8=tau_M2_dot_ratio8+tau;
		}
		else if ((ratio_now<=1.0)&&((ratio_now>=M2_ratio9))){
			tau_M2_dot_ratio9=tau_M2_dot_ratio9+tau;
		}
		
		
		
		double dv_circ = abs(v_circ - dr0[1]);
		
		detect = detect + tau*L2*r0[0];
		
		if ((r0[0] < d2)&&(r0[0] >= d1)){
			tau_d1 = tau_d1+tau;
		}
		else if ((r0[0] < d3)&&(r0[0] >= d2)){
			tau_d2 = tau_d2+tau;
		}
		else if ((r0[0] < d4)&&(r0[0] >= d3)){
			tau_d3 = tau_d3+tau;
		}
		else if ((r0[0] < d5)&&(r0[0] >= d4)){
			tau_d4 = tau_d4+tau;
		}
		else if ((r0[0] < d6)&&(r0[0] >= d5)){
			tau_d5 = tau_d5+tau;
		}
		else if ((r0[0] < d7)&&(r0[0] >= d6)){
			tau_d6 = tau_d6+tau;
		}
		else if ((r0[0] < d8)&&(r0[0] >= d7)){
			tau_d7 = tau_d7+tau;
		}
		else if ((r0[0] < d9)&&(r0[0] >= d8)){
			tau_d8 = tau_d8+tau;
		}
		else if ((r0[0] < d10)&&(r0[0] >= d9)){
			tau_d9 = tau_d9+tau;
		}
		else if ((r0[0] <= d11)&&(r0[0] >= d10)){
			tau_d10 = tau_d10+tau;
		}
		
		
		double logl= log10(L2/L1);
		if ((logl<L_ratio1)&&((logl>=-4.0))){
			tau_L0=tau_L0+tau;
		}
		else if ((logl<L_ratio2)&&(logl>=L_ratio1)){
			tau_L1=tau_L1+tau;
		}

		else if ((logl>=L_ratio2)&&(logl<L_ratio3)){
			tau_L2=tau_L2+tau;
		}
		else if ((logl>=L_ratio3)&&(logl<L_ratio4)){
			tau_L3=tau_L3+tau;
		}
		else if ((logl>=L_ratio4)&&(logl<L_ratio5)){
			tau_L4=tau_L4+tau;
		}
		else if ((logl>=L_ratio5)&&(logl<L_ratio6)){
			tau_L5=tau_L5+tau;
		}
		else if ((logl>=L_ratio6)&&(logl<L_ratio7)){
			tau_L6=tau_L6+tau;
		}
		else if ((logl>=L_ratio7)&&(logl<L_ratio8)){
			tau_L7=tau_L7+tau;
		}
		else if ((logl>=L_ratio8)&&(logl<L_ratio9)){
			tau_L8=tau_L8+tau;
		}
		else if ((logl>=L_ratio9)&&(logl< 3.0)){
			tau_L9=tau_L9+tau;
		}
		
		
		/*if (dv_circ <= 1.0E-5*abs(dr0[1]))*/
		/*{*/
			
			/*fprintf(fma, "%E\t%E\n", st/(3600.0*24*365), ma);*/
			/*double a_m = abs(r0[0]);*/
			/*e_mean = 0.0;*/
			/*fprintf(fa, "%E\t%E\n", a_m, ma);*/
			/*if (orb == 1.0)*/
			/*{*/
				/*e0=e_mean;*/
				/*a0=a_m;*/
			/*}*/
				
			/*fprintf(f12, "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", st/(3600.0*24*365), e_mean, tau/(3600.0*24*365), a_m, *(pointer101+29), *(pointer101+30), *(pointer101+20), pre_period, *(pointer101+37), *(pointer101+38), *(pointer101+39), *(pointer101+18), *(pointer101+19),*(pointer101+40));*/ /* time in yr */
			/*pw_DF_star_bulge=pw_DF_star_bulge+ *(pointer101+37)*(tau/(3600.0*24*365));*/
			/*pw_DF_star_disc=pw_DF_star_disc+ *(pointer101+38)*(tau/(3600.0*24*365));*/
			/*pw_DF_gas_disc=pw_DF_gas_disc+ *(pointer101+39)*(tau/(3600.0*24*365));*/
			
			/*pw_DF_star_bulge_phi=pw_DF_star_bulge_phi+ *(pointer101+29)*(tau/(3600.0*24*365));*/
			/*pw_DF_star_disc_phi=pw_DF_star_disc_phi+ *(pointer101+30)*(tau/(3600.0*24*365));*/
			/*pw_DF_gas_disc_phi=pw_DF_gas_disc_phi+ *(pointer101+20)*(tau/(3600.0*24*365));*/
				
			/*pw_DF_star_bulge_r=pw_DF_star_bulge_r+ *(pointer101+18)*(tau/(3600.0*24*365));*/
			/*pw_DF_star_disc_r=pw_DF_star_disc_r+ *(pointer101+19)*(tau/(3600.0*24*365));*/
			/*pw_DF_gas_disc_r=pw_DF_gas_disc_r+ *(pointer101+40)*(tau/(3600.0*24*365));*/
				
			/*pre_period = current_period;*/
			/*current_period=0.0;*/
			/*tw_DF_star_bulge=0.0;*/
			/*tw_DF_star_disc=0.0;*/
			/*tw_DF_gas_disc=0.0;*/
				
			/*tw_DF_star_bulge_phi=0.0;*/
			/*tw_DF_star_disc_phi=0.0;*/
			/*tw_DF_gas_disc_phi=0.0;*/
				
			/*tw_DF_star_bulge_r=0.0;*/
			/*tw_DF_star_disc_r=0.0;*/
			/*tw_DF_gas_disc_r=0.0;*/
		/*}*/
		
		
		if (sign_change< 0.0)
		{
			
			fill=fill+1.0;
			/*fprintf(fma, "%E\t%E\n", st/(3600.0*24*365), ma);*/
			if (pow(-1.0,fill)<=0.0)
			{
				a_mean[0]= abs(r0[0]);
				
			}
			else if (pow(-1.0,fill)>0.0)
			{
				
				a_mean[1]= abs(r0[0]);
				e_mean = abs(a_mean[0]-a_mean[1])/abs(a_mean[0]+a_mean[1]);
				
				/*fprintf(fa, "%E\t%E\n", abs(a_mean[0]+a_mean[1])/2.0, ma);*/
				if (fill==2.0)
				{
					e0=e_mean;
					a0=abs(a_mean[0]+a_mean[1])/2.0;
					/*printf("e0: %E\n", e0);*/
				}
				/*fprintf(f12, "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", st/(3600.0*24*365), e_mean, tau/(3600.0*24*365), abs(a_mean[0]+a_mean[1])/2.0, tw_DF_star_bulge_phi, tw_DF_star_disc_phi, tw_DF_gas_disc_phi, pre_period, tw_DF_star_bulge, tw_DF_star_disc, tw_DF_gas_disc, tw_DF_star_bulge_r, tw_DF_star_disc_r, tw_DF_gas_disc_r, r0[0]);*/ /* time in yr */
				pw_DF_star_bulge=pw_DF_star_bulge+ tw_DF_star_bulge*(pre_period/(3600.0*24*365));
				pw_DF_star_disc=pw_DF_star_disc+ tw_DF_star_disc*(pre_period/(3600.0*24*365));
				pw_DF_gas_disc=pw_DF_gas_disc+ tw_DF_gas_disc*(pre_period/(3600.0*24*365));
				
				pw_DF_star_bulge_phi=pw_DF_star_bulge_phi+ tw_DF_star_bulge_phi*(pre_period/(3600.0*24*365));
				pw_DF_star_disc_phi=pw_DF_star_disc_phi+ tw_DF_star_disc_phi*(pre_period/(3600.0*24*365));
				pw_DF_gas_disc_phi=pw_DF_gas_disc_phi+ tw_DF_gas_disc_phi*(pre_period/(3600.0*24*365));
				
				pw_DF_star_bulge_r=pw_DF_star_bulge_r+ tw_DF_star_bulge_r*(pre_period/(3600.0*24*365));
				pw_DF_star_disc_r=pw_DF_star_disc_r+ tw_DF_star_disc_r*(pre_period/(3600.0*24*365));
				pw_DF_gas_disc_r=pw_DF_gas_disc_r+ tw_DF_gas_disc_r*(pre_period/(3600.0*24*365));
				
				
				pre_period = current_period;
				current_period=0.0;
				tw_DF_star_bulge=0.0;
				tw_DF_star_disc=0.0;
				tw_DF_gas_disc=0.0;
				
				tw_DF_star_bulge_phi=0.0;
				tw_DF_star_disc_phi=0.0;
				tw_DF_gas_disc_phi=0.0;
				
				tw_DF_star_bulge_r=0.0;
				tw_DF_star_disc_r=0.0;
				tw_DF_gas_disc_r=0.0;
				
			}
			
		}
		
		v_circ = dr0[1];
		current_period= current_period+ tau; /* in s*/
		tw_DF_star_bulge_phi= tw_DF_star_bulge_phi+ (*(pointer101+29))*(tau/pre_period);
		tw_DF_star_disc_phi= tw_DF_star_disc_phi+ (*(pointer101+30))*(tau/pre_period);
		tw_DF_gas_disc_phi= tw_DF_gas_disc_phi+ (*(pointer101+20))*(tau/pre_period);
		
		tw_DF_star_bulge_r= tw_DF_star_bulge_r+ (*(pointer101+18))*(tau/pre_period);
		tw_DF_star_disc_r= tw_DF_star_disc_r+ (*(pointer101+19))*(tau/pre_period);
		tw_DF_gas_disc_r= tw_DF_gas_disc_r+ (*(pointer101+40))*(tau/pre_period);
		
		tw_DF_star_bulge= tw_DF_star_bulge+ (*(pointer101+37))*(tau/pre_period);
		tw_DF_star_disc= tw_DF_star_disc+ (*(pointer101+38))*(tau/pre_period);
		tw_DF_gas_disc= tw_DF_gas_disc+ (*(pointer101+39))*(tau/pre_period);
		
		
		/*fprintf(f12, "%E\t%E\t%E\n", st/(3600.0*24*365), e_mean, abs(a_mean[0]+a_mean[1])/2.0);*/
		
		if ((e_mean<ecc1)&&(e_mean>=0.0)){
			tau_e0=tau_e0+tau;
		}
		else if ((e_mean<ecc2)&&(e_mean>=ecc1)){
			tau_e1=tau_e1+tau;
		}
		else if ((e_mean>=ecc2)&&(e_mean<ecc3)){
			tau_e2=tau_e2+tau;
		}
		else if ((e_mean>=ecc3)&&(e_mean<ecc4)){
			tau_e3=tau_e3+tau;
		}
		else if ((e_mean>=ecc4)&&(e_mean<ecc5)){
			tau_e4=tau_e4+tau;
		}
		else if ((e_mean>=ecc5)&&(e_mean<ecc6)){
			tau_e5=tau_e5+tau;
		}
		else if ((e_mean>=ecc6)&&(e_mean<ecc7)){
			tau_e6=tau_e6+tau;
		}
		else if ((e_mean>=ecc7)&&(e_mean<ecc8)){
			tau_e7=tau_e7+tau;
		}
		else if ((e_mean>=ecc8)&&(e_mean<ecc9)){
			tau_e8=tau_e8+tau;
		}
		else if ((e_mean>=ecc9)&&(e_mean<1.0)){
			tau_e9=tau_e9+tau;
		}
		
		for (int i=0; i<10; i++){
			for (int j=0; j<10; j++){
				for (int k=0; k<10; k++){
					if ((r0[0]>=d123[i])&&(r0[0]<d123[i+1])&&(logl>=l123[j])&&(logl<l123[j+1])&&(e_mean>=e123[k])&&(e_mean<e123[k+1])){
						prob[i][j][k]=prob[i][j][k]+tau;
					}
				}
			}
		}
		
		for (int i=0; i<10; i++){
			for (int j=0; j<10; j++){
				if ((r0[0]>=d123[i])&&(r0[0]<=d123[i+1])&&(logl>=l123[j])&&(logl<l123[j+1])){
					prob_dl[i][j]=prob_dl[i][j]+tau;
				}
			}
		}
		
		for (int i=0; i<10; i++){
			for (int j=0; j<10; j++){
				if ((r0[0]>=d123[i])&&(r0[0]<=d123[i+1])&&(e_mean>=e123[j])&&(e_mean<e123[j+1])){
					prob_de[i][j]=prob_de[i][j]+tau;
				}
			}
		}
		
		for (int i=0; i<10; i++){
			for (int j=0; j<10; j++){
				if ((logl>=l123[i])&&(logl<l123[i+1])&&(e_mean>=e123[j])&&(e_mean<e123[j+1])){
					prob_le[i][j]=prob_le[i][j]+tau;
				}
			}
		}
		
		
		/*fprintf(f6, "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", x, y, r0[0], r0[1], dr0[0], dr0[1], *(pointer101+15), *(pointer101+16), *(pointer101+17), *(pointer101+18), *(pointer101+19), *(pointer101+20), st/(3600.0*24*365),*(pointer101+21), *(pointer101+22), *(pointer101+23),*(pointer101+24),*(pointer101+25),*(pointer101+26),*(pointer101+27),*(pointer101+28), e_mean, abs(a_mean[0]+a_mean[1])/2.0 ,*(pointer101+29), *(pointer101+30), sE, sL, *(pointer101+35), *(pointer101+41),*(pointer101+42));*/
		/*fprintf(f10, "%E\t%E\t%E\n", t0/(3600.0*24*365), r0[0], tau);*/
		/*fprintf(f7, "%E\t%E\t%E\t%E\t%E\t%E\n", r0[0], *(pointer101+3), *(pointer101+4), *(pointer101+5), *(pointer101+6), *(pointer101+7));*/
		/*fprintf(f8, "%E\t%E\t%E\t%E\n", r0[0], *(pointer101+8), *(pointer101+10), *(pointer101+11), *(pointer101+13));*/
		
		pointer101 = NULL;
	}
	
	/*char filename[22];*/
	/*sprintf(filename, "pl.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *fprob = fopen(filename, "a");*/
	/*for (int i=0; i<10; i++){*/
			/*for (int j=0; j<10; j++){*/
				/*for (int k=0; k<10; k++){*/
					/*fprintf(fprob, "%E\t%E\t%E\t%E\n", d123[i], l123[j+1], e123[k+1], prob[i][j][k]);*/
				/*}*/
			/*}*/
		/*}*/
		
	
	/*fclose(fprob);*/
	
	/*char filename[29];*/
	/*sprintf(filename, "emax_el.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *fprob = fopen(filename, "a");*/
	/*for (int i=0; i<10; i++){*/
			/*for (int j=0; j<10; j++){*/
				/*double prob_max=prob[i][j][0];*/
				/*int location_aha = 0;*/
				/*for (int k=1; k<10; k++){*/
					/*double nexte =prob[i][j][k];*/
					/*if (nexte > prob_max){*/
						/*prob_max=prob[i][j][k];*/
						/*location_aha = k;*/
					/*}*/
				/*}*/
				/*double e_max=e123[location_aha+1];*/
				/*printf("emax is: %E\n", e_max);*/
				/*fprintf(fprob, "%E\t%E\t%E\n", d123[i], l123[j+1], e_max);*/
			/*}*/
		/*}*/
		
	
	/*fclose(fprob);*/
	
	
	char filenamep1[25];
	sprintf(filenamep1, "dl_l2.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);
	FILE *fprob1 = fopen(filenamep1, "a");
	for (int i=0; i<10; i++){
			for (int j=0; j<10; j++){
				fprintf(fprob1, "%E\t%E\t%E\n", d123[i+1], l123[j+1], prob_dl[i][j]);
			}
		}
	char filenamep2[25];
	
	/*sprintf(filenamep2, "de_l2.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *fprob2 = fopen(filenamep2, "a");*/
	/*for (int i=0; i<10; i++){*/
			/*for (int j=0; j<10; j++){*/
				/*fprintf(fprob2, "%E\t%E\t%E\n", d123[i], e123[j+1], prob_de[i][j]);*/
			/*}*/
		/*}*/
	/*char filenamep3[25];*/
	/*sprintf(filenamep3, "le_l2.%d.%d.%d.%d.%d.txt", vnd, nnd, knd, jnd, ind);*/
	/*FILE *fprob3 = fopen(filenamep3, "a");*/
	/*for (int i=0; i<10; i++){*/
			/*for (int j=0; j<10; j++){*/
				/*fprintf(fprob3, "%E\t%E\t%E\n", l123[i+1], e123[j+1], prob_le[i][j]);*/
			/*}*/
		/*}*/
	
	fclose(fprob1);
	/*fclose(fprob2);*/
	/*fclose(fprob3);*/
	
	double tau_d_array[10] = {tau_d1, tau_d2,tau_d3,tau_d4,tau_d5,tau_d6,tau_d7,tau_d8,tau_d9,tau_d10};
	double tau_l_array[10] = {tau_L0,tau_L1,tau_L2,tau_L3,tau_L4,tau_L5,tau_L6,tau_L7,tau_L8,tau_L9 };
	double tau_e_array[10] = {tau_e0,tau_e1,tau_e2,tau_e3,tau_e4,tau_e5,tau_e6,tau_e7,tau_e8,tau_e9};
	
	double tau_max=tau_d1;
	int location = 0;
	for (int i=1; i<10; i++){
		double next = tau_d_array[i];
		if (next > tau_max){
			tau_max= tau_d_array[i];
			location= i;
			}
		
	}
	double d_ch = d123[location];
	
	double tau_l_max=tau_L0;
	int location_l = 0;
	for (int i=1; i<10; i++){
		double next_l = tau_l_array[i];
		if (next_l > tau_l_max){
			tau_l_max= tau_l_array[i];
			location_l= i;
			}
	}
	double l_ch = l123[location_l+1];
	
	double tau_e_max=tau_e0;
	int location_e = 0;
	for (int i=1; i<10; i++){
		double next_e = tau_e_array[i];
		if (next_e > tau_e_max){
			tau_e_max= tau_e_array[i];
			location_e= i;
			}
	}
	double e_ch = e123[location_e+1];
	
	double prob_dl_ch = prob_dl[location][location_l];
	
	double Lave_1; 
	double Lave_2; 
	double M2_ave;
	if (t0 ==0.0)
	{
		Lave_1 = 0.0;
		Lave_2 = 0.0;
		M2_ave = 0.0;
	}
	else
	{
		Lave_1 = E1_ave/t0;
		Lave_2 = E2_ave/t0;
		M2_ave = M2_delta/(t0/3600*24*365.0); /* solar mass per yr */
	}
	printf("L1_last: %E, L2_last: %E\n", L1, L2);
	printf("L1_ave: %E, L2_ave: %E\n", Lave_1, Lave_2);
	printf("M2_ave (solar mass/yr): %E\n", M2_ave);
	printf("Detectability (J*m): %E\n", detect);
    /*fclose(f6);*/
    /*fclose(f7);*/
    /*fclose(f8);*/
    /*fclose(f10);*/
	/*fclose(f12);*/
	/*fclose(fma);*/
	/*fclose(fa);*/
    static double out[79];
	out[0] = t0/(3600.0*24*365);
	out[1] = r0[0];
	out[2] = r0[1];
	out[3] = dr0[0];
	out[4] = dr0[1];
	/* this is energy */
	out[5] = sE;
	/* this is angular momentum */
	out[6] = sL;
	/*printf ("r0 (%E, %E), dr0 (%E, %E), E_total(J) (%E), L (%E)\n", r0[0], r0[1], dr0[0], dr0[1], out[5], out[6]);*/
	out[7] = st/(3600.0*24.0*365.0);
	out[8] = e_mean;
	out[9] = abs(a_mean[0]+a_mean[1])/2.0;
	out[10] = fg_1;
	out[11] = pw_DF_star_bulge_phi;
	out[12] = pw_DF_star_disc_phi;
	out[13] = pw_DF_gas_disc_phi;
	out[14] = pw_DF_star_bulge;
	out[15] = pw_DF_star_disc;
	out[16] = pw_DF_gas_disc;
	out[17] = pw_DF_star_bulge_r;
	out[18] = pw_DF_star_disc_r;
	out[19] = pw_DF_gas_disc_r;
	out[20] = e0;
	out[21] = a0;
	out[22] = flyout_error;
	out[23] = L1;
	out[24] = L2;
	out[25] = Lave_1;
	out[26] = Lave_2;
	out[27] = M2_ave;
	out[28] = M2;
	out[29] = detect;
	out[30] =tau_d1;
	out[31] =tau_d2;
	out[32] =tau_d3;
	out[33] =tau_d4;
	out[34] =tau_d5;
	out[35] =tau_d6;
	out[36] =tau_d7;
	out[37] =tau_d8;
	out[38] =tau_d9;
	out[39] =tau_d10;
	out[40] =tau_L0;
	out[41] =tau_L1;
	out[42] =tau_L2;
	out[43] =tau_L3;
	out[44] =tau_L4;
	out[45] =tau_L5;
	out[46] =tau_L6;
	out[47] =tau_L7;
	out[48] =tau_L8;
	out[49] =tau_L9;
	out[50] =tau_e0;
	out[51] =tau_e1;
	out[52] =tau_e2;
	out[53] =tau_e3;
	out[54] =tau_e4;
	out[55] =tau_e5;
	out[56] =tau_e6;
	out[57] =tau_e7;
	out[58] =tau_e8;
	out[59] =tau_e9;
	out[60] =(M2_dot_ratio_ave)/(t0/(3600.0*24*365));
	out[61] =tau_M2_dot_ratio0;
	out[62] =tau_M2_dot_ratio1;
	out[63] =tau_M2_dot_ratio2;
	out[64] =tau_M2_dot_ratio3;
	out[65] =tau_M2_dot_ratio4;
	out[66] =tau_M2_dot_ratio5;
	out[67] =tau_M2_dot_ratio6;
	out[68] =tau_M2_dot_ratio7;
	out[69] =tau_M2_dot_ratio8;
	out[70] =tau_M2_dot_ratio9;
	out[71] =M1_bondi;
	out[72] =M1_edd;
	out[73] =M2_delta/(t0/(3600.0*24*365));
	out[74] =M1_delta/(t0/(3600.0*24*365));
	out[75] = d_ch;
	out[76] = l_ch;
	out[77] = e_ch;
	out[78] = prob_dl_ch;
	
	return out;
	
	
}





/* Semi-analytical simulation of SMBH decay */
double *Decay(double z, double q, double gas,  double fg, double Mtot, double vg, int knd, int jnd, int nnd, int ind, int vnd)
{
	
	double M1 = Mtot * (1- q) ; /* in solar mass */
	double M2 = Mtot * (q) ; /* in solar mass */
	/*in s*/
	
	double time[1];
	double flyout_error1[1];
	double L1[1];
	double L2[1];
	double L1_ave[1];
	double L2_ave[1];
	double M2_ave[1];
	double M2_array[2] = {M2, 0.0};
	double detect_array[1] = {0.0};
	double pw_sb_phi[1];
	double pw_sd_phi[1];
	double pw_gd_phi[1];
	double pw_sb_r[1];
	double pw_sd_r[1];
	double pw_gd_r[1];
	double pw_sb[1];
	double pw_sd[1];
	double pw_gd[1];
	double stime[2]={0.0, 0.0};
	double encen[1]={0.0};
	double merger_status[2]={11.0, 10.0};
	double fg_1[1];
	double e0[1]={0.0};
	double lf[1];
	double a0;
	double e0_i;
	double a_mean;
	double tau_d[10];
	double tau_l[10];
	double tau_e[10];
	double tau_m2[10];
	double M2_dot_ra;
	double M1_dot, M1_edd;
	double ave_acc_m2, ave_acc_m1;
	double d_ch, l_ch, e_ch, prob_dl_ch;
	/*double rstart = astart*(1.0+e);*/ /* starts from the apo-center (farthest from the primary BH) */
	/*double rstart = astart*(1.0-e); *//*  starts from the peri-center (nearest from the primary BH) */
	
	/*v_phi in rad/s*/
	double *pointer_haha;
	pointer_haha = ddr(M1, M2, gas, fg, vg, 0.0 ,rstart, 0.0, 0.0, 0.0);
	double vi = *(pointer_haha+4);
	
	double af_array[5]={rstart, a1, a2, a3, a4};
	double *pointer_haha1;
	pointer_haha1 = ddr(M1, M2, gas, fg, vg, 0.0 ,af_array[1], 0.0, 0.0, 0.0);
	double vc1 = *(pointer_haha1+4);
	double lf1 = af_array[1]*3.086E+19*vc1;
	lf[0]=lf1;
	
	
	/*printf("vi: %E\n", vi);*/
	double go = *(pointer_haha+2);
	double go_d = *(pointer_haha+14);
	double fly_go = *(pointer_haha+12);
	pointer_haha = NULL;
	if (go == 0.0)
	{
		static double error[2] = {1.0, 1.0};
		return error;
	}
	
	if (go_d == 0.0)
	{
		static double error_d[2] = {1.0, 1.0};
		return error_d;
	}
	
	if (fly_go == 1.0)
	{
		static double fly_error[2] = {2.0, 2.0};
		return fly_error;
	
	}
	/* for an orbit with enceentricity e, the initial velocity at phi = 0 is: */
	/*double major_0 = astart/(1.0+e); *//* initial semi-major axis in unit of kpc, start from apocenter (farthest point)*/
	/*double vi_e = 0.5*vi *pow((1.0+e), 0.5);*/
	
	/*double vi_e = vi *pow(((1.0-e)/(1.0+e)), 0.5);*/ /* start from the apo-center */
	double vi_e = vi *pow(((1.0+e)), 0.5); /* start from the peri-center */
	double velocity_array[5][2] = {{0.0, vi_e},
									{0.0,0.0}, 
									{0.0,0.0}, 
									{0.0,0.0}, 
									{0.0,0.0}};
	
	
	double radius_array[5][2]= {{rstart, 0.0},
								 {0.0, 0.0}, 
								 {0.0, 0.0}, 
								 {0.0, 0.0}, 
								 {0.0, 0.0}};
	
	double energy_array[2];
	double L_array[2];
	int i;
	int sum=0;
	for (i=1; i<=1; ++i) { 
		sum = sum+1;
		double *pointer_1;
		int j = sum-1;
		/*printf("j: %d\n", j);*/
		
		/*double *pointer_aha_1;*/
		/*pointer_aha_1 = ddr(M1, M2, gas, fg,vg, 0.0 ,radius_array[j][0], 0.0, 0.0, 0.0);*/
		/*double time_stepsize=0.08*(*(pointer_aha_1+9));*/
		
		/*printf("r0, rf: %E, %E\n", radius_array[j][0], radius_array[j+1][0]);*/
		pointer_1= solve_eom(M1, M2_array[j], gas, fg, vg, 0, radius_array[j][0],radius_array[j][1], velocity_array[j][0],velocity_array[j][1], radius_array[j+1][0],radius_array[j+1][1], knd, jnd, nnd, ind, vnd, stime[j], af_array[j+1],af_array[j], lf[j]);
		time[j] = *pointer_1;
		stime[j+1] = *(pointer_1+7);
		flyout_error1[j] = *(pointer_1+22);
		/*printf("time[j]: %E\n", time[j]);*/
		radius_array[j+1][0]= *(pointer_1+1); 
		radius_array[j+1][1] = *(pointer_1+2);
		velocity_array[j+1][0] = *(pointer_1+3);
		velocity_array[j+1][1] = *(pointer_1+4);
		energy_array[j] = *(pointer_1+5);
		L_array[j] = *(pointer_1+6);
		encen[j] = *(pointer_1+8);
		a_mean = *(pointer_1+9);
		fg_1[j] = *(pointer_1+10);
		pw_sb_phi[j]=*(pointer_1+11);
		pw_sd_phi[j]=*(pointer_1+12);
		pw_gd_phi[j]=*(pointer_1+13);
		pw_sb[j]=*(pointer_1+14);
		pw_sd[j]=*(pointer_1+15);
		pw_gd[j]=*(pointer_1+16);
		pw_sb_r[j]=*(pointer_1+17);
		pw_sd_r[j]=*(pointer_1+18);
		pw_gd_r[j]=*(pointer_1+19);
		
		L1[j] = *(pointer_1+23);
		L2[j] = *(pointer_1+24);
		
		L1_ave[j] = *(pointer_1+25);
		L2_ave[j] = *(pointer_1+26);
		M2_ave[j] = *(pointer_1+27);
		e0[j]=*(pointer_1+20);
		M2_array[j+1] = *(pointer_1+28);
		detect_array[j] = *(pointer_1+29);
		
		tau_d[0] = *(pointer_1+30);
		tau_d[1] = *(pointer_1+31);
		tau_d[2] = *(pointer_1+32);
		tau_d[3] = *(pointer_1+33);
		tau_d[4] = *(pointer_1+34);
		tau_d[5] = *(pointer_1+35);
		tau_d[6] = *(pointer_1+36);
		tau_d[7] = *(pointer_1+37);
		tau_d[8] = *(pointer_1+38);
		tau_d[9] = *(pointer_1+39);
		
		tau_l[0] = *(pointer_1+40);
		tau_l[1] = *(pointer_1+41);
		tau_l[2] = *(pointer_1+42);
		tau_l[3] = *(pointer_1+43);
		tau_l[4] = *(pointer_1+44);
		tau_l[5] = *(pointer_1+45);
		tau_l[6] = *(pointer_1+46);
		tau_l[7] = *(pointer_1+47);
		tau_l[8] = *(pointer_1+48);
		tau_l[9] = *(pointer_1+49);
		
		tau_e[0] = *(pointer_1+50);
		tau_e[1] = *(pointer_1+51);
		tau_e[2] = *(pointer_1+52);
		tau_e[3] = *(pointer_1+53);
		tau_e[4] = *(pointer_1+54);
		tau_e[5] = *(pointer_1+55);
		tau_e[6] = *(pointer_1+56);
		tau_e[7] = *(pointer_1+57);
		tau_e[8] = *(pointer_1+58);
		tau_e[9] = *(pointer_1+59);
		M2_dot_ra = *(pointer_1+60);
		
		tau_m2[0] = *(pointer_1+61);
		tau_m2[1] = *(pointer_1+62);
		tau_m2[2] = *(pointer_1+63);
		tau_m2[3] = *(pointer_1+64);
		tau_m2[4] = *(pointer_1+65);
		tau_m2[5] = *(pointer_1+66);
		tau_m2[6] = *(pointer_1+67);
		tau_m2[7] = *(pointer_1+68);
		tau_m2[8] = *(pointer_1+69);
		tau_m2[9] = *(pointer_1+70);
		
		M1_dot = *(pointer_1+71);
		M1_edd = *(pointer_1+72);
		ave_acc_m2= *(pointer_1+73);
		ave_acc_m1= *(pointer_1+74);
		d_ch = *(pointer_1+75);
		l_ch= *(pointer_1+76);
		e_ch= *(pointer_1+77);
		prob_dl_ch= *(pointer_1+78);
		printf("M2_new (solar mass): %E\n", M2_array[j+1]);
		
		if (j==0)
		{
			a0= *(pointer_1+21);
			e0_i=*(pointer_1+20);
		}
        /*printf("L_array[j]:%g\n",L_array[j]);*/
		pointer_1 = NULL;
		if (stime[j+1] >= 1.4E+17){
			printf ("Stalled: This combination is not going to merger within Hubble time in the radius range:(%E-%E),(total mass: %E, gas density: %E, disc gas fraction: %E, mass ration :%E, gas speed: %E\n)", radius_array[j][0], radius_array[j+1][0], M1+M2, gas, fg, M2/(M1+M2), vg);
			merger_status[j+1] = 10.0; /*0 for stalled */
		}
		else {
			/*printf ("Stalled: This combination is not going to merger within Hubble time in the radius range:(%E-%E),(total mass: %E, gas density: %E, disc gas fraction: %E, mass ration :%E, gas speed: %E\n)", radius_array[j][0], radius_array[j+1][0], M1+M2, gas, fg, M2/(M1+M2), vg);*/
			merger_status[j+1] = 11.0; /*1 for continue to merge */
		}
		
		/*pointer_aha_1 = NULL;*/
		/*printf("radius vector is (%E, %E)", radius_array[j][0], radius_array[j][1]);*/
		/*printf("next radius vector is (%E, %E)", radius_array[j+1][0], radius_array[j+1][1]);*/
		/*printf("velocity vector is (%E, %E)", velocity_array[j][0],velocity_array[j][1]);*/
		/*printf("next velocity vector is (%E, %E)", velocity_array[j+1][0],velocity_array[j+1][1]);*/
	}
	static double phase_time[72];
	phase_time[0] = time[0];
	phase_time[1] = encen[0];
	phase_time[2] = fg_1[0];
	phase_time[3] = flyout_error1[0];
	phase_time[4] = pw_sb_phi[0];
	phase_time[5] = pw_sd_phi[0];
	phase_time[6] = pw_gd_phi[0];
	phase_time[7] = pw_sb[0];
	phase_time[8] = pw_sd[0];
	phase_time[9] = pw_gd[0];
	phase_time[10] = pw_sb_r[0];
	phase_time[11] = pw_sd_r[0];
	phase_time[12] = pw_gd_r[0];
	phase_time[13] = L1[0];
	phase_time[14] = L2[0];
	phase_time[15] = L1_ave[0];
	phase_time[16] = L2_ave[0];
	phase_time[17] = M2_ave[0];
	phase_time[18] = M2_array[1];
	phase_time[19] = detect_array[0];
	phase_time[20] = e0_i;
	phase_time[21] = a0;
	phase_time[22] = a_mean;
	phase_time[23] = tau_d[0];
	phase_time[24] = tau_d[1];
	phase_time[25] = tau_d[2];
	phase_time[26] = tau_d[3];
	phase_time[27] = tau_d[4];
	phase_time[28] = tau_d[5];
	phase_time[29] = tau_d[6];
	phase_time[30] = tau_d[7];
	phase_time[31] = tau_d[8];
	phase_time[32] = tau_d[9];
	
	phase_time[33] = tau_l[0];
	phase_time[34] = tau_l[1];
	phase_time[35] = tau_l[2];
	phase_time[36] = tau_l[3];
	phase_time[37] = tau_l[4];
	phase_time[38] = tau_l[5];
	phase_time[39] = tau_l[6];
	phase_time[40] = tau_l[7];
	phase_time[41] = tau_l[8];
	phase_time[42] = tau_l[9];
	
	phase_time[43] = tau_e[0];
	phase_time[44] = tau_e[1];
	phase_time[45] = tau_e[2];
	phase_time[46] = tau_e[3];
	phase_time[47] = tau_e[4];
	phase_time[48] = tau_e[5];
	phase_time[49] = tau_e[6];
	phase_time[50] = tau_e[7];
	phase_time[51] = tau_e[8];
	phase_time[52] = tau_e[9];
	phase_time[53] = M2_dot_ra;
	
	phase_time[54] = tau_m2[0];
	phase_time[55] = tau_m2[1];
	phase_time[56] = tau_m2[2];
	phase_time[57] = tau_m2[3];
	phase_time[58] = tau_m2[4];
	phase_time[59] = tau_m2[5];
	phase_time[60] = tau_m2[6];
	phase_time[61] = tau_m2[7];
	phase_time[62] = tau_m2[8];
	phase_time[63] = tau_m2[9];
	phase_time[64] = M1_dot;
	phase_time[65] = M1_edd;
	phase_time[66] =ave_acc_m2;
	phase_time[67] =ave_acc_m1;
	phase_time[68] = d_ch;
	phase_time[69] = l_ch;
	phase_time[70] = e_ch;
	phase_time[71] = prob_dl_ch;
	
	return phase_time;
	
}

double galaxy_number_density(double s, double z) 
{
	double x = log10(s);
	double f = 0.0035*pow((1.0+z),-2.20)*log(10.0)*pow(pow(10.0, (x-11.16-0.17*z+0.07*pow(z,2.0))),(1.0-1.18-0.082*z))* exp(-1*pow(10.0, (x-11.16-0.17*z+0.07*pow(z,2.0))));
	return f;
}
	
	
double sumintegral_last(double M_min, double M_max, double z)
{
	int n = 20;
	double lowbound = M_min;
	double upbound = M_max;
	double dx = (double) (upbound - lowbound)/ n;
	double cumsum = 0;
	for (int i = 1; i<n+1; i++)
	{
		double xi = lowbound + i*dx;
		double function_value = galaxy_number_density (xi, z);
		double rectangle_area = function_value*log10(dx);
		cumsum += rectangle_area;
	}
	return cumsum;
}


int main() {
	char filename11[14];
	sprintf(filename11, "Lnl36.txt");
	FILE *f11 = fopen(filename11, "a");
	double p1=1.0, p2=1.0, p3=1.0, p4= 1.0, p5=1.0, p6=1.0, p7= 1.0, p8= 1.0, p9= 1.0, p10= 1.0,p11= 1.0,p12= 1.0;
	double encentricity1;
	double merger_result;
	double fg_1_1, fg_1_2, fg_1_3, fg_1_4;
	double pw_sb, pw_sd, pw_gd;
	double pw_sb_phi, pw_sd_phi, pw_gd_phi;
	double pw_sb_r, pw_sd_r, pw_gd_r;
	double e0, a0;
	double M2_dot_ratio_average;
	double flyout=0.0;
	int i;
	int j, k, ni,nii, red = 0;
	int n_q = 3, n_gas = 3, n_fg = 3, n_Mtot = 3, n_vg = 3;
	double z = 1.0, q, fg, Mtot, vg;
	double H0 = 2.2683E-18; /* in 1/s *//*70.0 km/s/Mpc*/
	double OMEGA_Mh2= 0.1364;
	double OMEGA_M= 0.279;
	double OMEGA_LAMBDA= 0.721;
	double OMEGA_bh2= 0.02264;
	double OMEGA_B= 0.0463;
	double OMEGA_RAD= 8.4e-5;
	double H;/* in 1/s */
	H = H0*sqrt(OMEGA_RAD*(pow((1+z),4))+OMEGA_M*pow((1+z),3)+OMEGA_LAMBDA+(1.0-(OMEGA_M+OMEGA_LAMBDA))*pow((1+z),2));
	double dt1;
	double t_array[1];
	double f_dt1;
	double phase_time_1[n_q][1];
	double a_ratio[3] = {1.0, 1.0, 1.0}; /* the fraction of different q mass ratio */
	double n1_ratio[n_q]; /* the fraction of a_ratio*[dt1/(dt1+dt2+dt3+dt4)] for different q */
	double n2_ratio[n_q]; /* the fraction of a_ratio*[dt2/(dt1+dt2+dt3+dt4)] for different q */
	double n3_ratio[n_q]; /* the fraction of a_ratio*[dt3/(dt1+dt2+dt3+dt4)] for different q */
	double n4_ratio[n_q]; /* the fraction of a_ratio*[dt4/(dt1+dt2+dt3+dt4)] for different q */
	double b_ratio[3] = {1.0, 1.0, 1.0}; /* the fraction of different gas fraction */
	int size_gas = n_q*n_gas;
	double n1_ratio_b[n_gas][n_q]; /* the fraction of a_ratio*[dt1/(dt1+dt2+dt3+dt4)] for different q */
	double n2_ratio_b[n_gas][n_q]; /* the fraction of a_ratio*[dt2/(dt1+dt2+dt3+dt4)] for different q */
	double n3_ratio_b[n_gas][n_q]; /* the fraction of a_ratio*[dt3/(dt1+dt2+dt3+dt4)] for different q */
	double n4_ratio_b[n_gas][n_q]; /* the fraction of a_ratio*[dt4/(dt1+dt2+dt3+dt4)] for different q */
	double c_ratio[3] = {1.0, 1.0, 1.0}; /* the fraction of different total mass */
	double n1_ratio_c[n_Mtot][n_gas][n_q];
	double n2_ratio_c[n_Mtot][n_gas][n_q];
	double n3_ratio_c[n_Mtot][n_gas][n_q];
	double n4_ratio_c[n_Mtot][n_gas][n_q];
	double n_ratio_c[n_Mtot][n_gas][n_q];
	
	double d_ratio[3] = {1.0, 1.0, 1.0}; /* the fraction of different total mass */
	double n1_ratio_d[n_fg][n_Mtot][n_gas][n_q];
	double n2_ratio_d[n_fg][n_Mtot][n_gas][n_q];
	double n3_ratio_d[n_fg][n_Mtot][n_gas][n_q];
	double n4_ratio_d[n_fg][n_Mtot][n_gas][n_q];
	double n_ratio_d[n_fg][n_Mtot][n_gas][n_q];
	
	double e_ratio[3] = {1.0, 1.0, 1.0}; /* the fraction of different gas rotatioal velocity at 1.0 kpc */
	double n1_ratio_e[n_vg][n_fg][n_Mtot][n_gas][n_q];
	double n2_ratio_e[n_vg][n_fg][n_Mtot][n_gas][n_q];
	double n3_ratio_e[n_vg][n_fg][n_Mtot][n_gas][n_q];
	double n4_ratio_e[n_vg][n_fg][n_Mtot][n_gas][n_q];
	double n_ratio_e[n_vg][n_fg][n_Mtot][n_gas][n_q];
	
	int size_Mtot = n_vg*n_q*n_gas*n_Mtot*n_fg;
	double n1_ratio_total[size_Mtot]; /* the fraction of a_ratio*[dt1/(dt1+dt2+dt3+dt4)] for different q */
	double n2_ratio_total[size_Mtot]; /* the fraction of a_ratio*[dt2/(dt1+dt2+dt3+dt4)] for different q */
	double n3_ratio_total[size_Mtot]; /* the fraction of a_ratio*[dt3/(dt1+dt2+dt3+dt4)] for different q */
	double n4_ratio_total[size_Mtot]; /* the fraction of a_ratio*[dt4/(dt1+dt2+dt3+dt4)] for different q */
	double total = 0.0;
	double f1_total;
	double f2_total;
	double f3_total;
	double f4_total;
	double flyout_error2;
	/* Reduce duplicated calculation, and define dummy variables here: */
	double L1=0.0;
	double L2=0.0;
	double L3=0.0;
	double L4=0.0;
	double L5=0.0;
	
	double L1_r0=0.0;
	
	double L2_r1=0.0;

	double M2_r1=0.0;
	
	double L1_ave_r0=0.0;
	double L2_ave_r1=0.0;
	double M2_ave_r1=0.0;
	double M2_average;
	double L1_ave_tot = 0.0;
	double L2_ave_tot = 0.0;
	double M2_ave_tot = 0.0;
	double detect_r1=0.0;
	double rstart = 1.0;
	double M2_new_r1;
	int summmmm = 0;
	for(nii=1; nii<=n_vg; ++nii) {
		summmmm = summmmm + 1;
		int vnd = summmmm -1;
		
		int summmm = 0;
		for (ni=1; ni<=n_fg; ++ni) {
			summmm = summmm + 1;
			int nnd = summmm -1;
		
			int summm = 0;
			for (k=1; k<=n_Mtot; ++k) {
				summm = summm + 1;
				int knd = summm -1;
	
		
				int summ = 0;
				for (j=1; j<=n_gas; ++j) {
					summ = summ + 1;
					int jnd = summ -1;
		
			
					int sum = 0;
					for (i=1; i<=n_q; ++i) {
						sum = sum+1;
						int ind = sum-1;
					
						printf("vnd: %d\n", vnd);
						printf("nnd: %d\n", nnd);
						printf("knd: %d\n", knd);
						printf("jnd: %d\n", jnd);
						printf("ind: %d\n", ind);
						double *pointer_2;
				
						pointer_2= Decay(z, mass_ratio[ind], rho_gas[jnd], discgas_fraction[nnd], total_mass[knd], v_gas[vnd], knd, jnd, nnd, ind, vnd);
						flyout_error2 = *(pointer_2+3);
						if (*pointer_2 == 1.0)
						{
							printf ("ERROR: This combination is not a physical galaxy model: (gas velocity: %E, total mass: %E, rho_gas: %E, disc gas fraction: %E, mass ration :%E\n)", v_gas[vnd], total_mass[knd], rho_gas[jnd], discgas_fraction[nnd], mass_ratio[ind]);
							continue;
						}
				
						if (flyout_error2 ==0.0 )
						{
							flyout=0.0;
						}
						else if (flyout_error2 == 1.0)
						{
							printf ("ERROR: This combination is not going to merger: (gas velocity: %E, total mass: %E, rho_gas: %E, disc gas fraction: %E, mass ration :%E\n)",v_gas[vnd], total_mass[knd], rho_gas[jnd], discgas_fraction[nnd], mass_ratio[ind]);
							flyout=1.0;
							printf("flyout:%E", flyout);
						}
						phase_time_1[ind][0] = *pointer_2;
						encentricity1 = *(pointer_2+1);
						printf ("ef: %E\n", encentricity1);
						fg_1_1 = *(pointer_2+2);
						dt1 = phase_time_1[ind][0];
						double meh = dt1;
						/*if (dt4== 0.0)*/
						/*{*/
							/*merger_result = 0.0;*/
						/*}*/
						if (meh >= 1.4E+10)
						{
							merger_result = 0.0;
						}
						else
						{
							merger_result = 1.0;
						}
						
						pw_sb_phi = *(pointer_2+4) /(meh);
						pw_sd_phi = *(pointer_2+5) /(meh);
						pw_gd_phi = *(pointer_2+6) /(meh);
						
						pw_sb = *(pointer_2+7) /(meh);
						pw_sd = *(pointer_2+8) /(meh);
						pw_gd = *(pointer_2+9) /(meh);
						
						pw_sb_r = *(pointer_2+10) /(meh);
						pw_sd_r = *(pointer_2+11) /(meh);
						pw_gd_r = *(pointer_2+12) /(meh);
						
						L1_r0 = *(pointer_2+13) ;
						L2_r1 = *(pointer_2+14);
						
						L1_ave_r0 = *(pointer_2+15);
						L2_ave_r1 = *(pointer_2+16);
						
						M2_ave_r1 = *(pointer_2+17);
						M2_new_r1 = *(pointer_2+18);
						
						detect_r1 = *(pointer_2+19);
						
						e0 = *(pointer_2+20);
						a0 = *(pointer_2+21);
						double a_merger = *(pointer_2+22);
						
						double taud0=*(pointer_2+23);
						double taud1=*(pointer_2+24);
						double taud2=*(pointer_2+25);
						double taud3=*(pointer_2+26);
						double taud4=*(pointer_2+27);
						double taud5=*(pointer_2+28);
						double taud6=*(pointer_2+29);
						double taud7=*(pointer_2+30);
						double taud8=*(pointer_2+31);
						double taud9=*(pointer_2+32);
						
						double taul0=*(pointer_2+33);
						double taul1=*(pointer_2+34);
						double taul2=*(pointer_2+35);
						double taul3=*(pointer_2+36);
						double taul4=*(pointer_2+37);
						double taul5=*(pointer_2+38);
						double taul6=*(pointer_2+39);
						double taul7=*(pointer_2+40);
						double taul8=*(pointer_2+41);
						double taul9=*(pointer_2+42);
						
						double taue0=*(pointer_2+43);
						double taue1=*(pointer_2+44);
						double taue2=*(pointer_2+45);
						double taue3=*(pointer_2+46);
						double taue4=*(pointer_2+47);
						double taue5=*(pointer_2+48);
						double taue6=*(pointer_2+49);
						double taue7=*(pointer_2+50);
						double taue8=*(pointer_2+51);
						double taue9=*(pointer_2+52);
						
						M2_dot_ratio_average=*(pointer_2+53);
						
						double taum0=*(pointer_2+54);
						double taum1=*(pointer_2+55);
						double taum2=*(pointer_2+56);
						double taum3=*(pointer_2+57);
						double taum4=*(pointer_2+58);
						double taum5=*(pointer_2+59);
						double taum6=*(pointer_2+60);
						double taum7=*(pointer_2+61);
						double taum8=*(pointer_2+62);
						double taum9=*(pointer_2+63);
						
						double M1_dot =*(pointer_2+64);
						double M1_edd=*(pointer_2+65);
						double ave_acc_m2=*(pointer_2+66);
						double ave_acc_m1=*(pointer_2+67);
						double d_ch=*(pointer_2+68);
						double l_ch=*(pointer_2+69);
						double e_ch=*(pointer_2+70);
						double prob_dl_ch=*(pointer_2+71);
						
						
						/*printf("e0 written: %E\n", e0);*/
				
						n1_ratio[ind] = (a_ratio[ind]/(a_ratio[0]+a_ratio[1]+a_ratio[2]));
						n1_ratio_b[jnd][ind] = n1_ratio[ind]* (b_ratio[jnd]/(b_ratio[0]+b_ratio[1]+b_ratio[2]));
						n1_ratio_c[knd][jnd][ind] = n1_ratio_b[jnd][ind]* (c_ratio[knd]/(c_ratio[0]+c_ratio[1]+c_ratio[2]));
						n1_ratio_d[nnd][knd][jnd][ind] = n1_ratio_c[knd][jnd][ind]* (d_ratio[nnd]/(d_ratio[0]+d_ratio[1]+d_ratio[2]));
						n1_ratio_e[vnd][nnd][knd][jnd][ind] = n1_ratio_d[nnd][knd][jnd][ind]* (e_ratio[vnd]/(e_ratio[0]+e_ratio[1]+e_ratio[2]));
						
						n_ratio_e[vnd][nnd][knd][jnd][ind] = n1_ratio_e[vnd][nnd][knd][jnd][ind];
						total =total+ n_ratio_e[vnd][nnd][knd][jnd][ind];
						
						L1_ave_tot = L1_ave_tot + L1_ave_r0*n_ratio_e[vnd][nnd][knd][jnd][ind];
						printf("L1_ave_tot is: %E\n", L1_ave_tot);
						L2_ave_tot = L2_ave_tot + L2_ave_r1*n1_ratio_e[vnd][nnd][knd][jnd][ind];
						printf("L2_ave_tot is: %E\n", L2_ave_tot);
						M2_average= M2_ave_r1;
						M2_ave_tot = M2_ave_tot + M2_ave_r1*n1_ratio_e[vnd][nnd][knd][jnd][ind];
						printf("M2_ave_tot is (solar mass/yr): %E\n", M2_ave_tot);
						printf("########################################time: %E \n", dt1);
						printf("######################################## a_ merger: %E\n", a_merger);
						printf("################################################n(1)_ratio_d: %E\n", n1_ratio_e[vnd][nnd][knd][jnd][ind]);
						pointer_2 = NULL;
						fprintf(f11, "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", v_gas[vnd], discgas_fraction[nnd],total_mass[knd],rho_gas[jnd],mass_ratio[ind], n_ratio_e[vnd][nnd][knd][jnd][ind], encentricity1, merger_result, fg_1_1, pw_sb_phi, pw_sd_phi, pw_gd_phi, pw_sb, pw_sd, pw_gd, pw_sb_r, pw_sd_r, pw_gd_r, e0, a0, dt1, flyout_error2, L1_r0, L2_r1, L1_ave_r0, L2_ave_r1, n1_ratio_e[vnd][nnd][knd][jnd][ind], M2_average, M2_new_r1, detect_r1, taud0,taud1,taud2,taud3,taud4,taud5,taud6,taud7,taud8,taud9,taul0,taul1,taul2,taul3,taul4,taul5,taul6,taul7,taul8,taul9,taue0,taue1,taue2,taue3,taue4,taue5,taue6,taue7,taue8,taue9,M2_dot_ratio_average,taum0,taum1,taum2,taum3,taum4,taum5,taum6,taum7,taum8,taum9, M1_dot, M1_edd, ave_acc_m2,ave_acc_m1, d_ch,l_ch,e_ch,prob_dl_ch);
						L1 = L1 + (f_dt1*pow(10.0, (p1*log10(a1)+p2)))*pow(10.0, (p3*log10(mass_ratio[ind])+p4));
					}	
					L2 = L2+ L1*(b_ratio[jnd]/(b_ratio[0]+b_ratio[1]+b_ratio[2]))*pow(10.0, (p5*log10(rho_gas[jnd])+p6));
	
				}	
				L3 = L3+ L2*(c_ratio[knd]/(c_ratio[0]+c_ratio[1]+c_ratio[2]))*pow(10.0, (p7*log10(total_mass[knd])+p8));
			}
			L4 = L4+ L3*(d_ratio[nnd]/(d_ratio[0]+d_ratio[1]+d_ratio[2]))*pow(10.0, (p9*log10(discgas_fraction[nnd])+p10));
		}
		L5 = L5+ L4*(e_ratio[vnd]/(e_ratio[0]+e_ratio[1]+e_ratio[2]))*pow(10.0, (p11*log10(v_gas[vnd])+p12));
	}
	/* Calculating the number density of AGN mergers at this redshift using Hopkins merger rate: */
	/* For mass ratio < 10:1 including both major and minor mergers, above a minimum galaxy stellar mass M_min (equ 5, 9, 10 in Hopkins 2010) */
	/* dN/dt (dNdt) : Number of merger per galaxy (M > M_min) per Gyr */
	double dNdt, A, B, M_min, M_0 = 2*pow(10,10), M_max; /* M_0 in solar mass unit */
	M_min = (pow((2.0E+5*pow(10.0, ((log10(total_mass[0]*(1.0-mass_ratio[0]))-8.13)/(4.02)))),3)/ (10.0*6.67E-11*H))/2.0E+30;
	M_max = (pow((2.0E+5*pow(10.0, ((log10(total_mass[2]*(1.0-mass_ratio[2]))-8.13)/(4.02)))),3)/ (10.0*6.67E-11*H))/2.0E+30;
	double n_galaxy = sumintegral_last(M_min, M_max, z); /* number density of galaxy (M > M_min) */
	double observe_time = 10.0; /* in Gyr */
	double n_merger_total; /* total merger number density in observe time */
	A = 0.04*(1+pow((M_min/M_0), 0.8));
	B = 1.50- 0.25*log10(M_min/M_0);
	dNdt = A*(pow((1+z), B));
	n_merger_total = dNdt*n_galaxy*observe_time/2.0; /* divide by 2 to account for the double counting */
	
	/* L_total of all different M1 and fg considered in the problem and at all separation considered in the problem  in L_solar */
	/* L2_ave can take into account of the accretion of the secondary SMBH*/
	
	double L_total = (L1_ave_tot + L2_ave_tot)*n_merger_total;
	double n_pair_total; /* total AGN pair (with separation r0-r4) number density in observe time at z with fixed fg, Mtot */
	n_pair_total = n_merger_total*n1_ratio_e[0][0][0][0][0];
	printf ("AGN pair number density with q (%E), rho_gas (%E), disc gas fraction (%E), Mtot (%E), v_gas (%E) within separation (%E-%E) is: %E\n", mass_ratio[0], rho_gas[0], discgas_fraction[0], total_mass[0], v_gas[0], rstart, a1, n_pair_total);
	printf("total is: %E\n", total);
	printf("L_total from AGN pairs averaged over all parameters per unit volume is: %E\n", L_total);
	fclose(f11);
}