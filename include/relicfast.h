#ifndef __RELICFAST__
#define __RELICFAST__

	#include "common.h"
	#include "spectra.h"
	#include "lensing.h"

	#define _FALSE_ 0
	#define _TRUE_ 1


    #define Nkpoints 30 //number of p (momentum) points. Log-spaced, 30 is enough.
    #define option_initial_velocity 0 //	0 for sigma_dot/sigma, 1 for T_dot/T (@k_star=2/Mpc). They agree, but it's a good sanity check.

	
	//collapse code options:
	#define do_clustering_tag   _FALSE_   //whether to do clustering of extra species or not.
	#define debug_mode  0 	//whether to print out and save to file some tests. 0 (or _FALSE_) for none, 1 (or _TRUE_) for some, 2 or larger for A LOT.


	//and precision parameters:

	#define precision_scale 10 //from 1 to 10, precision in the scale-dependence of b, through \delta_crit. Recommended 2, above 4 barely any difference.
	#define precision_normalization 4  //from 1 to 10, precision in the normalization of the bias, it also makes the collapse ODE solver use more points, taking longer. Recommended 5, above 10 barely any difference.

	#define boost_initial_conditions 1.0 // to make initial conditions of \delta_short bigger by that amount. Code will complain if it cannot find \delta_short.
	#define delta_short_constant_centroid 0.015 //the centroid of the delta_shorts over which we search.
	#define delta_long_max (delta_short_constant_centroid/10.0) //the maximum delta_long we take for to calculate the bias. If you change z_i check that it still works.
	#define N_delta_long 2 //number of delta_longs for derivative. 2 is sufficient, since it's extremely linear.

	#define precision_clustering 1.0//from 0.5 to 2, makes the BKT clustering integral use more points. 1 yields % level precision for 1-eV neutrinos.
	#define Nrecusions_nu_clustering 1 * (do_clustering_tag > 0) //how many recursions of finding Mnu from R and viceversa we do.

	#define Mhalo_min_clustering   1e14//in Msuns
	#define mnu_min_clustering 0.3
		/*
		minimum neutrino mass to consider doing clustering, at M_halo = Mhalo_min_clustering Msun,
		we scale it as 1/sqrt(Mhalo); [1310.6459]
		multiplied by T_extra/Temp_nu for extra species (since it's m/T that matters).

		Feel free to change these parameters to avoid doing the clustering integral.
		If you only care about the scale dependence of b(k) you should be fine.
		With these parameters you should still preserve sub-percent-level precision in the (total) bias. The effect is mostly scale independent.
		*/


	#define Ninput 22 //how many inputs the code takes
	#define Ninput_int 5 //how many of the inputs are integers instead of doubles.







//other cosmological constants that do not change:
	#define omega_constant  4.4806e-7 //Boltzmann factor divided by rho crit * h^2
	#define kpivot 0.05 //in Mpc-1, chosen to match Planck2015

	#define m_SN_thresh_Neff  1.0 // in eV; threshold to consider sterile neutrinos as Neff or not (to subtract 1 from Neff).
	//if 3+SN have m_SN_thresh_Neff > m_SN, otherwise not.

	#define mass_constant_nu 93.14 //m_nu/\omega_nu in eV, good for T_nu = T_\gamma * 0.71599
	#define ratio_T_nu_T_gamma 0.71599 //recommended by CLASS for more accurately approximate neutrino decoupling, if set to (4/11)^(1/3) change ratio_T_nu_T_gamma to 94

	#define constant_change_Neff 1.0132 //this is by how much Neff changes when subtracting each neutrino. not exactly one (see CLASS).



//halo mass functions we implement:
	#define _MICE_ 1
	#define _wCDM_ 2
	#define _ST_ 3
	#define HMF_option _MICE_ //Option for what Halo Mass Function to use. 1:MICE, 2:wLCDM, 3:ST. Add your own in Bias.h






// Other code parameters:

	#define zi 200.0 //initial redshift for collapse. Same than LoVerde2014, we ignore non-linear evolution until then other than through CAMB.
	#define zf 0.0 //haloes should collapse by today--at the latest.


	//#define length_transfer 122 //Number of elements in the transfer function output from CLASS


	//#define Nz_transfer 100 //how many redshifts we take for transfer functions. For more than ~115 need to skip lines in CLASS file, so it's not supported.



	#define print_headers _TRUE_ //whether to print headers in output files	or not

	#define z_collapse_min 0.001 //smallest z_collapse before the code complains that it has to change the binning to linear instead of logarithmic.












///Mathematical constants and unit conversions:

	#define PI 3.1415926535

	#define Msuntokm 1.47578 //Msun to km
	#define kmtoMpc  3.24e-20 //km to Mpc
	#define MsuntoMpc 4.80e-20 //Msun to Mpc.
	#define c_light  299792 //speed of light in km/s
	#define pctoly  3.262 //parsec to lightyear

	#define KtoeV 8.617e-5 //we take k_B=1, so we convert K to eV and viceversa.

	#define  hbareVs 6.582e-16 //  hbar in eV * s, we will assume hbar=1 throughout, so we have to use to convert back sometimes.
	#define  hbareVkm  1.975e-10 // hbar in eV * km
	#define  hbareVMpc  6.398e-30 // hbar in eV * Mpc

	#define  eVtokg    1.783e-36 // 1 eV in kg
	#define  Msuntokg   1.9886e30 // 1 Msun in kg
	#define  eVtoMpc   4.304e-86//eVtokg/Msuntokg*MsuntoMpc; //here we assume G=1: geometrized units.




	struct relicfast{

//parameters read as input
		int relicfast_verbose;
		int length_transfer;
		int run_relicfast;           	
		int nu_hierarchy;
 		int tag_thermal_relic; //whether we include thermal relics
		int tag_sterile_nu; //whether we include sterile neutrinos


		double m_TR; //mass of thermal relic
		double omega_TR; //energy density of TR
		double m_SN; //mass of sterile nu (T is assumed to be T0_nu)

		double As; //amplitude of perturbations
		double ns; //scalar tilt
		double m_nu1; //mass of neutrino 1
		double m_nu2; //mass of neutrino 2
		double N_eff_input; //N_effective read, we will subtract 1 per each massive neutrino to get Neff.


//derived paramters:

		double N_eff; // from Neff_input subtracting neutrinos that are massive and so on.

		double omega_SN; //sterile neutrino energy density. m_SN/94.07

		double omega_extra; //energy density in TR or SN
		double Omega_extra; //omega_extra/h^2
		double m_extra; //mass of the extra species (TR or SN)


		int counter_massive_nus; //how many massive nus (0, 1, or 2, since 3rd nu comes in as a SN)

		double omega_nu1; //energy density in nu1
		double Omega_nu1;

		double omega_nu2; //energy density in nu2
		double Omega_nu2;

		double Omega_M; //CMB+b energy density (NOT INCLUDING RELICS/NUs!). Used for halo masses

		double T0_nu;//neutrino temperature at z=0 in K. Simpler to keep 1/94.07, there are more refined models.
		double T0_TR;// Fixed through cosmic abundance, assuming relativistic at decoupling and non-rel today
		double T0_extra;// either T0_TR or TO_nu (TR or SN)



		double Omega_nu_massless;//Omega_masslessneutrinos

		double Omega_R;// the part that is radiation (at all z)



		double Omega_L; //Omega_Lambda, we close the Friedmann Eq. assuming flat Universe.


		double H0_Mpc;// H0 in Mpc-1


//collapse and other parameters:

//read:

		int N_zcoll; //how many z_collapses
		int N_Mhalo; //how many halo masses
		int N_klong;	//how many k_longs.

		double z_collapse_bot; //minimum and maximum redshifts over which we calculate. Linearly spcaed.
		double z_collapse_top;

		//double Mhalo; //the halo mass for each iteration, in Msun.
		//double Mhalo_Mpc; //G*Mhalo, in Mpc.
		double Mhalo_min; //mhalo (in Msun) minimum and maximum over which we calculate. Logspaced.
		double Mhalo_max;

		double k_long_top; //Mpc-1, maximum and minimum k over which we calculate things.
		double k_long_bot;

//derived

		double *z_collapse_array; //array with z_collapses.
		double *Mhalo_array; //array with masses that we calculate for.
		double *klong_list_input; //for which k_longs we calculate. We use slightly different values to match CLASS/CAMB for higher precision.
		double ***b_L; //array of Lagrangian biases
		double ***b_E; //array of Eulerian biases
		double *bias_array; //array of Lagrangian biases for interpolation
		double ***Pmm; 
		double ***Pmh;
		double ***Phh;
		double *klong_list;
		double ****delta_long_collapse;
		double ****delta_short_collapse;
		double ****delta_short_crit;
		double *delta_long_list;
		double *zlist_transfer;
		int Nz_transfer;
		ErrorMsg error_message;

	};


/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


int relicfast_bias_at_z_Mmin_and_k(
			  struct relicfast * prf,
                          double z,
			  double k,
                          double *b
                          );


int relicfast_bias_at_z_M_and_k(
			  struct relicfast * prf,
                          double z,
                          double M,
			  double k,
                          double *b
                          );


int get_transfer(struct background * pba,
              struct perturbs * ppt,
              struct spectra * psp,
              struct relicfast *prf, 
	      double *kgrid, double **TF, double z);

int get_bulk_transfer(struct background * pba,
              struct perturbs * ppt,
              struct spectra * psp,
              struct relicfast *prf, 
	      double *kgrid,
	      double **TFm, double **TFgamma, double **TFnu_massless, double **TFnu1, double **TFnu2, double **TFextra);



double pressure_WDM(double mass, double Temp);

double density_WDM(double mass, double Temp);

double EoS_WDM(double mass, double Temp);


int get_bias(struct precision * ppr,
			struct perturbs * ppt,
			struct background * pba,
              		struct spectra * psp,
			struct relicfast *prf, double z_collapse);

//\delta_crit as a function of \delta_long]
int collapse(struct precision * ppr,
			struct perturbs * ppt,
			struct background * pba,
              		struct spectra * psp,
			struct relicfast * prf, double z_collapse);


//these are the functions that find z_collapse as a function of the initial conditions:
//for nothing massive (nu or extra)
double find_z_collapse_nothing
(struct background * pba, struct relicfast *prf, double Ri, double Rpi, double delta_long, 
double Tfm_klong, double *transfer_gamma_klong, double *transfer_nu_klong, double z_collapse, double Mhalo_Mpc
);

//for one massive species (nu or extra)
double find_z_collapse_1nu
(struct background * pba, struct relicfast *prf, double Ri, double Rpi, double delta_long, 
double Tfm_klong, double *transfer_gamma_klong, double *transfer_nu_klong, double *transfer_nu_massive_klong,
double zmin_EoS, double dz_EoS, long Nz_EoS, double *rholist_EoS, double *plist_EoS,
long Nz_solution, double *R_solution, double *Mnu_solution, double z_collapse, double Mhalo_Mpc
);

//for 2 massive species (nu1+nu2 OR nu1+extra)
double find_z_collapse_2nu
(struct background * pba, struct relicfast *prf, double Ri, double Rpi, double delta_long, 
double Tfm_klong, double *transfer_gamma_klong, double *transfer_nu_klong, double *transfer_nu1_klong, double *transfer_nu2_klong,
double zmin_EoS, double dz_EoS, long Nz_EoS, double *rholist1_EoS, double *plist1_EoS, double *rholist2_EoS, double *plist2_EoS,
long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution, double z_collapse
, double Mhalo_Mpc);

//for 2 massive neutrinos + extra
double find_z_collapse_3nu
(struct background * pba, struct relicfast *prf, double Ri, double Rpi, double delta_long,
double Tfm_klong, double *transfer_gamma_klong, double *transfer_nu_klong, double *transfer_nu1_klong, double *transfer_nu2_klong, double *transfer_extra_klong,
double zmin_EoS, double dz_EoS, long Nz_EoS, double *rholist1_EoS, double *plist1_EoS, double *rholist2_EoS, double *plist2_EoS, double *rholist_extra_EoS, double *plist_extra_EoS,
long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution, double *Mnu3_solution, double z_collapse, double Mhalo_Mpc
);

//this gets the solution for R(t) and finds the Mnu(t) within Rhalo. One for each species
//(could be made more efficient by solving together, but it is irrelevant for neutrinos, so we leave it as is)
int findmass_1nu (struct background * pba, struct relicfast *prf, double mass_0, double temp_0, double zmin_EoS, double dz_EoS, long Nz_EoS,
	double* rho1list_EoS, double* rho2list_EoS, double* rho3list_EoS, double* etalist,
	double logz_solution_min, double dlogz_solution, long Nz_solution, double* R_solution, double* Mnu_solution, double z_collapse, double Mhalo_Mpc);






//array functions:

double *allocate_1D_array(int n1);

int *allocate_1D_array_int(int n1);


double **allocate_2D_array(int n1, int n2);

void free_array_2D(double **array, int n1);


double ***allocate_3D_array(int n1, int n2, int n3);

void free_array_3D(double ***array, int n1,int n2);


double ****allocate_4D_array(int n1, int n2, int n3, int n4);

void free_array_4D(double ****array, int n1,int n2, int n3);


void do_check(int var_to_check);


long find_value(long length, double xtab[], double x);

long find_value_reverse(long length, double xtab[], double x);

void reverse(double *list,int N);


//mathemtical functions:

double interpol(double data[], double xtab[],int length, double x);

double interpol_2D(double **data, double xtab[], int lengthx, double ytab[], int lengthy, double x, double y);

double nintegrate(double data[],double xtab[], int length);

int triangle(double l1,double l2,double l3);

double interpol_cubic(double x0, double dx, double *ytab, int Nx, double x);

void solve_syst(double  **A, double  *X, double  *B, int N);

double det(double **M, int N);


//for progressbar:
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
void printProgress (double percentage);


//H/H0 and its log derivative for LCDM only
double E_LCDM(struct relicfast *prf, double z);
double dlogE_dz_LCDM(struct relicfast *prf, double z);


//time as a function of z with 1 additional rho.
void find_time(struct relicfast *prf, double *zlist_Fri, long Nz_Fri, double *tlist_Fri,
				 double zmin_EoS, double dz_EoS, long Nz_EoS, double *rholist_EoS);


//Compute growth factor (7.73 in Dodelson's book)
void find_growth(struct relicfast *prf, double *z_growth_array, double *Dp_growth_array, long N_growth);

double getsigma_M(struct relicfast *prf, double *k_transfer, double *transfer_function, double Mhalo_Mpc);
  //calculates sigma(M)
double getsigma_8(struct background *pba, struct relicfast *prf, double *k_transfer, double *transfer_function);
	//for R=8 Mpc/h. It's sigma_8=0.83 for Planck2015 best-fit

double getsigma_M_z(struct relicfast *prf, double z, double *z_transfer,  double *k_transfer,  double **transfer_function, double Mhalo_Mpc);
//calculates sigma(M) @z.


double getsigma_8_z(struct background *pba, struct relicfast *prf, double z, double *z_transfer, double *k_transfer, double **transfer_function);
//calculates sigma8 @z.


void findeta(struct background *pba, struct relicfast *prf, double *etalist, double zmin_EoS, double dz_EoS, long Nz_EoS,
            double* rholist_nu1_EoS, double* rholist_nu2_EoS, double* rholist_extra_EoS);
  //finds superconformal time etalist, at z_EoS list. Given LCDM + neutrinos and extra species.


int relicfast_init(
			struct precision * ppr,
			struct perturbs * ppt,
			struct background * pba,             		
              		struct spectra * psp,
			struct relicfast * prf
		   );



int relicfast_free(
			struct relicfast * prf
		   );


#ifdef __cplusplus
}
#endif
#endif
