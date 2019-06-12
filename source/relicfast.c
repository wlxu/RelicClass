
#include "relicfast.h"

int relicfast_bias_at_z_Mmin_and_k(
 struct relicfast * prf,
 double z,
 double k,
 double * b
 ) {
  class_call(relicfast_bias_at_z_M_and_k(prf,z,prf->Mhalo_min,k,b), prf->error_message, prf->error_message);

  return _SUCCESS_;

}



int relicfast_bias_at_z_M_and_k(
 struct relicfast * prf,
 double z,
 double M,
 double k,
 double * b
 ) {

  /** Summary: */

  /** - define local variables */


  int last_index;

  double * bias_at_z = NULL;
  double * bias_at_z_and_M = NULL;
  double * spline_z;
  double * spline_M_at_z;
  double * spline_k_at_z_and_M;

  //  double * b_result;


  /** - first step: check that k and z are in valid range */

  class_test((k < prf->k_long_bot) || (k > prf->k_long_top),
   prf->error_message,
   "k=%e out of bounds [%e:%e]",k,prf->k_long_bot,prf->k_long_top);

  
  class_test((z < prf->z_collapse_bot) || (z > prf->z_collapse_top),
   prf->error_message,
   "z=%e out of bounds [%e:%e]",z,prf->z_collapse_bot,prf->z_collapse_top);


  class_test((M < prf->Mhalo_min) || (M > prf->Mhalo_max),
   prf->error_message,
   "M=%e out of bounds [%e:%e]",M,prf->Mhalo_min,prf->Mhalo_max);


  //    class_alloc(b_result,
  //            sizeof(double),
  //            prf->error_message);


  class_alloc(bias_at_z,
    prf->N_klong*prf->N_Mhalo*sizeof(double),
    prf->error_message);

  class_alloc(bias_at_z_and_M,
    prf->N_klong*sizeof(double),
    prf->error_message);



  class_alloc(spline_z,
    sizeof(double)*prf->N_zcoll*prf->N_klong*prf->N_Mhalo,
    prf->error_message);
  class_alloc(spline_M_at_z,
    sizeof(double)*prf->N_klong*prf->N_Mhalo,
    prf->error_message);
  class_alloc(spline_k_at_z_and_M,
    sizeof(double)*prf->N_klong,
    prf->error_message);


  

    /*bias at z*/
  if(prf->N_zcoll>1){

    class_call(array_spline_table_lines(prf->z_collapse_array,
      prf->N_zcoll,
      prf->bias_array,
      prf->N_Mhalo*prf->N_klong,
      spline_z,
      _SPLINE_EST_DERIV_,
      prf->error_message),
    prf->error_message,
    prf->error_message);


    class_call(array_interpolate_spline(prf->z_collapse_array,
      prf->N_zcoll,
      prf->bias_array,
      spline_z,
      prf->N_Mhalo*prf->N_klong,
      z,
      &last_index,
      bias_at_z,
      prf->N_Mhalo*prf->N_klong,
      prf->error_message),
    prf->error_message,
    prf->error_message);
  }
  else{
   bias_at_z = prf->bias_array;
 }

    /*bias at z and M*/
 if(prf->N_Mhalo>1){
  class_call(array_spline_table_lines(prf->Mhalo_array,
    prf->N_Mhalo,
    bias_at_z,
    prf->N_klong,
    spline_M_at_z,
    _SPLINE_EST_DERIV_,
    prf->error_message),
  prf->error_message,
  prf->error_message);



  class_call(array_interpolate_spline(prf->Mhalo_array,
    prf->N_Mhalo,
    bias_at_z,
    spline_M_at_z,
    prf->N_klong,
    M,
    &last_index,
    bias_at_z_and_M,
    prf->N_klong,
    prf->error_message),
  prf->error_message,
  prf->error_message);
}
else{
  bias_at_z_and_M = bias_at_z;
}



    /*bias at z M and k*/
class_call(array_spline_table_lines(prf->klong_list,
  prf->N_klong,
  bias_at_z_and_M,
  1,
  spline_k_at_z_and_M,
  _SPLINE_EST_DERIV_,
  prf->error_message),
prf->error_message,
prf->error_message);



class_call(array_interpolate_spline(prf->klong_list,
  prf->N_klong,
  bias_at_z_and_M,
  spline_k_at_z_and_M,
  1,
  k,
  &last_index,
  b,
  1,
  prf->error_message),
prf->error_message,
prf->error_message);



    //    b = b_result[0];

free(bias_at_z);
free(bias_at_z_and_M);
free(spline_z);
free(spline_M_at_z);
free(spline_k_at_z_and_M);
    // free(b_result);



return _SUCCESS_;

}


  /** Summary: */

  /** Effective Relicfast Outputs -- for each choice of z_collapse, M_halo, and k_long
  /* Collapse Part: 
  /* delta_short_crit (extrapolated to z_collapse)
  /* delta_short_collapse (extrapolated to z_collapse)
  /* delta_long_collapse (extrapolated to z_collapse)
  /* delta_long_list 
  /* Bias Part: 
  /* lagrangian and Eulerian biases:
  /* b_L = d log(n)/d \delta_crit * d \delta_crit / d \delta_long 
  /* b_E = Phm/prf->Pmm.
  /* bias : Eulerian bias formatted for interpolation, formatted bias[z_index*N_klong*N_Mhalo + M_index*N_klong + k_index]
  /* Pmm, Pmh, Phh: Power spectra at different k_longs
  **/ 


int relicfast_init(	struct precision * ppr,
 struct perturbs * ppt,
 struct background * pba,
 struct spectra * psp,
 struct relicfast * prf
 ){


  if(prf->run_relicfast==_FALSE_){
   if(prf->relicfast_verbose>0){printf("Relicfast not requested. Module skipped\n");}
   return _SUCCESS_;
 }


  	#ifdef _OPENMP
		omp_set_nested(1); //we need to nest the two parallel loops since the limits on the inner for() depend on the outer one, so the scheduler does not know how to do collapse(2).
	#endif



  prf->H0_Mpc=100.*pba->h/(c_light); //H0 in Mpc-1





  class_alloc(prf->klong_list,prf->N_klong*sizeof(double), prf->error_message);
  
  int iM, iz, z_ind, M_ind, k_ind;
  prf->length_transfer = psp->ln_k_size;


  class_alloc(prf->b_L,prf->N_zcoll*sizeof(double**),prf->error_message);
  class_alloc(prf->b_E,prf->N_zcoll*sizeof(double**),prf->error_message);

  class_alloc(prf->Pmm,prf->N_zcoll*sizeof(double**),prf->error_message);
  class_alloc(prf->Pmh,prf->N_zcoll*sizeof(double**),prf->error_message);
  class_alloc(prf->Phh,prf->N_zcoll*sizeof(double**),prf->error_message);

  int al_ind, al_ind1, al_ind2;
  for(al_ind=0;al_ind<prf->N_zcoll;al_ind++){
   class_alloc(prf->b_L[al_ind],prf->N_Mhalo*sizeof(double*),prf->error_message);
   class_alloc(prf->b_E[al_ind],prf->N_Mhalo*sizeof(double*),prf->error_message);
   class_alloc(prf->Pmm[al_ind],prf->N_Mhalo*sizeof(double*),prf->error_message);
   class_alloc(prf->Pmh[al_ind],prf->N_Mhalo*sizeof(double*),prf->error_message);
   class_alloc(prf->Phh[al_ind],prf->N_Mhalo*sizeof(double*),prf->error_message);	
   for(al_ind1=0;al_ind1<prf->N_Mhalo;al_ind1++){
    class_alloc(prf->b_L[al_ind][al_ind1],prf->N_klong*sizeof(double),prf->error_message);
    class_alloc(prf->b_E[al_ind][al_ind1],prf->N_klong*sizeof(double),prf->error_message);
    class_alloc(prf->Pmm[al_ind][al_ind1],prf->N_klong*sizeof(double),prf->error_message);
    class_alloc(prf->Pmh[al_ind][al_ind1],prf->N_klong*sizeof(double),prf->error_message);
    class_alloc(prf->Phh[al_ind][al_ind1],prf->N_klong*sizeof(double),prf->error_message);
  }	
}

class_alloc(prf->bias_array, prf->N_zcoll*prf->N_Mhalo*prf->N_klong*sizeof(double),prf->error_message);
class_alloc(prf->delta_long_list,N_delta_long*sizeof(double),prf->error_message);

class_alloc(prf->delta_short_collapse,prf->N_zcoll*sizeof(double***),prf->error_message);
class_alloc(prf->delta_short_crit,prf->N_zcoll*sizeof(double***),prf->error_message);
class_alloc(prf->delta_long_collapse,prf->N_zcoll*sizeof(double***),prf->error_message);

for(al_ind=0;al_ind<prf->N_zcoll;al_ind++){
	class_alloc(prf->delta_short_collapse[al_ind],prf->N_Mhalo*sizeof(double**),prf->error_message);
	class_alloc(prf->delta_short_crit[al_ind],prf->N_Mhalo*sizeof(double**),prf->error_message);
	class_alloc(prf->delta_long_collapse[al_ind],prf->N_Mhalo*sizeof(double**),prf->error_message);
	
  for(al_ind1=0;al_ind1<prf->N_Mhalo;al_ind1++){
    class_alloc(prf->delta_short_collapse[al_ind][al_ind1],N_delta_long*sizeof(double*),prf->error_message);
    class_alloc(prf->delta_short_crit[al_ind][al_ind1],N_delta_long*sizeof(double*),prf->error_message);
    class_alloc(prf->delta_long_collapse[al_ind][al_ind1],N_delta_long*sizeof(double*),prf->error_message);


    for(al_ind2=0;al_ind2<N_delta_long;al_ind2++){
     class_alloc(prf->delta_short_collapse[al_ind][al_ind1][al_ind2],prf->N_klong*sizeof(double),prf->error_message);
     class_alloc(prf->delta_short_crit[al_ind][al_ind1][al_ind2],prf->N_klong*sizeof(double),prf->error_message);
     class_alloc(prf->delta_long_collapse[al_ind][al_ind1][al_ind2],prf->N_klong*sizeof(double),prf->error_message);

   }	
 }	
}	


#pragma omp parallel for 
for(iz=0;iz< prf->N_zcoll;iz++){
  int collapse_check = collapse(ppr, ppt, pba,psp, prf,prf->z_collapse_array[iz]);
  if (collapse_check == _FAILURE_){ printf("Problem in Relicfast: Collapse Step\n");} 
  int bias_check = get_bias(ppr, ppt, pba,psp, prf,prf->z_collapse_array[iz]);
  if (bias_check == _FAILURE_){ printf("Problem in Relicfast: Bias Step\n");} 
}

for(z_ind=0;z_ind<prf->N_zcoll;z_ind++){
	for(M_ind=0;M_ind<prf->N_Mhalo;M_ind++){
		for(k_ind=0;k_ind<prf->N_klong;k_ind++){
			prf->bias_array[z_ind*prf->N_klong*prf->N_Mhalo + M_ind*prf->N_klong + k_ind]=prf->b_L[z_ind][M_ind][k_ind];
    }
  }
}


return _SUCCESS_; 

}


int relicfast_free(
 struct relicfast * prf
 ) {
  if (prf->run_relicfast == _TRUE_){
   free(prf->Mhalo_array);
   free(prf->z_collapse_array);
   free(prf->klong_list_input);
   free(prf->klong_list);
   free_array_4D(prf->delta_long_collapse, prf->N_zcoll, prf->N_Mhalo, N_delta_long);
   free_array_4D(prf->delta_short_collapse,prf->N_zcoll, prf->N_Mhalo, N_delta_long);
   free_array_4D(prf->delta_short_crit,prf->N_zcoll, prf->N_Mhalo, N_delta_long);
   free(prf->delta_long_list);
   free_array_3D(prf->b_L,prf->N_zcoll, prf->N_Mhalo);
   free_array_3D(prf->b_E,prf->N_zcoll, prf->N_Mhalo);
   free_array_3D(prf->Pmm,prf->N_zcoll, prf->N_Mhalo);
   free_array_3D(prf->Pmh,prf->N_zcoll, prf->N_Mhalo);
   free_array_3D(prf->Phh,prf->N_zcoll, prf->N_Mhalo);
   free(prf->bias_array);
 }

 return _SUCCESS_;
}


int get_transfer(struct background * pba,
  struct perturbs * ppt,
  struct spectra * psp,
  struct relicfast *prf, 
  double *kgrid, double **TF, double z){



	int j;
	int size_data;

	double temp[psp->ln_k_size];

	double TF_grid_b[psp->ln_k_size];
	double TF_grid_c[psp->ln_k_size];
	double TF_grid_gamma[psp->ln_k_size];
	double TF_grid_extra[psp->ln_k_size];
	double TF_grid_nu_massless[psp->ln_k_size];
	double TF_grid_nu1[psp->ln_k_size];
	double TF_grid_nu2[psp->ln_k_size];
	double koverh[psp->ln_k_size]; 
	int Ncolum;

  int i;
  char titles[_MAXTITLESTRINGLENGTH_]={0};
  double *data;
  int index_extra;


  index_extra = 6 + (prf->m_nu1>0) + (prf->m_nu2>0) - 1;

  class_call(spectra_output_tk_titles(pba,ppt,class_format,titles),
   pba->error_message,
   prf->error_message);
  Ncolum = get_number_of_titles(titles);


  size_data = Ncolum*psp->ln_k_size;
  class_alloc(data, sizeof(double)*psp->ic_size[ppt->index_md_scalars]*size_data, prf->error_message);




  class_call(spectra_output_tk_data(pba,
    ppt,
    psp,
    class_format,
    z,
    Ncolum,
    data
    ),
  psp->error_message, prf->error_message);




  for(j=0; j<size_data/Ncolum; j++){
    koverh[j] = data[j*Ncolum + 0];
    TF_grid_gamma[j] =data[j*Ncolum + 1];
    TF_grid_b[j] =data[j*Ncolum + 2];
    TF_grid_c[j] = data[j*Ncolum + 3];
    TF_grid_nu_massless[j] = data[j*Ncolum + 4];
    if(prf->m_nu1>0){
      TF_grid_nu1[j] = data[j*Ncolum + 5];
      if(prf->m_nu2>0){
        TF_grid_nu2[j] = data[j*Ncolum + 6]; }
        else{TF_grid_nu2[j] = 0.0;}		
      }
      else if(prf->m_nu2>0){
        TF_grid_nu1[j] = 0.0;
        TF_grid_nu2[j] = data[j*Ncolum + 5]; 
      }
      else{
        TF_grid_nu1[j] = 0.0;
        TF_grid_nu2[j] = 0.0;
      }
      if(prf->m_extra>0){
        TF_grid_extra[j] = data[j*Ncolum + index_extra];
      }
      else{
        TF_grid_extra[j] = 0.0;
      }
    }

    for(j=0;j<size_data/Ncolum; ++j){
     kgrid[j] = pba->h * koverh[j];
     TF[0][j] = pba->Omega0_b/(pba->Omega0_b+ pba->Omega0_cdm) * TF_grid_b[j] + pba->Omega0_cdm/(pba->Omega0_b+pba->Omega0_cdm) * TF_grid_c[j];
     TF[1][j] = TF_grid_gamma[j];
     TF[2][j] = TF_grid_nu_massless[j];
     TF[3][j] = TF_grid_nu1[j];
     TF[4][j] = TF_grid_nu2[j];
     TF[5][j] = TF_grid_extra[j];
   }


   free(data);
   return _SUCCESS_;

 }


 int get_bulk_transfer(struct background * pba,
  struct perturbs * ppt,
  struct spectra * psp,
  struct relicfast *prf, 
  double *kgrid,
  double **TFm, double **TFgamma, double **TFnu_massless, double **TFnu1, double **TFnu2, double **TFextra){


   int jz, iz;
   double **TF_grid;
   int al_ind;
   class_alloc(TF_grid,6*sizeof(double*),prf->error_message);
   for(al_ind=0;al_ind<6;al_ind++){
    class_alloc(TF_grid[al_ind],psp->ln_k_size*sizeof(double),prf->error_message);
  }


  for(jz=0;jz<prf->Nz_transfer;jz++){
    class_call(get_transfer(pba,ppt,psp, prf, kgrid, TF_grid,prf->zlist_transfer[jz]),
     prf->error_message,
     prf->error_message);
    for (iz=0; iz<prf->length_transfer; iz++){
      TFm[jz][iz] = TF_grid[0][iz];
      TFgamma[jz][iz] = TF_grid[1][iz];
      TFnu_massless[jz][iz] = TF_grid[2][iz];
      TFnu1[jz][iz] = TF_grid[3][iz];
      TFnu2[jz][iz] = TF_grid[4][iz];
      TFextra[jz][iz] = TF_grid[5][iz];
    }
  }

  free_array_2D(TF_grid,6);


  return _SUCCESS_;

}


/**
 * This routine computes \delta_crit at collapse for Mhalo at zcollapse.
 *
 * @param prf  Input/output: pointer to relicfast structure
 * @return the error status
 */




/***COLLAPSE MODULE**/

int collapse(struct precision * ppr,
 struct perturbs * ppt,
 struct background * pba,
 struct spectra * psp,
 struct relicfast * prf, double z_collapse){

	int i,j;



//we establish if for the light relic and Mhalo we have to do clustering or not. Threshold given by mnu_min_clustering @ Mhalo=Mhalo_min_clustering.
	int i_clustering_index;
	int do_clustering;


/** Read the Transfer functions from the spectra module*/

	double *k_transfer_array; // k for transfer functions (common for all) (NOT over h)

	double **transfer_matter,**transfer_nu_massless,**transfer_nu1,**transfer_gamma;// transfer functions (not over k^2), for matter, massless nu, and photons.


	class_alloc(k_transfer_array,prf->length_transfer*sizeof(double),prf->error_message);
	class_alloc(transfer_matter,prf->Nz_transfer*sizeof(double*),prf->error_message);
	class_alloc(transfer_gamma,prf->Nz_transfer*sizeof(double*),prf->error_message);
	class_alloc(transfer_nu_massless,prf->Nz_transfer*sizeof(double*),prf->error_message);
	class_alloc(transfer_nu1,prf->Nz_transfer*sizeof(double*),prf->error_message);

	//have to create these arrays anyway, to read transfer files. They can be nonsense if m_nu2=0 or no extra species, is ok.
	double **transfer_nu2;
	class_alloc(transfer_nu2,prf->Nz_transfer*sizeof(double*),prf->error_message);

	double **transfer_extra;
	class_alloc(transfer_extra,prf->Nz_transfer*sizeof(double*),prf->error_message);

	int al_ind;
	for(al_ind=0;al_ind<prf->Nz_transfer;al_ind++){
		class_alloc(transfer_matter[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_gamma[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_nu_massless[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_nu1[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_nu2[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_extra[al_ind],prf->length_transfer*sizeof(double),prf->error_message);	
	}




  class_call(get_bulk_transfer(pba, ppt, psp, prf, k_transfer_array, transfer_matter, transfer_gamma, transfer_nu_massless, transfer_nu1, transfer_nu2, transfer_extra),
   prf->error_message,
   prf->error_message);



	const double kmin=k_transfer_array[0];
	const double kmax=k_transfer_array[prf->length_transfer-1];	

	




//find the closest value to k_input in k_transfer_array to avoid unnecesary interpolation.
	long jk;


	for(i=0;i< prf->N_klong;i++){

		jk=find_value(prf->length_transfer, k_transfer_array, prf->klong_list_input[i]);
		prf->klong_list[i]=k_transfer_array[jk];//	printf("j=%ld, out=%le, in=%le \n", jk, prf->klong_list[i],prf->klong_list_input[i]);
	}


//We compute T_dot to have two different initial conditions for comparison.
	double *d_transfer_array_zi; //dT/dz. The matter(c+b) transfer function and k/h obtained from CAMB code
	class_alloc(d_transfer_array_zi,prf->length_transfer*sizeof(double),prf->error_message);

	const int index_zip1 = prf->Nz_transfer-1;
	const int index_zip2 = prf->Nz_transfer-3;
	const int index_zi = prf->Nz_transfer-2;
	const double zip1=prf->zlist_transfer[index_zip1]; //we save zs before and after zi
	const double zip2=prf->zlist_transfer[index_zip2];




	for(j=0;j<prf->length_transfer;j++){
		d_transfer_array_zi[j] = (transfer_matter[index_zip1][j]-transfer_matter[index_zip2][j])/(zip1-zip2); //dT/dz@zi
	}


	double z; //redshift



	double **transfer_array_z0;
	class_alloc(transfer_array_z0,6*sizeof(double*),prf->error_message);


	double **transfer_array_zi;
	class_alloc(transfer_array_zi,6*sizeof(double*),prf->error_message);


	double **transfer_array_zip1,  **transfer_array_zip2; //for CDM+b, derivative


	class_alloc(transfer_array_zip1,6*sizeof(double*),prf->error_message);

	class_alloc(transfer_array_zip2,6*sizeof(double*),prf->error_message);


	double **transfer_array_z_collapse;

	class_alloc(transfer_array_z_collapse,6*sizeof(double*),prf->error_message);


	for( al_ind=0;al_ind<6;al_ind++){
		class_alloc(transfer_array_z0[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_array_zi[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_array_zip1[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_array_zip2[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_array_z_collapse[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
	}



  class_call(get_transfer(pba, ppt, psp, prf, k_transfer_array, transfer_array_z0, prf->zlist_transfer[0]),
   prf->error_message,
   prf->error_message);

  class_call(get_transfer(pba, ppt, psp, prf, k_transfer_array, transfer_array_zi, prf->zlist_transfer[prf->Nz_transfer-1]),
   prf->error_message,
   prf->error_message);


  class_call(get_transfer(pba, ppt, psp, prf, k_transfer_array, transfer_array_zip1, zip1),
   prf->error_message,
   prf->error_message);

  class_call(get_transfer(pba, ppt, psp, prf, k_transfer_array, transfer_array_zip2, zip2),
   prf->error_message,
   prf->error_message);
  	

//we also read the transfer function at zcollapse


	const int z_ind=find_value(prf->N_zcoll,prf->z_collapse_array, z_collapse);


  class_call(get_transfer(pba, ppt, psp, prf, k_transfer_array, transfer_array_z_collapse, z_collapse),
   prf->error_message,
   prf->error_message);





/** Interpolation tables for the equation of state*/

	const int Nz_EoS=1000; //how many redshifts we take for the EoS. 1000 is fine for z

//first for nu1
	double *zlist_EoS;
	double *plist_nu1_EoS, *rholist_nu1_EoS;
	double Temp_nu;

	class_alloc(zlist_EoS,Nz_EoS*sizeof(double),prf->error_message);
	class_alloc(plist_nu1_EoS,Nz_EoS*sizeof(double),prf->error_message);
	class_alloc(rholist_nu1_EoS,Nz_EoS*sizeof(double),prf->error_message);

	const double zmin_EoS=0.;//we go down to z=0, needed to take ratios with Omega.
	const double dz_EoS=(zi-zmin_EoS)/(Nz_EoS-1);


	//if(debug_mode > 0){
	//	lengthname=sprintf(filename,"tests/rho_nu1_%d.dat",prf->file_tag);
	//	fp=fopen(filename,"w");
	//	if(print_headers!=0){
	//		fprintf(fp, "z\t rho(z)/rho(0) w(z) \n");
	//	}
	//}

	for(i=0;i<Nz_EoS;i++){
		zlist_EoS[i]=zmin_EoS+i*dz_EoS; //z linearly from 0 to zi EQUALLY SPACED! Since we will use interpol_cubic.
		Temp_nu = prf->T0_nu * (1.+zlist_EoS[i]); //neutrino temperature in K
		rholist_nu1_EoS[i] = density_WDM(prf->m_nu1,Temp_nu);
		plist_nu1_EoS[i] = pressure_WDM(prf->m_nu1,Temp_nu);		//						printf("z=%.1le,  w=%.1le, Tnu/mnu=%.1le \n",zlist_EoS[i],plist_EoS[i]/rholist_EoS[i],Temp_nu*KtoeV/(m_nu1));//
		//if(debug_mode > 0){
		//	fprintf(fp, "%le %le %le \n",zlist_EoS[i], rholist_nu1_EoS[i]/rholist_nu1_EoS[0], plist_nu1_EoS[i]/rholist_nu1_EoS[i]);
		//}
	}
	//if(debug_mode > 0) fclose(fp);


	Temp_nu = prf->T0_nu * (1.+zi);
	double rho_nu1_zi=density_WDM(prf->m_nu1,Temp_nu);
	Temp_nu = prf->T0_nu;
	double rho_nu1_z0=density_WDM(prf->m_nu1,Temp_nu);

	double rho_nu1_ratio_zi=0.;
	if(prf->m_nu1>0){
		rho_nu1_ratio_zi=rho_nu1_zi/rho_nu1_z0;
	}

	if(debug_mode > 0){
		printf("kmin=%.3le and kmax=%.3le \n",kmin,kmax);
	}


//now for nu2.
	double *plist_nu2_EoS, *rholist_nu2_EoS;
	double rho_nu2_zi,rho_nu2_z0,rho_nu2_ratio_zi=0;

	class_alloc(plist_nu2_EoS,Nz_EoS*sizeof(double),prf->error_message);
	class_alloc(rholist_nu2_EoS,Nz_EoS*sizeof(double),prf->error_message);


	//if(debug_mode > 0){
	//	lengthname=sprintf(filename,"tests/rho_nu2_%d.dat",prf->file_tag);
	//	fp=fopen(filename,"w");
	//	if(print_headers!=0){
	//		fprintf(fp, "z\t rho(z)/rho(0) w(z) \n");
	//	}
	//}

	if(prf->Omega_nu2>0){
		for(i=0;i<Nz_EoS;i++){
			zlist_EoS[i]=zmin_EoS+i*dz_EoS; //z linearly from 0 to z_collapse EQUALLY SPACED! Since we will use interpol_cubic.
			Temp_nu = prf->T0_nu * (1.+zlist_EoS[i]); //neutrino temperature in K
			rholist_nu2_EoS[i] = density_WDM(prf->m_nu2,Temp_nu);
			plist_nu2_EoS[i] = pressure_WDM(prf->m_nu2,Temp_nu);		//						printf("z=%.1le,  w=%.1le, Tnu/mnu=%.1le \n",zlist_EoS[i],plist_EoS[i]/rholist_EoS[i],Temp_nu*KtoeV/(m_nu1));

			//if(debug_mode > 0){
			//	fprintf(fp, "%le %le %le \n",zlist_EoS[i], rholist_nu2_EoS[i]/rholist_nu2_EoS[0], plist_nu2_EoS[i]/rholist_nu2_EoS[i]);
			//}
		}

		Temp_nu = prf->T0_nu * (1.+zi);
		rho_nu2_zi=density_WDM(prf->m_nu2,Temp_nu);
		Temp_nu = prf->T0_nu;
		rho_nu2_z0=density_WDM(prf->m_nu2,Temp_nu);

		rho_nu2_ratio_zi=rho_nu2_zi/rho_nu2_z0;

	}
//	if(debug_mode > 0) fclose(fp);




//and for the extra species

	double *plist_extra_EoS, *rholist_extra_EoS;
	double rho_extra_zi,rho_extra_z0,rho_extra_ratio_zi=0;

	class_alloc(plist_extra_EoS,Nz_EoS*sizeof(double),prf->error_message);
	class_alloc(rholist_extra_EoS,Nz_EoS*sizeof(double),prf->error_message);

	double T0_extra, Temp_extra, m_extra;

	if(prf->tag_sterile_nu == _TRUE_){
		T0_extra=prf->T0_nu;
		m_extra=prf->m_SN;
	}
	else if(prf->tag_thermal_relic == _TRUE_){
		T0_extra=prf->T0_TR;
		m_extra=prf->m_TR;
	}
	else{
		m_extra=1.0;//just to make sure there are no divergencies.
		T0_extra=0.0;
	}


	if(prf->Omega_extra>0){
		for(i=0;i<Nz_EoS;i++){
			zlist_EoS[i]=zmin_EoS+i*dz_EoS; //z linearly from 0 to z_collapse EQUALLY SPACED! Since we will use interpol_cubic.
			Temp_extra = prf->T0_extra * (1.+zlist_EoS[i]); //neutrino temperature in K
			rholist_extra_EoS[i] = density_WDM(prf->m_extra,Temp_extra);
			plist_extra_EoS[i] = pressure_WDM(prf->m_extra,Temp_extra);		//						printf("z=%.1le,  w=%.1le, Tnu/mnu=%.1le \n",zlist_EoS[i],plist_EoS[i]/rholist_EoS[i],Temp_nu*KtoeV/(m_nu1));

		}

		Temp_extra = T0_extra * (1.+zi);
		rho_extra_zi=density_WDM(m_extra,Temp_extra);
		Temp_extra = T0_extra;
		rho_extra_z0=density_WDM(m_extra,Temp_extra);

		rho_extra_ratio_zi=rho_extra_zi/rho_extra_z0;

	}




	//if we ask for clustering do it only if at least one of the species is heavy enough to matter.
	if(do_clustering_tag>0){
		if((prf->m_nu1 > mnu_min_clustering*sqrt(Mhalo_min_clustering/prf->Mhalo_array[0])) ){
			do_clustering=1;
		}
		else if((prf->m_extra > prf->T0_extra/prf->T0_nu * mnu_min_clustering*sqrt(Mhalo_min_clustering/prf->Mhalo_array[0]))){
			do_clustering=1;
		}
		else{
			do_clustering=0;
		}
	}
	else{
		do_clustering=0;
	}


//we find the superconformal time eta for clustering. If required.
	double *etalist;
	class_alloc(etalist,Nz_EoS*sizeof(double),prf->error_message); //we find eta in the same z range we find EoS, where it matters.

	if(do_clustering_tag>0){
		findeta(pba,prf, etalist, zmin_EoS, dz_EoS, Nz_EoS,
      rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS);
	}//creates a cubic interpolator for eta at z_EoS list.



/** Setting the R(z) array that we will use to find extra species clustering*/


	const int Nz_solution=400;//number of z in solution of R(t). LOGSPACED.
	const double logz_solution_min=log(zi);
	const double logz_solution_max=log(z_collapse/2);//as in collapse code.
	const double dlogz_solution=(logz_solution_max-logz_solution_min)/(Nz_solution-1);

	double *Rhalo_solution; //radius of halo at different zs.
	class_alloc(Rhalo_solution,Nz_solution*sizeof(double),prf->error_message);
	double *Mnu1_solution; //neutrino 1 mass within R_halo, needed to find collapse.
	class_alloc(Mnu1_solution,Nz_solution*sizeof(double),prf->error_message);
	double *Mnu2_solution; //neutrino 1 mass within R_halo, needed to find collapse.
	class_alloc(Mnu2_solution,Nz_solution*sizeof(double),prf->error_message);
	double *Mextra_solution; //neutrino 1 mass within R_halo, needed to find collapse.
	class_alloc(Mextra_solution,Nz_solution*sizeof(double),prf->error_message);
	int mnu_ind;
	for (mnu_ind=0; mnu_ind<Nz_solution ; mnu_ind ++){ 
		Mnu1_solution[mnu_ind] = 0.0;
		Mnu2_solution[mnu_ind] = 0.0;
		Mextra_solution[mnu_ind] = 0.0;

	}
	//z array logspaced between zi and zcollapse/2, to match collapse code.




/** Setting Initial Conditions*/


	double OmL_i, OmM_i, OmR_i, Omnu1_i, Omnu2_i, Omextra_i;
	OmL_i=prf->Omega_L;
	OmM_i=prf->Omega_M*pow(1.+zi,3.);
	OmR_i=prf->Omega_R*pow(1.+zi,4.);
	Omnu1_i=prf->Omega_nu1*rho_nu1_ratio_zi;
	Omnu2_i=prf->Omega_nu2*rho_nu2_ratio_zi;
	Omextra_i=prf->Omega_extra*rho_extra_ratio_zi;
	const double Hi=prf->H0_Mpc*sqrt(OmL_i + OmM_i + OmR_i +
				Omnu1_i + Omnu2_i + Omextra_i);//H(zi) in Mpc-1

	if(debug_mode>0){
		printf("Hi=%.2le \n",Hi);
		printf("OmM_i=%.2le, OmR_i=%.2le, Omnu1=%.2le \n",OmM_i,OmR_i,Omnu1_i);
		printf("Omnu2_i=%.2le, Omextra=%.2le \n",Omnu2_i,Omextra_i);
	}


	const double sigma8_z0= getsigma_8(pba, prf, k_transfer_array, transfer_matter[0]);
	const double sigma8_collapse = getsigma_8(pba, prf, k_transfer_array, transfer_array_z_collapse[0]);
	const double sigma8_i = getsigma_8(pba, prf, k_transfer_array, transfer_matter[index_zi]);





//we set delta_short.
 	const double delta_short_min_constant=delta_short_constant_centroid/(2.0+boost_initial_conditions); //initial conditions, to do bipartition. these are the benchmark.
 	const double delta_short_max_constant=delta_short_constant_centroid*(2.0+boost_initial_conditions);



	const double tolerance=0.0002/fmax(fmin(precision_scale,10.),1.); //relative error before we call it converged. 3e-4 should guarantee 0.1% in bias.
	const double tolerance_ini= tolerance + (do_clustering>0)*tolerance*4.0; //first try does not have to be as precise (if we do clustering), since it is used to find Mnu collapse only.
	const double tolerance_z = z_collapse * 0.1; //this is just to make sure that we are not artificially converging to a "bisection" \delta_crit if our initial conditions are bad.



//we set delta_long. We assume Lambda (DE) does not cluster.
//for matter first, rest will depend on z.
	double *Ti_klong; //transfer function of CDM+b @ zi
	double *dTi_klong; //derivative of transfer function of CDM+b @ zi, wrt to z
	class_alloc(Ti_klong,prf->N_klong*sizeof(double),prf->error_message);
	class_alloc(dTi_klong,prf->N_klong*sizeof(double),prf->error_message);





//for usual species:
	double **transfer_gamma_klong;//transfer evaluated at klong, for quick interpolation.
	double **transfer_nu_massless_klong;
	double **transfer_nu1_klong;
	class_alloc(transfer_gamma_klong,prf->N_klong*sizeof(double*),prf->error_message);
	class_alloc(transfer_nu_massless_klong,prf->N_klong*sizeof(double*),prf->error_message);
	class_alloc(transfer_nu1_klong,prf->N_klong*sizeof(double*),prf->error_message);

	for(al_ind=0;al_ind<prf->N_klong;al_ind++){
		class_alloc(transfer_gamma_klong[al_ind],prf->Nz_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_nu_massless_klong[al_ind],prf->Nz_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_nu1_klong[al_ind],prf->Nz_transfer*sizeof(double),prf->error_message);
	}

	int i_klong;
//we set the initial cdm+b transfer functions, and the interpolators for all fluids.
		for(i_klong=0;i_klong<prf->N_klong;i_klong++){ //we find transfer function at initial redshift
			double k_long = prf->klong_list[i_klong];
			Ti_klong[i_klong]=interpol_2D(transfer_matter, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, zi, k_long);

			dTi_klong[i_klong]=(interpol_2D(transfer_matter, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, zip2, k_long)
       -interpol_2D(transfer_matter, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, zip1, k_long))/(zip2-zip1);

			for(i=0;i<prf->Nz_transfer;i++){
				z=prf->zlist_transfer[i];
				transfer_gamma_klong[i_klong][i]=interpol_2D(transfer_gamma, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, z, k_long);
				transfer_nu_massless_klong[i_klong][i]=interpol_2D(transfer_nu_massless, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, z, k_long);	//			printf("i_klong=%d, i=%ld, z=%.3le, Tf_gamma=%.1le, Tf_nu=%.1le \n",i_klong, i, z,transfer_gamma_klong[i_klong][i],transfer_nu_massless_klong[i_klong][i]);
				transfer_nu1_klong[i_klong][i]=interpol_2D(transfer_nu1, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, z, k_long);	//			printf("i_klong=%d, i=%ld, z=%.3le, Tf_gamma=%.1le, Tf_nu=%.1le \n",i_klong, i, z,transfer_gamma_klong[i_klong][i],transfer_nu1_klong[i_klong][i]);//			printf("k_long=%.1le, z=%.1le , T_g=%.1le \n",k_long, z, transfer_gamma_klong[i_klong][i]);
			}
		}


	// and for extra species. No need to add clause: if(m_nu2>0), since it will just be zeroes if so.
		double **transfer_nu2_klong;
		double **transfer_extra_klong;
		class_alloc(transfer_nu2_klong,prf->N_klong*sizeof(double*),prf->error_message);
		class_alloc(transfer_extra_klong,prf->N_klong*sizeof(double*),prf->error_message);

		for(al_ind=0;al_ind<prf->N_klong;al_ind++){
      class_alloc(transfer_nu2_klong[al_ind],prf->Nz_transfer*sizeof(double),prf->error_message);
      class_alloc(transfer_extra_klong[al_ind],prf->Nz_transfer*sizeof(double),prf->error_message);
    }

    for(i_klong=0;i_klong<prf->N_klong;i_klong++){
     double k_long = prf->klong_list[i_klong];
     for(i=0;i<prf->Nz_transfer;i++){
      z=prf->zlist_transfer[i];
				transfer_nu2_klong[i_klong][i]=interpol_2D(transfer_nu2, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, z, k_long);	//			printf("i_klong=%d, i=%ld, z=%.3le, Tf_gamma=%.1le, Tf_nu=%.1le \n",i_klong, i, z,transfer_gamma_klong[i_klong][i],transfer_nu1_klong[i_klong][i]);
				transfer_extra_klong[i_klong][i]=interpol_2D(transfer_extra, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, z, k_long);	//
			}
		}




	///////////////////////////////////////////////////////////////////////////////////
	////    we iterate over klong array defined above, as well as Mhalos					////
	/////////////////////////////////////////////////////////////////////////////////

		int iM;
		const int N_klong_calc [N_delta_long] = {1, prf->N_klong}; //what k_longs we calculate over. For delta_long=0 only one is needed, otherwise one per kmode.
		//N_delta_long is defined in common.h



		const double delta_long_min=0.; //the maximum value is specified in common.h for convenience to change it.
		const double delta_long_step=(delta_long_max-delta_long_min)/(N_delta_long-1.);

		for(i=0;i<N_delta_long;i++){
			prf->delta_long_list[i]=delta_long_min+delta_long_step*i;
		}



		






	double T_z_collapse_klong; //transfer function T(k_long) at zcollapse, to extrapolate \delta_long.



#pragma omp parallel for \
     private(iM) \
     schedule (static)
  for(iM=0;iM<prf->N_Mhalo;iM++){

   int m_ind = iM;
   double Mhalo = prf->Mhalo_array[iM];
	double Mhalo_Mpc = Mhalo * MsuntoMpc;//Mhalo in Mpc (*G)



//average R_i
	double Ribar= pow(prf->Omega_M*pow(1.+zi,3.)/2.*prf->H0_Mpc*prf->H0_Mpc/Mhalo_Mpc,-1./3); //R in Mpc.


	//sigma(M) at zi, z)collapse, and derivative (and z=0 just in case), for initial conditions of ODE. (option_initial_velocity=0)
	double sigmaM_0 = getsigma_M(prf,  k_transfer_array, transfer_matter[0], Mhalo_Mpc);
	double sigmaM_collapse = getsigma_M(prf,  k_transfer_array, transfer_array_z_collapse[0], Mhalo_Mpc);
	double sigmaM_i = getsigma_M(prf,  k_transfer_array, transfer_matter[index_zi], Mhalo_Mpc);
	double sigmaM_zip1 = getsigma_M(prf,  k_transfer_array, transfer_matter[index_zip1], Mhalo_Mpc);
	double sigmaM_zip2 = getsigma_M(prf,  k_transfer_array, transfer_matter[index_zip2], Mhalo_Mpc);
	double d_sigmaM_i = (sigmaM_zip2 - sigmaM_zip1)/(zip2-zip1);


	//if(debug_mode>=0) printf("sigmaM_i=%lf, sigmaM_colapse = %lf \n",sigmaM_i, sigmaM_collapse);



	//We also use the transfer functions (option_initial_velocity=1)
	double kstar = 1.0/pow(prf->Omega_M/2.*prf->H0_Mpc*prf->H0_Mpc/Mhalo_Mpc,-1./3); //Mpc-1; where we evaluate Tdot/T.

//we find the k that is closest to kstar, since that is easier than interpolating for no reason.
	int indexkstar = find_value(prf->length_transfer, k_transfer_array, kstar);
	double Ti_kstar=transfer_array_zi[0][indexkstar];
	double dTi_kstar=d_transfer_array_zi[indexkstar];
	double T_z_collapse_kstar=transfer_array_z_collapse[0][indexkstar]; //transfer function T(k_pivot) at zcollapse, to extrapolate \delta_crit.






	if(debug_mode>0){
		printf("IC_1:%.3le  IC_2:%.3le \n", d_sigmaM_i/sigmaM_i, dTi_kstar/Ti_kstar);
		printf("Ribar=%.3le \n", Ribar);

		printf("We will run some checks: \n");
		printf("@z=%.1le \t sigma_M= %.3le \n",0., sigmaM_0);
		printf("@z=%.1le \t sigma_8= %.3le \n",0., sigma8_z0);
		printf("@z=%.1le \t sigma_M= %.3le \n", z_collapse, sigmaM_collapse);
		printf("@z=%.1le \t sigma_8= %.3le \n", z_collapse, sigma8_collapse);
		printf("@z=%.1le \t sigma_M= %.3le \n",zi, sigmaM_i);
		printf("@z=%.1le  d sigma_M/dz= %.3le \n",zi, d_sigmaM_i);
		printf("@z=%.1le \t sigma_8= %.3le \n",zi, sigma8_i);
	}



	int i_delta_long;
	#pragma omp parallel for \
     private(i_delta_long) \
     schedule (static)
	for(i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
		#pragma omp parallel for \
     private(i_klong) \
     schedule (static)
		for(i_klong=0;i_klong<N_klong_calc[i_delta_long];i_klong++){



			double delta_long=prf->delta_long_list[i_delta_long]; //long-wavelength perturbation for CDM+b , for photons , and for massless neutrinos.


			double k_long=prf->klong_list[i_klong];

			double ddelta_long_dt=delta_long*(-Hi*(1.+zi))*dTi_klong[i_klong]/Ti_klong[i_klong]; //d delta_long/dt, calculated as dT/dt*1/T*delta_long.
			//delta_short (initial) to collapse for each delta_long and k_long.

			double delta_short_min=delta_short_min_constant; //we restart the initial conditions, and redo the procedure.
			double delta_short_max=delta_short_max_constant;
			double delta_short=0.0,ddelta_short_dt; //private copies for each run
			double Ri, Ridot, Rpi, z_coll_iteration=0.0; //initial conditions and z of collapse of each iteration.




			if(debug_mode >1 ){
				printf("deltashort_min=%.2le, deltashort_max=%.2le \n", delta_short_min, delta_short_max);
			}

			for(;delta_short_max/delta_short_min-1.>tolerance_ini;){


				delta_short=(delta_short_max+delta_short_min)/2.0; //bisection, we update the min or max value.

				if(option_initial_velocity==0){
			///Option 1: through sigma(M)
					ddelta_short_dt=delta_short*(-Hi*(1.+zi))*d_sigmaM_i/sigmaM_i; //since: 		ddelta_dt/delta = dsigma_dt/sigma = -H*(1+z) dsigma_dz/sigma
				}
				else if(option_initial_velocity==1) {
			///Option 2: through transfer function (they agree)
					ddelta_short_dt=delta_short*(-Hi*(1.+zi))*dTi_kstar/Ti_kstar; //since:		ddelta_dt/delta = dT_dt/T = -H*(1+z) dT_dz/sigma
				}
				else {
					printf("ERROR: option_initial_velocity has to be either 0 or 1 \n");
					exit(0);
				}


			//these are the i.c. we will feed the function find_z_collapse_XXX.
			Ri= Ribar*(1.0-1./3.*(delta_short+delta_long)); //R_i in Mpc.
			Ridot= Ri*Hi*(1.0-1./3.*(ddelta_short_dt+ddelta_long_dt)/Hi); //dR/dt in Mpc/Mpc
			Rpi= -Ridot/((1.0+zi)*Hi); //dR/dz in Mpc.


      if(debug_mode > 0){
       printf("Ri=%.3le, dRi/dz=%.3le, delta_short=%.2le \n", Ri, Rpi, delta_short);
//					printf("i=%d, T_cdm_i=%.2le, T_g=%.2le, T_nu=%.2le \n",i_klong, Ti_klong[i_klong], transfer_gamma_klong[i_klong][33], transfer_nu_massless_klong[i_klong][33]);
     }


//we now call the collapse routine relevant for each case:

 				if(prf->Omega_extra + prf->Omega_nu2 + prf->Omega_nu1 == 0){//nothing extra
           if (debug_mode>0) printf("Collapse with nothing extra \n");
           z_coll_iteration=find_z_collapse_nothing(pba, prf, Ri, Rpi, delta_long, 
             Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], z_collapse, Mhalo_Mpc);
         }
				else if(prf->Omega_extra + prf->Omega_nu2 == 0){//only m_nu1
					if (debug_mode>0) printf("Collapse with nu1 \n");
					z_coll_iteration=find_z_collapse_1nu(pba,prf, Ri, Rpi, delta_long, 
           Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
           transfer_nu1_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS,
           Nz_solution, Rhalo_solution, Mnu1_solution, z_collapse, Mhalo_Mpc);
				}
				else if (prf->Omega_extra == 0 && prf->Omega_nu2 > 0){//m_nu1 and m_nu2
					if (debug_mode>0) printf("Collapse with nu1 and nu2 \n");
					z_coll_iteration=find_z_collapse_2nu(pba,prf, Ri, Rpi, delta_long, 
           Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],
           zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS,
           Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution, z_collapse, Mhalo_Mpc);

				}
				else if(prf->Omega_nu1 + prf->Omega_nu2 == 0){//only extra
					if (debug_mode>0) printf("Collapse with extra \n");
					z_coll_iteration=find_z_collapse_1nu(pba,prf, Ri, Rpi, delta_long, 
           Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
           transfer_extra_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_extra_EoS, plist_extra_EoS,
           Nz_solution, Rhalo_solution, Mextra_solution, z_collapse, Mhalo_Mpc);
				}
				else if (prf->Omega_extra > 0 && prf->Omega_nu2 == 0){//m_nu1 and extra, useful for CAMB
					if (debug_mode>0) printf("Collapse with nu1 and extra \n");
					z_coll_iteration=find_z_collapse_2nu(pba,prf, Ri, Rpi, delta_long,
           Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],
           transfer_extra_klong[i_klong],
           zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_extra_EoS, plist_extra_EoS,
           Nz_solution, Rhalo_solution, Mnu1_solution, Mextra_solution, z_collapse, Mhalo_Mpc);

				}
				else{	//m_nu1, m_nu2 and extra
					if (debug_mode>0) printf("Collapse with nu1,nu2 & extra \n");
					z_coll_iteration=find_z_collapse_3nu(pba,prf, Ri, Rpi, delta_long,
           Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
           transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],transfer_extra_klong[i_klong],
           zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS, rholist_extra_EoS, plist_extra_EoS,
           Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution, Mextra_solution, z_collapse, Mhalo_Mpc);

				}



				if(z_coll_iteration>z_collapse){ //collapses too quickly
					delta_short_max=delta_short;
				}
				else{ 							//collapses too slowly
					delta_short_min=delta_short;
				}






			}	//end of delta_short loop

			if(debug_mode > 0){
       printf("delta_short=%le, delta_long=%le (with k_long=%le) \n",delta_short, delta_long, k_long);
       printf("z_coll_iteration= %le \n",z_coll_iteration);
     }


			if(abs(z_coll_iteration-z_collapse)>tolerance_z){ //collapses too quickly
				printf("Not converged to z_collapse. Initial conditions too narrow, make boost_initial_conditions bigger. \n");
				printf("z_coll_it=%.3le and z_collapse=%.3le \n",z_coll_iteration,z_collapse);
			}

			prf->delta_short_collapse[z_ind][m_ind][i_delta_long][i_klong]=delta_short;


		}	// end of k_long loop
	}	//end of delta_long loop



	// since the delta_long=0 does not depend on k we can just copy it for all ks.
	for(i_delta_long=0;i_delta_long<1;i_delta_long++){
		for(i_klong=1;i_klong<prf->N_klong;i_klong++){
			prf->delta_short_collapse[z_ind][m_ind][i_delta_long][i_klong]=prf->delta_short_collapse[z_ind][m_ind][i_delta_long][0];
		}
	}




	if (debug_mode > 1){
		for (i=0;i<Nz_solution;i++){
			printf("i=%d, z=%.1le, R(z)/Mpc=%.1le \n",i, logz_solution_min*exp(dlogz_solution*i),Rhalo_solution[i]);
		}
	}




  if(do_clustering == 1 && (prf->Omega_extra + prf->Omega_nu2 + prf->Omega_nu1 > 0)){

   for(i_clustering_index=0;i_clustering_index<Nrecusions_nu_clustering;i_clustering_index++){


    printf("We add clustering of neutrinos/extra species: \n");

		//we compute the collapse of neutrinos/other species_counter.
		//we do it one at a time, since no more than one usually matters.
  int cluster_check=0;

    if(prf->Omega_nu1>0 && (prf->m_nu1 > mnu_min_clustering*sqrt(Mhalo_min_clustering/Mhalo)) ){
     cluster_check = findmass_1nu(pba, prf, prf->m_nu1, prf->T0_nu, zmin_EoS, dz_EoS, Nz_EoS,
       rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS, etalist,
       logz_solution_min, dlogz_solution, Nz_solution, Rhalo_solution, Mnu1_solution, z_collapse, Mhalo_Mpc);
   }
   if(prf->Omega_nu2>0 && (prf->m_nu2 > mnu_min_clustering*sqrt(Mhalo_min_clustering/Mhalo)) ){
     cluster_check = findmass_1nu(pba, prf, prf->m_nu2, prf->T0_nu, zmin_EoS, dz_EoS, Nz_EoS,
       rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS, etalist,
       logz_solution_min, dlogz_solution, Nz_solution, Rhalo_solution, Mnu2_solution, z_collapse, Mhalo_Mpc);
   }
   if(prf->Omega_extra>0 && (prf->m_extra > prf->T0_extra/prf->T0_nu * mnu_min_clustering*sqrt(Mhalo_min_clustering/Mhalo))){
     cluster_check = findmass_1nu(pba, prf, prf->m_extra, prf->T0_extra, zmin_EoS, dz_EoS, Nz_EoS,
       rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS, etalist,
       logz_solution_min, dlogz_solution, Nz_solution, Rhalo_solution, Mextra_solution, z_collapse, Mhalo_Mpc);
   }

   if (cluster_check){printf("Problem with Nu clustering in RelicFast\n");}


		///and we do the same thing now that we have calculated the nu_collapse, but with more precision.

		#pragma omp parallel for
   for(i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			#pragma omp parallel for
     for(i_klong=0;i_klong<N_klong_calc[i_delta_long];i_klong++){

      double delta_long=prf->delta_long_list[i_delta_long];


      double ddelta_long_dt=delta_long*(-Hi*(1.+zi))*dTi_klong[i_klong]/Ti_klong[i_klong];

				double delta_short_min=delta_short_min_constant; //we restart the initial conditions, and redo the procedure.
				double delta_short_max=delta_short_max_constant;
				double delta_short,ddelta_short_dt; //private copies for each run
				double Ri, Ridot, Rpi, z_coll_iteration; //initial conditions and z of collapse of each iteration.


				for(;delta_short_max/delta_short_min-1.>tolerance;){


					delta_short=(delta_short_max+delta_short_min)/2.; //bisection, we update the min or max value.

					if(option_initial_velocity==0){
				///Option 1: through sigma(M)
						ddelta_short_dt=delta_short*(-Hi*(1.+zi))*d_sigmaM_i/sigmaM_i; //since: 		ddelta_dt/delta = dsigma_dt/sigma = -H*(1+z) dsigma_dz/sigma
					}
					else if(option_initial_velocity==1) {
				///Option 2: through transfer function (they agree)
						ddelta_short_dt=delta_short*(-Hi*(1.+zi))*dTi_kstar/Ti_kstar; //since:		ddelta_dt/delta = dT_dt/T = -H*(1+z) dT_dz/sigma
					}
					else {
						printf("option_initial_velocity has to be either 0 or 1 \n");
						exit(0);
					}


				//these are the i.c. we will feed the function find_z_collapse_long.
				Ri= Ribar*(1.0-1./3.*(delta_short+delta_long)); //R_i in Mpc.
				Ridot= Ri*Hi*(1.0-1./3.*(ddelta_short_dt+ddelta_long_dt)/Hi); //dR/dt in Mpc/Mpc
				Rpi= -Ridot/((1.0+zi)*Hi); //dR/dz in Mpc.



				if(debug_mode > 0){
					printf("Ri=%.3le, dRi/dz=%.3le, delta_short=%.2le \n", Ri, Rpi, delta_short);
//					printf("i=%d, T_cdm_i=%.2le, T_g=%.2le, T_nu=%.2le \n",i_klong, Ti_klong[i_klong], transfer_gamma_klong[i_klong][33], transfer_nu_massless_klong[i_klong][33]);
				}


		//we now call the collapse routine relevant for each case:

		 				if(prf->Omega_extra + prf->Omega_nu2 + prf->Omega_nu1 == 0){//nothing extra
             if (debug_mode>0) printf("Collapse with nothing extra \n");
             z_coll_iteration=find_z_collapse_nothing(pba, prf, Ri, Rpi, delta_long, 
               Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], z_collapse, Mhalo_Mpc);
           }
						else if(prf->Omega_extra + prf->Omega_nu2 == 0){//only m_nu1
							if (debug_mode>0) printf("Collapse with nu1 \n");
							z_coll_iteration=find_z_collapse_1nu(pba,prf, Ri, Rpi, delta_long, 
               Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
               transfer_nu1_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS,
               Nz_solution, Rhalo_solution, Mnu1_solution, z_collapse, Mhalo_Mpc);
						}
						else if (prf->Omega_extra == 0 && prf->Omega_nu2 > 0){//m_nu1 and m_nu2
							if (debug_mode>0) printf("Collapse with nu1 and nu2 \n");
							z_coll_iteration=find_z_collapse_2nu(pba,prf, Ri, Rpi, delta_long, 
               Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],
               zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS,
               Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution, z_collapse, Mhalo_Mpc);

						}
						else if(prf->Omega_nu1 + prf->Omega_nu2 == 0){//only extra
							if (debug_mode>0) printf("Collapse with extra \n");
							z_coll_iteration=find_z_collapse_1nu(pba,prf, Ri, Rpi, delta_long, 
               Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
               transfer_extra_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_extra_EoS, plist_extra_EoS,
               Nz_solution, Rhalo_solution, Mextra_solution, z_collapse, Mhalo_Mpc);
						}
						else if (prf->Omega_extra > 0 && prf->Omega_nu2 == 0){//m_nu1 and extra, useful for CAMB
							if (debug_mode>0) printf("Collapse with nu1 and extra \n");
							z_coll_iteration=find_z_collapse_2nu(pba,prf, Ri, Rpi, delta_long, 
               Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],
               transfer_extra_klong[i_klong],
               zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_extra_EoS, plist_extra_EoS,
               Nz_solution, Rhalo_solution, Mnu1_solution, Mextra_solution, z_collapse, Mhalo_Mpc);

						}
						else{	//m_nu1, m_nu2 and extra
							if (debug_mode>0) printf("Collapse with nu1,nu2 & extra \n");
							z_coll_iteration=find_z_collapse_3nu(pba,prf, Ri, Rpi, delta_long, 
               Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
               transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],transfer_extra_klong[i_klong],
               zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS, rholist_extra_EoS, plist_extra_EoS,
               Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution, Mextra_solution, z_collapse, Mhalo_Mpc);

						}



						if(z_coll_iteration>z_collapse){ //collapses too quickly
							delta_short_max=delta_short;
						}
						else{ 							//collapses too slowly
							delta_short_min=delta_short;
						}






					}	//end of delta_short loop

				//BLOOP



				if(abs(z_coll_iteration-z_collapse)>tolerance_z){ //collapses too quickly
					printf("Not converged to z_collapse (after doing Mnu_collapse). Initial conditions too narrow, make delta_short_max and delta_short_min wider. \n");
					printf("z_coll_it=%.3le and z_collapse=%.3le \n",z_coll_iteration,z_collapse);
				}

				prf->delta_short_collapse[z_ind][m_ind][i_delta_long][i_klong]=delta_short;


			}	// end of k_long loop
		}	//end of delta_long loop

		printf("\n Iteration of clustering #%d \n",i_clustering_index+1);

	} //end of do_clustering loop

} //end of do_clustering if statement.



// since the delta_long=0 does not depend on k we can just copy it for all ks (again).
for( i_delta_long=0;i_delta_long<1;i_delta_long++){
	for(i_klong=1;i_klong<prf->N_klong;i_klong++){
		prf->delta_short_collapse[z_ind][m_ind][i_delta_long][i_klong]=prf->delta_short_collapse[z_ind][m_ind][i_delta_long][0];
	}
}






//////////////////////////////////////////////////////////////////////////////
//// we now extrapolate  to find\delta_crit to the redshift of collapse	/////
////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for collapse(2)
for(i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
  for(i_klong=0;i_klong<prf->N_klong;i_klong++){
   double delta_long=prf->delta_long_list[i_delta_long];
   double k_long=prf->klong_list[i_klong];

   T_z_collapse_klong=interpol_2D(transfer_matter, prf->zlist_transfer, prf->Nz_transfer, k_transfer_array, prf->length_transfer, z_collapse, k_long);

			prf->delta_long_collapse[z_ind][m_ind][i_delta_long][i_klong]=delta_long*T_z_collapse_klong/Ti_klong[i_klong];	//this is right as well //			printf("klong=%.1le, prf->delta_long_collapse[z_ind][m_ind]/delta_long= %.3le \n", k_long, prf->delta_long_collapse[z_ind][m_ind][i_delta_long][i_klong]/delta_long);


			if(option_initial_velocity==0){
				prf->delta_short_crit[z_ind][m_ind][i_delta_long][i_klong]=prf->delta_short_collapse[z_ind][m_ind][i_delta_long][i_klong]*sigmaM_collapse/sigmaM_i;
			}
			else if(option_initial_velocity==1){
				prf->delta_short_crit[z_ind][m_ind][i_delta_long][i_klong]=prf->delta_short_collapse[z_ind][m_ind][i_delta_long][i_klong]*T_z_collapse_kstar/Ti_kstar;
			}
			else{
				printf("Error, select an option_initial_velocity \n");
			}
		}
	}

} //end of Mhalo loop.










//////////////////////////////////////////
//// we free the allocated memory	/////
////////////////////////////////////////

free_array_2D(transfer_array_z0,6);
free(k_transfer_array);


free_array_2D(transfer_array_zi,6);
free_array_2D(transfer_array_zip1,6);
free_array_2D(transfer_array_zip2,6);

free_array_2D(transfer_array_z_collapse,6);


free(d_transfer_array_zi);


free_array_2D(transfer_matter, prf->Nz_transfer);
free_array_2D(transfer_gamma, prf->Nz_transfer);
free_array_2D(transfer_nu_massless, prf->Nz_transfer);
free_array_2D(transfer_nu1, prf->Nz_transfer);
free_array_2D(transfer_nu2, prf->Nz_transfer);
free_array_2D(transfer_extra, prf->Nz_transfer);


free_array_2D(transfer_gamma_klong, prf->N_klong);
free_array_2D(transfer_nu_massless_klong, prf->N_klong);
free_array_2D(transfer_nu1_klong, prf->N_klong);
free_array_2D(transfer_nu2_klong, prf->N_klong);
free_array_2D(transfer_extra_klong, prf->N_klong);


free(Ti_klong);
free(dTi_klong);



free(zlist_EoS);
free(plist_nu1_EoS);
free(rholist_nu1_EoS);

free(plist_nu2_EoS);
free(rholist_nu2_EoS);

free(plist_extra_EoS);
free(rholist_extra_EoS);

free(etalist);



free(Rhalo_solution);
free(Mnu1_solution);
free(Mnu2_solution);
free(Mextra_solution);



return _SUCCESS_;
}






//now collapse functions:

double find_z_collapse_nothing
(struct background * pba, struct relicfast * prf, const double Ri, const double Rpi, const double delta_long, const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong, double z_collapse, double Mhalo_Mpc
  ){
// This function returns the z at which an overdensity collapses. global prflogical parameters assumed.
// In the presence of a long-wavelength PHOTON and MASSLESS NU perturbation.
//we read the instantaneous transfer function from CLASS/CAMB at each redshift.

	int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20


	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
  }

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;
	double dE, OmL, OmR, H, OmM; //z-dependent.
	double H2; //H for previous redshift, for derivatives.

	double OmG, Omnu_massless; //photon, massless and massive neutrino energy density.
	double OmRbar, OmGbar, Omnu_masslessbar; //average for all of them



	double T_gamma, T_nu;//transfer functions at each k and z.

	const double T_matter=Tfm_klong;



//we set the initial H
	OmGbar= pba->Omega0_g * pow(1.+zi,4.);
	Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+zi,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = prf->Omega_M * pow(1.+zi,3.);
	OmL = prf->Omega_L; //Omega_L(z), we take it as z-independent.
	H = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar); //H(zi)


//and the initial OmR
	T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	OmG= OmGbar * (1.+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	Omnu_massless= Omnu_masslessbar * (1.+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR2 = OmG + Omnu_massless;




	double z_next, dE2; //for the next values, to do Heun's method. We assume Omega_L does not change.

//symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	z_next = zi * exp(zstep_log);//first step
	OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
	Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = prf->Omega_M * pow(1.+z_next,3.);
	OmL = prf->Omega_L;
	H2  = prf->H0_Mpc * sqrt(OmL + OmM + OmRbar);// H(z_next)
	dE2 = (H2-H)/H2/(zi*zstep_log);

	double Rpp1, Rpp2; //d^2R(z)/dz^2


///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////
	for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
   R1 = R2;
   Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
   z=z_next;
   H=H2;
   OmR=OmR2;
   dE=dE2;

		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
		Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
		OmRbar = OmGbar + Omnu_masslessbar;
		OmM = prf->Omega_M * pow(1.+z_next,3.);


		T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);

		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
		OmR2 = OmG + Omnu_massless;


		H2  = prf->H0_Mpc * sqrt(OmL + OmM + OmRbar);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value at the end.



//these are too annoying to keep even with debug_mode, activate manually if you want to debug:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, prf->H0_Mpc*E_LCDM(prf, z), dE, dlogE_dz_LCDM(prf, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)

//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}

	return z;

}







double find_z_collapse_1nu
(struct background * pba, struct relicfast * prf, const double Ri, const double Rpi, const double delta_long, 
  const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong, double * const transfer_nu_massive_klong,
  const double zmin_EoS, const double dz_EoS, const long Nz_EoS, double * const rholist_EoS, double * const plist_EoS,
  const long Nz_solution, double * R_solution, double * Mnu_solution, double z_collapse, double Mhalo_Mpc
  ){
// This function returns the z at which an overdensity collapses. global prflogical parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 1 MASSIVE NU (OR OTHER EXTRA) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.
//we only modify R_solution and Mnu_solution, the rest are read-only.


	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20

	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
  }

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;



	const double T_matter=Tfm_klong;

	//the nu_massive can be any new species, just make sure that the Omegas are well defined here
  if((prf->Omega_extra>0) && (prf->m_nu1>0)){
    printf("Using wrong collapse function, XXX_1nu only accepts one species. \n");
  }
	double rhonu1_0=rholist_EoS[0];//rho at z=0. units are K^4




//we set the initial H. All these are average densities.
	double OmGbar= pba->Omega0_g * pow(1.+zi,4.); //photon
	double Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+zi,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = prf->Omega_M * pow(1.+zi,3.); //matter
	double OmL = prf->Omega_L; //Omega_Lambda(z), we take it as z-independent.
	double rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu1=rhonu1_z/rhonu1_0; //ratio to rho of zero
	double Omnu1bar_z =  prf->Omega_nu1 * rho_ratio_nu1 + prf->Omega_extra * rho_ratio_nu1; //Omega_nu1
	double H = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z); //H(zi)



//and the initial OmR (we assme c_s^2=w=1/3 for photons and massless nu)
	double T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//and initial perturbation in nu1	(massive) and pressure
	double pnu1_z = interpol_cubic(zmin_EoS, dz_EoS, plist_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu1_z = pnu1_z/rhonu1_z; //equation of state
	double T_nu1_z=interpol(transfer_nu_massive_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next





//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = zi * exp(zstep_log);//first step
	OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
	Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = prf->Omega_M * pow(1.+z_next,3.);
	OmL = prf->Omega_L;
	rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_EoS, Nz_EoS, z_next);
	rho_ratio_nu1=rhonu1_z/rhonu1_0;
	Omnu1bar_z =  prf->Omega_nu1 * rho_ratio_nu1 + prf->Omega_extra * rho_ratio_nu1;
	double H2 = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z); // H(z_next)
	double dE, dE2 = (H2-H)/H2/(zi*zstep_log);

//and the w and c_s^2 of nu at first step we calculate.
	pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist_EoS, Nz_EoS, z_next);
	double wnu1_z2=pnu1_z/rhonu1_z;
	double d_wnu1_z = (wnu1_z2-wnu1_z)/(z_next-zi);
	double csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/(3.0*(1+wnu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).


	double Rpp1, Rpp2; //d^2R(z)/dz^2

	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double Mhalo_nu1_Mpc=0; //M of the nu1 halo
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
  R_solution[i_solution]=R1;
  Mhalo_nu1_Mpc = Mnu_solution[i_solution];
  i_solution++;


///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////

  for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
   R1 = R2;
   Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
   z=z_next;
   H=H2;
   dE=dE2;
   OmR=OmR2;
   wnu1_z=wnu1_z2;

//we save R and read Mnu, only every few steps since it doesn't vary much.
   if (i % collapse_steps == 0){
    R_solution[i_solution]=R1;
    Mhalo_nu1_Mpc = Mnu_solution[i_solution];
    i_solution++;
    if(debug_mode >1){
      printf("z=%.1le, R=%.1le, Mnu/Mhalo=%1le \n\n", z, R_solution[i_solution], Mhalo_nu1_Mpc/Mhalo_Mpc);
    }
  }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + (1.+3.*wnu1_z + delta_nu1_z*(1+3.0*csq_ad_nu1_z))*Omnu1bar_z - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
		T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)

		Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
		T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)

		OmRbar = OmGbar + Omnu_masslessbar;
		OmR2 = OmG + Omnu_massless;

		OmM = prf->Omega_M * pow(1.+z_next,3.);



		rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_EoS, Nz_EoS, z_next);
		rho_ratio_nu1=rhonu1_z/rhonu1_0;
		Omnu1bar_z = prf->Omega_nu1 * rho_ratio_nu1 + prf->Omega_extra * rho_ratio_nu1;
		T_nu1_z=interpol(transfer_nu_massive_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
		pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist_EoS, Nz_EoS, z_next);
		wnu1_z2=pnu1_z/rhonu1_z;
		d_wnu1_z = (wnu1_z2-wnu1_z)/zstep_lin;
		csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/3.0*(1.0+z_next)/(1.0+wnu1_z2); //adiabatic sound speed squared.


		H2  = prf->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1bar_z);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value, although it's unnecessary.



		if(debug_mode > 1){ //check that derivatives are well defined, it should be for any reasonable value of npoints, but just in case:
			printf("w(z1)-w(z2)=%.3le , w(z2)=%.3le , dw/dz=%.3le \n", wnu1_z-wnu1_z2, wnu1_z2, d_wnu1_z);
			printf("H(z1)-H(z2)=%.3le , H(z2)=%.3le , dH/dz=%.3le \n", H-H2, H2, dE2*H2);
			printf("z=%.3le, cs^2/w(z2)=%.3le \n", z_next, csq_ad_nu1_z/wnu1_z2); //cs^2 and w are within 10% of each other.
   }


//these are too annoying to keep even with debug_mode, activate manually if you want:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, prf->H0_Mpc*E_LCDM(prf, z), dE, dlogE_dz_LCDM(prf, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

   R1tilde = R1 + zstep_lin * Rp1;
   Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
   Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
   - Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
   - Mhalo_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 + (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)


//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}



  return z;

}





double find_z_collapse_2nu
(struct background * pba, struct relicfast * prf, const double Ri, const double Rpi, const double delta_long, 
  const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong, double * const transfer_nu1_klong, double * const transfer_nu2_klong,
  const double zmin_EoS, const double dz_EoS, const long Nz_EoS, double * const rholist1_EoS, double * const plist1_EoS, double * const rholist2_EoS, double * const plist2_EoS,
  const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution, double z_collapse, double Mhalo_Mpc
  ){
// This function returns the z at which an overdensity collapses. global prflogical parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 2 MASSIVE NU (or 1 nu + 1 extra) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.




	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20

	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
  }

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;



	const double T_matter=Tfm_klong;

	//the nu1 is nu1, but nu2 can be any new species, just make sure that the Omegas are well defined here
  if((prf->Omega_extra>0) && (prf->m_nu2>0)){
    printf("Using wrong collapse function, XXX_2nu only accepts nu1 + one other species. \n");
  }
	double rhonu1_0=rholist1_EoS[0];//rho1 at z=0. units are K^4
  double rhonu2_0=rholist2_EoS[0];//rho2 at z=0. units are K^4




//we set the initial H. All these are average densities.
	double OmGbar= pba->Omega0_g * pow(1.+zi,4.); //photon
	double Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+zi,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = prf->Omega_M * pow(1.+zi,3.); //matter
	double OmL = prf->Omega_L; //Omega_Lambda(z), we take it as z-independent.

	double rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu1=rhonu1_z/rhonu1_0; //ratio to rho of zero
	double Omnu1bar_z =  prf->Omega_nu1 * rho_ratio_nu1; //Omega_nu1

	double rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu2=rhonu2_z/rhonu2_0; //ratio to rho of zero
	double Omnu2bar_z =  prf->Omega_nu2 * rho_ratio_nu2 + prf->Omega_extra * rho_ratio_nu2; //Omega_nu2 or omega extra, only one!

	double H = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z); //H(zi)


//and the initial OmR (we assme c_s^2=w=1/3 for photons and massless nu)
	double T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double OmG= OmGbar * (1.+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//and initial perturbation in nu1	(massive) and pressure
	double pnu1_z = interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu1_z = pnu1_z/rhonu1_z; //equation of state
	double T_nu1_z=interpol(transfer_nu1_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//and initial perturbation in nu2	(massive) and pressure
	double pnu2_z = interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu2_z = pnu2_z/rhonu2_z; //equation of state
	double T_nu2_z=interpol(transfer_nu2_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double delta_nu2_z = delta_long*T_nu2_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next



//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = zi * exp(zstep_log);//first step
	OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
	Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = prf->Omega_M * pow(1.+z_next,3.);
	OmL = prf->Omega_L;

	rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
	rho_ratio_nu1=rhonu1_z/rhonu1_0;
	Omnu1bar_z =  prf->Omega_nu1 * rho_ratio_nu1;

	rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
	rho_ratio_nu2=rhonu2_z/rhonu2_0;
	Omnu2bar_z =  prf->Omega_nu2 * rho_ratio_nu2 + prf->Omega_extra * rho_ratio_nu2;

	double H2 = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z); // H(z_next)
	double dE, dE2 = (H2-H)/H2/(zi*zstep_log);

//and the w and c_s^2 of nus at first step we calculate.
//for nu1
  pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
  double wnu1_z2=pnu1_z/rhonu1_z;
  double d_wnu1_z = (wnu1_z2-wnu1_z)/(z_next-zi);
		double csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/(3.0*(1+wnu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu2
		pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
		double wnu2_z2=pnu2_z/rhonu2_z;
		double d_wnu2_z = (wnu2_z2-wnu2_z)/(z_next-zi);
		double csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/(3.0*(1+wnu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).


	double Rpp1, Rpp2; //d^2R(z)/dz^2

	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double Mhalo_nu1_Mpc=0, Mhalo_nu2_Mpc=0; //M of the nu1 and nu2 haloes
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
  R_solution[i_solution]=R1;
  Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
  Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
  i_solution++;


///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////
  for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
   R1 = R2;
   Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
   z=z_next;
   H=H2;
   dE=dE2;
   OmR=OmR2;
   wnu1_z=wnu1_z2;
   wnu2_z=wnu2_z2;


//we save R and read Mnu, only every few steps since it doesn't vary much.
   if (i % collapse_steps == 0){
    R_solution[i_solution]=R1;
    Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
    Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
    i_solution++;
    if(debug_mode > 1){
      printf("z=%.1le, R=%.1le, m_nu1=%1le , m_nu2=%1le \n\n", z, R_solution[i_solution], Mhalo_nu1_Mpc, Mhalo_nu2_Mpc);
    }
  }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.


//these are the terms that go inside the sum, define them outside for clarity.
		double Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		double Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu2_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + Onu1 + Onu2 - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
		T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)

		Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
		T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)

		OmRbar = OmGbar + Omnu_masslessbar;
		OmR2 = OmG + Omnu_massless;

		OmM = prf->Omega_M * pow(1.+z_next,3.);


		rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
		rho_ratio_nu1=rhonu1_z/rhonu1_0;
		Omnu1bar_z = prf->Omega_nu1 * rho_ratio_nu1;
		T_nu1_z=interpol(transfer_nu1_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
		pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
		wnu1_z2=pnu1_z/rhonu1_z;
		d_wnu1_z = (wnu1_z2-wnu1_z)/zstep_lin;
		csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/3.0*(1.0+z_next)/(1.0+wnu1_z2); //adiabatic sound speed squared.

		rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
		rho_ratio_nu2=rhonu2_z/rhonu2_0;
		Omnu2bar_z = prf->Omega_nu2 * rho_ratio_nu2 + prf->Omega_extra * rho_ratio_nu2;
		T_nu2_z=interpol(transfer_nu2_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		delta_nu2_z = delta_long*T_nu2_z/T_matter; //long-wavlength neutrino perturbation
		pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
		wnu2_z2=pnu2_z/rhonu2_z;
		d_wnu2_z = (wnu2_z2-wnu2_z)/zstep_lin;
		csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/3.0*(1.0+z_next)/(1.0+wnu2_z2); //adiabatic sound speed squared.

//we update the Onu terms that go inside the integral.
		Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;

		H2  = prf->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value at the end.



//these are too annoying to keep even with debug_mode, activate manually if you want:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, prf->H0_Mpc*E_LCDM(prf, z), dE, dlogE_dz_LCDM(prf, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu2_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 + Onu1 + Onu2 - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)


//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}



  return z;

}



double find_z_collapse_3nu
(struct background * pba, struct relicfast * prf, const double Ri, const double Rpi, const double delta_long, 
  const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong, double * const transfer_nu1_klong, double * const transfer_nu2_klong, double * const transfer_nu3_klong,
  const double zmin_EoS, const double dz_EoS, const long Nz_EoS, double * const rholist1_EoS, double * const plist1_EoS, double * const rholist2_EoS, double * const plist2_EoS, double * const rholist3_EoS, double * const plist3_EoS,
  const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution, double *Mnu3_solution, double z_collapse, double Mhalo_Mpc
  ){
// This function returns the z at which an overdensity collapses. global prflogical parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 2 MASSIVE NU perturbation, and an EXTRA SPECIES.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.





	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20

	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
  }

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;



	const double T_matter=Tfm_klong;

	//the nu1 is nu1, but nu2 can be any new species, just make sure that the Omegas are well defined here
  if((prf->Omega_extra==0) || (prf->m_nu2==0) || (prf->m_nu1==0)){
    printf("Using wrong collapse function, XXX_3nu works with 3 nus (or 2+extra). \n");
  }
	double rhonu1_0=rholist1_EoS[0];//rho1 at z=0. units are K^4
	double rhonu2_0=rholist2_EoS[0];//rho2 at z=0. units are K^4
	double rhonu3_0=rholist3_EoS[0];//rho3 at z=0. units are K^4





//we set the initial H. All these are average densities.
	double OmGbar= pba->Omega0_g * pow(1.+zi,4.); //photon
	double Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+zi,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = prf->Omega_M * pow(1.+zi,3.); //matter
	double OmL = prf->Omega_L; //Omega_Lambda(z), we take it as z-independent.

	double rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu1=rhonu1_z/rhonu1_0; //ratio to rho of zero
	double Omnu1bar_z =  prf->Omega_nu1 * rho_ratio_nu1; //Omega_nu1

	double rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, zi); //rho of nu2 at z
	double rho_ratio_nu2=rhonu2_z/rhonu2_0; //ratio to rho of zero
	double Omnu2bar_z =  prf->Omega_nu2 * rho_ratio_nu2; //Omega_nu2

	double rhonu3_z=interpol_cubic(zmin_EoS, dz_EoS, rholist3_EoS, Nz_EoS, zi); //rho of nu3 at z
	double rho_ratio_nu3=rhonu3_z/rhonu3_0; //ratio to rho of zero
	double Omnu3bar_z = prf->Omega_extra * rho_ratio_nu3; //omeganu3

	double H = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z); //H(zi)


//and the initial OmR
	double T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double OmG= OmGbar * (1.+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//and initial perturbation in nu1	(massive) and pressure
	double pnu1_z = interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu1_z = pnu1_z/rhonu1_z; //equation of state
	double T_nu1_z=interpol(transfer_nu1_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//and initial perturbation in nu2	(massive) and pressure
	double pnu2_z = interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu2_z = pnu2_z/rhonu2_z; //equation of state
	double T_nu2_z=interpol(transfer_nu2_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double delta_nu2_z = delta_long*T_nu2_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//and initial perturbation in nu3	(massive) and pressure
	double pnu3_z = interpol_cubic(zmin_EoS, dz_EoS, plist3_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu3_z = pnu3_z/rhonu3_z; //equation of state
	double T_nu3_z=interpol(transfer_nu3_klong, prf->zlist_transfer, prf->Nz_transfer, zi);
	double delta_nu3_z = delta_long*T_nu3_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next


//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = zi * exp(zstep_log);//first step
	OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
	Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = prf->Omega_M * pow(1.+z_next,3.);
	OmL = prf->Omega_L;

	rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
	rho_ratio_nu1=rhonu1_z/rhonu1_0;
	Omnu1bar_z =  prf->Omega_nu1 * rho_ratio_nu1;

	rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
	rho_ratio_nu2=rhonu2_z/rhonu2_0;
	Omnu2bar_z =  prf->Omega_nu2 * rho_ratio_nu2;

	rhonu3_z=interpol_cubic(zmin_EoS, dz_EoS, rholist3_EoS, Nz_EoS, z_next);
	rho_ratio_nu3=rhonu3_z/rhonu3_0;
	Omnu3bar_z = prf->Omega_extra * rho_ratio_nu3;

	double H2 = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z); // H(z_next)
	double dE, dE2 = (H2-H)/H2/(zi*zstep_log);

//and the w and c_s^2 of nus at first step we calculate.
//for nu1
	pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
	double wnu1_z2=pnu1_z/rhonu1_z;
	double d_wnu1_z = (wnu1_z2-wnu1_z)/(z_next-zi);
	double csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/(3.0*(1+wnu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu2
	pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
	double wnu2_z2=pnu2_z/rhonu2_z;
	double d_wnu2_z = (wnu2_z2-wnu2_z)/(z_next-zi);
	double csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/(3.0*(1+wnu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu3
	pnu3_z=interpol_cubic(zmin_EoS, dz_EoS, plist3_EoS, Nz_EoS, z_next);
	double wnu3_z2=pnu3_z/rhonu3_z;
	double d_wnu3_z = (wnu3_z2-wnu3_z)/(z_next-zi);
	double csq_ad_nu3_z = wnu3_z2 + d_wnu3_z/(3.0*(1+wnu3_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).






	double Rpp1, Rpp2; //d^2R(z)/dz^2

	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double Mhalo_nu1_Mpc=0, Mhalo_nu2_Mpc=0, Mhalo_nu3_Mpc=0; //M of the nu1, nu2, and nu3 haloes
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
  R_solution[i_solution]=R1;
  Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
  Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
  Mhalo_nu3_Mpc = Mnu3_solution[i_solution];
  i_solution++;




///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////
  for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
   R1 = R2;
   Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
   z=z_next;
   H=H2;
   dE=dE2;
   OmR=OmR2;
   wnu1_z=wnu1_z2;
   wnu2_z=wnu2_z2;
   wnu3_z=wnu3_z2;

//we save R and read Mnu, only every few steps since it doesn't vary much.
   if (i % collapse_steps == 0){
    R_solution[i_solution]=R1;
    Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
    Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
    Mhalo_nu3_Mpc = Mnu3_solution[i_solution];
    i_solution++;
    if(debug_mode > 1){
      printf("z=%.1le, R=%.1le, m_nu1=%1le , m_nu2=%1le , Mnu3=%1le \n\n", z, R_solution[i_solution], Mhalo_nu1_Mpc, Mhalo_nu2_Mpc, Mhalo_nu3_Mpc);
    }
  }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

//these are the terms that go inside the sum, define them outside for clarity.
		double Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		double Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
		double Onu3 = (1.0+3.0*wnu3_z + delta_nu3_z*(1.0+3.0*csq_ad_nu3_z))*Omnu3bar_z;

	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu2_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu3_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + Onu1 + Onu2 + Onu3 - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= pba->Omega0_g * pow(1.+z_next,4.);
		T_gamma = interpol(transfer_gamma_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)

		Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_next,4.);
		T_nu = interpol(transfer_nu_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)

		OmRbar = OmGbar + Omnu_masslessbar;
		OmR2 = OmG + Omnu_massless;

		OmM = prf->Omega_M * pow(1.+z_next,3.);


		rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
		rho_ratio_nu1=rhonu1_z/rhonu1_0;
		Omnu1bar_z = prf->Omega_nu1 * rho_ratio_nu1;
		T_nu1_z=interpol(transfer_nu1_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
		pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
		wnu1_z2=pnu1_z/rhonu1_z;
		d_wnu1_z = (wnu1_z2-wnu1_z)/zstep_lin;
		csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/3.0*(1.0+z_next)/(1.0+wnu1_z2); //adiabatic sound speed squared.

		rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
		rho_ratio_nu2=rhonu2_z/rhonu2_0;
		Omnu2bar_z = prf->Omega_nu2 * rho_ratio_nu2;
		T_nu2_z=interpol(transfer_nu2_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		delta_nu2_z = delta_long*T_nu2_z/T_matter; //long-wavlength neutrino perturbation
		pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
		wnu2_z2=pnu2_z/rhonu2_z;
		d_wnu2_z = (wnu2_z2-wnu2_z)/zstep_lin;
		csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/3.0*(1.0+z_next)/(1.0+wnu2_z2); //adiabatic sound speed squared.

		rhonu3_z=interpol_cubic(zmin_EoS, dz_EoS, rholist3_EoS, Nz_EoS, z_next);
		rho_ratio_nu3=rhonu3_z/rhonu3_0;
		Omnu3bar_z = prf->Omega_extra * rho_ratio_nu3;
		T_nu3_z=interpol(transfer_nu3_klong, prf->zlist_transfer, prf->Nz_transfer, z_next);
		delta_nu3_z = delta_long*T_nu3_z/T_matter; //long-wavlength neutrino perturbation
		pnu3_z=interpol_cubic(zmin_EoS, dz_EoS, plist3_EoS, Nz_EoS, z_next);
		wnu3_z2=pnu3_z/rhonu3_z;
		d_wnu3_z = (wnu3_z2-wnu3_z)/zstep_lin;
		csq_ad_nu3_z = wnu3_z2 + d_wnu3_z/3.0*(1.0+z_next)/(1.0+wnu3_z2); //adiabatic sound speed squared.

//we update the Onu terms that go inside the integral.
		Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
		Onu3 = (1.0+3.0*wnu3_z + delta_nu3_z*(1.0+3.0*csq_ad_nu3_z))*Omnu3bar_z;


		H2  = prf->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value at the end.




//these are too annoying to keep even with debug_mode, activate manually if you want:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, prf->H0_Mpc*E_LCDM(prf, z), dE, dlogE_dz_LCDM(prf, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu2_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu3_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 + Onu1 + Onu2 + Onu3 - 2*OmL)*prf->H0_Mpc*prf->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)


//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}


  return z;

}






int findmass_1nu (struct background *pba, struct relicfast *prf, double mass_0, double temp_0, double zmin_EoS, double dz_EoS, long Nz_EoS,
	double* rho1list_EoS, double* rho2list_EoS, double* rho3list_EoS, double* etalist,
	double logz_solution_min, double dlogz_solution, long Nz_solution, double* R_solution, double* Mnu_solution, double z_collapse, double Mhalo_Mpc){
//this function finds Mnu(t) given the R(t) input. We use the BKT approximation for neutrino clustering.
//mass_0 and temp_0 are the mass and temperature (at z=0) for the particle (neutrino or extra), we specify it to distinguish which case we do.


  double precision_findmass = fmin(fmax(precision_clustering,0.5),2.0); //precision parameter to regulate speed vs. precision.

	long nzpoints=400*precision_findmass;//how many points we integrate over in z (or t)
	double zstep;//= (z-zi)/(nzpoints-1.); //linear binning in z, since we're interested in early times

  int Nx = (int) 12*precision_findmass; //how many positions
	int Nq = 10; //how many momenta, 10 is good.
	int Nmu = (int) 8*precision_findmass; //how many x \cdot p cosines.

	double ***f1;	//this is f1(x,p,mu) at different times
	class_alloc(f1,Nx*sizeof(double**),prf->error_message);

	int al_ind, al_ind1;
 for(al_ind=0;al_ind<Nx;al_ind++){
  class_alloc(f1[al_ind],Nq*sizeof(double*),prf->error_message);	
  for(al_ind1=0;al_ind1<Nq;al_ind1++){
   class_alloc(f1[al_ind][al_ind1],Nmu*sizeof(double),prf->error_message);

 }	
}


	double OmRbar, OmM, OmGbar, Omnu_masslessbar, Omnu1, Omnu2, Omnu3; //average for all of them, nu3=extra.
  double OmL=prf->Omega_L;



	double rho_ratio_nu1; //ratio of rho(z)/rho(z=0)
	double rho1_0=rho1list_EoS[0] + (prf->Omega_nu1 == 0) * 1;//rho at z=0. units are K^4. add 1 if not present to avoid NaNs
	double rho1_z;

	double rho_ratio_nu2; //ratio of rho(z)/rho(z=0)
	double rho2_0=rho2list_EoS[0] + (prf->Omega_nu2 == 0) * 1;//rho at z=0. units are K^4
	double rho2_z;

	double rho_ratio_nu3; //ratio of rho(z)/rho(z=0)
	double rho3_0=rho3list_EoS[0] + (prf->Omega_extra == 0) * 1;//rho at z=0. units are K^4
	double rho3_z;


  //all temperatures, mases, and momenta in eV unless otherwise specified (we'd say: _Mpc):
	double temp = temp_0*KtoeV;//at z=0 in eV, we will use comoving momenta
  double mass = mass_0;//
	//we do each species clustering every time.




  long i, iq, ix, imu, iz;



 // int lengthname=200;
 //	char *filename; //To open files.
 //	filename=(char *)malloc((lengthname+1)*sizeof(char));
//	FILE *fp; //to save Mnu(z) in case
//	if(debug_mode>0){
//	  lengthname=sprintf(filename,"tests/Mnu_z-%d.dat",prf->file_tag); //we save escape fraction in a file
//	  fp=fopen(filename,"w");
//	}


//we first find \eta(z), the superconformal time, \eta(z) = \int dz/(H(z)*(1+z))

  double H;


  double Mnusmooth;

	double eta_in, eta_out; //eta(t) and eta(t'), we integrate over t'.
	double z_in, z_out; //z insider integral and outside


  double factor = 2.0*mass/temp; //x, q, and z-independent factor, 2 because 2 degrees of freedom in f0. unitless.
  double factorx; //1/x^2
  double factorq; //df0/dq = e^q/T/(1+e^q/T)^2

  double x, q, mu;
  double parenthesis; //(x-xc)^2 (vectors).
  double R_halo, x_halo; //Rhalo(t) and xhalo=Rhalo(t)/(a(t));

  double integrand, xc;
  double deltaM, Mhalo_smooth;

  double xmin,xmax,dx;
  double qmin,qmax,dq_log; //q is logspaced around temp.
  double mumin,mumax,dmu;

  xmin=0.; //xmax will be z-dependent, since we go up to R_halo.

  mumin=-1.;//cosine always goes between -1 and 1, we do linear spacing.
  mumax=1.;
  dmu=(mumax-mumin)/(Nmu-1.);

  qmin=temp/20.;
  qmax=temp*10.;
  dq_log = log(qmax/qmin)/(Nq-1.); //logspaced is better.

  z_out=zi; //we set it first at z_initial~200 at the begining, for the first test.

  long Nz_solution_max = (log(z_collapse)-logz_solution_min)/dlogz_solution; //to make sure we don't continue after collapse.



	for (i=0;(i<Nz_solution) && (i<Nz_solution_max);i++){ //we stop when we cross z_collapse, which is before zf_code=z_collapse/2. We keep the same array than in collapse code for convenience.
		Mnu_solution[i] = 0.0;
		z_out = exp(logz_solution_min+dlogz_solution*i); //that's where R_solution is saved, and where we save Mnu_solution.
		eta_out=interpol_cubic(zmin_EoS, dz_EoS, etalist, Nz_EoS, z_out);

//    xmax = (1+z_out)*interpol_cubic(logz_solution_min, dlogz_solution, R_solution, Nz_solution, log(z_out)); //comoving radius of halo at z_out.
    xmax = (1+z_out)*R_solution[i];
    dx=(xmax-xmin)/(Nx-1.);

    if (debug_mode>1) 	printf("i=%ld, z=%.1le, eta=%.1le \n",i,z_out,eta_out);


    for(ix=1;ix<Nx;ix++){//comoving position, we avoid x=0.
      x=xmin+dx*ix;
      factorx=1.0/x/x;
      for(iq=0;iq<Nq;iq++){//comoving momentum
        q=qmin*exp(dq_log*iq);
        factorq = exp(q/temp)/pow(1.+exp(q/temp),2.); //assumming FERMIONS.
        for(imu=0;imu<Nmu;imu++){//cosine between them
          mu=mumin + dmu*imu;
          f1[ix][iq][imu] = 0.;//just in case.

          for (iz=0;iz<nzpoints-1;iz++){
                zstep = (zi-z_out)/(nzpoints-1.); //we integrate z_inside (z_in) from z_initial (z_i) to z_outside (z_out).
                z_in = z_out + zstep * iz;
                eta_in=interpol_cubic(zmin_EoS, dz_EoS, etalist, Nz_EoS, z_in);

                OmGbar= pba->Omega0_g * pow(1.+z_in,4.);
                Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z_in,4.);
                OmRbar = OmGbar + Omnu_masslessbar;
                OmM = prf->Omega_M * pow(1.+z_in,3.);

                rho1_z=interpol_cubic(zmin_EoS, dz_EoS, rho1list_EoS, Nz_EoS, z_in);
                rho_ratio_nu1=rho1_z/rho1_0;
                Omnu1 =  prf->Omega_nu1 * rho_ratio_nu1;

                rho2_z=interpol_cubic(zmin_EoS, dz_EoS, rho2list_EoS, Nz_EoS, z_in);
                rho_ratio_nu2=rho2_z/rho2_0;
                Omnu2 =  prf->Omega_nu2 * rho_ratio_nu2;

                rho3_z=interpol_cubic(zmin_EoS, dz_EoS, rho3list_EoS, Nz_EoS, z_in);
                rho_ratio_nu3=rho3_z/rho3_0;
                Omnu3 =  prf->Omega_extra * rho_ratio_nu3;

                H = prf->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1 + Omnu2 + Omnu3);


                xc=(eta_out-eta_in)*q/mass;
                R_halo = interpol_cubic(logz_solution_min, dlogz_solution, R_solution, Nz_solution, log(z_in));
                x_halo = (1+z_in) * R_halo;

                Mhalo_smooth = prf->H0_Mpc*prf->H0_Mpc*pow(R_halo,3.)*OmM/2.;
                deltaM = Mhalo_Mpc - Mhalo_smooth;//Mhalo-Mhalo_smooth. In Mpc.
                integrand = zstep/H * deltaM * (xc/x-mu);

                parenthesis = x*x + xc*xc - 2*x*xc*mu;

                if(parenthesis<x_halo*x_halo){
                  f1[ix][iq][imu]+= integrand*pow(x/x_halo,3.);
                }
                else{
                  f1[ix][iq][imu]+= integrand/pow(parenthesis/x/x,3/2.);
                }

//               printf("z_in=%.1le, deltaM=%.1le, xc=%.1le, x_halo=%.1le, parenthesis=%.1le \n",z_in,deltaM,xc,x_halo,parenthesis);

            }//z-loop
            f1[ix][iq][imu] *= factorq * factorx * factor;

            Mnu_solution[i] += (2*PI*x*x*dx) * (dmu) * (dq_log * q * q * q /(2.*PI*PI) ) * f1[ix][iq][imu]; //\int d^3x d^3q = dx x^2 (2PI) dmu dp p^2/(2PI^2).
        }//mu-loop
      }//q-loop
    }//x-loop

    Mnu_solution[i]*=mass*eVtoMpc/pow(hbareVMpc,3.); //to get Mnu in Mpc we need hbar^3 (that is why there is a 1/(2pi)^3 in p integral). There are 3 q in eV to convert to Mpc-1 (since x is in Mpc),
    //we have to use G to convert mass to Mpc -> eVtoMpc.



//smooth mass of ALL non-CDM components. Just for comparison.
    Mnusmooth = prf->H0_Mpc*prf->H0_Mpc*pow(R_solution[i],3.)*(prf->Omega_nu1+prf->Omega_nu2+prf->Omega_extra)/2.*pow(1+z_out,3.);



	//if(debug_mode>0){
	//	printf("z=%.1le, Mnu/Mhalo=%1le \n",z_out,Mnu_solution[i]/Mhalo_Mpc);
	//	fprintf(fp, "%le  %le  %le \n",z_out, Mnu_solution[i]/Mhalo_Mpc, (Mnu_solution[i]+Mnusmooth)/Mhalo_Mpc);
	//}



  }



  free_array_3D(f1,Nx,Nq);

  return _SUCCESS_;

}


/***BIAS MODULE**/


int get_bias(struct precision * ppr,
 struct perturbs * ppt,
 struct background * pba,
 struct spectra * psp,
 struct relicfast *prf, double z_collapse){





	int i_delta_long, i_klong, i, j;







//which z_collapse are we in:
	int z_ind=find_value(prf->N_zcoll,prf->z_collapse_array, z_collapse);



//we will loop over Mhalo and klong, to avoid reading transfer functions many time.
//there isn't a lot to win from looping over z inside here, as opposed to outside, since you still have to read Boltzmann output every time.


	/////////////////////////////////////////
	//// we read the transfer functions	////
	////////////////////////////////////////



		double *k_transfer_array, **transfer_array_z_coll; //for CDM+b, derivative



		class_alloc(k_transfer_array ,prf->length_transfer*sizeof(double),prf->error_message);
		class_alloc(transfer_array_z_coll,6*sizeof(double*),prf->error_message);
		int al_ind;
		for(al_ind=0;al_ind<6;al_ind++){
      class_alloc(transfer_array_z_coll[al_ind],prf->length_transfer*sizeof(double),prf->error_message);
    }
    class_call(get_transfer(pba, ppt, psp, prf, k_transfer_array, transfer_array_z_coll, z_collapse),
      prf->error_message,
      prf->error_message);


	//and we will compute sigma_M(z_coll).
    double sigmaM_coll;


		/////////////////////////////////////////
		//// these are the HMFs that we use	////
		////////////////////////////////////////

    double HMF;
	//choose your favorite halo mass function:
	//the input (called HMF) is simply dlogn/d delta_crit = d log f/d delta_crit.
	//if your HMF does not have an analytic derivative do numerical instead.



	//Option 1: MICE (0907.0019)
	//the fit is for f_MICE (\sigma), so we transform \sigma -> \sigma/\delta_crit and take the derivative.
	// HMF = -2 c delta_crit/(delta_ref^2 * sigma^2) + a/delta_crit/(1+b*(delta_ref * sigma/delta_crit)^a), where delta_ref=1.686, since they only quote it as a function of \sigma.


    double delta_ref=1.686;
    double a_MICE=1.37*pow(1+z_collapse,-0.15);
    double b_MICE=0.3*pow(1+z_collapse,-0.084);
    double c_MICE=1.036*pow(1+z_collapse,-0.024);

    double HMF_MICE;


	//Option 2: arXiv:1005.2239
	// HMF =d log(n)/d \delta_crit = (q-a (delta_crit/sigma)^2)/delta_crit - (2 p/delta_crit)/(1+(a delta_crit^2/sigma^2)^p)

    double q_hmf=1.795;
    double p_hmf=0.807;
		double a_hmf_0=0.788; // a_hmf (z) = a_hmf_0/(1+z)^a_hmf_exp
		double a_hmf_exp=0.01;
		double a_hmf=a_hmf_0/pow(1.+z_collapse,a_hmf_exp);

		double HMF_Batt;


	//Option 3: Sheth Tormen (astro-ph/9901122)
	// HMF = (1-a (delta_crit/sigma)^2)/delta_crit - (2 p/delta_crit)/(1+(a delta_crit^2/sigma^2)^p)
	// it's essentially the same as Option 2, but with no z-dependence and with q set to unity;

		double a_ST=0.707; //updated from original 0.707, from Ref.1005.2239 Table 3
		double p_ST=0.3;
		double HMF_ST;

		double delta_crit_sim=1.686; //usual value of delta_crit used in simulations to obtain the HMF.


		double k_long; //k_long for each iteration





		double derivative;

		int counter; //to get average of delta_crit.




		const int N_species = 1 + prf->counter_massive_nus + (prf->Omega_extra>0);//how many species (cdm+b) + others
		int species_counter=0;


		const double Omegatot = prf->Omega_M + prf->Omega_nu1 + prf->Omega_nu2 + prf->Omega_extra; //remember we call Omega_M=Omegab+Omegac

		double *fraction_species; //fractions of matter today in each species.
		class_alloc(fraction_species,N_species*sizeof(double),prf->error_message); // cdm, nu1, nu2, extra.

		fraction_species[0] = prf->Omega_M/Omegatot;//CDM+b
		if(prf->m_nu1>0){
			species_counter++;
			fraction_species[species_counter] = prf->Omega_nu1/Omegatot;
		}
		if(prf->m_nu2>0){
			species_counter++;
			fraction_species[species_counter] = prf->Omega_nu2/Omegatot;
		}
		if(prf->Omega_extra>0){
			species_counter++;
			fraction_species[species_counter] = prf->Omega_extra/Omegatot;
		}
		for(species_counter=0;species_counter<N_species;species_counter++){
			if (debug_mode>0) printf("fraction_%d=%.3le \n\n\n",species_counter,fraction_species[species_counter]);
		}
		species_counter=0; //we reset the counter.

		double *tf_species;
		class_alloc(tf_species,N_species*sizeof(double),prf->error_message); // cdm, nu1, nu2, extra.


		double **Power_spectra; //P_ij, for i and j species. At fixed k_long
		class_alloc(Power_spectra,N_species*sizeof(double*),prf->error_message);// cdm, nu1, nu2, extra.
		
		for(al_ind=0;al_ind<N_species;al_ind++){
      class_alloc(Power_spectra[al_ind],N_species*sizeof(double),prf->error_message);
    }


		double factor; //(2pi^2) and other factors.



    int iM;

    for(iM=0;iM<prf->N_Mhalo;iM++){

      double Mhalo = prf->Mhalo_array[iM];
      int m_ind = iM;
		double Mhalo_Mpc = Mhalo * MsuntoMpc;//Mhalo in Mpc (*G)




		sigmaM_coll=getsigma_M(prf, k_transfer_array, transfer_array_z_coll[0], Mhalo_Mpc);


	/////////////////////////////////////
	//// we now compute the bias.	////
	////////////////////////////////////




//we find the "average" value for delta_crit_sim for our collapse. Just to get a fair estimate of the amplitude of b_1^L
		delta_crit_sim=0.;
		counter=0;
		for(i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			for(i_klong=0;i_klong<prf->N_klong;i_klong++){
				counter++;
				delta_crit_sim+=(prf->delta_long_collapse[z_ind][m_ind][i_delta_long][i_klong]+prf->delta_short_crit[z_ind][m_ind][i_delta_long][i_klong]);
			}
		}
		delta_crit_sim/=(1.*counter);
		if(debug_mode>0){
			printf("delta_sim=1.686, delta_us=%.3le \n",delta_crit_sim);
		}


		HMF_MICE= -2. * c_MICE * delta_crit_sim/(delta_ref*delta_ref *sigmaM_coll*sigmaM_coll) +
    a_MICE/delta_crit_sim/(1. + b_MICE*pow(sigmaM_coll*delta_ref/delta_crit_sim,a_MICE));

    HMF_Batt= (q_hmf-a_hmf*pow(delta_crit_sim/sigmaM_coll,2.))/delta_crit_sim - (2*p_hmf/delta_crit_sim)/(1.+pow(a_hmf*pow(delta_crit_sim/sigmaM_coll,2.),p_hmf));

    HMF_ST= (1-a_ST*pow(delta_crit_sim/sigmaM_coll,2.))/delta_crit_sim - (2*p_ST/delta_crit_sim)/(1.+pow(a_ST*pow(delta_crit_sim/sigmaM_coll,2.),p_ST));


	//select a halo mass function, or add your own:
		if(HMF_option==1){ //MICE
			HMF=HMF_MICE;
		}
		else if(HMF_option==2){ //wLCDM
			HMF=HMF_Batt;
		}
		else if(HMF_option==3){ //Sheth-Tormen
			HMF=HMF_ST;
		}
		else {
			printf("Halo Mass Function not chosen, set HMF=HMF_MICE to 1, 2, or 3. MICE selected by default \n");
			HMF=HMF_MICE;
		}



	//	printf("HMF = %le \n",HMF);



		for(i_klong=0;i_klong<prf->N_klong;i_klong++){
			derivative=(prf->delta_short_crit[z_ind][m_ind][N_delta_long-1][i_klong]-prf->delta_short_crit[z_ind][m_ind][0][i_klong])/
					(prf->delta_long_collapse[z_ind][m_ind][N_delta_long-1][i_klong]-prf->delta_long_collapse[z_ind][m_ind][0][i_klong]); // ddelta_crit/ddelta_long
         prf->b_L[z_ind][m_ind][i_klong] = HMF * derivative;
       }


		// for(i_klong=0;i_klong<prf->length_transfer;i_klong++){
		// 	printf("%le \n",k_transfer_array[i_klong]);
		// }


       for(i_klong=0;i_klong<prf->N_klong;i_klong++){
         k_long=prf->klong_list[i_klong];
         factor=(2.*PI*PI)*prf->As*pow(k_long/kpivot,prf->ns-1.)/k_long/k_long/k_long;
         prf->Pmm[z_ind][m_ind][i_klong] = prf->Pmh[z_ind][m_ind][i_klong] = prf->Phh[z_ind][m_ind][i_klong] = 0.;

         tf_species[0] = interpol(transfer_array_z_coll[0],k_transfer_array,prf->length_transfer,k_long);

         if(prf->m_nu1>0){
           species_counter++;
           tf_species[species_counter] = interpol(transfer_array_z_coll[3],k_transfer_array,prf->length_transfer,k_long);
         }
         if(prf->m_nu2>0){
          species_counter++;
          tf_species[species_counter] = interpol(transfer_array_z_coll[4],k_transfer_array,prf->length_transfer,k_long);
        }
        if(prf->Omega_extra>0){
          species_counter++;
          tf_species[species_counter] = interpol(transfer_array_z_coll[5],k_transfer_array,prf->length_transfer,k_long);
        }
			species_counter=0; //we reset the counter.

			for(i=0;i<N_species;i++){
				for(j=0;j<N_species;j++){
					Power_spectra[i][j] = factor * tf_species[i] * tf_species[j]; // P_ij
					prf->Pmm[z_ind][m_ind][i_klong]+= Power_spectra[i][j] * fraction_species[i] * fraction_species[j]; //P_mm
	//				printf("k/h=%.3le, i=%d, j=%d, P=%.1le \n",k_long/h,i,j,Power_spectra[i][j]);
				}
				prf->Pmh[z_ind][m_ind][i_klong]+= (1.0+prf->b_L[z_ind][m_ind][i_klong]) * Power_spectra[i][0] * fraction_species[i]; //P_mhalo
			}
			prf->Phh[z_ind][m_ind][i_klong]+= (1.0+prf->b_L[z_ind][m_ind][i_klong]) * (1.0+prf->b_L[z_ind][m_ind][i_klong]) * Power_spectra[0][0]; //P_hh = (1+bL)^2Pcc


			prf->b_E[z_ind][m_ind][i_klong]=prf->Pmh[z_ind][m_ind][i_klong]/prf->Pmm[z_ind][m_ind][i_klong];

			if(N_species>1 && debug_mode>0){
				printf("k/h=%.1le,  Tx/Tc=%.3le, prf->b_E[z_ind][m_ind]=%.3le, prf->b_L[z_ind][m_ind]=%.3le \n",k_long/pba->h, Power_spectra[1][0]/Power_spectra[0][0], prf->Pmh[z_ind][m_ind][i_klong]/prf->Pmm[z_ind][m_ind][i_klong], prf->b_L[z_ind][m_ind][i_klong]);
			}

		}


}//end of Mhalo loop


//////////////////////////////////////////
//// we free the allocated memory	/////
////////////////////////////////////////






free(k_transfer_array);


free_array_2D(transfer_array_z_coll,6);


free_array_2D(Power_spectra, N_species);

free(fraction_species);
free(tf_species);

return _SUCCESS_;

}


/** PRESSURE MODULE */

double pressure_WDM(double mass, double Temp){
//calculates Pbar(mass,Temp) in Kelvin^4 as in Eq. (19) of Loverde14.
//for a WDM particle of mass mass (in eV) and temperature Temp (in K)

	double moverT=mass/(Temp*KtoeV); //mass in units of Temperature

	int j, Nj=Nkpoints;

	double deltaPmax=30.;
	double deltaPmin=10.;
	double pmin_log=log(1./deltaPmin);
	double pmax_log=log(1.*deltaPmax);//in units of Temp
	double pstep_log=(pmax_log-pmin_log)/(Nj-1.);

	double Pbar=0.; // P_nu = 2 (\int d^3p/(2pi^3) p^2/(3 sqrt(p^2+m^2)) * 1/(exp(p/Tnu)+1)).   //   	printf("%le %le %le \n",Toverm,pmin_log,pmin_log);

	double p, fermi, integrand;// p (in units of Temp), fermi factor f, integrand


	for(j=0;j<Nj-1;j++){
		p=exp(pmin_log+j*pstep_log);//in units of mnu
		fermi=1./(1.+exp(p)); //p in units of Tnu
		integrand= p*p/(3.*sqrt(p*p+moverT*moverT));
		Pbar+= pstep_log * p * p * p / (2.*PI*PI) * fermi * 2 * integrand; //in units of Temp^4, 2 because 2 degrees of freedom
//		printf("p/mnu=%le, fermi=%le, integrand= %le \n",p, fermi, integrand);
	}

	Pbar*=pow(Temp,4.); //in K^4

	return Pbar;
}


double density_WDM(double mass, double Temp){
//calculates rhobar(mass,Temp) in Kelvin^4 as in Eq. (19) of Loverde14.
//for a WDM particle of mass mass (in eV) and temperature Temp (in K)

	double moverT=mass/(Temp*KtoeV); //mnu in units of Tnu.


	int j, Nj=Nkpoints; //number of kpoints. Log-spaced, 30 seems enough.

	double deltaPmax=30.;
	double deltaPmin=10.;
	double pmin_log=log(1./deltaPmin);
	double pmax_log=log(1.*deltaPmax);//in units of Temp
	double pstep_log=(pmax_log-pmin_log)/(Nj-1.); //p will be logspaced in units of Tnu.

	double rhobar=0.; // rho_nu = 2 (\int d^3p/(2pi^3) sqrt(p^2+m^2) * 1/(exp(p/Tnu)+1)).

	double p, fermi, integrand;// p (in units of mnu), fermi factor f, integrand


	for(j=0;j<Nj-1;j++){
		p=exp(pmin_log+j*pstep_log);//in units of mnu
		fermi=1./(1.+exp(p)); //p in units of Tnu
		integrand= sqrt(p*p+moverT*moverT);
		rhobar+= pstep_log * p * p * p / (2.*PI*PI) * fermi * 2 * integrand; //in units of Temp^4, 2 because 2 degrees of freedom	//		printf("p/mnu=%le, fermi=%le, integrand= %le \n",p, fermi, integrand);
	}

	rhobar*=pow(Temp,4.); //in K^4

	return rhobar;
}




double EoS_WDM(double mass, double Temp){
//calculates w, Equation of state, of WDM with mass and Temp, (in eV and K, respectively)

	double w=pressure_WDM(mass,Temp)/density_WDM(mass,Temp);



	return w;
}



/** AUXILLIARY FUNCTIONS */






void free_array_2D(double **array, int n1){
  int i;
  for (i = 0; i < n1; i++){
    free(array[i]);
  }
  free(array);
}


void free_array_3D(double ***array, int n1, int n2){
  int i;  
  for (i = 0; i < n1; i++){
   free_array_2D(array[i],n2);
 }

 free(array);
}

void free_array_4D(double ****array, int n1, int n2, int n3){
  int i;  
  for (i = 0; i < n1; i++){
   free_array_3D(array[i], n2, n3);
 }

 free(array);
}


//////////////////////////////////////
/////	MATHEMATICAL FUNCTIONS:	/////
////////////////////////////////////


double interpol(double data[], double xtab[],int length, double x){
  //Linear interpolation algorithm

  double d1,res,v1,v2,x1,xmin,xmax;
  xmin=xtab[0];
  xmax=xtab[length-1];


  if((x<xmin)) {
    printf("interpol(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
    x=xmin;
  }

  if((x>xmax)) {
    printf("interpol(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
    x=xmax;
  }



  int i=floor(length/2),count;

  for(count=2; i<length-1 && (x>=xtab[i+1] || x<xtab[i]);++count){
   if (x>=xtab[i+1])
    i+=1;
//    		i=i+length/(count);
  else if (x<xtab[i])
   i-=1;
//    		i=i-length/(count);
 else
  printf("ERROR on interpol \n");
}

x1=xtab[i];
v1=data[i];


double step ;

//we avoid the edge by interpolating downwards if it's there.
if(i==length-1){
  i--;
}


step = xtab[i+1]-xtab[i];
v2=data[i+1];

d1=(v2-v1)/step;

res= v1 + d1*(x-x1);


    // printf("%e,%e, %f, %d \n",v1,v2,step,length);
    // printf("x=%e,v1=%e, res=%e \n",x1,v1,res);

return res;

}




double interpol_2D(double **data, double xtab[], int lengthx, double ytab[], int lengthy, double x, double y){
  //Linear interpolation algorithm in 2D
  //assumes both lists are ordered (upwards).

  double d1x,d1y,res,v1,v2x,v2y,x1,y1,xmin,xmax, ymin,ymax;


  xmin=xtab[0];
  xmax=xtab[lengthx-1];
  if((x<xmin)) {
    printf("interpol_2D(x,y) is out of range in x. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
    x=xmin;
  }

  if((x>xmax)) {
    printf("interpol_2D(x,y) is out of range in x. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
    x=xmax;
  }


  ymin=ytab[0];
  ymax=ytab[lengthy-1];
  if((y<ymin)) {
    printf("interpol_2D(x,y) is out of range in y. Range is y=(%f %f) and y=%le \n",ymin,ymax,y);
//        return 0;
    y=ymin;
  }

  if((y>ymax)) {
    printf("interpol_2D(x,y) is out of range in y. Range is y=(%f %f) and y=%le \n",ymin,ymax,y);
//        return 0;
    y=ymax;
  }

  int i, j, count;


  i=floor(lengthx/2);
  for(count=2; i<lengthx-1 && (x>=xtab[i+1] || x<xtab[i]);++count){
   if (x>=xtab[i+1])
    i+=1;
//    		i=i+length/(count);
  else if (x<xtab[i])
   i-=1;
//    		i=i-length/(count);
 else
  printf("ERROR on interpol_2D in x \n");
}


j=floor(lengthy/2);
for(count=2; j<lengthy-1 && (y>=ytab[j+1] || y<ytab[j]);++count){
 if (y>=ytab[j+1])
  j+=1;
//    		i=i+length/(count);
else if (y<ytab[j])
 j-=1;
//    		i=i-length/(count);
else
  printf("ERROR on interpol_2D in y \n");
}

v1=data[i][j];
x1=xtab[i];
y1=ytab[j];


double stepx, stepy;

//we avoid the edge by interpolating downwards if it's there.
if(i==lengthx-1){
  i--;
}

if(j==lengthy-1){
  j--;
}


stepx = xtab[i+1]-xtab[i];
v2x=data[i+1][j];
stepy = ytab[j+1]-ytab[j];
v2y=data[i][j+1];


d1x=(v2x-v1)/stepx;
d1y=(v2y-v1)/stepy;

res= v1 + d1x*(x-x1) +  d1y*(y-y1);

//    printf(" ix=%d (x-x1)=%.3le, iy=%d (y-y1)=%.3le \n", i ,(x-x1) ,j , (y-y1));


//	printf("%e,%e, %f, %d \n",v1,v2,step,length);
//	printf("x=%e, y=%e, v1=%e, res=%e \n",x1,y1,v1,res);

return res;

}






long find_value(long length, double xtab[], double x){
  //finds the position in the list xtab closest to x

  double xmin,xmax;
  xmin=xtab[0];
  xmax=xtab[length-1];


  if((x<xmin)) {
    printf("find_value(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
    return 0;
  }

  if((x>xmax)) {
    printf("find_value(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
    return length -1;
  }

    if(x==xmin){//same but with not warnings.
      return 0;
    }
    if(x==xmax){
      return length-1;
    }



    long i=floor(length/2),count;

    for(count=2; x>=xtab[i+1] || x<xtab[i];++count){
	//printf("finding value %f, currently at %f at index %d \n", x, xtab[i], i);
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
       i-=1;
//    		i=i-length/(count);
     else
      printf("ERROR on find_value \n");
  }

  if (xtab[i+1]-x < x - xtab[i]){
     i++; //if it's closer to the one above add one.
   }



   return i;

 }

 long find_value_reverse(long length, double xtab[], double x){
  //finds the position in the list xtab closest to x
  //here we assume the xlist is reversed in order (goes down)

  double xmin,xmax;
  xmin=xtab[length-1];
  xmax=xtab[0];


  if((x<xmin)) {
    printf("find_value_reverse(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
    x=xmin;
  }

  if((x>xmax)) {
    printf("find_value_reverse(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
    x=xmax;
  }



  long i=floor(length/2),count;

  for(count=2; x<=xtab[i+1] || x>xtab[i];++count){
   if (x<=xtab[i+1])
    i+=1;
//    		i=i+length/(count);
  else if (x>xtab[i])
   i-=1;
//    		i=i-length/(count);
 else
  printf("ERROR on find_value_reverse \n");
}

if (xtab[i+1]-x > x - xtab[i]){
     i++; //if it's closer to the one above add one.
   }



   return i;

 }




 double nintegrate(double data[],double xtab[], int length){

   double sum = 0;
	long double step; //To avoid possible overflows when two values are too close.

	long i;

  for (i = 0; i < length-2; ++i) {
   step = xtab[i+1]-xtab[i];
   sum +=  0.5*(data[i] + data[i+1]) * step;
//        printf("i= %d, sum=%lf \n",i,sum);
 }

 return sum;

}


int triangle(double l1,double l2,double l3){
	//Checking for triangle identities. Returns 0 if not valid and 1 if valid.
	int res=1;

	if((l1+l2-l3<0) || (l1-l2+l3<0) || (-l1+l2+l3<0)){
		res=0;
  }

  return res;
}


double interpol_cubic(double x0, double dx, double *ytab, int Nx, double x) {
  //cubic interpolation routine, assumes x equally spaced.
  long ix;
  double frac;

  if (Nx < 4) {
    fprintf(stderr, "Error: interpol_cubic: Table needs to be of dimension 4 at least\n");
    exit(1);
  }

  // Check if in range
  if (  (dx > 0 && (x<x0 || x>x0+dx*(Nx-1)))
    ||(dx < 0 && (x>x0 || x<x0+dx*(Nx-1))) ) {
   fprintf(stderr, "Error: interpol_cubic: value out of range. Range is (%le,%le) and x=%le. \n",x0,x0+dx*(Nx-1),x);
 //   exit(1);
}

  // Identify location to interpolate
ix = (long)floor((x-x0)/dx);
if (ix<1) ix=1;
if (ix>Nx-3) ix=Nx-3;
frac = (x-x0)/dx-ix;
ytab += ix-1;

  // Return value
return(
  -ytab[0]*frac*(1.-frac)*(2.-frac)/6.
  +ytab[1]*(1.+frac)*(1.-frac)*(2.-frac)/2.
  +ytab[2]*(1.+frac)*frac*(2.-frac)/2.
  -ytab[3]*(1.+frac)*frac*(1.-frac)/6.
  );
}







void reverse(double *list,int N){
//reverse a list. Pretty straightforward.
	double aux;
	int j;
	for(j=0;j<=floor(N/2);j++){
		aux=list[j];
		list[j]=list[N-1-j];
		list[N-1-j]=aux;
	}

}




//for progressbar:
//PBWIDTH and PBSTR defined in .h
void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}











//////////////////////////////////////
/////	COSMOLOGICAL FUNCTIONS:	/////
////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//unused, although they are useful for checks://
double E_LCDM(struct relicfast *prf, double z){
//computes E(z) = H(z)/H0 for given prflogical parameters.

	double E1 = sqrt(prf->Omega_L + prf->Omega_M * pow(1.+z,3.)+ prf->Omega_R * pow(1.+z,4.));

	return E1;

}
double dlogE_dz_LCDM(struct relicfast *prf, double z){
//computes dE(z)/dz * 1/E(z)  for given prflogical parameters.
// for H(z)= Sqrt[OmL + OmM (1 + z)^3 + OmR (1 + z)^4].
//H'(z)/H(z) = (3 OmM (1 + z)^2 + 4 OmR (1 + z)^3)/(2 (OmL + OmM (1 + z)^3 + OmR (1 + z)^4))
//saves the trouble of doing it numerically

	double dE1 = (3.*prf->Omega_M*(1+z)*(1+z) + 4.*prf->Omega_R*pow(1.+z,3.))/(2.*(prf->Omega_L + prf->Omega_M *pow(1.+z,3.) + prf->Omega_R *pow(1.+z,4.)));

	return dE1;

}
////////////////////////////////////////////////////////////////////////





void findeta(struct background *pba, struct relicfast *prf, double *etalist, double zmin_EoS, double dz_EoS, long Nz_EoS,
  double* rholist_nu1_EoS, double* rholist_nu2_EoS, double* rholist_extra_EoS){
  //finds superconformal time etalist, at z_EoS list. Given LCDM + neutrinos and extra species.
  //used for clustering of neutrinos in the BKT approximation.

  double z, H, rho_ratio_nu1, rho_ratio_nu2, rho_ratio_extra;
  double rho_nu1_z, rho_nu1_z0;
  double rho_nu2_z, rho_nu2_z0;
  double rho_extra_z, rho_extra_z0;

  double Omnuextra, OmM, OmL, OmGbar, OmRbar, Omnu_masslessbar;

  long i;

  rho_nu1_z0 = (prf->Omega_nu1>0) * rholist_nu1_EoS[0] + (prf->Omega_nu1==0) * 1; //we make sure that dividing by it does not blow anything up.
  rho_nu2_z0 = (prf->Omega_nu2>0) * rholist_nu2_EoS[0] + (prf->Omega_nu2==0) * 1;
  rho_extra_z0 = (prf->Omega_extra>0) * rholist_extra_EoS[0] + (prf->Omega_extra==0) * 1;

		etalist[0]=0.; //everything will be a function of eta-eta'. so initial value does not matter.

  OmL=prf->Omega_L; //z independent.

  for (i=1;i<Nz_EoS;i++){

   z=zmin_EoS + dz_EoS*i;

   OmGbar= pba->Omega0_g * pow(1.+z,4.);
   Omnu_masslessbar= prf->Omega_nu_massless * pow(1.+z,4.);
   OmRbar = OmGbar + Omnu_masslessbar;
   OmM = prf->Omega_M * pow(1.+z,3.);

   rho_nu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_nu1_EoS, Nz_EoS, z);
   rho_nu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_nu2_EoS, Nz_EoS, z);
   rho_extra_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_extra_EoS, Nz_EoS, z);
   rho_ratio_nu1=rho_nu1_z/rho_nu1_z0;
   rho_ratio_nu2=rho_nu2_z/rho_nu2_z0;
   rho_ratio_extra=rho_extra_z/rho_extra_z0;

   Omnuextra =  prf->Omega_nu1 * rho_ratio_nu1 + prf->Omega_nu2 * rho_ratio_nu2 + prf->Omega_extra * rho_ratio_extra;


			H = prf->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnuextra); //H(zi) in Mpc

			etalist[i]=etalist[i-1] - dz_EoS / H * (1.+z);//d eta = d t/a^2, - sign because eta grows at low z (it's super-conformal time)
//      printf("etalist[i]=%.1le , ratioextra=%.1le \n",etalist[i],rho_ratio_extra);
    } //cubic interpolator for eta.


  }



///////////////////////////////////////////////////////////////////////////////////////
//		Calculate sigma(M) and sigma_8 at some redshift
//
//	input Mhalo_Mpc and transfer-function data
///////////////////////////////////////////////////////////////////////////////////////
  double getsigma_M(struct relicfast *prf, double *k_transfer, double *transfer_function, double Mhalo_Mpc){

	int j, Nj=500; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[prf->length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = pow(prf->Omega_M/2.*prf->H0_Mpc*prf->H0_Mpc/Mhalo_Mpc,-1./3); //Omega_M=Omegac+Omegab. in Mpc.


	double k, window, transfer, Powcc;// k, W(kR), T(k), prf->Pmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol(transfer_function,k_transfer,prf->length_transfer,k);
		Powcc=(2.*PI*PI)*prf->As*pow(k/kpivot,prf->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);	//		printf("k=%le, Powcc=%le, transfer= %le \n",k,Powcc, transfer);
	}

	return sqrt(sigmaM_2);
}



double getsigma_8(struct background *pba, struct relicfast *prf, double *k_transfer, double *transfer_function){
	//for R=8 Mpc/h. It's sigma_8=0.83 for Planck2015 best-fit
	int j, Nj=500; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[prf->length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = 8./pba->h; //In Mpc/h.



	double k, window, transfer, Powcc;// k, W(kR), T(k), prf->Pmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol(transfer_function,k_transfer,prf->length_transfer,k);
		Powcc=(2.*PI*PI)*prf->As*pow(k/kpivot,prf->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);		//		
		//printf("k=%le, Powcc=%le, integrand= %le \n",k,Powcc, k*k * window*window * Powcc / (2.*PI*PI));
	}

	return sqrt(sigmaM_2);
}




double getsigma_M_z(struct relicfast *prf, double z, double *z_transfer, double *k_transfer, double **transfer_function, double Mhalo_Mpc){
//calculates sigma(M) @z.
	int j, Nj=1000; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[prf->length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = pow(prf->Omega_M/2.*prf->H0_Mpc*prf->H0_Mpc/Mhalo_Mpc,-1./3); //Omega_M=Omegac+Omegab. in Mpc.


	double k, window, transfer, Powcc;// k, W(kR), T(k), prf->Pmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol_2D(transfer_function,z_transfer, prf->Nz_transfer, k_transfer,prf->length_transfer,z,k);
		Powcc=(2.*PI*PI)*prf->As*pow(k/kpivot,prf->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);	//		printf("k=%le, Powcc=%le, transfer= %le \n",k,Powcc, transfer);
	}

	return sqrt(sigmaM_2);
}


double getsigma_8_z(struct background *pba, struct relicfast *prf, double z, double *z_transfer, double *k_transfer, double **transfer_function){
//calculates sigma(M) @z.
	int j, Nj=1000; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[prf->length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = 8./pba->h; //In Mpc/h.


	double k, window, transfer, Powcc;// k, W(kR), T(k), prf->PPmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol_2D(transfer_function,z_transfer, prf->Nz_transfer, k_transfer,prf->length_transfer,z ,k );
		Powcc=(2.*PI*PI)*prf->As*pow(k/kpivot,prf->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);	//		printf("k=%le, Powcc=%le, transfer= %le \n",k,Powcc, transfer);
	}

	return sqrt(sigmaM_2);
}
