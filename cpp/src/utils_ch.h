#ifndef __UTILS_CH_H__
#define __UTILS_CH_H__

#include <vector>
#include "chparams.h"

void compute_ch_nonlocal(const std::vector<double> &c,
			 std::vector<double> &dcdt,
			 const double t,
			 CHparamsVector& chpV,
			 SimInfo& info);

void compute_ch_nonlocal_stationary_boundaries(const std::vector<double> &c,
					       std::vector<double> &dcdt,
					       const double t,
					       CHparamsVector& chpV,
					       SimInfo& info);

void compute_ch_nonlocal_neumannBC(const std::vector<double> &c,
				   std::vector<double> &dcdt,
				   const double t,
				   CHparamsVector& chpV,
				   SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_bottom_dirichlet(const std::vector<double> &c,
							       std::vector<double> &dcdt,
							       const double t,
							       CHparamsVector& chpV,
							       SimInfo& info);

void compute_ch_nonlocal_mixedBC_neumann_with_top_dirichlet(const std::vector<double> &c,
							    std::vector<double> &dcdt,
							    const double t,
							    CHparamsVector& chpV,
							    SimInfo& info);

std::vector<double>& set_boundary_values_to_zero( std::vector<double> &dcdt ,
						  SimInfo& info );

std::vector<double>& apply_dirichlet_bc( std::vector<double>& c ,
					 SimInfo& info );

std::vector<double>& apply_neumann_bc( std::vector<double>& c ,
				       SimInfo& info );

std::vector<double>& apply_mixed_bc_neumann_with_bottom_dirichlet( std::vector<double>& c ,
								   SimInfo& info );

std::vector<double>& apply_mixed_bc_neumann_with_top_dirichlet( std::vector<double>& c ,
								SimInfo& info );

std::vector<double>& freeze_corners( std::vector<double>& dcdt ,
				     SimInfo& info );

double laplace_component(int i ,
                         const std::vector<double>& c ,
                         const std::vector<double>& u ,
                         const std::vector<double>& b );

CHparamsVector compute_chparams_using_temperature( CHparamsVector& chpV0 ,
						   SimInfo& info,
						   std::vector<double>& T );

double convert_temperature_to_flory_huggins( CHparamsVector& chpV ,
					     SimInfo& info,
					     const double T );

CHparamsVector compute_eps2_and_sigma_from_polymer_params( CHparamsVector& chpV0 ,
							   SimInfo& info,
							   std::vector<double>& T );


/*********************
Functions for implicit
**********************/

void compute_ch_nonlocal_implicit_dirichlet(const std::vector<double> &c,
                                            std::vector<double> &dcdt,
                                            const double t,
                                            CHparamsVector& chpV,
                                            SimInfo& info ,
                                            TS& ts ,
                                            Vec& F ,
                                            DM& da );

class ThirteenPointStencil {

  public:

    ThirteenPointStencil() { };
   ~ThirteenPointStencil() { };
    double c_i,c_im2,c_im1,c_ip1,c_ip2,c_jm2,c_jm1,c_jp1,c_jp2, c_ul,c_ur,c_bl,c_br;

};

ThirteenPointStencil get_stencil_dirichlet_bc( const double dirichlet_bc ,
                                               const double** uarray ,
                                               const int Mx ,
                                               const int My ,
                                               const int i ,
                                               const int j );


#endif
