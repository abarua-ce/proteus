#ifndef TADR_H
#define TADR_H
#include <cmath>
#include <iostream>
#include <valarray>
#include <random>
#include "CompKernel.h"
#include "ModelFactory.h"
#include "ArgumentsDict.h"
#include "xtensor-python/pyarray.hpp"
#define nnz nSpace

namespace py = pybind11;

#define POWER_SMOOTHNESS_INDICATOR 2
#define IS_BETAij_ONE 0

// Cell based methods:
//    * Galerkin (unstabilized)
//    * VMS(SUPG) with BDF1 or BDF2 time integration
//    * Explicit Taylor Galerkin with EV stabilization
// Edge based methods.
//    Low order via D. Kuzmin's
//    High order methods: Smoothness indicator with MC, EV commutator with MC, D.K with ML
//    Zalesak's FCT

namespace proteus
{
  enum class STABILIZATION : int { Galerkin=-1, VMS=0, TaylorGalerkinEV=1, EntropyViscosity=2, SmoothnessIndicator=3, Kuzmin=4};
  enum class ENTROPY : int { POWER=0, LOG=1};
  // Power entropy //
  inline double EPOWER(const double& phi, const double& phiL, const double& phiR)
  {
    return 1./2.*std::pow(fabs(phi),2.);
  }
  inline double DEPOWER(const double& phi, const double& phiL, const double& phiR)
  {
    return fabs(phi)*(phi>=0 ? 1 : -1);
  }
  // Log entropy // for level set from 0 to 1
  inline double ELOG(const double& phi, const double& phiL, const double& phiR)
  {
    return std::log(fabs((phi-phiL)*(phiR-phi))+1E-14);
  }
  inline double DELOG(const double& phi, const double& phiL, const double& phiR)
  {
    return (phiL+phiR-2*phi)*((phi-phiL)*(phiR-phi)>=0 ? 1 : -1)/(fabs((phi-phiL)*(phiR-phi))+1E-14);
  }
}

namespace proteus
{
  class TADR_base
  {
    //The base class defining the interface
  public:
    std::valarray<double> Rpos, Rneg;
    std::valarray<double> FluxCorrectionMatrix;
    std::valarray<double> TransportMatrix, DiffusionMatrix, TransposeTransportMatrix;
    std::valarray<double> psi, eta, global_entropy_residual, boundary_integral;
    std::valarray<double> maxVel,maxEntRes;
    virtual ~TADR_base(){}
    virtual void calculateResidual(arguments_dict& args)=0;
    virtual void calculateJacobian(arguments_dict& args)=0;
    virtual void FCTStep(arguments_dict& args)=0;
    virtual void Update_concentration_RWPT(arguments_dict& args)=0;
  };

  template<class CompKernelType,
           int nSpace,
           int nQuadraturePoints_element,
           int nDOF_mesh_trial_element,
           int nDOF_trial_element,
           int nDOF_test_element,
           int nQuadraturePoints_elementBoundary>
  class TADR : public TADR_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    TADR():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}

    inline
    void calculateCFL(const double& elementDiameter,
                      const double df[nSpace],
                      double& cfl)
    {
      double h,nrm_v;
      h = elementDiameter;
      nrm_v=0.0;
      for(int I=0;I<nSpace;I++)
        nrm_v+=df[I]*df[I];
      nrm_v = sqrt(nrm_v);
      cfl = nrm_v/h;
    }

    inline
// void evaluateCoefficients(const int rowptr[nSpace],
// 				                      const int colind[nnz],
//                               const double v[nSpace],
//                               const double* q_a,
//                               const double& u,
//                               double& m,
//                               double& dm,
//                               double f[nSpace],
//                               double df[nSpace], 
//                               double a[nnz],
//                               double da[nnz])
//     {
//       m = u;
//       dm= 1.0;
//       for (int I=0; I < nSpace; I++)
//         {
//           f[I] = v[I]*u;
//           df[I] = v[I];
//            for (int ii = rowptr[I]; ii < rowptr[I + 1]; ii++)
//         {
//             a[ii] = q_a[ii];
//             da[ii] = 0.0;
//         }

//         } 
//     }


void evaluateCoefficients(const int rowptr[nSpace], 
      const int colind[nnz],
      const double v[nSpace],
      const double alpha_L,
      const double alpha_T,
      const double Dm, // Molecular diffusion term
      //const double* q_a,
      const double& u,
      double& m,
      double& dm,
      double f[nSpace],
      double df[nSpace], 
      double a[nnz],
      double da[nnz])
{
    m = u;
    dm = 1.0;

// Compute velocity magnitude
    double v_mag = 0.0;
    for (int I = 0; I < nSpace; I++) {
    v_mag += v[I] * v[I];
    }
    v_mag = std::sqrt(v_mag);

// Compute unit velocity vector
    double v_unit[nSpace] = {0.0};
    if (v_mag > 1e-10) { // Avoid division by zero
      for (int I = 0; I < nSpace; I++) {
          v_unit[I] = v[I] / v_mag;
            }
      }

    for (int I = 0; I < nSpace; I++) {
    f[I] = v[I] * u;
    df[I] = v[I];

    for (int ii = rowptr[I]; ii < rowptr[I + 1]; ii++) {
    int J = colind[ii]; // Column index
    double deltaIJ = (I == J) ? 1.0 : 0.0;

    // Compute dispersion tensor elements, incorporating molecular diffusion (Dm)
    double dispersion_ij = 
    Dm * deltaIJ +  // Molecular diffusion (isotropic contribution)
    alpha_L * v_unit[I] * v_unit[J] * v_mag +  // Longitudinal dispersion
    alpha_T * v_mag * (deltaIJ - v_unit[I] * v_unit[J]);  // Transverse dispersion

    a[ii] = dispersion_ij;
    da[ii] = 0.0; // Assuming no dependency of dispersion on u
      }
      }
    }



inline
    void exteriorNumericalDiffusiveFlux(int* rowptr,
					int* colind,
					const int& isDOFBoundary,
					const int& isDiffusiveFluxBoundary,
					const double n[nSpace],
					double* bc_a,
					const double& bc_u,
					const double& bc_flux,
					double* a,
					const double grad_potential[nSpace],
					const double& u,
					const double& penalty,
					double& flux)
    {
      double diffusiveVelocityComponent_I;
      double penaltyFlux;
      double max_a;
      if (isDiffusiveFluxBoundary == 1)
	{
	  flux = bc_flux;
	}
      else if (isDOFBoundary == 1)
	{
	  flux = 0.0;
	  max_a = 0.0;
	  for (int I = 0; I < nSpace; I++)
	    {
	      diffusiveVelocityComponent_I = 0.0;
	      for (int m = rowptr[I]; m < rowptr[I+1]; m++)
		{
		  diffusiveVelocityComponent_I -= a[m] * grad_potential[colind[m]];
		  max_a = fmax(max_a, a[m]);
		}
	      flux += diffusiveVelocityComponent_I * n[I];
	    }
	  penaltyFlux = max_a * penalty * (u - bc_u);
	  flux += penaltyFlux;
	}
      else
	{
	  //std::cerr << "warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0" << std::endl;
	  flux = 0.0;
	}
    }

    inline
    double ExteriorNumericalDiffusiveFluxJacobian(int* rowptr,
						  int* colind,
						  const int& isDOFBoundary,
						  const int& isDiffusiveFluxBoundary,
						  const double n[nSpace],
						  double* a,
						  const double& v,
						  const double grad_v[nSpace],
						  const double& penalty)
    {
      double dvel_I;
      double tmp = 0.0;
      double max_a = 0.0;
      if ((isDiffusiveFluxBoundary == 0) && (isDOFBoundary == 1))
	{
	  for (int I = 0; I < nSpace; I++)
	    {
	      dvel_I = 0.0;
	      for (int m = rowptr[I]; m < rowptr[I + 1]; m++)
		{
		  dvel_I -= a[m] * grad_v[colind[m]];
		  max_a = fmax(max_a, a[m]);
		}
	      tmp += dvel_I * n[I];
	    }
	  tmp += max_a * penalty * v;
	}
      return tmp;
    }

    inline
    void calculateSubgridError_tau(const double& elementDiameter,
                                   const double& dmt,
                                   const double df[nSpace],
                                   double& cfl,
                                   double& tau)
    {
      //regular elements
      double h,nrm_v,oneByAbsdt;
      h = elementDiameter;
      nrm_v=0.0;
      for(int I=0;I<nSpace;I++)
        nrm_v+=df[I]*df[I];
      nrm_v = sqrt(nrm_v);
      cfl = nrm_v/h;
      oneByAbsdt =  fabs(dmt);
      tau = 1.0/(2.0*nrm_v/h + oneByAbsdt + 1.0e-8);
    }

    inline
    void calculateSubgridError_tau(const double&  Ct_sge,
                                   const double   G[nSpace*nSpace],
                                   const double&  A0,
                                   const double   Ai[nSpace],
                                   double& tau_v,
                                   double& cfl)
    {
      //metric-based tau for arbitrarily shaped elements
      double v_d_Gv=0.0;
      for(int I=0;I<nSpace;I++)
        {for (int J=0;J<nSpace;J++)
            v_d_Gv += Ai[I]*G[I*nSpace+J]*Ai[J];
          v_d_Gv += Ai[I]*G[I*nSpace+I];
          for(int J=0;J<nSpace;J++)
            {
              if(J!=I)
                v_d_Gv += 2.0*Ai[I]*G[I*nSpace+J];
            }
        }
      tau_v = 1.0/sqrt(Ct_sge*A0*A0 + v_d_Gv + 1.0e-8);
    }

    inline
    void calculateNumericalDiffusion(const double& shockCapturingDiffusion,
                                     const double& elementDiameter,
                                     const double& strong_residual,
                                     const double grad_u[nSpace],
                                     double& numDiff)
    {
      double h,
        num,
        den,
        n_grad_u;
      h = elementDiameter;
      n_grad_u = 0.0;
      for (int I=0;I<nSpace;I++)
        n_grad_u += grad_u[I]*grad_u[I];
      num = shockCapturingDiffusion*0.5*h*fabs(strong_residual);
      den = sqrt(n_grad_u) + 1.0e-8;
      numDiff = num/den;
    }

    inline
    void exteriorNumericalAdvectiveFlux(const int& isDOFBoundary_u,
                                        const int& isFluxBoundary_u,
                                        const double n[nSpace],
                                        const double& bc_u,
                                        const double& bc_flux_u,
                                        const double& u,
                                        const double velocity[nSpace],
                                        double& flux)
    {
    //   // Debug: Print input values
    // std:: cout << "Exterior Numerical Advective Flux Function"<< std::endl;
    // std::cout << "isDOFBoundary_u: " << isDOFBoundary_u << ", isFluxBoundary_u: " << isFluxBoundary_u << std::endl;
    // std::cout << "bc_u: " << bc_u << ", bc_flux_u: " << bc_flux_u << ", u: " << u << std::endl;
    // std::cout << "velocity: [" << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << "]" << std::endl;
    // std::cout << "Normal: [" << n[0] << ", " << n[1] << ", " << n[2] << "]\n";
    // std ::cout <<"End Line"<< "\n";

      double flow=0.0;
      for (int I=0; I < nSpace; I++)
        flow += n[I]*velocity[I];

      if (isDOFBoundary_u == 1)
        {
          if (flow >= 0.0)
            {
              flux = u*flow;
            }
          else
            {
              flux = bc_u*flow;
            }
        }
      else if (isFluxBoundary_u == 1)
        {
          flux = bc_flux_u;
        }
      else
        {
          if (flow >= 0.0)
            {
              flux = u*flow;
            }
          else
            {
              //std::cout<<"warning: TADR open boundary with no external trace, setting to zero for inflow"<<std::endl;
              flux = 0.0;
            }

        }
    }

    inline
    void exteriorNumericalAdvectiveFluxDerivative(const int& isDOFBoundary_u,
                                                  const int& isFluxBoundary_u,
                                                  const double n[nSpace],
                                                  const double velocity[nSpace],
                                                  double& dflux)
    {
      double flow=0.0;
      for (int I=0; I < nSpace; I++)
        flow += n[I]*velocity[I];
        
      dflux=0.0;//default to no flux
      if (isDOFBoundary_u == 1)
        {
          if (flow >= 0.0)
            {
              dflux = flow;
            }
          else
            {
              dflux = 0.0;
            }
        }
      else if (isFluxBoundary_u == 1)
        {
          dflux = 0.0;
        }
      else
        {
          if (flow >= 0.0)
            {
              dflux = flow;
            }
        }
    }
    inline

  void exteriorNumericalDiffusiveFluxDerivative(const int& isDOFBoundary,
                                                const int& isDiffusiveFluxBoundary,
                                                const int rowptr[nSpace],
                                                const int colind[nnz],
                                                const double n[nSpace],
                                                const double a[nnz],
                                                const double da[nnz],
                                                const double grad_psi[nSpace],
                                                const double grad_v[nSpace],
                                                const double v[nSpace],
                                                const double penalty,
                                                double& fluxJacobian)
{
    if (isDiffusiveFluxBoundary == 0 && isDOFBoundary == 1)
    {
        fluxJacobian = 0.0;
        for (int I = 0; I < nSpace; I++) {
            fluxJacobian += penalty * v[I];
            for(int m=rowptr[I]; m<rowptr[I+1]; m++)
        {
            fluxJacobian -= (a[m] * grad_v[colind[m]] + da[m] * v[I] * grad_psi[colind[m]]) * n[I];
        }
        }
    }
    else
    {
        fluxJacobian = 0.0;
    }
}

///////////////////////////ADD RWPT////////////////////////////////////

    inline double normal_random()
    {
      static thread_local std::mt19937 generator;
      static thread_local std::normal_distribution<double> distribution(0.0,1.0);
      return distribution(generator);
    }
    // Minimal lower-triangular Cholesky; assumes SPD and succeeds.
    inline void chol_lower(const double* A, double* L) {
      for (int i=0;i<nSpace;i++)
        for (int j=0;j<nSpace;j++)
          L[i*nSpace+j] = 0.0;

      for (int i=0;i<nSpace;i++) {
        for (int j=0;j<=i;j++) {
          double s = A[i*nSpace+j];
          for (int k=0;k<j;k++) s -= L[i*nSpace+k]*L[j*nSpace+k];
          if (i==j) L[i*nSpace+j] = std::sqrt(s);
          else      L[i*nSpace+j] = s / L[j*nSpace+j];
        }
      }
    }

  // --- small fixed-size inverse for nSpace=1/2/3 ---
inline bool invert_small_matrix(const double* A, double* Ainv)
{
  if constexpr (nSpace == 1) {
    const double det = A[0];
    if (std::fabs(det) < 1e-16) return false;
    Ainv[0] = 1.0 / det;
    return true;
  } else if constexpr (nSpace == 2) {
    const double a=A[0], b=A[1], c=A[2], d=A[3];
    const double det = a*d - b*c;
    if (std::fabs(det) < 1e-16) return false;
    const double invDet = 1.0/det;
    Ainv[0] =  d*invDet; Ainv[1] = -b*invDet;
    Ainv[2] = -c*invDet; Ainv[3] =  a*invDet;
    return true;
  } else { // nSpace==3
    const double a=A[0], b=A[1], c=A[2],
                 d=A[3], e=A[4], f=A[5],
                 g=A[6], h=A[7], i=A[8];
    const double A11 =  (e*i - f*h);
    const double A12 = -(d*i - f*g);
    const double A13 =  (d*h - e*g);
    const double det = a*A11 + b*A12 + c*A13;
    if (std::fabs(det) < 1e-16) return false;
    const double invDet = 1.0/det;
    Ainv[0] = A11*invDet;
    Ainv[1] = (c*h - b*i)*invDet;
    Ainv[2] = (b*f - c*e)*invDet;
    Ainv[3] = A12*invDet;
    Ainv[4] = (a*i - c*g)*invDet;
    Ainv[5] = (c*d - a*f)*invDet;
    Ainv[6] = A13*invDet;
    Ainv[7] = (b*g - a*h)*invDet;
    Ainv[8] = (a*e - b*d)*invDet;
    return true;
  }
}

// --- conservative P1 deposit: particles -> solution DOFs ---
inline void deposit_particles_to_nodes_conservative(
    const double*              p_X,                     // [pid * nSpace]
    const int                  pid,
    const double               m_p,                     // mass per particle
    const xt::pyarray<double>& mesh_dof,                // geometry coords [nNodes,3]
    const xt::pyarray<int>&    mesh_l2g,                // geom connectivity [nElem, nSpace+1]
    const xt::pyarray<int>&    u_l2g,                   // solution connectivity [nElem, nSpace+1]
    const int                  nElements_global,
    const int                  mesh_nDOF_trial_element, // == nSpace+1
    const int                  u_nDOF_trial_element,    // == nSpace+1
    xt::pyarray<double>&       b,                       // RHS mass per DOF (accumulate)
    long long&                 n_outside)
{
  if (mesh_nDOF_trial_element != nSpace+1 || u_nDOF_trial_element != nSpace+1) {
    std::cerr << "[RWPT][deposit] ERROR: need P1 simplex; got mesh="
              << mesh_nDOF_trial_element << " u=" << u_nDOF_trial_element
              << " while nSpace+1=" << (nSpace+1) << "\n";
    return;
  }

  const double tol = -1e-12; // allow tiny negative lambda due to FP error
  double total_mass_RHS = 0.0;

  for (int p=0; p<pid; ++p) {
    const double* xp = &p_X[p*nSpace];
    bool placed = false;

    for (int e=0; e<nElements_global && !placed; ++e) {
      const int* l2g_geom = &mesh_l2g.data()[e*mesh_nDOF_trial_element];

      // vertices X0..Xd
      double X0[nSpace];
      for (int i=0;i<nSpace;++i)
        X0[i] = mesh_dof.data()[l2g_geom[0]*3 + i];

      // A = [X1-X0, ..., Xd-X0]
      double A[nSpace*nSpace];
      for (int col=0; col<nSpace; ++col) {
        const int a = col+1;
        for (int i=0;i<nSpace;++i) {
          const double Xi = mesh_dof.data()[l2g_geom[a]*3 + i];
          A[i*nSpace + col] = Xi - X0[i];
        }
      }

      double Ainv[nSpace*nSpace];
      if (!invert_small_matrix(A, Ainv))
        continue; // degenerate element

      // alpha = A^{-1}(xp - X0)
      double rhs[nSpace];
      for (int i=0;i<nSpace;++i) rhs[i] = xp[i] - X0[i];
      double alpha[nSpace];
      for (int i=0;i<nSpace;++i) {
        double s=0; for (int j=0;j<nSpace;++j) s += Ainv[i*nSpace+j]*rhs[j];
        alpha[i] = s;
      }

      // barycentric lambdas
      double lambda[nSpace+1];
      double sumA=0; for (int i=0;i<nSpace;++i) sumA += alpha[i];
      lambda[0] = 1.0 - sumA;
      for (int i=0;i<nSpace;++i) lambda[i+1] = alpha[i];

      // inside test
      bool inside=true;
      for (int i=0;i<nSpace+1;++i) if (lambda[i] < tol) { inside=false; break; }
      if (!inside) continue;

      // deposit to solution DOFs
      const int* l2g_u = &u_l2g.data()[e*u_nDOF_trial_element];
      for (int j=0;j<u_nDOF_trial_element;++j) {
        const int I = l2g_u[j];
        const double add = m_p * lambda[j];
        b.data()[I] += add;
        total_mass_RHS += add;
      }
      placed = true;
    } // elements

    if (!placed) ++n_outside;
  } // particles

  std::cout << "[RWPT] deposit: pid=" << pid
            << " outside=" << n_outside
            << " total_mass_RHS=" << std::setprecision(16) << total_mass_RHS
            << " (should be pid * m_p = " << (static_cast<double>(pid)*m_p) << ")\n";
}


    //  inline void displacement(const double v[nSpace],
    //                           const double grad_theta[nSpace],
    //                           const double divD[nSpace],
    //                           const double theta,
    //                           const double alpha_L,
    //                           const double alpha_T,
    //                           const double Dm,
    //                           const double dt,
    //                           const double Z[nSpace],
    //                           double dx[nSpace])
    //   {
    //     // 1) Velocity magnitude and unit vector e
    //     double vmag = 0.0;
    //     for (int i=0;i<nSpace;i++) vmag += v[i]*v[i];
    //     vmag = std::sqrt(vmag);

    //     double e[nSpace];
    //     if (vmag > 0.0) {
    //       for (int i=0;i<nSpace;i++) e[i] = v[i]/vmag;
    //     } else {
    //       for (int i=0;i<nSpace;i++) e[i] = 0.0;
    //     }

    //     // 2) Dispersion tensor D = Dm I + vmag*( alpha_T I + (alpha_L-alpha_T) e e^T )
    //     double D[nSpace*nSpace];
    //     for (int I=0; I<nSpace; ++I) {
    //       for (int J=0; J<nSpace; ++J) {
    //         const double deltaIJ = (I==J) ? 1.0 : 0.0;
    //         D[I*nSpace+J] =
    //             Dm * deltaIJ
    //           + vmag * ( alpha_T*deltaIJ + (alpha_L - alpha_T)*e[I]*e[J] );
    //       }
    //     }

    //     // 3) Drift A = v + divD + (1/theta) D grad_theta
    //     double D_grad_theta[nSpace];
    //     for (int i=0;i<nSpace;i++) {
    //       double s = 0.0;
    //       for (int j=0;j<nSpace;j++) s += D[i*nSpace+j]*grad_theta[j];
    //       D_grad_theta[i] = s;
    //     }
    //     double A[nSpace];
    //     for (int i=0;i<nSpace;i++) {
    //       A[i] = v[i] + divD[i] + (D_grad_theta[i]/theta);
    //     }

    //     // 4) Cholesky of M = 2D -> L (lower)
    //     double M[nSpace*nSpace];
    //     for (int i=0;i<nSpace;i++)
    //       for (int j=0;j<nSpace;j++)
    //         M[i*nSpace+j] = 2.0 * D[i*nSpace+j];

    //     double L[nSpace*nSpace];
    //     chol_lower(M, L);

    //     // 5) Stochastic increment eta = L * Z * sqrt(dt)
    //     const double sdt = std::sqrt(dt);
    //     double eta[nSpace];
    //     for (int i=0;i<nSpace;i++) {
    //       double s = 0.0;
    //       for (int j=0;j<=i;j++) s += L[i*nSpace+j]*Z[j];
    //       eta[i] = s * sdt;
    //     }

    //     // 6) Final displacement
    //     for (int i=0;i<nSpace;i++) dx[i] = A[i]*dt + eta[i];
    //   }

inline void displacement(const double v[nSpace],
                         const double grad_theta[nSpace],
                         const double divD[nSpace],
                         const double theta,
                         const double alpha_L,
                         const double alpha_T,
                         const double Dm,
                         const double dt,
                         const double Z[nSpace],   // Z ~ N(0,I)
                         double dx[nSpace])
{
  // ---- guards -------------------------------------------------------------
  const double eps     = 1e-14;
  const double sdt     = (dt > 0.0) ? std::sqrt(dt) : 0.0;
  const double theta_s = (theta > eps) ? theta : 1.0;  // safe θ in drift

  // ---- velocity magnitude and unit direction e ----------------------------
  double vmag2 = 0.0;
  for (int i=0;i<nSpace;++i) vmag2 += v[i]*v[i];
  const double vmag = std::sqrt(vmag2);

  double e[nSpace] = {0.0};
  if (vmag > eps) {
    for (int i=0;i<nSpace;++i) e[i] = v[i]/vmag;
  } else {
    // no flow: pick arbitrary unit vector (x-axis)
    e[0] = 1.0;
    for (int i=1;i<nSpace;++i) e[i] = 0.0;
  }

  // ---- dispersion scalars (principal values) ------------------------------
  const double DL = std::max(0.0, Dm + alpha_L*vmag);  // along e
  const double DT = std::max(0.0, Dm + alpha_T*vmag);  // transverse

  // ---- drift term A = v + divD + (1/θ) D ∇θ --------------------------------
  double e_dot_grad = 0.0;
  for (int i=0;i<nSpace;++i) e_dot_grad += e[i]*grad_theta[i];

  double D_grad_theta[nSpace] = {0.0};
  for (int i=0;i<nSpace;++i) {
    D_grad_theta[i] = Dm*grad_theta[i]
                    + vmag*( alpha_T*grad_theta[i]
                           + (alpha_L - alpha_T)*e_dot_grad*e[i] );
  }

  double A[nSpace] = {0.0};
  for (int i=0;i<nSpace;++i)
    A[i] = v[i] + divD[i] + D_grad_theta[i]/theta_s;

  double t1[nSpace] = {0.0}, t2[nSpace] = {0.0};
  if (nSpace >= 2) {
    // pick a vector not colinear with e
    if (std::fabs(e[0]) < 0.9) { t1[0]=0; t1[1]= (nSpace>2? e[2]:0); t1[2]= (nSpace>2? -e[1]:0); }
    else                        { t1[0]= (nSpace>2? -e[2]:0); t1[1]=0; t1[2]= e[0]; }

    // normalize t1
    double n1 = 0.0; for (int i=0;i<nSpace;++i) n1 += t1[i]*t1[i];
    n1 = std::sqrt(std::max(n1,0.0));
    if (n1 > eps) for (int i=0;i<nSpace;++i) t1[i] /= n1; else { t1[0]=0; t1[1]=1; if(nSpace>2) t1[2]=0; }

    // t2 = e x t1 (only meaningful for 3D; in 2D we keep t2=0)
    if (nSpace == 3) {
      t2[0] = e[1]*t1[2] - e[2]*t1[1];
      t2[1] = e[2]*t1[0] - e[0]*t1[2];
      t2[2] = e[0]*t1[1] - e[1]*t1[0];
    }
  }
  // ---- stochastic displ. dW ~ N(0,2D dt) via principal axes ---------------
  const double sigmaL = std::sqrt(2.0*DL)*sdt;
  const double sigmaT = std::sqrt(2.0*DT)*sdt;

  // Compose random vector: Z[0] along e, Z[1] (and Z[2]) transverse
  double dRand[nSpace] = {0.0};
  for (int i=0;i<nSpace;++i) {
    dRand[i]  = sigmaL * Z[0] * e[i];
    if (nSpace >= 2) dRand[i] += sigmaT * Z[1] * t1[i];
    if (nSpace == 3) dRand[i] += sigmaT * Z[2] * t2[i];
  }

  // ---- final displacement --------------------------------------------------
  for (int i=0;i<nSpace;++i) dx[i] = A[i]*dt + dRand[i];

  // ---- last-resort NaN/Inf guard ------------------------------------------
  for (int i=0;i<nSpace;++i)
    if (!std::isfinite(dx[i])) dx[i] = 0.0;
}




 void calculateResidual(arguments_dict& args)
    {
      double dt = args.scalar<double>("dt");
      xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
      xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
      xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
      xt::pyarray<double>& mesh_velocity_dof = args.array<double>("mesh_velocity_dof");
      double MOVING_DOMAIN = args.scalar<double>("MOVING_DOMAIN");
      xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
      xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
      xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
      xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
      xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
      xt::pyarray<double>& u_grad_test_ref = args.array<double>("u_grad_test_ref");
      xt::pyarray<double>& mesh_trial_trace_ref = args.array<double>("mesh_trial_trace_ref");
      xt::pyarray<double>& mesh_grad_trial_trace_ref = args.array<double>("mesh_grad_trial_trace_ref");
      xt::pyarray<double>& dS_ref = args.array<double>("dS_ref");
      xt::pyarray<double>& u_trial_trace_ref = args.array<double>("u_trial_trace_ref");
      xt::pyarray<double>& u_grad_trial_trace_ref = args.array<double>("u_grad_trial_trace_ref");
      xt::pyarray<double>& u_test_trace_ref = args.array<double>("u_test_trace_ref");
      xt::pyarray<double>& u_grad_test_trace_ref = args.array<double>("u_grad_test_trace_ref");
      xt::pyarray<double>& normal_ref = args.array<double>("normal_ref");
      xt::pyarray<double>& boundaryJac_ref = args.array<double>("boundaryJac_ref");
      int nElements_global = args.scalar<int>("nElements_global");
      double useMetrics = args.scalar<double>("useMetrics");
      double alphaBDF = args.scalar<double>("alphaBDF");
      int lag_shockCapturing = args.scalar<int>("lag_shockCapturing");
      double shockCapturingDiffusion = args.scalar<double>("shockCapturingDiffusion");
      double sc_uref = args.scalar<double>("sc_uref");
      double sc_alpha = args.scalar<double>("sc_alpha");
      xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
      xt::pyarray<int>& r_l2g = args.array<int>("r_l2g");
      xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
      double degree_polynomial = args.scalar<double>("degree_polynomial");
      xt::pyarray<double>& u_dof = args.array<double>("u_dof");
      xt::pyarray<double>& u_dof_old = args.array<double>("u_dof_old");
      xt::pyarray<double>& velocity = args.array<double>("velocity");
      xt::pyarray<double>& q_m = args.array<double>("q_m");
      xt::pyarray<double>& q_u = args.array<double>("q_u");

      xt::pyarray<double>& q_a = args.array<double>("q_a");
      xt::pyarray<double>& q_r = args.array<double>("q_r");
      xt::pyarray<double>& ebqe_a = args.array<double>("ebq_a");
      

      xt::pyarray<double>& q_m_betaBDF = args.array<double>("q_m_betaBDF");
      xt::pyarray<double>& q_dV = args.array<double>("q_dV");
      xt::pyarray<double>& q_dV_last = args.array<double>("q_dV_last");
      xt::pyarray<double>& cfl = args.array<double>("cfl");
      xt::pyarray<double>& edge_based_cfl = args.array<double>("edge_based_cfl");
      xt::pyarray<double>& q_numDiff_u = args.array<double>("q_numDiff_u");
      xt::pyarray<double>& q_numDiff_u_last = args.array<double>("q_numDiff_u_last");
      int offset_u = args.scalar<int>("offset_u");
      int stride_u = args.scalar<int>("stride_u");
      xt::pyarray<double>& globalResidual = args.array<double>("globalResidual");
      int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
      xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
      xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
      xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
      xt::pyarray<double>& ebqe_velocity_ext = args.array<double>("ebqe_velocity_ext");
      xt::pyarray<int>& isDOFBoundary_u = args.array<int>("isDOFBoundary_u");
      xt::pyarray<double>& ebqe_bc_u_ext = args.array<double>("ebqe_bc_u_ext");
      xt::pyarray<int>& isFluxBoundary_u = args.array<int>("isFluxBoundary_u");
      xt::pyarray<double>& ebqe_bc_flux_u_ext = args.array<double>("ebqe_bc_flux_u_ext");
      xt::pyarray<double>& ebqe_bc_diffusiveFlux_u_ext = args.array<double>("ebqe_bc_diffusiveFlux_u_ext");
      
      double epsFact = args.scalar<double>("epsFact");
      xt::pyarray<double>& ebqe_u = args.array<double>("ebqe_u");
      xt::pyarray<double>& ebqe_flux = args.array<double>("ebqe_flux");
      int stage = args.scalar<int>("stage");
      xt::pyarray<double>&  uTilde_dof = args.array<double>("uTilde_dof");
      double cE = args.scalar<double>("cE");
      double cMax = args.scalar<double>("cMax");
      double cK = args.scalar<double>("cK");
      double uL = args.scalar<double>("uL");
      double uR = args.scalar<double>("uR");
      int numDOFs = args.scalar<int>("numDOFs");
      int NNZ = args.scalar<int>("NNZ");
      xt::pyarray<int>& csrRowIndeces_DofLoops = args.array<int>("csrRowIndeces_DofLoops");
      xt::pyarray<int>& csrColumnOffsets_DofLoops = args.array<int>("csrColumnOffsets_DofLoops");
      xt::pyarray<int>& csrRowIndeces_CellLoops = args.array<int>("csrRowIndeces_CellLoops");
      xt::pyarray<int>& csrColumnOffsets_CellLoops = args.array<int>("csrColumnOffsets_CellLoops");
      xt::pyarray<int>& csrColumnOffsets_eb_CellLoops = args.array<int>("csrColumnOffsets_eb_CellLoops");
      xt::pyarray<double>& ML = args.array<double>("ML");
      int LUMPED_MASS_MATRIX = args.scalar<int>("LUMPED_MASS_MATRIX");
      STABILIZATION STABILIZATION_TYPE = static_cast<STABILIZATION>(args.scalar<int>("STABILIZATION_TYPE"));
      ENTROPY ENTROPY_TYPE = static_cast<ENTROPY>(args.scalar<int>("ENTROPY_TYPE"));    
      //STABILIZATION STABILIZATION_TYPE{args.scalar<int>("STABILIZATION_TYPE")};
      //ENTROPY ENTROPY_TYPE{args.scalar<int>("ENTROPY_TYPE")};
      xt::pyarray<double>& uLow = args.array<double>("uLow");
      xt::pyarray<double>& dLow = args.array<double>("dLow");
      xt::pyarray<double>& dt_times_dH_minus_dL = args.array<double>("dt_times_dH_minus_dL");
      xt::pyarray<double>& min_u_bc = args.array<double>("min_u_bc");
      xt::pyarray<double>& max_u_bc = args.array<double>("max_u_bc");
      xt::pyarray<double>& quantDOFs = args.array<double>("quantDOFs");
      /////////////////////////////////////////////////////////////////////////
      xt::pyarray<int>& a_rowptr = args.array<int>("a_rowptr");
      xt::pyarray<int>& a_colind = args.array<int>("a_colind");
      //xt::pyarray<double>& D = args.array<double>("D");
      //initializeDToZero(D);
      ///////////////////////////////////////////
      xt::pyarray<int>& isDiffusiveFluxBoundary_u = args.array<int>("isDiffusiveFluxBoundary_u");
      xt::pyarray<int>& isAdvectiveFluxBoundary_u = args.array<int>("isAdvectiveFluxBoundary_u");
      xt::pyarray<double>& ebqe_bc_advectiveFlux_u_ext = args.array<double>("ebqe_bc_advectiveFlux_u_ext");
      xt::pyarray<double>& ebqe_penalty_ext = args.array<double>("ebqe_penalty_ext");
      //////////////////////////////////////////////////////////////////////////
      double physicalDiffusion = args.scalar<double>("physicalDiffusion");
      const double eb_adjoint_sigma = args.scalar<double>("eb_adjoint_sigma");
      // Extract alpha_L, alpha_T, and Dm from args
      const double alphaL_val = args.scalar<double>("alpha_L"); // Longitudinal dispersion coefficient
      const double alphaT_val = args.scalar<double>("alpha_T"); // Transverse dispersion coefficient
      const double Dm_val = args.scalar<double>("Dm");         // Molecular diffusion coefficient

      double meanEntropy = 0., meanOmega = 0., maxEntropy = -1E10, minEntropy = 1E10;
      maxVel.resize(nElements_global, 0.0);
      maxEntRes.resize(nElements_global, 0.0);
      double Ct_sge = 4.0;
      if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity or 
          STABILIZATION_TYPE==STABILIZATION::SmoothnessIndicator or 
          STABILIZATION_TYPE==STABILIZATION::Kuzmin)
        {
          TransportMatrix.resize(NNZ,0.0);
          DiffusionMatrix.resize(NNZ,0.0);
          TransposeTransportMatrix.resize(NNZ,0.0);
          // compute entropy and init global_entropy_residual and boundary_integral
          psi.resize(numDOFs,0.0);
          eta.resize(numDOFs,0.0);
          global_entropy_residual.resize(numDOFs,0.0);
          boundary_integral.resize(numDOFs,0.0);
          for (int i=0; i<numDOFs; i++)
            {
              // NODAL ENTROPY //
              if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity) //EV stab
                {
                  eta[i] = ENTROPY_TYPE == ENTROPY::POWER ? EPOWER(u_dof_old.data()[i],uL,uR) : ELOG(u_dof_old.data()[i],uL,uR);
                  global_entropy_residual[i]=0.;
                }
              boundary_integral[i]=0.;
            }
        }
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      //eN is the element index
      //eN_k is the quadrature point index for a scalar
      //eN_k_nSpace is the quadrature point index for a vector
      //eN_i is the element test function index
      //eN_j is the element trial function index
      //eN_k_j is the quadrature point index for a trial function
      //eN_k_i is the quadrature point index for a trial function
      for(int eN=0;eN<nElements_global;eN++)
        {
          //declare local storage for element residual and initialize
          double
            elementResidual_u[nDOF_test_element],
            element_entropy_residual[nDOF_test_element];
          double  elementTransport[nDOF_test_element][nDOF_trial_element];
          double  elementDiffusion[nDOF_test_element][nDOF_trial_element];
          double  elementTransposeTransport[nDOF_test_element][nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            {
              elementResidual_u[i]=0.0;
            }//i
          if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity or 
              STABILIZATION_TYPE==STABILIZATION::SmoothnessIndicator or 
              STABILIZATION_TYPE==STABILIZATION::Kuzmin)
            {
              for (int i=0;i<nDOF_test_element;i++)
                {
                  element_entropy_residual[i]=0.0;
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      elementTransport[i][j]=0.0;
                      elementDiffusion[i][j]=0.0;
                      elementTransposeTransport[i][j]=0.0;
                    }
                }
            }
          //loop over quadrature points and compute integrands
          for  (int k=0;k<nQuadraturePoints_element;k++)
            {
              //compute indeces and declare local storage
              int eN_k = eN*nQuadraturePoints_element+k,
                eN_k_nSpace = eN_k*nSpace,
                eN_nDOF_trial_element = eN*nDOF_trial_element;
                //int index_D = eN_k * a_rowptr.data()[nSpace];
              double
                entVisc_minus_artComp,
                u=0.0,un=0.0,
                grad_u[nSpace],grad_u_old[nSpace],grad_uTilde[nSpace],
                m=0.0,dm=0.0,mn=0.0,dmn=0.0,
                H=0.0,Hn=0.0,HTilde=0.0,
                f[nSpace],fn[nSpace],df[nSpace],dfn[nSpace],
                ////////////////////////////////////////////
                //a[nSpace], da[nSpace], an[nSpace], dan[nSpace],
                a[nnz], da[nnz], an[nnz], dan[nnz],
                m_t=0.0,dm_t=0.0,
                pdeResidual_u=0.0,
                Lstar_u[nDOF_test_element],
                subgridError_u=0.0,
                tau=0.0,tau0=0.0,tau1=0.0,
                numDiff0=0.0,numDiff1=0.0,
                jac[nSpace*nSpace],
                jacDet,
                jacInv[nSpace*nSpace],
                u_grad_trial[nDOF_trial_element*nSpace],
                u_test_dV[nDOF_trial_element],
                u_grad_test_dV[nDOF_test_element*nSpace],
                dV,x,y,z,xt,yt,zt,
                G[nSpace*nSpace],G_dd_G,tr_G,
                // for entropy residual
                aux_entropy_residual=0.0, DENTROPY_un, DENTROPY_uni;//norm_Rv;

              ck.calculateMapping_element(eN,
                                          k,
                                          mesh_dof.data(),
                                          mesh_l2g.data(),
                                          mesh_trial_ref.data(),
                                          mesh_grad_trial_ref.data(),
                                          jac,
                                          jacDet,
                                          jacInv,
                                          x,y,z);
              ck.calculateMappingVelocity_element(eN,
                                                  k,
                                                  mesh_velocity_dof.data(),
                                                  mesh_l2g.data(),
                                                  mesh_trial_ref.data(),
                                                  xt,yt,zt);
              //get the physical integration weight
              dV = fabs(jacDet)*dV_ref.data()[k];
              ck.calculateG(jacInv,G,G_dd_G,tr_G);
              //get the trial function gradients
              ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],
                                  jacInv,
                                  u_grad_trial);
              //get the solution
              ck.valFromDOF(u_dof.data(),
                            &u_l2g.data()[eN_nDOF_trial_element],
                            &u_trial_ref.data()[k*nDOF_trial_element],
                            u);
              ck.valFromDOF(u_dof_old.data(),
                            &u_l2g.data()[eN_nDOF_trial_element],
                            &u_trial_ref.data()[k*nDOF_trial_element],
                            un);
              //get the solution gradients
              ck.gradFromDOF(u_dof.data(),
                             &u_l2g.data()[eN_nDOF_trial_element],
                             u_grad_trial,
                             grad_u);
              ck.gradFromDOF(u_dof_old.data(),
                             &u_l2g.data()[eN_nDOF_trial_element],
                             u_grad_trial,
                             grad_u_old);
              ck.gradFromDOF(uTilde_dof.data(),
                             &u_l2g.data()[eN_nDOF_trial_element],
                             u_grad_trial,
                             grad_uTilde);
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;
                  for (int I=0;I<nSpace;I++)
                    {
                      u_grad_test_dV[j*nSpace+I] = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
                    }
                }

              //
              //
              //calculate pde coefficients at quadrature points

              //const double* q_a_ptr = &q_a[eN_k * a_rowptr[nSpace]];
              evaluateCoefficients(a_rowptr.data(),
				                           a_colind.data(),
                                   &velocity.data()[eN_k_nSpace],
                                   alphaL_val, alphaT_val, Dm_val, 
                                   //q_a_ptr, //&q_a.data()[eN_k*a_rowptr.data()[nSpace]],//[eN_k*nnz],
                                   u,
                                   m,
                                   dm,
                                   f,
                                   df,
                                   a,
                                   da);

            //D.data()[eN * nQuadraturePoints_element * nnz + k * nnz] calculates 
            //the correct starting index for D for the current element eN and 
            //quadrature point k, considering nnz as the number of non-zero entries for the sparse matrix. This ensures that the values for D are accessed correctly based on the element and quadrature point.

              //cek todo, this should be velocity_old
              evaluateCoefficients(a_rowptr.data(),
				                           a_colind.data(),
                                   &velocity.data()[eN_k_nSpace],
                                   alphaL_val, alphaT_val, Dm_val, 
                                   //q_a_ptr, //&q_a.data()[eN_k*a_rowptr.data()[nSpace]],//[eN_k*a_rowptr.data()[nSpace]],//[eN_k*nnz],
                                   un,
                                   mn,
                                   dmn,
                                   fn,
                                   dfn, 
                                   an, 
                                   dan);
              //an= &q_a.data()[eN_k * sd_rowptr.data()[nSpace]];

              //
              //moving mesh
              //
              double mesh_velocity[3];
              mesh_velocity[0] = xt;
              mesh_velocity[1] = yt;
              mesh_velocity[2] = zt;

              for (int I=0;I<nSpace;I++)
                {
                  f[I] -= MOVING_DOMAIN*m*mesh_velocity[I];
                  df[I] -= MOVING_DOMAIN*dm*mesh_velocity[I];
                  fn[I] -= MOVING_DOMAIN*mn*mesh_velocity[I];
                  dfn[I] -= MOVING_DOMAIN*dmn*mesh_velocity[I];
                }
              //
              //calculate time derivative at quadrature points
              //
              if (q_dV_last.data()[eN_k] <= -100)
                q_dV_last.data()[eN_k] = dV;
              q_dV.data()[eN_k] = dV;
              ck.bdf(alphaBDF,
                     q_m_betaBDF.data()[eN_k]*q_dV_last.data()[eN_k]/dV,//ensure prior mass integral is correct for  m_t with BDF1
                     m,
                     dm,
                     m_t,
                     dm_t);

              if (STABILIZATION_TYPE==STABILIZATION::TaylorGalerkinEV)
                {
                  double normVel=0., norm_grad_un=0.;
                  for (int I=0;I<nSpace;I++)
                    {
                      Hn += dfn[I]*grad_u_old[I];
                      HTilde += dfn[I]*grad_uTilde[I];
                      fn[I] = dfn[I]*un-MOVING_DOMAIN*m*mesh_velocity[I];//cek check this for moving domain
                      H += dfn[I]*grad_u[I];
                      normVel += dfn[I]*df[I];
                      norm_grad_un += grad_u_old[I]*grad_u_old[I];
                    }
                  normVel = std::sqrt(normVel);
                  norm_grad_un = std::sqrt(norm_grad_un)+1E-10;

                  // calculate CFL
                  calculateCFL(elementDiameter.data()[eN]/degree_polynomial,dfn,cfl.data()[eN_k]);


                  // compute max velocity at cell
                  maxVel[eN] = fmax(normVel,maxVel[eN]);

                  // Strong entropy residual
                  double entRes = (EPOWER(u,0,1)-EPOWER(un,0,1))/dt + 0.5*(DEPOWER(u,0,1)*H + DEPOWER(un,0,1)*Hn);
                  maxEntRes[eN] = fmax(maxEntRes[eN],fabs(entRes));

                  // Quantities for normalization factor //
                  meanEntropy += EPOWER(u,0,1)*dV;
                  meanOmega += dV;
                  maxEntropy = fmax(maxEntropy,EPOWER(u,0,1));
                  minEntropy = fmin(minEntropy,EPOWER(u,0,1));

                  // artificial compression
                  double hK=elementDiameter.data()[eN]/degree_polynomial;
                  entVisc_minus_artComp = fmax(1-cK*fmax(un*(1-un),0)/hK/norm_grad_un,0);
                }
              else if (STABILIZATION_TYPE==STABILIZATION::VMS)
                {
                  //
                  //calculate subgrid error (strong residual and adjoint)
                  //
                  //calculate strong residual
                  pdeResidual_u = ck.Mass_strong(m_t) + ck.Advection_strong(df,grad_u);
                  //calculate adjoint
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      int i_nSpace = i*nSpace;
                      Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
                    }
                  //calculate tau and tau*Res
                  calculateSubgridError_tau(elementDiameter.data()[eN],dm_t,df,cfl.data()[eN_k],tau0);
                  calculateSubgridError_tau(Ct_sge,
                                            G,
                                            dm_t,
                                            df,
                                            tau1,
                                            cfl.data()[eN_k]);
                  tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

                  subgridError_u = -tau*pdeResidual_u;
                  //
                  //calculate shock capturing diffusion
                  //
                  ck.calculateNumericalDiffusion(shockCapturingDiffusion,
                                                 elementDiameter.data()[eN],
                                                 pdeResidual_u,
                                                 grad_u,
                                                 numDiff0);
                  ck.calculateNumericalDiffusion(shockCapturingDiffusion,
                                                 sc_uref,
                                                 sc_alpha,
                                                 G,
                                                 G_dd_G,
                                                 pdeResidual_u,
                                                 grad_u,
                                                 numDiff1);
                  q_numDiff_u.data()[eN_k] = useMetrics*numDiff1+(1.0-useMetrics)*numDiff0;
                }
              else if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity)
              {
                for (int I=0;I<nSpace;I++)
                  aux_entropy_residual += dfn[I]*grad_u_old[I];
                DENTROPY_un = ENTROPY_TYPE==ENTROPY::POWER ? DEPOWER(un,uL,uR) : DELOG(un,uL,uR);
                calculateCFL(elementDiameter.data()[eN]/degree_polynomial,dfn,cfl.data()[eN_k]);
              }
              else
                calculateCFL(elementDiameter.data()[eN]/degree_polynomial,dfn,cfl.data()[eN_k]);

              for(int i=0;i<nDOF_test_element;i++)
                {
                  //int eN_k_i=eN_k*nDOF_test_element+i,
                  //eN_k_i_nSpace = eN_k_i*nSpace,
                  int i_nSpace=i*nSpace;
                  if (STABILIZATION_TYPE==STABILIZATION::TaylorGalerkinEV)
                    {
                      if (stage == 1)
                        elementResidual_u[i] +=
                          ck.Mass_weak(dt*m_t,u_test_dV[i]) +  // time derivative
                          1./3*dt*(ck.Advection_weak(fn,&u_grad_test_dV[i_nSpace]) +
                                   ck.Diffusion_weak(a_rowptr.data(),a_colind.data(),a,grad_u,&u_grad_test_dV[i_nSpace])+ 
                                   ck.NumericalDiffusion(physicalDiffusion, grad_u_old, &u_grad_test_dV[i_nSpace])) +
                          1./9*dt*dt*ck.NumericalDiffusion(Hn,dfn,&u_grad_test_dV[i_nSpace]) +
                          1./3*dt*entVisc_minus_artComp*ck.NumericalDiffusion(q_numDiff_u_last.data()[eN_k]+physicalDiffusion,
                                                                              grad_u_old,
                                                                              &u_grad_test_dV[i_nSpace]);
                      // TODO: Add part about moving mesh
                      else //stage 2
                        elementResidual_u[i] +=
                          ck.Mass_weak(dt*m_t,u_test_dV[i]) +  // time derivative
                          dt*(ck.Advection_weak(fn,&u_grad_test_dV[i_nSpace]) + 
                              ck.Diffusion_weak(a_rowptr.data(),a_colind.data(),an,grad_u,&u_grad_test_dV[i_nSpace])+
                              ck.NumericalDiffusion(physicalDiffusion, grad_u_old, &u_grad_test_dV[i_nSpace])) +
                          0.5*dt*dt*ck.NumericalDiffusion(HTilde,dfn,&u_grad_test_dV[i_nSpace]) +
                          dt*entVisc_minus_artComp*ck.NumericalDiffusion(q_numDiff_u_last.data()[eN_k]+physicalDiffusion,
                                                                         grad_u_old,
                                                                         &u_grad_test_dV[i_nSpace]);
                    }
                  else if (STABILIZATION_TYPE==STABILIZATION::VMS)
                    {
                      elementResidual_u[i] +=
                        ck.Mass_weak(m_t,u_test_dV[i]) +
                        ck.Advection_weak(f,&u_grad_test_dV[i_nSpace]) +
                        ck.Diffusion_weak(a_rowptr.data(),a_colind.data(),a,grad_u,&u_grad_test_dV[i_nSpace]) +    
                        ck.SubgridError(subgridError_u,Lstar_u[i]) +
                        ck.NumericalDiffusion(q_numDiff_u_last.data()[eN_k] + physicalDiffusion,//todo add full sparse diffusion terms
                                              grad_u,
                                              &u_grad_test_dV[i_nSpace]);
                    }
                  else if(STABILIZATION_TYPE==STABILIZATION::EntropyViscosity or 
                          STABILIZATION_TYPE==STABILIZATION::SmoothnessIndicator or 
                          STABILIZATION_TYPE==STABILIZATION::Kuzmin)
                    {
                      // VECTOR OF ENTROPY RESIDUAL //
                      int eN_i=eN*nDOF_test_element+i;
                      if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity) // EV stab
                        {
                          int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
                          DENTROPY_uni = ENTROPY_TYPE == ENTROPY::POWER ? DEPOWER(u_dof_old.data()[gi],uL,uR) : DELOG(u_dof_old.data()[gi],uL,uR);
                          element_entropy_residual[i] += (DENTROPY_un - DENTROPY_uni)*aux_entropy_residual*u_test_dV[i];
                        }
                      elementResidual_u[i] += (u-un)*u_test_dV[i];
                      ///////////////
                      // j-th LOOP // To construct transport matrices
                      ///////////////
                      for(int j=0;j<nDOF_trial_element;j++)
                        {
                          int j_nSpace = j*nSpace;
                          int i_nSpace = i*nSpace;
                          //cek todo, see if we really need elementTransposeTransport
                          // (can we just swap indices on local transport matrix?)
                          //or event TransposeTransportMatrix (can we just use pointers?)
                          elementTransport[i][j] +=
                            ck.AdvectionJacobian_weak(dfn,
                                                      u_trial_ref.data()[k*nDOF_trial_element+j],
                                                      &u_grad_test_dV[i_nSpace])
                                                      +
                            ck.SimpleDiffusionJacobian_weak(a_rowptr.data(),
										                                        a_colind.data(),
                                                            a,
                                                            &u_grad_trial[j_nSpace],
                                                            &u_grad_test_dV[i_nSpace]);
                                                      


                                               
                           elementDiffusion[i][j] += ck.NumericalDiffusionJacobian(physicalDiffusion,
                                                                                   &u_grad_trial[j_nSpace],
                                                                                   &u_grad_test_dV[i_nSpace]);
                          elementTransposeTransport[i][j] += ck.AdvectionJacobian_weak(dfn,
                                                                                       u_trial_ref.data()[k*nDOF_trial_element+i],
                                                                                       &u_grad_test_dV[j_nSpace])+
                                                             ck.SimpleDiffusionJacobian_weak(a_rowptr.data(),
                                                                                              a_colind.data(),
                                                                                              a,
                                                                                              &u_grad_trial[j_nSpace],
                                                                                              &u_grad_test_dV[i_nSpace]);                          
                                                                                       
                                                            
                                                          }
                    }
                  else
                    {
                      elementResidual_u[i] +=
                        ck.Mass_weak(m_t,u_test_dV[i]) +
                        ck.Advection_weak(f,&u_grad_test_dV[i_nSpace])+
                        ck.Diffusion_weak(a_rowptr.data(),a_colind.data(),a,grad_u,&u_grad_test_dV[i_nSpace]);

                        //std::cout << "elementResidual_u[" << i << "]: " << elementResidual_u[i] << std::endl;
                        //  +
                        

                        // ck.NumericalDiffusion(physicalDiffusion,//todo add full sparse diffusion terms
                        //                       grad_u,
                        //                       &u_grad_test_dV[i_nSpace]);
                    }
                }//i
              //
              //save solution for other models
              //
              q_u.data()[eN_k] = u;
              q_m.data()[eN_k] = m;
              }//k
          //
          //load element into global residual and save element residual
          //
          for(int i=0;i<nDOF_test_element;i++)
            {
              int eN_i=eN*nDOF_test_element+i;
              int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
              globalResidual.data()[gi] += elementResidual_u[i];
              if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity or 
                  STABILIZATION_TYPE==STABILIZATION::SmoothnessIndicator or 
                  STABILIZATION_TYPE==STABILIZATION::Kuzmin)
                {
                  
                  // distribute entropy_residual
                  if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity) // EV Stab
                    global_entropy_residual[gi] += element_entropy_residual[i];
                  // distribute transport matrices
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      int eN_i_j = eN_i*nDOF_trial_element+j;
                      TransportMatrix[csrRowIndeces_CellLoops.data()[eN_i] +
                                      csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransport[i][j];
                      DiffusionMatrix[csrRowIndeces_CellLoops.data()[eN_i] +
                                      csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementDiffusion[i][j];
                      TransposeTransportMatrix[csrRowIndeces_CellLoops.data()[eN_i] +
                                               csrColumnOffsets_CellLoops.data()[eN_i_j]]
                        += elementTransposeTransport[i][j];
                    }//j
                }//edge-based
            }//i
        }//eN
      //
      //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
        {
          double min_u_bc_local = 1E10, max_u_bc_local = -1E10;
          int ebN = exteriorElementBoundariesArray.data()[ebNE],
            eN  = elementBoundaryElementsArray.data()[ebN*2+0],
            ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0],
            eN_nDOF_trial_element = eN*nDOF_trial_element;
          double elementResidual_u[nDOF_test_element],
            fluxTransport[nDOF_test_element][nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            {
              elementResidual_u[i]=0.0;
              for (int j=0;j<nDOF_trial_element;j++)
                fluxTransport[i][j] = 0.0;
            }
          for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
            {
              int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
                ebNE_kb_nSpace = ebNE_kb*nSpace,
                ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                ebN_local_kb_nSpace = ebN_local_kb*nSpace;
              double u_ext=0.0,
                grad_u_ext[nSpace],
                m_ext=0.0,
                dm_ext=0.0,
                f_ext[nSpace],
                df_ext[nSpace],
                /////////////////////

                a_ext[nnz],
		            da_ext[nnz],

                bc_a_ext[nnz],
		            bc_da_ext[nnz],

                /////////////////////////////
                flux_ext=0.0,
                dflux_u_u_ext=0.0,
                bc_u_ext=0.0,
                bc_m_ext=0.0,
                bc_dm_ext=0.0,

                flux_diff_ext=0.0,
                difffluxjacobian_ext=0.0,
                bc_f_ext[nSpace],
                bc_df_ext[nSpace],
                jac_ext[nSpace*nSpace],
                jacDet_ext,
                jacInv_ext[nSpace*nSpace],
                boundaryJac[nSpace*(nSpace-1)],
                metricTensor[(nSpace-1)*(nSpace-1)],
                metricTensorDetSqrt,
                dS,
                u_test_dS[nDOF_test_element],
                u_grad_trial_trace[nDOF_trial_element*nSpace],
                u_grad_test_dS[nDOF_trial_element*nSpace],
                normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                
                G[nSpace*nSpace],G_dd_G,tr_G;
                
              //
              //calculate the solution and gradients at quadrature points
              //
              //compute information about mapping from reference element to physical element
              ck.calculateMapping_elementBoundary(eN,
                                                  ebN_local,
                                                  kb,
                                                  ebN_local_kb,
                                                  mesh_dof.data(),
                                                  mesh_l2g.data(),
                                                  mesh_trial_trace_ref.data(),
                                                  mesh_grad_trial_trace_ref.data(),
                                                  boundaryJac_ref.data(),
                                                  jac_ext,
                                                  jacDet_ext,
                                                  jacInv_ext,
                                                  boundaryJac,
                                                  metricTensor,
                                                  metricTensorDetSqrt,
                                                  normal_ref.data(),
                                                  normal,
                                                  x_ext,y_ext,z_ext);
              ck.calculateMappingVelocity_elementBoundary(eN,
                                                          ebN_local,
                                                          kb,
                                                          ebN_local_kb,
                                                          mesh_velocity_dof.data(),
                                                          mesh_l2g.data(),
                                                          mesh_trial_trace_ref.data(),
                                                          xt_ext,yt_ext,zt_ext,
                                                          normal,
                                                          boundaryJac,
                                                          metricTensor,
                                                          integralScaling);
              dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];
              //get the metric tensor
              //cek todo use symmetry
              ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
              //compute shape and solution information
              //shape
              ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],
                                  jacInv_ext,
                                  u_grad_trial_trace);
              //solution and gradients
              if (STABILIZATION_TYPE==STABILIZATION::TaylorGalerkinEV) //explicit
                {
                  ck.valFromDOF(u_dof_old.data(),
                                &u_l2g.data()[eN_nDOF_trial_element],
                                &u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],
                                u_ext);
                  ck.gradFromDOF(u_dof_old.data(),
                                 &u_l2g.data()[eN_nDOF_trial_element],
                                 u_grad_trial_trace,
                                 grad_u_ext);
                }
              else
                {
                  ck.valFromDOF(u_dof.data(),
                                &u_l2g.data()[eN_nDOF_trial_element],
                                &u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],
                                u_ext);
                  ck.gradFromDOF(u_dof.data(),
                                 &u_l2g.data()[eN_nDOF_trial_element],
                                 u_grad_trial_trace,
                                 grad_u_ext);
                }
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                  for (int I=0;I<nSpace;I++)
		                u_grad_test_dS[j*nSpace+I] = u_grad_trial_trace[j*nSpace+I]*dS;//cek hack, using trial
	
                }
              //
              //load the boundary values
              //
              bc_u_ext = isDOFBoundary_u.data()[ebNE_kb]*ebqe_bc_u_ext.data()[ebNE_kb]+
                          (1-isDOFBoundary_u.data()[ebNE_kb])*u_ext;

                      //
              //
              //calculate the pde coefficients using the solution and the boundary values for the solution
              //
              const double* qb_a_ptr = &ebqe_a[ebNE_kb * a_rowptr[nSpace]];
              evaluateCoefficients(a_rowptr.data(),
				                           a_colind.data(),
                                   &ebqe_velocity_ext.data()[ebNE_kb_nSpace],
                                   alphaL_val, alphaT_val, Dm_val, 
                                   //qb_a_ptr, //&q_a.data()[ebNE * nQuadraturePoints_element * nnz + kb * nnz],//[ebNE_kb*a_rowptr.data()[nSpace]],//[ebNE_kb*nnz],
                                   u_ext,
                                   m_ext,
                                   dm_ext,
                                   f_ext,
                                   df_ext,
                                   a_ext,
                                   da_ext);
              
              evaluateCoefficients(a_rowptr.data(),
				                           a_colind.data(),
                                   &ebqe_velocity_ext.data()[ebNE_kb_nSpace],
                                   alphaL_val, alphaT_val, Dm_val, 
                                   //qb_a_ptr,//&q_a.data()[ebNE * nQuadraturePoints_element * nnz + kb * nnz],//[ebNE_kb*a_rowptr.data()[nSpace]],//[ebNE_kb*nnz],
                                   bc_u_ext,
                                   bc_m_ext,
                                   bc_dm_ext,
                                   bc_f_ext,
                                   bc_df_ext,
                                   bc_a_ext,
                                   bc_da_ext);         
              //moving mesh
              //
              double mesh_velocity[3];
              mesh_velocity[0] = xt_ext;
              mesh_velocity[1] = yt_ext;
              mesh_velocity[2] = zt_ext;

              for (int I=0;I<nSpace;I++)
                {
                  f_ext[I] -= MOVING_DOMAIN*m_ext*mesh_velocity[I];
                  df_ext[I] -= MOVING_DOMAIN*dm_ext*mesh_velocity[I];
                  bc_f_ext[I] -= MOVING_DOMAIN*bc_m_ext*mesh_velocity[I];
                  bc_df_ext[I] -= MOVING_DOMAIN*bc_dm_ext*mesh_velocity[I];
                }
              //
              //calculate the numerical fluxes
              //
              exteriorNumericalDiffusiveFlux(a_rowptr.data(),
                                             a_colind.data(),
                                             isDOFBoundary_u.data()[ebNE_kb],
                                             isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                             normal,
                                             a_ext,
                                             bc_u_ext,
                                             ebqe_bc_diffusiveFlux_u_ext.data()[ebNE_kb],
                                             a_ext,
                                             grad_u_ext,
                                             u_ext,
                                             ebqe_penalty_ext.data()[ebNE_kb],
                                             flux_diff_ext);

              exteriorNumericalAdvectiveFlux(isDOFBoundary_u.data()[ebNE_kb],
                                             isFluxBoundary_u.data()[ebNE_kb],
                                             normal,
                                             bc_u_ext,
                                             ebqe_bc_flux_u_ext.data()[ebNE_kb],
                                             u_ext,
                                             df_ext,
                                             flux_ext);
            
                  
    
              
               //std::cout<<"Advection EXT"<<flux_ext<<std::endl;
               //std::cout<<"Diffusion  Ext"<<flux_diff_ext<<std::endl;
                    
              flux_ext += flux_diff_ext;  
              //std::cout<<"Combine EXT"<<flux_ext<<std::endl;
              ebqe_flux.data()[ebNE_kb] = flux_ext;
              //save for other models? cek need to be consistent with numerical flux
              if(flux_ext >= 0.0)
                ebqe_u.data()[ebNE_kb] = u_ext;
              else
                ebqe_u.data()[ebNE_kb] = bc_u_ext;              
              if (STABILIZATION_TYPE==STABILIZATION::TaylorGalerkinEV)
                {
                  if (stage == 1)
                    flux_ext *= 1./3*dt;
                  else
                    flux_ext *= dt;
                }

              //
              //update residuals
              //
              //cek todo, these are brought in from EV residual and are not correct
              //the element residual should be updated and the global residual and transport updated
              //after the closure of the quadrature loop
              for (int i=0;i<nDOF_test_element;i++)
                {
                  if (STABILIZATION_TYPE == STABILIZATION::Galerkin or 
                      STABILIZATION_TYPE == STABILIZATION::VMS or 
                      STABILIZATION_TYPE == STABILIZATION::TaylorGalerkinEV) 
                      {
                    elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i])+
                                            ck.ExteriorElementBoundaryDiffusionAdjoint(isDOFBoundary_u.data()[ebNE_kb],
                                                                                            isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                                                            eb_adjoint_sigma,
                                                                                            u_ext,
                                                                                            bc_u_ext,
                                                                                            normal,
                                                                                            a_rowptr.data(),
                                                                                            a_colind.data(),
                                                                                            a_ext,
                                                                                            &u_grad_test_dS[i*nSpace]);
                      }
                  else if (STABILIZATION_TYPE == STABILIZATION::EntropyViscosity or 
                      STABILIZATION_TYPE == STABILIZATION::SmoothnessIndicator or 
                      STABILIZATION_TYPE == STABILIZATION::Kuzmin)
                    {                 
                      exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u.data()[ebNE_kb],
                                                           isFluxBoundary_u.data()[ebNE_kb],
                                                           normal,
                                                           df_ext,
                                                           dflux_u_u_ext);  
                      exteriorNumericalDiffusiveFluxDerivative(isDOFBoundary_u.data()[ebNE_kb],
                                                           isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                           a_rowptr.data(),
                                                           a_colind.data(),
                                                           normal,
                                                           a_ext,
                                                           da_ext,
                                                           grad_u_ext,
                                                           &u_grad_trial_trace[nSpace],
                                                           &u_trial_trace_ref.data()[ebN_local_kb*nSpace],
                                                           ebqe_penalty_ext.data()[ebNE_kb],
                                                           difffluxjacobian_ext);                 
                      
                      if (dflux_u_u_ext> 0.0)
                      { 
                        int ebN_local_kb_i = ebN_local_kb*nDOF_test_element+i;
                        for (int j=0;j<nDOF_trial_element;j++)
                          fluxTransport[j][i] += (dflux_u_u_ext + difffluxjacobian_ext)*
                            u_trial_trace_ref.data()[ebN_local_kb_i]*
                            u_test_dS[j];
                           
                      }
                      else
                        elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i])+
                                                ck.ExteriorElementBoundaryDiffusionAdjoint(isDOFBoundary_u.data()[ebNE_kb],
                                                                                            isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                                                            eb_adjoint_sigma,
                                                                                            u_ext,
                                                                                            bc_u_ext,
                                                                                            normal,
                                                                                            a_rowptr.data(),
                                                                                            a_colind.data(),
                                                                                            a_ext,
                                                                                            &u_grad_test_dS[i*nSpace]);                                                                                        
                    }                   
                }//i
              // local min/max at boundary
              min_u_bc_local = fmin(ebqe_u.data()[ebNE_kb], min_u_bc_local);
              max_u_bc_local = fmax(ebqe_u.data()[ebNE_kb], max_u_bc_local);
            }//kb
          //
          //update the element and global residual storage
          //
          for (int i=0;i<nDOF_test_element;i++)
            {
              int eN_i = eN*nDOF_test_element+i;
              int gi = offset_u+stride_u*u_l2g.data()[eN_i]; //global i-th index
              if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity or STABILIZATION_TYPE==STABILIZATION::SmoothnessIndicator or STABILIZATION_TYPE==STABILIZATION::Kuzmin)
                {
                  globalResidual.data()[gi] += dt*elementResidual_u[i];
                  boundary_integral[gi] += elementResidual_u[i];
                  min_u_bc[gi] = fmin(min_u_bc_local,min_u_bc[gi]);
                  max_u_bc[gi] = fmax(max_u_bc_local,max_u_bc[gi]);
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
                      TransportMatrix[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_eb_CellLoops.data()[ebN_i_j]]
                        += fluxTransport[i][j];
                      TransposeTransportMatrix[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_eb_CellLoops.data()[ebN_i_j]]
                        += fluxTransport[j][i];
                    }//j
                }
              else
                {
                  globalResidual.data()[offset_u+stride_u*r_l2g.data()[eN_i]] += elementResidual_u[i];
                }
            }//i
        }//ebNE
      if (STABILIZATION_TYPE==STABILIZATION::TaylorGalerkinEV)
        {
          meanEntropy /= meanOmega;
          double norm_factor = fmax(fabs(maxEntropy - meanEntropy), fabs(meanEntropy-minEntropy));
          for(int eN=0;eN<nElements_global;eN++)
            {
              double hK=elementDiameter.data()[eN]/degree_polynomial;
              double linear_viscosity = cMax*hK*maxVel[eN];
              double entropy_viscosity = cE*hK*hK*maxEntRes[eN]/norm_factor;
              for  (int k=0;k<nQuadraturePoints_element;k++)
                {
                  int eN_k = eN*nQuadraturePoints_element+k;
                  q_numDiff_u.data()[eN_k] = fmin(linear_viscosity,entropy_viscosity);
                }
            }
        }
      //edge based stabilization
      else if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity or 
               STABILIZATION_TYPE==STABILIZATION::SmoothnessIndicator or 
               STABILIZATION_TYPE==STABILIZATION::Kuzmin)
        {
          /////////////////////////////////////////////////////////////////
          // COMPUTE SMOOTHNESS INDICATOR and NORMALIZE ENTROPY RESIDUAL //
          /////////////////////////////////////////////////////////////////
          // NOTE: see NCLS.h for a different but equivalent implementation of this.
          //cek todo: can these loops over numDOFs be collapsed?
          int ij = 0;
          for (int i=0; i<numDOFs; i++)
            {
              double etaMaxi, etaMini;
              if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity) //EV
                {
                  // For eta min and max
                  etaMaxi = fabs(eta[i]);
                  etaMini = fabs(eta[i]);
                }
              // for smoothness indicator //
              double alpha_numerator = 0., alpha_denominator = 0.;
              for (int offset=csrRowIndeces_DofLoops.data()[i]; offset<csrRowIndeces_DofLoops.data()[i+1]; offset++)
                { // First loop in j (sparsity pattern)
                  int j = csrColumnOffsets_DofLoops.data()[offset];
                  if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity) //EV Stabilization
                    {
                      // COMPUTE ETA MIN AND ETA MAX //
                      etaMaxi = fmax(etaMaxi,fabs(eta[j]));
                      etaMini = fmin(etaMini,fabs(eta[j]));
                    }
                  alpha_numerator += u_dof_old.data()[i]-u_dof_old.data()[j];
                  alpha_denominator += fabs(u_dof_old.data()[i]-u_dof_old.data()[j]);
                  //update ij
                  ij+=1;
                }
              if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity) //EV Stab
                {
                  // Normalize entropy residual
                  global_entropy_residual[i] *= etaMini == etaMaxi ? 0. : 2*cE/(etaMaxi-etaMini);
                  quantDOFs[i] = fabs(global_entropy_residual[i]);
                }

              double alphai = alpha_numerator/(alpha_denominator+1E-15);
              quantDOFs[i] = alphai;


              if (POWER_SMOOTHNESS_INDICATOR==0)
                psi[i] = 1.0;
              else
                psi[i] = std::pow(alphai,POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
            }
          /////////////////////////////////////////////
          // ** LOOP IN DOFs FOR EDGE BASED TERMS ** //
          /////////////////////////////////////////////
          ij=0;
          for (int i=0; i<numDOFs; i++)
            {
              double solni = u_dof_old.data()[i]; // solution at time tn for the ith DOF
              double ith_dissipative_term = 0;
              double ith_low_order_dissipative_term = 0;
              double ith_flux_term = 0;
              double dLii = 0.;

              // loop over the sparsity pattern of the i-th DOF
              for (int offset=csrRowIndeces_DofLoops.data()[i]; offset<csrRowIndeces_DofLoops.data()[i+1]; offset++)
                {
                  int j = csrColumnOffsets_DofLoops.data()[offset];
                  double solnj = u_dof_old.data()[j]; // solution at time tn for the jth DOF
                  double dLowij, dLij, dEVij, dHij;

                  ith_flux_term += (TransportMatrix[ij]+ DiffusionMatrix[ij])*solnj;
                 
                  
                  if (i != j) //NOTE: there is really no need to check for i!=j (see formula for ith_dissipative_term)
                    {
                      // artificial compression
                      double solij = 0.5*(solni+solnj);
                      double Compij = cK*fmax(solij*(1.0-solij),0.0)/(fabs(solni-solnj)+1E-14);
                      // first-order dissipative operator
                      dLowij = fmax(fabs(TransportMatrix[ij]),fabs(TransposeTransportMatrix[ij]));
                      //std::cout << dLowij;
                      //dLij = fmax(0.,fmax(psi[i]*TransportMatrix[ij], // Approach by S. Badia
                      //              psi[j]*TransposeTransportMatrix[ij]));
                      dLij = dLowij*fmax(psi[i],psi[j]); // Approach by JLG & BP
                      
                      if (STABILIZATION_TYPE==STABILIZATION::EntropyViscosity) //EV Stab
                        {
                          // high-order (entropy viscosity) dissipative operator
                          dEVij = fmax(fabs(global_entropy_residual[i]),fabs(global_entropy_residual[j]));
                          dHij = fmin(dLowij,dEVij) * fmax(1.0-Compij,0.0); // artificial compression
                        }
                      else // smoothness based indicator
                        {
                          dHij = dLij * fmax(1.0-Compij,0.0); // artificial compression
                        }
                      //dissipative terms
                      ith_dissipative_term += dHij*(solnj-solni);
                      ith_low_order_dissipative_term += dLowij*(solnj-solni);
                      //dHij - dLij. This matrix is needed during FCT step
                      dt_times_dH_minus_dL[ij] = dt*(dHij - dLowij);
                      //std::cout << dLij;

                      dLii -= dLij;
                      dLow[ij] = dLowij;
                      
                    }
                  else //i==j
                    {
                      // NOTE: this is incorrect. Indeed, dLii = -sum_{j!=i}(dLij) and similarly for dCii.
                      // However, it is irrelevant since during the FCT step we do (dL-dC)*(solnj-solni)
                      dt_times_dH_minus_dL[ij]=0;
                      dLow[ij]=0;
                    }
                  //update ij
                  ij+=1;
                }
              double mi = ML.data()[i];
              // compute edge_based_cfl
              edge_based_cfl.data()[i] = 2.*fabs(dLii)/mi;
              // Debugging output to check values
              //std::cout << "Element: " << i << ", dLii: " << fabs(dLii) << ", mi: " << mi << std::endl;

              uLow[i] = u_dof_old.data()[i] - dt/mi*(ith_flux_term
                                                     + boundary_integral[i]
                                                     - ith_low_order_dissipative_term);

              // update residual
              if (LUMPED_MASS_MATRIX==1)
                globalResidual.data()[i] = u_dof_old.data()[i] - dt/mi*(ith_flux_term
                                                                        + boundary_integral[i]
                                                                        - ith_dissipative_term);
              else
                globalResidual.data()[i] += dt*(ith_flux_term - ith_dissipative_term);//cek todo: shouldn't this have boundaryIntegral?
            }//i
        }//edge-based
    }
  

    void calculateJacobian(arguments_dict& args)
    {
      xt::pyarray<double>& mesh_trial_ref = args.array<double>("mesh_trial_ref");
      xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
      xt::pyarray<double>& mesh_dof = args.array<double>("mesh_dof");
      xt::pyarray<double>& mesh_velocity_dof = args.array<double>("mesh_velocity_dof");
      double MOVING_DOMAIN = args.scalar<double>("MOVING_DOMAIN");
      xt::pyarray<int>& mesh_l2g = args.array<int>("mesh_l2g");
      xt::pyarray<double>& dV_ref = args.array<double>("dV_ref");
      xt::pyarray<double>& u_trial_ref = args.array<double>("u_trial_ref");
      xt::pyarray<double>& u_grad_trial_ref = args.array<double>("u_grad_trial_ref");
      xt::pyarray<double>& u_test_ref = args.array<double>("u_test_ref");
      xt::pyarray<double>& u_grad_test_ref = args.array<double>("u_grad_test_ref");
      xt::pyarray<double>& mesh_trial_trace_ref = args.array<double>("mesh_trial_trace_ref");
      xt::pyarray<double>& mesh_grad_trial_trace_ref = args.array<double>("mesh_grad_trial_trace_ref");
      xt::pyarray<double>& dS_ref = args.array<double>("dS_ref");
      xt::pyarray<double>& u_trial_trace_ref = args.array<double>("u_trial_trace_ref");
      xt::pyarray<double>& u_grad_trial_trace_ref = args.array<double>("u_grad_trial_trace_ref");
      xt::pyarray<double>& u_test_trace_ref = args.array<double>("u_test_trace_ref");
      xt::pyarray<double>& u_grad_test_trace_ref = args.array<double>("u_grad_test_trace_ref");
      xt::pyarray<double>& normal_ref = args.array<double>("normal_ref");
      xt::pyarray<double>& boundaryJac_ref = args.array<double>("boundaryJac_ref");
      int nElements_global = args.scalar<int>("nElements_global");
      double useMetrics = args.scalar<double>("useMetrics");
      double alphaBDF = args.scalar<double>("alphaBDF");
      int lag_shockCapturing = args.scalar<int>("lag_shockCapturing");
      double shockCapturingDiffusion = args.scalar<double>("shockCapturingDiffusion");
      xt::pyarray<int>& u_l2g = args.array<int>("u_l2g");
      xt::pyarray<int>& r_l2g = args.array<int>("r_l2g");
      xt::pyarray<double>& elementDiameter = args.array<double>("elementDiameter");
      xt::pyarray<double>& u_dof = args.array<double>("u_dof");
      xt::pyarray<double>& velocity = args.array<double>("velocity");
      xt::pyarray<double>& q_m_betaBDF = args.array<double>("q_m_betaBDF");
      xt::pyarray<double>& cfl = args.array<double>("cfl");
      xt::pyarray<double>& q_numDiff_u_last = args.array<double>("q_numDiff_u_last");
      xt::pyarray<int>& csrRowIndeces_u_u = args.array<int>("csrRowIndeces_u_u");
      xt::pyarray<int>& csrColumnOffsets_u_u = args.array<int>("csrColumnOffsets_u_u");
      xt::pyarray<double>& globalJacobian = args.array<double>("globalJacobian");
      int nExteriorElementBoundaries_global = args.scalar<int>("nExteriorElementBoundaries_global");
      xt::pyarray<int>& exteriorElementBoundariesArray = args.array<int>("exteriorElementBoundariesArray");
      xt::pyarray<int>& elementBoundaryElementsArray = args.array<int>("elementBoundaryElementsArray");
      xt::pyarray<int>& elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
      xt::pyarray<double>& ebqe_velocity_ext = args.array<double>("ebqe_velocity_ext");
      xt::pyarray<int>& isDOFBoundary_u = args.array<int>("isDOFBoundary_u");
      xt::pyarray<double>& ebqe_bc_u_ext = args.array<double>("ebqe_bc_u_ext");
      xt::pyarray<int>& isFluxBoundary_u = args.array<int>("isFluxBoundary_u");
      xt::pyarray<double>& ebqe_bc_flux_u_ext = args.array<double>("ebqe_bc_flux_u_ext");
      xt::pyarray<int>& csrColumnOffsets_eb_u_u = args.array<int>("csrColumnOffsets_eb_u_u");
      STABILIZATION STABILIZATION_TYPE = static_cast<STABILIZATION>(args.scalar<int>("STABILIZATION_TYPE"));
//      ENTROPY ENTROPY_TYPE = static_cast<ENTROPY>(args.scalar<int>("ENTROPY_TYPE"));    
//      STABILIZATION STABILIZATION_TYPE{args.scalar<int>("STABILIZATION_TYPE")};
      double physicalDiffusion = args.scalar<double>("physicalDiffusion");
      xt::pyarray<double>& q_a = args.array<double>("q_a");
      xt::pyarray<double>& ebqe_a = args.array<double>("ebq_a");
      
      // Extract alpha_L, alpha_T, and Dm from args
      const double alphaL_val = args.scalar<double>("alpha_L"); // Longitudinal dispersion coefficient
      const double alphaT_val = args.scalar<double>("alpha_T"); // Transverse dispersion coefficient
      const double Dm_val = args.scalar<double>("Dm");         // Molecular diffusion coefficient


      double Ct_sge = 4.0;



      /////////////////////////////////////////////////////////////////////////
      xt::pyarray<int>& a_rowptr = args.array<int>("a_rowptr");
      xt::pyarray<int>& a_colind = args.array<int>("a_colind");
      //xt::pyarray<double>& D = args.array<double>("D");
      //////////////////////////////////////////////////////////////////////////
      xt::pyarray<int>& isDiffusiveFluxBoundary_u = args.array<int>("isDiffusiveFluxBoundary_u");
      xt::pyarray<double>& ebqe_penalty_ext = args.array<double>("ebqe_penalty_ext");
      

      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
        {
          double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
          for (int i=0;i<nDOF_test_element;i++)
            for (int j=0;j<nDOF_trial_element;j++)
              {
                elementJacobian_u_u[i][j]=0.0;
              }
          for  (int k=0;k<nQuadraturePoints_element;k++)
            {
              int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
                eN_k_nSpace = eN_k*nSpace,
                eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

              //declare local storage
              double u=0.0,
                grad_u[nSpace],
                m=0.0,dm=0.0,
                f[nSpace],df[nSpace],
                a[nnz],da[nnz],
                
                m_t=0.0,dm_t=0.0,
                dpdeResidual_u_u[nDOF_trial_element],
                Lstar_u[nDOF_test_element],
                dsubgridError_u_u[nDOF_trial_element],
                tau=0.0,tau0=0.0,tau1=0.0,
                jac[nSpace*nSpace],
                jacDet,
                jacInv[nSpace*nSpace],
                u_grad_trial[nDOF_trial_element*nSpace],
                dV,
                u_test_dV[nDOF_test_element],
                u_grad_test_dV[nDOF_test_element*nSpace],
                x,y,z,xt,yt,zt,
                G[nSpace*nSpace],G_dd_G,tr_G;
              //
              //calculate solution and gradients at quadrature points
              //
              //get jacobian, etc for mapping reference element
              ck.calculateMapping_element(eN,
                                          k,
                                          mesh_dof.data(),
                                          mesh_l2g.data(),
                                          mesh_trial_ref.data(),
                                          mesh_grad_trial_ref.data(),
                                          jac,
                                          jacDet,
                                          jacInv,
                                          x,y,z);
              ck.calculateMappingVelocity_element(eN,
                                                  k,
                                                  mesh_velocity_dof.data(),
                                                  mesh_l2g.data(),
                                                  mesh_trial_ref.data(),
                                                  xt,yt,zt);
              //get the physical integration weight
              dV = fabs(jacDet)*dV_ref.data()[k];
              ck.calculateG(jacInv,G,G_dd_G,tr_G);
              //get the trial function gradients
              ck.gradTrialFromRef(&u_grad_trial_ref.data()[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
              //get the solution
              ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_ref.data()[k*nDOF_trial_element],u);
              //get the solution gradients
              ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial,grad_u);
              //precalculate test function products with integration weights
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  u_test_dV[j] = u_test_ref.data()[k*nDOF_trial_element+j]*dV;
                  for (int I=0;I<nSpace;I++)
                    {
                      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
                    }
                }
              //
              //calculate pde coefficients and derivatives at quadrature points
              //
              
              //const double* q_a_ptr = &q_a[eN_k * a_rowptr[nSpace]];
              evaluateCoefficients(a_rowptr.data(),
				                           a_colind.data(),
                                   &velocity.data()[eN_k_nSpace],
                                   alphaL_val, alphaT_val, Dm_val, 
                                   //q_a_ptr, //&q_a.data()[eN * nQuadraturePoints_element * nnz + k * nnz],//[eN_k*a_rowptr.data()[nSpace]],//[eN_k*nnz],
                                   u,
                                   m,
                                   dm,
                                   f,
                                   df,
                                   a,
                                   da);
              //
              //moving mesh
              //
              double mesh_velocity[3];
              mesh_velocity[0] = xt;
              mesh_velocity[1] = yt;
              mesh_velocity[2] = zt;
            
              for(int I=0;I<nSpace;I++)
                {
                  f[I] -= MOVING_DOMAIN*m*mesh_velocity[I];
                  df[I] -= MOVING_DOMAIN*dm*mesh_velocity[I];
                }
              //
              //calculate time derivatives
              //
              ck.bdf(alphaBDF,
                     q_m_betaBDF.data()[eN_k],//since m_t isn't used, we don't have to correct mass
                     m,
                     dm,
                     m_t,
                     dm_t);
              if (STABILIZATION_TYPE == STABILIZATION::VMS)
                {
                  //
                  //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
                  //
                  //calculate the adjoint times the test functions
                  for (int i=0;i<nDOF_test_element;i++)
                    {
                      int i_nSpace = i*nSpace;
                      Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
                    }
                  //calculate the Jacobian of strong residual
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      int j_nSpace = j*nSpace;
                      dpdeResidual_u_u[j]= ck.MassJacobian_strong(dm_t,u_trial_ref.data()[k*nDOF_trial_element+j]) +
                        ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);
                    }
                  //tau and tau*Res
                  calculateSubgridError_tau(elementDiameter.data()[eN],
                                            dm_t,
                                            df,
                                            cfl.data()[eN_k],
                                            tau0);

                  calculateSubgridError_tau(Ct_sge,
                                            G,
                                            dm_t,
                                            df,
                                            tau1,
                                            cfl.data()[eN_k]);
                  tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

                  for(int j=0;j<nDOF_trial_element;j++)
                    dsubgridError_u_u[j] = -tau*dpdeResidual_u_u[j];
                }
              for(int i=0;i<nDOF_test_element;i++)
                {
                  for(int j=0;j<nDOF_trial_element;j++)
                    {
                      int j_nSpace = j*nSpace;
                      int i_nSpace = i*nSpace;
                      if (STABILIZATION_TYPE==STABILIZATION::Galerkin)
                        {
                          elementJacobian_u_u[i][j] +=
                            ck.MassJacobian_weak(dm_t,
                                                 u_trial_ref.data()[k*nDOF_trial_element+j],
                                                 u_test_dV[i]) +
                            ck.AdvectionJacobian_weak(df,
                                                      u_trial_ref.data()[k*nDOF_trial_element+j],
                                                      &u_grad_test_dV[i_nSpace]) +
                            ck.DiffusionJacobian_weak(a_rowptr.data(),a_colind.data(),a,da,
						                                          grad_u,&u_grad_test_dV[i_nSpace],1.0,
						                                          u_trial_ref.data()[k*nDOF_trial_element+j],&u_grad_trial[j_nSpace])
                                                      +
                            ck.NumericalDiffusionJacobian(physicalDiffusion,
                                                          &u_grad_trial[j_nSpace],
                                                          &u_grad_test_dV[i_nSpace]); //implicit
                        }
                      else if (STABILIZATION_TYPE==STABILIZATION::VMS)
                        {
                          elementJacobian_u_u[i][j] +=
                            ck.MassJacobian_weak(dm_t,
                                                 u_trial_ref.data()[k*nDOF_trial_element+j],
                                                 u_test_dV[i]) +
                            ck.AdvectionJacobian_weak(df,
                                                      u_trial_ref.data()[k*nDOF_trial_element+j],
                                                      &u_grad_test_dV[i_nSpace]) +
                            ck.DiffusionJacobian_weak(a_rowptr.data(),a_colind.data(),a,da,
                                                      grad_u,&u_grad_test_dV[i_nSpace],1.0,
                                                      u_trial_ref.data()[k*nDOF_trial_element+j],&u_grad_trial[j_nSpace])+

                            ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) +
                            ck.NumericalDiffusionJacobian(q_numDiff_u_last.data()[eN_k] + physicalDiffusion,
                                                          &u_grad_trial[j_nSpace],
                                                          &u_grad_test_dV[i_nSpace]); //implicit
                        }
                      else if (STABILIZATION_TYPE==STABILIZATION::TaylorGalerkinEV or 
                               STABILIZATION_TYPE==STABILIZATION::EntropyViscosity or
                               STABILIZATION_TYPE==STABILIZATION::SmoothnessIndicator or 
                               STABILIZATION_TYPE==STABILIZATION::Kuzmin)
                        {
                          elementJacobian_u_u[i][j] +=
                            ck.MassJacobian_weak(1.0,
                                                 u_trial_ref.data()[k*nDOF_trial_element+j],
                                                 u_test_dV[i]);
                        }
                    }//j
                }//i
            }//k
          //
          //load into element Jacobian into global Jacobian
          //
          for (int i=0;i<nDOF_test_element;i++)
            {
              int eN_i = eN*nDOF_test_element+i;
              for (int j=0;j<nDOF_trial_element;j++)
                {
                  int eN_i_j = eN_i*nDOF_trial_element+j;
                  globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += elementJacobian_u_u[i][j];
                }//j
            }//i
        }//elements
      //
      //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
      //
      if (STABILIZATION_TYPE==STABILIZATION::VMS or STABILIZATION_TYPE==STABILIZATION::Galerkin)
        {
          for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
            {
              int ebN = exteriorElementBoundariesArray.data()[ebNE];
              int eN  = elementBoundaryElementsArray.data()[ebN*2+0],
                ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN*2+0],
                eN_nDOF_trial_element = eN*nDOF_trial_element;
              double fluxJacobian_u_u[nDOF_test_element][nDOF_trial_element];
              for (int i=0;i<nDOF_test_element;i++)
                for (int j=0;j<nDOF_trial_element;j++)
                  {
                    fluxJacobian_u_u[i][j]=0.0;
                  }
              for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
                {
                  int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
                    ebNE_kb_nSpace = ebNE_kb*nSpace,
                    ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
                    ebN_local_kb_nSpace = ebN_local_kb*nSpace;
                  double u_ext=0.0,
                    grad_u_ext[nSpace],
                    m_ext=0.0,
                    dm_ext=0.0,
                    f_ext[nSpace],
                    df_ext[nSpace],

                    a_ext[nnz],
                    da_ext[nnz],
                    bc_a_ext[nnz],
                    bc_da_ext[nnz],

                    dflux_u_u_ext=0.0,
                    difffluxjacobian_ext=0.0,
                    bc_u_ext=0.0,
                    //bc_grad_u_ext[nSpace],
                    bc_m_ext=0.0,
                    bc_dm_ext=0.0,
                    bc_f_ext[nSpace],
                    bc_df_ext[nSpace],
                    //////////////
                    diffusiveFluxJacobian_u_u[nDOF_trial_element],
                    ////////////////////////////
                    jac_ext[nSpace*nSpace],
                    jacDet_ext,
                    jacInv_ext[nSpace*nSpace],
                    boundaryJac[nSpace*(nSpace-1)],
                    metricTensor[(nSpace-1)*(nSpace-1)],
                    metricTensorDetSqrt,
                    dS,
                    u_test_dS[nDOF_test_element],
                    u_grad_trial_trace[nDOF_trial_element*nSpace],
                    u_grad_test_dS[nDOF_trial_element*nSpace],
                    normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
                    //
                    G[nSpace*nSpace],G_dd_G,tr_G;

                  //
                  //calculate the solution and gradients at quadrature points
                  //
                  ck.calculateMapping_elementBoundary(eN,
                                                      ebN_local,
                                                      kb,
                                                      ebN_local_kb,
                                                      mesh_dof.data(),
                                                      mesh_l2g.data(),
                                                      mesh_trial_trace_ref.data(),
                                                      mesh_grad_trial_trace_ref.data(),
                                                      boundaryJac_ref.data(),
                                                      jac_ext,
                                                      jacDet_ext,
                                                      jacInv_ext,
                                                      boundaryJac,
                                                      metricTensor,
                                                      metricTensorDetSqrt,
                                                      normal_ref.data(),
                                                      normal,
                                                      x_ext,y_ext,z_ext);
                  ck.calculateMappingVelocity_elementBoundary(eN,
                                                              ebN_local,
                                                              kb,
                                                              ebN_local_kb,
                                                              mesh_velocity_dof.data(),
                                                              mesh_l2g.data(),
                                                              mesh_trial_trace_ref.data(),
                                                              xt_ext,yt_ext,zt_ext,
                                                              normal,
                                                              boundaryJac,
                                                              metricTensor,
                                                              integralScaling);
                  dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref.data()[kb];
                  ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
                  //compute shape and solution information
                  //shape
                  ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
                  //solution and gradients
                  ck.valFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],&u_trial_trace_ref.data()[ebN_local_kb*nDOF_test_element],u_ext);
                  ck.gradFromDOF(u_dof.data(),&u_l2g.data()[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
                  //precalculate test function products with integration weights
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb*nDOF_test_element+j]*dS;
                      for (int I=0;I<nSpace;I++)
                      {
                        u_grad_test_dS[j*nSpace+I]= u_grad_trial_trace[j*nSpace+I]*dS;
                      }
                    }
                  //
                  //load the boundary values
                  //
                  bc_u_ext = isDOFBoundary_u.data()[ebNE_kb]*ebqe_bc_u_ext.data()[ebNE_kb]+(1-isDOFBoundary_u.data()[ebNE_kb])*u_ext;
                  //
                  //
                  //calculate the internal and external trace of the pde coefficients
                  //

                  const double* qb_a_ptr = &ebqe_a[ebNE_kb * a_rowptr[nSpace]];
                  evaluateCoefficients(a_rowptr.data(),
                                       a_colind.data(),
                                       &ebqe_velocity_ext.data()[ebNE_kb_nSpace],
                                       alphaL_val, alphaT_val, Dm_val, 
                                       //qb_a_ptr,//&q_a.data()[ebNE * nQuadraturePoints_element * nnz + kb * nnz],//[ebNE_kb*a_rowptr.data()[nSpace]],//[ebNE_kb* nnz],
                                       u_ext,
                                       m_ext,
                                       dm_ext,
                                       f_ext,
                                       df_ext,
                                       a_ext,
                                       da_ext
                                       );

                  evaluateCoefficients(a_rowptr.data(),
                                       a_colind.data(),
                                       &ebqe_velocity_ext.data()[ebNE_kb_nSpace],
                                       alphaL_val, alphaT_val, Dm_val, 
                                       //qb_a_ptr,//&q_a.data()[ebNE * nQuadraturePoints_element * nnz + kb * nnz],//[ebNE_kb*a_rowptr.data()[nSpace]],//[ebNE_kb* nnz],
                                       bc_u_ext,
                                       bc_m_ext,
                                       bc_dm_ext,
                                       bc_f_ext,
                                       bc_df_ext,
                                       bc_a_ext,
                                       bc_da_ext);
                  //
                  //moving domain
                  //
                  double mesh_velocity[3];
                  mesh_velocity[0] = xt_ext;
                  mesh_velocity[1] = yt_ext;
                  mesh_velocity[2] = zt_ext;
                  for (int I=0;I<nSpace;I++)
                    {
                      f_ext[I] -= MOVING_DOMAIN*m_ext*mesh_velocity[I];
                      df_ext[I] -= MOVING_DOMAIN*dm_ext*mesh_velocity[I];
                      bc_f_ext[I] -= MOVING_DOMAIN*bc_m_ext*mesh_velocity[I];
                      bc_df_ext[I] -= MOVING_DOMAIN*bc_dm_ext*mesh_velocity[I];
                    }
                  //
                  //calculate the numerical fluxes
                  //
                  exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u.data()[ebNE_kb],
                                                           isFluxBoundary_u.data()[ebNE_kb],
                                                           normal,
                                                           df_ext,
                                                           dflux_u_u_ext);
                  exteriorNumericalDiffusiveFluxDerivative(isDOFBoundary_u.data()[ebNE_kb],
                                                           isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                           a_rowptr.data(),
                                                           a_colind.data(),
                                                           normal,
                                                           a_ext,
                                                           da_ext,
                                                           grad_u_ext,
                                                           &u_grad_trial_trace[nSpace],
                                                           &u_trial_trace_ref.data()[ebN_local_kb*nSpace],
                                                           ebqe_penalty_ext.data()[ebNE_kb],
                                                           difffluxjacobian_ext);
                  
                  
                  //
                  //calculate the flux jacobian
                  //
                  for (int i=0;i<nDOF_test_element;i++)
                    for (int j=0;j<nDOF_trial_element;j++)
                      {
                        int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
  
                        fluxJacobian_u_u[i][j]+=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref.data()[ebN_local_kb_j])*u_test_dS[i]+
                                                ExteriorNumericalDiffusiveFluxJacobian(a_rowptr.data(),
                                                                                       a_colind.data(),
                                                                                       isDOFBoundary_u.data()[ebNE_kb],
                                                                                       isDiffusiveFluxBoundary_u.data()[ebNE_kb],
                                                                                       normal,
                                                                                       a_ext,
                                                                                       u_trial_trace_ref.data()[ebN_local_kb_j],
                                                                                       &u_grad_trial_trace[j*nSpace],
                                                                                       ebqe_penalty_ext.data()[ebNE_kb])*u_test_dS[i];
                      }//j

              //
              //update the global Jacobian from the flux Jacobian
              //
              for (int i=0;i<nDOF_test_element;i++)
                {
                  int eN_i = eN*nDOF_test_element+i;
                  for (int j=0;j<nDOF_trial_element;j++)
                    {
                      int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
                      globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_eb_u_u.data()[ebN_i_j]] += fluxJacobian_u_u[i][j];
//                                                                                                                         
                    }//j
                }//i
              }//kb
            }//ebNE
        }//VMS and Galerkin
    }//computeJacobian

  void FCTStep(arguments_dict& args)
  {
    double dt = args.scalar<double>("dt");
    int NNZ = args.scalar<int>("NNZ");
    int numDOFs = args.scalar<int>("numDOFs");
    xt::pyarray<double>& lumped_mass_matrix = args.array<double>("lumped_mass_matrix");
    xt::pyarray<double>& soln = args.array<double>("soln");
    xt::pyarray<double>& solH = args.array<double>("solH");
    xt::pyarray<double>& uLow = args.array<double>("uLow");
    xt::pyarray<double>& dLow = args.array<double>("dLow");
    xt::pyarray<double>& limited_solution = args.array<double>("limited_solution");
    xt::pyarray<int>& csrRowIndeces_DofLoops = args.array<int>("csrRowIndeces_DofLoops");
    xt::pyarray<int>& csrColumnOffsets_DofLoops = args.array<int>("csrColumnOffsets_DofLoops");
    xt::pyarray<double>& MassMatrix = args.array<double>("MassMatrix");
    xt::pyarray<double>& dt_times_dH_minus_dL = args.array<double>("dt_times_dH_minus_dL");
    xt::pyarray<double>& min_u_bc = args.array<double>("min_u_bc");
    xt::pyarray<double>& max_u_bc = args.array<double>("max_u_bc");
    int LUMPED_MASS_MATRIX = args.scalar<int>("LUMPED_MASS_MATRIX");
//    STABILIZATION STABILIZATION_TYPE{args.scalar<int>("STABILIZATION_TYPE")};
    STABILIZATION STABILIZATION_TYPE = static_cast<STABILIZATION>(args.scalar<int>("STABILIZATION_TYPE"));
//      ENTROPY ENTROPY_TYPE = static_cast<ENTROPY>(args.scalar<int>("ENTROPY_TYPE"));    
    Rpos.resize(numDOFs,0.0);
    Rneg.resize(numDOFs,0.0);
    FluxCorrectionMatrix.resize(NNZ,0.0);
    int ij=0;
    //loop over nodes (i)
    for (int i=0; i<numDOFs; i++)
      {
        //read some vectors
        double solHi = solH.data()[i];
        double solni = soln.data()[i];
        double mi = lumped_mass_matrix.data()[i];
        double uLowi = uLow.data()[i];
        double uDotLowi = (uLowi - solni)/dt;
        double mini=min_u_bc.data()[i], maxi=max_u_bc.data()[i]; // init min/max with value at BCs (NOTE: if no boundary then min=1E10, max=-1E10)
        double Pposi=0, Pnegi=0;
        // Loop over neighbors (j)
        for (int offset=csrRowIndeces_DofLoops.data()[i]; offset<csrRowIndeces_DofLoops.data()[i+1]; offset++)
          {
            assert(offset == ij); // (CSR matrix consistency)
            int j = csrColumnOffsets_DofLoops.data()[offset];
            double solnj = soln.data()[j];
            ////////////////////////
            // COMPUTE THE BOUNDS //
            ////////////////////////
            mini = fmin(mini,solnj);
            maxi = fmax(maxi,solnj);
            double uLowj = uLow.data()[j];
            double uDotLowj = (uLowj - solnj)/dt;
            // i-th row of flux correction matrix
            if (STABILIZATION_TYPE == STABILIZATION::Kuzmin)
              {
                FluxCorrectionMatrix[ij] = dt*(MassMatrix.data()[ij]*(uDotLowi-uDotLowj)
                                               + dLow.data()[ij]*(uLowi-uLowj));
              }
            else
              {
                double ML_minus_MC =
                  (LUMPED_MASS_MATRIX == 1 ? 0. : (i==j ? 1. : 0.)*mi - MassMatrix.data()[ij]);
                FluxCorrectionMatrix[ij] = ML_minus_MC * (solH.data()[j]-solnj - (solHi-solni))
                  + dt_times_dH_minus_dL.data()[ij]*(solnj-solni);
              }
            Pposi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] > 0) ? 1. : 0.);
            Pnegi += FluxCorrectionMatrix[ij]*((FluxCorrectionMatrix[ij] < 0) ? 1. : 0.);
            ij+=1;
          }//j
        double Qposi = mi*(maxi-uLow.data()[i]);
        double Qnegi = mi*(mini-uLow.data()[i]);
        Rpos[i] = ((Pposi==0) ? 1. : fmin(1.0,Qposi/Pposi));
        Rneg[i] = ((Pnegi==0) ? 1. : fmin(1.0,Qnegi/Pnegi));
      }//i
    ij=0;
    for (int i=0; i<numDOFs; i++)
      {
        double ith_Limiter_times_FluxCorrectionMatrix = 0.;
        double Rposi = Rpos[i], Rnegi = Rneg[i];
        for (int offset=csrRowIndeces_DofLoops.data()[i]; offset<csrRowIndeces_DofLoops.data()[i+1]; offset++)
          {
            assert(offset == ij); // (CSR matrix consistency
            int j = csrColumnOffsets_DofLoops.data()[offset];
            double Lij = 1;
            Lij = ((FluxCorrectionMatrix[ij]>0) ? fmin(Rposi,Rneg[j]) : fmin(Rnegi,Rpos[j]));
            ith_Limiter_times_FluxCorrectionMatrix += Lij * FluxCorrectionMatrix[ij];
            ij+=1;
          }
        limited_solution.data()[i] = uLow.data()[i] + 1./lumped_mass_matrix.data()[i]*ith_Limiter_times_FluxCorrectionMatrix;
      }
    }//FCTStep


    // NOTE: expects q_x to be a C-contiguous buffer with trailing dimension == nSpace,
// e.g. Python side: q['x'].shape == (nElements_global*nQuadraturePoints_element, nSpace)
// or (nElements_global, nQuadraturePoints_element, nSpace) with order='C'.



void Update_concentration_RWPT(arguments_dict& args)
{
  const double dt = args.scalar<double>("dt");

  // --- mesh / FE ---
  xt::pyarray<double>& mesh_trial_ref      = args.array<double>("mesh_trial_ref");
  xt::pyarray<double>& mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
  xt::pyarray<double>& mesh_dof            = args.array<double>("mesh_dof");
  xt::pyarray<int>&    mesh_l2g            = args.array<int>("mesh_l2g");
  xt::pyarray<double>& dV_ref              = args.array<double>("dV_ref");
  xt::pyarray<double>& u_trial_ref         = args.array<double>("u_trial_ref");
  xt::pyarray<double>& u_grad_trial_ref    = args.array<double>("u_grad_trial_ref");
  xt::pyarray<double>& u_dof               = args.array<double>("u_dof");
  xt::pyarray<int>&    u_l2g               = args.array<int>("u_l2g");
  const int            nElements_global    = args.scalar<int>("nElements_global");

  // --- qp fields ---
  xt::pyarray<double>& q_u        = args.array<double>("q_u");        // concentration C at qp
  xt::pyarray<double>& q_porosity = args.array<double>("q_porosity"); // porosity at qp
  xt::pyarray<double>& velocity   = args.array<double>("velocity");   // velocity at qp

  // --- optional inputs for grad porosity (currently unused, grad_phi set to 0) ---
  xt::pyarray<int>&    phi_l2g    = args.array<int>("phi_l2g");

  // --- parameters ---
  const double mass_per_particle = args.scalar<double>("mass_per_particle");
  const double alphaL_val        = args.scalar<double>("alpha_L");
  const double alphaT_val        = args.scalar<double>("alpha_T");
  const double Dm_val            = args.scalar<double>("Dm");

  // --- per-qp outputs (pre-allocated) ---
  xt::pyarray<double>& qp_dV    = args.array<double>("qp_dV");
  xt::pyarray<double>& qp_mass  = args.array<double>("qp_mass");
  xt::pyarray<int>&    qp_spawn = args.array<int>("qp_spawn");
  xt::pyarray<double>& q_x      = args.array<double>("q_x");   // flattened in C order: [nQP_total * nSpace]

  // --- projection data (pre-allocated) ---
  xt::pyarray<double>& lumped_mass_matrix = args.array<double>("lumped_mass_matrix"); // ML (size: numDOFs)
  xt::pyarray<double>& u_new              = args.array<double>("u_new");              // accumulator (size: numDOFs)
  const int            numDOFs            = args.scalar<int>("numDOFs");

  // ---------- DEBUG: input u_dof stats ----------
  {
    double umin = 1e300, umax = -1e300;
    for (int I = 0; I < numDOFs; ++I) {
      const double v = u_dof.data()[I];
      if (v < umin) umin = v;
      if (v > umax) umax = v;
    }
    const int showN = (numDOFs < 5 ? numDOFs : 5);
    std::cerr << "[RWPT] u_dof(in) min=" << umin << " max=" << umax << " head:";
    for (int i = 0; i < showN; ++i) std::cerr << " " << u_dof.data()[i];
    std::cerr << " | mpp=" << mass_per_particle << "\n";
  }

  // clear qp outputs
  std::fill(qp_dV.begin(),    qp_dV.end(),    0.0);
  std::fill(qp_mass.begin(),  qp_mass.end(),  0.0);
  std::fill(qp_spawn.begin(), qp_spawn.end(), 0);

  // ---------------- PASS 1: compute qp stats & spawn counts, record qp coords ----------------
  long long N_total = 0; // total particles to spawn this step

  // DEBUG accumulators for PASS 1
  long long qps_with_spawn = 0;
  int       max_n_particles = 0;
  double    max_u_qp = 0.0, max_m_qp_init = 0.0, sum_m_qp_init = 0.0;
  double    lambda_sum = 0.0;          // new: sum of all λ over qps
  long long qps_ge1    = 0;            // new: qps where λ ≥ 1
  int       printed_qp = 0;            // print only a few qp samples to avoid spam
  const int PRINT_QP   = 5;

  for (int eN = 0; eN < nElements_global; ++eN)
  {
    for (int k = 0; k < nQuadraturePoints_element; ++k)
    {
      const int eN_k                  = eN * nQuadraturePoints_element + k;
      const int eN_k_nSpace           = eN_k * nSpace;
      const int eN_nDOF_trial_element = eN * nDOF_trial_element;

      // element mapping + qp coordinates
      double jac[nSpace*nSpace], jacInv[nSpace*nSpace], jacDet, x=0.0, y=0.0, z=0.0;
      ck.calculateMapping_element(eN, k,
                                  mesh_dof.data(), mesh_l2g.data(),
                                  mesh_trial_ref.data(), mesh_grad_trial_ref.data(),
                                  jac, jacDet, jacInv, x,y,z);

      // store qp coordinates in q_x (C order, trailing dim = nSpace)
      const int qFlat = eN*nQuadraturePoints_element + k;
      double* qx = &q_x.data()[ qFlat * nSpace ];
      qx[0] = x;
      if (nSpace > 1) qx[1] = y;
      if (nSpace > 2) qx[2] = z;

      // qp volume weight: ensure positive (dV_ref may carry sign on some refs)
      const double dV = std::fabs(jacDet * dV_ref.data()[k]);
      qp_dV.data()[eN_k] = dV;

      // interpolate concentration at qp
      double u_qp = 0.0;
      ck.valFromDOF(u_dof.data(),
                    &u_l2g.data()[eN_nDOF_trial_element],
                    &u_trial_ref.data()[k*nDOF_trial_element],
                    u_qp);

      // include porosity in initial mass (consistent with PASS 4)
      const double phi        = q_porosity.data()[eN_k];
      const double m_qp_init  = u_qp * phi * dV;   // mass at this qp
      qp_mass.data()[eN_k]    = m_qp_init;

      // particles to spawn at this qp
      const double lambda = (mass_per_particle > 0.0) ? (m_qp_init / mass_per_particle) : 0.0;
      lambda_sum += lambda;
      if (lambda >= 1.0) qps_ge1++;

      // ---- MODIFIED: stochastic rounding instead of pure floor ----
      const double base = std::floor(std::max(lambda, 0.0));
      const double frac = std::max(lambda - base, 0.0);
      const double u01  = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX); // [0,1)
      const int    n_particles = static_cast<int>(base + (u01 < frac ? 1 : 0));
      // -------------------------------------------------------------

      qp_spawn.data()[eN_k] = n_particles;
      N_total              += static_cast<long long>(n_particles);

      // ---------- DEBUG: track and sample-print ----------
      if (n_particles > 0) qps_with_spawn++;
      if (n_particles > max_n_particles) max_n_particles = n_particles;
      if (u_qp > max_u_qp) max_u_qp = u_qp;
      if (m_qp_init > max_m_qp_init) max_m_qp_init = m_qp_init;
      sum_m_qp_init += m_qp_init;

      // print first few qp details (element 0 preferred) to see n_particles, u_qp, dV, m_qp_init
      if (printed_qp < PRINT_QP && eN == 0) {
        std::cout << "[RWPT] PASS1 sample qp #" << printed_qp
                  << " (e=" << eN << ", k=" << k << "): "
                  << "u_qp=" << u_qp << " dV=" << dV
                  << " phi=" << phi << " m_qp_init=" << m_qp_init
                  << " lambda=" << lambda << " n_particles=" << n_particles << "\n";
        printed_qp++;
      }
    }
  }

  // ---------- DEBUG: PASS 1 summary ----------
  std::cout << "[RWPT] PASS1 summary: "
            << "N_total=" << N_total
            << " qps_with_spawn=" << qps_with_spawn
            << " qps_ge1=" << qps_ge1
            << " lambda_sum=" << lambda_sum
            << " max(n_particles)=" << max_n_particles
            << " max(u_qp)=" << max_u_qp
            << " max(m_qp_init)=" << max_m_qp_init
            << " sum(m_qp_init)=" << sum_m_qp_init << "\n";

  // ---------------- PASS 2: spawn & move particles (local, no args for p_X/p_qp) --------------
  const long long N_cap = (N_total > 0) ? N_total : 0;
  std::vector<double> p_X(static_cast<size_t>(N_cap) * nSpace, 0.0); // positions after move
  std::vector<int>    p_qp(static_cast<size_t>(N_cap), 0);           // source qp index
  int pid = 0;

  for (int eN = 0; eN < nElements_global; ++eN)
  {
    for (int k = 0; k < nQuadraturePoints_element; ++k)
    {
      const int eN_k        = eN * nQuadraturePoints_element + k;
      const int eN_k_nSpace = eN_k * nSpace;

      const int n_particles = qp_spawn.data()[eN_k];
      if (n_particles <= 0) continue;

      // qp coordinates
      double xq[3] = {0.0, 0.0, 0.0};
      xq[0] = q_x.data()[eN_k_nSpace + 0];
      if (nSpace > 1) xq[1] = q_x.data()[eN_k_nSpace + 1];
      if (nSpace > 2) xq[2] = q_x.data()[eN_k_nSpace + 2];

      // qp velocity
      double v_qp[3] = {0.0, 0.0, 0.0};
      for (int I = 0; I < nSpace; ++I)
        v_qp[I] = velocity.data()[eN_k_nSpace + I];

      // porosity and (currently zero) grad_phi + divD
      const double phi = q_porosity.data()[eN_k];
      double grad_phi[3] = {0.0, 0.0, 0.0};
      double divD[3]     = {0.0, 0.0, 0.0};

      for (int p = 0; p < n_particles; ++p)
      {
        // draw Z ~ N(0, I)
        double Z[3] = {0.0, 0.0, 0.0};
        for (int I = 0; I < nSpace; ++I)
          Z[I] = normal_random();

        // displacement over dt
        double dxp[3] = {0.0, 0.0, 0.0};
        displacement(v_qp, grad_phi, divD, phi,
                     alphaL_val, alphaT_val, Dm_val,
                     dt, Z, dxp);

        // new particle position
        double xp[3] = {xq[0], xq[1], xq[2]};
        for (int I = 0; I < nSpace; ++I)
          xp[I] += dxp[I];

        // store
        for (int I = 0; I < nSpace; ++I)
          p_X[static_cast<size_t>(pid)*nSpace + I] = xp[I];
        p_qp[static_cast<size_t>(pid)] = eN_k;

        ++pid;
      }
    }
  }

  // ---------- DEBUG: PASS 2 summary ----------
  std::cout << "[RWPT] PASS2 summary: particles_created(pid)=" << pid
            << " (should equal N_total=" << N_total << ")\n";

  std::fill(u_new.begin(), u_new.end(), 0.0);

// 3b) Deposit each particle's mass m_p to the element nodes with P1 shape values
  long long n_outside = 0;
    deposit_particles_to_nodes_conservative(
    p_X.data(),                // particle positions [pid * nSpace]
    pid,                       // number of alive particles
    mass_per_particle,         // m_p
    mesh_dof,                  // node coordinates (shape: [nNodes, 3], uses first nSpace)
    mesh_l2g,                  // element->node connectivity
    u_l2g,                    // element->u_dof connectivity
    nElements_global,          // number of elements
    nSpace+1,
    nDOF_trial_element,        // must be nSpace+1 for P1 simplex
    u_new,                     // RHS accumulator b_I  (mass per node)
    n_outside                  // counter for particles outside the mesh
);
if (pid > 0) {
  double Mtot_RHS = 0.0;
  for (int I = 0; I < numDOFs; ++I) Mtot_RHS += u_new.data()[I];
  std::cout << "[RWPT] deposit: pid=" << pid
            << " outside=" << n_outside
            << " total_mass_RHS=" << std::setprecision(16) << Mtot_RHS
            << " (should be pid * m_p = " << (static_cast<double>(pid)*mass_per_particle) << ")\n";
}
// Convert node mass to concentration using lumped mass matrix
for (int I = 0; I < numDOFs; ++I) {
  const double MLI = std::max(lumped_mass_matrix.data()[I], 1e-16);
  u_dof.data()[I] = u_new.data()[I] / MLI;
}

// (Optional) Rebuild q_u and qp_mass from u_dof for diagnostics/output
for (int eN = 0; eN < nElements_global; ++eN) {
  const int* l2g_e = &u_l2g.data()[eN*nDOF_trial_element];
  for (int k = 0; k < nQuadraturePoints_element; ++k) {
    const int q = eN*nQuadraturePoints_element + k;
    const double* Nk = &u_trial_ref.data()[k*nDOF_trial_element];

    double val = 0.0;
    for (int j = 0; j < nDOF_trial_element; ++j)
      val += Nk[j] * u_dof.data()[ l2g_e[j] ];

    q_u.data()[q] = val;
    const double phi = q_porosity.data()[q];
    const double dV  = qp_dV.data()[q];
    qp_mass.data()[q] = val * phi * dV;
  }
}

// Final stats
{
  double fmin = 1e300, fmax = -1e300;
  for (int I = 0; I < numDOFs; ++I) {
    const double v = u_dof.data()[I];
    if (v < fmin) fmin = v;
    if (v > fmax) fmax = v;
  }
  const int showN = (numDOFs < 5 ? numDOFs : 5);
  std::cerr << "[RWPT] u_dof(out) min=" << fmin << " max=" << fmax << " head:";
  for (int i = 0; i < showN; ++i) std::cerr << " " << u_dof.data()[i];
  std::cerr << "\n";
}
}


  // // ---------------- PASS 3: map particles to nearest global quadrature point -------------------
  // const int nQP_total = nElements_global * nQuadraturePoints_element;

  // for (int q = 0; q < std::min(5, nQP_total); ++q) {
  //   std::cout << "[RWPT] q_x[" << q << "] = (";
  //   for (int I=0; I<nSpace; ++I) {
  //     if (I) std::cout << ", ";
  //     std::cout << q_x.data()[q*nSpace + I];
  //   }
  //   std::cout << ")\n";
  // }
  // for (int p = 0; p < std::min(5, pid); ++p) {
  //   std::cout << "[RWPT] p_X[" << p << "] = (";
  //   for (int I=0; I<nSpace; ++I) {
  //     if (I) std::cout << ", ";
  //     std::cout << p_X[static_cast<size_t>(p)*nSpace + I];
  //   }
  //   std::cout << ")\n";
  // }

  // #pragma omp parallel for schedule(static)
  // for (int p = 0; p < pid; ++p)
  // {
  //   const double* xpp = &p_X[static_cast<size_t>(p)*nSpace];

  //   double best_d2 = 1.0e300;
  //   int    best_q  = 0;

  //   for (int q = 0; q < nQP_total; ++q)
  //   {
  //     const double* qx = &q_x.data()[q*nSpace];

  //     double d2 = 0.0;
  //     for (int I = 0; I < nSpace; ++I)
  //     {
  //       const double d = xpp[I] - qx[I];
  //       d2 += d*d;
  //     }
  //     if (d2 <= best_d2)
  //     {
  //       best_d2 = d2;
  //       best_q  = q;
  //     }
  //   }
  //   p_qp[static_cast<size_t>(p)] = best_q;
  // }

//   // ---------------- PASS 4: counts → mass → concentration at each qp --------------------------
//   xt::pyarray<int>& qp_counts = args.array<int>("qp_counts"); // pre-allocated: size nQP_total
//   std::fill(qp_counts.begin(), qp_counts.end(), 0);

//   #pragma omp parallel for schedule(static)
//   for (int p = 0; p < pid; ++p)
//   {
//     const int q = p_qp[static_cast<size_t>(p)];
//     if (0 <= q && q < nQP_total)
//     {
//       #pragma omp atomic
//       qp_counts.data()[q] += 1;
//     }
//   }

//   long long sum_counts = 0;
//   int       max_count  = 0;
//   for (int q = 0; q < nQP_total; ++q) {
//     const int c = qp_counts.data()[q];
//     sum_counts += c;
//     if (c > max_count) max_count = c;
//   }
//   std::cout << "[RWPT] PASS4 counts: sum_counts=" << sum_counts
//             << " max_count=" << max_count << "\n";

//   #pragma omp parallel for schedule(static)
//   for (int q = 0; q < nQP_total; ++q)
//   {
//     const double phi = q_porosity.data()[q];
//     const double dV  = qp_dV.data()[q];
//     const double m_q = static_cast<double>(qp_counts.data()[q]) * mass_per_particle;

//     qp_mass.data()[q] = m_q;

//     const double denom = phi * dV;
//     q_u.data()[q] = (denom > 1.0e-16) ? (m_q / denom) : 0.0;
//   }

//   {
//     const int PRINT_Q = 5;
//     for (int q = 0; q < std::min(nQP_total, PRINT_Q); ++q) {
//       std::cout << "[RWPT] PASS4 sample q=" << q
//                 << " count=" << qp_counts.data()[q]
//                 << " m_q=" << qp_mass.data()[q]
//                 << " q_u=" << q_u.data()[q] << "\n";
//     }
//   }

//   // ---------------- scatter back to DOFs (lumped projection with ML) --------------------------
//   std::fill(u_new.begin(), u_new.end(), 0.0);

//   for (int eN = 0; eN < nElements_global; ++eN)
//   {
//     const int* l2g_e = &u_l2g.data()[eN*nDOF_trial_element];

//     for (int k = 0; k < nQuadraturePoints_element; ++k)
//     {
//       const int eN_k     = eN*nQuadraturePoints_element + k;
//       const double w     = std::fabs(qp_dV.data()[eN_k]);            // ensure positive weight
//       const double val   = q_u.data()[eN_k];                         // updated concentration at qp
//       const double* Nk   = &u_trial_ref.data()[k*nDOF_trial_element];

//       for (int j = 0; j < nDOF_trial_element; ++j)
//       {
//         const int    I  = l2g_e[j];
//         const double Nj = Nk[j];
//         u_new.data()[I] += Nj * (val * w);                           // b_i += ∫ N_i * u_h
//       }
//     }
//   }

//   // ---------- DEBUG: u_new (RHS) stats ----------
//   {
//     double bmin = 1e300, bmax = -1e300;
//     for (int I = 0; I < numDOFs; ++I) {
//       const double v = u_new.data()[I];
//       if (v < bmin) bmin = v;
//       if (v > bmax) bmax = v;
//     }
//     const int showN = (numDOFs < 5 ? numDOFs : 5);
//     std::cerr << "[RWPT] u_new(RHS) min=" << bmin << " max=" << bmax << " head:";
//     for (int i = 0; i < showN; ++i) std::cerr << " " << u_new.data()[i];
//     std::cerr << "\n";
//   }

//   // finalize DOFs with precomputed lumped mass (same philosophy as FCTStep)
//   for (int I = 0; I < numDOFs; ++I)
//   {
//     const double mI = (lumped_mass_matrix.data()[I] > 1.0e-16)
//                       ? lumped_mass_matrix.data()[I] : 1.0e-16;
//     u_dof.data()[I] = u_new.data()[I] / mI;
//   }

//   // ---------- DEBUG: u_dof (final) stats ----------
//   {
//     double fmin = 1e300, fmax = -1e300;
//     for (int I = 0; I < numDOFs; ++I) {
//       const double v = u_dof.data()[I];
//       if (v < fmin) fmin = v;
//       if (v > fmax) fmax = v;
//     }
//     const int showN = (numDOFs < 5 ? numDOFs : 5);
//     std::cerr << "[RWPT] u_dof(out) min=" << fmin << " max=" << fmax << " head:";
//     for (int i = 0; i < showN; ++i) std::cerr << " " << u_dof.data()[i];
//     std::cerr << "\n";
//   }
// }



};//TADR


inline TADR_base* newTADR(int nSpaceIn,
                          int nQuadraturePoints_elementIn,
                          int nDOF_mesh_trial_elementIn,
                          int nDOF_trial_elementIn,
                          int nDOF_test_elementIn,
                          int nQuadraturePoints_elementBoundaryIn,
                          int CompKernelFlag)
{
  if (nSpaceIn == 1)
    return proteus::chooseAndAllocateDiscretization1D<TADR_base,TADR,CompKernel>(nSpaceIn,
                                                                                 nQuadraturePoints_elementIn,
                                                                                 nDOF_mesh_trial_elementIn,
                                                                                 nDOF_trial_elementIn,
                                                                                 nDOF_test_elementIn,
                                                                                 nQuadraturePoints_elementBoundaryIn,
                                                                                 CompKernelFlag);
  else if (nSpaceIn == 2)
    return proteus::chooseAndAllocateDiscretization2D<TADR_base,TADR,CompKernel>(nSpaceIn,
                                                                                 nQuadraturePoints_elementIn,
                                                                                 nDOF_mesh_trial_elementIn,
                                                                                 nDOF_trial_elementIn,
                                                                                 nDOF_test_elementIn,
                                                                                 nQuadraturePoints_elementBoundaryIn,
                                                                                 CompKernelFlag);
  else
    return proteus::chooseAndAllocateDiscretization<TADR_base,TADR,CompKernel>(nSpaceIn,
                                                                               nQuadraturePoints_elementIn,
                                                                               nDOF_mesh_trial_elementIn,
                                                                               nDOF_trial_elementIn,
                                                                               nDOF_test_elementIn,
                                                                               nQuadraturePoints_elementBoundaryIn,
                                                                               CompKernelFlag);
}
}//proteus
#endif
