#ifndef CO2_H
#define CO2_H
#include <cmath>
#include <iostream>
#include <valarray>
#include "CompKernel.h"
#include "ModelFactory.h"
#include "../mprans/ArgumentsDict.h"
#include "xtensor-python/pyarray.hpp"
#define nnz nSpace

namespace py = pybind11;
#define POWER_SMOOTHNESS_INDICATOR 2
#define IS_BETAij_ONE              0
#define GLOBAL_FCT                 0
namespace proteus
{
enum class STABILIZATION : int {
  Galerkin         = 0,
  EV_Stab          = 1,
  EntropyViscosity = 2,
  Implicit_FCT     = 3
};
enum class PSK : int {
  VG_PSK = 0,   
  BC_PSK  = 1   
};



namespace mphase_co2
{
//cek todo: revisit entry for mass transport form
// Power entropy //
inline double ENTROPY(const double &phi, const double &phiL, const double &phiR)
{
  return 1. / 2. * std::pow(fabs(phi), 2.);
}
inline double DENTROPY(const double &phi, const double &phiL, const double &phiR)
{
  return fabs(phi) * (phi >= 0 ? 1 : -1);
}
// Log entropy // for level set from 0 to 1
inline double ENTROPY_LOG(const double &phi, const double &phiL, const double &phiR)
{
  return std::log(fabs((phi - phiL) * (phiR - phi)) + 1E-14);
}
inline double DENTROPY_LOG(const double &phi, const double &phiL, const double &phiR)
{
  return (phiL + phiR - 2 * phi) * ((phi - phiL) * (phiR - phi) >= 0 ? 1 : -1) / (fabs((phi - phiL) * (phiR - phi)) + 1E-14);
}
} // namespace mphase_co2
} // namespace proteus
namespace proteus
{
namespace mphase_co2
{
class Mphase_co2_base {
  //The base class defining the interface
public:
  virtual ~Mphase_co2_base() { double anb_seepage_flux = 1e-16; }
  virtual void calculateResidual(arguments_dict &args)                   = 0;
  virtual void calculateJacobian(arguments_dict &args)                   = 0;
  virtual void invert(arguments_dict &args)                              = 0;
  virtual void FCTStep(arguments_dict &args)                             = 0;
  virtual void kth_FCT_step(arguments_dict &args)                        = 0;
  virtual void calculateResidual_entropy_viscosity(arguments_dict &args) = 0;
  virtual void calculateMassMatrix(arguments_dict &args)                 = 0;
};

template <class CompKernelType, int nSpace, int nQuadraturePoints_element, int nDOF_mesh_trial_element, int nDOF_trial_element, int nDOF_test_element, int nQuadraturePoints_elementBoundary>
class mphase_co2 : public Mphase_co2_base {
public:
  const int      nDOF_test_X_trial_element;
  CompKernelType ck;
  mphase_co2() : nDOF_test_X_trial_element(nDOF_test_element * nDOF_trial_element), ck() { }

inline void evaluateCoefficients(const int rowptr[nSpace], const int colind[nnz],
                                 const int k, // phase flag (k=0 water, k=1 air)
                                 const double rho_water, const double rho_air,        // phase density
                                 const double beta_water, const double beta_air,      // phase beta
                                 const double gravity[nSpace],
                                 const double alpha, const double n_vg,
                                 const double thetaR, const double thetaSR,
                                 const double KWs[nnz],
                                 const double &u_water,  const double &u_air,         // both heads
                                 const double Y_water,   const double Y_air,          // mass fractions
                                 double &m, double &dm,
                                 double f[nSpace], double df[nSpace],
                                 double a[nnz], double da[nnz], double as[nnz],
                                 double &kr, double &dkr,
                                 const PSK PSK_TYPE,
                                 double &Swater_out , double &Sair_out,               // Sw,Sg outputs
                                 const double BC_entry_head, const double BC_lambda)
{
  const int nSpace2 = nSpace * nSpace;
  double    psiC;
  double    pcBar;
  double    pcBar_n;
  double    pcBar_nM1;
  double    pcBar_nM2;
  double    onePlus_pcBar_n;
  double    sBar;
  double    sqrt_sBar;
  double    DsBar_DpsiC;
  double    thetaW;
  double    DthetaW_DpsiC;
  double    vBar;
  double    vBar2;
  double    DvBar_DpsiC;
  double    KWr;
  double    DKWr_DpsiC;
  // double    rho2 = rho * rho; // no single-phase rho 
  double    thetaS;
  double    rhom;
  double    drhom;
  double    m_vg;
  double    pcBarStar;
  double    sqrt_sBarStar;

  // -------------------- NEW: select phase, build psic --------------------
  const double rho_k  = (k==0) ? rho_water  : rho_air;     // phase density
  const double beta_k = (k==0) ? beta_water : beta_air;    // phase compressibility
  const double u_k    = (k==0) ? u_water    : u_air;       // phase head
  const double Y_k    = (k==0) ? Y_water    : Y_air;       // phase mass fraction

  thetaS = thetaR + thetaSR;                               // (as before)
  // Capillary head psiC = (rho_a/rho_w) * u_a − u_w
  psiC = (rho_air / rho_water) * u_air - u_water;          // NEW

  // Chain rule factors dpsiC/du
  const double dpsic_duw = -1.0;                           // NEW
  const double dpsic_dua = (rho_air / rho_water);          // NEW
  const double dpsic_duk = (k==0) ? dpsic_duw : dpsic_dua; // NEW

  // Use selected phase density for ρ^2 terms
  const double rho2 = rho_k * rho_k;                       // NEW

  // -------------------- PSK branches (VG/BC) to get Se and dSe/dψc --------------------
  if (PSK_TYPE == PSK::VG_PSK) {
    psiC   = psiC;                     // (keep symbol; already set above)
    m_vg   = 1.0 - 1.0 / n_vg;
    if (psiC > 0.0) {
      pcBar     = alpha * psiC;
      pcBarStar = (pcBar < 1.0e-8) ? 1.0e-8 : pcBar;
      pcBar_nM2       = pow(pcBarStar, n_vg - 2.0);
      pcBar_nM1       = pcBar_nM2 * pcBar;
      pcBar_n         = pcBar_nM1 * pcBar;
      onePlus_pcBar_n = 1.0 + pcBar_n;
      sBar = pow(onePlus_pcBar_n, -m_vg);
      DsBar_DpsiC = alpha * (1.0 - n_vg) * (sBar / onePlus_pcBar_n) * pcBar_nM1;
      vBar        = 1.0 - pcBar_nM1 * sBar;
      vBar2       = vBar * vBar;
      DvBar_DpsiC = -alpha * (n_vg - 1.0) * pcBar_nM2 * sBar - pcBar_nM1 * DsBar_DpsiC;

      thetaW        = thetaSR * sBar + thetaR;
      DthetaW_DpsiC = thetaSR * DsBar_DpsiC;

      // Wetting rel-perm krw, dkrw/dψ (keep your symbols KWr and DKWr_DpsiC)
      sqrt_sBar     = sqrt(sBar);
      sqrt_sBarStar = (sqrt_sBar < 1.0e-8) ? 1.0e-8 : sqrt_sBar;
      KWr           = sqrt_sBar * vBar2;
      DKWr_DpsiC    = (0.5 / sqrt_sBarStar) * DsBar_DpsiC * vBar2
                    + 2.0 * sqrt_sBar * vBar * DvBar_DpsiC;
    } else {
      thetaW        = thetaS;
      DthetaW_DpsiC = 0.0;
      KWr           = 1.0;
      DKWr_DpsiC    = 0.0;
      sBar          = 1.0;
      DsBar_DpsiC   = 0.0;
    }
  } else {
    // -------------------- Brooks–Corey (Mualem) branch (kept) --------------------
    //psiC   = psiC; // already set above
    const double hb  = (BC_entry_head > 1e-12) ? BC_entry_head : 1e-12;
    const double lam = (BC_lambda     > 1e-12) ? BC_lambda     : 1e-12;

    if (psiC > 0.0) {
      if (psiC >= hb) {
        sBar        = pow(hb / psiC, lam);       // Se
        DsBar_DpsiC = -lam * sBar / psiC;        // dSe/dψ

        thetaW        = thetaR + thetaSR * sBar; // θw
        DthetaW_DpsiC = thetaSR * DsBar_DpsiC;   // dθw/dψ

        // Wetting rel-perm (Mualem–BC): krw = Se^{3 + 2/λ}
        const double expo = 3.0 + 2.0 / lam;
        KWr           = pow(sBar, expo);
        DKWr_DpsiC    = (expo * pow(sBar, expo - 1.0)) * DsBar_DpsiC;
      } else {
        // Saturated plateau
        sBar          = 1.0;
        DsBar_DpsiC   = 0.0;
        thetaW        = thetaS;
        DthetaW_DpsiC = 0.0;
        KWr           = 1.0;
        DKWr_DpsiC    = 0.0;
      }
    } else {
      sBar          = 1.0;
      DsBar_DpsiC   = 0.0;
      thetaW        = thetaS;
      DthetaW_DpsiC = 0.0;
      KWr           = 1.0;
      DKWr_DpsiC    = 0.0;
    }
  }

  // -------------------- NEW: pick phase content and kra/krw --------------------
  // Saturations (Sw + Sg = 1)
  const double Sw = thetaW / std::max(thetaS, 1.0e-12);    // NEW
  const double Sg = 1.0 - Sw;                              // NEW
  Swater_out = Sw;                                         // NEW
  Sair_out   = Sg;                                         // NEW

  // Phase content theta_k and derivative w.r.t psiC
  double theta_phase, DthetaPhase_DpsiC;                   // NEW
  if (k==0) { // water
    theta_phase        = thetaW;
    DthetaPhase_DpsiC  = DthetaW_DpsiC;
  } else {    // air: theta_a = theta_s − theta_w
    theta_phase        = thetaS - thetaW;
    DthetaPhase_DpsiC  = -DthetaW_DpsiC;
  }

  // Select relative permeability for phase k
  double kr_phase, DkrPhase_DpsiC;                         // NEW
  if (k==0) {
    // Use the wetting values already computed
    kr_phase       = KWr;
    DkrPhase_DpsiC = DKWr_DpsiC;
  } else {
    // Non-wetting rel-perm (kra) and its derivative
    const double Se_cl = std::fmin(std::fmax(sBar, 1.0e-12), 1.0 - 1.0e-12);
    if (PSK_TYPE == PSK::VG_PSK) {
      const double mvg     = 1.0 - 1.0/n_vg;
      const double Se_pow  = std::pow(Se_cl, 1.0/mvg);
      const double t2      = 1.0 - Se_pow;

      kr_phase = std::sqrt(1.0 - Se_cl) * std::pow(t2, 2.0*mvg);

      const double dt1_dSe = -0.5 / std::fmax(std::sqrt(1.0 - Se_cl), 1.0e-8);
      const double dt2_dSe = -(1.0/mvg) * std::pow(Se_cl, 1.0/mvg - 1.0);
      const double dkr_dSe = dt1_dSe * std::pow(t2, 2.0*mvg)
                           + std::sqrt(1.0 - Se_cl) * (2.0*mvg) * std::pow(t2, 2.0*mvg - 1.0) * dt2_dSe;
      DkrPhase_DpsiC = dkr_dSe * DsBar_DpsiC;
    } else {
      const double lam     = (BC_lambda > 1.0e-12) ? BC_lambda : 1.0e-12;
      const double expo_a  = 2.0 + 1.0/lam;
      kr_phase             = std::pow(1.0 - Se_cl, expo_a);
      const double dkr_dSe = -expo_a * std::pow(1.0 - Se_cl, expo_a - 1.0);
      DkrPhase_DpsiC       = dkr_dSe * DsBar_DpsiC;
    }
  }
  rhom  = rho_k * std::exp(beta_k * u_k);                   // CHANGED
  drhom = beta_k * rhom;                                    // CHANGED

  const double dtheta_duk = DthetaPhase_DpsiC * dpsic_duk;  // NEW
  const double dkr_duk    = DkrPhase_DpsiC    * dpsic_duk;  // NEW

  // Mass term with optional mass fraction Y_k (set Y_k=1.0 for pure phase mass)
  m  = rhom * theta_phase * Y_k;                            // CHANGED
  dm = rhom * (dtheta_duk * Y_k) + drhom * theta_phase * Y_k;// CHANGED

  // -------------------- CHANGED: assemble using selected-phase kr and derivative wrt u_k --------------------
  for (int I = 0; I < nSpace; I++) {
    f[I]  = 0.0;
    df[I] = 0.0;
    for (int ii = rowptr[I]; ii < rowptr[I + 1]; ii++) {
      const int J = colind[ii];

      // Gravity/advection-like term
      f[I]  += rho2 * kr_phase * KWs[ii] * gravity[J];     // rho2 is selected-phase
      df[I] += rho2 * dkr_duk  * KWs[ii] * gravity[J];     // derivative wrt selected u_k

      // Diffusion/mobility-like tensors
      a[ii]  = rho_k * kr_phase * KWs[ii];                  // 
      da[ii] = rho_k * dkr_duk  * KWs[ii];                  // 

      as[ii] = rho_k * KWs[ii];                             // 
    }
  }

  // Return selected-phase kr and its derivative wrt selected variable
  kr  = kr_phase;                                           // NEW
  dkr = dkr_duk;                                            // used to be -DKWr_DpsiC
}


// Two-phase inverse: recover the selected phase head (k=0 water, k=1 air)
// from the stored mass m and the other phase head, using PSK inverse.
//
// Notes:
// - Uses m ≈ ρ_k * θ_k * Y_k  (same as forward) and ignores slight
//   compressibility in the inversion (like your original did).
// - Requires both heads by reference; updates only the selected one.
// - If you still need the original single-phase inverse, keep it under a
//   different name (e.g., evaluateInverseCoefficients_1ph).

inline void evaluateInverseCoefficients_2ph(const int rowptr[nSpace], const int colind[nnz],
                                            const int k, // NEW: phase flag (0=water, 1=air)
                                            const double rho_water, const double rho_air,      // NEW
                                            const double beta_water, const double beta_air,    // (unused here, parity with forward)
                                            const double gravity[nSpace],
                                            const double alpha, const double n_vg,
                                            const double thetaR, const double thetaSR,
                                            const double KWs[nnz],
                                            double &u_water, double &u_air,                    // NEW: both heads (in/out)
                                            const double &m,                                   // mass of selected phase (or component)
                                            const double Y_water, const double Y_air,          // NEW: mass fractions
                                            const PSK PSK_TYPE,                                // VG or BC
                                            const double BC_entry_head, const double BC_lambda)
{
    // --- constants & guards ---
  const double eps     = 1.0e-12;
  const double thetaS  = thetaR + thetaSR;

  // --- select phase inputs ---
  const double rho_k = (k==0) ? rho_water : rho_air;     // NEW
  const double Y_k   = (k==0) ? Y_water   : Y_air;       // NEW

  // --- invert m ≈ ρ_k * θ_k * Y_k  →  θ_k  (ignore compressibility in inverse, like original) ---
  const double denom = std::fmax(rho_k * std::fmax(Y_k, eps), eps);          // NEW
  double theta_k     = m / denom;                                            // NEW

  // --- map to θ_w and S_e ---
  double thetaW = (k==0) ? theta_k : (thetaS - theta_k);                     // NEW
  // original gating (keep your guard band)
  if (thetaW > 1.01*thetaR && thetaW < thetaS)                               // unchanged threshold style
  {
    // S_e = (θ_w − θ_R)/θ_SR  in (0,1)
    const double Se = std::fmin(std::fmax((thetaW - thetaR)/std::fmax(thetaSR, eps), eps),
                                1.0 - eps);                                   // NEW clamp

    // -------- inverse PSK to get ψ_c --------
    double psic = 0.0;
    if (PSK_TYPE == PSK::VG_PSK) {
      // Van Genuchten inverse:  ψ_c = ( (Se^{-1/m}-1)^{1/n} ) / α
      const double m_vg  = 1.0 - 1.0/n_vg;
      const double pc_n  = std::pow(Se, -1.0/m_vg) - 1.0;
      const double pc    = std::pow(std::fmax(pc_n, 0.0), 1.0/n_vg); // guard tiny negatives
      psic = pc / std::fmax(alpha, eps);
    } else {
      // Brooks–Corey inverse:  ψ_c = h_b * Se^{-1/λ}
      const double hb  = (BC_entry_head > eps) ? BC_entry_head : eps;
      const double lam = (BC_lambda     > eps) ? BC_lambda     : eps;
      psic = hb * std::pow(Se, -1.0/lam);
    }

    // -------- recover selected head from ψ_c and the other head --------
    // ψ_c = (ρ_a/ρ_w) u_a − u_w
    if (k==0) {
      // water head
      u_water = (rho_air/rho_water) * u_air - psic;                           // NEW
    } else {
      // air head
      u_air   = (psic + u_water) * (rho_water/rho_air);                       // NEW
    }
  }
  // else: outside physical inversion band ⇒ leave u_water/u_air as passed in.
}






inline void evaluateInverseCoefficients(const int rowptr[nSpace], const int colind[nnz],
                                        const double rho, const double beta,
                                        const double gravity[nSpace],
                                        const double alpha, const double n_vg,
                                        const double thetaR, const double thetaSR,
                                        const double KWs[nnz],
                                        double &u,                
                                        const double &m, const double &dm,
                                        const double f[nSpace], const double df[nSpace],
                                        const double a[nnz], const double da[nnz],
                                        const PSK PSK_TYPE,       // NEW: selector
                                        const double BC_entry_head,  // NEW: h_b (>0)
                                        const double BC_lambda)      // NEW: lambda (>0)
{
  double psiC, pcBar, pcBar_n, sBar, thetaW, thetaS, m_vg;
  m_vg   = 1.0 - 1.0 / n_vg;
  thetaS = thetaR + thetaSR;
  thetaW = m / rho;

  // same gating as your original
  if (thetaW > 1.01*thetaR && thetaW < thetaS) {
    sBar = (thetaW - thetaR) / thetaSR;          // S_e

    if (PSK_TYPE == PSK::VG_PSK) {
      // ------- original VG inverse (unchanged) -------
      pcBar_n = pow(sBar, -1.0 / m_vg) - 1.0;
      pcBar   = pow(pcBar_n, 1.0 / n_vg);
      psiC    = pcBar / alpha;
      u       = -psiC;
    } else {
      // ------- Brooks–Corey inverse -------
      const double hb  = (BC_entry_head > 1e-12) ? BC_entry_head : 1e-12;
      const double lam = (BC_lambda     > 1e-12) ? BC_lambda     : 1.0e-12;
      // psi_c = h_b * S_e^{-1/λ}   (for 0<S_e<1)
      const double Se_clamped = std:fmax(std::fmin(sBar, 1.0-1e-12), 1e-12);
      psiC = hb * std::pow(Se_clamped, -1.0/lam);
      u    = -psiC;
    }
  }
}

  inline void calculateCFL(const double &elementDiameter, const double df[nSpace], double &cfl)
  {
    double h, nrm_v;
    h     = elementDiameter;
    nrm_v = 0.0;
    for (int I = 0; I < nSpace; I++) nrm_v += df[I] * df[I];
    nrm_v = sqrt(nrm_v);
    cfl   = nrm_v / h;
  }

  inline void calculateSubgridError_tau(const double &elementDiameter, const double &dmt, const double dH[nSpace], double &cfl, double &tau)
  {
    double h, nrm_v, oneByAbsdt;
    h     = elementDiameter;
    nrm_v = 0.0;
    for (int I = 0; I < nSpace; I++) nrm_v += dH[I] * dH[I];
    nrm_v      = sqrt(nrm_v);
    cfl        = nrm_v / h;
    oneByAbsdt = fabs(dmt);
    tau        = 1.0 / (2.0 * nrm_v / h + oneByAbsdt + 1.0e-8);
  }

  inline void calculateSubgridError_tau(const double &Ct_sge, const double G[nSpace * nSpace], const double &A0, const double Ai[nSpace], double &tau_v, double &cfl)
  {
    double v_d_Gv = 0.0;
    for (int I = 0; I < nSpace; I++)
      for (int J = 0; J < nSpace; J++) v_d_Gv += Ai[I] * G[I * nSpace + J] * Ai[J];
    tau_v = 1.0 / sqrt(Ct_sge * A0 * A0 + v_d_Gv);
  }

  inline void calculateNumericalDiffusion(const double &shockCapturingDiffusion, const double &elementDiameter, const double &strong_residual, const double grad_u[nSpace], double &numDiff)
  {
    double h, num, den, n_grad_u;
    h        = elementDiameter;
    n_grad_u = 0.0;
    for (int I = 0; I < nSpace; I++) n_grad_u += grad_u[I] * grad_u[I];
    num     = shockCapturingDiffusion * 0.5 * h * fabs(strong_residual);
    den     = sqrt(n_grad_u) + 1.0e-8;
    numDiff = num / den;
  }

  inline void exteriorNumericalFlux(const double &bc_flux, int rowptr[nSpace], int colind[nnz], int isSeepageFace, int &isDOFBoundary, double n[nSpace], double bc_u, double K[nnz], double grad_psi[nSpace], double u, double K_rho_g[nSpace], double penalty, double &flux)
  {
    double v_I, bc_u_seepage = 0.0;
    if (isSeepageFace || isDOFBoundary) {
      flux = 0.0;
      for (int I = 0; I < nSpace; I++) {
        //gravity
        v_I = K_rho_g[I];
        //pressure head
        for (int m = rowptr[I]; m < rowptr[I + 1]; m++) { v_I -= K[m] * grad_psi[colind[m]]; }
        flux += v_I * n[I];
      }
      if (isSeepageFace) bc_u = bc_u_seepage;
      flux += penalty * (u - bc_u);
      //flux -= penalty * bc_u;
      if (isSeepageFace) {
        if (flux > 0.0) {
          isDOFBoundary = 1;
          bc_u          = bc_u_seepage;
        } else {
          isDOFBoundary = 0;
          flux          = 0.0;
        }
      }
    } else flux = bc_flux;
  }

  void exteriorNumericalFluxJacobian(const int rowptr[nSpace], const int colind[nnz], const int isDOFBoundary, const double n[nSpace], const double K[nnz], const double dK[nnz], const double grad_psi[nSpace], const double grad_v[nSpace], const double dK_rho_g[nSpace], const double v, const double penalty, double &fluxJacobian)
  {
    if (isDOFBoundary) {
      fluxJacobian = 0.0;
      for (int I = 0; I < nSpace; I++) {
        //gravity
        fluxJacobian += dK_rho_g[I] * v * n[I];
        //pressure head
        for (int m = rowptr[I]; m < rowptr[I + 1]; m++) { fluxJacobian -= (K[m] * grad_v[colind[m]] + dK[m] * v * grad_psi[colind[m]]) * n[I]; }
      }
      //Dirichlet penalty
      fluxJacobian += penalty * v;
    } else fluxJacobian = 0.0;
  }

inline void exteriorNumericalFlux2(const double &bc_flux, int rowptr[nSpace], int colind[nnz], int isSeepageFace, int &isDOFBoundary, double n[nSpace], double bc_u, double K[nnz], double grad_psi[nSpace], double u, double K_rho_g[nSpace], double penalty, double &flux, double &bflux)
  {
    double v_I, bc_u_seepage = 0.0;
    if (isSeepageFace || isDOFBoundary) {
      flux = 0.0;
      bflux = 0.0;
      for (int I = 0; I < nSpace; I++) {
        //gravity
        v_I = K_rho_g[I];
        //pressure head
        for (int m = rowptr[I]; m < rowptr[I + 1]; m++) { v_I -= K[m] * grad_psi[colind[m]]; }
        flux += v_I * n[I];
      }
      if (isSeepageFace) bc_u = bc_u_seepage;
      flux += penalty * (u - bc_u);
      bflux += penalty * (u - bc_u);
      if (isSeepageFace) {
        if (flux > 0.0) {
          isDOFBoundary = 1;
        } else {
          isDOFBoundary = 0;
          flux          = 0.0;
          bflux         = 0.0;
        }
      }
    } else {
      flux = bc_flux;
      bflux = bc_flux;
    }
  }

  void exteriorNumericalFluxJacobian2(const int rowptr[nSpace], const int colind[nnz], const int isDOFBoundary, const double n[nSpace],  const double Ks[nnz], const double K[nnz], const double dK[nnz], const double grad_psi[nSpace], const double grad_v[nSpace], const double dK_rho_g[nSpace], const double v, const double penalty, double &fluxJacobian, double &bfluxJacobian)
  {
    if (isDOFBoundary) {
      fluxJacobian = 0.0;
      bfluxJacobian = 0.0;
      for (int I = 0; I < nSpace; I++) {
        for (int m = rowptr[I]; m < rowptr[I + 1]; m++) { 
          fluxJacobian -= Ks[m] * grad_v[colind[m]] * n[I]; 
        }
      }
      //Dirichlet penalty
      bfluxJacobian = penalty * v;
    } else {
      fluxJacobian = 0.0;
      bfluxJacobian = 0.0;
    }
  }

  double seepagefluxcalculator(double anb_seepage_flux, int isSeepageFace, double dS, double flux_ext)
  {
    if (isSeepageFace) { anb_seepage_flux += flux_ext * dS; }
    return anb_seepage_flux;
  }
// First make a loop in the calculateResidual function for two phase air and water 
  void calculateResidual(arguments_dict &args)
  {
    xt::pyarray<double> &mesh_trial_ref                             = args.array<double>("mesh_trial_ref");
    xt::pyarray<double> &mesh_grad_trial_ref                        = args.array<double>("mesh_grad_trial_ref");
    xt::pyarray<double> &mesh_dof                                   = args.array<double>("mesh_dof");
    xt::pyarray<double> &mesh_velocity_dof                          = args.array<double>("mesh_velocity_dof");
    double               MOVING_DOMAIN                              = args.scalar<double>("MOVING_DOMAIN");
    xt::pyarray<int>    &mesh_l2g                                   = args.array<int>("mesh_l2g");
    xt::pyarray<double> &dV_ref                                     = args.array<double>("dV_ref");
    xt::pyarray<double> &u_trial_ref                                = args.array<double>("u_trial_ref");
    xt::pyarray<double> &u_grad_trial_ref                           = args.array<double>("u_grad_trial_ref");
    xt::pyarray<double> &u_test_ref                                 = args.array<double>("u_test_ref");
    xt::pyarray<double> &u_grad_test_ref                            = args.array<double>("u_grad_test_ref");
    xt::pyarray<double> &mesh_trial_trace_ref                       = args.array<double>("mesh_trial_trace_ref");
    xt::pyarray<double> &mesh_grad_trial_trace_ref                  = args.array<double>("mesh_grad_trial_trace_ref");
    xt::pyarray<double> &dS_ref                                     = args.array<double>("dS_ref");
    xt::pyarray<double> &u_trial_trace_ref                          = args.array<double>("u_trial_trace_ref");
    xt::pyarray<double> &u_grad_trial_trace_ref                     = args.array<double>("u_grad_trial_trace_ref");
    xt::pyarray<double> &u_test_trace_ref                           = args.array<double>("u_test_trace_ref");
    xt::pyarray<double> &u_grad_test_trace_ref                      = args.array<double>("u_grad_test_trace_ref");
    xt::pyarray<double> &normal_ref                                 = args.array<double>("normal_ref");
    xt::pyarray<double> &boundaryJac_ref                            = args.array<double>("boundaryJac_ref");
    int                  nElements_global                           = args.scalar<int>("nElements_global");
    xt::pyarray<double> &ebqe_penalty_ext                           = args.array<double>("ebqe_penalty_ext");
    xt::pyarray<int>    &elementMaterialTypes                       = args.array<int>("elementMaterialTypes");
    xt::pyarray<int>    &isSeepageFace                              = args.array<int>("isSeepageFace");
    xt::pyarray<int>    &a_rowptr                                   = args.array<int>("a_rowptr");
    xt::pyarray<int>    &a_colind                                   = args.array<int>("a_colind");
    xt::pyarray<double> &gravity                                    = args.array<double>("gravity");
    
    xt::pyarray<double> &thetaR                                     = args.array<double>("thetaR");
    xt::pyarray<double> &thetaSR                                    = args.array<double>("thetaSR");
    xt::pyarray<double> &KWs                                        = args.array<double>("KWs");
    double               useMetrics                                 = args.scalar<double>("useMetrics");
    double               alphaBDF                                   = args.scalar<double>("alphaBDF");
    int                  lag_shockCapturing                         = args.scalar<int>("lag_shockCapturing");
    double               shockCapturingDiffusion                    = args.scalar<double>("shockCapturingDiffusion");
    double               sc_uref                                    = args.scalar<double>("sc_uref");
    double               sc_alpha                                   = args.scalar<double>("sc_alpha");
    xt::pyarray<double> &cfl                                        = args.array<double>("cfl");
    xt::pyarray<double> &elementDiameter                            = args.array<double>("elementDiameter");

        // PARAMETERS FOR EDGE VISCOSITY
    int               numDOFs                       = args.scalar<int>("numDOFs");
    int               NNZ                           = args.scalar<int>("NNZ");
    xt::pyarray<int> &csrRowIndeces_DofLoops        = args.array<int>("csrRowIndeces_DofLoops");
    xt::pyarray<int> &csrColumnOffsets_DofLoops     = args.array<int>("csrColumnOffsets_DofLoops");
    xt::pyarray<int> &csrRowIndeces_CellLoops       = args.array<int>("csrRowIndeces_CellLoops");
    xt::pyarray<int> &csrColumnOffsets_CellLoops    = args.array<int>("csrColumnOffsets_CellLoops");
    xt::pyarray<int> &csrColumnOffsets_eb_CellLoops = args.array<int>("csrColumnOffsets_eb_CellLoops");
    // C matrices
    xt::pyarray<double> &Cx         = args.array<double>("Cx");
    xt::pyarray<double> &Cy         = args.array<double>("Cy");
    xt::pyarray<double> &Cz         = args.array<double>("Cz");
    xt::pyarray<double> &CTx        = args.array<double>("CTx");
    xt::pyarray<double> &CTy        = args.array<double>("CTy");
    xt::pyarray<double> &CTz        = args.array<double>("CTz");
    xt::pyarray<double> &ML         = args.array<double>("ML");
    xt::pyarray<double> &delta_x_ij = args.array<double>("delta_x_ij");
        // VMS
    double VMS = args.scalar<double>("VMS");
    // PARAMETERS FOR EDGE BASED STABILIZATION
    double cE = args.scalar<double>("cE");
    double cK = args.scalar<double>("cK");
    // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
    int LUMPED_MASS_MATRIX = args.scalar<int>("LUMPED_MASS_MATRIX");
    STABILIZATION STABILIZATION_TYPE{static_cast<STABILIZATION>(args.scalar<int>("STABILIZATION_TYPE"))};
    //////////////////////////////////////For Brooks- Corey/////////////////
    PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
    double BC_entry_head = args.scalar<double>("BC_entry_head");
    double BC_lambda     = args.scalar<double>("BC_lambda");
    int                  nExteriorElementBoundaries_global          = args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int>    &exteriorElementBoundariesArray             = args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int>    &elementBoundaryElementsArray               = args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int>    &elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
    double               epsFact                                    = args.scalar<double>("epsFact");    
    ///////////////////////////////For both phases/////////////////////////////////////////////////////////
    double               rho_water                                  = args.scalar<double>("rho_water");
    double               rho_air                                    = args.scalar<double>("rho_air");  
    double               beta_water                                 = args.scalar<double>("beta_water");
    double               beta_air                                   = args.scalar<double>("beta_air");   
        // mass fraction
    double Y_air                                                    = args.scalar<double>("Y_air");     // e.g. 1.0
    double Y_water                                                  = args.scalar<double>("Y_water");   // e.g. 1.0
    xt::pyarray<int>    &u_l2g                                      = args.array<int>("u_l2g");
    //xt::pyarray<int>    &u_l2g_air                                  = args.array<int>("u_l2g_air");   
    xt::pyarray<double> &u_dof_water                                = args.array<double>("u_dof_water");
    xt::pyarray<double> &u_dof_air                                  = args.array<double>("u_dof_air");
    xt::pyarray<double> &u_dof_old_water                            = args.array<double>("u_dof_old_water");
    xt::pyarray<double> &u_dof_old_air                              = args.array<double>("u_dof_old_air");
    xt::pyarray<double> &velocity_water                             = args.array<double>("velocity_water");
    xt::pyarray<double> &velocity_air                               = args.array<double>("velocity_air");
    xt::pyarray<double> &q_m_water                                  = args.array<double>("q_m_water");
    xt::pyarray<double> &q_m_air                                    = args.array<double>("q_m_air");
    xt::pyarray<double> &q_u_water                                  = args.array<double>("q_u_water");
    xt::pyarray<double> &q_u_air                                    = args.array<double>("q_u_air");
    xt::pyarray<double> &q_dV_water                                 = args.array<double>("q_dV_water");
    xt::pyarray<double> &q_dV_air                                   = args.array<double>("q_dV_air");
    xt::pyarray<double> &q_m_betaBDF_water                          = args.array<double>("q_m_betaBDF_water");
    xt::pyarray<double> &q_m_betaBDF_air                            = args.array<double>("q_m_betaBDF_air");
    xt::pyarray<double> &q_numDiff_u_water                          = args.array<double>("q_numDiff_u_water");    
    xt::pyarray<double> &q_numDiff_u_air                            = args.array<double>("q_numDiff_u_air");
    xt::pyarray<double> &q_numDiff_u_last_water                     = args.array<double>("q_numDiff_u_last_water");
    xt::pyarray<double> &q_numDiff_u_last_air                       = args.array<double>("q_numDiff_u_last_air");
    int                  offset_u_water                             = args.scalar<int>("offset_u_water");
    int                  offset_u_air                               = args.scalar<int>("offset_u_air");
    int                  stride_u_water                             = args.scalar<int>("stride_u_water");
    int                  stride_u_air                               = args.scalar<int>("stride_u_air");
    xt::pyarray<double> &globalResidual_water                       = args.array<double>("globalResidual_water");
    xt::pyarray<double> &globalResidual_air                         = args.array<double>("globalResidual_air");
    xt::pyarray<double> &ebqe_velocity_ext_water                    = args.array<double>("ebqe_velocity_ext_water");
    xt::pyarray<double> &ebqe_velocity_ext_air                      = args.array<double>("ebqe_velocity_ext_air");
    xt::pyarray<int>    &isDOFBoundary_u_water                      = args.array<int>("isDOFBoundary_u_water");
    xt::pyarray<int>    &isDOFBoundary_u_air                        = args.array<int>("isDOFBoundary_u_air");
    xt::pyarray<double> &ebqe_bc_u_ext_water                        = args.array<double>("ebqe_bc_u_ext_water");
    xt::pyarray<double> &ebqe_bc_u_ext_air                          = args.array<double>("ebqe_bc_u_ext_air");
    xt::pyarray<int>    &isFluxBoundary_u_water                     = args.array<int>("isFluxBoundary_u_water");
    xt::pyarray<int>    &isFluxBoundary_u_air                       = args.array<int>("isFluxBoundary_u_air");
    xt::pyarray<double> &ebqe_bc_flux_ext_water                     = args.array<double>("ebqe_bc_flux_ext_water");
    xt::pyarray<double> &ebqe_bc_flux_ext_air                       = args.array<double>("ebqe_bc_flux_ext_air");
    xt::pyarray<double> &ebqe_phi_water                             = args.array<double>("ebqe_phi_water");
    xt::pyarray<double> &ebqe_phi_air                               = args.array<double>("ebqe_phi_air");
    xt::pyarray<double> &ebqe_u_water                               = args.array<double>("ebqe_u_water");
    xt::pyarray<double> &ebqe_u_air                                 = args.array<double>("ebqe_u_air");
    xt::pyarray<double> &ebqe_flux_water                            = args.array<double>("ebqe_flux_water");
    xt::pyarray<double> &ebqe_flux_air                              = args.array<double>("ebqe_flux_air");


    // PARAMETERS FOR LOG BASED ENTROPY FUNCTION
    double uL_water = args.scalar<double>("uL_water");
    double uL_air = args.scalar<double>("uL_air");
    double uR_water = args.scalar<double>("uR_water");
    double uR_air = args.scalar<double>("uR_air");
    


 
    

    ///////////////////////////////////////////////////////////////////
    int ENTROPY_TYPE = args.scalar<int>("ENTROPY_TYPE");
    // FOR FCT
    xt::pyarray<double> &dLow                 = args.array<double>("dLow");
    xt::pyarray<double> &fluxMatrix           = args.array<double>("fluxMatrix");
    // AUX QUANTITIES OF INTEREST
    xt::pyarray<double> &quantDOFs = args.array<double>("quantDOFs");

    assert(a_rowptr.data()[nSpace] == nnz);
    assert(a_rowptr.data()[nSpace] == nSpace);
    //cek should this be read in?
    double Ct_sge = 4.0;

    xt::pyarray<double> &anb_seepage_flux_n = args.array<double>("anb_seepage_flux_n");

    //double anb_seepage_flux=0.0;
    double &anb_seepage_flux(args.scalar<double>("anb_seepage_flux"));
    anb_seepage_flux = 0.0;

    //loop over elements to compute volume integrals and load them into element and global residual
    //
    //eN is the element index
    //eN_k is the quadrature point index for a scalar
    //eN_k_nSpace is the quadrature point index for a vector
    //eN_i is the element test function index
    //eN_j is the element trial function index
    //eN_k_j is the quadrature point index for a trial function
    //eN_k_i is the quadrature point index for a trial function
    for (int eN = 0; eN < nElements_global; eN++) {
      //declare local storage for element residual and initialize
      double elementResidual_u_water[nDOF_test_element];
      double elementResidual_u_air[nDOF_test_element];
      for (int i = 0; i < nDOF_test_element; i++) { elementResidual_u_water[i] = 0.0; } //i
      for (int i = 0; i < nDOF_test_element; i++) { elementResidual_u_air[i] = 0.0; } //i
      //loop over quadrature points and compute integrands
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        //compute indeces and declare local storage
        int eN_k = eN * nQuadraturePoints_element + k, eN_k_nSpace = eN_k * nSpace, eN_nDOF_trial_element = eN * nDOF_trial_element;
        double dV, x,y,z, xt,yt,zt;
        double jac[nSpace*nSpace], jacDet, jacInv[nSpace*nSpace];
        double G[nSpace*nSpace], G_dd_G, tr_G;
        double u_grad_trial[nDOF_trial_element*nSpace];
        double u_test_dV[nDOF_trial_element];
        double u_grad_test_dV[nDOF_test_element*nSpace];

        // --- phase water ---
        double u_water=0.0, grad_u_water[nSpace];
        double m_water=0.0, dm_water=0.0, m_t_water=0.0, dm_t_water=0.0;
        double f_water[nSpace], df_water[nSpace];
        double a_water[nnz],   da_water[nnz],   as_water[nnz];
        double Kr_water=0.0, dKr_water=0.0;
        double pdeRes_water=0.0, Lstar_water[nDOF_test_element], subgrid_water=0.0;
        double tau_water=0.0, tau0_water=0.0, tau1_water=0.0, numDiff0_water=0.0, numDiff1_water=0.0;

        // fill water, assemble to elementResidual_u_water[...]

        // --- phase air ---
        double u_air=0.0, grad_u_air[nSpace];
        double m_air=0.0, dm_air=0.0, m_t_air=0.0, dm_t_air=0.0;
        double f_air[nSpace], df_air[nSpace];
        double a_air[nnz],   da_air[nnz],   as_air[nnz];
        double Kr_air=0.0, dKr_air=0.0;
        double pdeRes_air=0.0, Lstar_air[nDOF_test_element], subgrid_air=0.0;
        double tau_air=0.0, tau0_air=0.0, tau1_air=0.0, numDiff0_air=0.0, numDiff1_air=0.0;

        // double u_water = 0.0, grad_u_water[nSpace], grad_u_old_water[nSpace], m_water = 0.0, dm_water = 0.0, f_water[nSpace], df_water[nSpace], a_water[nnz], da_water[nnz], as_water[nnz], m_t_water = 0.0, dm_t_water = 0.0, pdeResidual_u_water = 0.0, Lstar_u_water[nDOF_test_element], subgridError_u_water = 0.0, tau_water = 0.0, tau0_water = 0.0, tau1_water = 0.0, numDiff0_water = 0.0, numDiff1_water = 0.0, jac_water[nSpace * nSpace], jacDet_water, jacInv_water[nSpace * nSpace], u_grad_trial_water[nDOF_trial_element * nSpace], u_test_dV_water[nDOF_trial_element], u_grad_test_dV_water[nDOF_test_element * nSpace], dV, x, y, z, xt, yt, zt, G[nSpace * nSpace], G_dd_G, tr_G, norm_Rv;
        // double u_air = 0.0, grad_u_air[nSpace], grad_u_old_air[nSpace], m_air = 0.0, dm_air = 0.0, f_air[nSpace], df_air[nSpace], a_air[nnz], da_air[nnz], as_air[nnz], m_t_air = 0.0, dm_t_air = 0.0, pdeResidual_u_air = 0.0, Lstar_u_air[nDOF_test_element], subgridError_u_air = 0.0, tau_air = 0.0, tau0_air = 0.0, tau1_air = 0.0, numDiff0_air = 0.0, numDiff1_air = 0.0, jac_air[nSpace * nSpace], jacDet_air, jacInv_air[nSpace * nSpace], u_grad_trial_air[nDOF_trial_element * nSpace], u_test_dV_air[nDOF_trial_element], u_grad_test_dV_air[nDOF_test_element * nSpace], dV, x, y, z, xt, yt, zt, G[nSpace * nSpace], G_dd_G, tr_G, norm_Rv;

        //
        //compute solution and gradients at quadrature points
        //
        ck.calculateMapping_element(eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y, z);
        ck.calculateMappingVelocity_element(eN, k, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), xt, yt, zt);
        //get the physical integration weight
        dV                = fabs(jacDet) * dV_ref.data()[k];
        q_dV.data()[eN_k] = dV;
        ck.calculateG(jacInv, G, G_dd_G, tr_G);
        //get the trial function gradients
        ck.gradTrialFromRef(&u_grad_trial_ref.data()[k * nDOF_trial_element * nSpace], jacInv, u_grad_trial);
        
        //get the solution for both phase
        ck.valFromDOF(u_dof_water.data(),                             // DOFs for WATER
                      &u_l2g.data()[eN_nDOF_trial_element],     // map for WATER
                      &u_trial_ref.data()[k * nDOF_trial_element],
                      u_water);
        ck.valFromDOF(u_dof_air.data(),                             // DOFs for AIR
                      &u_l2g.data()[eN_nDOF_trial_element],     // map for AIR
                      &u_trial_ref.data()[k * nDOF_trial_element],
                      u_air);
                        
        //get the solution gradient for both phases

        ck.gradFromDOF(u_dof_water.data(),
                       &u_l2g.data()[eN_nDOF_trial_element], 
                       u_grad_trial, 
                       grad_u_water);

        ck.gradFromDOF(u_dof_air.data(),
                       &u_l2g.data()[eN_nDOF_trial_element], 
                       u_grad_trial, 
                       grad_u_air);
        //get the solution
        //ck.valFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], u);
        //get the solution gradients
  //      ck.gradFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], u_grad_trial, grad_u);
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) {
          u_test_dV[j] = u_test_ref.data()[k * nDOF_trial_element + j] * dV;
          for (int I = 0; I < nSpace; I++) {
            u_grad_test_dV[j * nSpace + I] = u_grad_trial[j * nSpace + I] * dV; //cek warning won't work for Petrov-Galerkin
          }
        }
        //
        //calculate pde coefficients at quadrature points
        //
        //double Kr_water, dKr_water;
        double Sw=0.0, Sg=0.0;

       // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
       //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u, m, dm, f, df, a, da, as, Kr, dKr);
        
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                              0,//phase k=0 
                              rho_water, beta_water, 
                              gravity.data(),
                              alpha.data()[elementMaterialTypes.data()[eN]],
                              n.data()[elementMaterialTypes.data()[eN]],
                              thetaR.data()[elementMaterialTypes.data()[eN]], thetaSR.data()[elementMaterialTypes.data()[eN]],
                              &KWs.data()[elementMaterialTypes.data()[eN] * nnz],
                              u_water,  m_water, dm_water, 
                              f_water, df_water, 
                              a_water, da_water, 
                              as_water, 
                              Kr_water, dKr_water, 
                              PSK_TYPE,   // 0: VG_PSK, 1: BC_PSK
                              Sw, Sg,       
                              BC_entry_head,     // h_e (only used if BC_PSK)
                              BC_lambda);        // λ   (only used if BC_PSK)
                             //

        
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                              1,//phase k=0 
                              rho_air, beta_air, 
                              gravity.data(),
                              alpha.data()[elementMaterialTypes.data()[eN]],
                              n.data()[elementMaterialTypes.data()[eN]],
                              thetaR.data()[elementMaterialTypes.data()[eN]], thetaSR.data()[elementMaterialTypes.data()[eN]],
                              &KWs.data()[elementMaterialTypes.data()[eN] * nnz],
                              u_air,  m_air, dm_air, 
                              f_air, df_air, 
                              a_air, da_air, 
                              as_air, 
                              Kr_air, dKr_air, 
                              PSK_TYPE,   // 0: VG_PSK, 1: BC_PSK
                              Sw, Sg,       
                              BC_entry_head,     // h_e (only used if BC_PSK)
                              BC_lambda);        // λ   (only used if BC_PSK)
                             //
        //calculate time derivative at quadrature points
        //
        ck.bdf(alphaBDF, q_m_betaBDF_water.data()[eN_k], m_water, dm_water, m_t_water, dm_t_water);
        ck.bdf(alphaBDF, q_m_betaBDF_air.data()[eN_k], m_air, dm_air, m_t_air, dm_t_air);
        
        //
        //calculate subgrid error (strong residual and adjoint)
        //
        //calculate strong residual
        pdeResidual_u_water = ck.Mass_strong(m_t_water) + ck.Advection_strong(df_water, grad_u_water);      
        pdeResidual_u_air = ck.Mass_strong(m_t_air) + ck.Advection_strong(df_air, grad_u_air);

        //calculate adjoint
        for (int i = 0; i < nDOF_test_element; i++) {
          int i_nSpace = i * nSpace;
          Lstar_u_water[i]   = ck.Advection_adjoint(df_water, &u_grad_test_dV[i_nSpace]);
          Lstar_u_air[i]   = ck.Advection_adjoint(df_air, &u_grad_test_dV[i_nSpace]);
        }
        //calculate tau and tau*Res
        //phase water
        calculateSubgridError_tau(elementDiameter[eN], dm_t_water, df_water, cfl[eN_k], tau0_water);
        calculateSubgridError_tau(Ct_sge, G, dm_t_water, df_water, tau1_water, cfl[eN_k]);
        tau_water = useMetrics * tau1_water + (1.0 - useMetrics) * tau0_water;
        subgridError_u_water = -tau_water * pdeResidual_u_water;
        
                //phase air
        calculateSubgridError_tau(elementDiameter[eN], dm_t_air, df_air, cfl[eN_k], tau0_air);
        calculateSubgridError_tau(Ct_sge, G, dm_t_air, df_air, tau1_air, cfl[eN_k]);
        tau_air = useMetrics * tau1_air + (1.0 - useMetrics) * tau0_air;
        subgridError_u_air = -tau_air * pdeResidual_u_air;
        
        //subgridError_u = -tau * pdeResidual_u;
        //
        //calculate shock capturing diffusion
        //
        //phase water
        ck.calculateNumericalDiffusion(shockCapturingDiffusion, elementDiameter[eN], pdeResidual_u_water, grad_u_water, numDiff0_water);
        ck.calculateNumericalDiffusion(shockCapturingDiffusion, sc_uref, sc_alpha, G, G_dd_G, pdeResidual_u_water, grad_u_water, numDiff1_water);
        q_numDiff_u[eN_k] = useMetrics * numDiff1_water + (1.0 - useMetrics) * numDiff0_water;
        //phase air
        ck.calculateNumericalDiffusion(shockCapturingDiffusion, elementDiameter[eN], pdeResidual_u_air, grad_u_air, numDiff0_air);
        ck.calculateNumericalDiffusion(shockCapturingDiffusion, sc_uref, sc_alpha, G, G_dd_G, pdeResidual_u_air, grad_u_air, numDiff1_air);
        q_numDiff_u[eN_k] = useMetrics * numDiff1_air + (1.0 - useMetrics) * numDiff0_air;        
        
        //
        //update element residual
        //
        for (int i = 0; i < nDOF_test_element; i++) {
          int eN_k_i = eN_k * nDOF_test_element + i, eN_k_i_nSpace = eN_k_i * nSpace, i_nSpace = i * nSpace;
          elementResidual_u_water[i] += ck.Mass_weak(m_t_water, u_test_dV[i]) + ck.Advection_weak(f_water, &u_grad_test_dV[i_nSpace]) + ck.Diffusion_weak(a_rowptr.data(), a_colind.data(), a_water, grad_u_water, &u_grad_test_dV[i_nSpace]) + VMS * ck.SubgridError(subgridError_u_water, Lstar_u_water[i]) + VMS * ck.NumericalDiffusion(q_numDiff_u_last_water[eN_k], grad_u_water, &u_grad_test_dV[i_nSpace]);
          elementResidual_u_air[i] += ck.Mass_weak(m_t_air, u_test_dV[i]) + ck.Advection_weak(f_air, &u_grad_test_dV[i_nSpace]) + ck.Diffusion_weak(a_rowptr.data(), a_colind.data(), a_air, grad_u_air, &u_grad_test_dV[i_nSpace]) + VMS * ck.SubgridError(subgridError_u_air, Lstar_u_air[i]) + VMS * ck.NumericalDiffusion(q_numDiff_u_last_air[eN_k], grad_u_air, &u_grad_test_dV[i_nSpace]);

        } //i
        //phase water
        q_m_water.data()[eN_k] = m_water;
        q_u_water.data()[eN_k] = u_water;
        //phase air
        q_m_air.data()[eN_k] = m_air;
        q_u.data()[eN_k] = u_air;
      }
      //
      //load element into global residual and save element residual
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        globalResidual_water.data()[offset_u_water + stride_u_water * u_l2g.data()[eN_i]] += elementResidual_u_water[i];
        globalResidual_air.data()[offset_u_air + stride_u_air * u_l2g.data()[eN_i]] += elementResidual_u_air[i];
      } //i
    } //elements
    //
    //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
    //
    //ebNE is the Exterior element boundary INdex
    //ebN is the element boundary INdex
    //eN is the element index
    for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) {
      int    ebN = exteriorElementBoundariesArray.data()[ebNE], eN = elementBoundaryElementsArray.data()[ebN * 2 + 0], ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN * 2 + 0], eN_nDOF_trial_element = eN * nDOF_trial_element;
      double elementResidual_u_water[nDOF_test_element];
      double elementResidual_u_air[nDOF_test_element];
      for (int i = 0; i < nDOF_test_element; i++) { elementResidual_u_water[i] = 0.0; }
      for (int i = 0; i < nDOF_test_element; i++) { elementResidual_u_air[i] = 0.0; }
      for (int kb = 0; kb < nQuadraturePoints_elementBoundary; kb++) {
        int    ebNE_kb = ebNE * nQuadraturePoints_elementBoundary + kb, ebNE_kb_nSpace = ebNE_kb * nSpace, ebN_local_kb = ebN_local * nQuadraturePoints_elementBoundary + kb, ebN_local_kb_nSpace = ebN_local_kb * nSpace;
        double jac_ext[nSpace * nSpace], jacDet_ext_water, jacInv_ext_water[nSpace * nSpace], boundaryJac_water[nSpace * (nSpace - 1)];       
        double u_ext_water = 0.0, grad_u_ext_water[nSpace], m_ext_water = 0.0, dm_ext_water = 0.0, 
               f_ext_water[nSpace], df_ext_water[nSpace], a_ext_water[nnz], da_ext_water[nnz], as_ext_water[nnz], flux_ext_water = 0.0,
               //anb_seepage_flux=0.0, // for flux calculation
               bc_u_ext_water = 0.0, bc_grad_u_ext_water[nSpace], bc_m_ext_water = 0.0, bc_dm_ext_water = 0.0, 
               bc_f_ext_water[nSpace], bc_df_ext_water[nSpace], bc_a_ext_water[nnz], bc_da_ext_water[nnz], bc_as_ext_water[nnz];
        
        double u_ext_air = 0.0, grad_u_ext_air[nSpace], m_ext_air = 0.0, dm_ext_air = 0.0, 
               f_ext_air[nSpace], df_ext_air[nSpace], a_ext_air[nnz], da_ext_air[nnz], as_ext_air[nnz], flux_ext_air = 0.0,
               //anb_seepage_flux=0.0, // for flux calculation
               bc_u_ext_air = 0.0, bc_grad_u_ext_air[nSpace], bc_m_ext_air = 0.0, bc_dm_ext_air = 0.0, 
               bc_f_ext_air[nSpace], bc_df_ext_air[nSpace], bc_a_ext_air[nnz], bc_da_ext_air[nnz], bc_as_ext_air[nnz];         
        
          //
        //calculate the solution and gradients at quadrature points
        //
        //compute information about mapping from reference element to physical element
        ck.calculateMapping_elementBoundary(eN, ebN_local, kb, ebN_local_kb, mesh_dof.data(), mesh_l2g.data(), mesh_trial_trace_ref.data(), mesh_grad_trial_trace_ref.data(), boundaryJac_ref.data(), jac_ext, jacDet_ext, jacInv_ext, boundaryJac, metricTensor, metricTensorDetSqrt,
                                            normal_ref.data(), normal, x_ext, y_ext, z_ext);
        ck.calculateMappingVelocity_elementBoundary(eN, ebN_local, kb, ebN_local_kb, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_trace_ref.data(), xt_ext, yt_ext, zt_ext, normal, boundaryJac, metricTensor, integralScaling);
        dS = ((1.0 - MOVING_DOMAIN) * metricTensorDetSqrt + MOVING_DOMAIN * integralScaling) * dS_ref.data()[kb];
        //get the metric tensor
        //cek todo use symmetry
        ck.calculateG(jacInv_ext, G, G_dd_G, tr_G);
        //compute shape and solution information
        //shape
        ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace * nDOF_trial_element], jacInv_ext, u_grad_trial_trace);
        //solution 
        ck.valFromDOF(u_dof_water.data(), 
                      &u_l2g.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_water);
        ck.valFromDOF(u_dof_air.data(), 
                      &u_l2g.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_air);
        
        //Gradient
        ck.gradFromDOF(u_dof_water.data(), 
                       &u_l2g.data()[eN_nDOF_trial_element], 
                       u_grad_trial_trace, 
                       grad_u_ext_water);
        ck.gradFromDOF(u_dof_air.data(), 
                       &u_l2g.data()[eN_nDOF_trial_element], 
                       u_grad_trial_trace, 
                       grad_u_ext_air);
               
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) { u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb * nDOF_test_element + j] * dS; }
        //
        //load the boundary values
        //
        bc_u_ext_water = isDOFBoundary_u_water.data()[ebNE_kb] * ebqe_bc_u_ext_water.data()[ebNE_kb] + (1 - isDOFBoundary_u.data()[ebNE_kb]) * u_ext_water;
        bc_u_ext_air = isDOFBoundary_u_air.data()[ebNE_kb] * ebqe_bc_u_ext_air.data()[ebNE_kb] + (1 - isDOFBoundary_u.data()[ebNE_kb]) * u_ext;
      
        //
        //calculate the pde coefficients using the solution and the boundary values for the solution
        //
        double Kr, dKr;
        // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
        //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u_ext, m_ext, dm_ext, f_ext, df_ext, a_ext, da_ext, as_ext, Kr, dKr);
        // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
        //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], bc_u_ext, bc_m_ext, bc_dm_ext, bc_f_ext, bc_df_ext, bc_a_ext, bc_da_ext, bc_as_ext, Kr, dKr);
        

        // int num_phase= 2;

        // for (int i=0; i<num_phase, i++){

      //}
      double Sw_ext=0.0;
      double Sg_ext =0.0;

        //phase water
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             0, 
                             rho_water, beta_water, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_ext_water, m_ext_water, dm_ext_water, 
                             f_ext_water, df_ext_water, 
                             a_ext_water, da_ext_water, 
                             as_ext_water, 
                             Kr, dKr,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK));
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), 
                            0,
                            rho_water, beta_water, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            bc_u_ext_water, bc_m_ext_water, bc_dm_ext_water, 
                            bc_f_ext_water, bc_df_ext_water, 
                            bc_a_ext_water, bc_da_ext_water, bc_as_ext_water, 
                            Kr, dKr,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw_ext, Sg_ext,
                            BC_entry_head,     // h_e (only used if BC_PSK)
                            BC_lambda);        // λ   (only used if BC_PSK));
        //phase air        
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             1, 
                             rho_air, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_ext_air, m_ext_air, dm_ext_air, 
                             f_ext_air, df_ext_air, 
                             a_ext_air, da_ext_air, 
                             as_ext_air, 
                             Kr, dKr,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw, Sg
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK));
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), 
                            1,
                            rho_air, beta_air, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            bc_u_ext_air, bc_m_ext_air, bc_dm_ext_air, 
                            bc_f_ext_air, bc_df_ext_air, 
                            bc_a_ext_air, bc_da_ext_air, bc_as_ext_air, 
                            Kr, dKr,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw, Sg,
                            BC_entry_head,     // h_e (only used if BC_PSK)
                            BC_lambda);        // λ   (only used if BC_PSK));
                            //
        //calculate the numerical fluxes
        //
        exteriorNumericalFlux(ebqe_bc_flux_ext_water[ebNE_kb], a_rowptr.data(), a_colind.data(),
                              isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                              isDOFBoundary_u_water.data()[ebNE_kb], normal, 
                              bc_u_ext_water, a_ext_water, grad_u_ext_water, 
                              u_ext_water, f_ext_water,
                              ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                              flux_ext_water);
        
        ebqe_flux_water.data()[ebNE_kb] = flux_ext_water;
        
        exteriorNumericalFlux(ebqe_bc_flux_ext_air[ebNE_kb], a_rowptr.data(), a_colind.data(),
                              isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                              isDOFBoundary_u_air.data()[ebNE_kb], normal, 
                              bc_u_ext_air, a_ext_air, grad_u_ext_air, 
                              u_ext_air, f_ext_air,
                              ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                              flux_ext_air);
        
        ebqe_flux_air.data()[ebNE_kb] = flux_ext_air;

        //anb_seepage_flux             = seepagefluxcalculator(anb_seepage_flux, isSeepageFace.data()[ebNE], dS, flux_ext);
        //anb_seepage_flux_n.data()[0] = anb_seepage_flux;
        ebqe_u.data_water()[ebNE_kb]     = u_ext_water;
        ebqe_u.data_air()[ebNE_kb]       = u_ext_air;
        
        //
        //update residuals
        //
        for (int i = 0; i < nDOF_test_element; i++) {
          elementResidual_u_water[i] += ck.ExteriorElementBoundaryFlux(flux_ext_water, u_test_dS[i]);
          elementResidual_u_air[i] += ck.ExteriorElementBoundaryFlux(flux_ext_air, u_test_dS[i]);
          
        } //i
      } //kb

      //
      //update the element and global residual storage
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        globalResidual_water.data()[offset_u_water + stride_u_water * u_l2g.data()[eN_i]] += elementResidual_u_water[i];
        globalResidual_air.data()[offset_u_air + stride_u_air * u_l2g.data()[eN_i]] += elementResidual_u_air[i];        
      } //i
    } //ebNE
  }

  void calculateJacobian(arguments_dict &args)
  {
    xt::pyarray<double> &mesh_trial_ref            = args.array<double>("mesh_trial_ref");
    xt::pyarray<double> &mesh_grad_trial_ref       = args.array<double>("mesh_grad_trial_ref");
    xt::pyarray<double> &mesh_dof                  = args.array<double>("mesh_dof");
    xt::pyarray<double> &mesh_velocity_dof         = args.array<double>("mesh_velocity_dof");
    double               MOVING_DOMAIN             = args.scalar<double>("MOVING_DOMAIN");
    xt::pyarray<int>    &mesh_l2g                  = args.array<int>("mesh_l2g");
    xt::pyarray<double> &dV_ref                    = args.array<double>("dV_ref");
    xt::pyarray<double> &u_trial_ref               = args.array<double>("u_trial_ref");
    xt::pyarray<double> &u_grad_trial_ref          = args.array<double>("u_grad_trial_ref");
    xt::pyarray<double> &u_test_ref                = args.array<double>("u_test_ref");
    xt::pyarray<double> &u_grad_test_ref           = args.array<double>("u_grad_test_ref");
    xt::pyarray<double> &mesh_trial_trace_ref      = args.array<double>("mesh_trial_trace_ref");
    xt::pyarray<double> &mesh_grad_trial_trace_ref = args.array<double>("mesh_grad_trial_trace_ref");
    xt::pyarray<double> &dS_ref                    = args.array<double>("dS_ref");
    xt::pyarray<double> &u_trial_trace_ref         = args.array<double>("u_trial_trace_ref");
    xt::pyarray<double> &u_grad_trial_trace_ref    = args.array<double>("u_grad_trial_trace_ref");
    xt::pyarray<double> &u_test_trace_ref          = args.array<double>("u_test_trace_ref");
    xt::pyarray<double> &u_grad_test_trace_ref     = args.array<double>("u_grad_test_trace_ref");
    xt::pyarray<double> &normal_ref                = args.array<double>("normal_ref");
    xt::pyarray<double> &boundaryJac_ref           = args.array<double>("boundaryJac_ref");
    int                  nElements_global          = args.scalar<int>("nElements_global");
    xt::pyarray<double> &ebqe_penalty_ext          = args.array<double>("ebqe_penalty_ext");
    xt::pyarray<int>    &elementMaterialTypes      = args.array<int>("elementMaterialTypes");
    xt::pyarray<int>    &isSeepageFace             = args.array<int>("isSeepageFace");
    xt::pyarray<int>    &a_rowptr                  = args.array<int>("a_rowptr");
    xt::pyarray<int>    &a_colind                  = args.array<int>("a_colind");
    xt::pyarray<double> &gravity                   = args.array<double>("gravity");
    xt::pyarray<double> &alpha                     = args.array<double>("alpha");
    xt::pyarray<double> &n                         = args.array<double>("n");
    xt::pyarray<double> &thetaR                    = args.array<double>("thetaR");
    xt::pyarray<double> &thetaSR                   = args.array<double>("thetaSR");
    xt::pyarray<double> &KWs                       = args.array<double>("KWs");
    double               useMetrics                = args.scalar<double>("useMetrics");
    double               alphaBDF                  = args.scalar<double>("alphaBDF");
    int                  lag_shockCapturing        = args.scalar<int>("lag_shockCapturing");
    double               shockCapturingDiffusion   = args.scalar<double>("shockCapturingDiffusion");
    // VMS
    double               VMS                       = args.scalar<double>("VMS");
    xt::pyarray<double> &elementDiameter           = args.array<double>("elementDiameter");
    xt::pyarray<double> &cfl                       = args.array<double>("cfl");
     xt::pyarray<int>   &csrRowIndeces_u_u         = args.array<int>("csrRowIndeces_u_u");
    xt::pyarray<int>    &csrColumnOffsets_u_u      = args.array<int>("csrColumnOffsets_u_u");
    int                  nExteriorElementBoundaries_global          = args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int>    &exteriorElementBoundariesArray             = args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int>    &elementBoundaryElementsArray               = args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int>    &elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
    xt::pyarray<int>    &csrColumnOffsets_eb_u_u                    = args.array<int>("csrColumnOffsets_eb_u_u");
    int                  LUMPED_MASS_MATRIX                         = args.scalar<int>("LUMPED_MASS_MATRIX");
   
    //////////////////////////Two phase//////////////////////////////////////////////////
    double               rho_water                                  = args.scalar<double>("rho_water");
    double               rho_air                                    = args.scalar<double>("rho_air");
    double               beta_water                                 = args.scalar<double>("beta_water");
    double               beta_air                                   = args.scalar<double>("beta_air");
    xt::pyarray<int>    &u_l2g                                      = args.array<int>("u_l2g");
   // xt::pyarray<int>    &u_l2g_air                                  = args.array<int>("u_l2g_air");   
    xt::pyarray<double> &u_dof_water                                = args.array<double>("u_dof_water");
    xt::pyarray<double> &u_dof_air                                  = args.array<double>("u_dof_air");
    xt::pyarray<double> &velocity_water                             = args.array<double>("velocity_water");
    xt::pyarray<double> &velocity_air                               = args.array<double>("velocity_air");
    xt::pyarray<double> &q_m_betaBDF_water                          = args.array<double>("q_m_betaBDF_water");
    xt::pyarray<double> &q_m_betaBDF_air                            = args.array<double>("q_m_betaBDF_air");
    xt::pyarray<double> &q_numDiff_u_water                          = args.array<double>("q_numDiff_u_water");
    xt::pyarray<double> &q_numDiff_u_air                            = args.array<double>("q_numDiff_u_air");
    xt::pyarray<double> &q_numDiff_u_last_water                     = args.array<double>("q_numDiff_u_last_water");
    xt::pyarray<double> &q_numDiff_u_last_air                       = args.array<double>("q_numDiff_u_last_air");
    xt::pyarray<double> &globalJacobian_water                       = args.array<double>("globalJacobian_water");
    xt::pyarray<double> &globalJacobian_air                         = args.array<double>("globalJacobian_air");
    xt::pyarray<double> &ebqe_velocity_ext_water                    = args.array<double>("ebqe_velocity_ext_water");
    xt::pyarray<double> &ebqe_velocity_ext_air                      = args.array<double>("ebqe_velocity_ext_air");
    xt::pyarray<int>    &isDOFBoundary_u_water                      = args.array<int>("isDOFBoundary_u_water");
    xt::pyarray<int>    &isDOFBoundary_u_air                        = args.array<int>("isDOFBoundary_u_air");
    xt::pyarray<double> &ebqe_bc_u_ext_water                        = args.array<double>("ebqe_bc_u_ext_water");
    xt::pyarray<double> &ebqe_bc_u_ext_air                          = args.array<double>("ebqe_bc_u_ext_air");
    xt::pyarray<int>    &isFluxBoundary_u_water                     = args.array<int>("isFluxBoundary_u_water");  
    xt::pyarray<int>    &isFluxBoundary_u_air                       = args.array<int>("isFluxBoundary_u_air");
    xt::pyarray<double> &ebqe_bc_flux_ext_water                     = args.array<double>("ebqe_bc_flux_ext_water");
    xt::pyarray<double> &ebqe_bc_flux_ext_air                       = args.array<double>("ebqe_bc_flux_ext_air");
    
    assert(a_rowptr.data()[nSpace] == nnz);
    assert(a_rowptr.data()[nSpace] == nSpace);
    double Ct_sge = 4.0;
    ////////////////////Added for Brooks Correy///////////////////////////////
    PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
    double BC_entry_head = args.scalar<double>("BC_entry_head");
    double BC_lambda     = args.scalar<double>("BC_lambda");

    //
    //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
    //
    for (int eN = 0; eN < nElements_global; eN++) {
      double elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++) {
        for (int j = 0; j < nDOF_trial_element; j++) { elementJacobian_u_u[i][j] = 0.0; }
      }
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        int eN_k                  = eN * nQuadraturePoints_element + k, //index to a scalar at a quadrature point
          eN_k_nSpace             = eN_k * nSpace,
          eN_nDOF_trial_element   = eN * nDOF_trial_element; //index to a vector at a quadrature point

        //declare local storage
        double  dV, x, y, z, xt, yt, zt;
        double  jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace];
        double  G[nSpace * nSpace], G_dd_G, tr_G;
        double  u_grad_trial[nDOF_trial_element * nSpace], u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element * nSpace]
        double  u_water = 0.0, grad_u_water[nSpace], m_water = 0.0, dm_water = 0.0, 
                f_water[nSpace], df_water[nSpace], a_water[nnz], da_water[nnz], as_water[nnz], m_t_water = 0.0, dm_t_water = 0.0, 
                dpdeResidual_u_u_water[nDOF_trial_element], Lstar_u_water[nDOF_test_element], dsubgridError_u_u_water[nDOF_trial_element], 
                tau_water = 0.0, tau0_water = 0.0, tau1_water = 0.0;
        double  u_air = 0.0, grad_u_air[nSpace], m_air = 0.0, dm_air = 0.0, 
                f_air[nSpace], df_air[nSpace], a_air[nnz], da_air[nnz], as_air[nnz], m_t_air = 0.0, dm_t_air = 0.0, 
                dpdeResidual_u_u_air[nDOF_trial_element], Lstar_u_air[nDOF_test_element], dsubgridError_u_u_air[nDOF_trial_element], 
                tau_air = 0.0, tau0_air = 0.0, tau1_air = 0.0;
        //
        //calculate solution and gradients at quadrature points
        //
        //get jacobian, etc for mapping reference element
        ck.calculateMapping_element(eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y, z);
        ck.calculateMappingVelocity_element(eN, k, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), xt, yt, zt);
        //get the physical integration weight
        dV = fabs(jacDet) * dV_ref.data()[k];
        ck.calculateG(jacInv, G, G_dd_G, tr_G);
        //get the trial function gradients
        ck.gradTrialFromRef(&u_grad_trial_ref.data()[k * nDOF_trial_element * nSpace], jacInv, u_grad_trial);
        //get the solution for both phases
        ck.valFromDOF(u_dof_water.data(), 
                      &u_l2g.data()[eN_nDOF_trial_element], 
                      &u_trial_ref.data()[k * nDOF_trial_element], 
                      u_water);
        
        ck.valFromDOF(u_dof_air.data(), 
                      &u_l2g.data()[eN_nDOF_trial_element], 
                      &u_trial_ref.data()[k * nDOF_trial_element], 
                      u_air);
        
        //get the solution gradients for both phases
        ck.gradFromDOF(u_dof_water.data(),
                       &u_l2g.data()[eN_nDOF_trial_element], 
                       u_grad_trial, 
                       grad_u_water);

        ck.gradFromDOF(u_dof_air.data(),
                       &u_l2g.data()[eN_nDOF_trial_element], 
                       u_grad_trial, 
                       grad_u_air);

                       //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) {
          u_test_dV[j] = u_test_ref.data()[k * nDOF_trial_element + j] * dV;
          for (int I = 0; I < nSpace; I++) {
            u_grad_test_dV[j * nSpace + I] = u_grad_trial[j * nSpace + I] * dV; //cek warning won't work for Petrov-Galerkin
          }
        }
        //
        //calculate pde coefficients and derivatives at quadrature points
        //
        double Kr_water, dKr_water;
        double Kr_air, dKr_air;
        // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
        //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u, m, dm, f, df, a, da, as, Kr, dKr);
        
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), 
                             0,
                             rho_water, beta_water, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_water, m_water, dm_water, 
                             f_water, df_water, 
                             a_water, da_water, 
                             as_water, Kr_water, dKr_water,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_water, Sg_water,
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK));

        evaluateCoefficients(a_rowptr.data(), a_colind.data(), 
                             1,
                             rho_air, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_air, m_air, dm_air, 
                             f_air, df_air, 
                             a_air, da_air, 
                             as_air, Kr_air, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_air, Sg_air,
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK));
        
                             //
        //calculate time derivatives
        //
        ck.bdf(alphaBDF, q_m_betaBDF_water.data()[eN_k], m_water, dm_water, m_t_water, dm_t_water);
        ck.bdf(alphaBDF, q_m_betaBDF_air.data()[eN_k], m_air, dm_air, m_t_air, dm_t_air);
        //
        //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
        //
        //calculate the adjoint times the test functions
        for (int i = 0; i < nDOF_test_element; i++) {
          int i_nSpace = i * nSpace;
          Lstar_u_water[i]   = ck.Advection_adjoint(df_water, &u_grad_test_dV[i_nSpace]);
          Lstar_u_air[i]   = ck.Advection_adjoint(df_air, &u_grad_test_dV[i_nSpace]); 
        }
        //calculate the Jacobian of strong residual
        for (int j = 0; j < nDOF_trial_element; j++) {
          int j_nSpace        = j * nSpace;
          dpdeResidual_u_u_water[j] = ck.MassJacobian_strong(dm_t_water, u_trial_ref[k * nDOF_trial_element + j]) + ck.AdvectionJacobian_strong(df_water, &u_grad_trial[j_nSpace]);
          dpdeResidual_u_u_air[j] = ck.MassJacobian_strong(dm_t_air, u_trial_ref[k * nDOF_trial_element + j]) + ck.AdvectionJacobian_strong(df_air, &u_grad_trial[j_nSpace]);        
        }
        //tau and tau*Res for phase 0 :: water

        calculateSubgridError_tau(elementDiameter[eN], dm_t_water, df_water, cfl[eN_k], tau0_water);
        calculateSubgridError_tau(Ct_sge, G, dm_t_water, df_water, tau1_water, cfl[eN_k]);
        tau_water = useMetrics * tau1_water + (1.0 - useMetrics) * tau0_water;
        
        //tau and tau*Res for phase 1 :: air

        calculateSubgridError_tau(elementDiameter[eN], dm_t_air, df_air, cfl[eN_k], tau0_air);
        calculateSubgridError_tau(Ct_sge, G, dm_t_air, df_air, tau1_air, cfl[eN_k]);
        tau_air = useMetrics * tau1_air + (1.0 - useMetrics) * tau0_air;
        
        
        for (int j = 0; j < nDOF_trial_element; j++) {
          dsubgridError_u_u_water[j] = -tau_water * dpdeResidual_u_u_water[j];
          dsubgridError_u_u_air[j] = -tau_air * dpdeResidual_u_u_air[j];
        }
        for (int i = 0; i < nDOF_test_element; i++) {
          for (int j = 0; j < nDOF_trial_element; j++) {
            int j_nSpace = j * nSpace;
            int i_nSpace = i * nSpace;
            elementJacobian_u_u_water[i][j] += ck.MassJacobian_weak(dm_t_water, u_trial_ref.data()[k * nDOF_trial_element + j], u_test_dV[i]) + 
                                               ck.AdvectionJacobian_weak(df_water, u_trial_ref.data()[k * nDOF_trial_element + j], &u_grad_test_dV[i_nSpace]) +
                                               ck.DiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), 
                                               a_water, da_water, grad_u_water, &u_grad_test_dV[i_nSpace], 1.0, u_trial_ref.data()[k * nDOF_trial_element + j], 
                                               &u_grad_trial[j_nSpace]) + VMS * ck.SubgridErrorJacobian(dsubgridError_u_u_water[j], Lstar_u_water[i]) + 
                                               VMS * ck.NumericalDiffusionJacobian(q_numDiff_u_last_water[eN_k], &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
          
            elementJacobian_u_u_air[i][j] += ck.MassJacobian_weak(dm_t_air, u_trial_ref.data()[k * nDOF_trial_element + j], u_test_dV[i]) + 
                                               ck.AdvectionJacobian_weak(df_air, u_trial_ref.data()[k * nDOF_trial_element + j], &u_grad_test_dV[i_nSpace]) +
                                               ck.DiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), 
                                               a_air, da_air, grad_u_air, &u_grad_test_dV[i_nSpace], 1.0, u_trial_ref.data()[k * nDOF_trial_element + j], 
                                               &u_grad_trial[j_nSpace]) + VMS * ck.SubgridErrorJacobian(dsubgridError_u_u_air[j], Lstar_u_air[i]) + 
                                               VMS * ck.NumericalDiffusionJacobian(q_numDiff_u_last_air[eN_k], &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);         
                                              } //j
        } //i
      } //k
      //
      //load into element Jacobian into global Jacobian
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        for (int j = 0; j < nDOF_trial_element; j++) {
          int eN_i_j = eN_i * nDOF_trial_element + j;
          globalJacobian_water.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u_water[i][j];
          globalJacobian_air.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u_air[i][j];

        } //j
      } //i
    } //elements
    //
    //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
    //
    for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) {
      int ebN = exteriorElementBoundariesArray.data()[ebNE];
      int eN = elementBoundaryElementsArray.data()[ebN * 2 + 0], ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN * 2 + 0], eN_nDOF_trial_element = eN * nDOF_trial_element;
      for (int kb = 0; kb < nQuadraturePoints_elementBoundary; kb++) {
        int ebNE_kb = ebNE * nQuadraturePoints_elementBoundary + kb, ebNE_kb_nSpace = ebNE_kb * nSpace, ebN_local_kb = ebN_local * nQuadraturePoints_elementBoundary + kb, ebN_local_kb_nSpace = ebN_local_kb * nSpace;
        double jac_ext[nSpace * nSpace], jacDet_ext, jacInv_ext[nSpace * nSpace], boundaryJac[nSpace * (nSpace - 1)], 
                metricTensor[(nSpace - 1) * (nSpace - 1)], metricTensorDetSqrt, dS, u_test_dS[nDOF_test_element], u_grad_trial_trace[nDOF_trial_element * nSpace], normal[3], 
                x_ext, y_ext, z_ext, xt_ext, yt_ext, zt_ext, integralScaling, G[nSpace * nSpace], G_dd_G, tr_G;
        double u_ext_water = 0.0, grad_u_ext_water[nSpace], m_ext_water = 0.0, dm_ext_water = 0.0, f_ext_water[nSpace], df_ext_water[nSpace], 
                a_ext_water[nnz], da_ext_water[nnz], as_ext_water[nnz], dflux_u_u_ext_water = 0.0, bc_u_ext_water = 0.0,
               //bc_grad_u_ext[nSpace],
                bc_m_ext_water = 0.0, bc_dm_ext_water = 0.0, bc_f_ext_water[nSpace], bc_df_ext_water[nSpace], bc_a_ext_water[nnz], bc_da_ext_water[nnz], 
                bc_as_ext_water[nnz], fluxJacobian_u_u_water[nDOF_trial_element];

        double u_ext_air = 0.0, grad_u_ext_air[nSpace], m_ext_air = 0.0, dm_ext_air = 0.0, f_ext_air[nSpace], df_ext_air[nSpace], 
                a_ext_air[nnz], da_ext_air[nnz], as_ext_air[nnz], dflux_u_u_ext_air = 0.0, bc_u_ext_air = 0.0,
               //bc_grad_u_ext[nSpace],
                bc_m_ext_air = 0.0, bc_dm_ext_air = 0.0, bc_f_ext_air[nSpace], bc_df_ext_air[nSpace], bc_a_ext_air[nnz], bc_da_ext_air[nnz], 
                bc_as_ext_air[nnz], fluxJacobian_u_u_air[nDOF_trial_element];
        //
        //calculate the solution and gradients at quadrature points
        //
        ck.calculateMapping_elementBoundary(eN, ebN_local, kb, ebN_local_kb, mesh_dof.data(), mesh_l2g.data(), mesh_trial_trace_ref.data(), mesh_grad_trial_trace_ref.data(), boundaryJac_ref.data(), jac_ext, jacDet_ext, jacInv_ext, boundaryJac, metricTensor, metricTensorDetSqrt,
                                            normal_ref.data(), normal, x_ext, y_ext, z_ext);
        ck.calculateMappingVelocity_elementBoundary(eN, ebN_local, kb, ebN_local_kb, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_trace_ref.data(), xt_ext, yt_ext, zt_ext, normal, boundaryJac, metricTensor, integralScaling);
        dS = ((1.0 - MOVING_DOMAIN) * metricTensorDetSqrt + MOVING_DOMAIN * integralScaling) * dS_ref.data()[kb];
        ck.calculateG(jacInv_ext, G, G_dd_G, tr_G);
        //compute shape and solution information
        //shape
        ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace * nDOF_trial_element], jacInv_ext, u_grad_trial_trace);
        //solution for both phases
        ck.valFromDOF(u_dof_water.data(), 
                      &u_l2g.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_water);
        
        ck.valFromDOF(u_dof_air.data(), 
                      &u_l2g.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_air);


        // Gradient for both phases
        ck.gradFromDOF(u_dof_water.data(), 
                        &u_l2g.data()[eN_nDOF_trial_element], 
                        u_grad_trial_trace, 
                        grad_u_ext_water);
        ck.gradFromDOF(u_dof_air.data(), 
                        &u_l2g.data()[eN_nDOF_trial_element], 
                        u_grad_trial_trace, 
                        grad_u_ext_air);

      
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) { u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb * nDOF_test_element + j] * dS; }
        //
        //load the boundary values
        //
        bc_u_ext_water = isDOFBoundary_u_water.data()[ebNE_kb] * ebqe_bc_u_ext_water.data()[ebNE_kb] + (1 - isDOFBoundary_u_water.data()[ebNE_kb]) * u_ext_water;
        bc_u_ext_air = isDOFBoundary_u_air.data()[ebNE_kb] * ebqe_bc_u_ext_air.data()[ebNE_kb] + (1 - isDOFBoundary_u_air.data()[ebNE_kb]) * u_ext_air;
        
        //
        //calculate the internal and external trace of the pde coefficients
        //
        double Kr_water, dKr_water;
        double Kr_air, dKr_air;
        double Sw=0.0, Sg=0.0;
        
        // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
        //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u_ext, m_ext, dm_ext, f_ext, df_ext, a_ext, da_ext, as_ext, Kr, dKr);
        // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
        //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], bc_u_ext, bc_m_ext, bc_dm_ext, bc_f_ext, bc_df_ext, bc_a_ext, bc_da_ext, bc_as_ext, Kr, dKr);
        
        // int num_phase= 2;

        // //Need to change variables for both phases

        // for (int i=0; i<num_phase, i++){

        //phase water
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             0, 
                             rho_water, beta_water, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_ext_water, m_ext_water, dm_ext_water, 
                             f_ext_water, df_ext_water, 
                             a_ext_water, da_ext_water, 
                             as_ext_water, 
                             Kr, dKr,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK));
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), 
                            0,
                            rho_water, beta_water, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            bc_u_ext_water, bc_m_ext_water, bc_dm_ext_water, 
                            bc_f_ext_water, bc_df_ext_water, 
                            bc_a_ext_water, bc_da_ext_water, bc_as_ext_water, 
                            Kr, dKr,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw_ext, Sg_ext,
                            BC_entry_head,     // h_e (only used if BC_PSK)
                            BC_lambda);        // λ   (only used if BC_PSK));
        //phase air        
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             1, 
                             rho_air, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_ext_air, m_ext_air, dm_ext_air, 
                             f_ext_air, df_ext_air, 
                             a_ext_air, da_ext_air, 
                             as_ext_air, 
                             Kr, dKr,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw, Sg
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK));
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), 
                            1,
                            rho_air, beta_air, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            bc_u_ext_air, bc_m_ext_air, bc_dm_ext_air, 
                            bc_f_ext_air, bc_df_ext_air, 
                            bc_a_ext_air, bc_da_ext_air, bc_as_ext_air, 
                            Kr, dKr,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw, Sg,
                            BC_entry_head,     // h_e (only used if BC_PSK)
                            BC_lambda);        // λ   (only used if BC_PSK));
                            //

        // }
        
                             //
        //calculate the flux jacobian
        //
        for (int j = 0; j < nDOF_trial_element; j++) {
          exteriorNumericalFluxJacobian(a_rowptr.data(), a_colind.data(), 
                                        isDOFBoundary_u_water.data()[ebNE_kb], 
                                        normal, 
                                        a_ext_water, da_ext_water, grad_u_ext_water, 
                                        &u_grad_trial_trace[j * nSpace], 
                                        df_ext_water, 
                                        u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u_water[j]);
          exteriorNumericalFluxJacobian(a_rowptr.data(), a_colind.data(), 
                                        isDOFBoundary_u_air.data()[ebNE_kb], 
                                        normal, 
                                        a_ext_air, da_ext_air, grad_u_ext_air, 
                                        &u_grad_trial_trace[j * nSpace], 
                                        df_ext_air, 
                                        u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u_air[j]);


        } //j
        //
        //update the global Jacobian from the flux Jacobian
        //
        for (int i = 0; i < nDOF_test_element; i++) {
          int eN_i = eN * nDOF_test_element + i;
          for (int j = 0; j < nDOF_trial_element; j++) {
            int ebN_i_j = ebN * 4 * nDOF_test_X_trial_element + i * nDOF_trial_element + j;
            globalJacobian_water.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_water[j] * u_test_dS[i];
            globalJacobian_air.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_air[j] * u_test_dS[i];
          } //j
        } //i
      } //kb
    } //ebNE
  } //computeJacobian

  void FCTStep(arguments_dict &args)
  {
    xt::pyarray<double> &bc_mask                   = args.array<double>("bc_mask");
    int                  NNZ                       = args.scalar<int>("NNZ");     //number on non-zero entries on sparsity pattern
    int                  numDOFs                   = args.scalar<int>("numDOFs"); //number of DOFs
    double               dt                        = args.scalar<double>("dt");
    xt::pyarray<double> &ML                        = args.array<double>("ML"); //lumped mass matrix (as vector)
    xt::pyarray<double> &mn                        = args.array<double>("mn");               //DOFs of solution at time tn
    xt::pyarray<double> &mHigh                     = args.array<double>("mHigh");               //DOFs of high order solution at tnp1
    xt::pyarray<double> &mLow                      = args.array<double>("mLow");
    xt::pyarray<double> &mDotHigh                     = args.array<double>("mDotHigh");               //DOFs of high order solution at tnp1
    xt::pyarray<double> &mDotLow                      = args.array<double>("mDotLow");
    xt::pyarray<double> &limited_solution          = args.array<double>("limited_solution");
    xt::pyarray<int>    &csrRowIndeces_DofLoops    = args.array<int>("csrRowIndeces_DofLoops");    //csr row indeces
    xt::pyarray<int>    &csrColumnOffsets_DofLoops = args.array<int>("csrColumnOffsets_DofLoops"); //csr column offsets
    xt::pyarray<double> &MC                        = args.array<double>("MC");             //mass matrix
    xt::pyarray<double> &dt_times_fH_minus_fL      = args.array<double>("dt_times_fH_minus_fL");   //low minus high order dissipative matrices
    xt::pyarray<double> &min_m_bc                  = args.array<double>("min_m_bc");               //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
    xt::pyarray<double> &max_m_bc                  = args.array<double>("max_m_bc");
    xt::pyarray<double> &fluxCorrection                  = args.array<double>("fluxCorrection");
    //flags
    int                  LUMPED_MASS_MATRIX        = args.scalar<int>("LUMPED_MASS_MATRIX");
    int                  MONOLITHIC                = args.scalar<int>("MONOLITHIC");
    double               Rpos[numDOFs], Rneg[numDOFs];
    double               FluxCorrectionMatrix[NNZ];
    double               mDot[numDOFs];

    //////////////////
    // LOOP in DOFs //
    //////////////////
    int ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      mDot[i] = (mLow.data()[i] - mn.data()[i])/dt;
      //cek todo: add boundary data--these are just initialized
      //will need to pass p_bc at DOF and calc M
      double mini=min_m_bc.data()[i], maxi=max_m_bc.data()[i];
      //we're doing local FCT
      //if (GLOBAL_FCT == 1) {
      //  mini = 0.;
      //  maxi = 1.;
      //}

      double Pposi = 0, Pnegi = 0;
      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops.data()[offset];
        ////////////////////////
        // COMPUTE THE BOUNDS //
        ////////////////////////
        if (GLOBAL_FCT == 0) {
          if (MONOLITHIC == 0) {
            mini = fmin(mini, mLow[j]);
            maxi = fmax(maxi, mLow[j]);
          } else {
            mini = fmin(mini, mn.data()[j]);
            maxi = fmax(maxi, mn.data()[j]);
          }
        }
        // i-th row of flux correction matrix
        //double I_plus_ML_minus_MC = (i == j ? 1. : 0.) * (1. + ML.data()[i]) - MC.data()[ij];
        //mDot[i] += I_plus_ML_minus_MC * (mHigh.data()[j] - mn.data()[j]) / ML.data()[i];
        mDot[j] = (mLow.data()[j] - mn.data()[j])/dt;
        if (MONOLITHIC == 0) {
          FluxCorrectionMatrix[ij] = (LUMPED_MASS_MATRIX == 1 ? 0. : 1.) * dt * MC.data()[ij] * (mDotLow.data()[i] - mDotLow.data()[j]) + dt_times_fH_minus_fL.data()[ij];
        } else {
          FluxCorrectionMatrix[ij] = dt_times_fH_minus_fL.data()[ij];
        }
        ///////////////////////
        // COMPUTE P VECTORS //
        ///////////////////////
        Pposi += FluxCorrectionMatrix[ij] * ((FluxCorrectionMatrix[ij] > 0) ? 1. : 0.);
        Pnegi += FluxCorrectionMatrix[ij] * ((FluxCorrectionMatrix[ij] < 0) ? 1. : 0.);

        //update ij
        ij += 1;
      }
      ///////////////////////
      // COMPUTE Q VECTORS //
      ///////////////////////
      double gamma;
      double Qposi;
      double Qnegi;
      if (MONOLITHIC == 0) {
        Qposi = ML.data()[i] * (maxi - mLow[i]);
        Qnegi = ML.data()[i] * (mini - mLow[i]);
      } else {
        //cek todo: don't think this is right for Richards
        gamma = 10.0 * ML.data()[i];
        Qposi = fmin(0.5 * ML.data()[i] * (1.0 - mn.data()[i]), gamma * (maxi - mn[i]));
        Qnegi = fmax(0.5 * ML.data()[i] * (0.0 - mn.data()[i]), gamma * (mini - mn[i]));
      }
      ///////////////////////
      // COMPUTE R VECTORS //
      ///////////////////////
      Rpos[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
      Rneg[i] = ((Pnegi == 0) ? 1. : fmin(1.0, Qnegi / Pnegi));
    } // i DOFs

    //////////////////////
    // COMPUTE LIMITERS //
    //////////////////////
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      double ith_Limiter_times_FluxCorrectionMatrix = 0.;
      double alpha_fA, alpha_dot, beta_ij = 1.0;
      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops.data()[offset];
        alpha_fA     = ((FluxCorrectionMatrix[ij] > 0) ? fmin(Rpos[i], Rneg[j]) : fmin(Rneg[i], Rpos[j])) * FluxCorrectionMatrix[ij];
        alpha_dot    = fmin(1.0, beta_ij * fabs(alpha_fA) / MC.data()[ij] / fmax(1.0e-8, fabs(mDot[i] - mDot[j])));
        if (MONOLITHIC == 0) {
          ith_Limiter_times_FluxCorrectionMatrix += alpha_fA;
        } else {
          ith_Limiter_times_FluxCorrectionMatrix += alpha_fA + (LUMPED_MASS_MATRIX == 1 ? 0. : 1.) * dt * alpha_dot * MC.data()[ij] * (mDot[i] - mDot[j]);
        }
        ij += 1;
      }

      fluxCorrection.data()[i] = -ith_Limiter_times_FluxCorrectionMatrix*bc_mask[i]/dt;
      limited_solution.data()[i] = mLow[i] + 1. / ML.data()[i] * ith_Limiter_times_FluxCorrectionMatrix * bc_mask[i];

      //cek todo: double check that the below is not necesary. The limted_solution should already be within the bounds
      //Calculate the min and max mass bounds
      //double mMin = rho * thetaR.data()[elementMaterialTypes.data()[0]];
      //double mMax = rho * (thetaR.data()[elementMaterialTypes.data()[0]] + thetaSR.data()[elementMaterialTypes.data()[0]]);

      // Check if the limited mass is within bounds
      //if (limited_mass < mMin || limited_mass > mMax) {
      //  limited_solution.data()[i] = solL[i]; // Fallback to lower-order solution
      //} else {
      //  limited_solution.data()[i] = limited_mass; // Assign the limited mass
      //}
    }
  }

  void kth_FCT_step(arguments_dict &args)
  {
    int                  NNZ                       = args.scalar<int>("NNZ");     //number on non-zero entries on sparsity pattern
    int                  numDOFs                   = args.scalar<int>("numDOFs"); //number of DOFs
    int                  num_fct_iter              = args.scalar<int>("num_fct_iter");
    double               dt                        = args.scalar<double>("dt");
    xt::pyarray<double> &lumped_mass_matrix        = args.array<double>("lumped_mass_matrix"); //lumped mass matrix (as vector)
    xt::pyarray<double> &soln                      = args.array<double>("soln");               //DOFs of solution at time tn
    xt::pyarray<double> &pn                        = args.array<double>("pn");                 //DOFs of solution at time tn
    xt::pyarray<double> &solH                      = args.array<double>("solH");               //DOFs of high order solution at tnp1
    xt::pyarray<double> &uLow                      = args.array<double>("uLow");
    xt::pyarray<double> &uDotLow                   = args.array<double>("uDotLow");
    xt::pyarray<double> &dLow                      = args.array<double>("dLow");
    xt::pyarray<double> &solLim                    = args.array<double>("limited_solution");
    xt::pyarray<double> &MC                        = args.array<double>("MC");
    xt::pyarray<double> &ML                        = args.array<double>("ML");
    xt::pyarray<double> &FluxMatrix                = args.array<double>("FluxMatrix");
    xt::pyarray<double> &limitedFlux               = args.array<double>("limited_Flux");
    xt::pyarray<int>    &csrRowIndeces_DofLoops    = args.array<int>("csrRowIndeces_DofLoops");    //csr row indeces
    xt::pyarray<int>    &csrColumnOffsets_DofLoops = args.array<int>("csrColumnOffsets_DofLoops"); //csr column offsets
    xt::pyarray<double> &MassMatrix                = args.array<double>("MassMatrix");             //mass matrix
    xt::pyarray<double> &dt_times_fH_minus_fL      = args.array<double>("dt_times_fH_minus_fL");   //low minus high order dissipative matrices
    xt::pyarray<double> &min_m_bc                  = args.array<double>("min_m_bc");               //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
    xt::pyarray<double> &max_m_bc                  = args.array<double>("max_m_bc");
    int                  LUMPED_MASS_MATRIX        = args.scalar<int>("LUMPED_MASS_MATRIX");
    int                  MONOLITHIC                = args.scalar<int>("MONOLITHIC");
    double               Rpos[numDOFs], Rneg[numDOFs];
    int                  ij = 0;

    //////////////////////////////////////////////////////
    // ********** COMPUTE LOW ORDER SOLUTION ********** //
    //////////////////////////////////////////////////////
    if (num_fct_iter == 0) { // No FCT for global bounds
      for (int i = 0; i < numDOFs; i++) { solLim.data()[i] = uLow.data()[i]; }
    } else // do FCT iterations (with global bounds) on low order solution
    {
      for (int iter = 0; iter < num_fct_iter; iter++) {
        ij = 0;
        for (int i = 0; i < numDOFs; i++) {
          double maxi = 1.0, Pposi = 0;
          for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
            int j = csrColumnOffsets_DofLoops.data()[offset];
            // compute Flux correction
            double Fluxij = FluxMatrix.data()[ij] - limitedFlux.data()[ij];
            Pposi += Fluxij * ((Fluxij > 0) ? 1. : 0.);
            // update ij
            ij += 1;
          }
          // compute Q vectors
          double mi      = ML.data()[i];
          double solLimi = solLim.data()[i];
          double Qposi   = mi * (maxi - solLimi);
          // compute R vectors
          Rpos[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
        }
        ij = 0;
        for (int i = 0; i < numDOFs; i++) {
          double ith_Limiter_times_FluxCorrectionMatrix = 0.;
          double Rposi                                  = Rpos[i];
          for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
            int j = csrColumnOffsets_DofLoops.data()[offset];
            // Flux Correction
            double Fluxij = FluxMatrix.data()[ij] - limitedFlux.data()[ij];
            // compute limiter
            double Lij = 1.0;
            Lij        = (Fluxij > 0 ? Rposi : Rpos[j]);
            // compute limited flux
            ith_Limiter_times_FluxCorrectionMatrix += Lij * Fluxij;

            // update limited flux
            limitedFlux.data()[ij] = Lij * Fluxij;

            //update FluxMatrix
            FluxMatrix.data()[ij] = Fluxij;

            //update ij
            ij += 1;
          }
          //update limited solution
          double mi = ML.data()[i];
        }
      }
    }

    // ***************************************** //
    // ********** HIGH ORDER SOLUTION ********** //
    // ***************************************** //
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      double mini = soln.data()[i], maxi = soln.data()[i];
      double Pposi = 0, Pnegi = 0.;
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops.data()[offset];
        // compute local bounds //
        mini = fmin(mini, soln.data()[j]);
        maxi = fmax(maxi, soln.data()[j]);
        // compute P vectors //
        double fij = (MC.data()[ij] * (uDotLow.data()[i] - uDotLow.data()[j]) / dt + dLow.data()[ij] * (uLow.data()[i] - uLow.data()[j]));
        Pposi += fij * (fij > 0 ? 1. : 0.);
        Pnegi += fij * (fij < 0 ? 1. : 0.);
        //update ij
        ij += 1;
      }
      // compute Q vectors //
      double mi    = ML.data()[i];
      double Qposi = mi * (maxi - solLim.data()[i]);
      double Qnegi = mi * (mini - solLim.data()[i]);
      // compute R vectors //
      Rpos[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
      Rneg[i] = ((Pnegi == 0) ? 1. : fmin(1.0, Qnegi / Pnegi));
    }

    // COMPUTE LIMITERS //
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      double ith_limited_flux_correction = 0;
      double Rposi                       = Rpos[i];
      double Rnegi                       = Rneg[i];
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops.data()[offset];
        // compute flux correction
        double fij = (MC.data()[ij] * (uDotLow.data()[i] - uDotLow.data()[j]) / dt + dLow.data()[ij] * (uLow.data()[i] - uLow.data()[j]));

        // compute limiters
        double Lij = 1.0;
        Lij        = fij > 0 ? fmin(Rposi, Rneg[j]) : fmin(Rnegi, Rpos[j]);
        // compute ith_limited_flux_correction
        ith_limited_flux_correction += Lij * fij;
        ij += 1;
      }
      double mi = ML.data()[i];
      solLim[i] += 1. / mi * ith_limited_flux_correction;
    }
  }

  void calculateResidual_entropy_viscosity(arguments_dict &args)
  {
    xt::pyarray<double> &globalJacobian            = args.array<double>("globalJacobian");
    double               Theta                     = args.scalar<double>("Theta");
    xt::pyarray<double> &bc_mask                   = args.array<double>("bc_mask");
    double               dt                        = args.scalar<double>("dt");
    xt::pyarray<double> &mesh_trial_ref            = args.array<double>("mesh_trial_ref");
    xt::pyarray<double> &mesh_grad_trial_ref       = args.array<double>("mesh_grad_trial_ref");
    xt::pyarray<double> &mesh_dof                  = args.array<double>("mesh_dof");
    xt::pyarray<double> &mesh_velocity_dof         = args.array<double>("mesh_velocity_dof");
    double               MOVING_DOMAIN             = args.scalar<double>("MOVING_DOMAIN");
    xt::pyarray<int>    &mesh_l2g                  = args.array<int>("mesh_l2g");
    xt::pyarray<double> &dV_ref                    = args.array<double>("dV_ref");
    xt::pyarray<double> &u_trial_ref               = args.array<double>("u_trial_ref");
    xt::pyarray<double> &u_grad_trial_ref          = args.array<double>("u_grad_trial_ref");
    xt::pyarray<double> &u_test_ref                = args.array<double>("u_test_ref");
    xt::pyarray<double> &u_grad_test_ref           = args.array<double>("u_grad_test_ref");
    xt::pyarray<double> &mesh_trial_trace_ref      = args.array<double>("mesh_trial_trace_ref");
    xt::pyarray<double> &mesh_grad_trial_trace_ref = args.array<double>("mesh_grad_trial_trace_ref");
    xt::pyarray<double> &dS_ref                    = args.array<double>("dS_ref");
    xt::pyarray<double> &u_trial_trace_ref         = args.array<double>("u_trial_trace_ref");

    xt::pyarray<double> &u_grad_trial_trace_ref                     = args.array<double>("u_grad_trial_trace_ref");
    xt::pyarray<double> &u_test_trace_ref                           = args.array<double>("u_test_trace_ref");
    xt::pyarray<double> &u_grad_test_trace_ref                      = args.array<double>("u_grad_test_trace_ref");
    xt::pyarray<double> &normal_ref                                 = args.array<double>("normal_ref");
    xt::pyarray<double> &boundaryJac_ref                            = args.array<double>("boundaryJac_ref");
    int                  nElements_global                           = args.scalar<int>("nElements_global");
    xt::pyarray<double> &ebqe_penalty_ext                           = args.array<double>("ebqe_penalty_ext");
    xt::pyarray<int>    &elementMaterialTypes                       = args.array<int>("elementMaterialTypes");
    xt::pyarray<int>    &isSeepageFace                              = args.array<int>("isSeepageFace");
    xt::pyarray<int>    &a_rowptr                                   = args.array<int>("a_rowptr");
    xt::pyarray<int>    &a_colind                                   = args.array<int>("a_colind");
    double               rho                                        = args.scalar<double>("rho");
    double               beta                                       = args.scalar<double>("beta");
    xt::pyarray<double> &gravity                                    = args.array<double>("gravity");
    xt::pyarray<double> &alpha                                      = args.array<double>("alpha");
    xt::pyarray<double> &n                                          = args.array<double>("n");
    xt::pyarray<double> &thetaR                                     = args.array<double>("thetaR");
    xt::pyarray<double> &thetaSR                                    = args.array<double>("thetaSR");
    xt::pyarray<double> &KWs                                        = args.array<double>("KWs");
    double               useMetrics                                 = args.scalar<double>("useMetrics");
    double               alphaBDF                                   = args.scalar<double>("alphaBDF");
    int                  lag_shockCapturing                         = args.scalar<int>("lag_shockCapturing");
    double               shockCapturingDiffusion                    = args.scalar<double>("shockCapturingDiffusion");
    double               sc_uref                                    = args.scalar<double>("sc_uref");
    double               sc_alpha                                   = args.scalar<double>("sc_alpha");
    xt::pyarray<int>    &u_l2g                                      = args.array<int>("u_l2g");
    xt::pyarray<int>    &r_l2g                                      = args.array<int>("r_l2g");
    xt::pyarray<double> &elementDiameter                            = args.array<double>("elementDiameter");
    int                  degree_polynomial                          = args.scalar<int>("degree_polynomial");
    xt::pyarray<double> &u_dof                                      = args.array<double>("u_dof");
    xt::pyarray<double> &u_dof_old                                  = args.array<double>("u_dof_old");
    xt::pyarray<double> &velocity                                   = args.array<double>("velocity");
    xt::pyarray<double> &q_m                                        = args.array<double>("q_m");
    xt::pyarray<double> &q_u                                        = args.array<double>("q_u");
    xt::pyarray<double> &q_dV                                       = args.array<double>("q_dV");
    xt::pyarray<double> &q_m_betaBDF                                = args.array<double>("q_m_betaBDF");
    xt::pyarray<double> &cfl                                        = args.array<double>("cfl");
    xt::pyarray<double> &q_numDiff_u                                = args.array<double>("q_numDiff_u");
    xt::pyarray<double> &q_numDiff_u_last                           = args.array<double>("q_numDiff_u_last");
    int                  offset_u                                   = args.scalar<int>("offset_u");
    int                  stride_u                                   = args.scalar<int>("stride_u");
    xt::pyarray<double> &globalResidual                             = args.array<double>("globalResidual");
    int                  nExteriorElementBoundaries_global          = args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int>    &exteriorElementBoundariesArray             = args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int>    &elementBoundaryElementsArray               = args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int>    &elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
    xt::pyarray<double> &ebqe_velocity_ext                          = args.array<double>("ebqe_velocity_ext");
    xt::pyarray<int>    &isDOFBoundary_u                            = args.array<int>("isDOFBoundary_u");
    xt::pyarray<double> &ebqe_bc_u_ext                              = args.array<double>("ebqe_bc_u_ext");
    xt::pyarray<int>    &isFluxBoundary_u                           = args.array<int>("isFluxBoundary_u");
    xt::pyarray<double> &ebqe_bc_flux_ext                           = args.array<double>("ebqe_bc_flux_ext");
    xt::pyarray<double> &ebqe_phi                                   = args.array<double>("ebqe_phi");
    double               epsFact                                    = args.scalar<double>("epsFact");
    xt::pyarray<double> &ebqe_u                                     = args.array<double>("ebqe_u");
    xt::pyarray<double> &ebqe_flux                                  = args.array<double>("ebqe_flux");
    // PARAMETERS FOR EDGE BASED STABILIZATION
    double cE = args.scalar<double>("cE");
    double cK = args.scalar<double>("cK");
    // PARAMETERS FOR LOG BASED ENTROPY FUNCTION
    double uL = args.scalar<double>("uL");
    double uR = args.scalar<double>("uR");
    // PARAMETERS FOR EDGE VISCOSITY
    int               numDOFs                       = args.scalar<int>("numDOFs");
    int               NNZ                           = args.scalar<int>("NNZ");
    xt::pyarray<int> &csrRowIndeces_DofLoops        = args.array<int>("csrRowIndeces_DofLoops");
    xt::pyarray<int> &csrColumnOffsets_DofLoops     = args.array<int>("csrColumnOffsets_DofLoops");
    xt::pyarray<int> &csrRowIndeces_CellLoops       = args.array<int>("csrRowIndeces_CellLoops");
    xt::pyarray<int> &csrColumnOffsets_CellLoops    = args.array<int>("csrColumnOffsets_CellLoops");
    xt::pyarray<int> &csrColumnOffsets_eb_CellLoops = args.array<int>("csrColumnOffsets_eb_CellLoops");
    // C matrices
    xt::pyarray<double> &Cx  = args.array<double>("Cx");
    xt::pyarray<double> &Cy  = args.array<double>("Cy");
    xt::pyarray<double> &Cz  = args.array<double>("Cz");
    xt::pyarray<double> &CTx = args.array<double>("CTx");
    xt::pyarray<double> &CTy = args.array<double>("CTy");
    xt::pyarray<double> &CTz = args.array<double>("CTz");
    xt::pyarray<double> &ML  = args.array<double>("ML");
    xt::pyarray<double> &MC  = args.array<double>("MC");

    xt::pyarray<double> &delta_x_ij = args.array<double>("delta_x_ij");
    // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
    int LUMPED_MASS_MATRIX = args.scalar<int>("LUMPED_MASS_MATRIX");
    STABILIZATION STABILIZATION_TYPE{static_cast<STABILIZATION>(args.scalar<int>("STABILIZATION_TYPE"))};
    
    
    //////////////////////////////////////For Brooks- Corey/////////////////
    PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
    double BC_entry_head = args.scalar<double>("BC_entry_head");
    double BC_lambda     = args.scalar<double>("BC_lambda");


    int ENTROPY_TYPE = args.scalar<int>("ENTROPY_TYPE");
    // FOR FCT
    xt::pyarray<double> &dLow                 = args.array<double>("dLow");
    xt::pyarray<double> &fluxMatrix           = args.array<double>("fluxMatrix");
    xt::pyarray<double> &mDotLow              = args.array<double>("mDotLow");
    xt::pyarray<double> &mDotHigh              = args.array<double>("mDotHigh");
    xt::pyarray<double> &mLow                 = args.array<double>("mLow");
    xt::pyarray<double> &dt_times_fH_minus_fL = args.array<double>("dt_times_fH_minus_fL");
    xt::pyarray<double> &min_m_bc             = args.array<double>("min_m_bc");
    xt::pyarray<double> &max_m_bc             = args.array<double>("max_m_bc");
    // AUX QUANTITIES OF INTEREST
    xt::pyarray<double> &quantDOFs = args.array<double>("quantDOFs");
    xt::pyarray<double> &mn        = args.array<double>("mn");
    xt::pyarray<double> &fluxCorrection        = args.array<double>("fluxCorrection");
    xt::pyarray<double> &limited_solution          = args.array<double>("limited_solution");

    xt::pyarray<double> &anb_seepage_flux_n = args.array<double>("anb_seepage_flux_n");
    xt::pyarray<double> &q_velocity = args.array<double>("q_velocity");
    double &anb_seepage_flux(args.scalar<double>("anb_seepage_flux"));
    anb_seepage_flux = 0.0;
    xt::pyarray<int>    &csrRowIndeces_u_u                          = args.array<int>("csrRowIndeces_u_u");
    xt::pyarray<int>    &csrColumnOffsets_u_u                       = args.array<int>("csrColumnOffsets_u_u");
    xt::pyarray<int>    &csrColumnOffsets_eb_u_u                    = args.array<int>("csrColumnOffsets_eb_u_u");
    
    double Rpos[numDOFs], Rneg[numDOFs];
    //double FluxCorrectionMatrix[NNZ];
    // NOTE: This function follows a different (but equivalent) implementation of the smoothness based indicator than NCLS.h
    // Allocate space for the transport matrices
    // This is used for first order KUZMIN'S METHOD
    double                TransportMatrix[NNZ], TransportMatrixConsistent[NNZ];
    double                TransportMatrixn[NNZ], TransportMatrixConsistentn[NNZ];
    std::valarray<double> u_free_dof(numDOFs);
    std::valarray<double> u_free_dof_old(numDOFs);
    std::valarray<double> ML2(numDOFs);

    for (int eN = 0; eN < nElements_global; eN++)
      for (int j = 0; j < nDOF_trial_element; j++) {
        int eN_nDOF_trial_element                               = eN * nDOF_trial_element;
        u_free_dof[r_l2g.data()[eN_nDOF_trial_element + j]]     = u_dof.data()[u_l2g.data()[eN_nDOF_trial_element + j]];
        u_free_dof_old[r_l2g.data()[eN_nDOF_trial_element + j]] = u_dof_old.data()[u_l2g.data()[eN_nDOF_trial_element + j]];
      }
    for (int i = 0; i < NNZ; i++) {
      TransportMatrix[i]            = 0.;
      TransportMatrixConsistent[i]  = 0.;
      TransportMatrixn[i]           = 0.;
      TransportMatrixConsistentn[i] = 0.;
    }

    // compute entropy and init global_entropy_residual and boundary_integral
    double psi[numDOFs], eta[numDOFs], global_entropy_residual[numDOFs], boundary_integral[numDOFs];
    for (int i = 0; i < numDOFs; i++) {
      // NODAL ENTROPY //
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV stab
      {
        double solni = 1.0 * u_free_dof_old[i];
        eta[i]                      = ENTROPY_TYPE == 1 ? ENTROPY(solni, uL, uR) : ENTROPY_LOG(solni, uL, uR);
        global_entropy_residual[i]  = 0.;
      }
      boundary_integral[i] = 0.;
      ML2[i]               = 0.0;
    }

    //////////////////////////////////////////////
    // ** LOOP IN CELLS FOR CELL BASED TERMS ** //
    //////////////////////////////////////////////
    // HERE WE COMPUTE:
    //    * Time derivative term. u_t
    //    * cell based CFL (for reference)
    //    * Entropy residual
    //    * Transport matrices
    for (int eN = 0; eN < nElements_global; eN++) {
      //declare local storage for local contributions and initialize
      double elementResidual_u[nDOF_test_element], element_entropy_residual[nDOF_test_element], Phi[nDOF_trial_element], Phi_n[nDOF_trial_element];
      double elementTransport[nDOF_test_element][nDOF_trial_element], elementTransportConsistent[nDOF_test_element][nDOF_trial_element];
      double elementTransportn[nDOF_test_element][nDOF_trial_element], elementTransportConsistentn[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++) {
        Phi[i]   = u_dof[i];
        Phi_n[i] = u_dof_old[i];
        for (int I = 0; I < nSpace; I++) {
          Phi[i] -= rho * mesh_dof[i * 3 + I] * gravity[I];
          Phi_n[i] -= rho * mesh_dof[i * 3 + I] * gravity[I];
        }
        elementResidual_u[i]        = 0.0;
        element_entropy_residual[i] = 0.0;
        for (int j = 0; j < nDOF_trial_element; j++) {
          elementTransport[i][j]            = 0.0;
          elementTransportConsistent[i][j]  = 0.0;
          elementTransportn[i][j]           = 0.0;
          elementTransportConsistentn[i][j] = 0.0;
        }
      }
      //loop over quadrature points and compute integrands
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        //compute indeces and declare local storage
        int eN_k = eN * nQuadraturePoints_element + k, eN_k_nSpace = eN_k * nSpace, eN_nDOF_trial_element = eN * nDOF_trial_element;
        double
          // for entropy residual
          aux_entropy_residual = 0.,
          DENTROPY_un, DENTROPY_uni,
          //for mass matrix contributions
          u = 0.0, un = 0.0, grad_u[nSpace], grad_un[nSpace], velocity_loc[nSpace], u_test_dV[nDOF_trial_element], u_grad_trial[nDOF_trial_element * nSpace], u_grad_test_dV[nDOF_test_element * nSpace],
          //for general use
          jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace], dV, x, y, z, xt, yt, zt, m, dm, f[nSpace], df[nSpace], a[nnz], da[nnz], as[nnz], mn, dmn, fn[nSpace], dfn[nSpace], an[nnz], dan[nnz], asn[nnz];
        //get the physical integration weight
        ck.calculateMapping_element(eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y, z);
        ck.calculateMappingVelocity_element(eN, k, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), xt, yt, zt);
        dV = fabs(jacDet) * dV_ref.data()[k];
        //get the solution (of Newton's solver). To compute time derivative term
        ck.valFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], u);
        //get the solution at quad point at tn and tnm1 for entropy viscosity
        ck.valFromDOF(u_dof_old.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], un);
        //get the solution gradients at tn for entropy viscosity
        ck.gradTrialFromRef(&u_grad_trial_ref.data()[k * nDOF_trial_element * nSpace], jacInv, u_grad_trial);
        //precalculate test function products with integration weights for mass matrix terms
        for (int I = 0; I < nSpace; I++) {
          grad_u[I]  = 0.0;
          grad_un[I] = 0.0;
        }

        for (int j = 0; j < nDOF_trial_element; j++) {
          u_test_dV[j] = u_test_ref.data()[k * nDOF_trial_element + j] * dV;
          for (int I = 0; I < nSpace; I++) {
            grad_un[I] += Phi_n[j] * u_grad_trial[j * nSpace + I];//note: grad u is grad phi
            grad_u[I] += Phi[j] * u_grad_trial[j * nSpace + I];
            u_grad_test_dV[j * nSpace + I] = u_grad_trial[j * nSpace + I] * dV; //cek warning won't work for Petrov-Galerkin
          }
        }
        //
        //calculate pde coefficients at quadrature points
        //
        double Kr, dKr, Krn, dKrn;
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes[eN]], n.data()[elementMaterialTypes[eN]], thetaR.data()[elementMaterialTypes[eN]], thetaSR.data()[elementMaterialTypes[eN]],
                             &KWs.data()[elementMaterialTypes[eN] * nnz], un, mn, dmn, fn, dfn, an, dan, asn, Krn, dKrn,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK)););
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes[eN]], n.data()[elementMaterialTypes[eN]], thetaR.data()[elementMaterialTypes[eN]], thetaSR.data()[elementMaterialTypes[eN]],
                             &KWs.data()[elementMaterialTypes[eN] * nnz], u, m, dm, f, df, a, da, as, Kr, dKr,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);        // lambda  (only used if BC_PSK)););

        // Darcy velocity calculation
        for (int I = 0; I < nSpace; I++) { velocity_loc[I] = 0.0; }
        for (int I = 0; I < nSpace; I++) {
          for (int J = 0; J < nSpace; J++) { velocity_loc[I] -= Kr * KWs.data()[elementMaterialTypes[eN] * nSpace * nSpace + I * nSpace + J] * grad_u[J]; }
        }
        for (int I = 0; I < nSpace; I++) { q_velocity.data()[eN_k_nSpace + I] = velocity_loc[I]; }

        //
        //moving mesh
        //
        double mesh_velocity[3];
        mesh_velocity[0] = xt;
        mesh_velocity[1] = yt;
        mesh_velocity[2] = zt;
        //relative velocity at tn
        for (int I = 0; I < nSpace; I++) {
          f[I] -= MOVING_DOMAIN * m * mesh_velocity[I];
          velocity_loc[I] = df[I] * (2.0 * dm * dm / (dm * dm + fmax(1.0e-16, dm * dm)));
        }
        //////////////////////////////
        // CALCULATE CELL BASED CFL //
        //////////////////////////////
        calculateCFL(elementDiameter.data()[eN] / degree_polynomial, velocity_loc, cfl.data()[eN_k]);

        //////////////////////////////////////////////
        // CALCULATE ENTROPY RESIDUAL AT QUAD POINT //
        //////////////////////////////////////////////
        if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) // EV stab
        {
          for (int I = 0; I < nSpace; I++) aux_entropy_residual += velocity_loc[I] * grad_un[I];
          DENTROPY_un = ENTROPY_TYPE == 1 ? DENTROPY(un, uL, uR) : DENTROPY_LOG(un, uL, uR);
        }
        //////////////
        // ith-LOOP //
        //////////////
        for (int i = 0; i < nDOF_test_element; i++) {
          // VECTOR OF ENTROPY RESIDUAL //
          int eN_i = eN * nDOF_test_element + i;
          ML2[u_l2g.data()[eN_i]] += u_test_dV[i];
          if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) // EV stab
          {
            int    gi                 = offset_u + stride_u * u_l2g.data()[eN_i]; //global i-th index
            double uni = u_dof_old.data()[gi];
            DENTROPY_uni              = ENTROPY_TYPE == 1 ? DENTROPY(uni, uL, uR) : DENTROPY_LOG(uni, uL, uR);
            element_entropy_residual[i] += (DENTROPY_un - DENTROPY_uni) * aux_entropy_residual * u_test_dV[i];
          }
          elementResidual_u[i] += m * u_test_dV[i];
          ///////////////
          // j-th LOOP // To construct transport matrices
          ///////////////
          for (int j = 0; j < nDOF_trial_element; j++) {
            int j_nSpace = j * nSpace;
            int i_nSpace = i * nSpace;
            elementTransport[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), as, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportConsistent[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), a, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportn[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), asn, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportConsistentn[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), an, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
          }
        } //i
        //save solution for other models
        q_u.data()[eN_k] = u;
        q_m.data()[eN_k] = m;
      }
      /////////////////
      // DISTRIBUTE // load cell based element into global residual
      ////////////////
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        int gi   = offset_u + stride_u * r_l2g.data()[eN_i]; //global i-th index
        // distribute entropy_residual
        if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) // EV Stab
          global_entropy_residual[gi] += element_entropy_residual[i];
        // distribute transport matrices
        for (int j = 0; j < nDOF_trial_element; j++) {
          int eN_i_j = eN_i * nDOF_trial_element + j;
          TransportMatrix[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransport[i][j];
          TransportMatrixConsistent[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportConsistent[i][j];
          TransportMatrixn[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportn[i][j];
          TransportMatrixConsistentn[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportConsistentn[i][j];
        } //j
      } //i
    } //elementsxw

        //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
    //
    //ebNE is the Exterior element boundary INdex
    //ebN is the element boundary INdex
    //eN is the element index
    for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) {
      int    ebN = exteriorElementBoundariesArray.data()[ebNE], eN = elementBoundaryElementsArray.data()[ebN * 2 + 0], ebN_local = elementBoundaryLocalElementBoundariesArray.data()[ebN * 2 + 0], eN_nDOF_trial_element = eN * nDOF_trial_element;
      double elementResidual_u[nDOF_test_element];
      for (int i = 0; i < nDOF_test_element; i++) { elementResidual_u[i] = 0.0; }
      for (int kb = 0; kb < nQuadraturePoints_elementBoundary; kb++) {
        int    ebNE_kb = ebNE * nQuadraturePoints_elementBoundary + kb, ebNE_kb_nSpace = ebNE_kb * nSpace, ebN_local_kb = ebN_local * nQuadraturePoints_elementBoundary + kb, ebN_local_kb_nSpace = ebN_local_kb * nSpace;
        double u_ext = 0.0, un_ext, grad_u_ext[nSpace], m_ext = 0.0, dm_ext = 0.0, f_ext[nSpace], df_ext[nSpace], a_ext[nnz], da_ext[nnz], as_ext[nnz], 
        mn_ext = 0.0, dmn_ext = 0.0, fn_ext[nSpace], dfn_ext[nSpace], an_ext[nnz], dan_ext[nnz], asn_ext[nnz], flux_ext = 0.0, bflux_ext = 0.0,
               //anb_seepage_flux=0.0, // for flux calculation
          bc_u_ext = 0.0, bc_grad_u_ext[nSpace], bc_m_ext = 0.0, bc_dm_ext = 0.0, bc_f_ext[nSpace], bc_df_ext[nSpace], bc_a_ext[nnz], bc_da_ext[nnz], bc_as_ext[nnz], jac_ext[nSpace * nSpace], jacDet_ext, jacInv_ext[nSpace * nSpace], boundaryJac[nSpace * (nSpace - 1)], metricTensor[(nSpace - 1) * (nSpace - 1)], metricTensorDetSqrt, dS, u_test_dS[nDOF_test_element], u_grad_trial_trace[nDOF_trial_element * nSpace], normal[3], x_ext, y_ext, z_ext, xt_ext, yt_ext, zt_ext, integralScaling, G[nSpace * nSpace], G_dd_G, tr_G, fluxJacobian_u_u[nDOF_trial_element], bfluxJacobian_u_u[nDOF_trial_element], fluxJacobian_un_un[nDOF_trial_element];
        //
        //calculate the solution and gradients at quadrature points
        //
        //compute information about mapping from reference element to physical element
        ck.calculateMapping_elementBoundary(eN, ebN_local, kb, ebN_local_kb, mesh_dof.data(), mesh_l2g.data(), mesh_trial_trace_ref.data(), mesh_grad_trial_trace_ref.data(), boundaryJac_ref.data(), jac_ext, jacDet_ext, jacInv_ext, boundaryJac, metricTensor, metricTensorDetSqrt,
                                            normal_ref.data(), normal, x_ext, y_ext, z_ext);
        ck.calculateMappingVelocity_elementBoundary(eN, ebN_local, kb, ebN_local_kb, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_trace_ref.data(), xt_ext, yt_ext, zt_ext, normal, boundaryJac, metricTensor, integralScaling);
        dS = ((1.0 - MOVING_DOMAIN) * metricTensorDetSqrt + MOVING_DOMAIN * integralScaling) * dS_ref.data()[kb];
        //get the metric tensor
        //cek todo use symmetry
        ck.calculateG(jacInv_ext, G, G_dd_G, tr_G);
        //compute shape and solution information
        //shape
        ck.gradTrialFromRef(&u_grad_trial_trace_ref.data()[ebN_local_kb_nSpace * nDOF_trial_element], jacInv_ext, u_grad_trial_trace);
        //solution and gradient
        ck.valFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], u_ext);
        ck.valFromDOF(u_dof_old.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], un_ext);
        ck.gradFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], u_grad_trial_trace, grad_u_ext);
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) { u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb * nDOF_test_element + j] * dS; }
        //
        //load the boundary values
        //
        bc_u_ext = isDOFBoundary_u.data()[ebNE_kb] * ebqe_bc_u_ext.data()[ebNE_kb] + (1 - isDOFBoundary_u.data()[ebNE_kb]) * u_ext;
        //
        //calculate the pde coefficients using the solution and the boundary values for the solution
        //
        double bc_Kr, bc_dKr,bc_Kr_ext, bc_dKr_ext, bc_Krn, bc_dKrn;
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u_ext, m_ext, dm_ext, f_ext, df_ext, a_ext, da_ext, as_ext, bc_Kr, bc_dKr,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], un_ext, mn_ext, dmn_ext, fn_ext, dfn_ext, an_ext, dan_ext, asn_ext, bc_Krn, bc_dKrn,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], bc_u_ext, bc_m_ext, bc_dm_ext, bc_f_ext, bc_df_ext, bc_a_ext, bc_da_ext, bc_as_ext, bc_Kr_ext, bc_dKr_ext,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);
        //
        //calculate the numerical fluxes
        //
        bool useConsistentFlux=false;
        if (useConsistentFlux) {
          exteriorNumericalFlux(ebqe_bc_flux_ext[ebNE_kb], a_rowptr.data(), a_colind.data(),
                                isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                                isDOFBoundary_u.data()[ebNE_kb], normal, bc_u_ext, a_ext, grad_u_ext, u_ext, f_ext,
                                ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                                flux_ext);
        } else {
          exteriorNumericalFlux2(ebqe_bc_flux_ext[ebNE_kb], a_rowptr.data(), a_colind.data(),
                              isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                              isDOFBoundary_u.data()[ebNE_kb], normal, bc_u_ext, a_ext, grad_u_ext, u_ext, f_ext,
                              ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                              flux_ext, bflux_ext);
        }
        ebqe_flux.data()[ebNE_kb] = flux_ext;

        anb_seepage_flux             = seepagefluxcalculator(anb_seepage_flux, isSeepageFace.data()[ebNE], dS, flux_ext);
        anb_seepage_flux_n.data()[0] = anb_seepage_flux;
        ebqe_u.data()[ebNE_kb]       = u_ext;
        //
        //update residuals
        //
        for (int i = 0; i < nDOF_test_element; i++) {
          if (useConsistentFlux) {
            elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext, u_test_dS[i]);
          } else {
            elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(bflux_ext, u_test_dS[i]);
          }
        } //i
        for (int j = 0; j < nDOF_trial_element; j++) {
          if (useConsistentFlux) {
          exteriorNumericalFluxJacobian(a_rowptr.data(), a_colind.data(), isDOFBoundary_u.data()[ebNE_kb], normal, a_ext, da_ext, grad_u_ext, &u_grad_trial_trace[j * nSpace], df_ext, u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u[j]);
          } else {
            exteriorNumericalFluxJacobian2(a_rowptr.data(), a_colind.data(), isDOFBoundary_u.data()[ebNE_kb], normal, as_ext, a_ext, da_ext, grad_u_ext, &u_grad_trial_trace[j * nSpace], df_ext, u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u[j],bfluxJacobian_u_u[j]);
          }
          //probably need isDOFBoundary_un here
          //exteriorNumericalFluxJacobian(a_rowptr.data(), a_colind.data(), isDOFBoundary_u.data()[ebNE_kb], normal, asn_ext, dan_ext, grad_u_ext, &u_grad_trial_trace[j * nSpace], dfn_ext, u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
          //                              ebqe_penalty_ext.data()[ebNE_kb], //penalty,
          //                              fluxJacobian_un_un[j]);
        } //j
        //
        //update the element and global residual storage
        //
        for (int i = 0; i < nDOF_test_element; i++) {
          int eN_i = eN * nDOF_test_element + i;
          for (int j = 0; j < nDOF_trial_element; j++) {
            int ebN_i_j = ebN * 4 * nDOF_test_X_trial_element + i * nDOF_trial_element + j;
            if (useConsistentFlux) {
              globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j] * u_test_dS[i];
            } else {
              globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += bfluxJacobian_u_u[j] * u_test_dS[i];
              TransportMatrix[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j] * u_test_dS[i];
              TransportMatrixConsistent[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j] * u_test_dS[i];
              TransportMatrixn[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_un_un[j] * u_test_dS[i];
              TransportMatrixConsistentn[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_un_un[j] * u_test_dS[i];
            }
          } //j
        } //i
      } //kb
      for (int i = 0; i < nDOF_test_element; i++) {
          int eN_i = eN * nDOF_test_element + i;
          globalResidual.data()[offset_u + stride_u * u_l2g.data()[eN_i]] += elementResidual_u[i];
      }//i
    } //ebNE
    /////////////////////////////////////////////////////////////////
    // COMPUTE SMOOTHNESS INDICATOR and NORMALIZE ENTROPY RESIDUAL //
    /////////////////////////////////////////////////////////////////
    // NOTE: see NCLS.h for a different but equivalent implementation of this.
    int ij = 0;
    double cflux[numDOFs];
    for (int i = 0; i < numDOFs; i++) {
      double gi[nSpace], Cij[nSpace], xi[nSpace], etaMaxi, etaMini;
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stabilization
      {
        // For eta min and max
        etaMaxi = fabs(eta[i]);
        etaMini = fabs(eta[i]);
      }
      double solni = u_free_dof_old[i];
      // initialize gi and compute xi
      for (int I = 0; I < nSpace; I++) {
        gi[I] = 0.;
        xi[I] = mesh_dof.data()[i * 3 + I];
      }
      // for smoothness indicator //
      double alpha_numerator_pos = 0., alpha_numerator_neg = 0., alpha_denominator_pos = 0., alpha_denominator_neg = 0.;
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) { // First loop in j (sparsity pattern)
        int j = csrColumnOffsets_DofLoops.data()[offset];
        if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stabilization
        {
          // COMPUTE ETA MIN AND ETA MAX //
          etaMaxi = fmax(etaMaxi, fabs(eta[j]));
          etaMini = fmin(etaMini, fabs(eta[j]));
        }
        double solnj = u_free_dof_old[j];
        // Update Cij matrices
        Cij[0] = Cx[ij];
#if nSpace == 2
        Cij[1] = Cy[ij];
#endif
#if nSpace == 3
        Cij[2] = Cz[ij];
#endif
        // COMPUTE gi VECTOR. gi=1/mi*sum_j(Cij*solj)
        for (int I = 0; I < nSpace; I++) gi[I] += Cij[I] * solnj;

        // COMPUTE numerator and denominator of smoothness indicator
        double alpha_num = solni - solnj;
        if (alpha_num >= 0.) {
          alpha_numerator_pos += alpha_num;
          alpha_denominator_pos += alpha_num;
        } else {
          alpha_numerator_neg += alpha_num;
          alpha_denominator_neg += fabs(alpha_num);
        }
        //update ij
        ij += 1;
      }
      // scale g vector by lumped mass matrix
      //double mass_matrix_error = abs(ML.data()[i] - ML2[i]);
      //if (mass_matrix_error > 1.0e-16) std::cout << mass_matrix_error<<" ML " << ML.data()[i] << '\t' << ML2[i] << std::endl;
      for (int I = 0; I < nSpace; I++) gi[I] /= ML.data()[i];
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stab
      {
        // Normalizae entropy residual
        global_entropy_residual[i] *= etaMini == etaMaxi ? 0. : 2 * cE / (etaMaxi - etaMini);
        quantDOFs.data()[i] = fabs(global_entropy_residual[i]);
      }

      // Now that I have the gi vectors, I can use them for the current i-th DOF
      double SumPos = 0., SumNeg = 0.;
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) { // second loop in j (sparsity pattern)
        int j = csrColumnOffsets_DofLoops.data()[offset];
        // compute xj
        double xj[nSpace];
        for (int I = 0; I < nSpace; I++) xj[I] = mesh_dof.data()[j * 3 + I];
        // compute gi*(xi-xj)
        double gi_times_x = 0.;
        for (int I = 0; I < nSpace; I++) {
          gi_times_x += gi[I] * delta_x_ij.data()[offset * 3 + I];
        }
        // compute the positive and negative part of gi*(xi-xj)
        SumPos += gi_times_x > 0 ? gi_times_x : 0;
        SumNeg += gi_times_x < 0 ? gi_times_x : 0;
      }
      double sigmaPosi  = fmin(1., (fabs(SumNeg) + 1E-15) / (SumPos + 1E-15));
      double sigmaNegi  = fmin(1., (SumPos + 1E-15) / (fabs(SumNeg) + 1E-15));
      double alpha_numi = fabs(sigmaPosi * alpha_numerator_pos + sigmaNegi * alpha_numerator_neg);
      double alpha_deni = sigmaPosi * alpha_denominator_pos + sigmaNegi * alpha_denominator_neg;
      if (IS_BETAij_ONE == 1) {
        alpha_numi = fabs(alpha_numerator_pos + alpha_numerator_neg);
        alpha_deni = alpha_denominator_pos + alpha_denominator_neg;
      }
      double alphai       = alpha_numi / (alpha_deni + 1E-15);
      quantDOFs.data()[i] = alphai;

      if (POWER_SMOOTHNESS_INDICATOR == 0) psi[i] = 1.0;
      else psi[i] = std::pow(alphai, POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
    }
    /////////////////////////////////////////////
    // ** LOOP IN DOFs FOR EDGE BASED TERMS ** //
    /////////////////////////////////////////////
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      int    ii;
      double sum_abs_dt_times_fH_minus_fL = 0.0, phi_i  = u_free_dof[i], phin_i = u_free_dof_old[i], MLi = ML.data()[i];
      double Kr, dKr, Krn, dKrn;
      double J_ii = 0.0;
      double ith_dissipative_term           = 0;
      double ith_low_order_dissipative_term = 0;
      double ith_flux_term                  = 0;
      double ith_consistent_flux_term       = 0;
      double dLii                           = 0.;
      double m, dm, f[nSpace], df[nSpace], a[nnz], da[nnz], as[nnz];
      double dmn, fn[nSpace], dfn[nSpace], an[nnz], dan[nnz], asn[nnz];

      for (int I = 0; I < nSpace; I++) {
        phi_i -= rho * gravity.data()[I] * mesh_dof.data()[i * 3 + I];
        phin_i -= rho * gravity.data()[I] * mesh_dof.data()[i * 3 + I];
      }
      // loop over the sparsity pattern of the i-th DOF
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops.data()[offset];
        if (i == j) ii = ij;
        double phi_j  = u_free_dof[j], phin_j = u_free_dof_old[j];

        for (int I = 0; I < nSpace; I++) {
          phi_j -= rho * gravity.data()[I] * mesh_dof.data()[j * 3 + I];
          phin_j -= rho * gravity.data()[I] * mesh_dof.data()[j * 3 + I];
        }

        double dLowij, dLij, dEVij, dHij, fH, fL, fA=0.0;
        fH = -Theta * TransportMatrixConsistent[ij] * (phi_j - phi_i) - (1 - Theta) * TransportMatrixConsistentn[ij] * (phin_j - phin_i);
        ith_consistent_flux_term += fH;
        fA = fH;
        if (-TransportMatrix[ij] * (phi_j - phi_i) <= 0.0) {
          evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(),
                               alpha.data()[elementMaterialTypes.data()[0]], //cek hack, only for 1 material
                               n.data()[elementMaterialTypes.data()[0]], thetaR.data()[elementMaterialTypes.data()[0]], thetaSR.data()[elementMaterialTypes.data()[0]], &KWs.data()[elementMaterialTypes.data()[0] * nnz], u_free_dof[i], m, dm, f, df, a, da, as, Kr, dKr,
                               PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                               BC_entry_head,     // h_e (only used if BC_PSK)
                               BC_lambda);
          fL = Theta * Kr * fmax(0.0, -TransportMatrix[ij]) * (phi_j - phi_i);
          if (i != j) {
            globalJacobian.data()[ij] -= Theta * Kr * fmax(0.0, -TransportMatrix[ij]);
            J_ii -= -Theta * Kr * fmax(0.0, -TransportMatrix[ij]) + Theta * dKr * fmax(0.0, -TransportMatrix[ij]) * (phi_j - phi_i);
          }
          ith_flux_term += fL;
          fA -= fL;
        } else {
          evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(),
                               alpha.data()[elementMaterialTypes.data()[0]], //cek hack, only for 1 material
                               n.data()[elementMaterialTypes.data()[0]], thetaR.data()[elementMaterialTypes.data()[0]], thetaSR.data()[elementMaterialTypes.data()[0]], &KWs.data()[elementMaterialTypes.data()[0] * nnz], u_free_dof[j], m, dm, f, df, a, da, as, Kr, dKr,
                               PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                               BC_entry_head,     // h_e (only used if BC_PSK)
                               BC_lambda);
          fL = Theta * Kr * fmax(0.0, -TransportMatrix[ij]) * (phi_j - phi_i);
          if (i != j) {
            globalJacobian.data()[ij] -= Theta * Kr * fmax(0.0, -TransportMatrix[ij]) + Theta * dKr * fmax(0.0, -TransportMatrix[ij]) * (phi_j - phi_i);
            J_ii -= -Theta * Kr * fmax(0.0, -TransportMatrix[ij]);
          }
          ith_flux_term += fL;
          fA -= fL;
        }
        if (-TransportMatrixn[ij] * (phin_j - phin_i) <= 0.0) {
          evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(),
                               alpha.data()[elementMaterialTypes.data()[0]], //cek hack, only for 1 material
                               n.data()[elementMaterialTypes.data()[0]], thetaR.data()[elementMaterialTypes.data()[0]], thetaSR.data()[elementMaterialTypes.data()[0]], &KWs.data()[elementMaterialTypes.data()[0] * nnz], u_free_dof_old[i], m, dm, f, df, a, da, as, Kr, dKr,
                               PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                               BC_entry_head,     // h_e (only used if BC_PSK)
                               BC_lambda);
          fL = (1 - Theta) * Kr * fmax(0.0, -TransportMatrixn[ij]) * (phin_j - phin_i);
          ith_flux_term += fL;
          fA -= fL;
        } else {
          evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(),
                               alpha.data()[elementMaterialTypes.data()[0]], //cek hack, only for 1 material
                               n.data()[elementMaterialTypes.data()[0]], thetaR.data()[elementMaterialTypes.data()[0]], thetaSR.data()[elementMaterialTypes.data()[0]], &KWs.data()[elementMaterialTypes.data()[0] * nnz], u_free_dof_old[j], m, dm, f, df, a, da, as, Kr, dKr,
                               PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                               BC_entry_head,     // h_e (only used if BC_PSK)
                               BC_lambda);
          fL = (1 - Theta) * Kr * fmax(0.0, -TransportMatrixn[ij]) * (phin_j - phin_i);
          ith_flux_term += fL;
          fA -= fL;
        }
        dt_times_fH_minus_fL.data()[ij] = dt * fA;
        ij += 1;
      }

      mDotLow.data()[i] = ith_flux_term/MLi;
      cflux[i] = ith_consistent_flux_term;
      evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]], //cek hack, only for 1 material
                           n.data()[elementMaterialTypes.data()[0]], thetaR.data()[elementMaterialTypes.data()[0]], thetaSR.data()[elementMaterialTypes.data()[0]], &KWs.data()[elementMaterialTypes.data()[0] * nnz], u_free_dof[i], m, dm, f, df, a, da, as, Kr, dKr,
                           PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                           BC_entry_head,     // h_e (only used if BC_PSK)
                           BC_lambda);
      evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]], //cek hack, only for 1 material
                           n.data()[elementMaterialTypes.data()[0]], thetaR.data()[elementMaterialTypes.data()[0]], thetaSR.data()[elementMaterialTypes.data()[0]], &KWs.data()[elementMaterialTypes.data()[0] * nnz], u_free_dof_old[i], mn.data()[i], dmn, fn, dfn, an, dan, asn, Krn, dKrn,
                           PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                           BC_entry_head,     // h_e (only used if BC_PSK)
                           BC_lambda);
      mLow.data()[i] = m;
      globalResidual.data()[i] += bc_mask.data()[i] * (MLi * (m - mn.data()[i]) / dt - ith_flux_term);
      globalJacobian.data()[ii] += bc_mask.data()[i] * (MLi * dm / dt + J_ii) + (1.0 - bc_mask.data()[i]);
    }
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      mDotHigh[i] = cflux[i];
      for (int offset = csrRowIndeces_DofLoops.data()[i]; offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops.data()[offset];
        mDotHigh[i] -= MC.data()[ij]*cflux[j]/ML.data()[j];
        ij +=1;
      }
      mDotHigh[i] = (cflux[i] + mDotHigh[i])/ML.data()[i];
    }
    if (STABILIZATION_TYPE == STABILIZATION::Implicit_FCT) {
      FCTStep(args);
      for (int i = 0; i < numDOFs; i++) {
        globalResidual.data()[i] += fluxCorrection.data()[i];
      }
    }

  }

  void invert(arguments_dict &args)
  {
    xt::pyarray<int>    &a_rowptr             = args.array<int>("a_rowptr");
    xt::pyarray<int>    &a_colind             = args.array<int>("a_colind");
    double               rho                  = args.scalar<double>("rho");
    double               beta                 = args.scalar<double>("beta");
    xt::pyarray<double> &gravity              = args.array<double>("gravity");
    xt::pyarray<double> &alpha                = args.array<double>("alpha");
    xt::pyarray<double> &n                    = args.array<double>("n");
    xt::pyarray<double> &thetaR               = args.array<double>("thetaR");
    xt::pyarray<double> &thetaSR              = args.array<double>("thetaSR");
    xt::pyarray<double> &KWs                  = args.array<double>("KWs");
    xt::pyarray<int>    &elementMaterialTypes = args.array<int>("elementMaterialTypes");
    int                  numDOFs              = args.scalar<int>("numDOFs");
    xt::pyarray<double> &mIn = args.array<double>("limited_solution");
    xt::pyarray<double> &pOut = args.array<double>("u_dof");

        //////////////////////////////////////For Brooks- Corey/////////////////
    PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
    double BC_entry_head = args.scalar<double>("BC_entry_head");
    double BC_lambda     = args.scalar<double>("BC_lambda");


    for (int i = 0; i < numDOFs; i++) {
      double dm, f[nSpace], df[nSpace], a[nnz], da[nnz];
      double mMin = rho * thetaR.data()[elementMaterialTypes.data()[0]];
      double mMax = rho * (thetaR.data()[elementMaterialTypes.data()[0]] + thetaSR.data()[elementMaterialTypes.data()[0]]);

      if (mIn.data()[i] < mMin - 0.001 || mIn.data()[i] > mMax + 0.001) { std::cout << "mass out of bounds " << mMin << '\t' << mIn.data()[i] << '\t' << mMax << std::endl; }

      evaluateInverseCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[0]], n.data()[elementMaterialTypes.data()[0]], thetaR.data()[elementMaterialTypes.data()[0]],
                                  thetaSR.data()[elementMaterialTypes.data()[0]], &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                                  pOut.data()[i], mIn.data()[i],
                                  dm, f, df, a, da,
                                  PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                                  BC_entry_head,     // h_e (only used if BC_PSK)
                                  BC_lambda);
    }
  }

  void calculateMassMatrix(arguments_dict &args)
  {
    //element
    double               dt                  = args.scalar<double>("dt");
    xt::pyarray<double> &mesh_trial_ref      = args.array<double>("mesh_trial_ref");
    xt::pyarray<double> &mesh_grad_trial_ref = args.array<double>("mesh_grad_trial_ref");
    xt::pyarray<double> &mesh_dof            = args.array<double>("mesh_dof");
    xt::pyarray<double> &mesh_velocity_dof   = args.array<double>("mesh_velocity_dof");
    double               MOVING_DOMAIN       = args.scalar<double>("MOVING_DOMAIN");
    xt::pyarray<int>    &mesh_l2g            = args.array<int>("mesh_l2g");
    xt::pyarray<double> &dV_ref              = args.array<double>("dV_ref");
    xt::pyarray<double> &u_trial_ref         = args.array<double>("u_trial_ref");
    xt::pyarray<double> &u_grad_trial_ref    = args.array<double>("u_grad_trial_ref");
    xt::pyarray<double> &u_test_ref          = args.array<double>("u_test_ref");
    xt::pyarray<double> &u_grad_test_ref     = args.array<double>("u_grad_test_ref");
    //element boundary
    xt::pyarray<double> &mesh_trial_trace_ref      = args.array<double>("mesh_trial_trace_ref");
    xt::pyarray<double> &mesh_grad_trial_trace_ref = args.array<double>("mesh_grad_trial_trace_ref");
    xt::pyarray<double> &dS_ref                    = args.array<double>("dS_ref");
    xt::pyarray<double> &u_trial_trace_ref         = args.array<double>("u_trial_trace_ref");
    xt::pyarray<double> &u_grad_trial_trace_ref    = args.array<double>("u_grad_trial_trace_ref");
    xt::pyarray<double> &u_test_trace_ref          = args.array<double>("u_test_trace_ref");
    xt::pyarray<double> &u_grad_test_trace_ref     = args.array<double>("u_grad_test_trace_ref");
    xt::pyarray<double> &normal_ref                = args.array<double>("normal_ref");
    xt::pyarray<double> &boundaryJac_ref           = args.array<double>("boundaryJac_ref");
    //physics
    int nElements_global = args.scalar<int>("nElements_global");
    //new
    xt::pyarray<double> &ebqe_penalty_ext     = args.array<double>("ebqe_penalty_ext");
    xt::pyarray<int>    &elementMaterialTypes = args.array<int>("elementMaterialTypes");
    xt::pyarray<int>    &isSeepageFace        = args.array<int>("isSeepageFace");
    xt::pyarray<int>    &a_rowptr             = args.array<int>("a_rowptr");
    xt::pyarray<int>    &a_colind             = args.array<int>("a_colind");
    double               rho                  = args.scalar<double>("rho");
    double               beta                 = args.scalar<double>("beta");
    xt::pyarray<double> &gravity              = args.array<double>("gravity");
    xt::pyarray<double> &alpha                = args.array<double>("alpha");
    xt::pyarray<double> &n                    = args.array<double>("n");
    xt::pyarray<double> &thetaR               = args.array<double>("thetaR");
    xt::pyarray<double> &thetaSR              = args.array<double>("thetaSR");
    xt::pyarray<double> &KWs                  = args.array<double>("KWs");
    //end new
    double               useMetrics                                 = args.scalar<double>("useMetrics");
    double               alphaBDF                                   = args.scalar<double>("alphaBDF");
    int                  lag_shockCapturing                         = args.scalar<int>("lag_shockCapturing");
    double               shockCapturingDiffusion                    = args.scalar<double>("shockCapturingDiffusion");
    xt::pyarray<int>    &u_l2g                                      = args.array<int>("u_l2g");
    xt::pyarray<int>    &r_l2g                                      = args.array<int>("r_l2g");
    xt::pyarray<double> &elementDiameter                            = args.array<double>("elementDiameter");
    int                  degree_polynomial                          = args.scalar<int>("degree_polynomial");
    xt::pyarray<double> &u_dof                                      = args.array<double>("u_dof");
    xt::pyarray<double> &velocity                                   = args.array<double>("velocity");
    xt::pyarray<double> &q_m_betaBDF                                = args.array<double>("q_m_betaBDF");
    xt::pyarray<double> &cfl                                        = args.array<double>("cfl");
    xt::pyarray<double> &q_numDiff_u_last                           = args.array<double>("q_numDiff_u_last");
    xt::pyarray<int>    &csrRowIndeces_u_u                          = args.array<int>("csrRowIndeces_u_u");
    xt::pyarray<int>    &csrColumnOffsets_u_u                       = args.array<int>("csrColumnOffsets_u_u");
    xt::pyarray<double> &globalJacobian                             = args.array<double>("globalJacobian");
    xt::pyarray<double> &delta_x_ij                                 = args.array<double>("delta_x_ij");
    int                  nExteriorElementBoundaries_global          = args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int>    &exteriorElementBoundariesArray             = args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int>    &elementBoundaryElementsArray               = args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int>    &elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
    xt::pyarray<double> &ebqe_velocity_ext                          = args.array<double>("ebqe_velocity_ext");
    xt::pyarray<int>    &isDOFBoundary_u                            = args.array<int>("isDOFBoundary_u");
    xt::pyarray<double> &ebqe_bc_u_ext                              = args.array<double>("ebqe_bc_u_ext");
    xt::pyarray<int>    &isFluxBoundary_u                           = args.array<int>("isFluxBoundary_u");
    xt::pyarray<double> &ebqe_bc_flux_ext                           = args.array<double>("ebqe_bc_flux_ext");
    xt::pyarray<int>    &csrColumnOffsets_eb_u_u                    = args.array<int>("csrColumnOffsets_eb_u_u");
        //////////////////////////////////////For Brooks- Corey/////////////////
    PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
    double BC_entry_head = args.scalar<double>("BC_entry_head");
    double BC_lambda     = args.scalar<double>("BC_lambda");

    
    
    int                  LUMPED_MASS_MATRIX                         = args.scalar<int>("LUMPED_MASS_MATRIX");
    double Ct_sge = 4.0;
    //
    //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
    //
    for (int eN = 0; eN < nElements_global; eN++) {
      double elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
      for (int i = 0; i < nDOF_test_element; i++)
        for (int j = 0; j < nDOF_trial_element; j++) { elementJacobian_u_u[i][j] = 0.0; }
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        int eN_k                  = eN * nQuadraturePoints_element + k, //index to a scalar at a quadrature point
          eN_k_nSpace             = eN_k * nSpace,
            eN_nDOF_trial_element = eN * nDOF_trial_element; //index to a vector at a quadrature point

        //declare local storage
        double u = 0.0, grad_u[nSpace], m = 0.0, dm = 0.0, f[nSpace], df[nSpace], a[nnz], da[nnz], as[nnz], m_t = 0.0, dm_t = 0.0, dpdeResidual_u_u[nDOF_trial_element], Lstar_u[nDOF_test_element], dsubgridError_u_u[nDOF_trial_element], tau = 0.0, tau0 = 0.0, tau1 = 0.0, jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace], u_grad_trial[nDOF_trial_element * nSpace], dV, u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element * nSpace], x, y, z, xt, yt, zt,
          G[nSpace * nSpace], G_dd_G, tr_G;

        //get jacobian, etc for mapping reference element
        ck.calculateMapping_element(eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y, z);
        ck.calculateMappingVelocity_element(eN, k, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), xt, yt, zt);
        //get the physical integration weight
        dV = fabs(jacDet) * dV_ref.data()[k];
        ck.calculateG(jacInv, G, G_dd_G, tr_G);
        //get the trial function gradients
        ck.gradTrialFromRef(&u_grad_trial_ref.data()[k * nDOF_trial_element * nSpace], jacInv, u_grad_trial);
        //get the solution
        ck.valFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], u);
        //get the solution gradients
        ck.gradFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], u_grad_trial, grad_u);
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) {
          u_test_dV[j] = u_test_ref.data()[k * nDOF_trial_element + j] * dV;
          for (int I = 0; I < nSpace; I++) {
            u_grad_test_dV[j * nSpace + I] = u_grad_trial[j * nSpace + I] * dV; //cek warning won't work for Petrov-Galerkin
          }
        }
        //
        //calculate pde coefficients and derivatives at quadrature points
        //
        double Kr, dKr;
        evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u, m, dm, f, df, a, da, as, Kr, dKr,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             BC_entry_head,     // h_e (only used if BC_PSK)
                             BC_lambda);
        //
        //moving mesh
        //
        double mesh_velocity[3];
        mesh_velocity[0] = xt;
        mesh_velocity[1] = yt;
        mesh_velocity[2] = zt;
        for (int I = 0; I < nSpace; I++) {
          f[I] -= MOVING_DOMAIN * m * mesh_velocity[I];
          df[I] -= MOVING_DOMAIN * dm * mesh_velocity[I];
        }
        //
        //calculate time derivatives
        //
        //cek hack
        dm = 1.0;
        ck.bdf(alphaBDF,
               q_m_betaBDF.data()[eN_k], //since m_t isn't used, we don't have to correct mass
               m, dm, m_t, dm_t);
        //
        //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
        //
        //calculate the adjoint times the test functions
        for (int i = 0; i < nDOF_test_element; i++) {
          int i_nSpace = i * nSpace;
          Lstar_u[i]   = ck.Advection_adjoint(df, &u_grad_test_dV[i_nSpace]);
        }
        //calculate the Jacobian of strong residual
        for (int j = 0; j < nDOF_trial_element; j++) {
          int j_nSpace        = j * nSpace;
          dpdeResidual_u_u[j] = ck.MassJacobian_strong(dm_t, u_trial_ref.data()[k * nDOF_trial_element + j]) + ck.AdvectionJacobian_strong(df, &u_grad_trial[j_nSpace]);
        }
        //tau and tau*Res
        calculateSubgridError_tau(elementDiameter.data()[eN], dm_t, df, cfl.data()[eN_k], tau0);

        calculateSubgridError_tau(Ct_sge, G, dm_t, df, tau1, cfl.data()[eN_k]);
        tau = useMetrics * tau1 + (1.0 - useMetrics) * tau0;

        for (int j = 0; j < nDOF_trial_element; j++) dsubgridError_u_u[j] = -tau * dpdeResidual_u_u[j];
        for (int i = 0; i < nDOF_test_element; i++) {
          for (int j = 0; j < nDOF_trial_element; j++) {
            if (LUMPED_MASS_MATRIX == 1) {
              if (i == j) elementJacobian_u_u[i][j] += u_test_dV[i];
            } else {
              int j_nSpace = j * nSpace;
              int i_nSpace = i * nSpace;
              dm_t = 1.0; //we are solving for continuum density explicitly
              elementJacobian_u_u[i][j] += ck.MassJacobian_weak(dm_t, u_trial_ref.data()[k * nDOF_trial_element + j], u_test_dV[i]);
            }
          } //j
        } //i
      } //k
      //
      //load into element Jacobian into global Jacobian
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        int I    = u_l2g.data()[eN_i];
        for (int j = 0; j < nDOF_trial_element; j++) {
          int eN_i_j = eN_i * nDOF_trial_element + j;
          int J      = u_l2g.data()[eN * nDOF_trial_element + j];
          //globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += elementJacobian_u_u[i][j];
          delta_x_ij.data()[3 * (csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]) + 0] = mesh_dof.data()[I * 3 + 0] - mesh_dof.data()[J * 3 + 0];
          delta_x_ij.data()[3 * (csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]) + 1] = mesh_dof.data()[I * 3 + 1] - mesh_dof.data()[J * 3 + 1];
          delta_x_ij.data()[3 * (csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]) + 2] = mesh_dof.data()[I * 3 + 2] - mesh_dof.data()[J * 3 + 2];
        } //j
      } //i
    } //elements
  } //computeMassMatrix
}; //Richards

inline Mphase_co2_base* newmphase_co2(int nSpaceIn, int nQuadraturePoints_elementIn, int nDOF_mesh_trial_elementIn, int nDOF_trial_elementIn, int nDOF_test_elementIn, int nQuadraturePoints_elementBoundaryIn, int CompKernelFlag)
{
  if (nSpaceIn == 1)
    return proteus::chooseAndAllocateDiscretization1D<Mphase_co2_base, mphase_co2, CompKernel>(nSpaceIn, nQuadraturePoints_elementIn, nDOF_mesh_trial_elementIn, nDOF_trial_elementIn, nDOF_test_elementIn, nQuadraturePoints_elementBoundaryIn, CompKernelFlag);
  else if (nSpaceIn == 2)
    return proteus::chooseAndAllocateDiscretization2D<Mphase_co2_base, mphase_co2, CompKernel>(nSpaceIn, nQuadraturePoints_elementIn, nDOF_mesh_trial_elementIn, nDOF_trial_elementIn, nDOF_test_elementIn, nQuadraturePoints_elementBoundaryIn, CompKernelFlag);
  else {
    assert(nSpaceIn == 3);
    return proteus::chooseAndAllocateDiscretization<Mphase_co2_base, mphase_co2, CompKernel>(nSpaceIn, nQuadraturePoints_elementIn, nDOF_mesh_trial_elementIn, nDOF_trial_elementIn, nDOF_test_elementIn, nQuadraturePoints_elementBoundaryIn, CompKernelFlag);
  }
}
} // namespace mphase_co2
} // namespace proteus
#endif
