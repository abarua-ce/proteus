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
//  virtual void kth_FCT_step(arguments_dict &args)                        = 0;
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
                                 const double rho_water, const double rho_air,
                                 const double beta_water , const double beta_air, 
                                 const double gravity[nSpace],
                                 const double alpha, 
                                 const double n_vg,
                                 const double thetaR, const double thetaSR,
                                 const double KWs[nnz],
                                 const double &u_water,  const double &u_air,   // heads
                                 const double Y_water,   const double Y_air,
                                 double &m_water, double &m_air,
                                 double &dm_water, double &dm_air,
                                 double f_water[nSpace], double f_air[nSpace],
                                 double df_water[nSpace], double df_air[nSpace],
                                 double a_water[nnz], double a_air[nnz],
                                 double da_water[nnz], double da_air[nnz],
                                 double as_water[nnz], double as_air[nnz],
                                 double &kr_water, double &dkr_water,
                                 double &kr_air,   double &dkr_air,
                                 const PSK PSK_TYPE,
                                 double &Swater_out, double &Sair_out,
                                 const double BC_entry_head, const double BC_lambda)
{ const int nSpace2 = nSpace * nSpace;
  double psiC;
  double pcBar, pcBar_n, pcBar_nM1, pcBar_nM2, onePlus_pcBar_n;
  double sBar, sqrt_sBar, DsBar_DpsiC;
  double thetaW, DthetaW_DpsiC;
  double vBar, vBar2, DvBar_DpsiC;
  double KWr, DKWr_DpsiC;
  double thetaS, m_vg, pcBarStar, sqrt_sBarStar;
  
  double Se      = 1.0;
  double dSe_dpsic = 0.0;
  thetaS = thetaR + thetaSR;
  psiC = (rho_air / rho_water) * u_air - u_water;
 
  const double dpsic_duw = -1.0;
  const double dpsic_dua = (rho_air / rho_water);

  if (PSK_TYPE == PSK::VG_PSK)
  {
    m_vg = 1.0 - 1.0 / n_vg;

    if (psiC > 0.0) {
      pcBar     = alpha * psiC;
      pcBarStar = (pcBar < 1.0e-8) ? 1.0e-8 : pcBar;

      pcBar_nM2       = std::pow(pcBarStar, n_vg - 2.0);
      pcBar_nM1       = pcBar_nM2 * pcBar;
      pcBar_n         = pcBar_nM1 * pcBar;
      onePlus_pcBar_n = 1.0 + pcBar_n;

      sBar        = std::pow(onePlus_pcBar_n, -m_vg);
      DsBar_DpsiC = alpha * (1.0 - n_vg) * (sBar / onePlus_pcBar_n) * pcBar_nM1;

      vBar        = 1.0 - pcBar_nM1 * sBar;
      vBar2       = vBar * vBar;
      DvBar_DpsiC = -alpha * (n_vg - 1.0) * pcBar_nM2 * sBar - pcBar_nM1 * DsBar_DpsiC;

      thetaW        = thetaSR * sBar + thetaR;
      DthetaW_DpsiC = thetaSR * DsBar_DpsiC;

      sqrt_sBar     = std::sqrt(std::max(sBar, 1.0e-16));
      sqrt_sBarStar = (sqrt_sBar < 1.0e-8) ? 1.0e-8 : sqrt_sBar;

      KWr           = sqrt_sBar * vBar2;                                                                
      DKWr_DpsiC    = (0.5 / sqrt_sBarStar) * DsBar_DpsiC * vBar2 + 2.0 * sqrt_sBar * vBar * DvBar_DpsiC;

      Se            = sBar;
      dSe_dpsic     = DsBar_DpsiC;
    }
    else {
      thetaW        = thetaS;
      DthetaW_DpsiC = 0.0;
      KWr           = 1.0;
      DKWr_DpsiC    = 0.0;

      Se            = 1.0;
      dSe_dpsic     = 0.0;
    }
  }
  else // Brooks–Corey
  {
    const double hb  = (BC_entry_head > 1.0e-12) ? BC_entry_head : 1.0e-12;
    const double lam = (BC_lambda     > 1.0e-12) ? BC_lambda     : 1.0e-12;

    // ---- Se(ψc) and dSe/dψc (UNCHANGED) ----
    if (psiC > 0.0) {
      if (psiC >= hb) {
        Se            = std::pow(hb / psiC, lam);
        dSe_dpsic     = -lam * Se / psiC;

        thetaW        = thetaR + thetaSR * Se;
        DthetaW_DpsiC = thetaSR * dSe_dpsic;

        // ---- k_rℓ : Brooks–Corey wetting rel-perm (same exponent as before) ----
        // CHANGED: annotate formula explicitly: kr_w = Se^{(2+3λ)/λ} = Se^{3 + 2/λ}
        const double expo_w = 3.0 + 2.0 / lam;
        KWr           = std::pow(Se, expo_w);
        DKWr_DpsiC    = expo_w * std::pow(Se, expo_w - 1.0) * dSe_dpsic;

        // ---- k_rg : Brooks–Corey non-wetting per user's formula ----
        // CHANGED: use kr_g = (1 - Se)^2 * (1 - Se^{(2+λ)/λ})
        const double Se_cl   = std::min(std::max(Se, 1.0e-12), 1.0 - 1.0e-12);
        const double expo_bw = (2.0 + lam) / lam;            // = 2/λ + 1
        const double one_m_Se= 1.0 - Se_cl;
        const double Se_bw   = std::pow(Se_cl, expo_bw);
        double Kra, dKra_dSe;

        Kra      = (one_m_Se * one_m_Se) * (1.0 - Se_bw);
        // d/dSe [ (1-Se)^2 * (1 - Se^{bw}) ]
        dKra_dSe = -2.0 * one_m_Se * (1.0 - Se_bw)
                   - (one_m_Se * one_m_Se) * (expo_bw * std::pow(Se_cl, expo_bw - 1.0));

        // pack using existing names (used later)
        // NOTE: dKra_dpsic = dKra_dSe * dSe/dψc, and OWN-variable mapping is applied below
        // We reuse 'dKra_dpsic' slot via local scope after this block
        // To keep structure, store to temporary file-scope variables:
        // We'll compute dKra_dpsic below outside the if's after Se known
        // For now stash Kra into 'kr_air' via locals:
        kr_air = Kra;
        dkr_air = dKra_dSe * dSe_dpsic;  // d(kr_g)/dψc (OWN mapping applied later)

      } else {
        // saturated
        Se            = 1.0;  dSe_dpsic = 0.0;
        thetaW        = thetaS; DthetaW_DpsiC = 0.0;
        KWr           = 1.0;  DKWr_DpsiC = 0.0;

        // CHANGED: gas rel-perm when Se=1 -> kr_g = 0
        kr_air        = 0.0;  dkr_air    = 0.0;
      }
    } else {
      Se            = 1.0;  dSe_dpsic = 0.0;
      thetaW        = thetaS; DthetaW_DpsiC = 0.0;
      KWr           = 1.0;  DKWr_DpsiC = 0.0;

      // CHANGED: gas rel-perm when Se=1 -> kr_g = 0
      kr_air        = 0.0;  dkr_air    = 0.0;
    }

    // If we are in the general case above (ψc ≥ hb), kr_air and dkr_air are already set.
    // Otherwise they are zero from the branches.
  }

  // --- Non-wetting rel-perm & derivative for VG path or if BC computed locally ---
  // For VG branch we already had Kra via vg-formula below; for BC we put it in kr_air/dkr_air above.
  const double Se_cl = std::min(std::max(Se, 1.0e-12), 1.0 - 1.0e-12);
  double Kra, dKra_dSe;

  if (PSK_TYPE == PSK::VG_PSK) {
    const double mvg    = 1.0 - 1.0 / n_vg;
    const double Se_pow = std::pow(Se_cl, 1.0 / mvg);
    const double t1     = std::sqrt(1.0 - Se_cl);
    const double t2     = 1.0 - Se_pow;

    Kra      = t1 * std::pow(t2, 2.0 * mvg);

    const double dt1_dSe = -0.5 / std::max(t1, 1.0e-8);
    const double dt2_dSe = -(1.0 / mvg) * std::pow(Se_cl, 1.0 / mvg - 1.0);
    dKra_dSe = dt1_dSe * std::pow(t2, 2.0 * mvg)
             + t1 * (2.0 * mvg) * std::pow(t2, 2.0 * mvg - 1.0) * dt2_dSe;
  }
  else { 
    // BC branch: if we already computed kr_air/dkr_air above, reuse it; otherwise compute here
    if (psiC > 0.0 && psiC >= ((BC_entry_head > 1.0e-12) ? BC_entry_head : 1.0e-12)) {
      // Already set in the BC block
      Kra = kr_air;
      dKra_dSe = (dkr_air / std::max(dSe_dpsic, 1e-30)); // back out dKra/dSe safely
    } else {
      Kra      = 0.0;
      dKra_dSe = 0.0;
    }
  }
  const double dKra_dpsic = dKra_dSe * dSe_dpsic;

  // ---------- Phase saturations ----------
  const double Sw = thetaW / std::max(thetaS, 1.0e-12);
  const double Sg = 1.0 - Sw;
  Swater_out = Sw;
  Sair_out   = Sg;

  const double thetaA        = thetaS - thetaW;
  const double dthetaA_dpsic = -DthetaW_DpsiC;

  // ---------- Slight compressibility (per phase) ----------
  const double rho_w  = rho_water * std::exp(beta_water * u_water);
  const double drho_w = beta_water * rho_w;
  const double rho_a  = rho_air   * std::exp(beta_air * u_air);
  const double drho_a = beta_air  * rho_a;

  // chain rule to own primary variable
  const double dthetaW_duw = DthetaW_DpsiC * dpsic_duw;
  const double dthetaW_dua = DthetaW_DpsiC * dpsic_dua;
  const double dthetaA_duw = dthetaA_dpsic * dpsic_duw; // = -dthetaW_duw
  const double dthetaA_dua = dthetaA_dpsic * dpsic_dua; // = -dthetaW_dua

  const double dKWr_duw = DKWr_DpsiC * dpsic_duw;
  const double dKWr_dua = DKWr_DpsiC * dpsic_dua;
  const double dKra_duw = dKra_dpsic * dpsic_duw;
  const double dKra_dua = dKra_dpsic * dpsic_dua;

  // ---------- Mass and its derivative (with optional mass fractions Y_k) ----------
  m_water  = rho_w * thetaW * Y_water;
  m_air    = rho_a * thetaA * Y_air;

  dm_water = drho_w * thetaW * Y_water + rho_w * dthetaW_duw * Y_water; // wrt u_water
  dm_air   = drho_a * thetaA * Y_air   + rho_a * dthetaA_dua * Y_air;   // wrt u_air

  // ---------- Flux-like terms (keep original loops/shape) ----------
  for (int I = 0; I < nSpace; ++I) {
    f_water[I] = 0.0; df_water[I] = 0.0;
    f_air[I]   = 0.0; df_air[I]   = 0.0;
  }
  for (int ii = 0; ii < nnz; ++ii) {
    a_water[ii] = 0.0;  da_water[ii] = 0.0;  as_water[ii] = 0.0;
    a_air[ii]   = 0.0;  da_air[ii]   = 0.0;  as_air[ii]   = 0.0;
  }

  const double rho2_w = rho_w * rho_w;
  const double rho2_a = rho_a * rho_a;

  for (int I = 0; I < nSpace; ++I) {
    for (int ii = rowptr[I]; ii < rowptr[I + 1]; ++ii) {
      const int J = colind[ii];

      // wetting (water)
      f_water[I]  += rho2_w * KWr * KWs[ii] * gravity[J];
      df_water[I] += rho2_w * dKWr_duw * KWs[ii] * gravity[J];

      a_water[ii]  = rho_w * KWr * KWs[ii];
      da_water[ii] = rho_w * dKWr_duw * KWs[ii];

      as_water[ii] = rho_w * KWs[ii];

      // non-wetting (air)
      f_air[I]  += rho2_a * Kra * KWs[ii] * gravity[J];
      df_air[I] += rho2_a * dKra_dua * KWs[ii] * gravity[J];

      a_air[ii]  = rho_a * Kra * KWs[ii];
      da_air[ii] = rho_a * dKra_dua * KWs[ii];

      as_air[ii] = rho_a * KWs[ii];
    }
  }

  // ---------- Report rel-perms and OWN-variable derivatives ----------
  kr_water = KWr;  dkr_water = dKWr_duw;
  kr_air   = Kra;  dkr_air   = dKra_dua;
}



inline void evaluateInverseCoefficients_2ph(const int rowptr[nSpace], const int colind[nnz],
                                            const double rho_water, const double rho_air,
                                            const double beta_water, const double beta_air,   // kept for parity (unused)
                                            const double gravity[nSpace],
                                            const double alpha, const double n_vg,
                                            const double thetaR, const double thetaSR,
                                            const double KWs[nnz],
                                            double &u_water, double &u_air,                   // IN/OUT: both heads
                                            const double &m_water, const double &m_air,       // IN: masses for both phases
                                            const double Y_water, const double Y_air,         // mass fractions
                                            const PSK PSK_TYPE,
                                            const double BC_entry_head, 
                                            const double BC_lambda)
{
  const double eps    = 1.0e-12;
  const double thetaS = thetaR + thetaSR;

  // --- Invert masses to phase contents (ignore compressibility in inverse, like before) ---
  const double denom_w = std::fmax(rho_water * std::fmax(Y_water, eps), eps);
  const double denom_a = std::fmax(rho_air   * std::fmax(Y_air,   eps), eps);

  const double thetaW_from_w = m_water / denom_w;           
  const double thetaA_from_a = m_air   / denom_a;            
  const double thetaW_from_a = thetaS - thetaA_from_a;      

  // choose a stable θw candidate, prefer water mass if reasonable
  const bool w_ok = (thetaW_from_w > thetaR - 1e-8) && (thetaW_from_w < thetaS + 1e-8);
  const bool a_ok = (thetaW_from_a > thetaR - 1e-8) && (thetaW_from_a < thetaS + 1e-8);

  double thetaW = w_ok ? thetaW_from_w : (a_ok ? thetaW_from_a : thetaW_from_w);
  thetaW = std::fmin(std::fmax(thetaW, thetaR + eps), thetaS - eps);

  // --- Effective saturation ---
  double Se = (thetaW - thetaR) / std::fmax(thetaSR, eps);
  Se = std::fmin(std::fmax(Se, eps), 1.0 - eps);

  // --- Invert PSK to capillary head ψc ---
  double psiC = 0.0;
  if (PSK_TYPE == PSK::VG_PSK) {
    const double m_vg = 1.0 - 1.0 / n_vg;
    const double pc_n = std::pow(Se, -1.0 / m_vg) - 1.0;
    const double pc   = std::pow(std::fmax(pc_n, 0.0), 1.0 / n_vg);
    psiC = pc / std::fmax(alpha, eps);
  } else {
    const double hb  = (BC_entry_head > eps) ? BC_entry_head : eps;
    const double lam = (BC_lambda     > eps) ? BC_lambda     : eps;
    psiC = hb * std::pow(Se, -1.0 / lam);
  }

  double uw_new = (rho_air / rho_water) * u_air - psiC;     // update water from previous air
  double ua_new = (psiC + uw_new) * (rho_water / rho_air);  // then air from updated water

  u_water = uw_new;
  u_air   = ua_new;


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
      const double Se_clamped = std::fmax(std::fmin(sBar, 1.0-1e-12), 1e-12);
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
    //double BC_entry_head = args.scalar<double>("BC_entry_head");
    //double BC_lambda     = args.scalar<double>("BC_lambda");

    xt::pyarray<double> &BC_entry_head = args.array<double>("BC_entry_head"); // size: nMaterials
    xt::pyarray<double> &BC_lambda     = args.array<double>("BC_lambda");     // size: nMaterials





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
    xt::pyarray<int>    &u_l2g_water                                = args.array<int>("u_l2g_water");
    xt::pyarray<int>    &u_l2g_air                                  = args.array<int>("u_l2g_air");   
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
    xt::pyarray<double> &q_dV                                       = args.array<double>("q_dV");
    
 //   xt::pyarray<double> &q_dV_water                                 = args.array<double>("q_dV_water");
 //   xt::pyarray<double> &q_dV_air                                   = args.array<double>("q_dV_air");
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
    xt::pyarray<double> &globalResidual                       = args.array<double>("globalResidual");
    // xt::pyarray<double> &globalResidual                       = args.array<double>("globalResidual");
    // xt::pyarray<double> &globalResidual                         = args.array<double>("globalResidual");
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
    //xt::pyarray<double> &ebqe_phi_water                             = args.array<double>("ebqe_phi_water");
    //xt::pyarray<double> &ebqe_phi_air                               = args.array<double>("ebqe_phi_air");
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
        double tau_water=0.0, tau0_water=0.0, tau1_water=0.0, numDiff0_water=0.0, numDiff1_water=0.0;
        double pdeResidual_u_water = 0.0, Lstar_u_water[nDOF_test_element], subgridError_u_water = 0.0;

        // fill water, assemble to elementResidual_u_water[...]

        // --- phase air ---
        double u_air=0.0, grad_u_air[nSpace];
        double m_air=0.0, dm_air=0.0, m_t_air=0.0, dm_t_air=0.0;
        double f_air[nSpace], df_air[nSpace];
        double a_air[nnz],   da_air[nnz],   as_air[nnz];
        double Kr_air=0.0, dKr_air=0.0;
        double pdeResidual_u_air = 0.0, Lstar_u_air[nDOF_test_element], subgridError_u_air = 0.0;
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
                      &u_l2g_water.data()[eN_nDOF_trial_element],     // map for WATER
                      &u_trial_ref.data()[k * nDOF_trial_element],
                      u_water);
        ck.valFromDOF(u_dof_air.data(),                             // DOFs for AIR
                      &u_l2g_air.data()[eN_nDOF_trial_element],     // map for AIR
                      &u_trial_ref.data()[k * nDOF_trial_element],
                      u_air);
                        
        //get the solution gradient for both phases

        ck.gradFromDOF(u_dof_water.data(),
                       &u_l2g_water.data()[eN_nDOF_trial_element], 
                       u_grad_trial, 
                       grad_u_water);

        ck.gradFromDOF(u_dof_air.data(),
                       &u_l2g_air.data()[eN_nDOF_trial_element], 
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
                            rho_water, rho_air, 
                            beta_water, beta_air, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            u_water, u_air, 
                            Y_water, Y_air,
                            m_water, m_air, 
                            dm_water, dm_air,
                            f_water, f_air,
                            df_water, df_air,
                            a_water, a_air,
                            da_water, da_air, 
                            as_water, as_air,
                            Kr_water, dKr_water,
                            Kr_air, dKr_air,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw, Sg,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));       
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
        q_numDiff_u_water[eN_k] = useMetrics * numDiff1_water + (1.0 - useMetrics) * numDiff0_water;
        //phase air
        ck.calculateNumericalDiffusion(shockCapturingDiffusion, elementDiameter[eN], pdeResidual_u_air, grad_u_air, numDiff0_air);
        ck.calculateNumericalDiffusion(shockCapturingDiffusion, sc_uref, sc_alpha, G, G_dd_G, pdeResidual_u_air, grad_u_air, numDiff1_air);
        q_numDiff_u_air[eN_k] = useMetrics * numDiff1_air + (1.0 - useMetrics) * numDiff0_air;        
        
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
        q_u_air.data()[eN_k] = u_air;
      }
      //
      //load element into global residual and save element residual
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        globalResidual.data()[offset_u_water + stride_u_water * u_l2g_water.data()[eN_i]] += elementResidual_u_water[i];
        globalResidual.data()[offset_u_air + stride_u_air * u_l2g_air.data()[eN_i]] += elementResidual_u_air[i];
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
        double jac_ext[nSpace * nSpace], jacDet_ext, jacInv_ext[nSpace * nSpace], boundaryJac[nSpace * (nSpace - 1)], metricTensor[(nSpace - 1) * (nSpace - 1)], metricTensorDetSqrt;// metricTensor, metricTensorDetSqrt;   
        double normal[3], x_ext , y_ext, z_ext,  xt_ext, yt_ext, zt_ext, integralScaling;    
        double dS, u_test_dS[nDOF_test_element], u_grad_trial_trace[nDOF_trial_element * nSpace], G[nSpace * nSpace], G_dd_G, tr_G;
               //fluxJacobian_u_u[nDOF_trial_element], bfluxJacobian_u_u[nDOF_trial_element], fluxJacobian_un_un[nDOF_trial_element];
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
        ck.calculateMapping_elementBoundary(eN, ebN_local, kb, ebN_local_kb, mesh_dof.data(), mesh_l2g.data(), mesh_trial_trace_ref.data(), mesh_grad_trial_trace_ref.data(), 
                                            boundaryJac_ref.data(), jac_ext, jacDet_ext, jacInv_ext, boundaryJac, metricTensor, metricTensorDetSqrt,
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
                      &u_l2g_water.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_water);
        ck.valFromDOF(u_dof_air.data(), 
                      &u_l2g_air.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_air);
        
        //Gradient
        ck.gradFromDOF(u_dof_water.data(), 
                       &u_l2g_water.data()[eN_nDOF_trial_element], 
                       u_grad_trial_trace, 
                       grad_u_ext_water);
        ck.gradFromDOF(u_dof_air.data(), 
                       &u_l2g_air.data()[eN_nDOF_trial_element], 
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
        //calculate the pde coefficients using the solution and the boundary values for the solution
        //
        double Kr_water, dKr_water, Kr_air, dKr_air;
        double Sw_ext=0.0;
        double Sg_ext =0.0;
        double bc_Kr_water, bc_dKr_water ,bc_Kr_ext_water, bc_dKr_ext_water, bc_Krn_water, bc_dKrn_water;
        double bc_Kr_air, bc_dKr_air ,bc_Kr_ext_air, bc_dKr_ext_air, bc_Krn_air, bc_dKrn_air;
        
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_ext_water, u_ext_air, 
                             Y_water, Y_air,
                             m_ext_water, m_ext_air, 
                             dm_ext_water, dm_ext_air,
                             f_ext_water, f_ext_air,
                             df_ext_water, df_ext_air,
                             a_ext_water, a_ext_air,
                             da_ext_water, da_ext_air, 
                             as_ext_water, as_ext_air,
                             Kr_water, dKr_water,
                             Kr_air, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                             BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                             BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));

        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             bc_u_ext_water, bc_u_ext_air,
                             Y_water, Y_air, 
                             bc_m_ext_water, bc_m_ext_air, 
                             bc_dm_ext_water, bc_dm_ext_air, 
                             bc_f_ext_water, bc_f_ext_air, 
                             bc_df_ext_water, bc_df_ext_air, 
                             bc_a_ext_water, bc_a_ext_air, 
                             bc_da_ext_water, bc_da_ext_air, 
                             bc_as_ext_water, bc_as_ext_air,
                             Kr_water, dKr_water,
                             Kr_air, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                             BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                             BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));
     // λ   (only used if BC_PSK));
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
        ebqe_u_water.data()[ebNE_kb]     = u_ext_water;
        ebqe_u_air.data()[ebNE_kb]       = u_ext_air;
        
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
        globalResidual.data()[offset_u_water + stride_u_water * u_l2g_water.data()[eN_i]] += elementResidual_u_water[i];
        globalResidual.data()[offset_u_air + stride_u_air * u_l2g_water.data()[eN_i]] += elementResidual_u_air[i];        
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
    double Y_air                                                    = args.scalar<double>("Y_air");     // e.g. 1.0
    double Y_water                                                  = args.scalar<double>("Y_water");   // e.g. 1.0   
    xt::pyarray<int>    &u_l2g_water                                = args.array<int>("u_l2g_water");
    xt::pyarray<int>    &u_l2g_air                                  = args.array<int>("u_l2g_air");   
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
    xt::pyarray<double> &globalJacobian                       = args.array<double>("globalJacobian");
    //xt::pyarray<double> &globalJacobian                         = args.array<double>("globalJacobian");
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
    // double BC_entry_head = args.scalar<double>("BC_entry_head");
    // double BC_lambda     = args.scalar<double>("BC_lambda");
    xt::pyarray<double> &BC_entry_head = args.array<double>("BC_entry_head"); // size: nMaterials
    xt::pyarray<double> &BC_lambda     = args.array<double>("BC_lambda");     // size: nMaterials

    //
    //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
    //
    for (int eN = 0; eN < nElements_global; eN++) {
      double elementJacobian_u_u_water[nDOF_test_element][nDOF_trial_element];
      double elementJacobian_u_u_air[nDOF_test_element][nDOF_trial_element];

      for (int i = 0; i < nDOF_test_element; i++) {
        for (int j = 0; j < nDOF_trial_element; j++) { 
          elementJacobian_u_u_water[i][j] = 0.0;
          elementJacobian_u_u_air[i][j] = 0.0;
         }
      }
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        int eN_k                  = eN * nQuadraturePoints_element + k, //index to a scalar at a quadrature point
          eN_k_nSpace             = eN_k * nSpace,
          eN_nDOF_trial_element   = eN * nDOF_trial_element; //index to a vector at a quadrature point

        //declare local storage
        double  dV, x, y, z, xt, yt, zt;
        double  jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace];
        double  G[nSpace * nSpace], G_dd_G, tr_G;
        double  u_grad_trial[nDOF_trial_element * nSpace], u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element * nSpace];
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
                      &u_l2g_water.data()[eN_nDOF_trial_element], 
                      &u_trial_ref.data()[k * nDOF_trial_element], 
                      u_water);
        
        ck.valFromDOF(u_dof_air.data(), 
                      &u_l2g_air.data()[eN_nDOF_trial_element], 
                      &u_trial_ref.data()[k * nDOF_trial_element], 
                      u_air);
        
        //get the solution gradients for both phases
        ck.gradFromDOF(u_dof_water.data(),
                       &u_l2g_water.data()[eN_nDOF_trial_element], 
                       u_grad_trial, 
                       grad_u_water);

        ck.gradFromDOF(u_dof_air.data(),
                       &u_l2g_air.data()[eN_nDOF_trial_element], 
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
        double Sw=0.0, Sg=0.0;
        // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
        //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u, m, dm, f, df, a, da, as, Kr, dKr);
        
         evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_water, u_air, 
                             Y_water, Y_air,
                             m_water, m_air, 
                             dm_water, dm_air,
                             f_water, f_air,
                             df_water, df_air,
                             a_water, a_air,
                             da_water, da_air, 
                             as_water, as_air,
                             Kr_water, dKr_water, 
                             Kr_air, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw, Sg,
                             BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                             BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));


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
          globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u_water[i][j];
          globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u_air[i][j];

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
                      &u_l2g_water.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_water);
        
        ck.valFromDOF(u_dof_air.data(), 
                      &u_l2g_air.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_air);


        // Gradient for both phases
        ck.gradFromDOF(u_dof_water.data(), 
                        &u_l2g_water.data()[eN_nDOF_trial_element], 
                        u_grad_trial_trace, 
                        grad_u_ext_water);
        ck.gradFromDOF(u_dof_air.data(), 
                        &u_l2g_air.data()[eN_nDOF_trial_element], 
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
        double Sw_ext=0.0, Sg_ext=0.0;
        
 
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_ext_water, u_ext_air, 
                             Y_water, Y_air,
                             m_ext_water, m_ext_air, 
                             dm_ext_water, dm_ext_air,
                             f_ext_water, f_ext_air,
                             df_ext_water, df_ext_air,
                             a_ext_water, a_ext_air,
                             da_ext_water, da_ext_air, 
                             as_ext_water, as_ext_air,
                             Kr_water, dKr_water,
                             Kr_air, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));

        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             bc_u_ext_water, bc_u_ext_air,
                             Y_water, Y_air, 
                             bc_m_ext_water, bc_m_ext_air, 
                             bc_dm_ext_water, bc_dm_ext_air, 
                             bc_f_ext_water, bc_f_ext_air, 
                             bc_df_ext_water, bc_df_ext_air, 
                             bc_a_ext_water, bc_a_ext_air, 
                             bc_da_ext_water, bc_da_ext_air, 
                             bc_as_ext_water, bc_as_ext_air,
                             Kr_water, dKr_water,
                             Kr_air, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));


        
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
            globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_water[j] * u_test_dS[i];
            globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_air[j] * u_test_dS[i];
          } //j
        } //i
      } //kb
    } //ebNE
  } //computeJacobian

  void FCTStep(arguments_dict &args)
  {
    xt::pyarray<double> &bc_mask_water              = args.array<double>("bc_mask_water");
    xt::pyarray<double> &bc_mask_air                = args.array<double>("bc_mask_air");
    
    int                  NNZ                       = args.scalar<int>("NNZ");     //number on non-zero entries on sparsity pattern
    int                  numDOFs                   = args.scalar<int>("numDOFs"); //number of DOFs
    double               dt                        = args.scalar<double>("dt");
    xt::pyarray<int>    &csrRowIndeces_DofLoops_water    = args.array<int>("csrRowIndeces_DofLoops_water");    //csr row indeces
    xt::pyarray<int>    &csrColumnOffsets_DofLoops_water = args.array<int>("csrColumnOffsets_DofLoops_water"); //csr column offsets
    xt::pyarray<int>    &csrRowIndeces_DofLoops_air    = args.array<int>("csrRowIndeces_DofLoops_water");    //csr row indeces
    xt::pyarray<int>    &csrColumnOffsets_DofLoops_air = args.array<int>("csrColumnOffsets_DofLoops_air"); //csr column offsets


    //flags
    int                  LUMPED_MASS_MATRIX        = args.scalar<int>("LUMPED_MASS_MATRIX");
    int                  MONOLITHIC                = args.scalar<int>("MONOLITHIC");

    //////=====================Water Phase=====================////////////
    xt::pyarray<double> &ML_water                  = args.array<double>("ML_water"); //lumped mass matrix (as vector)
    xt::pyarray<double> &mn_water                  = args.array<double>("mn_water");               //DOFs of solution at time tn
    xt::pyarray<double> &mHigh_water               = args.array<double>("mHigh_water");               //DOFs of high order solution at tnp1
    xt::pyarray<double> &mLow_water                      = args.array<double>("mLow_water");
    xt::pyarray<double> &mDotHigh_water                  = args.array<double>("mDotHigh_water");               //DOFs of high order solution at tnp1
    xt::pyarray<double> &mDotLow_water                   = args.array<double>("mDotLow_water");
    xt::pyarray<double> &limited_solution_water          = args.array<double>("limited_solution_water");
    xt::pyarray<double> &MC_water                        = args.array<double>("MC_water");             //mass matrix
    xt::pyarray<double> &dt_times_fH_minus_fL_water      = args.array<double>("dt_times_fH_minus_fL_water");   //low minus high order dissipative matrices
    xt::pyarray<double> &min_m_bc_water                  = args.array<double>("min_m_bc_water");               //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
    xt::pyarray<double> &max_m_bc_water                  = args.array<double>("max_m_bc_water");
    xt::pyarray<double> &fluxCorrection_water             = args.array<double>("fluxCorrection_water");
   
    double               Rpos_water[numDOFs], Rneg_water[numDOFs];
    double               FluxCorrectionMatrix_water[NNZ] ;
    double               mDot_water[numDOFs];  

    //////////////////
    // LOOP in DOFs //
    //////////////////
    int ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      mDot_water[i] = (mLow_water.data()[i] - mn_water.data()[i])/dt;
      //cek todo: add boundary data--these are just initialized
      //will need to pass p_bc at DOF and calc M
      double mini=min_m_bc_water.data()[i], maxi=max_m_bc_water.data()[i];
      //we're doing local FCT
      //if (GLOBAL_FCT == 1) {
      //  mini = 0.;
      //  maxi = 1.;
      //}

      double Pposi = 0, Pnegi = 0;
      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops_water.data()[i]; offset < csrRowIndeces_DofLoops_water.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops_water.data()[offset];
        ////////////////////////
        // COMPUTE THE BOUNDS //
        ////////////////////////
        if (GLOBAL_FCT == 0) {
          if (MONOLITHIC == 0) {
            mini = fmin(mini, mLow_water[j]);
            maxi = fmax(maxi, mLow_water[j]);
          } else {
            mini = fmin(mini, mn_water.data()[j]);
            maxi = fmax(maxi, mn_water.data()[j]);
          }
        }
        // i-th row of flux correction matrix
        //double I_plus_ML_minus_MC = (i == j ? 1. : 0.) * (1. + ML.data()[i]) - MC.data()[ij];
        //mDot[i] += I_plus_ML_minus_MC * (mHigh.data()[j] - mn.data()[j]) / ML.data()[i];
        mDot_water[j] = (mLow_water.data()[j] - mn_water.data()[j])/dt;
        if (MONOLITHIC == 0) {
          FluxCorrectionMatrix_water[ij] = (LUMPED_MASS_MATRIX == 1 ? 0. : 1.) * dt * MC_water.data()[ij] * (mDotLow_water.data()[i] - mDotLow_water.data()[j]) + dt_times_fH_minus_fL_water.data()[ij];
        } else {
          FluxCorrectionMatrix_water[ij] = dt_times_fH_minus_fL_water.data()[ij];
        }
        ///////////////////////
        // COMPUTE P VECTORS //
        ///////////////////////
        Pposi += FluxCorrectionMatrix_water[ij] * ((FluxCorrectionMatrix_water[ij] > 0) ? 1. : 0.);
        Pnegi += FluxCorrectionMatrix_water[ij] * ((FluxCorrectionMatrix_water[ij] < 0) ? 1. : 0.);

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
        Qposi = ML_water.data()[i] * (maxi - mLow_water[i]);
        Qnegi = ML_water.data()[i] * (mini - mLow_water[i]);
      } else {
        //cek todo: don't think this is right for Richards
        gamma = 10.0 * ML_water.data()[i];
        Qposi = fmin(0.5 * ML_water.data()[i] * (1.0 - mn_water.data()[i]), gamma * (maxi - mn_water[i]));
        Qnegi = fmax(0.5 * ML_water.data()[i] * (0.0 - mn_water.data()[i]), gamma * (mini - mn_water[i]));
      }
      ///////////////////////
      // COMPUTE R VECTORS //
      ///////////////////////
      Rpos_water[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
      Rneg_water[i] = ((Pnegi == 0) ? 1. : fmin(1.0, Qnegi / Pnegi));
    } // i DOFs

    //////////////////////
    // COMPUTE LIMITERS //
    //////////////////////
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      double ith_Limiter_times_FluxCorrectionMatrix = 0.;
      double alpha_fA, alpha_dot, beta_ij = 1.0;
      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops_water.data()[i]; offset < csrRowIndeces_DofLoops_water.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops_water.data()[offset];
        alpha_fA     = ((FluxCorrectionMatrix_water[ij] > 0) ? fmin(Rpos_water[i], Rneg_water[j]) : fmin(Rneg_water[i], Rpos_water[j])) * FluxCorrectionMatrix_water[ij];
        alpha_dot    = fmin(1.0, beta_ij * fabs(alpha_fA) / MC_water.data()[ij] / fmax(1.0e-8, fabs(mDot_water[i] - mDot_water[j])));
        if (MONOLITHIC == 0) {
          ith_Limiter_times_FluxCorrectionMatrix += alpha_fA;
        } else {
          ith_Limiter_times_FluxCorrectionMatrix += alpha_fA + (LUMPED_MASS_MATRIX == 1 ? 0. : 1.) * dt * alpha_dot * MC_water.data()[ij] * (mDot_water[i] - mDot_water[j]);
        }
        ij += 1;

      
      }

      fluxCorrection_water.data()[i] = -ith_Limiter_times_FluxCorrectionMatrix*bc_mask_water[i]/dt;
      limited_solution_water.data()[i] = mLow_water[i] + 1. / ML_water.data()[i] * ith_Limiter_times_FluxCorrectionMatrix * bc_mask_water[i];

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


    //////=====================Water Phase=====================////////////
    xt::pyarray<double> &ML_air                  = args.array<double>("ML_air"); //lumped mass matrix (as vector)
    xt::pyarray<double> &mn_air                  = args.array<double>("mn_air");               //DOFs of solution at time tn
    xt::pyarray<double> &mHigh_air               = args.array<double>("mHigh_air");               //DOFs of high order solution at tnp1
    xt::pyarray<double> &mLow_air                      = args.array<double>("mLow_air");
    xt::pyarray<double> &mDotHigh_air                  = args.array<double>("mDotHigh_air");               //DOFs of high order solution at tnp1
    xt::pyarray<double> &mDotLow_air                   = args.array<double>("mDotLow_air");
    xt::pyarray<double> &limited_solution_air          = args.array<double>("limited_solution_air");
    xt::pyarray<double> &MC_air                        = args.array<double>("MC_air");             //mass matrix
    xt::pyarray<double> &dt_times_fH_minus_fL_air      = args.array<double>("dt_times_fH_minus_fL_air");   //low minus high order dissipative matrices
    xt::pyarray<double> &min_m_bc_air                  = args.array<double>("min_m_bc_air");               //min/max value at BCs. If DOF is not at boundary then min=1E10, max=-1E10
    xt::pyarray<double> &max_m_bc_air                  = args.array<double>("max_m_bc_air");
    xt::pyarray<double> &fluxCorrection_air             = args.array<double>("fluxCorrection_air");
   
    double               Rpos_air[numDOFs], Rneg_air[numDOFs];
    double               FluxCorrectionMatrix_air[NNZ];
    double               mDot_air[numDOFs];
//    double               mDot_air[numDOFs];
    

    //////////////////
    // LOOP in DOFs //
    //////////////////
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      mDot_air[i] = (mLow_air.data()[i] - mn_air.data()[i])/dt;
      //cek todo: add boundary data--these are just initialized
      //will need to pass p_bc at DOF and calc M
      double mini=min_m_bc_air.data()[i], maxi=max_m_bc_air.data()[i];
      //we're doing local FCT
      //if (GLOBAL_FCT == 1) {
      //  mini = 0.;
      //  maxi = 1.;
      //}

      double Pposi = 0, Pnegi = 0;
      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops_air.data()[i]; offset < csrRowIndeces_DofLoops_air.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops_air.data()[offset];
        ////////////////////////
        // COMPUTE THE BOUNDS //
        ////////////////////////
        if (GLOBAL_FCT == 0) {
          if (MONOLITHIC == 0) {
            mini = fmin(mini, mLow_air[j]);
            maxi = fmax(maxi, mLow_air[j]);
          } else {
            mini = fmin(mini, mn_air.data()[j]);
            maxi = fmax(maxi, mn_air.data()[j]);
          }
        }
        // i-th row of flux correction matrix
        //double I_plus_ML_minus_MC = (i == j ? 1. : 0.) * (1. + ML.data()[i]) - MC.data()[ij];
        //mDot[i] += I_plus_ML_minus_MC * (mHigh.data()[j] - mn.data()[j]) / ML.data()[i];
        mDot_air[j] = (mLow_air.data()[j] - mn_air.data()[j])/dt;
        if (MONOLITHIC == 0) {
          FluxCorrectionMatrix_air[ij] = (LUMPED_MASS_MATRIX == 1 ? 0. : 1.) * dt * MC_air.data()[ij] * (mDotLow_air.data()[i] - mDotLow_air.data()[j]) + dt_times_fH_minus_fL_air.data()[ij];
        } else {
          FluxCorrectionMatrix_air[ij] = dt_times_fH_minus_fL_air.data()[ij];
        }
        ///////////////////////
        // COMPUTE P VECTORS //
        ///////////////////////
        Pposi += FluxCorrectionMatrix_air[ij] * ((FluxCorrectionMatrix_air[ij] > 0) ? 1. : 0.);
        Pnegi += FluxCorrectionMatrix_air[ij] * ((FluxCorrectionMatrix_air[ij] < 0) ? 1. : 0.);

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
        Qposi = ML_air.data()[i] * (maxi - mLow_air[i]);
        Qnegi = ML_air.data()[i] * (mini - mLow_air[i]);
      } else {
        //cek todo: don't think this is right for Richards
        gamma = 10.0 * ML_air.data()[i];
        Qposi = fmin(0.5 * ML_air.data()[i] * (1.0 - mn_air.data()[i]), gamma * (maxi - mn_air[i]));
        Qnegi = fmax(0.5 * ML_air.data()[i] * (0.0 - mn_air.data()[i]), gamma * (mini - mn_air[i]));
      }
      ///////////////////////
      // COMPUTE R VECTORS //
      ///////////////////////
      Rpos_air[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
      Rneg_air[i] = ((Pnegi == 0) ? 1. : fmin(1.0, Qnegi / Pnegi));
    } // i DOFs

    //////////////////////
    // COMPUTE LIMITERS //
    //////////////////////
    ij = 0;
    for (int i = 0; i < numDOFs; i++) {
      double ith_Limiter_times_FluxCorrectionMatrix_air = 0.;
      double alpha_fA, alpha_dot, beta_ij = 1.0;
      // LOOP OVER THE SPARSITY PATTERN (j-LOOP)//
      for (int offset = csrRowIndeces_DofLoops_air.data()[i]; offset < csrRowIndeces_DofLoops_air.data()[i + 1]; offset++) {
        int j = csrColumnOffsets_DofLoops_air.data()[offset];
        alpha_fA     = ((FluxCorrectionMatrix_air[ij] > 0) ? fmin(Rpos_air[i], Rneg_air[j]) : fmin(Rneg_air[i], Rpos_air[j])) * FluxCorrectionMatrix_air[ij];
        alpha_dot    = fmin(1.0, beta_ij * fabs(alpha_fA) / MC_air.data()[ij] / fmax(1.0e-8, fabs(mDot_air[i] - mDot_air[j])));
        if (MONOLITHIC == 0) {
          ith_Limiter_times_FluxCorrectionMatrix_air += alpha_fA;
        } else {
          ith_Limiter_times_FluxCorrectionMatrix_air += alpha_fA + (LUMPED_MASS_MATRIX == 1 ? 0. : 1.) * dt * alpha_dot * MC_air.data()[ij] * (mDot_air[i] - mDot_air[j]);
        }
        ij += 1;

      }     
      fluxCorrection_air.data()[i] = -ith_Limiter_times_FluxCorrectionMatrix_air *bc_mask_air[i]/dt;
      limited_solution_air.data()[i] = mLow_air[i] + 1. / ML_air.data()[i] * ith_Limiter_times_FluxCorrectionMatrix_air * bc_mask_air[i];
  }
}
// void kth_FCT_step(arguments_dict &args)
// {
//   int                  NNZ                       = args.scalar<int>("NNZ");     //number on non-zero entries on sparsity pattern
//   int                  numDOFs                   = args.scalar<int>("numDOFs"); //number of DOFs
//   int                  num_fct_iter              = args.scalar<int>("num_fct_iter");
//   double               dt                        = args.scalar<double>("dt");
//   xt::pyarray<int>    &csrRowIndeces_DofLoops    = args.array<int>("csrRowIndeces_DofLoops");    //csr row indeces
//   xt::pyarray<int>    &csrColumnOffsets_DofLoops = args.array<int>("csrColumnOffsets_DofLoops"); //csr column offsets
//   int                  LUMPED_MASS_MATRIX        = args.scalar<int>("LUMPED_MASS_MATRIX");
//   int                  MONOLITHIC                = args.scalar<int>("MONOLITHIC");

//   // ============================
//   // ========= WATER ============
//   // ============================
//   xt::pyarray<double> &lumped_mass_matrix_water = args.array<double>("lumped_mass_matrix_water");
//   xt::pyarray<double> &soln_water               = args.array<double>("soln_water");
//   xt::pyarray<double> &pn_water                 = args.array<double>("pn_water");
//   xt::pyarray<double> &solH_water               = args.array<double>("solH_water");
//   xt::pyarray<double> &uLow_water               = args.array<double>("uLow_water");
//   xt::pyarray<double> &uDotLow_water            = args.array<double>("uDotLow_water");
//   xt::pyarray<double> &dLow_water               = args.array<double>("dLow_water");
//   xt::pyarray<double> &solLim_water             = args.array<double>("limited_solution_water");
//   xt::pyarray<double> &MC_water                 = args.array<double>("MC_water");
//   xt::pyarray<double> &ML_water                 = args.array<double>("ML_water");
//   xt::pyarray<double> &FluxMatrix_water         = args.array<double>("FluxMatrix_water");
//   xt::pyarray<double> &limitedFlux_water        = args.array<double>("limited_Flux_water");
//   xt::pyarray<double> &MassMatrix_water         = args.array<double>("MassMatrix_water");
//   xt::pyarray<double> &dt_times_fH_minus_fL_water = args.array<double>("dt_times_fH_minus_fL_water");
//   xt::pyarray<double> &min_m_bc_water           = args.array<double>("min_m_bc_water");
//   xt::pyarray<double> &max_m_bc_water           = args.array<double>("max_m_bc_water");

//   double Rpos_water[numDOFs], Rneg_water[numDOFs];
//   int    ij = 0;

//   //////////////////////////////////////////////////////
//   // ********** COMPUTE LOW ORDER SOLUTION ********** //
//   //////////////////////////////////////////////////////
//   if (num_fct_iter == 0) { // No FCT for global bounds
//     for (int i = 0; i < numDOFs; i++) { solLim_water.data()[i] = uLow_water.data()[i]; }
//   } else { // do FCT iterations (with global bounds) on low order solution
//     for (int iter = 0; iter < num_fct_iter; iter++) {
//       ij = 0;
//       for (int i = 0; i < numDOFs; i++) {
//         double maxi = 1.0, Pposi = 0;
//         for (int offset = csrRowIndeces_DofLoops.data()[i];
//              offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//           int j = csrColumnOffsets_DofLoops.data()[offset];
//           // compute Flux correction
//           double Fluxij = FluxMatrix_water.data()[ij] - limitedFlux_water.data()[ij];
//           Pposi += Fluxij * ((Fluxij > 0) ? 1. : 0.);
//           // update ij
//           ij += 1;
//         }
//         // compute Q vectors
//         double mi      = ML_water.data()[i];
//         double solLimi = solLim_water.data()[i];
//         double Qposi   = mi * (maxi - solLimi);
//         // compute R vectors
//         Rpos_water[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
//       }
//       ij = 0;
//       for (int i = 0; i < numDOFs; i++) {
//         double ith_Limiter_times_FluxCorrectionMatrix = 0.;
//         double Rposi                                  = Rpos_water[i];
//         for (int offset = csrRowIndeces_DofLoops.data()[i];
//              offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//           int j = csrColumnOffsets_DofLoops.data()[offset];
//           // Flux Correction
//           double Fluxij = FluxMatrix_water.data()[ij] - limitedFlux_water.data()[ij];
//           // compute limiter
//           double Lij = 1.0;
//           Lij        = (Fluxij > 0 ? Rposi : Rpos_water[j]);
//           // compute limited flux
//           ith_Limiter_times_FluxCorrectionMatrix += Lij * Fluxij;

//           // update limited flux
//           limitedFlux_water.data()[ij] = Lij * Fluxij;

//           // update FluxMatrix
//           FluxMatrix_water.data()[ij] = Fluxij;

//           // update ij
//           ij += 1;
//         }
//         // update limited solution (same structure)
//         double mi = ML_water.data()[i];
//         // (no additional op here in your original loop body)
//       }
//     }
//   }

//   // ***************************************** //
//   // ********** HIGH ORDER SOLUTION ********** //
//   // ***************************************** //
//   ij = 0;
//   for (int i = 0; i < numDOFs; i++) {
//     double mini = soln_water.data()[i], maxi = soln_water.data()[i];
//     double Pposi = 0, Pnegi = 0.;
//     for (int offset = csrRowIndeces_DofLoops.data()[i];
//          offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//       int j = csrColumnOffsets_DofLoops.data()[offset];
//       // compute local bounds //
//       mini = fmin(mini, soln_water.data()[j]);
//       maxi = fmax(maxi, soln_water.data()[j]);
//       // compute P vectors //
//       double fij = (MC_water.data()[ij] * (uDotLow_water.data()[i] - uDotLow_water.data()[j]) / dt
//                     + dLow_water.data()[ij] * (uLow_water.data()[i] - uLow_water.data()[j]));
//       Pposi += fij * (fij > 0 ? 1. : 0.);
//       Pnegi += fij * (fij < 0 ? 1. : 0.);
//       // update ij
//       ij += 1;
//     }
//     // compute Q vectors //
//     double mi    = ML_water.data()[i];
//     double Qposi = mi * (maxi - solLim_water.data()[i]);
//     double Qnegi = mi * (mini - solLim_water.data()[i]);
//     // compute R vectors //
//     Rpos_water[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
//     Rneg_water[i] = ((Pnegi == 0) ? 1. : fmin(1.0, Qnegi / Pnegi));
//   }

//   // COMPUTE LIMITERS //
//   ij = 0;
//   for (int i = 0; i < numDOFs; i++) {
//     double ith_limited_flux_correction = 0;
//     double Rposi                       = Rpos_water[i];
//     double Rnegi                       = Rneg_water[i];
//     for (int offset = csrRowIndeces_DofLoops.data()[i];
//          offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//       int j = csrColumnOffsets_DofLoops.data()[offset];
//       // compute flux correction
//       double fij = (MC_water.data()[ij] * (uDotLow_water.data()[i] - uDotLow_water.data()[j]) / dt
//                     + dLow_water.data()[ij] * (uLow_water.data()[i] - uLow_water.data()[j]));

//       // compute limiters
//       double Lij = 1.0;
//       Lij        = fij > 0 ? fmin(Rposi, Rneg_water[j]) : fmin(Rnegi, Rpos_water[j]);
//       // compute ith_limited_flux_correction
//       ith_limited_flux_correction += Lij * fij;
//       ij += 1;
//     }
//     double mi = ML_water.data()[i];
//     solLim_water[i] += 1. / mi * ith_limited_flux_correction;
//   }

//   // ============================
//   // =========== AIR ============
//   // ============================
//   xt::pyarray<double> &lumped_mass_matrix_air = args.array<double>("lumped_mass_matrix_air");
//   xt::pyarray<double> &soln_air               = args.array<double>("soln_air");
//   xt::pyarray<double> &pn_air                 = args.array<double>("pn_air");
//   xt::pyarray<double> &solH_air               = args.array<double>("solH_air");
//   xt::pyarray<double> &uLow_air               = args.array<double>("uLow_air");
//   xt::pyarray<double> &uDotLow_air            = args.array<double>("uDotLow_air");
//   xt::pyarray<double> &dLow_air               = args.array<double>("dLow_air");
//   xt::pyarray<double> &solLim_air             = args.array<double>("limited_solution_air");
//   xt::pyarray<double> &MC_air                 = args.array<double>("MC_air");
//   xt::pyarray<double> &ML_air                 = args.array<double>("ML_air");
//   xt::pyarray<double> &FluxMatrix_air         = args.array<double>("FluxMatrix_air");
//   xt::pyarray<double> &limitedFlux_air        = args.array<double>("limited_Flux_air");
//   xt::pyarray<double> &MassMatrix_air         = args.array<double>("MassMatrix_air");
//   xt::pyarray<double> &dt_times_fH_minus_fL_air = args.array<double>("dt_times_fH_minus_fL_air");
//   xt::pyarray<double> &min_m_bc_air           = args.array<double>("min_m_bc_air");
//   xt::pyarray<double> &max_m_bc_air           = args.array<double>("max_m_bc_air");

//   double Rpos_air[numDOFs], Rneg_air[numDOFs];
//   // reuse ij var; reset before use
//   ij = 0;

//   //////////////////////////////////////////////////////
//   // ********** COMPUTE LOW ORDER SOLUTION ********** //
//   //////////////////////////////////////////////////////
//   if (num_fct_iter == 0) { // No FCT for global bounds
//     for (int i = 0; i < numDOFs; i++) { solLim_air.data()[i] = uLow_air.data()[i]; }
//   } else { // do FCT iterations (with global bounds) on low order solution
//     for (int iter = 0; iter < num_fct_iter; iter++) {
//       ij = 0;
//       for (int i = 0; i < numDOFs; i++) {
//         double maxi = 1.0, Pposi = 0;
//         for (int offset = csrRowIndeces_DofLoops.data()[i];
//              offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//           int j = csrColumnOffsets_DofLoops.data()[offset];
//           // compute Flux correction
//           double Fluxij = FluxMatrix_air.data()[ij] - limitedFlux_air.data()[ij];
//           Pposi += Fluxij * ((Fluxij > 0) ? 1. : 0.);
//           // update ij
//           ij += 1;
//         }
//         // compute Q vectors
//         double mi      = ML_air.data()[i];
//         double solLimi = solLim_air.data()[i];
//         double Qposi   = mi * (maxi - solLimi);
//         // compute R vectors
//         Rpos_air[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
//       }
//       ij = 0;
//       for (int i = 0; i < numDOFs; i++) {
//         double ith_Limiter_times_FluxCorrectionMatrix = 0.;
//         double Rposi                                  = Rpos_air[i];
//         for (int offset = csrRowIndeces_DofLoops.data()[i];
//              offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//           int j = csrColumnOffsets_DofLoops.data()[offset];
//           // Flux Correction
//           double Fluxij = FluxMatrix_air.data()[ij] - limitedFlux_air.data()[ij];
//           // compute limiter
//           double Lij = 1.0;
//           Lij        = (Fluxij > 0 ? Rposi : Rpos_air[j]);
//           // compute limited flux
//           ith_Limiter_times_FluxCorrectionMatrix += Lij * Fluxij;

//           // update limited flux
//           limitedFlux_air.data()[ij] = Lij * Fluxij;

//           // update FluxMatrix
//           FluxMatrix_air.data()[ij] = Fluxij;

//           // update ij
//           ij += 1;
//         }
//         // update limited solution (same structure)
//         double mi = ML_air.data()[i];
//         // (no additional op here in your original loop body)
//       }
//     }
//   }

//   // ***************************************** //
//   // ********** HIGH ORDER SOLUTION ********** //
//   // ***************************************** //
//   ij = 0;
//   for (int i = 0; i < numDOFs; i++) {
//     double mini = soln_air.data()[i], maxi = soln_air.data()[i];
//     double Pposi = 0, Pnegi = 0.;
//     for (int offset = csrRowIndeces_DofLoops.data()[i];
//          offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//       int j = csrColumnOffsets_DofLoops.data()[offset];
//       // compute local bounds //
//       mini = fmin(mini, soln_air.data()[j]);
//       maxi = fmax(maxi, soln_air.data()[j]);
//       // compute P vectors //
//       double fij = (MC_air.data()[ij] * (uDotLow_air.data()[i] - uDotLow_air.data()[j]) / dt
//                     + dLow_air.data()[ij] * (uLow_air.data()[i] - uLow_air.data()[j]));
//       Pposi += fij * (fij > 0 ? 1. : 0.);
//       Pnegi += fij * (fij < 0 ? 1. : 0.);
//       // update ij
//       ij += 1;
//     }
//     // compute Q vectors //
//     double mi    = ML_air.data()[i];
//     double Qposi = mi * (maxi - solLim_air.data()[i]);
//     double Qnegi = mi * (mini - solLim_air.data()[i]);
//     // compute R vectors //
//     Rpos_air[i] = ((Pposi == 0) ? 1. : fmin(1.0, Qposi / Pposi));
//     Rneg_air[i] = ((Pnegi == 0) ? 1. : fmin(1.0, Qnegi / Pnegi));
//   }

//   // COMPUTE LIMITERS //
//   ij = 0;
//   for (int i = 0; i < numDOFs; i++) {
//     double ith_limited_flux_correction = 0;
//     double Rposi                       = Rpos_air[i];
//     double Rnegi                       = Rneg_air[i];
//     for (int offset = csrRowIndeces_DofLoops.data()[i];
//          offset < csrRowIndeces_DofLoops.data()[i + 1]; offset++) {
//       int j = csrColumnOffsets_DofLoops.data()[offset];
//       // compute flux correction
//       double fij = (MC_air.data()[ij] * (uDotLow_air.data()[i] - uDotLow_air.data()[j]) / dt
//                     + dLow_air.data()[ij] * (uLow_air.data()[i] - uLow_air.data()[j]));

//       // compute limiters
//       double Lij = 1.0;
//       Lij        = fij > 0 ? fmin(Rposi, Rneg_air[j]) : fmin(Rnegi, Rpos_air[j]);
//       // compute ith_limited_flux_correction
//       ith_limited_flux_correction += Lij * fij;
//       ij += 1;
//     }
//     double mi = ML_air.data()[i];
//     solLim_air[i] += 1. / mi * ith_limited_flux_correction;
//   }
// }

  void calculateResidual_entropy_viscosity(arguments_dict &args)
  {
    xt::pyarray<double> &globalJacobian = args.array<double>("globalJacobian");
  //  xt::pyarray<double> &globalJacobian   = args.array<double>("globalJacobian");
    
    double               Theta                     = args.scalar<double>("Theta");
    xt::pyarray<double> &bc_mask_water             = args.array<double>("bc_mask_water");
    xt::pyarray<double> &bc_mask_air               = args.array<double>("bc_mask_air");
    
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
    xt::pyarray<double> &ebqe_phi_water                                   = args.array<double>("ebqe_phi_water");
    xt::pyarray<double> &ebqe_phi_air                                   = args.array<double>("ebqe_phi_air");
    
    double               epsFact                                    = args.scalar<double>("epsFact");
    xt::pyarray<double> &cfl                                        = args.array<double>("cfl");
    



    /////////////////////////Two-phase/////////////////////////////////////////
    double               rho_water                                  = args.scalar<double>("rho_water");
    double               rho_air                                    = args.scalar<double>("rho_air");  
    double               beta_water                                 = args.scalar<double>("beta_water");
    double               beta_air                                   = args.scalar<double>("beta_air");    
    double Y_air                                                    = args.scalar<double>("Y_air");     // e.g. 1.0
    double Y_water                                                  = args.scalar<double>("Y_water");   // e.g. 1.0
    
    //double               rho                                        = args.scalar<double>("rho");
    //double               beta                                       = args.scalar<double>("beta");
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
    xt::pyarray<int>    &u_l2g_water                                = args.array<int>("u_l2g_water");
    xt::pyarray<int>    &r_l2g_water                                = args.array<int>("r_l2g_water");

    xt::pyarray<int>    &u_l2g_air                                  = args.array<int>("u_l2g_air");
    xt::pyarray<int>    &r_l2g_air                                  = args.array<int>("r_l2g_air");

    xt::pyarray<double> &elementDiameter                            = args.array<double>("elementDiameter");
    int                  degree_polynomial                          = args.scalar<int>("degree_polynomial");
    int                  nExteriorElementBoundaries_global          = args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int>    &exteriorElementBoundariesArray             = args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int>    &elementBoundaryElementsArray               = args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int>    &elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
    
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
    xt::pyarray<double> &q_dV                                       = args.array<double>("q_dV");
    
    // xt::pyarray<double> &q_dV_water                                 = args.array<double>("q_dV_water");
    // xt::pyarray<double> &q_dV_air                                   = args.array<double>("q_dV_air");
    // xt::pyarray<double> &q_m_betaBDF_water                          = args.array<double>("q_m_betaBDF_water");
    xt::pyarray<double> &q_m_betaBDF_air                            = args.array<double>("q_m_betaBDF_air");
    xt::pyarray<double> &q_numDiff_u_water                          = args.array<double>("q_numDiff_u_water");    
    xt::pyarray<double> &q_numDiff_u_air                            = args.array<double>("q_numDiff_u_air");
    xt::pyarray<double> &q_numDiff_u_last_water                     = args.array<double>("q_numDiff_u_last_water");
    xt::pyarray<double> &q_numDiff_u_last_air                       = args.array<double>("q_numDiff_u_last_air");
    //int                  offset_u                                   = args.scalar<int>("offset_u");
    //int                  stride_u                                   = args.scalar<int>("stride_u");
    
    int                  offset_u_water                             = args.scalar<int>("offset_u_water");
    int                  offset_u_air                               = args.scalar<int>("offset_u_air");
    int                  stride_u_water                             = args.scalar<int>("stride_u_water");
    int                  stride_u_air                               = args.scalar<int>("stride_u_air");
    
    
    xt::pyarray<double> &globalResidual                       = args.array<double>("globalResidual");
    //xt::pyarray<double> &globalResidual                         = args.array<double>("globalResidual");
    
    
    
    xt::pyarray<double> &ebqe_velocity_ext_water                    = args.array<double>("ebqe_velocity_ext_water");
    xt::pyarray<double> &ebqe_velocity_ext_air                      = args.array<double>("ebqe_velocity_ext_air");
    xt::pyarray<int>    &isDOFBoundary_u_water                      = args.array<int>("isDOFBoundary_u_water");
    xt::pyarray<int>    &isDOFBoundary_u_air                        = args.array<int>("isDOFBoundary_u_air");
    xt::pyarray<double> &ebqe_bc_u_ext_water                        = args.array<double>("ebqe_bc_u_ext_water");
    xt::pyarray<double> &ebqe_bc_u_ext_air                          = args.array<double>("ebqe_bc_u_ext_air");
    
    
    // xt::pyarray<double> &globalResidual                             = args.array<double>("globalResidual");
    
    // xt::pyarray<double> &ebqe_velocity_ext                          = args.array<double>("ebqe_velocity_ext");
    // xt::pyarray<int>    &isDOFBoundary_u                            = args.array<int>("isDOFBoundary_u");
    // xt::pyarray<double> &ebqe_bc_u_ext                              = args.array<double>("ebqe_bc_u_ext");
    // xt::pyarray<int>    &isFluxBoundary_u                           = args.array<int>("isFluxBoundary_u");
    // xt::pyarray<double> &ebqe_bc_flux_ext                           = args.array<double>("ebqe_bc_flux_ext");
    
    xt::pyarray<double> &ebqe_u_water                               = args.array<double>("ebqe_u_water");
    xt::pyarray<double> &ebqe_u_air                                 = args.array<double>("ebqe_u_air");
    xt::pyarray<double> &ebqe_flux_water                            = args.array<double>("ebqe_flux_water");
    xt::pyarray<double> &ebqe_flux_air                              = args.array<double>("ebqe_flux_air");

    //xt::pyarray<double> &ebqe_u                                     = args.array<double>("ebqe_u");
    //xt::pyarray<double> &ebqe_flux                                  = args.array<double>("ebqe_flux");
    // PARAMETERS FOR EDGE BASED STABILIZATION
    double cE = args.scalar<double>("cE");
    double cK = args.scalar<double>("cK");
    // PARAMETERS FOR LOG BASED ENTROPY FUNCTION
    double uL_water = args.scalar<double>("uL_water");
    double uL_air = args.scalar<double>("uL_air");
    double uR_water = args.scalar<double>("uR_water");
    double uR_air = args.scalar<double>("uR_air");
    // PARAMETERS FOR EDGE VISCOSITY
    int               numDOFs                       = args.scalar<int>("numDOFs");
    int               NNZ                           = args.scalar<int>("NNZ");
    xt::pyarray<int> &csrRowIndeces_DofLoops_water        = args.array<int>("csrRowIndeces_DofLoops_water");
    xt::pyarray<int> &csrColumnOffsets_DofLoops_water     = args.array<int>("csrColumnOffsets_DofLoops_water");

    xt::pyarray<int> &csrRowIndeces_DofLoops_air        = args.array<int>("csrRowIndeces_DofLoops_air");
    xt::pyarray<int> &csrColumnOffsets_DofLoops_air     = args.array<int>("csrColumnOffsets_DofLoops_air");
    
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

    xt::pyarray<double> &ML_water  = args.array<double>("ML_water");
    xt::pyarray<double> &ML_air  = args.array<double>("ML_air");
    xt::pyarray<double> &MC_water  = args.array<double>("MC_water");
    xt::pyarray<double> &MC_air  = args.array<double>("MC_air");
    xt::pyarray<double> &delta_x_ij_water = args.array<double>("delta_x_ij_water");
    xt::pyarray<double> &delta_x_ij_air = args.array<double>("delta_x_ij_air");
    // PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
    int LUMPED_MASS_MATRIX = args.scalar<int>("LUMPED_MASS_MATRIX");
    STABILIZATION STABILIZATION_TYPE{static_cast<STABILIZATION>(args.scalar<int>("STABILIZATION_TYPE"))};
    
    
    //////////////////////////////////////For Brooks- Corey/////////////////
    PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
    // double BC_entry_head = args.scalar<double>("BC_entry_head");
    // double BC_lambda     = args.scalar<double>("BC_lambda");
    xt::pyarray<double> &BC_entry_head = args.array<double>("BC_entry_head"); // size: nMaterials
    xt::pyarray<double> &BC_lambda     = args.array<double>("BC_lambda");     // size: nMaterials



    int ENTROPY_TYPE = args.scalar<int>("ENTROPY_TYPE");
    // FOR FCT
    xt::pyarray<double> &dLow_water            = args.array<double>("dLow_water");
    xt::pyarray<double> &dLow_air                 = args.array<double>("dLow_air");

    xt::pyarray<double> &fluxMatrix_water           = args.array<double>("fluxMatrix_water");
    xt::pyarray<double> &fluxMatrix_air           = args.array<double>("fluxMatrix_air");

    xt::pyarray<double> &mDotLow_water         = args.array<double>("mDotLow_water");
    xt::pyarray<double> &mDotHigh_water              = args.array<double>("mDotHigh_water");
    xt::pyarray<double> &mLow_water                 = args.array<double>("mLow_water");
    xt::pyarray<double> &dt_times_fH_minus_fL_water = args.array<double>("dt_times_fH_minus_fL_water");
    xt::pyarray<double> &min_m_bc_water             = args.array<double>("min_m_bc_water");
    xt::pyarray<double> &max_m_bc_water             = args.array<double>("max_m_bc_water");

    xt::pyarray<double> &mDotLow_air         = args.array<double>("mDotLow_air");
    xt::pyarray<double> &mDotHigh_air              = args.array<double>("mDotHigh_air");
    xt::pyarray<double> &mLow_air                 = args.array<double>("mLow_air");
    xt::pyarray<double> &dt_times_fH_minus_fL_air = args.array<double>("dt_times_fH_minus_fL_air");
    xt::pyarray<double> &min_m_bc_air             = args.array<double>("min_m_bc_air");
    xt::pyarray<double> &max_m_bc_air             = args.array<double>("max_m_bc_air");


    // AUX QUANTITIES OF INTEREST
    xt::pyarray<double> &quantDOFs_water = args.array<double>("quantDOFs_water");
    xt::pyarray<double> &quantDOFs_air = args.array<double>("quantDOFs_air");
    
    xt::pyarray<double> &mn        = args.array<double>("mn");

    xt::pyarray<double> &fluxCorrection_water        = args.array<double>("fluxCorrection_water");
    xt::pyarray<double> &fluxCorrection_air        = args.array<double>("fluxCorrection_air");
    
    xt::pyarray<double> &limited_solution_water          = args.array<double>("limited_solution_water");
    xt::pyarray<double> &limited_solution_air          = args.array<double>("limited_solution_air");

    xt::pyarray<double> &anb_seepage_flux_n = args.array<double>("anb_seepage_flux_n");
    xt::pyarray<double> &q_velocity_water = args.array<double>("q_velocity_water");
    xt::pyarray<double> &q_velocity_air = args.array<double>("q_velocity_air");
    
    double &anb_seepage_flux(args.scalar<double>("anb_seepage_flux"));
    anb_seepage_flux = 0.0;
    xt::pyarray<int>    &csrRowIndeces_u_u                          = args.array<int>("csrRowIndeces_u_u");
    xt::pyarray<int>    &csrColumnOffsets_u_u                       = args.array<int>("csrColumnOffsets_u_u");
    xt::pyarray<int>    &csrColumnOffsets_eb_u_u                    = args.array<int>("csrColumnOffsets_eb_u_u");

    xt::pyarray<double> &ebqe_bc_flux_ext_water                     = args.array<double>("ebqe_bc_flux_ext_water");
    xt::pyarray<double> &ebqe_bc_flux_ext_air                       = args.array<double>("ebqe_bc_flux_ext_air");
    
    double Rpos[numDOFs], Rneg[numDOFs];
    //double FluxCorrectionMatrix[NNZ];
    // NOTE: This function follows a different (but equivalent) implementation of the smoothness based indicator than NCLS.h
    // Allocate space for the transport matrices
    // This is used for first order KUZMIN'S METHOD
    double                TransportMatrix_water[NNZ], TransportMatrixConsistent_water[NNZ];
    double                TransportMatrixn_water[NNZ], TransportMatrixConsistentn_water[NNZ];
    std::valarray<double> u_free_dof_water(numDOFs);
    std::valarray<double> u_free_dof_old_water(numDOFs);
    std::valarray<double> ML2_water(numDOFs);

    double                TransportMatrix_air[NNZ], TransportMatrixConsistent_air[NNZ];
    double                TransportMatrixn_air[NNZ], TransportMatrixConsistentn_air[NNZ];
    std::valarray<double> u_free_dof_air(numDOFs);
    std::valarray<double> u_free_dof_old_air(numDOFs);
    std::valarray<double> ML2_air(numDOFs);

    //phase :0 and phase 1
    for (int eN = 0; eN < nElements_global; eN++)
      for (int j = 0; j < nDOF_trial_element; j++) {
        int eN_nDOF_trial_element                               = eN * nDOF_trial_element;
        u_free_dof_water[r_l2g_water.data()[eN_nDOF_trial_element + j]]  = u_dof_water.data()[u_l2g_water.data()[eN_nDOF_trial_element + j]];
        u_free_dof_air[r_l2g_air.data()[eN_nDOF_trial_element + j]]    = u_dof_air.data()[u_l2g_air.data()[eN_nDOF_trial_element + j]];
        
        u_free_dof_old_water[r_l2g_water.data()[eN_nDOF_trial_element + j]] = u_dof_old_water.data()[u_l2g_water.data()[eN_nDOF_trial_element + j]];        
        u_free_dof_old_air[r_l2g_air.data()[eN_nDOF_trial_element + j]] = u_dof_old_air.data()[u_l2g_air.data()[eN_nDOF_trial_element + j]];
      }
    for (int i = 0; i < NNZ; i++) {
      TransportMatrix_water[i]            = 0.;
      TransportMatrixConsistent_water[i]  = 0.;
      TransportMatrixn_water[i]           = 0.;
      TransportMatrixConsistentn_water[i] = 0.;
    }

    for (int i = 0; i < NNZ; i++) {
      TransportMatrix_air[i]            = 0.;
      TransportMatrixConsistent_air[i]  = 0.;
      TransportMatrixn_air[i]           = 0.;
      TransportMatrixConsistentn_air[i] = 0.;
    }

    // compute entropy and init global_entropy_residual and boundary_integral
    double psi_water[numDOFs], eta_water[numDOFs], global_entropy_residual_water[numDOFs], boundary_integral_water[numDOFs];
    double psi_air[numDOFs], eta_air[numDOFs], global_entropy_residual_air[numDOFs], boundary_integral_air[numDOFs];
    
    for (int i = 0; i < numDOFs; i++) {
      // NODAL ENTROPY //
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV stab
      {
        double solni_water = 1.0 * u_free_dof_old_water[i];
        double solni_air  = 1.0 * u_free_dof_old_air[i];
        eta_water[i]                      = ENTROPY_TYPE == 1 ? ENTROPY(solni_water, uL_water, uR_water) : ENTROPY_LOG(solni_water, uL_water, uR_water);
        eta_air[i]                        = ENTROPY_TYPE == 1 ? ENTROPY(solni_air, uL_air, uR_air) : ENTROPY_LOG(solni_air, uL_air, uR_air);

        global_entropy_residual_water[i]  = 0.;
        global_entropy_residual_air[i]  = 0.;
      }
      boundary_integral_water[i] = 0.;
      boundary_integral_air[i] = 0.;
      
      ML2_water[i]               = 0.0;
      ML2_air[i]               = 0.0;
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
      //phase 0::water
      double elementResidual_u_water[nDOF_test_element], element_entropy_residual_water[nDOF_test_element], Phi_water[nDOF_trial_element], Phi_n_water[nDOF_trial_element];
      double elementTransport_water[nDOF_test_element][nDOF_trial_element], elementTransportConsistent_water[nDOF_test_element][nDOF_trial_element];
      double elementTransportn_water[nDOF_test_element][nDOF_trial_element], elementTransportConsistentn_water[nDOF_test_element][nDOF_trial_element];
      //phase 1:: air
      double elementResidual_u_air[nDOF_test_element], element_entropy_residual_air[nDOF_test_element], Phi_air[nDOF_trial_element], Phi_n_air[nDOF_trial_element];
      double elementTransport_air[nDOF_test_element][nDOF_trial_element], elementTransportConsistent_air[nDOF_test_element][nDOF_trial_element];
      double elementTransportn_air[nDOF_test_element][nDOF_trial_element], elementTransportConsistentn_air[nDOF_test_element][nDOF_trial_element];
      //Later: Need to look for mesh_dof
      //phase 0::water
      for (int i = 0; i < nDOF_test_element; i++) {
        Phi_water[i]   = u_dof_water[i];
        Phi_n_water[i] = u_dof_old_water[i];
        Phi_air[i]   = u_dof_air[i];
        Phi_n_air[i] = u_dof_old_air[i];
        for (int I = 0; I < nSpace; I++) {
          Phi_water[i] -= rho_water * mesh_dof[i * 3 + I] * gravity[I];
          Phi_n_water[i] -= rho_water * mesh_dof[i * 3 + I] * gravity[I];
          Phi_air[i] -= rho_air * mesh_dof[i * 3 + I] * gravity[I];
          Phi_n_air[i] -= rho_air * mesh_dof[i * 3 + I] * gravity[I];
        }

        elementResidual_u_water[i] = 0.0, element_entropy_residual_water[i] = 0.0;
        elementResidual_u_air[i] = 0.0, element_entropy_residual_air[i] = 0.0;

        for (int j = 0; j < nDOF_trial_element; j++) {
         //phase 0:: water
          elementTransport_water[i][j]            = 0.0;
          elementTransportConsistent_water[i][j]  = 0.0;
          elementTransportn_water[i][j]           = 0.0;
          elementTransportConsistentn_water[i][j] = 0.0;
          //phase 1:: air
          elementTransport_air[i][j]            = 0.0;
          elementTransportConsistent_air[i][j]  = 0.0;
          elementTransportn_air[i][j]           = 0.0;
          elementTransportConsistentn_air[i][j] = 0.0;
        }
      }
      //loop over quadrature points and compute integrands
      for (int k = 0; k < nQuadraturePoints_element; k++) {
        //compute indeces and declare local storage
        int eN_k = eN * nQuadraturePoints_element + k, eN_k_nSpace = eN_k * nSpace, eN_nDOF_trial_element = eN * nDOF_trial_element;
        double
          // for entropy residual phase 0::water
          aux_entropy_residual_water = 0.,
          DENTROPY_un_water, DENTROPY_uni_water,
          // for entropy residual phase 1::air
          aux_entropy_residual_air = 0.,
          DENTROPY_un_air, DENTROPY_uni_air,
          u_test_dV[nDOF_trial_element], u_grad_trial[nDOF_trial_element * nSpace], u_grad_test_dV[nDOF_test_element * nSpace],
          //for mass matrix contributions water
          u_water = 0.0, un_water = 0.0, grad_u_water[nSpace], grad_un_water[nSpace], velocity_loc_water[nSpace],
          u_air = 0.0, un_air = 0.0, grad_u_air[nSpace], grad_un_air[nSpace], velocity_loc_air[nSpace],
          //for general use
          jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace], dV, x, y, z, xt, yt, zt, 
          m_water, dm_water, f_water[nSpace], df_water[nSpace], a_water[nnz], da_water[nnz], as_water[nnz], mn_water, dmn_water, fn_water[nSpace], dfn_water[nSpace], an_water[nnz], dan_water[nnz], 
          asn_water[nnz],
          m_air, dm_air, f_air[nSpace], df_air[nSpace], a_air[nnz], da_air[nnz], as_air[nnz], mn_air, dmn_air, fn_air[nSpace], dfn_air[nSpace], an_air[nnz], dan_air[nnz], 
          asn_air[nnz];
        //get the physical integration weight
        ck.calculateMapping_element(eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y, z);
        ck.calculateMappingVelocity_element(eN, k, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), xt, yt, zt);
        dV = fabs(jacDet) * dV_ref.data()[k];

        //get the solution (of Newton's solver). To compute time derivative term
        ck.valFromDOF(u_dof_water.data(), &u_l2g_water.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], u_water);
        ck.valFromDOF(u_dof_air.data(), &u_l2g_air.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], u_air);
        
        
        //get the solution at quad point at tn and tnm1 for entropy viscosity
        ck.valFromDOF(u_dof_old_water.data(), &u_l2g_water.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], un_water);
        ck.valFromDOF(u_dof_old_air.data(), &u_l2g_air.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], un_air);
        

        //get the solution gradients at tn for entropy viscosity
        ck.gradTrialFromRef(&u_grad_trial_ref.data()[k * nDOF_trial_element * nSpace], jacInv, u_grad_trial);
        //precalculate test function products with integration weights for mass matrix terms
        for (int I = 0; I < nSpace; I++) {
          grad_u_water[I]  = 0.0;
          grad_un_water[I] = 0.0;
          grad_u_air[I]  = 0.0;
          grad_un_air[I] = 0.0;          
        }

        for (int j = 0; j < nDOF_trial_element; j++) {
          u_test_dV[j] = u_test_ref.data()[k * nDOF_trial_element + j] * dV;
          for (int I = 0; I < nSpace; I++) {
            grad_un_water[I] += Phi_n_water[j] * u_grad_trial[j * nSpace + I];//note: grad u is grad phi
            grad_u_water[I] += Phi_water[j] * u_grad_trial[j * nSpace + I];

            grad_un_air[I] += Phi_n_air[j] * u_grad_trial[j * nSpace + I];//note: grad u is grad phi
            grad_u_air[I] += Phi_air[j] * u_grad_trial[j * nSpace + I];

            u_grad_test_dV[j * nSpace + I] = u_grad_trial[j * nSpace + I] * dV; //cek warning won't work for Petrov-Galerkin
          }
        }
        //
        //calculate pde coefficients at quadrature points
        //
        double Kr_water, dKr_water, Krn_water, dKrn_water;
        double Kr_air, dKr_air, Krn_air, dKrn_air;
        double Sw=0.0, Sg=0.0;


        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                            rho_water, rho_air, 
                            beta_water, beta_air, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            un_water, un_air, 
                            Y_water, Y_air,
                            mn_water, mn_air, 
                            dmn_water, dmn_air,
                            fn_water, fn_air,
                            dfn_water, dfn_air,
                            an_water, an_air,
                            dan_water, dan_air, 
                            asn_water, asn_air,
                            Krn_water, dKrn_water,
                            Krn_air, dKrn_air,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw, Sg,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));

        
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                            rho_water, rho_air, 
                            beta_water, beta_air, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            u_water, u_air, 
                            Y_water, Y_air,
                            m_water, m_air, 
                            dm_water, dm_air,
                            f_water, f_air,
                            df_water, df_air,
                            a_water, a_air,
                            da_water, da_air, 
                            as_water, as_air,
                            Kr_water, 
                            dKr_water,Kr_air, dKr_air,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw, Sg,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));

        // Darcy velocity calculation
        for (int I = 0; I < nSpace; I++) { 
          velocity_loc_water[I] = 0.0;
          velocity_loc_air[I] = 0.0;
         }
        for (int I = 0; I < nSpace; I++) {
          for (int J = 0; J < nSpace; J++) { 
            velocity_loc_water[I] -= Kr_water * KWs.data()[elementMaterialTypes[eN] * nSpace * nSpace + I * nSpace + J] * grad_u_water[J];
            velocity_loc_air[I] -= Kr_air * KWs.data()[elementMaterialTypes[eN] * nSpace * nSpace + I * nSpace + J] * grad_u_air[J]; 
          }
        }
        for (int I = 0; I < nSpace; I++) { 
          q_velocity_water.data()[eN_k_nSpace + I] = velocity_loc_water[I];
          q_velocity_air.data()[eN_k_nSpace + I] = velocity_loc_air[I];
        
        }

        //
        //moving mesh
        //
        double mesh_velocity[3];
        mesh_velocity[0] = xt;
        mesh_velocity[1] = yt;
        mesh_velocity[2] = zt;
        //relative velocity at tn
        for (int I = 0; I < nSpace; I++) {
          f_water[I] -= MOVING_DOMAIN * m_water * mesh_velocity[I];
          f_air[I] -= MOVING_DOMAIN * m_air * mesh_velocity[I];        
          velocity_loc_water[I] = df_water[I] * (2.0 * dm_water * dm_water / (dm_water * dm_water + fmax(1.0e-16, dm_water * dm_water)));
          velocity_loc_air[I] = df_air[I] * (2.0 * dm_air * dm_air / (dm_air * dm_air + fmax(1.0e-16, dm_air * dm_air)));
        }
        //////////////////////////////
        // CALCULATE CELL BASED CFL //
        //////////////////////////////
        //////////////////////////////
      
       const double h_eff = elementDiameter.data()[eN] / degree_polynomial;
       double cfl_w = 0.0, cfl_a = 0.0;
       calculateCFL(h_eff, velocity_loc_water, cfl_w);
       calculateCFL(h_eff, velocity_loc_air,  cfl_a);
       cfl.data()[eN_k] = fmax(cfl_w, cfl_a);
      

       // calculateCFL(elementDiameter.data()[eN] / degree_polynomial, velocity_loc, cfl.data()[eN_k]);

        //////////////////////////////////////////////
        // CALCULATE ENTROPY RESIDUAL AT QUAD POINT //
        //////////////////////////////////////////////
        if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) // EV stab
        {
          for (int I = 0; I < nSpace; I++) {
            aux_entropy_residual_water += velocity_loc_water[I] * grad_un_water[I];
            aux_entropy_residual_air += velocity_loc_air[I] * grad_un_air[I];
          }
          DENTROPY_un_water = ENTROPY_TYPE == 1 ? DENTROPY(un_water, uL_water, uR_water) : DENTROPY_LOG(un_water, uL_water, uR_water);
          DENTROPY_un_air = ENTROPY_TYPE == 1 ? DENTROPY(un_air, uL_air, uR_air) : DENTROPY_LOG(un_air, uL_air, uR_air);        
        }
        //////////////
        // ith-LOOP //
        //////////////
        for (int i = 0; i < nDOF_test_element; i++) {
          // VECTOR OF ENTROPY RESIDUAL //
          int eN_i = eN * nDOF_test_element + i;
          ML2_water[u_l2g_water.data()[eN_i]] += u_test_dV[i];
          ML2_air[u_l2g_air.data()[eN_i]] += u_test_dV[i];
          if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) // EV stab
          {
            int    gi_water = offset_u_water + stride_u_water * u_l2g_water.data()[eN_i]; //global i-th index
            int    gi_air   = offset_u_air + stride_u_air * u_l2g_air.data()[eN_i]; //global i-th index
            
            double uni_water = u_dof_old_water.data()[gi_water];
            double uni_air = u_dof_old_air.data()[gi_air];
            DENTROPY_uni_water              = ENTROPY_TYPE == 1 ? DENTROPY(uni_water, uL_water, uR_water) : DENTROPY_LOG(uni_water, uL_water, uR_water);           
            DENTROPY_uni_air              = ENTROPY_TYPE == 1 ? DENTROPY(uni_air, uL_air, uR_air) : DENTROPY_LOG(uni_air, uL_air, uR_air);
            element_entropy_residual_water[i] += (DENTROPY_un_water - DENTROPY_uni_water) * aux_entropy_residual_water * u_test_dV[i];
            element_entropy_residual_air[i] += (DENTROPY_un_air - DENTROPY_uni_air) * aux_entropy_residual_air * u_test_dV[i];
          }
          elementResidual_u_water[i] += m_water * u_test_dV[i];
          elementResidual_u_air[i] += m_air * u_test_dV[i];
          
          ///////////////
          // j-th LOOP // To construct transport matrices
          ///////////////
          // phase 0:: water
          for (int j = 0; j < nDOF_trial_element; j++) {
            int j_nSpace = j * nSpace;
            int i_nSpace = i * nSpace;
            elementTransport_water[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), as_water, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportConsistent_water[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), a_water, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportn_water[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), asn_water, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportConsistentn_water[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), an_water, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
          }

          for (int j = 0; j < nDOF_trial_element; j++) {
            int j_nSpace = j * nSpace;
            int i_nSpace = i * nSpace;
            elementTransport_air[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), as_air, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportConsistent_air[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), a_air, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportn_air[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), asn_air, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
            elementTransportConsistentn_air[i][j] += ck.SimpleDiffusionJacobian_weak(a_rowptr.data(), a_colind.data(), an_air, &u_grad_trial[j_nSpace], &u_grad_test_dV[i_nSpace]);
          }
        } //i
        //save solution for other models
        q_u_water.data()[eN_k] = u_water;
        q_m_water.data()[eN_k] = m_water;

        q_u_air.data()[eN_k] = u_air;
        q_m_air.data()[eN_k] = m_air;
      }
      /////////////////
      // DISTRIBUTE // load cell based element into global residual
      ////////////////
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        int gi_water   = offset_u_water + stride_u_water * r_l2g_water.data()[eN_i]; //global i-th index
        int gi_air   = offset_u_air + stride_u_air * r_l2g_air.data()[eN_i]; //global i-th index
        // distribute entropy_residual
        if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) // EV Stab
          global_entropy_residual_water[gi_water] += element_entropy_residual_water[i];
          global_entropy_residual_air[gi_air] += element_entropy_residual_air[i];
        // distribute transport matrices
        for (int j = 0; j < nDOF_trial_element; j++) {
          int eN_i_j = eN_i * nDOF_trial_element + j;
        ///phase 0::water
          TransportMatrix_water[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransport_water[i][j];
          TransportMatrixConsistent_water[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportConsistent_water[i][j];
          TransportMatrixn_water[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportn_water[i][j];
          TransportMatrixConsistentn_water[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportConsistentn_water[i][j];
        ///phase 1: air
          TransportMatrix_air[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransport_air[i][j];
          TransportMatrixConsistent_air[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportConsistent_air[i][j];
          TransportMatrixn_air[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportn_air[i][j];
          TransportMatrixConsistentn_air[csrRowIndeces_CellLoops.data()[eN_i] + csrColumnOffsets_CellLoops.data()[eN_i_j]] += elementTransportConsistentn_air[i][j];                
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
      double elementResidual_u_water[nDOF_test_element] , elementResidual_u_air[nDOF_test_element] ;
      for (int i = 0; i < nDOF_test_element; i++) { 
        elementResidual_u_water[i] = 0.0;
        elementResidual_u_air[i] = 0.0;
       }

      for (int kb = 0; kb < nQuadraturePoints_elementBoundary; kb++) {
        int    ebNE_kb = ebNE * nQuadraturePoints_elementBoundary + kb, ebNE_kb_nSpace = ebNE_kb * nSpace, ebN_local_kb = ebN_local * nQuadraturePoints_elementBoundary + kb, ebN_local_kb_nSpace = ebN_local_kb * nSpace;
        double jac_ext[nSpace * nSpace], jacDet_ext, jacInv_ext[nSpace * nSpace], boundaryJac[nSpace * (nSpace - 1)], 
                metricTensor[(nSpace - 1) * (nSpace - 1)], metricTensorDetSqrt, dS, u_test_dS[nDOF_test_element], u_grad_trial_trace[nDOF_trial_element * nSpace], normal[3], 
                x_ext, y_ext, z_ext, xt_ext, yt_ext, zt_ext, integralScaling, G[nSpace * nSpace], G_dd_G, tr_G;
        double u_ext_water = 0.0, grad_u_ext_water[nSpace], m_ext_water = 0.0, dm_ext_water = 0.0, 
               f_ext_water[nSpace], df_ext_water[nSpace], a_ext_water[nnz], da_ext_water[nnz], as_ext_water[nnz], flux_ext_water = 0.0,
               //anb_seepage_flux=0.0, // for flux calculation
               bc_u_ext_water = 0.0, bc_grad_u_ext_water[nSpace], bc_m_ext_water = 0.0, bc_dm_ext_water = 0.0, 
               bc_f_ext_water[nSpace], bc_df_ext_water[nSpace], bc_a_ext_water[nnz], bc_da_ext_water[nnz], bc_as_ext_water[nnz];
        
        double un_ext_water, mn_ext_water = 0.0, dmn_ext_water = 0.0, fn_ext_water[nSpace], dfn_ext_water[nSpace], an_ext_water[nnz], dan_ext_water[nnz], asn_ext_water[nnz], bflux_ext_water = 0.0;
        double un_ext_air , mn_ext_air = 0.0, dmn_ext_air = 0.0, fn_ext_air[nSpace], dfn_ext_air[nSpace], an_ext_air[nnz], dan_ext_air[nnz], asn_ext_air[nnz],  bflux_ext_air = 0.0;
        
        double u_ext_air = 0.0, grad_u_ext_air[nSpace], m_ext_air = 0.0, dm_ext_air = 0.0, 
               f_ext_air[nSpace], df_ext_air[nSpace], a_ext_air[nnz], da_ext_air[nnz], as_ext_air[nnz], flux_ext_air = 0.0,
               //anb_seepage_flux=0.0, // for flux calculation
               bc_u_ext_air = 0.0, bc_grad_u_ext_air[nSpace], bc_m_ext_air = 0.0, bc_dm_ext_air = 0.0, 
               bc_f_ext_air[nSpace], bc_df_ext_air[nSpace], bc_a_ext_air[nnz], bc_da_ext_air[nnz], bc_as_ext_air[nnz];         
        
        
        double fluxJacobian_u_u_water[nDOF_trial_element], bfluxJacobian_u_u_water[nDOF_trial_element], fluxJacobian_un_un_water[nDOF_trial_element];

        double fluxJacobian_u_u_air[nDOF_trial_element], bfluxJacobian_u_u_air[nDOF_trial_element], fluxJacobian_un_un_air[nDOF_trial_element];
        
        
        
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
        // ck.valFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], u_ext);

        ck.valFromDOF(u_dof_water.data(), 
                      &u_l2g_water.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_water);
        ck.valFromDOF(u_dof_air.data(), 
                      &u_l2g_air.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      u_ext_air);

        
        ck.valFromDOF(u_dof_old_water.data(), 
                      &u_l2g_water.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      un_ext_water);
        ck.valFromDOF(u_dof_old_air.data(), 
                      &u_l2g_air.data()[eN_nDOF_trial_element], 
                      &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], 
                      un_ext_air);
        //Gradient
        ck.gradFromDOF(u_dof_water.data(), 
                       &u_l2g_water.data()[eN_nDOF_trial_element], 
                       u_grad_trial_trace, 
                       grad_u_ext_water);
        ck.gradFromDOF(u_dof_air.data(), 
                       &u_l2g_air.data()[eN_nDOF_trial_element], 
                       u_grad_trial_trace, 
                       grad_u_ext_air);
               



                      //        ck.valFromDOF(u_dof_old.data(), &u_l2g.data()[eN_nDOF_trial_element], &u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element], un_ext);
       // ck.gradFromDOF(u_dof.data(), &u_l2g.data()[eN_nDOF_trial_element], u_grad_trial_trace, grad_u_ext);
        //precalculate test function products with integration weights
        for (int j = 0; j < nDOF_trial_element; j++) { u_test_dS[j] = u_test_trace_ref.data()[ebN_local_kb * nDOF_test_element + j] * dS; }
        //
        //load the boundary values
        //
        bc_u_ext_water = isDOFBoundary_u_water.data()[ebNE_kb] * ebqe_bc_u_ext_water.data()[ebNE_kb] + (1 - isDOFBoundary_u_water.data()[ebNE_kb]) * u_ext_water;
        bc_u_ext_air = isDOFBoundary_u_air.data()[ebNE_kb] * ebqe_bc_u_ext_air.data()[ebNE_kb] + (1 - isDOFBoundary_u_air.data()[ebNE_kb]) * u_ext_air;
      
//        bc_u_ext = isDOFBoundary_u.data()[ebNE_kb] * ebqe_bc_u_ext.data()[ebNE_kb] + (1 - isDOFBoundary_u.data()[ebNE_kb]) * u_ext;
        //
        //calculate the pde coefficients using the solution and the boundary values for the solution
        //

        double Kr_water, dKr_water,Kr_ext_water, dKr_ext_water, Krn_water, dKrn_water;
        double Kr_air, dKr_air,Kr_ext_air, dKr_ext_air, Krn_air, dKrn_air;


        double bc_Kr_water, bc_dKr_water,bc_Kr_ext_water, bc_dKr_ext_water, bc_Krn_water, bc_dKrn_water;
        double bc_Kr_air, bc_dKr_air,bc_Kr_ext_air, bc_dKr_ext_air, bc_Krn_air, bc_dKrn_air;
        double Sw_ext, Sg_ext, Swn_ext, Sgn_ext;

        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air,
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             u_ext_water, u_ext_air, 
                             Y_water, Y_air,
                             m_ext_water, m_ext_air, 
                             dm_ext_water, dm_ext_air,
                             f_ext_water, f_ext_air,
                             df_ext_water, df_ext_air,
                             a_ext_water, a_ext_air,
                             da_ext_water, da_ext_air, 
                             as_ext_water, as_ext_air,
                             Kr_water, Kr_air,
                             dKr_water, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));



        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air, 
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             un_ext_water, un_ext_air, 
                             Y_water, Y_air,
                             mn_ext_water, mn_ext_air, 
                             dmn_ext_water, dmn_ext_air,
                             fn_ext_water, fn_ext_air,
                             dfn_ext_water, dfn_ext_air,
                             an_ext_water, an_ext_air,
                             dan_ext_water, dan_ext_air, 
                             asn_ext_water, asn_ext_air,
                             Krn_water, Krn_air,
                             dKrn_water, dKrn_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Swn_ext, Sgn_ext,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));



        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                             rho_water, rho_air, 
                             beta_water, beta_air ,
                             gravity.data(), 
                             alpha.data()[elementMaterialTypes.data()[eN]], 
                             n.data()[elementMaterialTypes.data()[eN]], 
                             thetaR.data()[elementMaterialTypes.data()[eN]],
                             thetaSR.data()[elementMaterialTypes.data()[eN]], 
                             &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                             bc_u_ext_water, bc_u_ext_air,
                             Y_water, Y_air, 
                             bc_m_ext_water, bc_m_ext_air, 
                             bc_dm_ext_water, bc_dm_ext_air, 
                             bc_f_ext_water, bc_f_ext_air, 
                             bc_df_ext_water, bc_df_ext_air, 
                             bc_a_ext_water, bc_a_ext_air, 
                             bc_da_ext_water, bc_da_ext_air, 
                             bc_as_ext_water, bc_as_ext_air,
                             Kr_water, Kr_air,
                             dKr_water, dKr_air,
                             PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                             Sw_ext, Sg_ext,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]);        // lambda  (only used if BC_PSK));


         //
        //calculate the numerical fluxes
        //
        bool useConsistentFlux=false;
        if (useConsistentFlux) {
          exteriorNumericalFlux(ebqe_bc_flux_ext_water[ebNE_kb], a_rowptr.data(), a_colind.data(),
                                isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                                isDOFBoundary_u_water.data()[ebNE_kb], normal, bc_u_ext_water, a_ext_water, grad_u_ext_water,
                                u_ext_water, f_ext_water,
                                ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                                flux_ext_water);
          exteriorNumericalFlux(ebqe_bc_flux_ext_air[ebNE_kb], a_rowptr.data(), a_colind.data(),
                                isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                                isDOFBoundary_u_air.data()[ebNE_kb], normal, bc_u_ext_air, a_ext_air, grad_u_ext_air,
                                u_ext_air, f_ext_air,
                                ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                                flux_ext_air);
        } else {
          exteriorNumericalFlux2(ebqe_bc_flux_ext_water[ebNE_kb], a_rowptr.data(), a_colind.data(),
                              isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                              isDOFBoundary_u_water.data()[ebNE_kb], normal, bc_u_ext_water, a_ext_water, grad_u_ext_water, 
                              u_ext_water, f_ext_water,
                              ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                              flux_ext_water, bflux_ext_water);
          exteriorNumericalFlux2(ebqe_bc_flux_ext_air[ebNE_kb], a_rowptr.data(), a_colind.data(),
                              isSeepageFace.data()[ebNE], //tricky, this is a face flag not face quad
                              isDOFBoundary_u_air.data()[ebNE_kb], normal, bc_u_ext_air, a_ext_air, grad_u_ext_air, 
                              u_ext_air, f_ext_air,
                              ebqe_penalty_ext.data()[ebNE_kb], // penalty,
                              flux_ext_air, bflux_ext_air);
                            }
        ebqe_flux_water.data()[ebNE_kb] = flux_ext_water;
        ebqe_flux_air.data()[ebNE_kb] = flux_ext_air;
      
        anb_seepage_flux                 = seepagefluxcalculator(anb_seepage_flux, isSeepageFace.data()[ebNE], dS, flux_ext_water);
        anb_seepage_flux_n.data()[0]     = anb_seepage_flux;
        ebqe_u_water.data()[ebNE_kb]     = u_ext_water;
        ebqe_u_air.data()[ebNE_kb]       = u_ext_air;
        
        //
        //update residuals
        //
        for (int i = 0; i < nDOF_test_element; i++) {
          if (useConsistentFlux) {
            elementResidual_u_water[i] += ck.ExteriorElementBoundaryFlux(flux_ext_water, u_test_dS[i]);
            elementResidual_u_air[i] += ck.ExteriorElementBoundaryFlux(flux_ext_air, u_test_dS[i]);

          } else {
            elementResidual_u_water[i] += ck.ExteriorElementBoundaryFlux(bflux_ext_water, u_test_dS[i]);
            elementResidual_u_air[i] += ck.ExteriorElementBoundaryFlux(bflux_ext_air, u_test_dS[i]);

          }
        } //i
        for (int j = 0; j < nDOF_trial_element; j++) {
          if (useConsistentFlux) {
          exteriorNumericalFluxJacobian(a_rowptr.data(), a_colind.data(), isDOFBoundary_u_water.data()[ebNE_kb], normal, a_ext_water, da_ext_water, grad_u_ext_water, &u_grad_trial_trace[j * nSpace], 
                                        df_ext_water, u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u_water[j]);
          exteriorNumericalFluxJacobian(a_rowptr.data(), a_colind.data(), isDOFBoundary_u_air.data()[ebNE_kb], normal, a_ext_air, da_ext_air, grad_u_ext_air, &u_grad_trial_trace[j * nSpace], 
                                        df_ext_air, u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u_air[j]);
          } 
          else {
            exteriorNumericalFluxJacobian2(a_rowptr.data(), a_colind.data(), isDOFBoundary_u_water.data()[ebNE_kb], normal, as_ext_water, a_ext_water, da_ext_water, grad_u_ext_water, &u_grad_trial_trace[j * nSpace], 
                                        df_ext_water, u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u_water[j],bfluxJacobian_u_u_water[j]);
            exteriorNumericalFluxJacobian2(a_rowptr.data(), a_colind.data(), isDOFBoundary_u_air.data()[ebNE_kb], normal, as_ext_air, a_ext_air, da_ext_air, grad_u_ext_air, &u_grad_trial_trace[j * nSpace], 
                                        df_ext_air, u_trial_trace_ref.data()[ebN_local_kb * nDOF_test_element + j],
                                        ebqe_penalty_ext.data()[ebNE_kb], //penalty,
                                        fluxJacobian_u_u_air[j],bfluxJacobian_u_u_air[j]);

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
              globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_water[j] * u_test_dS[i];
              globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_air[j] * u_test_dS[i];
            } else {
              globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += bfluxJacobian_u_u_water[j] * u_test_dS[i];
              globalJacobian.data()[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += bfluxJacobian_u_u_air[j] * u_test_dS[i];
              //phase 0:: water
              TransportMatrix_water[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_water[j] * u_test_dS[i];
              TransportMatrixConsistent_water[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_water[j] * u_test_dS[i];
              TransportMatrixn_water[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_un_un_water[j] * u_test_dS[i];
              TransportMatrixConsistentn_water[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_un_un_water[j] * u_test_dS[i];

              //phase 1:: air
              TransportMatrix_air[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_air[j] * u_test_dS[i];
              TransportMatrixConsistent_air[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_air[j] * u_test_dS[i];
              TransportMatrixn_air[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_un_un_air[j] * u_test_dS[i];
              TransportMatrixConsistentn_air[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_un_un_air[j] * u_test_dS[i];
            }
          } //j
        } //i
      } //kb
      for (int i = 0; i < nDOF_test_element; i++) {
          int eN_i = eN * nDOF_test_element + i;
          globalResidual.data()[offset_u_water + stride_u_water * u_l2g_water.data()[eN_i]] += elementResidual_u_water[i];
          globalResidual.data()[offset_u_air + stride_u_air * u_l2g_air.data()[eN_i]] += elementResidual_u_air[i];
          
      }//i
    } //ebNE



    /////////////////////////////////////////////////////////////////
    // COMPUTE SMOOTHNESS INDICATOR and NORMALIZE ENTROPY RESIDUAL //
    /////////////////////////////////////////////////////////////////
    // NOTE: see NCLS.h for a different but equivalent implementation of this.
    
    
    int ij = 0;
    double cflux_water[numDOFs];
    for (int i = 0; i < numDOFs; i++) {
      double gi_water[nSpace], Cij[nSpace], xi[nSpace], etaMaxi, etaMini;
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stabilization
      {
        // For eta min and max
        etaMaxi = fabs(eta_water[i]);
        etaMini = fabs(eta_water[i]);
      }
      double solni = u_free_dof_old_water[i];
      // initialize gi and compute xi
      for (int I = 0; I < nSpace; I++) {
        gi_water[I] = 0.;
        xi[I] = mesh_dof.data()[i * 3 + I];
      }
      // for smoothness indicator //
      double alpha_numerator_pos = 0., alpha_numerator_neg = 0., alpha_denominator_pos = 0., alpha_denominator_neg = 0.;
      for (int offset = csrRowIndeces_DofLoops_water.data()[i]; offset < csrRowIndeces_DofLoops_water.data()[i + 1]; offset++) { // First loop in j (sparsity pattern)
        int j = csrColumnOffsets_DofLoops_water.data()[offset];
        if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stabilization
        {
          // COMPUTE ETA MIN AND ETA MAX //
          etaMaxi = fmax(etaMaxi, fabs(eta_water[j]));
          etaMini = fmin(etaMini, fabs(eta_water[j]));
        }
        double solnj = u_free_dof_old_water[j];
        // Update Cij matrices
        Cij[0] = Cx[ij];
#if nSpace == 2
        Cij[1] = Cy[ij];
#endif
#if nSpace == 3
        Cij[2] = Cz[ij];
#endif
        // COMPUTE gi VECTOR. gi=1/mi*sum_j(Cij*solj)
        for (int I = 0; I < nSpace; I++) gi_water[I] += Cij[I] * solnj;

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
      for (int I = 0; I < nSpace; I++) gi_water[I] /= ML_water.data()[i];
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stab
      {
        // Normalizae entropy residual
        global_entropy_residual_water[i] *= etaMini == etaMaxi ? 0. : 2 * cE / (etaMaxi - etaMini);
        quantDOFs_water.data()[i] = fabs(global_entropy_residual_water[i]);
      }

      // Now that I have the gi vectors, I can use them for the current i-th DOF
      double SumPos = 0., SumNeg = 0.;
      for (int offset = csrRowIndeces_DofLoops_water.data()[i]; offset < csrRowIndeces_DofLoops_water.data()[i + 1]; offset++) { // second loop in j (sparsity pattern)
        int j = csrColumnOffsets_DofLoops_water.data()[offset];
        // compute xj
        double xj[nSpace];
        for (int I = 0; I < nSpace; I++) xj[I] = mesh_dof.data()[j * 3 + I];
        // compute gi*(xi-xj)
        double gi_times_x = 0.;
        for (int I = 0; I < nSpace; I++) {
          gi_times_x += gi_water[I] * delta_x_ij_water.data()[offset * 3 + I];
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
      quantDOFs_water.data()[i] = alphai;

      if (POWER_SMOOTHNESS_INDICATOR == 0) psi_water[i] = 1.0;
      else psi_water[i] = std::pow(alphai, POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
    }
// =============== air ================= //

    ij = 0;
    double cflux_air[numDOFs];
    for (int i = 0; i < numDOFs; i++) {
      double gi_air[nSpace], Cij[nSpace], xi[nSpace], etaMaxi, etaMini;
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stabilization
      {
        // For eta min and max
        etaMaxi = fabs(eta_air[i]);
        etaMini = fabs(eta_air[i]);
      }
      double solni = u_free_dof_old_air[i];
      // initialize gi and compute xi
      for (int I = 0; I < nSpace; I++) {
        gi_air[I] = 0.;
        xi[I] = mesh_dof.data()[i * 3 + I];
      }
      // for smoothness indicator //
      double alpha_numerator_pos = 0., alpha_numerator_neg = 0., alpha_denominator_pos = 0., alpha_denominator_neg = 0.;
      for (int offset = csrRowIndeces_DofLoops_air.data()[i]; offset < csrRowIndeces_DofLoops_air.data()[i + 1]; offset++) { // First loop in j (sparsity pattern)
        int j = csrColumnOffsets_DofLoops_air.data()[offset];
        if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stabilization
        {
          // COMPUTE ETA MIN AND ETA MAX //
          etaMaxi = fmax(etaMaxi, fabs(eta_air[j]));
          etaMini = fmin(etaMini, fabs(eta_air[j]));
        }
        double solnj = u_free_dof_old_air[j];
        // Update Cij matrices
        Cij[0] = Cx[ij];
#if nSpace == 2
        Cij[1] = Cy[ij];
#endif
#if nSpace == 3
        Cij[2] = Cz[ij];
#endif
        // COMPUTE gi VECTOR. gi=1/mi*sum_j(Cij*solj)
        for (int I = 0; I < nSpace; I++) gi_air[I] += Cij[I] * solnj;

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
      for (int I = 0; I < nSpace; I++) gi_air[I] /= ML_air.data()[i];
      if (STABILIZATION_TYPE == STABILIZATION::EV_Stab) //EV Stab
      {
        // Normalizae entropy residual
        global_entropy_residual_air[i] *= etaMini == etaMaxi ? 0. : 2 * cE / (etaMaxi - etaMini);
        quantDOFs_air.data()[i] = fabs(global_entropy_residual_air[i]);
      }

      // Now that I have the gi vectors, I can use them for the current i-th DOF
      double SumPos = 0., SumNeg = 0.;
      for (int offset = csrRowIndeces_DofLoops_air.data()[i]; offset < csrRowIndeces_DofLoops_air.data()[i + 1]; offset++) { // second loop in j (sparsity pattern)
        int j = csrColumnOffsets_DofLoops_air.data()[offset];
        // compute xj
        double xj[nSpace];
        for (int I = 0; I < nSpace; I++) xj[I] = mesh_dof.data()[j * 3 + I];
        // compute gi*(xi-xj)
        double gi_times_x = 0.;
        for (int I = 0; I < nSpace; I++) {
          gi_times_x += gi_air[I] * delta_x_ij_air.data()[offset * 3 + I];
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
      quantDOFs_air.data()[i] = alphai;

      if (POWER_SMOOTHNESS_INDICATOR == 0) psi_air[i] = 1.0;
      else psi_air[i] = std::pow(alphai, POWER_SMOOTHNESS_INDICATOR); //NOTE: they use alpha^2 in the paper
    }


/////////////////////////////////////////////
// ** LOOP IN DOFs FOR EDGE BASED TERMS ** //
// =============== WATER ================= //
/////////////////////////////////////////////
ij = 0;
for (int i = 0; i < numDOFs; i++) {
  int    ii;
  double sum_abs_dt_times_fH_minus_fL = 0.0, phi_i  = u_free_dof_water[i], phin_i = u_free_dof_old_water[i], MLi = ML_water.data()[i];
  double Kr, dKr, Krn, dKrn;
  double J_ii = 0.0;
  double ith_dissipative_term           = 0;
  double ith_low_order_dissipative_term = 0;
  double ith_flux_term                  = 0;
  double ith_consistent_flux_term       = 0;
  double dLii                           = 0.;
  // scratch (kept; not all used because bi-phase call returns both)
  double kr_water, dkr_water, krn_water, dkrn_water;
  double kr_air, dkr_air, krn_air, dkrn_air;
  double m_water, dm_water, f_water[nSpace], df_water[nSpace], a_water[nnz], da_water[nnz], as_water[nnz];
  double m_air, dm_air, f_air[nSpace], df_air[nSpace], a_air[nnz], da_air[nnz], as_air[nnz];
  

  double mn_water, dmn_water, fn_water[nSpace], dfn_water[nSpace], an_water[nnz], dan_water[nnz], asn_water[nnz];
  double mn_air, dmn_air, fn_air[nSpace], dfn_air[nSpace], an_air[nnz], dan_air[nnz], asn_air[nnz];
  //double Kr, dKr, Krn, dKrn;
  for (int I = 0; I < nSpace; I++) {
    phi_i  -= rho_water * gravity.data()[I] * mesh_dof.data()[i * 3 + I];
    phin_i -= rho_water * gravity.data()[I] * mesh_dof.data()[i * 3 + I];
  }

  // loop over the sparsity pattern of the i-th DOF
  for (int offset = csrRowIndeces_DofLoops_water.data()[i]; offset < csrRowIndeces_DofLoops_water.data()[i + 1]; offset++) {
    int j = csrColumnOffsets_DofLoops_water.data()[offset];
    if (i == j) ii = ij;
    double phi_j  = u_free_dof_water[j], phin_j = u_free_dof_old_water[j];

    for (int I = 0; I < nSpace; I++) {
      phi_j  -= rho_water * gravity.data()[I] * mesh_dof.data()[j * 3 + I];
      phin_j -= rho_water * gravity.data()[I] * mesh_dof.data()[j * 3 + I];
    }

    double dLowij, dLij, dEVij, dHij, fH, fL, fA=0.0;
    fH = -Theta * TransportMatrixConsistent_water[ij] * (phi_j - phi_i) - (1 - Theta) * TransportMatrixConsistentn_water[ij] * (phin_j - phin_i);
    ith_consistent_flux_term += fH;
    fA = fH;

    // ---- low-order at t^{n+1} (donor based on sign) ----
    if (-TransportMatrix_water[ij] * (phi_j - phi_i) <= 0.0) {
      // donor = i  (use heads at i for BOTH phases; take water kr)
      double Sw, Sg;
      //double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      //double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           /*u_w*/ u_free_dof_water[i], /*u_a*/ u_free_dof_air[i],
                           /*Yw,Ya*/ 1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, 
                           df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg, 
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

      Kr  = kr_water; dKr = dkr_water;
      fL = Theta * Kr * fmax(0.0, -TransportMatrix_water[ij]) * (phi_j - phi_i);
      if (i != j) {
        globalJacobian.data()[ij] -= Theta * Kr * fmax(0.0, -TransportMatrix_water[ij]);
        J_ii -= -Theta * Kr * fmax(0.0, -TransportMatrix_water[ij]) + Theta * dKr * fmax(0.0, -TransportMatrix_water[ij]) * (phi_j - phi_i);
      }
      ith_flux_term += fL;
      fA -= fL;
    } else {
      // donor = j
      double Sw, Sg;
      //double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      //double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           u_free_dof_water[j], u_free_dof_air[j],
                           1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, 
                           df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg, 
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

      Kr  = kr_water; dKr = dkr_water;
      fL = Theta * Kr * fmax(0.0, -TransportMatrix_water[ij]) * (phi_j - phi_i);
      if (i != j) {
        globalJacobian.data()[ij] -= Theta * Kr * fmax(0.0, -TransportMatrix_water[ij]) + Theta * dKr * fmax(0.0, -TransportMatrix_water[ij]) * (phi_j - phi_i);
        J_ii -= -Theta * Kr * fmax(0.0, -TransportMatrix_water[ij]);
      }
      ith_flux_term += fL;
      fA -= fL;
    }

    // ---- low-order at t^n (old) ----
    if (-TransportMatrixn_water[ij] * (phin_j - phin_i) <= 0.0) {
      double Sw, Sg;
      //double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      //double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           u_free_dof_old_water[i], u_free_dof_old_air[i],
                           1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, 
                           df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg, 
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));
                           
      Krn = kr_water;
      fL  = (1 - Theta) * Krn * fmax(0.0, -TransportMatrixn_water[ij]) * (phin_j - phin_i);
      ith_flux_term += fL;
      fA -= fL;
    } else {
      double Sw, Sg;
      //double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      //double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           u_free_dof_old_water[j], u_free_dof_old_air[j],
                           1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg,
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

      Krn = kr_water;
      fL  = (1 - Theta) * Krn * fmax(0.0, -TransportMatrixn_water[ij]) * (phin_j - phin_i);
      ith_flux_term += fL;
      fA -= fL;
    }

    dt_times_fH_minus_fL_water.data()[ij] = dt * fA;
    ij += 1;
  }

  mDotLow_water.data()[i] = ith_flux_term/MLi;
  cflux_water[i] = ith_consistent_flux_term;

  // mass/current at i (take water outputs)
  {
    double Sw, Sg;
    //double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
    //double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

    evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                         rho_water, rho_air, 
                         beta_water, beta_air, 
                         gravity.data(),
                         alpha.data()[elementMaterialTypes.data()[0]],
                         n.data()[elementMaterialTypes.data()[0]],
                         thetaR.data()[elementMaterialTypes.data()[0]],
                         thetaSR.data()[elementMaterialTypes.data()[0]],
                         &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                         u_free_dof_water[i], u_free_dof_air[i],
                         1.0, 1.0,
                         m_water, m_air, 
                         dm_water, dm_air,
                         f_water, f_air, 
                         df_water, df_air,
                         a_water, a_air, 
                         da_water, da_air, 
                         as_water, as_air,
                         kr_water, dkr_water, 
                         kr_air, dkr_air,
                         PSK_TYPE, 
                         Sw, Sg, 
                         BC_entry_head.data()[elementMaterialTypes.data()[0]],
                         BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

    mLow_water.data()[i] = m_water;
    // mass/old
//    double mn_w, mn_a, dmn_w, dmn_a;
    evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                         rho_water, rho_air, 
                         beta_water, beta_air, 
                         gravity.data(),
                         alpha.data()[elementMaterialTypes.data()[0]],
                         n.data()[elementMaterialTypes.data()[0]],
                         thetaR.data()[elementMaterialTypes.data()[0]],
                         thetaSR.data()[elementMaterialTypes.data()[0]],
                         &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                         u_free_dof_old_water[i], u_free_dof_old_air[i],
                         1.0, 1.0,
                         mn_water, mn_air, 
                         dmn_water, dmn_air,
                         fn_water, fn_air, 
                         dfn_water, dfn_air,
                         an_water, an_air, 
                         dan_water, dan_air, 
                         asn_water, asn_air,
                         krn_water, dkrn_water, 
                         krn_air, dkrn_air,
                         PSK_TYPE, 
                         Sw, Sg, 
                         BC_entry_head.data()[elementMaterialTypes.data()[0]],
                         BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));
                         

    globalResidual.data()[i] += bc_mask_water.data()[i] * (MLi * (m_water - mn_water) / dt - ith_flux_term);
    globalJacobian.data()[ii] += bc_mask_water.data()[i] * (MLi * dm_water / dt + J_ii) + (1.0 - bc_mask_water.data()[i]);
  }
}

// high-order rate (water)
ij = 0;
for (int i = 0; i < numDOFs; i++) {
  mDotHigh_water[i] = cflux_water[i];
  for (int offset = csrRowIndeces_DofLoops_water.data()[i]; offset < csrRowIndeces_DofLoops_water.data()[i + 1]; offset++) {
    int j = csrColumnOffsets_DofLoops_water.data()[offset];
    mDotHigh_water[i] -= MC_water.data()[ij]*cflux_water[j]/ML_water.data()[j];
    ij +=1;
  }
  mDotHigh_water[i] = (cflux_water[i] + mDotHigh_water[i])/ML_water.data()[i];
}
if (STABILIZATION_TYPE == STABILIZATION::Implicit_FCT) {
  FCTStep(args);
  for (int i = 0; i < numDOFs; i++) {
    globalResidual.data()[i] += fluxCorrection_water.data()[i];
  }
}

/////////////////////////////////////////////
// ** LOOP IN DOFs FOR EDGE BASED TERMS ** //
// ================ AIR ================== //
/////////////////////////////////////////////
ij = 0;
for (int i = 0; i < numDOFs; i++) {
  int    ii;
  double sum_abs_dt_times_fH_minus_fL = 0.0, phi_i  = u_free_dof_air[i], phin_i = u_free_dof_old_air[i], MLi = ML_air.data()[i];
  double Kr, dKr, Krn, dKrn;
  double J_ii = 0.0;
  double ith_dissipative_term           = 0;
  double ith_low_order_dissipative_term = 0;
  double ith_flux_term                  = 0;
  double ith_consistent_flux_term       = 0;
  double dLii                           = 0.;

  double kr_water, dkr_water, krn_water, dkrn_water;
  double kr_air, dkr_air, krn_air, dkrn_air;

  double m_water, dm_water, f_water[nSpace], df_water[nSpace], a_water[nnz], da_water[nnz], as_water[nnz];
  double m_air, dm_air, f_air[nSpace], df_air[nSpace], a_air[nnz], da_air[nnz], as_air[nnz];
  double mn_water, dmn_water, fn_water[nSpace], dfn_water[nSpace], an_water[nnz], dan_water[nnz], asn_water[nnz];
  double mn_air, dmn_air, fn_air[nSpace], dfn_air[nSpace], an_air[nnz], dan_air[nnz], asn_air[nnz];



  // double m, dm, f[nSpace], df[nSpace], a[nnz], da[nnz], as[nnz];
  // double dmn, fn[nSpace], dfn[nSpace], an[nnz], dan[nnz], asn[nnz];

  for (int I = 0; I < nSpace; I++) {
    phi_i  -= rho_air * gravity.data()[I] * mesh_dof.data()[i * 3 + I];
    phin_i -= rho_air * gravity.data()[I] * mesh_dof.data()[i * 3 + I];
  }

  for (int offset = csrRowIndeces_DofLoops_air.data()[i]; offset < csrRowIndeces_DofLoops_air.data()[i + 1]; offset++) {
    int j = csrColumnOffsets_DofLoops_air.data()[offset];
    if (i == j) ii = ij;
    double phi_j  = u_free_dof_air[j], phin_j = u_free_dof_old_air[j];

    for (int I = 0; I < nSpace; I++) {
      phi_j  -= rho_air * gravity.data()[I] * mesh_dof.data()[j * 3 + I];
      phin_j -= rho_air * gravity.data()[I] * mesh_dof.data()[j * 3 + I];
    }

    double dLowij, dLij, dEVij, dHij, fH, fL, fA=0.0;
    fH = -Theta * TransportMatrixConsistent_air[ij] * (phi_j - phi_i) - (1 - Theta) * TransportMatrixConsistentn_air[ij] * (phin_j - phin_i);
    ith_consistent_flux_term += fH;
    fA = fH;

    // t^{n+1}
    if (-TransportMatrix_air[ij] * (phi_j - phi_i) <= 0.0) {
      double Sw, Sg;
      //double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      //double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           /*u_w*/ u_free_dof_water[i], /*u_a*/ u_free_dof_air[i],
                           /*Yw,Ya*/ 1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, 
                           df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg, 
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

      Kr  = kr_air; dKr = dkr_air;

      fL = Theta * Kr * fmax(0.0, -TransportMatrix_air[ij]) * (phi_j - phi_i);
      if (i != j) {
        globalJacobian.data()[ij] -= Theta * Kr * fmax(0.0, -TransportMatrix_air[ij]);
        J_ii -= -Theta * Kr * fmax(0.0, -TransportMatrix_air[ij]) + Theta * dKr * fmax(0.0, -TransportMatrix_air[ij]) * (phi_j - phi_i);
      }
      ith_flux_term += fL;
      fA -= fL;
    } else {
      double Sw, Sg;
      // double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      // double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           u_free_dof_water[j], u_free_dof_air[j],
                           1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, 
                           df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg, 
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

      Kr  = kr_air; dKr = dkr_air;

      fL = Theta * Kr * fmax(0.0, -TransportMatrix_air[ij]) * (phi_j - phi_i);
      if (i != j) {
        globalJacobian.data()[ij] -= Theta * Kr * fmax(0.0, -TransportMatrix_air[ij]) + Theta * dKr * fmax(0.0, -TransportMatrix_air[ij]) * (phi_j - phi_i);
        J_ii -= -Theta * Kr * fmax(0.0, -TransportMatrix_air[ij]);
      }
      ith_flux_term += fL;
      fA -= fL;
    }

    // t^n
    if (-TransportMatrixn_air[ij] * (phin_j - phin_i) <= 0.0) {
      double  Sw, Sg;
      //double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      //double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           u_free_dof_old_water[i], u_free_dof_old_air[i],
                           1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, 
                           df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg, 
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

      Krn = kr_air;
      fL  = (1 - Theta) * Krn * fmax(0.0, -TransportMatrixn_air[ij]) * (phin_j - phin_i);
      ith_flux_term += fL;
      fA -= fL;
    } else {
      double  Sw, Sg;
      // double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
      // double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

      evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                           rho_water, rho_air, 
                           beta_water, beta_air, 
                           gravity.data(),
                           alpha.data()[elementMaterialTypes.data()[0]],
                           n.data()[elementMaterialTypes.data()[0]],
                           thetaR.data()[elementMaterialTypes.data()[0]],
                           thetaSR.data()[elementMaterialTypes.data()[0]],
                           &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                           u_free_dof_old_water[j], u_free_dof_old_air[j],
                           1.0, 1.0,
                           m_water, m_air, 
                           dm_water, dm_air,
                           f_water, f_air, df_water, df_air,
                           a_water, a_air, 
                           da_water, da_air, 
                           as_water, as_air,
                           kr_water, dkr_water, 
                           kr_air, dkr_air,
                           PSK_TYPE, 
                           Sw, Sg,
                           BC_entry_head.data()[elementMaterialTypes.data()[0]],
                           BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));


      Krn = kr_air;
      fL  = (1 - Theta) * Krn * fmax(0.0, -TransportMatrixn_air[ij]) * (phin_j - phin_i);
      ith_flux_term += fL;
      fA -= fL;
    }

    dt_times_fH_minus_fL_air.data()[ij] = dt * fA;
    ij += 1;
  }

  mDotLow_air.data()[i] = ith_flux_term/MLi;
  cflux_air[i] = ith_consistent_flux_term;

  // mass/current at i (take air outputs)
  {
    double Sw, Sg;
    // double f_w[nSpace], f_a[nSpace], df_w[nSpace], df_a[nSpace];
    // double a_w[nnz], a_a[nnz], da_w[nnz], da_a[nnz], as_w[nnz], as_a[nnz];

    evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                         rho_water, rho_air, 
                         beta_water, beta_air, 
                         gravity.data(),
                         alpha.data()[elementMaterialTypes.data()[0]],
                         n.data()[elementMaterialTypes.data()[0]],
                         thetaR.data()[elementMaterialTypes.data()[0]],
                         thetaSR.data()[elementMaterialTypes.data()[0]],
                         &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                         u_free_dof_water[i], u_free_dof_air[i],
                         1.0, 1.0,
                         m_water, m_air, 
                         dm_water, dm_air,
                         f_water, f_air, 
                         df_water, df_air,
                         a_water, a_air, 
                         da_water, da_air, 
                         as_water, as_air,
                         kr_water, dkr_water, 
                         kr_air, dkr_air,
                         PSK_TYPE, 
                         Sw, Sg, 
                         BC_entry_head.data()[elementMaterialTypes.data()[0]],
                         BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));

    mLow_air.data()[i] = m_air;

    // mass/old
//    double mn_w, mn_a, dmn_w, dmn_a;
    evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                         rho_water, rho_air, 
                         beta_water, beta_air, 
                         gravity.data(),
                         alpha.data()[elementMaterialTypes.data()[0]],
                         n.data()[elementMaterialTypes.data()[0]],
                         thetaR.data()[elementMaterialTypes.data()[0]],
                         thetaSR.data()[elementMaterialTypes.data()[0]],
                         &KWs.data()[elementMaterialTypes.data()[0] * nnz],
                         u_free_dof_old_water[i], u_free_dof_old_air[i],
                         1.0, 1.0,
                         mn_water, mn_air, 
                         dmn_water, dmn_air,
                         fn_water, fn_air, 
                         dfn_water, dfn_air,
                         an_water, an_air, 
                         dan_water, dan_air, 
                         asn_water, asn_air,
                         krn_water, dkrn_water, 
                         krn_air, dkrn_air,
                         PSK_TYPE, 
                         Sw, Sg, 
                         BC_entry_head.data()[elementMaterialTypes.data()[0]],
                         BC_lambda.data()[elementMaterialTypes.data()[0]]);        // lambda  (only used if BC_PSK));


    globalResidual.data()[i] += bc_mask_air.data()[i] * (MLi * (m_air - mn_air) / dt - ith_flux_term);
    globalJacobian.data()[ii] += bc_mask_air.data()[i] * (MLi * dm_air / dt + J_ii) + (1.0 - bc_mask_air.data()[i]);
  }
}

// high-order rate (air)
ij = 0;
for (int i = 0; i < numDOFs; i++) {
  mDotHigh_air[i] = cflux_air[i];
  for (int offset = csrRowIndeces_DofLoops_air.data()[i]; offset < csrRowIndeces_DofLoops_air.data()[i + 1]; offset++) {
    int j = csrColumnOffsets_DofLoops_air.data()[offset];
    mDotHigh_air[i] -= MC_air.data()[ij]*cflux_air[j]/ML_air.data()[j];
    ij +=1;
  }
  mDotHigh_air[i] = (cflux_air[i] + mDotHigh_air[i])/ML_air.data()[i];
}
if (STABILIZATION_TYPE == STABILIZATION::Implicit_FCT) {
  FCTStep(args);
  for (int i = 0; i < numDOFs; i++) {
    globalResidual.data()[i] += fluxCorrection_air.data()[i];
  }
}
}
  
void invert(arguments_dict &args)
{
  // --- common tables/params (unchanged shape) ---
  xt::pyarray<int>    &a_rowptr             = args.array<int>("a_rowptr");
  xt::pyarray<int>    &a_colind             = args.array<int>("a_colind");
  xt::pyarray<double> &gravity              = args.array<double>("gravity");
  xt::pyarray<double> &alpha                = args.array<double>("alpha");
  xt::pyarray<double> &n                    = args.array<double>("n");
  xt::pyarray<double> &thetaR               = args.array<double>("thetaR");
  xt::pyarray<double> &thetaSR              = args.array<double>("thetaSR");
  xt::pyarray<double> &KWs                  = args.array<double>("KWs");
  xt::pyarray<int>    &elementMaterialTypes = args.array<int>("elementMaterialTypes");
  int                  numDOFs              = args.scalar<int>("numDOFs");

  // --- two-phase props & IO (NEW) ---
  double rho_water  = args.scalar<double>("rho_water");
  double rho_air    = args.scalar<double>("rho_air");
  double beta_water = args.scalar<double>("beta_water");
  double beta_air   = args.scalar<double>("beta_air");
  double Y_water    = args.scalar<double>("Y_water");  // set to 1.0 if not using mass fractions
  double Y_air      = args.scalar<double>("Y_air");

  xt::pyarray<double> &mIn_water = args.array<double>("limited_solution_water"); // per-DOF water mass
  xt::pyarray<double> &mIn_air   = args.array<double>("limited_solution_air");   // per-DOF air   mass
  xt::pyarray<double> &pOut_water= args.array<double>("u_dof_water");            // per-DOF water head (in/out)
  xt::pyarray<double> &pOut_air  = args.array<double>("u_dof_air");              // per-DOF air   head (in/out)

  // --- PSK (unchanged) ---
  PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
  xt::pyarray<double> &BC_entry_head = args.array<double>("BC_entry_head"); // size: nMaterials
  xt::pyarray<double> &BC_lambda     = args.array<double>("BC_lambda");     // size: nMaterials
  // double BC_entry_head = args.scalar<double>("BC_entry_head");
  // double BC_lambda     = args.scalar<double>("BC_lambda");

  // single material hack preserved
  //const int matId = elementMaterialTypes.data()[0];

  for (int i = 0; i < numDOFs; i++) {
    // keep scratch locals to match prior structure (unused by new inverse)
    double dm, f[nSpace], df[nSpace], a[nnz], da[nnz];

    // ---- mass sanity checks (mirrors your previous pattern) ----
    const double mMin_w = rho_water * thetaR.data()[elementMaterialTypes.data()[0]];
    const double mMax_w = rho_water * (thetaR.data()[elementMaterialTypes.data()[0]] + thetaSR.data()[elementMaterialTypes.data()[0]]);
    if (mIn_water.data()[i] < mMin_w - 1.0e-3 || mIn_water.data()[i] > mMax_w + 1.0e-3) {
      std::cout << "WATER mass out of bounds " << mMin_w << '\t' << mIn_water.data()[i]
                << '\t' << mMax_w << std::endl;
    }

    const double mMin_a = rho_air * thetaR.data()[elementMaterialTypes.data()[0]]; // use symmetric bounds; set to 0 if you prefer
    const double mMax_a = rho_air * (thetaR.data()[elementMaterialTypes.data()[0]] + thetaSR.data()[elementMaterialTypes.data()[0]]);
    if (mIn_air.data()[i] < mMin_a - 1.0e-3 || mIn_air.data()[i] > mMax_a + 1.0e-3) {
      std::cout << "AIR mass out of bounds " << mMin_a << '\t' << mIn_air.data()[i]
                << '\t' << mMax_a << std::endl;
    }

    // heads in/out (both needed)
    double uw = pOut_water.data()[i];
    double ua = pOut_air.data()[i];

    // ---- NEW: single call updates BOTH phases from BOTH masses ----
    evaluateInverseCoefficients_2ph(
      a_rowptr.data(), a_colind.data(),
      rho_water, rho_air,
      beta_water, beta_air,
      gravity.data(),
      alpha.data()[elementMaterialTypes.data()[0]],
      n.data()[elementMaterialTypes.data()[0]],
      thetaR.data()[elementMaterialTypes.data()[0]], 
      thetaSR.data()[elementMaterialTypes.data()[0]],
      &KWs.data()[elementMaterialTypes.data()[0] * nnz],
      uw, ua,                       // in/out
      mIn_water.data()[i], 
      mIn_air.data()[i],
      Y_water, Y_air,
      PSK_TYPE,
      BC_entry_head[elementMaterialTypes.data()[0]], 
      BC_lambda[elementMaterialTypes.data()[0]]
    );

    // write back heads
    pOut_water.data()[i] = uw;
    pOut_air.data()[i]   = ua;
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
    double               rho_water            = args.scalar<double>("rho_water");
    double               rho_air              = args.scalar<double>("rho_air");
    double               beta_water           = args.scalar<double>("beta_water");    
    double               beta_air             = args.scalar<double>("beta_air");
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
    xt::pyarray<int>    &u_l2g_water                                      = args.array<int>("u_l2g_water");
    xt::pyarray<int>    &r_l2g_water                                      = args.array<int>("r_l2g_water");
    xt::pyarray<int>    &u_l2g_air                                      = args.array<int>("u_l2g_air");
    xt::pyarray<int>    &r_l2g_air                                      = args.array<int>("r_l2g_air");
    
    xt::pyarray<double> &elementDiameter                            = args.array<double>("elementDiameter");
    int                  degree_polynomial                          = args.scalar<int>("degree_polynomial");
    xt::pyarray<double> &u_dof_water                                = args.array<double>("u_dof_water");
    xt::pyarray<double> &u_dof_air                                  = args.array<double>("u_dof_air");    
    xt::pyarray<double> &velocity_water                             = args.array<double>("velocity_water");
    xt::pyarray<double> &velocity_air                               = args.array<double>("velocity_air");
    xt::pyarray<double> &q_m_betaBDF_water                                = args.array<double>("q_m_betaBDF_water");    
    xt::pyarray<double> &q_m_betaBDF_air                                = args.array<double>("q_m_betaBDF_air");
    xt::pyarray<double> &cfl                                        = args.array<double>("cfl");
    xt::pyarray<double> &q_numDiff_u_last_water                     = args.array<double>("q_numDiff_u_last_water");
    xt::pyarray<double> &q_numDiff_u_last_air                       = args.array<double>("q_numDiff_u_last_air");
    
    xt::pyarray<int>    &csrRowIndeces_w_w                          = args.array<int>("csrRowIndeces_w_w");
    xt::pyarray<int>    &csrColumnOffsets_w_w                       = args.array<int>("csrColumnOffsets_w_w");
    xt::pyarray<int>    &csrRowIndeces_a_a                          = args.array<int>("csrRowIndeces_a_a");
    xt::pyarray<int>    &csrColumnOffsets_a_a                       = args.array<int>("csrColumnOffsets_a_a");
  

  //  xt::pyarray<int>    &csrRowIndeces_u_u                          = args.array<int>("csrRowIndeces_u_u");
  //  xt::pyarray<int>    &csrColumnOffsets_u_u                       = args.array<int>("csrColumnOffsets_u_u");
    xt::pyarray<double> &globalJacobian                             = args.array<double>("globalJacobian");
//    xt::pyarray<double> &globalJacobian                             = args.array<double>("globalJacobian");
    
    xt::pyarray<double> &delta_x_ij                                 = args.array<double>("delta_x_ij");
    int                   nExteriorElementBoundaries_global          = args.scalar<int>("nExteriorElementBoundaries_global");
    xt::pyarray<int>    &exteriorElementBoundariesArray             = args.array<int>("exteriorElementBoundariesArray");
    xt::pyarray<int>    &elementBoundaryElementsArray               = args.array<int>("elementBoundaryElementsArray");
    xt::pyarray<int>    &elementBoundaryLocalElementBoundariesArray = args.array<int>("elementBoundaryLocalElementBoundariesArray");
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
    xt::pyarray<int>    &csrColumnOffsets_eb_u_u                    = args.array<int>("csrColumnOffsets_eb_u_u");
        //////////////////////////////////////For Brooks- Corey/////////////////
    PSK PSK_TYPE{static_cast<PSK>(args.scalar<int>("PSK_MODEL"))};
    xt::pyarray<double> &BC_entry_head = args.array<double>("BC_entry_head"); // size: nMaterials
    xt::pyarray<double> &BC_lambda     = args.array<double>("BC_lambda");     // size: nMaterials
    double Y_air                                                    = args.scalar<double>("Y_air");     // e.g. 1.0
    double Y_water                                                  = args.scalar<double>("Y_water");   // e.g. 1.0
   
    // double BC_entry_head = args.scalar<double>("BC_entry_head");
    // double BC_lambda     = args.scalar<double>("BC_lambda");  
    int                  LUMPED_MASS_MATRIX                         = args.scalar<int>("LUMPED_MASS_MATRIX");
    double Ct_sge = 4.0;
    //
    //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
    //
    for (int eN = 0; eN < nElements_global; eN++) {
      double elementJacobian_u_u_water[nDOF_test_element][nDOF_trial_element];
      double elementJacobian_u_u_air[nDOF_test_element][nDOF_trial_element];
      
      for (int i = 0; i < nDOF_test_element; i++)
        for (int j = 0; j < nDOF_trial_element; j++) { 
      elementJacobian_u_u_water[i][j] = 0.0;
      elementJacobian_u_u_air[i][j] = 0.0;
     }
     for (int k = 0; k < nQuadraturePoints_element; k++) {
        int eN_k                  = eN * nQuadraturePoints_element + k, //index to a scalar at a quadrature point
          eN_k_nSpace             = eN_k * nSpace,
            eN_nDOF_trial_element = eN * nDOF_trial_element; //index to a vector at a quadrature point

        //declare local storage
        double jac[nSpace * nSpace], jacDet, jacInv[nSpace * nSpace], u_grad_trial[nDOF_trial_element * nSpace], dV, u_test_dV[nDOF_test_element], u_grad_test_dV[nDOF_test_element * nSpace], x, y, z, xt, yt, zt,
          G[nSpace * nSpace], G_dd_G, tr_G;
        double u_water = 0.0, grad_u_water[nSpace], m_water = 0.0, dm_water = 0.0, f_water[nSpace], df_water[nSpace], 
        a_water[nnz], da_water[nnz], as_water[nnz], m_t_water = 0.0, dm_t_water = 0.0, dpdeResidual_u_u_water[nDOF_trial_element], Lstar_u_water[nDOF_test_element], 
        dsubgridError_u_u_water[nDOF_trial_element], tau_water = 0.0, tau0_water = 0.0, tau1_water = 0.0;

        double u_air = 0.0, grad_u_air[nSpace], m_air = 0.0, dm_air = 0.0, f_air[nSpace], df_air[nSpace], 
        a_air[nnz], da_air[nnz], as_air[nnz], m_t_air = 0.0, dm_t_air = 0.0, dpdeResidual_u_u_air[nDOF_trial_element], Lstar_u_air[nDOF_test_element], 
        dsubgridError_u_u_air[nDOF_trial_element], tau_air = 0.0, tau0_air = 0.0, tau1_air = 0.0 ;
        
        
        //get jacobian, etc for mapping reference element
        ck.calculateMapping_element(eN, k, mesh_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), mesh_grad_trial_ref.data(), jac, jacDet, jacInv, x, y, z);
        ck.calculateMappingVelocity_element(eN, k, mesh_velocity_dof.data(), mesh_l2g.data(), mesh_trial_ref.data(), xt, yt, zt);
        //get the physical integration weight
        dV = fabs(jacDet) * dV_ref.data()[k];
        ck.calculateG(jacInv, G, G_dd_G, tr_G);
        //get the trial function gradients
        ck.gradTrialFromRef(&u_grad_trial_ref.data()[k * nDOF_trial_element * nSpace], jacInv, u_grad_trial);
        //get the solution
        ck.valFromDOF(u_dof_water.data(), &u_l2g_water.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], u_water);
        ck.valFromDOF(u_dof_air.data(), &u_l2g_air.data()[eN_nDOF_trial_element], &u_trial_ref.data()[k * nDOF_trial_element], u_air);
        

        //get the solution gradients
        ck.gradFromDOF(u_dof_water.data(), &u_l2g_water.data()[eN_nDOF_trial_element], u_grad_trial, grad_u_water);
        ck.gradFromDOF(u_dof_air.data(), &u_l2g_air.data()[eN_nDOF_trial_element], u_grad_trial, grad_u_air);
        
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
        double Sw=0.0, Sg= 0.0;    
        evaluateCoefficients(a_rowptr.data(), a_colind.data(),
                            rho_water, rho_air, 
                            beta_water, beta_air, 
                            gravity.data(), 
                            alpha.data()[elementMaterialTypes.data()[eN]], 
                            n.data()[elementMaterialTypes.data()[eN]], 
                            thetaR.data()[elementMaterialTypes.data()[eN]],
                            thetaSR.data()[elementMaterialTypes.data()[eN]], 
                            &KWs.data()[elementMaterialTypes.data()[eN] * nnz], 
                            u_water, u_air, 
                            Y_water, Y_air,
                            m_water, m_air, 
                            dm_water, dm_air,
                            f_water, f_air,
                            df_water, df_air,
                            a_water, a_air,
                            da_water, da_air, 
                            as_water, as_air,
                            Kr_water, dKr_water,
                            Kr_air, dKr_air,
                            PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
                            Sw, Sg,
                            BC_entry_head.data()[elementMaterialTypes.data()[eN]],
                            BC_lambda.data()[elementMaterialTypes.data()[eN]]); 

        // evaluateCoefficients(a_rowptr.data(), a_colind.data(), rho, beta, gravity.data(), alpha.data()[elementMaterialTypes.data()[eN]], n.data()[elementMaterialTypes.data()[eN]], thetaR.data()[elementMaterialTypes.data()[eN]],
        //                      thetaSR.data()[elementMaterialTypes.data()[eN]], &KWs.data()[elementMaterialTypes.data()[eN] * nnz], u, m, dm, f, df, a, da, as, Kr, dKr,
        //                      PSK_TYPE,          // 0: VG_PSK, 1: BC_PSK
        //                      BC_entry_head,     // h_e (only used if BC_PSK)
        //                      BC_lambda);
        //
        //moving mesh
        //
        double mesh_velocity[3];
        mesh_velocity[0] = xt;
        mesh_velocity[1] = yt;
        mesh_velocity[2] = zt;
        for (int I = 0; I < nSpace; I++) {
          f_water[I] -= MOVING_DOMAIN * m_water * mesh_velocity[I];
          df_water[I] -= MOVING_DOMAIN * dm_water * mesh_velocity[I];
          f_air[I] -= MOVING_DOMAIN * m_air * mesh_velocity[I];
          df_air[I] -= MOVING_DOMAIN * dm_air * mesh_velocity[I];

        }
        //
        //calculate time derivatives
        //
        //cek hack
        dm_water = 1.0;
        dm_air = 1.0;

        ck.bdf(alphaBDF,
               q_m_betaBDF_water.data()[eN_k], //since m_t isn't used, we don't have to correct mass
               m_water, dm_water, m_t_water, dm_t_water);
        
        ck.bdf(alphaBDF,
               q_m_betaBDF_air.data()[eN_k], //since m_t isn't used, we don't have to correct mass
               m_air, dm_air, m_t_air, dm_t_air);
        
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
          dpdeResidual_u_u_water[j] = ck.MassJacobian_strong(dm_t_water, u_trial_ref.data()[k * nDOF_trial_element + j]) + 
                                      ck.AdvectionJacobian_strong(df_water, &u_grad_trial[j_nSpace]);
          dpdeResidual_u_u_air[j]   = ck.MassJacobian_strong(dm_t_air, u_trial_ref.data()[k * nDOF_trial_element + j]) + 
                                      ck.AdvectionJacobian_strong(df_air, &u_grad_trial[j_nSpace]);
        }
        //tau and tau*Res
        calculateSubgridError_tau(elementDiameter.data()[eN], dm_t_water, df_water, cfl.data()[eN_k], tau0_water);
        calculateSubgridError_tau(Ct_sge, G, dm_t_water, df_water, tau1_water, cfl.data()[eN_k]);
        tau_water = useMetrics * tau1_water + (1.0 - useMetrics) * tau0_water;

        calculateSubgridError_tau(elementDiameter.data()[eN], dm_t_air, df_air, cfl.data()[eN_k], tau0_air);
        calculateSubgridError_tau(Ct_sge, G, dm_t_air, df_air, tau1_air, cfl.data()[eN_k]);
        tau_air = useMetrics * tau1_air + (1.0 - useMetrics) * tau0_air;

        for (int j = 0; j < nDOF_trial_element; j++) 
        {
          dsubgridError_u_u_water[j] = -tau_water * dpdeResidual_u_u_water[j];
          dsubgridError_u_u_air[j] = -tau_air * dpdeResidual_u_u_air[j];
        }
          for (int i = 0; i < nDOF_test_element; i++) {
          for (int j = 0; j < nDOF_trial_element; j++) {
            if (LUMPED_MASS_MATRIX == 1) {
              if (i == j) 
              elementJacobian_u_u_water[i][j] += u_test_dV[i];
              elementJacobian_u_u_air[i][j] += u_test_dV[i];              
            } else {
              int j_nSpace = j * nSpace;
              int i_nSpace = i * nSpace;
              dm_t_water = 1.0; //we are solving for continuum density explicitly
              dm_t_air = 1.0;
              elementJacobian_u_u_water[i][j] += ck.MassJacobian_weak(dm_t_water, u_trial_ref.data()[k * nDOF_trial_element + j], u_test_dV[i]);
              elementJacobian_u_u_air[i][j] += ck.MassJacobian_weak(dm_t_air, u_trial_ref.data()[k * nDOF_trial_element + j], u_test_dV[i]);
            }
          } //j
        } //i
      } //k
      //
      //load into element Jacobian into global Jacobian
      //
      for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        int I    = u_l2g_water.data()[eN_i];
        for (int j = 0; j < nDOF_trial_element; j++) {
          int eN_i_j = eN_i * nDOF_trial_element + j;
          int J      = u_l2g_water.data()[eN * nDOF_trial_element + j];
          //globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += elementJacobian_u_u[i][j];
          delta_x_ij.data()[3 * (csrRowIndeces_w_w.data()[eN_i] + csrColumnOffsets_w_w.data()[eN_i_j]) + 0] = mesh_dof.data()[I * 3 + 0] - mesh_dof.data()[J * 3 + 0];
          delta_x_ij.data()[3 * (csrRowIndeces_w_w.data()[eN_i] + csrColumnOffsets_w_w.data()[eN_i_j]) + 1] = mesh_dof.data()[I * 3 + 1] - mesh_dof.data()[J * 3 + 1];
          delta_x_ij.data()[3 * (csrRowIndeces_w_w.data()[eN_i] + csrColumnOffsets_w_w.data()[eN_i_j]) + 2] = mesh_dof.data()[I * 3 + 2] - mesh_dof.data()[J * 3 + 2];
        } //j
      } //i

        for (int i = 0; i < nDOF_test_element; i++) {
        int eN_i = eN * nDOF_test_element + i;
        int I    = u_l2g_air.data()[eN_i];
        for (int j = 0; j < nDOF_trial_element; j++) {
          int eN_i_j = eN_i * nDOF_trial_element + j;
          int J      = u_l2g_air.data()[eN * nDOF_trial_element + j];
          //globalJacobian.data()[csrRowIndeces_u_u.data()[eN_i] + csrColumnOffsets_u_u.data()[eN_i_j]] += elementJacobian_u_u[i][j];
          delta_x_ij.data()[3 * (csrRowIndeces_a_a.data()[eN_i] + csrColumnOffsets_a_a.data()[eN_i_j]) + 0] = mesh_dof.data()[I * 3 + 0] - mesh_dof.data()[J * 3 + 0];
          delta_x_ij.data()[3 * (csrRowIndeces_a_a.data()[eN_i] + csrColumnOffsets_a_a.data()[eN_i_j]) + 1] = mesh_dof.data()[I * 3 + 1] - mesh_dof.data()[J * 3 + 1];
          delta_x_ij.data()[3 * (csrRowIndeces_a_a.data()[eN_i] + csrColumnOffsets_a_a.data()[eN_i_j]) + 2] = mesh_dof.data()[I * 3 + 2] - mesh_dof.data()[J * 3 + 2];
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
