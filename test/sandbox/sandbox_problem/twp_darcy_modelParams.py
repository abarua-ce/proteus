# Set the model parameters to be given to the impes model for two-phase flow.

#choose in problem-specific class now
#model = 'VGM'#'simp'#'VGM'  # Choose from {simp,BCB,BCM,VGB,VGM}
#for now have to have PorousMedia set here
PorousMedia = 'clay-2'  # 'sand' or 'clay-2'

# Model Parameters
rhow = 997.0  # density wetting (Water 997 kg/m^3)
rhon = 1.205  # density nonwetting (Air 1.205 kg/m^3)
muw  = 1.002e-3  # viscosity wetting (Water := 1.002e-3 kg/m s)
mun  = 1.81e-5  # viscosity nonwetting (Air := 1.81e-5 kg/ m s)

gmag = 9.8  # magnitude of gravity (Earth 9.8 m/s^2)
m_per_d_to_m_per_s = 1.1574074e-5

# Water density parameters
beta_w = rhow * gmag * 4.524e-10
density_w_exponential = {
    'model': 'Exponential',
    'rho_0': rhow,
    'psi_0': 0.0,
    'beta': beta_w
}

# Air density parameters
Temp = 288.15  # Temperature
R = 8.314  # Ideal Gas Constant (N/mol K)
M = 2.896e-2  # Molar Mass (Air := 2.896e-2)
p_o = 1.013e5  # (N/m^2)
W = 2.9e-4  # molar weight [kg/mol]
m_to_Newton_per_m2 = rhow * gmag  # mwf should be rhow * gmag

density_n_ideal = {
    'model': 'IdealGas',
    'T': Temp,
    'W': W,
    'R': R,
    'headToPressure': m_to_Newton_per_m2,
    'rho_0': rhon,
    'psi_0': 0.0
}

a = 1.0e-6

# Porous medium data
porousMediumDatabase = {
    'sand': {},
    'loam': {},
    'clay': {},
    'clay-2': {},
    'clay-3': {}
}

# Populate the database
porousMediumDatabase['sand'].update({
    'sw_min': 0.308, 'sw_max': 1.0, 'omega': 0.301, 'Ks': 5.040,
    'mvg_alpha': 5.470, 'mvg_n': 4.264, 'bc_pd': 0.120, 'bc_lambda': 3.264
})

porousMediumDatabase['loam'].update({
    'sw_min': 0.308, 'sw_max': 1.0, 'omega': 0.301, 'Ks': 5.040,
    'mvg_alpha': 5.470, 'mvg_n': 4.264, 'bc_pd': 0.120, 'bc_lambda': 3.264
})

porousMediumDatabase['clay'].update({
    'sw_min': 0.232, 'sw_max': 1.0, 'omega': 0.410, 'Ks': 0.062,
    'mvg_alpha': 1.900, 'mvg_n': 1.310, 'bc_pd': 0.044, 'bc_lambda': 0.310
})

porousMediumDatabase['clay-2'].update({
    'sw_min': 0.277, 'sw_max': 1.0, 'omega': 0.368, 'Ks': 7.970,
    'mvg_alpha': 3.350, 'mvg_n': 2.00, 'bc_pd': 1.0 / 3.350, 'bc_lambda': 1.00
})

porousMediumDatabase['clay-3'].update({
    'sw_min': 0.277, 'sw_max': 1.0, 'omega': 0.268, 'Ks': 7.970e-2,
    'mvg_alpha': 3.350, 'mvg_n': 2.00, 'bc_pd': 1.0 / 3.350, 'bc_lambda': 1.00
})

# Add missing values
for medium in porousMediumDatabase.keys():
    porousMediumDatabase[medium]['mvg_m'] = 1.0 - 1.0 / porousMediumDatabase[medium]['mvg_n']
    porousMediumDatabase[medium]['Ksw'] = porousMediumDatabase[medium]['Ks'] * m_per_d_to_m_per_s

# Check PorousMedia in the database using `in`
assert PorousMedia in porousMediumDatabase, f"PorousMedia={PorousMedia} not found in database keys= {list(porousMediumDatabase.keys())}"

# Assign default values based on PorousMedia
sw_min = porousMediumDatabase[PorousMedia]['sw_min']
sw_max = porousMediumDatabase[PorousMedia]['sw_max']
omega = porousMediumDatabase[PorousMedia]['omega']
Ks = porousMediumDatabase[PorousMedia]['Ks']
mvg_alpha = porousMediumDatabase[PorousMedia]['mvg_alpha']
mvg_n = porousMediumDatabase[PorousMedia]['mvg_n']
bc_pd = porousMediumDatabase[PorousMedia]['bc_pd']
bc_lambda = porousMediumDatabase[PorousMedia]['bc_lambda']
mvg_m = porousMediumDatabase[PorousMedia]['mvg_m']
Ksw = porousMediumDatabase[PorousMedia]['Ksw']

# Helper functions
def seVGM(psic, alVG, nVG, mVG):
    if psic <= 0:
        return 1.0
    tmp1 = pow(alVG * psic, nVG)
    tmp2 = pow(1. + tmp1, -mVG)
    return min(max(tmp2, 0.0), 1.0)

def pcVGM(se, alVG, nVG, mVG):
    if se >= 1.0:
        return 0.0
    tmp1 = pow(se, -1. / mVG)
    tmp2 = pow(tmp1 - 1.0, 1.0 / nVG) / alVG
    return tmp2

def seBCB(psic, pdBC, lamBC):
    if psic <= pdBC:
        return 1.0
    tmp1 = pow(pdBC / psic, lamBC)
    return min(max(tmp1, 0.0), 1.0)

def pcBCB(se, pdBC, lamBC):
    if se >= 1.0:
        return 0.0
    tmp1 = pow(se, -1.0 / lamBC)
    return pdBC * tmp1
