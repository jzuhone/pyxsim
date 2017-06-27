import numpy as np
from yt.frontends.stream.api import load_particles
from yt.units.yt_array import uconcatenate, YTArray, \
    YTQuantity
from yt.utilities.physical_ratios import keV_per_erg
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline
from six import string_types
from pyxsim.photon_list import PhotonList
from pyxsim.source_models import PowerLawSourceModel
from pyxsim.utils import mylog
from soxs.utils import parse_prng

"""
Papers referenced throughout this code:

Gilfanov, M. 2004, MNRAS, 349, 146

Fragos, T., Lehmer, B., Tremmel, M., et al. 2013, ApJ, 764, 41

Mineo, S., Gilfanov, M., & Sunyaev, R. 2012, MNRAS, 419, 2095 

"""

# Begin Fragos et al. 2013 arrays derived from Figure 2

# Luminosity as function of age arrays

t_lm = np.array([7.07038000e-03,   7.42810000e-03,   7.90964000e-03,
                 8.42238000e-03,   8.92825000e-03,   9.42208000e-03,
                 1.01688000e-02,   1.07796000e-02,   1.14271000e-02,
                 1.24438000e-02,   1.35510000e-02,   1.46905000e-02,
                 1.58546000e-02,   1.73426000e-02,   1.86328000e-02,
                 2.06573000e-02,   2.26973000e-02,   2.55040000e-02,
                 2.75246000e-02,   3.02433000e-02,   3.20596000e-02,
                 3.46001000e-02,   3.66781000e-02,   3.95844000e-02,
                 4.27204000e-02,   4.75741000e-02,   5.32194000e-02,
                 5.64152000e-02,   6.06131000e-02,   6.48309000e-02,
                 6.78033000e-02,   7.28446000e-02,   7.79105000e-02,
                 8.44592000e-02,   9.19716000e-02,   1.00151000e-01,
                 1.07600000e-01,   1.17169000e-01,   1.28739000e-01,
                 1.45309000e-01,   1.60375000e-01,   1.78599000e-01,
                 1.97117000e-01,   2.16583000e-01,   2.35846000e-01,
                 2.53392000e-01,   2.69811000e-01,   2.92495000e-01,
                 3.19943000e-01,   3.45294000e-01,   3.76011000e-01,
                 4.02184000e-01,   4.32109000e-01,   4.66355000e-01,
                 5.05547000e-01,   5.43143000e-01,   6.04851000e-01,
                 6.67557000e-01,   7.36767000e-01,   8.39077000e-01,
                 9.38607000e-01,   1.03591000e+00,   1.16400000e+00,
                 1.26184000e+00,   1.35566000e+00,   1.44995000e+00,
                 1.59307000e+00,   1.72693000e+00,   1.84702000e+00,
                 1.99325000e+00,   2.12232000e+00,   2.28011000e+00,
                 2.46054000e+00,   2.59645000e+00,   2.75216000e+00,
                 2.87822000e+00,   3.05081000e+00,   3.21929000e+00,
                 3.32181000e+00,   3.48957000e+00,   3.60070000e+00,
                 3.73205000e+00,   3.90304000e+00,   4.17444000e+00,
                 4.42481000e+00,   4.71124000e+00,   5.01629000e+00,
                 5.41341000e+00,   5.71239000e+00,   6.10952000e+00,
                 6.36079000e+00,   6.56327000e+00,   6.80264000e+00,
                 6.98759000e+00,   7.17781000e+00,   7.64251000e+00,
                 8.28517000e+00,   8.39777000e+00,   8.51198000e+00,
                 8.82334000e+00,   9.22836000e+00,   9.69535000e+00,
                 1.00041000e+01,   1.05096000e+01,   1.09911000e+01,
                 1.13411000e+01,   1.18081000e+01,   1.19676000e+01,
                 1.21835000e+01,   1.24031000e+01,   1.27410000e+01,
                 1.27976000e+01,   1.28546000e+01])

L_lm = np.array([ 39.8482,  39.8886,  39.9326,  39.973 ,  40.017 ,  40.05  ,
                  40.1051,  40.1528,  40.1931,  40.2445,  40.2885,  40.3215,
                  40.3582,  40.3876,  40.4022,  40.4168,  40.4241,  40.4277,
                  40.4424,  40.4791,  40.5121,  40.5524,  40.5891,  40.6221,
                  40.6331,  40.622 ,  40.6587,  40.6881,  40.7211,  40.732 ,
                  40.7137,  40.6732,  40.6328,  40.6218,  40.6364,  40.6364,
                  40.629 ,  40.6216,  40.6252,  40.6252,  40.6251,  40.6324,
                  40.636 ,  40.647 ,  40.6543,  40.6689,  40.6763,  40.6909,
                  40.7055,  40.7312,  40.7605,  40.8009,  40.8266,  40.8779,
                  40.8485,  40.8264,  40.8117,  40.797 ,  40.7896,  40.7748,
                  40.7601,  40.7453,  40.7343,  40.7269,  40.6864,  40.6497,
                  40.6056,  40.5578,  40.5137,  40.4623,  40.4182,  40.3704,
                  40.2712,  40.2124,  40.139 ,  40.0655,  39.9884,  39.9112,
                  39.8561,  39.779 ,  39.7276,  39.6688,  39.61  ,  39.5586,
                  39.4961,  39.4227,  39.3712,  39.3124,  39.2426,  39.1728,
                  39.103 ,  39.0333,  38.9635,  38.857 ,  38.7945,  38.7321,
                  38.7651,  38.8275,  38.9009,  38.9707,  39.0221,  39.0735,
                  39.0147,  38.9706,  38.9118,  38.8567,  38.8494,  38.8016,
                  38.7245,  38.6364,  38.5923,  38.5299,  38.4711])

t_hm = np.array([ 0.00702082,  0.00744214,  0.00785336,  0.00821324,  0.00855123,
                  0.00898331,  0.00935293,  0.00982547,  0.010322  ,  0.0107949 ,
                  0.0113912 ,  0.0120205 ,  0.0125713 ,  0.0129717 ,  0.013505  ,
                  0.0138105 ,  0.0143141 ,  0.014836  ,  0.0153769 ,  0.0157954 ,
                  0.0162982 ,  0.0169682 ,  0.0173523 ,  0.0180663 ,  0.0192364 ,
                  0.0210408 ,  0.022707  ,  0.0246155 ,  0.0268046 ,  0.0289276 ,
                  0.0309399 ,  0.0333904 ,  0.0360352 ,  0.0383699 ,  0.0410401 ,
                  0.0442912 ,  0.0475857 ,  0.0513545 ,  0.0537085 ,  0.0561698 ,
                  0.0582198 ,  0.0606155 ,  0.0639636 ,  0.0678013 ,  0.070591  ,
                  0.0741574 ,  0.0789607 ,  0.0844536 ,  0.0899232 ,  0.0927863 ,
                  0.0953134 ,  0.0974692 ,  0.102394  ,  0.110009  ,  0.124722  ,
                  0.138891  ,  0.151239  ,  0.163949  ,  0.174566  ,  0.181748  ,
                  0.188377  ,  0.196125  ,  0.20511   ,  0.215468  ,  0.225341  ,
                  0.238856  ,  0.253181  ,  0.264785  ,  0.276911  ,  0.285722  ,
                  0.293496  ,  0.301485  ,  0.313878  ,  0.323867  ,  0.335682  ,
                  0.344823  ,  0.354214  ,  0.365496  ,  0.377128  ,  0.382213  ,
                  0.390865  ,  0.394366  ,  0.403288  ])

L_hm = np.array([42.4149,  42.3929,  42.3561,  42.3084,  42.2680,  42.2239,
                 42.1761,  42.1247,  42.0843,  42.0292,  41.9741,  41.9189,
                 41.8712,  41.8161,  41.7279,  41.6582,  41.5847,  41.5002,
                 41.4158,  41.3460,  41.2725,  41.1807,  41.1330,  41.0852,
                 41.0558,  41.0227,  40.9933,  40.9786,  40.9748,  40.9601,
                 40.9454,  40.9307,  40.9270,  40.9233,  40.9416,  40.9452,
                 40.9452,  40.9305,  40.8974,  40.8533,  40.8166,  40.7725,
                 40.7100,  40.6733,  40.6255,  40.5741,  40.5484,  40.5336,
                 40.5005,  40.4381,  40.3940,  40.3206,  40.2765,  40.2544,
                 40.2176,  40.1845,  40.1515,  40.1184,  40.0779,  40.0302,
                 39.9604,  39.8943,  39.8282,  39.7510,  39.6996,  39.6335,
                 39.5673,  39.5269,  39.4425,  39.3543,  39.2735,  39.2074,
                 39.1119,  39.0311,  38.9723,  38.9209,  38.8768,  38.8291,
                 38.7483,  38.6712,  38.6198,  38.5537,  38.4875])

# Conversion factor table for the two metallicity bins as function of age

t_zlo = np.array([  7.02778000e-03,   7.58466000e-03,   8.11267000e-03,
                    8.87417000e-03,   9.57748000e-03,   1.03365000e-02,
                    1.10561000e-02,   1.19322000e-02,   1.28780000e-02,
                    1.40241000e-02,   1.48664000e-02,   1.56890000e-02,
                    1.66316000e-02,   1.77094000e-02,   1.93713000e-02,
                    2.09997000e-02,   2.25618000e-02,   2.47893000e-02,
                    2.77295000e-02,   3.04671000e-02,   3.31760000e-02,
                    3.71109000e-02,   4.09578000e-02,   4.50014000e-02,
                    4.94441000e-02,   5.53092000e-02,   6.29893000e-02,
                    7.07782000e-02,   7.88196000e-02,   8.85651000e-02,
                    9.68742000e-02,   1.08366000e-01,   1.20139000e-01,
                    1.34392000e-01,   1.50335000e-01,   1.66668000e-01,
                    1.82307000e-01,   2.05773000e-01,   2.29156000e-01,
                    2.54054000e-01,   2.86758000e-01,   3.19346000e-01,
                    3.44638000e-01,   3.71931000e-01,   4.03191000e-01,
                    4.51011000e-01,   4.86721000e-01,   5.34769000e-01,
                    6.00896000e-01,   6.78237000e-01,   7.41869000e-01,
                    8.15124000e-01,   9.03677000e-01,   1.03380000e+00,
                    1.15643000e+00,   1.28205000e+00,   1.43413000e+00,
                    1.57571000e+00,   1.75473000e+00,   1.95411000e+00,
                    2.13749000e+00,   2.33809000e+00,   2.62728000e+00,
                    2.87384000e+00,   3.15770000e+00,   3.36236000e+00,
                    3.69448000e+00,   3.98715000e+00,   4.40052000e+00,
                    4.90049000e+00,   5.48170000e+00,   5.96907000e+00,
                    6.35567000e+00,   6.95192000e+00,   7.40241000e+00,
                    7.91716000e+00,   8.13277000e+00,   8.35428000e+00,
                    8.50503000e+00,   8.54300000e+00,   8.73639000e+00,
                    9.17777000e+00,   9.77214000e+00,   1.05459000e+01,
                    1.11286000e+01,   1.14326000e+01,   1.18511000e+01,
                    1.21749000e+01,   1.25076000e+01,   1.28491000e+01,
                    1.32589000e+01,   1.38664000e+01])

t_zhi = np.array([  7.02148000e-03,   7.61183000e-03,   8.32627000e-03,
                    9.14870000e-03,   9.91802000e-03,   1.08489000e-02,
                    1.19206000e-02,   1.30395000e-02,   1.39474000e-02,
                    1.47854000e-02,   1.59575000e-02,   1.70683000e-02,
                    1.84210000e-02,   1.97034000e-02,   2.16496000e-02,
                    2.40021000e-02,   2.68496000e-02,   2.89766000e-02,
                    3.27061000e-02,   3.60968000e-02,   3.98390000e-02,
                    4.35778000e-02,   4.76667000e-02,   5.03012000e-02,
                    5.42843000e-02,   5.83217000e-02,   6.29415000e-02,
                    6.70211000e-02,   7.16862000e-02,   7.66756000e-02,
                    8.27509000e-02,   9.01118000e-02,   9.81270000e-02,
                    1.06374000e-01,   1.16876000e-01,   1.28417000e-01,
                    1.39837000e-01,   1.56425000e-01,   1.70336000e-01,
                    1.88839000e-01,   2.14100000e-01,   2.37360000e-01,
                    2.57314000e-01,   2.87840000e-01,   3.13443000e-01,
                    3.39789000e-01,   3.71668000e-01,   4.23279000e-01,
                    4.52737000e-01,   4.97450000e-01,   5.53983000e-01,
                    6.16926000e-01,   6.93216000e-01,   8.00184000e-01,
                    9.19528000e-01,   1.04723000e+00,   1.15582000e+00,
                    1.29293000e+00,   1.45933000e+00,   1.62516000e+00,
                    1.83434000e+00,   2.01549000e+00,   2.23448000e+00,
                    2.43326000e+00,   2.64974000e+00,   2.82147000e+00,
                    3.01791000e+00,   3.21354000e+00,   3.39135000e+00,
                    3.59505000e+00,   3.77692000e+00,   4.05790000e+00,
                    4.41878000e+00,   4.41878000e+00,   4.70496000e+00,
                    4.87670000e+00,   5.03191000e+00,   5.21547000e+00,
                    5.40576000e+00,   5.83388000e+00,   6.35286000e+00,
                    6.91811000e+00,   7.33357000e+00,   7.98547000e+00,
                    8.35132000e+00,   8.61721000e+00,   8.81209000e+00,
                    9.29885000e+00,   9.72481000e+00,   1.05896000e+01,
                    1.13780000e+01,   1.21161000e+01,   1.26723000e+01,
                    1.34336000e+01,   1.38609000e+01])

f_zhi = np.array([ 0.0429098,  0.0462239,  0.0536475,  0.0622628,  0.071194 ,
                   0.0838693,  0.0988003,  0.118139 ,  0.149953 ,  0.202036 ,
                   0.268171 ,  0.325494 ,  0.395064 ,  0.479512 ,  0.556518 ,
                   0.608478 ,  0.626755 ,  0.645628 ,  0.635913 ,  0.607968 ,
                   0.572648 ,  0.626132 ,  0.626016 ,  0.580973 ,  0.507927 ,
                   0.485633 ,  0.485557 ,  0.54701  ,  0.63489  ,  0.715235 ,
                   0.793807 ,  0.867954 ,  0.921141 ,  0.829698 ,  0.781503 ,
                   0.781352 ,  0.735979 ,  0.714189 ,  0.714064 ,  0.652811 ,
                   0.614851 ,  0.614721 ,  0.623854 ,  0.633081 ,  0.671875 ,
                   0.642373 ,  0.614156 ,  0.578439 ,  0.632492 ,  0.691559 ,
                   0.756121 ,  0.73374  ,  0.744586 ,  0.766901 ,  0.80176  ,
                   0.801546 ,  0.825645 ,  0.837857 ,  0.801    ,  0.825067 ,
                   0.849834 ,  0.888544 ,  0.942958 ,  1.06225  ,  1.23286  ,
                   1.36833  ,  1.71111  ,  2.01586  ,  2.55878  ,  3.1058   ,
                   3.55154  ,  3.71348  ,  3.65788  ,  3.65788  ,  3.19806  ,
                   2.83821  ,  2.07498  ,  1.63441  ,  1.32636  ,  1.23087  ,
                   1.38658  ,  1.70819  ,  1.95332  ,  1.585    ,  1.22996  ,
                   0.954477 ,  0.697819 ,  0.549634 ,  0.413984 ,  0.395803 ,
                   0.557661 ,  0.821672 ,  1.02755  ,  1.32388  ,  0.885034 ])


f_zlo = np.array([  5.38123,   6.15319,   7.46849,   8.93035,  11.0019 ,  13.5541 ,
                    15.9679 ,  19.094  ,  24.5995 ,  32.6512 ,  39.6314 ,  49.5606 ,
                    64.8124 ,  70.8696 ,  75.2117 ,  75.1992 ,  75.1881 ,  69.7724 ,
                    62.842  ,  56.602  ,  50.9825 ,  45.2388 ,  40.1433 ,  36.7004 ,
                    33.0562 ,  31.6027 ,  28.8902 ,  28.4557 ,  26.8021 ,  24.8704 ,
                    23.7779 ,  23.4206 ,  23.4156 ,  23.4102 ,  23.7565 ,  23.7515 ,
                    24.1039 ,  25.579  ,  26.7434 ,  27.5472 ,  30.1181 ,  32.9298 ,
                    31.0182 ,  27.9392 ,  25.9275 ,  22.666  ,  19.2339 ,  16.815  ,
                    16.5621 ,  16.558  ,  15.8306 ,  15.8276 ,  15.3593 ,  15.1277 ,
                    14.4626 ,  13.6222 ,  12.8305 ,  11.5565 ,  10.7238 ,  10.7214 ,
                    11.5492 ,  12.6279 ,  14.4382 ,  15.7867 ,  18.0507 ,  20.3352 ,
                    23.6009 ,  24.3115 ,  23.2432 ,  22.2214 ,  19.426  ,  16.9833 ,
                    15.0714 ,  13.986  ,  14.624  ,  12.0452 ,   9.92198,   8.2958 ,
                    6.24874,   5.30334,   4.36855,   3.54507,   3.05356,   2.55283,
                    2.23196,   2.83325,   4.23753,   5.62521,   7.46732,   9.20045,
                    9.19985,   7.46574])

# Two metallicity bins (solar units)

Z_lo = 0.1
Z_hi = 1.5

# End Fragos et al. 2013 arrays

# Function to calculate the scale factor for a power
# law with F = K*E**-alpha (K in units of ct/s/keV)
def get_scale_factor(ind, emin, emax):
    if ind == 2.0:
        k = np.log(emax/emin)
    else:
        k = (emax**(2.0-ind)-emin**(2.0-ind))/(2.0-ind)
    return keV_per_erg/k

# Function to convert between two different energy
# bands for a single power law
def convert_bands(ind, emin_a, emax_a, emin_b, emax_b):
    if ind == 2.0:
        k = np.log(emax_a/emin_a)
        k /=  np.log(emax_b/emin_b)
    else:
        k = (emax_a**(2.0-ind)-emin_a**(2.0-ind))
        k /= (emax_b**(2.0-ind)-emin_b**(2.0-ind))
    return k

# Spectral indices for both types of XRBs

alpha_lmxb = 1.56
alpha_hmxb = 2.0

# Energy bands for luminosities in XRB 
# distribution functions

emin_lmxb = 0.5
emax_lmxb = 8.0

emin_hmxb = 2.0
emax_hmxb = 10.0

emin_lum = 2.0
emax_lum = 10.0

bolometric_correction = 0.3

# Range of luminosities common to both types of XRBs

Lmin = 1.0e-3*bolometric_correction
Lcut = 500.0*bolometric_correction
nbins = 1000
Lbins = np.logspace(np.log10(Lmin), np.log10(Lcut), nbins+1)
logLbins = np.log10(Lbins)
logLmid = 0.5*(logLbins[1:]+logLbins[:-1])

# LMXB distribution function from Gilfanov 2004

alpha1 = 1.0
alpha2 = 2.0
alpha3 = 5.0

# The luminosities from Gilfanov 2004 are
# in the 0.5-8 keV band. Here we convert
# to 2-10 keV assuming our spectral index
# so that both LMXBs and HMXBs have the 
# same band

kappa = convert_bands(alpha_lmxb, emin_hmxb, emax_hmxb,
                      emin_lmxb, emax_lmxb)

Lb1 = 0.2*kappa
Lb2 = 5.0*kappa

K1 = 1.0
K2 = K1*(Lb1/Lb2)**alpha2
K3 = K2*(Lb2/Lcut)**alpha3

C1 = Lb1*K1
C2 = Lb2*K2/(1.0-alpha2)
C3 = Lcut*K3/(1.0-alpha3)

D2 = C2*(Lb1/Lb2)**(1.0-alpha2)
D3 = C3*(Lb2/Lcut)**(1.0-alpha3)

I1 = C1*np.log(Lb1/Lmin)
I2 = C2 - D2 + I1
I3 = C3 - D3 + I2

def lmxb_cdf(L):
    if L < Lb1:
        N = C1*np.log(L/Lmin)
    elif Lb1 <= L < Lb2:
        N = C2*(L/Lb2)**(1.0-alpha2) - D2 + I1
    elif Lb2 <= L < Lcut:
        N = C3*(L/Lcut)**(1.0-alpha3) - D3 + I2
    else:
        N = I3
    return N

def lmxb_pdf(L):
    if L < Lb1:
        dNdL = K1*(L/Lb1)**-alpha1
    elif Lb1 <= L < Lb2:
        dNdL = K2*(L/Lb2)**-alpha2
    elif Lb2 <= L < Lcut:
        dNdL = K3*(L/Lcut)**-alpha3
    else:
        dNdL = 0.0
    return dNdL

# HMXB distribution function from Mineo et al. 2012

gamma1 = 1.58
gamma2 = 2.73
Lb = 110.0
A = Lb**(gamma2-gamma1)

E1 = 1.0/(1.0-gamma1)
E2 = A/(1.0-gamma2)

F1 = E1*Lmin**(1.0-gamma1)
F2 = E2*Lb**(1.0-gamma2)

J1 = E1*Lb**(1.0-gamma1) - F1
J2 = E2*Lcut**(1.0-gamma2) - F2 + J1

def hmxb_cdf(L):
    if L < Lb:
        N = E1*L**(1.0-gamma1) - F1
    elif Lb <= L < Lcut:
        N = E2*L**(1.0-gamma2) - F2 + J1
    else:
        N = J2
    return N

def hmxb_pdf(L):
    if L < Lb:
        dNdL = L**-gamma1
    elif Lb <= L < Lcut:
        dNdL = A*L**-gamma2
    else:
        dNdL = 0.0
    return dNdL

# Function to compute number of XRBs from
# distribution function

def make_xrbs(Ls, Lfunc, Nfunc, prng):
    dL = lambda L: L*Lfunc(L)
    Ltot = quad(dL, Lmin, Lcut)[0]
    K = Ls.v / (Ltot*1.0e38)
    Ntot = K*Nfunc(Lcut)
    N = prng.poisson(lam=Ntot)
    return N, Ntot

def make_xrb_particles(data_source, metallicity_field, age_field,
                       scale_length, output_lums=None, prng=None):
    r"""
    This routine generates an in-memory dataset composed of X-ray binary particles
    from an input data source containing star particles. 

    Parameters
    ----------
    data_source : :class:`~yt.data_objects.data_containers.YTSelectionContainer`
        The yt data source to obtain the data from, such as a sphere, box, disk, 
        etc.
    metallicity_field : string or (type, name) field tuple
        The stellar metallicity field.
    age_field : string or (type, name) field tuple
        The stellar age field. Must be in some kind of time units. 
    prng : integer or :class:`~numpy.random.RandomState` object 
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is to use the :mod:`numpy.random` module.
    """
    prng = parse_prng(prng)

    ds = data_source.ds

    ptype = data_source._determine_fields(metallicity_field)[0][0]

    t = data_source[age_field].to("Gyr").d
    Z = data_source[metallicity_field].to("Zsun").d
    m = data_source[(ptype, "particle_mass")].to("Msun")

    npart = t.size

    scale_field = None
    if isinstance(scale_length, tuple):
        if isinstance(scale_length[0], string_types):
            scale_field = scale_length
    elif isinstance(scale_length, string_types):
        scale_field = (ptype, scale_length)

    if scale_field is None:
        if isinstance(scale_length, tuple):
            scale = YTArray([scale_length[0]]*npart, scale_length[1])
        elif isinstance(scale_length, YTQuantity):
            scale = YTArray([scale_length]*npart)
        else:
            scale = YTArray([scale_length[0]]*npart, "kpc")
    else:
        scale = data_source[scale_length]

    scale = scale.to('kpc').d

    l_l = np.interp(t, t_lm, L_lm, left=0.0, right=0.0)
    l_h = np.interp(t, t_hm, L_hm, left=0.0, right=0.0)

    l_l[l_l > 0.0] = 10**l_l[l_l > 0.0]
    l_h[l_h > 0.0] = 10**l_h[l_h > 0.0]

    l_l = YTArray(l_l, "erg/s/(1.0e10*Msun)")
    l_h = YTArray(l_h, "erg/s/(1.0e10*Msun)")

    f_zh = np.interp(t, t_zhi, f_zhi, left=0.0, right=0.0)
    f_zl = np.interp(t, t_zlo, f_zlo, left=0.0, right=0.0)

    l_l[Z >= Z_hi] *= f_zh[Z >= Z_hi]
    l_l[Z <= Z_lo] *= f_zl[Z <= Z_lo]

    l_h[Z >= Z_hi] *= f_zh[Z >= Z_hi]
    l_h[Z <= Z_lo] *= f_zl[Z <= Z_lo]

    l_l *= m
    l_h *= m

    # These are bolometric luminosities. Now convert them
    # to the emin_lum-emax_lum keV band assuming the 
    # "low-hard" state conversion factor uncorrected for
    # absorption in the 2-10 keV band from Table 2 of 
    # Fragos et al 2013.

    l_l *= bolometric_correction
    l_h *= bolometric_correction

    l_l.convert_to_units("erg/s")
    l_h.convert_to_units("erg/s")

    N_l, lam_l = make_xrbs(l_l, lmxb_pdf, lmxb_cdf, prng)
    N_h, lam_h = make_xrbs(l_h, hmxb_pdf, hmxb_cdf, prng)

    if output_lums is not None:
        np.savetxt("%s_lmxb.dat" % output_lums, 
                   np.transpose([l_l/bolometric_correction, lam_l, N_l]),
                   delimiter="\t")
        np.savetxt("%s_hmxb.dat" % output_lums, 
                   np.transpose([l_h/bolometric_correction, lam_h, N_h]),
                   delimiter="\t")

    N_all = (N_l+N_h).sum()

    mylog.info("Number of low-mass X-ray binaries: %s" % N_l.sum())
    mylog.info("Number of high-mass X-ray binaries: %s" % N_h.sum())

    if N_all == 0:
        raise RuntimeError("There are no X-ray binaries to generate!")

    # Compute conversion factors from luminosity to count rate

    lmxb_factor = get_scale_factor(alpha_lmxb, emin_lum, emax_lum)
    hmxb_factor = get_scale_factor(alpha_hmxb, emin_lum, emax_lum)

    xp = []
    yp = []
    zp = []
    vxp = []
    vyp = []
    vzp = []
    lp = []
    rp = []
    ap = []

    if N_l.sum() > 0:

        F_l = np.zeros(nbins+1)
        for i in range(1, nbins+1):
            F_l[i] = lmxb_cdf(Lbins[i]) 
        F_l /= F_l[-1]
        invcdf_l = InterpolatedUnivariateSpline(F_l, logLbins)

        for i, n in enumerate(N_l):
            if n > 0:
                randvec = prng.uniform(size=n)
                l = YTArray(10**invcdf_l(randvec)*1.0e38, "erg/s")
                r = YTArray(l.v*lmxb_factor, "photons/s/keV")
                # Now convert output luminosities back to bolometric
                l /= bolometric_correction
                x = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                y = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                z = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                x += data_source[ptype, "particle_position_x"][i].to("kpc")
                y += data_source[ptype, "particle_position_y"][i].to("kpc")
                z += data_source[ptype, "particle_position_z"][i].to("kpc")
                vx = YTArray([data_source[ptype, "particle_velocity_x"][i]]*n).to('km/s')
                vy = YTArray([data_source[ptype, "particle_velocity_y"][i]]*n).to('km/s')
                vz = YTArray([data_source[ptype, "particle_velocity_z"][i]]*n).to('km/s')
                xp.append(x)
                yp.append(y)
                zp.append(z)
                vxp.append(vx)
                vyp.append(vy)
                vzp.append(vz)
                lp.append(l)
                rp.append(r)
                ap.append(np.array([alpha_lmxb]*n))

    if N_h.sum() > 0:

        F_h = np.zeros(nbins+1)
        for i in range(1, nbins+1):
            F_h[i] = hmxb_cdf(Lbins[i])
        F_h /= F_h[-1]
        invcdf_h = InterpolatedUnivariateSpline(F_h, logLbins)

        for i, n in enumerate(N_h):
            if n > 0:
                randvec = prng.uniform(size=n)
                l = YTArray(10**invcdf_h(randvec)*1.0e38, "erg/s")
                r = YTArray(l.v*hmxb_factor, "photons/s/keV")
                # Now convert output luminosities back to bolometric
                l /= bolometric_correction
                x = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                y = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                z = YTArray(prng.normal(scale=scale[i], size=n), "kpc")
                x += data_source[ptype, "particle_position_x"][i].to("kpc")
                y += data_source[ptype, "particle_position_y"][i].to("kpc")
                z += data_source[ptype, "particle_position_z"][i].to("kpc")
                vx = YTArray([data_source[ptype, "particle_velocity_x"][i]]*n).to('km/s')
                vy = YTArray([data_source[ptype, "particle_velocity_y"][i]]*n).to('km/s')
                vz = YTArray([data_source[ptype, "particle_velocity_z"][i]]*n).to('km/s')
                xp.append(x)
                yp.append(y)
                zp.append(z)
                vxp.append(vx)
                vyp.append(vy)
                vzp.append(vz)
                lp.append(l)
                rp.append(r)
                ap.append(np.array([alpha_hmxb]*n))

    xp = uconcatenate(xp)
    yp = uconcatenate(yp)
    zp = uconcatenate(zp)
    vxp = uconcatenate(vxp)
    vyp = uconcatenate(vyp)
    vzp = uconcatenate(vzp)
    lp = uconcatenate(lp)
    rp = uconcatenate(rp)
    ap = uconcatenate(ap)

    data = {("io", "particle_position_x"): xp,
            ("io", "particle_position_y"): yp,
            ("io", "particle_position_z"): zp,
            ("io", "particle_velocity_x"): vxp,
            ("io", "particle_velocity_y"): vyp,
            ("io", "particle_velocity_z"): vzp,
            ("io", "particle_luminosity"): lp,
            ("io", "particle_count_rate"): rp,
            ("io", "particle_index"): ap}

    dle = ds.domain_left_edge.to("kpc").v
    dre = ds.domain_right_edge.to("kpc").v

    bbox = np.array([[dle[i], dre[i]] for i in range(3)])

    new_ds = load_particles(data, bbox=bbox, length_unit="kpc", 
                            time_unit="Myr", mass_unit="Msun", 
                            velocity_unit="km/s")

    return new_ds

def make_xrb_photons(ds, redshift, area, exp_time, emin, emax, 
                     center="c", cosmology=None, prng=None):
    dd = ds.all_data()
    e0 = (1.0, "keV")
    prng = parse_prng(prng)
    xrb_model = PowerLawSourceModel(e0, emin, emax, 
                                    ("io", "particle_count_rate"),
                                    ("io", "particle_index"), prng=prng)
    photons = PhotonList.from_data_source(dd, redshift, area, exp_time,
                                          xrb_model, center=center,
                                          cosmology=cosmology)
    return photons