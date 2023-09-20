
import math

#  Subsurface :: deco.c :: WV_PRESSURE
#  Inspired gas loading equations depend on the partial pressure of inert gas in the alveolar.
#  - RQ_frac  = (1 - RQ) / RQ
#  - WV       = PP_H20 - RQ_frac * PP_CO2
#  - P_alv    = Q * (P_amb - WV)
#  Output:
#  - WV       effective water vapour pressure
#  - P_alv    alveolar partial pressure of inert gas
#  Inputs:
#  - P_amb    ambient pressure
#  - PP_H2O   water vapour partial pressure = ~0.0627 bar
#  - PP_CO2   carbon dioxide partial pressure = ~0.0534 bar
#  - RQ       respiratory quotient (O2 consumption / CO2 production)
#  - Q        fraction of inert gas
#  
#  Buhlmann ignored the contribution of CO2 (i.e. Rq = 1.), whereas Schreiner adopted Rq = 0.8
#  WV_Buhlmann  = PP_H2O = 0.0627 bar
#  WV_Schreiner = 0.0627 - (1 - 0.8) / RQ * 0.0534 = 0.0493 bar
#  Buhlmann calculations use the Buhlmann value, VPM-B calculations use the Schreiner value.

PP = { 'H2O': 0.0627, 'CO2': 0.0534, }
RQ = { 'schreiner': 0.8, 'usnavy': 0.9, 'buhlmann': 1., }
# Schreiner is the most conservative, Buhlmann is the most aggressive
# https://www.shearwater.com/wp-content/uploads/2012/08/Introductory-Deco-Lessons.pdf

def palv( P_ambient, inert_gas_frac, RQ_name = 'buhlmann' ) :
    assert RQ_name in RQ
    RQ_coef = RQ [RQ_name]
    RQ_frac = (1. - RQ_coef) / RQ_coef
    P_eff_vapor = PP ['H2O'] - RQ_frac * PP ['CO2']
    return inert_gas_frac * (P_ambient - P_eff_vapor)

# return k: constant for tissue compartment (min^-1)
# Th : tissue compartment half-time in minutes
# for 5-minute compartment it's 0.8452
def kay( Th ) :
    assert Th > 0. if isinstance( Th, float ) else all( Th > 0. )
    return math.log( 2 ) / Th

# return rate of pressure change in bar/min
# d0 : start pressure, bar
# dt : end pressure, bar
# t : time, min
# Q : fraction of inert gas (same Q as in palv()
# Formerly called arr
def dP_dt( d0 = 1.0, dt = 1.0, t = 1, Q = 0.79 ) :
    assert float( t ) > 0.0
    dP = (float( dt ) - float( d0 )) / float( t )
    rc = dP * float( Q )
    return round( rc, 4 )

# Schreiner equation
# Palv + R * (t - 1/k) - (Palv - Pi - R/k) * e^(-k * t)
#
# returns pressure in tissue compartment after time t at depth Pa & dP
#
# Pi: initial pressure of inert gas in tissue (bar)
# Palv: initial pressure of inert gas in the lungs (bar, output of palv())
# t: time (minutes)\n",
# R: rate of pressure change (output of dP_dt()),
# k: gas decay constant (output of kay()
#
# (Intermediate variables b/c I was playing with rounding)
#

def schreiner( Pi = 0.7451, Palv = 0.7451, t = 1, R = 0, k = 0.1386, verbose = False ) :

    assert k != 0. if isinstance( k, float ) else all( k != 0. )
    x1 = R * (t - 1. / k)
    x2 = Palv - Pi - R / k
    x3 = math.e ** ( -k * t )
    rc = Palv + x1 - x2 * x3
    if verbose : print( "x1: %f, x2: %f, x3: %f, rc: %f\n" % (x1, x2, x3, rc,) )
    return rc

# M-value: workman to buhlmann
# P is ambient pressure in bar
# returns pair ( a, b )
#
# TODO: add GF?
#
def m_w2b( M0 = 2.9624, dM = 1.7928, P = 1 ) :
    assert (dM > 1.) if isinstance (dM, float) else all (dM > 1.)
    a = M0 - dM * P
    b = 1. / dM
    return (a, b)

# M-value: buhlmann to workman
# returns pair ( M0, dM )
#
def m_b2w( a = 1.1696, b = 0.5578, P = 1 ) :
    assert (b > 0.) if isinstance (b, float) else all (b > 0.)
    M0 = a + P / b
    dM = 1. / b
    return (M0, dM)

class MValueConversion (object):

    @staticmethod
    def getWorkmanFromBuhlmann(
            zero_intercept_A     : float,
            inverse_slope_B      : float,
            sea_level_pressure_P : float,
        ) \
        -> (float, float):
            slope_dM = 1. / inverse_slope_B
            one_intercept_M0 = zero_intercept_A + slope_dM * sea_level_pressure_P
            return (one_intercept_M0, slope_dM)

    @staticmethod
    def getBuhlmannFromWorkman(
        one_intercept_M0     : float,
        slope_dM             : float,
        sea_level_pressure_P : float,
    ) \
    -> (float, float):
        inverse_slope_B = 1. / slope_dM
        zero_intercept_A = one_intercept_M0 - slope_dM * sea_level_pressure_P
        return (zero_intercept_A, inverse_slope_B)

# no-stop time by Schreiner
#
# Palv: initial pressure of inert gas in the lungs (bar, output of palv())
# t: time (minutes)
# R: rate of pressure change (output of dP_dt()),
# k: gas decay constant (output of kay()
#  -- same as schreiner()
# M0: surfacing M-value as per Workman
#
def ndl( Palv = 0.7451, M0 = 2.9624, t = 0, R = 0, k = 0.1386, verbose = False ) :
    # (M0 - Palv - R * (t - 1/k)) * math.e ** (k * t) + Palv - R / k
    assert k != 0. if isinstance( k, float ) else all( k != 0. )
    x1 = M0 - Palv - R * (t - 1. / k)
    x2 = math.e ** (k * t)
    rc = x1 * x2 + Palv - R / k
    return rc

# Buhlman formula with GF and Helium
# returns safe ascent ceiling
#
# Pn is N2 pressure in tissue compartment
# Phe is He pressure in tissue compartment
# an is N2 a coefficient
# bn is N2 b coefficient
# ahe is He a coefficient
# bhe is He b coefficient
# gf is current gradient factor
#
def buhlmann( Pn, an, bn, Phe = 0, ahe = 0, bhe = 0, gf = 1 ) :

    P = Pn + Phe
    assert P != 0. if isinstance( P, float ) else all( P != 0. )
    a = (an * Pn + ahe * Phe) / P
    b = (bn * Pn + bhe * Phe) / P
    num = P - a * gf
    den = gf / b + 1. - gf
    assert den != 0. if isinstance( den, float ) else all( den != 0. )
    rc = num / den
    return rc
