
import math

# return alveolar inert gas pressure
# with P amb = 1 bar, fraction of inert gas = 0.79, and RQ = 0.9
# this should return 0.79 - 0.0567 = 0.7451 or 0.7452 dep. on where you round it
#
def palv( Pamb = 1, Q = 0.79, RQ = 0.9 ) :
    assert float( RQ ) != 0.0
    vw = float( Pamb ) - 0.0627 + (1.0 - float( RQ )) / float( RQ ) * 0.0534
    return round( vw * float( Q ), 4 )

# return k: constant for tissue compartment (min^-1)
# Th : tissue compartment half-time in minutes
# for 5-minute compartment it's 0.8452
#
def kay( Th = 5 ) :
    assert float( Th ) > 0.0
    return round( math.log( 2 ) / float( Th ), 4 )

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

    assert float( k ) != 0.0
    x1 = float( R ) * (float( t ) - 1.0 / float( k ))
    x2 = float( Palv ) - float( Pi ) - float( R ) / float( k )
    x3 = math.e ** (float( -k ) * float( t ))
    rc = round( float( Palv ) + x1 - x2 * x3, 4 )
    if verbose : print( "x1: %f, x2: %f, x3: %f, rc: %f\n" % (x1, x2, x3, rc,) )
    return round( rc, 4 )

# M-value: workman to buhlmann
# P is ambient pressure in bar
# returns pair ( a, b )
#
# TODO: add GF?
#
def m_w2b( M0 = 2.9624, dM = 1.7928, P = 1 ) :
    assert float( dM ) >= 1.0
    a = float( M0 ) - float( dM ) * float( P )
    b = 1.0 / float( dM )
    return (round( a, 4 ), round( b, 4 ))

# M-value: buhlmann to workman
# returns pair ( M0, dM )
#
def m_b2w( a = 1.1696, b = 0.5578, P = 1 ) :
    assert float( b ) > 0.0
    M0 = float( a ) + float( P ) / float( b )
    dM = 1.0 / float( b )
    return (round( M0, 4 ), round( dM, 4 ))

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
    assert float( k ) != 0.0
    x1 = float( M0 ) - float( Palv ) - float( R ) * (float( t ) - 1.0 / float( k ))
    x2 = math.e ** (float( k ) * float( t ))
    rc = x1 * x2 + float( Palv ) - float( R ) / float( k )
    if verbose : print( "x1: %f, x2: %f, rc: %f\n" % (x1, x2, rc,) )
    return round( rc, 4 )

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

    P = float( Pn ) + float( Phe )
    assert float( P ) != 0.0
    a = (float( an ) * float( Pn ) + float( ahe ) * float( Phe )) / P
    b = (float( bn ) * float( Pn ) + float( bhe ) * float( Phe )) / P
    num = float( P ) - float( a ) * float( gf )
    den = float( gf ) / float( b ) + 1.0 - float( gf )
    assert den != 0.0
    rc = num / den
    return round( rc, 4 )
