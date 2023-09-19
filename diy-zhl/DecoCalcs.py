#!/usr/bin/python -u
#
# (K) Copy Rites Reversed: reuse what you like (but give credit)
#
# Credits:
#
# Mark Powell's "Deco for Divers"
# Erik Baker's papers: "Decolessons", "Understanding M-values", and "Clearing up the confusion about Deep Stops" in particular
# Buhlmann's "Decompression - Decompression Sickness", English edition
# Several open-source implementations, most notably Subsurface software (and people, Robert in particular)
# Plenty of other on-line sources, e.g. Stuart Morrison's "DIY Decompresion"
#
# The goal here is "by the book" implementation to use for learning this stuff.
#

import re, sys
import pprint

import Constants, Equations, Utilities

class TimeParser( object ):

    _Pattern_ = re.compile( r"(?:(\d{1,2}):)?(\d{1,2}):(\d{1,2})" )

    # helper function: takes human-readable time string like "1:30" and returns minutes: 1.5
    @classmethod
    def parse( cls, t = "0:0" ) :
        if t is None : return 0
        m = cls._Pattern_.search( str( t ).strip() )
        if not m : raise Exception( "Invalid time string %s" % (t,) )
        rc = float( m.group( 2 ) ) + float( m.group( 3 ) ) / 60.
        if m.group( 1 ) is not None :
            rc += 60. * float( m.group( 1 ) )
        return round( rc, 1 )

class Dive( object ) :

    def __init__( self, GFHi, use_4min_not_5min = False ) :

        # air, sea level, USN RQ.
        self._T = 0
        self._S = 1.0       # Surface pressure (const)
        self._P = self._S   # Current pressure (var)
        self._Q = 0.79
        self._RQ = 0.9
        self._GFHi = GFHi

        self.TC = Utilities.TableValues.fetchZHL( 16, 'N', use_4min_not_5min = use_4min_not_5min )
        self.TC['P'] = Equations.palv( Pamb = self._P, Q = self._Q, RQ = self._RQ )
        self.TC['C'] = Equations.buhlmann( self.TC['P'], self.TC['a'], self.TC['b'], gf = self._GFHi )

    def updateSchreinerBuhlmann( self, Palv, R, t ):
        kay = Equations.kay( self.TC['t'] )
        pres = Equations.schreiner( Pi = self.TC['P'], Palv = Palv, R = R, t = t, k = kay )
        ceil = Equations.buhlmann( Pn = pres, an = self.TC['a'], bn = self.TC['b'], gf = self._GFHi )
        ndl = Equations.ndl( Palv = Palv, t = t, R = R, M0 = self.TC['M0'], k = kay )
        return pres, ceil, ndl

    def _calcNDLTimeTC_( self, i, start = 0 ):
        Palv = Equations.palv( Pamb = self._P, Q = self._Q, RQ = self._RQ )
        calc_ceil = lambda t: self._calcNDLCeilingTC_( i, Palv, 0., t )
        ceil_below_water = lambda c: c > 0.1
        return Utilities.BinarySearch.hop( calc_ceil, ceil_below_water, start, step = 1, eps = 1e-3 )

    # new_depth is new depth in 0.1 bar / depth in meters
    # timestr is time as [hours:]minutes:seconds string. *it is the total elapsed* time
    def segment( self, new_depth = 0.0, newtimestr = "1:0" ) :
        assert new_depth >= 0.0
        # if new_depth == 0. : newP = self._P
        newP = round( self._S + new_depth / 10, 1 )
        t = TimeParser.parse( newtimestr ) - self._T

        Palv = Equations.palv( Pamb = self._P, Q = self._Q, RQ = self._RQ )
        R = Equations.dP_dt( d0 = self._P, dt = newP, t = t, Q = self._Q )
        self.TC['P'], self.TC['C'], self.TC['L'] = self.updateSchreinerBuhlmann( Palv, R, t )

        self._P = newP
        self._T += t
