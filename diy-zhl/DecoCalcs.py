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

import Constants, Equations

class TimeParser( object ):

    _Pattern_ = re.compile( r"(?:(\d{1,2}):)?(\d{1,2}):(\d{1,2})" )

    # helper function: takes human-readable time string like "1:30" and returns minutes: 1.5
    @classmethod
    def parse( cls, t = "0:0" ) :
        if t is None : return 0
        m = cls._Pattern_.search( str( t ).strip() )
        if not m : raise Exception( "Invalid time string %s" % (t,) )
        rc = 0.0
        if m.group( 1 ) is not None :
            rc = float( m.group( 1 ) ) * 60.0
        rc += float( m.group( 2 ) )
        rc += float( m.group( 3 ) ) / 60.0
        return round( rc, 1 )

class Dive( object ) :
    
    def __init__( self, GFHi, use_4min_not_5min = False, verbose = False ) :

        self._verbose = bool( verbose )
        
        # air, sea level, USN RQ.
        self._T = 0
        self._S = 1.0       # Surface pressure (const)
        self._P = self._S   # Current pressure (var)
        self._Q = 0.79
        self._RQ = 0.9
        self._GFHi = GFHi
        self._TCs = []

        # starting Pt (same for all TCs)
        sp = Equations.palv( Pamb = self._P, Q = self._Q, RQ = self._RQ )

        ZHL = [ Constants.ZHL16N_5m, Constants.ZHL16N_4m, ] [use_4min_not_5min]
        for tc in ZHL.keys() :
            self._TCs.append( { 
                "t" : ZHL[tc]["t"],
                "a" : ZHL[tc]["a"]["C"],
                "b" : ZHL[tc]["b"],
                "P" : sp,
            } )

        # init. ceiling
        for i in range( len( self._TCs ) ) :
            self._TCs[i]["C"] = self.runBuhlmanTC( i )
            
        if self._verbose : pprint.pprint( self._TCs )
            
    @property
    def compartments( self ) :
        return list( map( lambda TC: TC[ 't' ], self._TCs ) )
    
    @property
    def loadings( self ) :
        return list( map( lambda TC: TC[ 'P' ], self._TCs ) )

    def runBuhlmanTC( self, i ):
        return Equations.buhlmann(
            Pn = self._TCs[i]["P"],
            an = self._TCs[i]["a"],
            bn = self._TCs[i]["b"],
            gf = self._GFHi,
        )

    # newdepth is new depth in 0.1 bar
    # timestr is time as [hours:]minutes:seconds string. *it is the total elapsed* time
    def segment( self, newdepth = 0.0, newtimestr = "1:0" ) :
        assert float( newdepth ) >= 0.0
        if float( newdepth ) == 0.0 :
            newP = self._P
        else :
            newP = round( self._S + float( newdepth ) / 10, 1 )
        t = TimeParser.parse( newtimestr ) - self._T
    
        for i in range( len( self._TCs ) ) :
            Palv = Equations.palv( Pamb = self._P, Q = self._Q, RQ = self._RQ )
            p = Equations.schreiner(
                Pi = self._TCs[i]["P"], Palv = Palv, t = t,
                R = Equations.arr( d0 = self._P, dt = newP, t = t, Q = self._Q ),
                k = Equations.kay( Th = self._TCs[i]["t"] ),
            )
            # ndl(
            #    Palv = Palv, t = t,
            #    R = arr( d0 = self._P, dt = newP, t = t, Q = self._Q ),
            #    k = kay( Th = self._TCs[i]["t"] )
            #)
            self._TCs[i]["P"] = p
            self._TCs[i]["C"] = self.runBuhlmanTC( i )

        self._P = newP
        self._T += t
    
        if self._verbose :
            sys.stdout.write( "* At time %f, P %f:\n" % (self._T, self._P,) )
            pprint.pprint( self._TCs )
