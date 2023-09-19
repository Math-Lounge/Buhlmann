
import numpy as np
import pandas as pd

import Constants, Equations

class TableValues (object):

    @staticmethod
    def _fetch_ (model, version, sub_dict_a):
        values = getattr (Constants, model)
        df = pd.DataFrame.from_dict (values).T
        if sub_dict_a:
            df ['a'] = df ['a'].apply (lambda x: x [version])
        for col in df.columns: df [col] = df [col].astype (float)
        return df

    @classmethod
    def fetchZHL (cls, compartments, gas, **kwargs):
        assert gas in [ 'He', 'N', ]
        assert compartments in [ 12, 16, ]
        model = f'ZHL{compartments}{gas}'
        sub_dict_a = (compartments == 16)
        version = [ 'B', 'C', ] [sub_dict_a]
        df = cls._fetch_ (model, version, sub_dict_a)
        if (model == 'ZHL16N'):
            use_4m_not_5m = kwargs.pop ('use_4m_not_5m', True)
            drop = [ 1, 1.1, ] [use_4m_not_5m]
            df = df [~df.index.isin ([drop])]
            df.index = np.floor (df.index).astype (int)
        df ['M0'], df ['dM'] = Equations.MValueConversion.getWorkmanFromBuhlmann (df ['a'], df ['b'], 1.)
        return df

    @classmethod
    def fetchWorkman (cls):
        return cls._fetch_ ('WORKMAN', None, False)

    @classmethod
    def fetchDSAT (cls):
        return cls._fetch_ ('DSAT', None, False)

class TimeSeriesFrame (object):

    def __init__ (self, dive):
        self.dive = dive
        self.tissues = pd.DataFrame ()

    @staticmethod
    def snapFromDiveObj (dive):
        df = dive.TC.drop ([ 'a', 'b', ], axis = 1).copy ()
        df ['compart_num'] = 1 + np.arange (len (dive.TC))
        df ['time_min'] = dive._T
        return df.rename (columns = { 't': 'half_life_min', 'P': 'pressure_bar', 'C': 'ceiling', })

    def diveUpdate (self):
        df = self.snapFromDiveObj (self.dive)
        self.tissues = pd.concat ([ self.tissues, df, ], ignore_index = True)

class BinarySearch (object):

    @staticmethod
    def _increasingSearch_ (func, cond, start, step, eps, max_iter):
        last_X, last_Y = start, np.NaN
        for i in range (max_iter):
            curr_X = start + 2 ** i * step
            curr_Y = func (curr_X)
            if abs (curr_Y - last_Y) < eps:
                raise RuntimeError (f'Failed to converge = {curr_Y}')
            if cond (curr_Y): return (last_X, curr_X, True)
            last_X, last_Y = curr_X, curr_Y
        return (last_X, curr_X, False)

    @staticmethod
    def _decreasingSearch_ (func, cond, LHS, RHS, step):
        assert not cond (LHS), f'LHS value {LHS} should be false'
        assert     cond (RHS), f'RHS value {RHS} should be true'
        while (RHS - LHS) > step:
            mid_X = (LHS + RHS) / 2
            curr_Y = func (mid_X)
            if cond (curr_Y): RHS = mid_X
            else            : LHS = mid_X
        return LHS # last value before ceiling violated

    @classmethod
    def hop (cls, func, cond, start, step = 1, eps = 1e-4, max_iter = 20):
        last_X, curr_X, success = cls._increasingSearch_ (func, cond, start, step, eps, max_iter)
        if not success: raise RuntimeWarning (f'Exceeded max_iter = {max_iter}')
        # Why is this necessary? This checks if the first value exceeded ceiling
        if abs (last_X - curr_X) <= step: return last_X
        return cls._decreasingSearch_ (func, cond, last_X, curr_X, step)
