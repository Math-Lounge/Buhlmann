
import numpy as np
import pandas as pd

import Constants

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
        df = pd.DataFrame (dive._TCs).drop ([ 'a', 'b', ], axis = 1)
        df ['compart_num'] = 1 + np.arange (len (dive._TCs))
        df ['time_min'] = dive._T
        return df.rename (columns = { 't': 'half_life_min', 'P': 'pressure_bar', 'C': 'ceiling', })

    def diveUpdate (self):
        df = self.snapFromDiveObj (self.dive)
        self.tissues = pd.concat ([ self.tissues, df, ], ignore_index = True)
