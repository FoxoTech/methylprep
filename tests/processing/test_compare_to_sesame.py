import methylprep
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt # debug
LOCAL = Path('docs/example_data/mouse')


class TestSesame():

    def test_vs_sesame_mouse(self):
        sample = '204879580038_R06C02'
        sesame_files = [
            #'sesame_mouse_raw.csv',
            'sesame_mouse_infer.csv',
            #'sesame_mouse_noob.csv',
            #'sesame_mouse_poobah.csv',
            'sesame_mouse_dye.csv',
            'sesame_mouse_betas.csv',
        ]

        for filename in sesame_files:
            attrib = filename.split('.')[0]
            df = pd.read_csv(Path(LOCAL,filename)).set_index('Unnamed: 0').sort_index()
            setattr(self, attrib, df)
        methylprep_files = [
            'meth_values.pkl',
            'unmeth_values.pkl',
            'noob_meth_values.pkl',
            'noob_unmeth_values.pkl',
            'beta_values.pkl',
        ]
        for filename in methylprep_files:
            attrib = filename.split('.')[0]
            df = pd.read_pickle(Path(LOCAL,filename))
            setattr(self, attrib, df)

        # 1 compare meth/unmeth vs sesame raw
        # self.sesame_mouse_raw # M and U
        mraw = pd.concat([
            self.meth_values.rename(columns={sample: 'M'}),
            self.unmeth_values.rename(columns={sample: 'U'}),
        ], axis='columns').sort_index()
        combined = mraw.join(self.sesame_mouse_infer, lsuffix='m', rsuffix='s').astype(float)
        mean_meth_diff = (combined.Mm - combined.Ms).mean()
        mean_unmeth_diff = (combined.Um - combined.Us).mean()
        print(mean_meth_diff, mean_unmeth_diff)
        if mean_meth_diff > 0.1 or mean_unmeth_diff > 0.1: # actual 0.0043 0.0671, but histogram shows EXACT match
            raise AssertionError(f"raw inferred-switched intensities don't match: meth {mean_meth_diff} unmeth {mean_unmeth_diff}")

        # 2 compare NOOB meth/unmeth vs sesame raw
        # self.sesame_mouse_raw # M and U
        mdye = pd.concat([
            self.noob_meth_values.rename(columns={sample: 'M'}),
            self.noob_unmeth_values.rename(columns={sample: 'U'}),
        ], axis='columns').sort_index()
        combined = mdye.join(self.sesame_mouse_dye, lsuffix='m', rsuffix='s').astype(float)
        mean_meth_diff = (combined.Mm - combined.Ms).mean()
        mean_unmeth_diff = (combined.Um - combined.Us).mean()
        print(mean_meth_diff, mean_unmeth_diff)
        import pdb;pdb.set_trace()
        if mean_meth_diff > 0.1 or mean_unmeth_diff > 0.1:
            raise AssertionError(f"noob-dye intensities don't match: meth {mean_meth_diff} unmeth {mean_unmeth_diff}")
