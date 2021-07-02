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
        do_run_pipeline = False
        for filename in methylprep_files:
            if not Path(LOCAL,filename).exists():
                do_run_pipeline = True
                print(f"MUST re-run pipeline on {LOCAL} because files are missing.")
                break
        if do_run_pipeline:
            methylprep.make_pipeline(LOCAL, steps=['all'], exports=['all'], make_sample_sheet=True, save_uncorrected=True)
            # same as CLI -d . --all
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
        print("sesame_mouse_infer vs raw mprep meth/unmeth", mean_meth_diff, mean_unmeth_diff)
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
        print("noob meth/unmeth vs sesame mouse_noob_dye meth/unmeth",mean_meth_diff, mean_unmeth_diff)
        #import pdb;pdb.set_trace()
        #if mean_meth_diff > 0.1 or mean_unmeth_diff > 0.1:
        #    raise AssertionError(f"noob-dye intensities don't match: meth {mean_meth_diff} unmeth {mean_unmeth_diff}")

        # 3 betas
        #ses = pd.read_csv('sesame_mouse_betas.csv').set_index('Unnamed: 0')
        #oses = pd.read_csv('open_sesame_mouse_betas.csv').set_index('Unnamed: 0')
        #combined = pd.concat([ses, oses], axis='columns').rename(columns={'X204879580038_R06C02':'open sesame','x':'sesame'})
        #mprep = pd.read_pickle('beta_values.pkl')
        #combined = pd.concat([ses, mprep], axis='columns').rename(columns={'X204879580038_R06C02':'ses','204879580038_R06C02':'mprep'})
        combined = pd.concat([self.sesame_mouse_betas, self.beta_values], axis='columns').rename(columns={'X204879580038_R06C02':'ses','204879580038_R06C02':'mprep'})
        #import matplotlib.pyplot as plt
        #(combined.ses - combined.mprep).hist(bins=200, range=[-0.05, 0.05])
        #plt.show()
        beta_mean_diff = (combined.ses - combined.mprep).mean()
        print("beta mean diff vs sesame", beta_mean_diff) # actual: 0.00762
        if beta_mean_diff > 0.02:
            raise AssertionError(f"beta_mean_diff exceeds 0.02: (actual: {beta_mean_diff})")
