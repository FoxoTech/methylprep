# analog of sesame's SigSet class
import pandas as pd
import numpy as np
import methylprep
from pathlib import Path

def sesame_convert(filename, LOCAL_PATH = Path('/Volumes/LEGX/GSE69852/idats/'), drop_rs=True):
    data = pd.read_csv(Path(LOCAL_PATH,filename)).rename(columns={'Unnamed: 0':'IlmnID','M':'meth','U':'unmeth'}).set_index('IlmnID').sort_index().round()
    if drop_rs:
        data = data[~data.index.str.startswith('rs')] # dropping SNP probes, which are included in sesame output
    return data

class Sesame():
    """ used for testing, comparing methylprep with sesame output. """
    def __init__(self, LOCAL_PATH, SESAME_PATH=None, run_tests=True, array_type=methylprep.ArrayType.ILLUMINA_450K ):
        self.LOCAL_PATH = LOCAL_PATH
        self.SESAME_PATH = SESAME_PATH
        self.sample_sheet = methylprep.files.get_sample_sheet(LOCAL_PATH)
        self.manifest = methylprep.Manifest(array_type)

        # raw datasets are the red and green idat probe values, indexed by illumina ids (not cpg000000)
        self.raw_datasets = methylprep.models.raw_dataset.get_raw_datasets(self.sample_sheet)

        # meth datasets have idat values linked to IlmnID probe ids.
        # storing as a dict pair of objects
        meth_datasets = []
        for raw_dataset in self.raw_datasets:
            methylated = methylprep.models.MethylationDataset.methylated(raw_dataset, self.manifest)
            unmethylated = methylprep.models.MethylationDataset.unmethylated(raw_dataset, self.manifest)
            meth_datasets.append({'meth':methylated, 'unmeth':unmethylated})

        # only testing first sample
        raw_meth = meth_datasets[0]['meth'].data_frame.sort_index()
        raw_unmeth = meth_datasets[0]['unmeth'].data_frame.sort_index()

        self.raw_IG = pd.DataFrame(data={
            'meth':raw_meth[(raw_meth['Infinium_Design_Type'] == 'I') & (raw_meth['Color_Channel'] == 'Grn')]['mean_value'].astype(int),
            'unmeth':raw_unmeth[(raw_unmeth['Infinium_Design_Type'] == 'I') & (raw_unmeth['Color_Channel'] == 'Grn')]['mean_value'].astype(int)
            })
        self.raw_IR = pd.DataFrame(data={
            'meth': raw_meth[(raw_meth['Infinium_Design_Type'] == 'I') & (raw_meth['Color_Channel'] == 'Red')]['mean_value'].astype(int),
            'unmeth': raw_unmeth[(raw_unmeth['Infinium_Design_Type'] == 'I') & (raw_unmeth['Color_Channel'] == 'Red')]['mean_value'].astype(int)
            })
        self.raw_II = pd.DataFrame(data={
            'meth': raw_meth[(raw_meth['Infinium_Design_Type'] == 'II')]['mean_value'].astype(int),
            'unmeth': raw_unmeth[(raw_unmeth['Infinium_Design_Type'] == 'II')]['mean_value'].astype(int)
            })

        if run_tests:
            self.compare_with_sesame()

        # minimal processing here, skip all steps. just create a container for testing
        # no noob, no dye, no poobah yet -- future, maybe pass this in as containers, instead
        self.containers = []
        self.snps = []
        for raw_dataset in self.raw_datasets:
            container = methylprep.processing.pipeline.SampleDataContainer(raw_dataset,
                self.manifest,
                retain_uncorrected_probe_intensities=False,
                bit='float32',
                pval=False,
                poobah_decimals=3,
                poobah_sig=0.05,
                do_noob=False,
                quality_mask=False,
                switch_probes=False,
                correct_dye_bias=False,
                debug=False)
            container.process_all()
            self.containers.append(container)

            snp = methylprep.processing.postprocess.one_sample_control_snp(container)
            snp = snp[~snp['snp_meth'].isna()][['snp_meth','snp_unmeth']].rename({'snp_meth': 'meth', 'snp_unmeth':'unmeth'})
            self.snps.append(snp)

        if run_tests:
            self.test_dye_bias()


    def compare_with_sesame(self):
        ses_raw_II = sesame_convert('raw_II.csv', self.SESAME_PATH)
        ses_raw_IG = sesame_convert('raw_IG.csv', self.SESAME_PATH)
        ses_raw_IR = sesame_convert('raw_IR.csv', self.SESAME_PATH)
        print('IG matches:',self.raw_IG.equals(ses_raw_IG))
        print('IR matches:',self.raw_IR.equals(ses_raw_IR))
        print('II matches:',self.raw_II.equals(ses_raw_II))

    def test_dye_bias(self):
        for container in self.containers:
            container.betas_from_raw = pd.DataFrame(container._SampleDataContainer__data_frame['beta_value'].sort_index())

        # run DYE STEP
        self.dye_debug = methylprep.processing.dye_bias.nonlinear_dye_bias_correction(self.containers[0], debug=True)

        self.test_IR1 = pd.read_csv(Path(self.LOCAL_PATH,'test_IR1.csv'))['x']
        self.test_IG1 = pd.read_csv(Path(self.LOCAL_PATH,'test_IG1.csv'))['x']

        self.test_IR2 = pd.read_csv(Path(self.LOCAL_PATH,'test_IR2.csv'))['x']
        self.test_IG2 = pd.read_csv(Path(self.LOCAL_PATH,'test_IG2.csv'))['x']

        self.ses_II = sesame_convert('ses_II.csv', self.SESAME_PATH)
        self.ses_IG = sesame_convert('ses_IG.csv', self.SESAME_PATH)
        self.ses_IR = sesame_convert('ses_IR.csv', self.SESAME_PATH)

        raw_IGrs = sesame_convert('raw_IG.csv', self.SESAME_PATH, drop_rs=False)
        raw_IRrs = sesame_convert('raw_IR.csv', self.SESAME_PATH, drop_rs=False)

        IR0 = pd.concat( [self.containers[0].IR, self.containers[0].snp_IR.rename(columns={'meth':'noob_meth', 'unmeth':'noob_unmeth'})] ).sort_index()
        IG0 = pd.concat( [self.containers[0].IG, self.containers[0].snp_IG.rename(columns={'meth':'noob_meth', 'unmeth':'noob_unmeth'})] ).sort_index()

        print('IG+rs matches:', IR0.equals(raw_IGrs))
        print('IR+rs matches:', IG0.equals(raw_IRrs))

        try:
            print('IR1 match:', self.test_IR1.equals(pd.Series(self.dye_debug['IR1']).astype(int)))
            self.compare(self.test_IR1, self.dye_debug['IR1'])
            print('IG1 match:', self.test_IG1.equals(pd.Series(self.dye_debug['IG1']).astype(int)))
            self.compare(self.test_IG1, self.dye_debug['IG1'])

            self.matrixIR1 = pd.read_csv(Path(LOCAL,'test_matrixIR1.csv'))
            print('matrix(IR1) match:', pd.Series(self.dye_debug['IR1']).astype(int) .equals( self.matrixIR1['V1'] ))
            self.compare( pd.Series(self.dye_debug['IR1']).astype(int), self.matrixIR1['V1'], atol=1)

            self.asvectorIG0 = pd.read_csv(Path(LOCAL,'test_asvectorIG0.csv'))
            print('asvectorIG0 doesnt match, unless sort values and reset index', pd.Series(self.dye_debug['IG1']).astype(int) .equals(   self.asvectorIG0['x'].sort_values() ) )
            print('asvectorIG0 vs IG1 match:', pd.Series(self.dye_debug['IG1']).astype(int) .equals(   self.asvectorIG0['x'].sort_values().reset_index(drop=True) ))
            self.compare(list( pd.Series(self.dye_debug['IG1']).astype(int) ),
                         list( self.asvectorIG0['x'].sort_values() ),
                         atol=1)

            print('IR2 match:', self.test_IR2.astype(int).equals(pd.Series(self.dye_debug['IR2']).astype(int)))
            self.compare(self.test_IR2, self.dye_debug['IR2'])
            print('IG2 match:', self.test_IG2.astype(int).equals(pd.Series(self.dye_debug['IG2']).astype(int)))
            self.compare(self.test_IG2, self.dye_debug['IG2'])

            print('** transformed values, mean_diff **')
            print('II meth', (self.ses_II['meth'] - self.dye_debug['IIm'].sort_index()).mean())
            (self.ses_II['meth'] - self.dye_debug['IIm'].sort_index()).plot.hist(bins=100)
            plt.show()
            print('II unmeth', (self.ses_II['unmeth'] - self.dye_debug['IIu'].sort_index()).mean())
            (self.ses_II['unmeth'] - self.dye_debug['IIu'].sort_index()).plot.hist(bins=100)
            plt.show()
            print('IR meth', (self.ses_IR['meth'] - self.dye_debug['IR'].sort_index()).mean())
            (self.ses_IR['meth'] - self.dye_debug['IR'].sort_index()).plot.hist(bins=100)
            plt.show()
            print('IR unmeth', (self.ses_IR['unmeth'] - self.dye_debug['IRu'].sort_index()).mean())
            (self.ses_IR['unmeth'] - self.dye_debug['IRu'].sort_index()).plot.hist(bins=100)
            plt.show()
            print('IG meth', (self.ses_IG['meth'] - self.dye_debug['IGm'].sort_index()).mean())
            (self.ses_IG['meth'] - self.dye_debug['IGm'].sort_index()).plot.hist(bins=100)
            plt.show()
            print('IG unmeth', (self.ses_IG['unmeth'] - self.dye_debug['IG'].sort_index()).mean())
            (self.ses_IG['unmeth'] - self.dye_debug['IG'].sort_index()).plot.hist(bins=100)
            plt.show()
        except Exception as e:
            print(e)
            pass

        input_dataframe = self.containers[0]._SampleDataContainer__data_frame
        self.containers[0]._postprocess(input_dataframe, methylprep.processing.postprocess.calculate_beta_value, 'beta_value', offset=0)
        self.pre_dye = self.containers[0].betas_from_raw # presorted df
        self.post_dye = pd.DataFrame(self.containers[0]._SampleDataContainer__data_frame['beta_value'].sort_index())
        print('Effect of dye bias on beta values:', (self.post_dye - self.pre_dye).mean())
        (self.post_dye - self.pre_dye).plot.hist(bins=100)
        plt.show()

        ## post - pre split by probe types
        IG = self.containers[0].IG.index
        IR = self.containers[0].IR.index
        II = self.containers[0].II.index
        pre_IG = self.pre_dye[self.pre_dye.index.isin(IG)]
        pre_IR = self.pre_dye[self.pre_dye.index.isin(IR)]
        pre_II = self.pre_dye[self.pre_dye.index.isin(II)]
        (self.post_dye[self.post_dye.index.isin(IG)] - pre_IG).plot.hist(bins=100, title= f"post - pre IG {(self.post_dye[self.post_dye.index.isin(IG)] - pre_IG).mean()}")
        (self.post_dye[self.post_dye.index.isin(IR)] - pre_IR).plot.hist(bins=100, title= f"post - pre IR {(self.post_dye[self.post_dye.index.isin(IR)] - pre_IR).mean()}")
        (self.post_dye[self.post_dye.index.isin(II)] - pre_II).plot.hist(bins=100, title= f"post - pre II {(self.post_dye[self.post_dye.index.isin(II)] - pre_II).mean()}")
        plt.show()


        print("sesame betas (raw) vs mprep betas (raw)")
        self.ses_raw_betas = pd.read_csv(Path(LOCAL,'ses_raw_betas.csv')).rename(columns={'Unnamed: 0':'IlmnID','x':'beta_value'}).set_index('IlmnID').sort_index().round(5)
        self.ses_raw_betas = self.ses_raw_betas[~self.ses_raw_betas.index.str.startswith('rs')]
        if self.ses_raw_betas.shape == self.pre_dye.shape:
            (self.ses_raw_betas - self.pre_dye).plot.hist(bins=100)
            print("sesame betas (dye) vs mprep betas (dye)")
        else:
            print('shape doesnt match pre-dye')
        self.ses_dye_betas = pd.read_csv(Path(LOCAL,'ses_dye_betas.csv')).rename(columns={'Unnamed: 0':'IlmnID','x':'beta_value'}).set_index('IlmnID').sort_index().round(5)
        self.ses_dye_betas = self.ses_dye_betas[~self.ses_dye_betas.index.str.startswith('rs')]
        (self.ses_dye_betas - self.post_dye.sort_index()).plot.hist(bins=100)


    @staticmethod
    def compare(L1, L2, atol=1):
        try:
            L1 = L1.sort_index()
            L2 = L2.sort_index()
        except:
            pass
        print( round(100*(np.isclose(L1, L2, atol=atol, equal_nan=True)).sum()/ len(L1),4) )
