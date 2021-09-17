#PATH = '/Volumes/LEGX/Barnes/44668_MURMETVEP/204617710009'
#PATH = '/Volumes/LEGX/Barnes/48230_MURMETVEP/361821/204879580038'
PATH = 'docs/example_data/mouse'
import methylprep
from pathlib import Path
import pandas as pd

def test_mouse_array_for_bugs():
    print('* loading mouse manifest')
    manifest = methylprep.files.Manifest(methylprep.models.ArrayType('mouse'))
    if not isinstance(manifest.mouse_data_frame, pd.DataFrame):
        raise AssertionError("mouse_data_frame")
    print('* loading one idat pair of files')
    #green_filepath = Path(PATH, '204617710009_R06C02_Grn.idat')
    red_filepath = Path(PATH, '204879580038_R06C02_Red.idat')
    #red_filepath = Path(PATH, '204617710009_R06C02_Red.idat')
    green_filepath = Path(PATH, '204879580038_R06C02_Grn.idat')
    green_idat = methylprep.files.IdatDataset(green_filepath, channel=methylprep.models.Channel.GREEN)
    red_idat = methylprep.files.IdatDataset(red_filepath, channel=methylprep.models.Channel.RED)
    print(f"* GREEN --> {green_filepath.name} -- {green_idat.probe_means.shape}")
    print(f"* RED --> {red_filepath.name} -- {red_idat.probe_means.shape}")
    sample = 1
    sigset = methylprep.models.SigSet(sample, green_idat, red_idat, manifest, debug=True)
    """ on 2021-08-23, got this:
        II 228271
        IG 17697
        IR 47231
        oobG 47231
        oobR 17697
        methylated 293196
        unmethylated 293197
        snp_methylated 1485
        snp_unmethylated 1486
        ibG 245968
        ibR 275502
    """
    print("Checking sigset for any duplicated probes...")
    for subset, parts in sigset.subsets.items():
        print(f"----- {subset}: {getattr(sigset,subset).shape[0]}, {getattr(sigset,subset).index.duplicated().sum()}")

    print("manifest starts:")
    print(f"{sigset.man.shape} {sigset.ctl_man.shape} {sigset.snp_man.shape}") # later: man: (291713, 10) + 1486 = 293199
    df = methylprep.run_pipeline(PATH, debug=True)
    #print('* raw_dataset')
    #raw_dataset = methylprep.processing.raw_dataset.RawDataset(sample, green_idat, red_idat)
    #print('* meth_dataset.unmethylated')
    #unmethylated = methylprep.models.MethylationDataset.unmethylated(raw_dataset, manifest)
    #print('* meth_dataset.methylated')
    #methylated = methylprep.models.MethylationDataset.methylated(raw_dataset, manifest)
    #return green_idat, red_idat
    #grn, red = test_noob_df_same_size()
    #manifest = methylprep.files.Manifest(methylprep.models.ArrayType('450k'))
    #green_filepath = Path(PATH, '9247377085_R04C02_Grn.idat')
    #red_filepath = Path(PATH, '9247377085_R04C02_Red.idat')
    #print(f"* GREEN --> {green_filepath.name}")
    #print(f"* RED --> {red_filepath.name}")
    #green_idat = methylprep.files.IdatDataset(green_filepath, channel=methylprep.models.Channel.GREEN)
    #red_idat = methylprep.files.IdatDataset(red_filepath, channel=methylprep.models.Channel.RED)
    #sample = 1
    #sigset = methylprep.models.SigSet(sample, green_idat, red_idat, manifest, debug=True)
