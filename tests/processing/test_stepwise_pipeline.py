#PATH = '/Volumes/LEGX/Barnes/44668_MURMETVEP/204617710009'
#PATH = '/Volumes/LEGX/Barnes/48230_MURMETVEP/361821/204879580038'
PATH = '/docs/example_data/mouse_noob'
import methylprep
from pathlib import Path

def _test_noob_df_same_size():
    print('* loading mouse manifest')
    manifest = methylprep.files.Manifest(methylprep.models.ArrayType('mouse'))
    print('* loading one idat pair of files')
    #green_filepath = Path(PATH, '204617710009_R06C02_Grn.idat')
    green_filepath = Path(PATH, '204879580038_R06C02_Red.idat')
    #red_filepath = Path(PATH, '204617710009_R06C02_Red.idat')
    red_filepath = Path(PATH, '204879580038_R06C02_Grn.idat')
    print(f"* GREEN --> {green_filepath.name}")
    print(f"* RED --> {red_filepath.name}")
    green_idat = methylprep.files.IdatDataset(green_filepath, channel=methylprep.models.Channel.GREEN)
    red_idat = methylprep.files.IdatDataset(red_filepath, channel=methylprep.models.Channel.RED)
    sample = 1
    print('* raw_dataset')
    raw_dataset = methylprep.processing.raw_dataset.RawDataset(sample, green_idat, red_idat)

    print('* meth_dataset.unmethylated')
    unmethylated = methylprep.models.MethylationDataset.unmethylated(raw_dataset, manifest)

    print('* meth_dataset.methylated')
    methylated = methylprep.models.MethylationDataset.methylated(raw_dataset, manifest)
    return green_idat, red_idat

#grn, red = test_noob_df_same_size()
