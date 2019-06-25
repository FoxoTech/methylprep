# App
from methpype.files import manifests
from methpype.models import ArrayType


class TestManifestConstants():
    def test_has_correct_path_values(self):
        assert manifests.MANIFEST_DIR_NAME == '.methpype_manifest_files'
        assert manifests.MANIFEST_DIR_PATH == '~/.methpype_manifest_files'
        assert manifests.MANIFEST_REMOTE_PATH == 'https://s3.amazonaws.com/array-manifest-files/'

    def test_array_type_filenames(self):
        array_type_filenames = manifests.ARRAY_TYPE_MANIFEST_FILENAMES

        assert array_type_filenames[ArrayType.ILLUMINA_450K] == 'HumanMethylation450_15017482_v1-2.CoreColumns.csv.gz'
        assert array_type_filenames[ArrayType.ILLUMINA_EPIC] == 'MethylationEPIC_v-1-0_B4.CoreColumns.csv.gz'
        assert array_type_filenames[ArrayType.ILLUMINA_EPIC_PLUS] == 'CombinedManifestEPIC.manifest.CoreColumns.csv.gz'
