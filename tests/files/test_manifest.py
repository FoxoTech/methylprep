# App
from methylprep.files import manifests
from methylprep.models import ArrayType
from pathlib import Path


class TestManifestConstants():
    def test_has_correct_path_values(self):
        assert manifests.MANIFEST_DIR_NAME == '.methylprep_manifest_files'
        assert manifests.MANIFEST_DIR_PATH == '~/.methylprep_manifest_files'
        assert manifests.MANIFEST_REMOTE_PATH == 'https://s3.amazonaws.com/array-manifest-files/'

    def test_array_type_filenames(self):
        array_type_filenames = manifests.ARRAY_TYPE_MANIFEST_FILENAMES

        assert array_type_filenames[ArrayType.ILLUMINA_450K] == 'HumanMethylation450_15017482_v1-2.CoreColumns.csv.gz'
        assert array_type_filenames[ArrayType.ILLUMINA_EPIC] == 'MethylationEPIC_v-1-0_B4.CoreColumns.csv.gz'
        assert array_type_filenames[ArrayType.ILLUMINA_EPIC_PLUS] == 'CombinedManifestEPIC.manifest.CoreColumns.csv.gz'

    def _control_probe_selection(self):
        """ checks that control sections of each manifest match the number of probes the model expects.
        this TEST DOES NOT RUN because it would require downloading all 4 manifest files each time. """
        array_type_filenames = manifests.ARRAY_TYPE_MANIFEST_FILENAMES # array_type : manifest_filename

        files = {
            array_type : Path(manifests.MANIFEST_DIR_PATH, manifest_filename).expanduser()
            for array_type, manifest_filename in array_type_filenames.items()
            }
        for array_type, filepath in files.items():
            man = manifests.Manifest(array_type, filepath)
            if ArrayType(array_type).num_controls != man._Manifest__control_data_frame.shape[0]:
                raise AssertionError(f'Control probes found ({man._Manifest__control_data_frame.shape[0]}) in file ({filepath}) does not match expected number: {ArrayType(array_type).num_controls}')
