# App
from methylprep.files import manifests
from methylprep.models import ArrayType
from pathlib import Path
from methylprep.utils.files import download_file

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

    def test_download_manifest_dummy_file(self):
        """ will download a tiny file from the array-manifest-files s3 bucket, to test the SSL connection on all platforms.
        The dummy file is not a proper manifest CSV, so doesn't test format.
        download_file now defaults to non-SSL if SSL fails, with warning to user."""
        test_filename = 'unittest.txt'
        test_s3_bucket = 'https://array-manifest-files.s3.amazonaws.com'  # 's3://array-manifest-files'
        dest_dir = ''
        # use the .download_file() method in files.py to test the download step specifically. this is called by Manifests() class.
        download_file(test_filename, test_s3_bucket, dest_dir, overwrite=False)
        # in testing mode, this should not exist, and should get deleted right after each successful test.
        if not Path(dest_dir,test_filename).is_file():
            raise AssertionError()
        Path(dest_dir,test_filename).unlink() # deletes file.
