# App
from methylprep.files import manifests
from methylprep.models import ArrayType
from pathlib import Path
from methylprep.utils.files import download_file
import pytest

class TestManifestConstants():
    def test_has_correct_path_values(self):
        assert manifests.MANIFEST_DIR_NAME == '.methylprep_manifest_files'
        assert manifests.MANIFEST_DIR_PATH == '~/.methylprep_manifest_files'
        assert manifests.MANIFEST_REMOTE_PATH == 'https://s3.amazonaws.com/array-manifest-files/'

    def __test_control_probe_selection(self):
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


    def test_non_manifest_file_fails(self):
        with pytest.raises(ValueError) as e:
            TEST_FILEPATH = Path('docs/example_data/mouse/open_sesame_mouse_betas_subdata.csv')
            manifests.Manifest(ArrayType('27k'), TEST_FILEPATH)
            # Usecols do not match columns
        with pytest.raises(EOFError) as e:
            TEST_FILEPATH = Path('docs/example_data/epic_plus/samplesheet.csv')
            manifests.Manifest(ArrayType('mouse'), TEST_FILEPATH)

    def test_manifest_get_probe_details_errors(self):
        man = manifests.Manifest(manifests.ArrayType('mouse'))
        with pytest.raises(Exception) as e:
            man.get_probe_details('III', manifests.Channel('RED'))
        with pytest.raises(Exception) as e:
            man.get_probe_details(manifests.ProbeType('II'), 'RED')
        with pytest.raises(ValueError) as e:
            man.get_probe_details(manifests.ProbeType('II'), manifests.Channel('BLACK'))
        if man.get_probe_details(manifests.ProbeType('I'), manifests.Channel('Grn')).shape != (17469, 10):
            raise ValueError(f"get_probe_details (used in infer channel) shape mismatch: IG {man.get_probe_details(manifests.ProbeType('I'), manifests.Channel('Grn')).shape}")
        if man.get_probe_details(manifests.ProbeType('I'), manifests.Channel('Red')).shape != (46545, 10):
            raise ValueError(f"get_probe_details (used in infer channel) shape mismatch: IR {man.get_probe_details(manifests.ProbeType('I'), manifests.Channel('Red')).shape}")
        if man.get_probe_details(manifests.ProbeType('II'), None).shape != (227699, 10):
            raise ValueError(f"get_probe_details (used in infer channel) shape mismatch: II {man.get_probe_details(manifests.ProbeType('II'), None).shape}")
        if man.get_probe_details(manifests.ProbeType('II'), manifests.Channel('Grn')).shape != (0, 10):
            raise ValueError(f"get_probe_details (used in infer channel) shape mismatch: II-G {man.get_probe_details(manifests.ProbeType('II'), manifests.Channel('Grn')).shape}")
        if man.get_probe_details(manifests.ProbeType('II'), manifests.Channel('Red')).shape != (2, 10):
            raise ValueError(f"get_probe_details (used in infer channel) shape mismatch: II-R {man.get_probe_details(manifests.ProbeType('II'), manifests.Channel('Grn')).shape}")
