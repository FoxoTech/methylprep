import pytest
import methylprep
from pathlib import Path
PATH = 'docs/example_data/58496/samplesheet.csv'

def test_sample_sheet_meta_data():
    """ TEST: meta_data and samplesheet objects match, when using alt column names in a samplesheet """
    sample_sheet = methylprep.files.SampleSheet(PATH, PATH)
    #sample_sheet.fields SHOULD BE a complete mapping of original and renamed_fields
    cols = list(sample_sheet.fields.values()) + ['Sample_ID']
    with pytest.raises(ValueError):
        sample_sheet.build_meta_data()
    samples = sample_sheet.get_samples()
    meta_data = sample_sheet.build_meta_data(samples)
    assert sample_sheet.fields == {'Unnamed: 0': 'Unnamed_0', 'Sentrix_ID': 'Sentrix_ID',
    'Sentrix_Position': 'Sentrix_Position', 'Sample_Group': 'Sample_Group',
    'Sample_Name': 'Sample_Name', 'Sample_Plate': 'Sample_Plate',
    'Sample_Type': 'Sample_Type', 'Sub_Type': 'Sub_Type', 'Sample_Well': 'Sample_Well', 'Pool_ID': 'Pool_ID',
    'GSM_ID': 'GSM_ID', 'Control': 'Control'}
