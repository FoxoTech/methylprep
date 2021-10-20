# -*- coding: utf-8 -*-
from pathlib import Path
import gzip
TESTPATH = 'docs/example_data/read_geo'
#app
import methylprep
import methylcheck

class TestReadGeo():
    # a more complete version of this test is in methylcheck, where the function is maintained then copied here.
    # these files are small sub-sets of the real GEO files, created to ensure read_geo() can parse various file structures without commiting huge files to the repo.
    unit_test_files = [
        'GSE111165_test.csv',
        'GSE72120_test.txt',
        'GSE133355_processed_test.xlsx', # multiline header here, and extra columns that are not sample betas.
    ]
    unit_test_file_shapes = {
        'GSE111165_test.csv': (200, 101),
        'GSE72120_test.txt': (200, 72),
        'GSE133355_processed_test.xlsx': (6735, 44),
    }

    def test_read_csv(self):
        files = [file for file in self.unit_test_files if ('.csv' in Path(TESTPATH,file).suffixes)]
        for infile in files:
            df = methylcheck.read_geo(Path(TESTPATH,infile), verbose=False)
            if not hasattr(df,'shape'):
                raise AssertionError(f"[CSV] {infile.name} failed to return a dataframe")
            if df.shape != self.unit_test_file_shapes[Path(infile).name]:
                raise ValueError(f"[CSV] {infile.name} shape did not match ({df.shape} vs {unit_test_file_shapes[Path(infile).name]})")

    def test_read_xlsx(self):
        files = [file for file in self.unit_test_files if Path(TESTPATH,file).suffix == '.xlsx']
        for infile in files:
            df = methylcheck.read_geo(Path(TESTPATH,infile), verbose=False)
            if not hasattr(df,'shape'):
                raise AssertionError(f"[XLSX] {infile} failed to return a dataframe")

    def test_read_txt(self):
        files = [file for file in self.unit_test_files if Path(TESTPATH,file).suffix == '.txt']
        for infile in files:
            df = methylcheck.read_geo(Path(TESTPATH,infile), verbose=False)
            if not hasattr(df,'shape'):
                raise AssertionError(f"[TXT] {infile} failed to return a dataframe")
