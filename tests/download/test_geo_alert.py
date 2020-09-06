import methylprep
import pandas as pd
from pathlib import Path

def test_geo_alert(keys="blood spleen"):
    df = methylprep.download.search(keys)
    if type(df) is not type(pd.DataFrame()):
        raise AssertionError()
    if Path(f'geo_alert {keys}.csv').exists():
        Path(f'geo_alert {keys}.csv').unlink()
    else:
        print(f'file not found')
