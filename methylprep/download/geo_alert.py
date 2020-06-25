# Lib
import logging
from pathlib import Path, PurePath
from urllib.request import urlopen
from bs4 import BeautifulSoup
import dateutil # python-dateutil, non-built-in
import datetime
import pandas as pd
from tqdm import tqdm
# App
from .process_data import confirm_dataset_contains_idats, get_attachment_info

def search(keyword):
    """ CLI/cron function to check for new datasets.
    set up as a weekly cron.
    uses a local storage file to compare with old datasets in <pattern>_meta.csv.
    saves the dates of each dataset from GEO; calculates any new ones as new rows. updates csv.

options:
    pass in -k keyword
    verbose (True|False) --- reports to page; saves csv too

returns:
    saves a CSV to disk and returns a dataframe of results"""
    # FIRST: get any previous search results
    filename = f'geo_alert {keyword}.csv'
    if Path(filename).exists():
        prev_data = pd.read_csv(Path(filename))
        print(f'Previous search: {len(prev_data)} results')
    else:
        prev_data = None

    ROOT = 'http://www.ncbi.nlm.nih.gov/'
    # NOTE: query_page was limited to 20 reesults, so had to use summary table page instead
    query_page = 'http://www.ncbi.nlm.nih.gov/gds/?term=GPL13534'
    summary_page = 'http://www.ncbi.nlm.nih.gov/geo/browse/?view=series&display=500&zsort=date&search='  #methylation'
    if keyword:
        for word in keyword.split():
            query_page += f'+{word}'
            summary_page += f'+{word}'
    query_page += '+AND+%22gse%22%5BEntry+Type%5D' # limits to datasets
    print(summary_page)
    summary_html = urlopen(summary_page).read()
    #query_html = urlopen(query_page).read()
    try:
        soup = BeautifulSoup(summary_html, 'html.parser')
        total_series = soup.find(id='count')
        print(total_series.text)
        table = soup.find(id="geo_data")
        data = []
        for row in tqdm(table.find_all('tr'), desc='Checking for idats'):
            gse = row.find('a') # always first row
            gse = gse.text if gse else 'None'
            if gse == 'None':
                continue # header
            title = row.find('td', {'class': 'title'})
            title = title.text.strip() if title else 'No title'
            samples = [i.text for i in row.find_all('a') if 'samples' in i.get('href','')]
            samples = samples[0] if len(samples) > 0 else 0
            try:
                samples = int(samples)
            except:
                if samples == 'Samples':
                    continue
            date = row.find('td', {'class': 'date'})
            if date:
                date = dateutil.parser.parse(date.text)
                date = date.strftime('%Y-%m-%d') # YEAR-month-day
            else:
                date = ''
            ROW = {
                'title': title,
                'GSE': gse,
                'url': ROOT + row.find('a').get('href', ''),
                'samples': samples,
                'date': date,
                'platform': None,
                'idats': None,
            }
            geo_acc_page = f"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
            html = urlopen(geo_acc_page).read()
            PLATFORMS = ['GPL21145', 'GPL13534', 'GPL23976', 'GPL8490', 'GPL16304']
            for platform in PLATFORMS:
                if platform in str(html):
                    ROW['platform'] = platform
                    break
            usable = confirm_dataset_contains_idats(gse) # True|False
            info = get_attachment_info(gse) # list of dicts for each file with name, link, size
            # prefill, so csv is square; always 3 file spots saved
            for i in range(3):
                ROW[f'file_name_{i}'] = ''
                ROW[f'file_size_{i}'] = ''
                ROW[f'file_link_{i}'] = ''
            for i, filedata in enumerate(info):
                ROW[f'file_name_{i}'] = filedata['name']
                ROW[f'file_size_{i}'] = filedata['size']
                ROW[f'file_link_{i}'] = filedata['link']
                if i == 2:
                    break
            ROW['idats'] = '1' if usable else '0'
            #soup = BeautifulSoup(query_html, 'html.parser')
            data.append(ROW)
    except Exception as e:
        print(e)
    if type(prev_data) != type(None):
        new_results = len(data) - len(prev_data)
        print(f'{new_results} new results found.')
        # do stuff here with only the new results.
        # need a hook API here (i.e. slackbot) - call another CLI function
    # overwrite results either way.
    df = pd.DataFrame(data)
    df.to_csv(filename)
    print(filename,'written')
    return df
