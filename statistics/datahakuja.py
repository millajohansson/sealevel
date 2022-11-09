# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 12:43:32 2022

@author: johanssm
"""

# Scripts for getting sea level data from databases

import requests
import numpy as np
import sys

sys.path.append('../tools/') 

from mareos import fmisid_to_chr, mareo_names_ver2, N2000mar
from mw_tools import read_mw2

# Get annual and monthly mean/max/min statistics from Smartmet Server,
# convert to N2000, and save both MW and N2000 statistics to text files
def get_smartmet_statistics(yr1 = 1887, yr2 = 2021, outdir = 'd:/tulos/'):

    fmisids = [134254,100669,132310,134253,134252,134225,134224,134266,134251,134223,134250,100540,134248,100539]
    
    url = 'http://smartmet.fmi.fi/timeseries'
    payload = {
        "producer": "observations_fmi",
         "starttime": str(yr1) + "-01-01T00:00:00",
         "endtime": str(yr2 + 1) + "-01-01T00:00:00",
         "timestep": "data",
         "format": "json"
    }
    
    mw = read_mw2()

    # Get monthly and annual statistics separately    
    for mon in (True, False):

        payload["param"] = "fmisid,utctime," + \
            ("WLEV_P1M_AVG,WLEV_P1M_MAX,WLEV_P1M_MIN" if mon else \
             "WLEV_P1Y_AVG,WLEV_P1Y_MAX,WLEV_P1Y_MIN")

        # for all tide gauges        
        for sid in fmisids:
            mar = fmisid_to_chr[sid]

            payload['fmisid'] = sid
            r = requests.get(url, params=payload)
            rj = r.json()

            txt_out = 'Year' + ('\tMonth' if mon else '') + '\tMean\tMax\tMin\n'
            txt_out_n2 = txt_out
            for row in rj:
                if row['utctime'][(6 if mon else 4):] != (('' if mon else '01') + '01T000000'):
                    continue
                if row['fmisid'] == sid:

                    yr = int(row['utctime'][0:4])                    

                    mmean = row['WLEV_P1M_AVG' if mon else 'WLEV_P1Y_AVG']
                    mmax = row['WLEV_P1M_MAX' if mon else 'WLEV_P1Y_MAX'] 
                    mmin = row['WLEV_P1M_MIN' if mon else 'WLEV_P1Y_MIN'] 

                    if mmean == None:
                        mmean = np.nan
                    if mmax == None:
                        mmax = np.nan
                    if mmin == None:
                        mmin = np.nan

                    txt_out += row['utctime'][0:4] + \
                        (('\t' + row['utctime'][4:6]) if mon else '')
                    txt_out_n2 += row['utctime'][0:4] + \
                        (('\t' + row['utctime'][4:6]) if mon else '')
                    
                    txt_out += '\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(mmean, mmax, mmin)
                    txt_out_n2 += '\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(mmean + mw[yr][mar] - N2000mar[mar], \
                            mmax + mw[yr][mar] - N2000mar[mar], \
                            mmin + mw[yr][mar] - N2000mar[mar])
        
            with open(outdir + mareo_names_ver2[fmisid_to_chr[sid]].replace('ä','a').replace('ö','o') \
                  + ('_monthly' if mon else '_annual') + 'stats_MW.txt', 'w') as f:
                f.write(txt_out)
            with open(outdir + mareo_names_ver2[fmisid_to_chr[sid]].replace('ä','a').replace('ö','o') \
                  + ('_monthly' if mon else '_annual') + 'stats_N2000.txt', 'w') as f:
                f.write(txt_out_n2)
