# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:41:26 2020

@author: johanssm
"""

# Scripts for calculating sea level statistics for the 30-year reference 
# period 1991-2020.
#
# 03/2021 Milla Johansson

# The WMO criteria considering calculated/corrected/uncertain values
# have been implemented in the code, but their operation not completely 
# checked as the WMO criteria were not applied in the final calculations.
# Note: WMO criteria considering missing values have not been implemented
# at all, as there were no missing values in 1991-2020 sea level data.


import numpy as np
from time import strptime
from calendar import timegm
import matplotlib.pyplot as plt
from matplotlib import rc

####### Define some constants related to tide gauge data: ###############

# Tide gauge names
mareo_names = {'a': 'Kemi', 'o': 'Oulu', 'b': 'Raahe', 'p': 'Pietarsaari', \
           'v': 'Vaasa', 's': 'Kaskinen', 'm': 'Mantyluoto', 'r': 'Rauma', \
           't': 'Turku', 'd': 'Foglo', 'h': 'Hanko', 'e': 'Helsinki', \
           'f': 'Hamina', 'c': 'Porvoo'}

mareo_names_long = {'a': 'Kemi Ajos', 'o': 'Oulu Toppila', 'b': 'Raahe Lapaluoto', 
                    'p': 'Pietarsaari Leppäluoto', \
                    'v': 'Vaasa Vaskiluoto', 's': 'Kaskinen Ådskär', 
                    'm': 'Pori Mäntyluoto Kallo', 'r': 'Rauma Ulko-Petäjäs', \
                    't': 'Turku Ruissalo Saaronniemi', 'd': 'Föglö Degerby', 
                    'h': 'Hanko Pikku Kolalahti', 'e': 'Helsinki Kaivopuisto', \
                    'f': 'Hamina Pitäjänsaari', 'c': 'Porvoo Emäsalo Vaarlahti'}

mareo_names_tables = {'a': 'Kemi',
                 'o': 'Oulu - Uleåborg',
                 'b': 'Raahe - Brahestad',
                 'p': 'Pietarsaari - Jakobstad',
                 'v': 'Vaasa - Vasa',
                 's': 'Kaskinen - Kaskö',
                 'm': 'Mäntyluoto',
                 'r': 'Rauma - Raumo',
                 't': 'Turku - Åbo',
                 'd': 'Föglö',
                 'h': 'Hanko - Hangö',
                 'e': 'Helsinki - Helsingfors',
                 'f': 'Hamina - Fredrikshamn'}

    
# Station id's
fmisids = {'f': 134254, 'c': 100669, 'e': 132310, 'h': 134253, \
       'd': 134252, 't': 134225, 'r': 134224, 'm': 134266, \
       's': 134251, 'v': 134223, 'p': 134250, 'b': 100540, \
       'o': 134248, 'a': 100539}

# N2000 vs. reference level at tide gauges
N2000mar = {'a': 1264, 'o': 1288, 'b': 1242, 'p': 1202, 'v': 1233, 's': 1270, \
            'm': 1331, 'r': 1386, 't': 1512, 'd': 1522, 'h': 1609, 'e': 1642, \
            'c': 1651, 'f': 1674}


############ File readers for different data file types: ###############

def read_from_file(fname):    
# Read an ascii file linewise
    try:
        with open(fname, 'r') as f:
            data0 = f.readlines()
    except FileNotFoundError:
            print('File not found: ' + fname)
            data0 = ''
    return data0
    
    
def read_hourlydata(fname):
# Read hourly sea level data produced by 'sealevel'-program.
# Returns: {'Tstamp': timestamps, 'Sl': sea level, 'Int': interpolation flag}
# Interpolation flag = True for interpolated/calculated/adjusted values 
# (which are marked with a comma in the original data file.)

    def convert_timestamp(yr, mon, day, hr):
    
        # Calculate the beginning of the day in seconds (UTC).
        tstamp = timegm(strptime(yr + '-' + '{:02d}'.format(int(mon)) \
                    + '-' + '{:02d}'.format(int(day)), '%Y-%m-%d'))
        #  Add hours. (This also handles the timestamps with hour 24.)
        tstamp += int(hr) * 3600.0
        # Convert EET to UTC.
        tstamp -= 7200.0
        
        return tstamp
    
    
    # Read from file
    data0 = read_from_file(fname)

    cols = ('Tstamp', 'Sl', 'Int')
    
    data = {}
    for col in cols:
        data[col] = []

    for row in data0:
        # Ignore empty rows and rows not starting with a year.
        if len(row) < 5:
            continue
        if row[0] not in ('1', '2'):
            continue
        
        rspl = row.split()
        # Only process rows with 5 columns.
        if len(rspl) != 5:
            continue

        # The time stamps are given in four columns, and in EET = UTC + 2.
        # Convert to UTC seconds.
        data['Tstamp'].append(convert_timestamp(rspl[0], rspl[1], rspl[2],
                                                rspl[3]))
        # Interpolated/calculated value: set 'Int' to True:
        data['Int'].append(1 if rspl[4][-1] == ',' else 0)
        # Actual sea level value (in mm, N2000):
        data['Sl'].append(float(rspl[4].strip(',')))
                        
    for d in data:
        data[d] = np.array(data[d])
            
    return(data)            


def read_cldbdata(fname, fmisid):
# Read sea level data, csv files retrieved with SQL from CLDB.
# This reads hourly, monthly or annual data.
# The measurand id tells what the data are. E.g. 74 = hourly, 
# (313, 314, 315) = monthly, (316, 317, 318) = annual.
#
# File content: Timestamp, Station_id, Measurand_id, sea level MW, 
#   sea level N2000, Source flag, Quality flag.
#
# Returns: {'Tstamp': timestamps, 'SlMW': sea level MW, 
#        'SlN2000': sea level N2000, 'Src': source flag,
#       'Qual': quality flag}

    
    # Read from file
    data0 = read_from_file(fname)

    cols = ('Tstamp', 'SlMW', 'SlN2000', 'Src', 'Qual', 'Measid')
    
    data = {}
    for col in cols:
        data[col] = []

    for row in data0:
        # Ignore empty rows and rows not starting with a year.
        if len(row) < 5:
            continue
        if row[0] not in ('1', '2'):
            continue
        
        rspl = row.split(',')
        # Only process rows with 7 columns.
        if len(rspl) != 7:
            continue

        # Check station id. 
        if float(rspl[1]) != fmisid:
            print('Wrong fmisid in file: ' + rspl[1] + ' vs. ' + str(fmisid) + '. Skipping line.')
            continue

        # Convert timestamp to UTC seconds.
        data['Tstamp'].append(timegm(strptime(rspl[0], '%Y-%m-%dT%H:%M:%S')))

        # Actual sea level value (in mm, N2000):
        data['SlMW'].append(float(rspl[3]))

        # Actual sea level value (in mm, MW):
        data['SlN2000'].append(float(rspl[4]))

        # Measurand id:
        data['Measid'].append(float(rspl[2])) 
           
        # Source and quality flag:
        data['Src'].append(int(rspl[5]))
        try:
            data['Qual'].append(int(rspl[6]))
        except ValueError:
            data['Qual'].append(-999)
                        
    for d in data:
        data[d] = np.array(data[d])
            
    return(data)            


def read_txtdata(fname, aari = False):
# Read annual/monthly sea level statistics from "Tekstikanta" file.
# Returns: {'Yr': year, 'Mon': month, 'Mean': mean,
#       'Max': maximum, 'Min': minimum}
# aari = False => Read annual/monthly means,
# aari = True => Read annual/monthly extremes.

    
    # Read from file
    data0 = read_from_file(fname)

    cols = ('Yr', 'Mon', 'Mean', 'Max', 'Min')


    tt = str.maketrans(',*', '  ')    # Muunnetaan pilkut ja tahdet välilyönneiksi
        
    data = {}
    for col in cols:
        data[col] = []

    i = 0
    step = 2 if aari else 1
    
    while i < len(data0):

        # Ignore empty rows and rows not starting with a year.
        if len(data0[i]) < 5:
            i += 1
            continue
        if data0[i][0] not in ('1', '2'):
            i += 1
            continue

        if aari:
            row = data0[i].strip() + ' ' + data0[i + 1].strip()
        else:
            row = data0[i].strip()

        row1 = row.translate(tt).strip()

        rspl = row1.split()

        if not aari:
            if len(rspl) != 15:
                i += 1
                continue
        else:            
            if len(rspl) != 28:
                i += 1
                continue

        for mon in np.arange(1, 14, 1):
            data['Yr'].append(float(rspl[0]))
            data['Mon'].append(mon)
            if aari:
                data['Max'].append(float(rspl[mon]))
                data['Min'].append(float(rspl[mon + 14]))
            else:
                data['Mean'].append(float(rspl[mon]))
                
        i += step

                                    
    for d in data:
        data[d] = np.array(data[d])
            
    return(data)            



def read_statfile(fname, y1 = 1991, y2 = 2020, ismonthly = False):
# Read annual/monthly statistics file produced by the function process_stats().
# File content: Year, (Month), Mean, Max, Min.
    
    # Read from file
    data0 = read_from_file(fname)

    cols = ('Yr', 'Mon', 'Mean', 'Max', 'Min')
    
    data = {}
    for col in cols:
        data[col] = []

    for row in data0:
        # Ignore empty rows and rows not starting with a year.
        if len(row) < 5:
            continue
        if row[0] not in ('1', '2'):
            continue
        
        rspl = row.split()
        # Only process rows with 4/5 columns.
        if len(rspl) != (5 if ismonthly else 4):
            continue

        yr = float(rspl[0])
        if (yr < y1) or (yr > y2):
            continue
        data['Yr'].append(yr)
        if ismonthly:
            data['Mon'].append(float(rspl[1]))
        data['Mean'].append(float(rspl[-3]))
        data['Max'].append(float(rspl[-2]))
        data['Min'].append(float(rspl[-1]))
                                
    for d in data:
        data[d] = np.array(data[d])
            
    return(data)            


def read_resultfile(fname):
# Read final statistics file produced by the function annual_monthly_stats().
# File content: Month	Lowest mean	Mean	Highest mean	Std(mean)	N(mean)	
#    Lowest max	Average max	Highest max	Std(max)	N(max)	
#   Lowest min	Average min	Highest min	Std(min)	N(min)
#   Year(max)     Year(min)

    # Read from file
    data0 = read_from_file(fname)

    cols = data0[0].strip().split('\t')
    
    data = {}
    for col in cols:
        data[col] = []

    for row in data0[1:]:
        # Ignore empty rows.
        if len(row) < 5:
            continue
        
        rspl = row.split()
        # Only process rows with right number of columns.
        if len(rspl) != len(cols):
            continue

        for i in range(len(cols)):
            data[cols[i]].append(float(rspl[i]))
                                
    for d in data:
        data[d] = np.array(data[d])
            
    return(data)            


def read_distribution(fname, monthly = True):
# Read probability distribution file produced with Sealevel-program.
# Reads either monthly distribution for one station (monthly = True), 
# or the annual distribution for all stations.

    # Read from file
    data0 = read_from_file(fname)

    levels = {}

    for row in data0:
        # Ignore empty rows.
        if len(row) < 5:
            continue
        
        rspl = row.strip().split()

        # Only process rows with right number of columns.
        if len(rspl) != (13 if monthly else 15):
            continue

        if rspl[0] in ('MAX', 'MIN', 'MHW', 'MLW'):
            levels[rspl[0]] = [float(x) for x in rspl[1:]]
        elif rspl[0] in ('Nro', 'Nimi', 'Alk.'):
            continue
        else:
            levels[float(rspl[0].replace(',', '.'))] = [float(x) for x in rspl[1:]]
                                            
    return(levels) 


def read_distribution2(fname):
# Read probability distribution file produced by pick_distribution_levels()

    # Read from file
    data0 = read_from_file(fname)
    
    
    data = {}
    
    for row in data0[1:]:
        # Ignore empty rows.
        if len(row) < 5:
            continue
        
        rspl = row.split(';')
        # Only process rows with right number of columns.
        if len(rspl) != 15:
            continue

        try:
            pros = float(rspl[1])
        except ValueError:
            pros = rspl[1]

        data[pros] = [float(x) for x in rspl[2:]]

                                        
    return(data)            



def read_monthlycrit(fname):
# Read the file produced by monthly_criteria(), which lists those months
# that don't fulfill the WMO criteria.    
    
    res = {}
    for mar in 'aobpvsmrtdhef':
        res[mar] = {'yr': [], 'mon': []}
    
    txt = read_from_file(fname)
        
    for row in txt:
        if len(row) < 5:
            continue
        if row[0:6] == 'Months':
            continue
        rspl = row.split()
        for mar in 'aobpvsmrtdhef':
            if mareo_names[mar] == rspl[0]:
                break
        res[mar]['yr'].append(float(rspl[1]))
        res[mar]['mon'].append(float(rspl[2]))

    return res


def read_mw(fname):
# Read theoretical mean sea level (MW) values from file.
# The file should be 'mwtaulu_1887-2021.txt' produced by Sealevel-program.

    mareos = 'aobpvsmrtdhefc'
    mw = {'Year': []}
    for mar in mareos:
        mw[mar] = []

    mwdata = read_from_file(fname)
    
    for row in mwdata:
        if len(row) < 5:
            continue
        if row[0] not in ('1', '2'):
            continue
        
        rspl =row.split()
        if len(rspl) != 15:
            print('read_mw: faulty datarow in ' + fname + ', skipping.')
            continue
        
        mw['Year'].append(float(rspl[0]))
        for i in range(14):
            mw[mareos[i]].append(float(rspl[i + 1]))

    for k in mw:
        mw[k] = np.array(mw[k])

    return mw



############# Functions for calculating statistics etc. ##################


def qc_statistics(do_monthly = False, use_txt = False, startyear = 1991,
                  ddir = 'd:/2021/vertailukausi30/Data/Raw/',
                  ddir2 = 'd:/2021/vertailukausi30/Data_other/Raw_txt/',
                  outdir = 'd:/tulos/'):
# Count the amount of missing, interpolated, corrected etc. values
# for each month from the hourly sea level data.
# If do_monthly = True, calculate also monthly statistics (not just annual).
# The monthly calculation is very slow (> 40 min)!
# If use_txt = True, also calculate statistics from "Tekstikanta" files 
# (otherwise, only CLDB data is used).

    startmar = {'a': 1922, 'o': 1922, 'b': 1922, 'p': 1922, 'v': 1922, 's': 1926,
                'm': 1925, 'r': 1933, 't': 1922, 'd': 1923, 'h': 1887, 'e': 1904, 'f': 1928}

    if startyear != 0:
        t0 = timegm(strptime(str(startyear) + '-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
    t1 = timegm(strptime('2021-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
    if use_txt:
        # Tiedostoissa viimeinen arvo klo 24:00 Suomen aikaa, UTC 23 puuttuu.
        t1 = timegm(strptime('2020-12-31T23:00:00', '%Y-%m-%dT%H:%M:%S'))
        

    qflags = (-999, 1, 2, 3, 4)
    sflags = (1, 2, 3, 4, 5, 6, 7)
    iflags = (0, 1)

    stats = {}
    
    staul = 'Station\tTotal\tGood\tLower_q\tUnknown\tMissing\n'

    for mar in 'aobpvsmrtdhef':
        
        if startyear == 0:
            t0 = timegm(strptime(str(startmar[mar]) + '-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
                
        data = read_cldbdata(ddir + 'sealevels' + \
                                    str(fmisids[mar]) + '.txt',
                                    fmisids[mar])
        ok = (data['Tstamp'] >= t0) & (data['Tstamp'] < t1)
        for d in data:
            data[d] = data[d][ok]
        ntotal = len(data['SlMW'])

        if use_txt:
#            data2 = read_hourlydata(ddir2 + mareo_names[mar] + str(startyear) + '_2020.txt')
            data2 = read_hourlydata(ddir2 + mareo_names[mar] + '_all_n2000.txt')
            ok2 = (data2['Tstamp'] >= t0) & (data2['Tstamp'] < t1)
            for d in data2:
                data2[d] = data2[d][ok2]
            ntotal2 = len(data2['Sl'])

            if ntotal != ntotal2:
                print('Mareo ' + mar + ', different ntotals (' + str(ntotal) + \
                      ' vs ' + str(ntotal2) + '). Skipping.')
                continue
            
            print(data2['Tstamp'][0] - data['Tstamp'][0])
            print(data2['Tstamp'][190000] - data['Tstamp'][190000])
        ismissing = np.isnan(data['SlMW'])
        
        print(mareo_names[mar])

        stats[mar] = {}
        nmiss = np.sum(ismissing)
#        nsum = len(data['SlMW'])
        nsum = nmiss

        if use_txt:
            tul = 'Src\tQual\tInt\tCount\n'
            
            for s in sflags:
                for q in qflags:
                    for i in iflags:
                        statsm = np.sum((data['Qual'] == q) & \
                                                  (data['Src'] == s) & \
                                                  (data2['Int'] == i) & \
                                                  (~ismissing))
                        nsum += statsm
                        if (statsm > 0):
                            tul += '{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(s, q, i, statsm)

            tul += '-9\t-9\t-9\t{:.0f}\n'.format(nmiss)
            tul += 'All\tAll\tAll\t{:.0f}\n'.format(nsum)

            print(tul)   
            with open(outdir + mareo_names[mar] + '_statistics.txt', 'w') as f:
                f.write(tul)


        stats[mar]['good'] = np.sum((data['Qual'] == 1) | (data['Qual'] == 2))
        stats[mar]['lowq'] = np.sum((data['Qual'] == 3) | (data['Qual'] == 4))
        stats[mar]['unknown'] = np.sum(data['Qual'] == -999)
        stats[mar]['missing'] = nmiss
        stats[mar]['total'] = nsum
        
        staul += mareo_names[mar] + \
            '\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(stats[mar]['total'], \
                stats[mar]['good'], stats[mar]['lowq'], 
                stats[mar]['unknown'], stats[mar]['missing'])

            
        if do_monthly:
            mon_tul = 'Year\tMonth\tTotal\tGood\tLower_q\tUnknown\tMissing\n'
            for yr in np.arange(startyear, 2021):
                maxmon = (11 if yr == 2020 else 12)
                maxdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
                if yr in (1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020):
                    maxdays[1] = 29
                for mon in np.arange(1, maxmon + 1):
                    t0m = timegm(strptime('{:d}-{:02d}-01T00:00:00'.format(yr, mon), '%Y-%m-%dT%H:%M:%S'))
                    t1m = timegm(strptime('{:d}-{:02d}-{:02d}T23:59:59'.format(yr, mon, maxdays[mon-1]), '%Y-%m-%dT%H:%M:%S'))
                    ok = (data['Tstamp'] >= t0m) & (data['Tstamp'] < t1m)
                    mstats = {
                        'good': np.sum((data['Qual'][ok] == 1) | (data['Qual'][ok] == 2)),
                        'lowq': np.sum((data['Qual'][ok] == 3) | (data['Qual'][ok] == 4)),
                        'unknown': np.sum(data['Qual'][ok] == -999),
                        'missing': 0,
                        'total': sum(ok)
                    }
                    if (mstats['lowq'] > 0) or (mstats['unknown'] > 0) or (mstats['missing'] > 0):
                        mon_tul += ('{:d}\t{:02d}\t{:.0f}\t{:.0f}\t' + \
                                    '{:.0f}\t{:.0f}\t{:.0f}\n').format(yr, mon, \
                                        mstats['total'], mstats['good'], mstats['lowq'], \
                                        mstats['unknown'], mstats['missing'])

            print(mareo_names[mar])
            print(mon_tul)
            with open(outdir + mareo_names[mar] + '_monthly_stats.txt', 'w') as f:
                f.write(mon_tul)
            
    print(staul)
    with open(outdir + 'all_statistics.txt', 'w') as f:
        f.write(staul)

    

def monthly_criteria(limit = 240, 
                     qdir = 'd:/2021/vertailukausi30/Data/Quality/',
                     outdir = 'd:/tulos/'):
# Read the missing/interpolated statistics produced by qc_statistics()
# and pick those years/months that don't fulfill the WMO criteria.
#
# Note that there are no missing values in 1991-2020 sea level data.
# Thus, the only criterion applied is: "less than 10 days of 
# interpolated/corrected/calculated values".

    result = 'Months with more than ' + str(limit) + ' hours of lower_q or unknown data:\n'
    for mar in 'aobpvsmrtdhef':
        mstat = read_from_file(qdir + mareo_names[mar] + '_monthly_stats.txt')
        
        for row in mstat[1:]:
            rdata = [float(x) for x in row.split()]
            if rdata[4] + rdata[5] > limit:
                result += mareo_names[mar] + ' ' + row
                
    with open(outdir + 'months_failing' + str(limit) + '.txt', 'w') as f:
        f.write(result)




def process_stats(from_hourly = False, use_eet = False, from_txt = False,
                  startyear = 1991, endyear = 2020, use_wmo = False,
                  ddir = 'd:/2021/Vertailukausi30/Data/Raw/',
                  ddir2 = 'd:/2021/vertailukausi30/Data_other/Raw_txt/',
                  tdir = 'd:/2021/vertailukausi30/Data_other/Tekstikanta/',
                  outdir = 'd:/tulos/'):
# Calculate or pick monthly and annual statistics (mean, max, min)
# from different data sources.
# from_hourly: calculate from hourly (CLDB) data
# use_eet: calculate using EET, not UTC (to correspond to the old data)
# from_txt: use "Tekstikanta" annual/monthly statistics
# use_wmo: discard months which don't fulfil the WMO criteria
    
    def stats_from_hourly(data, t1, t2):
        # Calculate statistics (mean, max, min) from hourly data
        ii = (data['Tstamp'] >= t1) & (data['Tstamp'] < t2)

        if any(ii):
            dmn = [np.nanmean(data[sl][ii]),]
            dmax = [np.nanmax(data[sl][ii]),]
            dmin = [np.nanmin(data[sl][ii]),]
        else:
            dmn = [float('nan'),]
            dmax = [float('nan'),]
            dmin = [float('nan'),]

        return dmn[0], dmax[0], dmin[0]


    def stats_from_cldb(data, t1, t2, ismonthly, sl):
        # Pick statistics (mean, max, min) from CLDB monthly/annual data

        if ismonthly:
            measids = (313, 314, 315)
        else:
            measids = (316, 317, 318)

        dmn = data[sl][(data['Tstamp'] == t1) & (data['Measid'] == measids[0])]
        dmax = data[sl][(data['Tstamp'] == t1) & (data['Measid'] == measids[1])]
        dmin = data[sl][(data['Tstamp'] == t1) & (data['Measid'] == measids[2])]
        if len(dmn) == 0:
            dmn = [float('nan'),]
        if len(dmax) == 0:
            dmax = [float('nan'),]
        if len(dmin) == 0:
            dmin = [float('nan'),]

        return dmn[0], dmax[0], dmin[0]


    def fills_criteria(crit, mar, yr, mon):
    # Determine whether a month meets the WMO criteria.
    # (Uses the file produced by monthly_criteria())

        for i in range(len(crit[mar]['yr'])):
            if (crit[mar]['yr'][i] == yr) and (crit[mar]['mon'][i] == mon):
                return False
        
        return True


    # Read MW height system values used to convert Tekstikanta data from
    # reference level to MW
    if from_txt:
        mw = read_mw('d:/2021/Vertailukausi30/Data/mwtaulu_1887-2021.txt')

    if from_txt and from_hourly:
        print('Calculating from hourly Tekstikanta data not implemented.')
        return
    if not from_hourly:
        use_eet = False

    htext = '_hourly' if from_hourly else ''
    ttext = '_tekstikanta' if from_txt else ''
    
    delta = 7200.0 if (use_eet and from_hourly) else 0.0

    # Read information on months fulfilling WMO criteria
    critfile = 'd:/2021/Vertailukausi30/Data/months_failing240.txt'
    crit = read_monthlycrit(critfile)

    # Calculate for all tide gauges
    for mar in 'aobpvsmrtdhef':

        if from_hourly:
            # Read hourly data
            if startyear >= 1991:
                hdata = read_cldbdata(ddir + 'sealevels' + \
                                    str(fmisids[mar]) + '.txt',
                                    fmisids[mar])
            else:
                hdata = read_hourlydata(ddir2 + mareo_names[mar] + '_all_n2000.txt')
                hdata['SlN2000'] = hdata['Sl']
                hdata['SlMW'] = float('nan')* hdata['Sl']
        else:
            if from_txt:
                # Read Tekstikanta statistics files
                adata = read_txtdata(tdir + mareo_names[mar] + '_aariarvot_ref.txt', 
                                 aari = True)
                kdata = read_txtdata(tdir + mareo_names[mar] + '_keskiarvot_ref.txt', 
                                 aari = False)
            else:
                # Read CLDB statistics files
                ydata = read_cldbdata(ddir + 'annual' + \
                                str(fmisids[mar]) + '.txt',
                                fmisids[mar])

                mdata = read_cldbdata(ddir + 'monthly' + \
                                str(fmisids[mar]) + '.txt',
                                fmisids[mar])
                
        for hsyst in ('MW', 'N2000'):
            sl = 'Sl' + hsyst
            yout = '#Yr\tMean\tMax\tMin\n'
            mout = '#Yr\tMonth\tMean\tMax\tMin\n'
            for yr in np.arange(startyear, endyear + 1, 1):
                
                if from_txt:
                    for mon in np.arange(1, (13 if use_wmo else 14), 1):
                        ai = (adata['Yr'] == yr) & (adata['Mon'] == mon)
                        ki = (kdata['Yr'] == yr) & (kdata['Mon'] == mon)
                        if any(ai):
                            amean = kdata['Mean'][ki]
                            amax = adata['Max'][ai]
                            amin = adata['Min'][ai]
                            amean = amean[0]
                            amax = amax[0]
                            amin = amin[0]
                        else:
                            amean = float('nan')
                            amax = float('nan')
                            amin = float('nan')

                        # In Tekstikanta, missing values = -999
                        if amean < 0: amean = float('nan')
                        if amax < 0: amax = float('nan')
                        if amin < 0: amin = float('nan')

                        # Convert height system from reference level
                        if hsyst == 'N2000':
                            amean -= N2000mar[mar]
                            amax -= N2000mar[mar]
                            amin -= N2000mar[mar]
                        elif hsyst == 'MW':
                            iy = np.where(mw['Year'] == yr)[0][0]
                            amean -= mw[mar][iy]
                            amax -= mw[mar][iy]
                            amin -= mw[mar][iy]

                        # Discard months not fulfilling WMO criteria
                        if use_wmo and \
                            (not fills_criteria(crit, mar, yr, mon)):
                            amean = float('nan')
                            amax = float('nan')
                            amin = float('nan')

                        if mon == 13:
                            # Annual statistics
                            yout += '{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(yr, amean, amax, amin)
                        else:
                            # Monthly statistics
                            mout += '{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(yr, mon, \
                                                            amean, amax, amin)

                else:
                    ty = timegm(strptime(str(yr) + '-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
                    ty2 = timegm(strptime(str(yr + 1) + '-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))

                    # Pick/calculate annual statistics only if WMO criteria 
                    # are not applied. (They cannot be straightforward 
                    # applied to annual statistics.)
                    if not use_wmo:
                        if from_hourly:
                            ymn, ymax, ymin = stats_from_hourly(hdata, 
                                               t1 = ty - delta, t2 = ty2 - delta)
                        else:
                            ymn, ymax, ymin = stats_from_cldb(ydata, 
                                               t1 = ty - delta, t2 = ty2 - delta, 
                                               ismonthly = False, sl = sl)

                        yout += '{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(yr, ymn, ymax, ymin)
    
                    for mon in np.arange(1, 13, 1):
                        tm = timegm(strptime('{:.0f}-{:.0f}-01T00:00:00'.format(yr, mon), 
                                             '%Y-%m-%dT%H:%M:%S'))
                        tm2 = timegm(strptime('{:.0f}-{:.0f}-01T00:00:00'.format((yr + 1 if mon == 12 else yr), (1 if mon == 12 else mon + 1)), 
                                              '%Y-%m-%dT%H:%M:%S'))
    
                        if use_wmo and (not fills_criteria(crit, mar, yr, mon)):
                            mmn = float('nan')
                            mmax = float('nan')
                            mmin = float('nan')
                        else:
                            if from_hourly:
                                mmn, mmax, mmin = stats_from_hourly(hdata, 
                                               t1 = tm - delta, t2 = tm2 - delta)
                            else:
                                mmn, mmax, mmin = stats_from_cldb(mdata, 
                                               t1 = tm - delta, t2 = tm2 - delta, 
                                               ismonthly = True, sl = sl)

                        mout += '{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\n'.format(yr, mon, \
                                                            mmn, mmax, mmin)

            # Save results    
            if not use_wmo:
                with open(outdir + mareo_names[mar] + '_annualstats' + \
                      htext + '_' + hsyst + ('_eet' if use_eet else '') + \
                          ttext + '.txt', 'w') as f:
                    f.write(yout)

            with open(outdir + mareo_names[mar] + '_monthlystats' + \
                      htext + '_' + hsyst + ('_eet' if use_eet else '') + \
                          ttext + ('_wmo' if use_wmo else '') + '.txt', 'w') as f:
                f.write(mout)



def compare_stats(data1 = 'hourly_eet', data2 = 'stat_txt',
                  sdir = 'd:/2021/vertailukausi30/Data/Stats/',
                  outdir = 'd:/tulos/', skip_mw = False):
# Compare annual/monthly statistics calculated from different source data.
# Plot time series of differences.
# Possible datasets (data1, data2):
#   hourly_utc: calculated from hourly data, UTC
#   hourly_eet: calculated from hourly data, EET
#   stat_txt: statistics from Tekstikanta
#   stat_cldb: statistics from CLDB

    def read_stats(dtype, sdir, mar, ym, hsyst):
        if dtype == 'hourly_utc':
            ystat = read_statfile(sdir + mareo_names[mar] + '_' + \
                      ('monthly' if ym else 'annual') + 'stats_hourly_' + \
                      hsyst + '.txt', ismonthly = ym)
        elif dtype == 'hourly_eet':
            ystat = read_statfile(sdir + mareo_names[mar] + '_' + \
                      ('monthly' if ym else 'annual') + 'stats_hourly_' + \
                      hsyst + '_eet.txt', ismonthly = ym)
        elif dtype == 'stat_cldb':
            ystat = read_statfile(sdir + mareo_names[mar] + '_' + \
                      ('monthly' if ym else 'annual') + 'stats_' + \
                      hsyst + '.txt', ismonthly = ym)
        elif dtype == 'stat_txt':
            ystat = read_statfile(sdir + mareo_names[mar] + '_' + \
                      ('monthly' if ym else 'annual') + 'stats_' + \
                      hsyst + '_tekstikanta.txt', ismonthly = ym)

        return ystat


    sout = ''
    for ym in (False, True):
        for mar in 'aobpvsmrtdhef':
            for hsyst in ('MW', 'N2000'):
                if skip_mw and (hsyst == 'MW'):
                    continue
                sout += mareo_names[mar] + ' ' + hsyst + (' monthly' if ym else ' annual') + '\n'

                ystat1 = read_stats(data1, sdir, mar, ym, hsyst)
                ystat2 = read_stats(data2, sdir, mar, ym, hsyst)

                fig, ax = plt.subplots(figsize = [15, 8])
                for dd in ('Mean', 'Max', 'Min'):
                    ddiff = ystat1[dd] - ystat2[dd]
                    sout += dd + ' diff.: {:.0f} - {:.0f}\n'.format(np.nanmin(ddiff), 
                                                                  np.nanmax(ddiff))
                    ax.plot(ystat1['Yr'] + (((ystat1['Mon'] - 1) / 12) if ym else 0),
                            ddiff, label = dd)
                    
                ax.legend()
                fig.savefig(outdir + mareo_names[mar] + \
                            ('_monthly' if ym else '_annual') + \
                                '_' + hsyst + '_diffs_' + \
                                data1 + '_vs_' + data2 + '.png', dpi = 300)
                plt.close(fig)

    print(sout)


    
def stats(data, mincount = 0):
# Calculate basic statistics (mean, max, min, std, count, indmax, indmin)
# from data.
# If count of non-nan values < mincount, return nan
    if np.sum(~np.isnan(data)) >= mincount:
        return np.nanmean(data), np.nanmax(data), np.nanmin(data), \
            np.nanstd(data), np.sum(~np.isnan(data)), \
            data == np.nanmax(data), data == np.nanmin(data)
            
    else:
        return [float('nan'), float('nan'), float('nan'), \
            float('nan'), np.sum(~np.isnan(data)), [], []]


def is_leap(yr):
# Determine if yr is a leap year.
    return ( ((yr % 4) == 0) and (((yr % 100) != 0) or ((yr % 400) == 0)))

                                                                       
def calculate_annual(mstat):
# Calculate annual statistics from monthly statistics
# (Not used currently.)
    mdays = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

    y1 = mstat['Yr'][0]
    y2 = mstat['Yr'][-1]
    
    ystat = {'Yr': [], 'Mon': [], 'Max': [], 'Min': [], 'Mean': []}
    
    for yr in np.arange(y1, y2):
        mdays[1] = (29 if is_leap(yr) else 28)
        iy = (mstat['Yr'] == yr)
        ymdata = mstat['Mean'][iy]
        if len(ymdata) != 12:
            print('Virhe kuukausidatassa.')
            return {}
        
        ystat['Yr'].append(yr)
        ymn = 0
        for m in range(12):
            ymn += (mdays[m] * ymdata[m]) / np.sum(mdays)            

        ystat['Mean'].append(np.nanmean(mstat['Mean'][iy]))
        ystat['Max'].append(np.nanmax(mstat['Max'][iy]))
        ystat['Min'].append(np.nanmin(mstat['Min'][iy]))

    for y in ystat:
        ystat[y] = np.array(ystat[y])
        
    return ystat                                            


def annual_monthly_stats(y1 = 1991, y2 = 2020, use_wmo = False, 
                         use_txt = False,
                         sdir = 'd:/2021/vertailukausi30/Data/Stats/',
                         outdir = 'd:/tulos/'):
# Calculate the 30-year annual/monthly statistics from the files
# produced by process_stats()

    if use_txt:
        ttxt = '_tekstikanta'
    else:
        ttxt = ''

    ystr = '_' + str(y1) + '_' + str(y2)
    
    for mar in 'aobpvsmrtdhef':
        msid = fmisids[mar]
        for hsyst in ('MW', 'N2000'):
            if not use_wmo:
                ystat = read_statfile(sdir + mareo_names[mar] + '_' + \
                                 'annualstats' + \
                      '_' + hsyst + ttxt + '.txt', 
                      y1 = y1, y2 = y2, ismonthly = False)
            mstat = read_statfile(sdir + mareo_names[mar] + '_' + \
                                      'monthlystats' + \
                      '_' + hsyst + ttxt + ('_wmo' if use_wmo else '') \
                          + '.txt', y1 = y1, y2 = y2, ismonthly = True)
        
            sout = 'Month\tLowest mean\tMean\tHighest mean\tStd(mean)\tN(mean)' + \
            '\tLowest max\tAverage max\tHighest max\tStd(max)\tN(max)' + \
            '\tLowest min\tAverage min\tHighest min\tStd(min)\tN(min)\tYear(max)\tYear(min)\n'

            sout_kantaan = 'FMISID;Month;Lowest_mean;Mean;Highest_mean' + \
            ';Average_max;Highest_max;Lowest_min;Average_min;Year_max;Year_min\n'

            y_mnstat = [[], [], [], [], [], [], []]
            y_maxstat = [[], [], [], [], [], [], []]
            y_minstat = [[], [], [], [], [], [], []]
            
            for mon in np.arange(1, 13, 1):

                ik = (mstat['Mon'] == mon)
                mns = mstat['Mean'][ik] / 10.0
                maxs = mstat['Max'][ik] / 10.0
                mins = mstat['Min'][ik] / 10.0
                yrs = mstat['Yr'][ik]
                    
                mnstat = stats(mns, 24)
                maxstat = stats(maxs, 24)
                maxyear = yrs[maxstat[5]][0] if len(maxstat[5]) > 0 else float('nan')
                minstat = stats(mins, 24)
                minyear = yrs[minstat[6]][0] if len(minstat[6]) > 0 else float('nan')
                if use_wmo:
                    if mnstat[4] < 24:
                    # WMO: ääriarvon normaalia ei lasketa jos
                    # keskiarvoparametria ei ole olemassa väh. 24 vuodelta
                        for i in range(4):
                            maxstat[i] = float('nan')
                            minstat[i] = float('nan')
                                                            
                sout += ('{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}' + \
                         '\t{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}\t' + \
                             '{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}\t{:.0f}\t{:.0f}\n').format(mon, mnstat[2], mnstat[0], \
                        mnstat[1], mnstat[3], mnstat[4], \
                        maxstat[2], maxstat[0], maxstat[1], \
                        maxstat[3], maxstat[4], \
                        minstat[2], minstat[0], minstat[1], \
                        minstat[3], minstat[4], maxyear, minyear)

                sout_kantaan += ('{:.0f};{:.0f};{:.0f};{:.0f};{:.0f}' + \
                         ';{:.0f};{:.0f};' + \
                             '{:.0f};{:.0f};{:.0f};{:.0f}\n').format(msid, mon, \
                        mnstat[2] * 10.0, mnstat[0] * 10.0, mnstat[1] * 10.0, \
                        maxstat[0] * 10.0, maxstat[1] * 10.0, \
                        minstat[2] * 10.0, minstat[0] * 10.0, \
                        maxyear, minyear)

                                                                                                
                for i in range(5):
                    y_mnstat[i].append(mnstat[i])
                    y_maxstat[i].append(maxstat[i])
                    y_minstat[i].append(minstat[i])

            if use_wmo:
            # Calculate annual stats:
                sout += ('{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}' + \
                     '\t{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}\t' + \
                     '{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}\t{:.0f}\t{:.0f}\n').format(13, \
                        float('nan'),
                np.mean(np.array(y_mnstat[0])),
                float('nan'),
                float('nan'),
                float('nan'),
                float('nan'),
                float('nan'),
                np.max(np.array(y_maxstat[1])),
                float('nan'),
                float('nan'),
                np.min(np.array(y_minstat[2])),
                float('nan'),
                float('nan'),
                float('nan'),
                float('nan'),
                float('nan'),
                float('nan'))
            else:
                mns = ystat['Mean'] / 10.0
                maxs = ystat['Max'] / 10.0
                mins = ystat['Min'] / 10.0
                    
                mnstat = stats(mns, 24)
                maxstat = stats(maxs, 24)
                minstat = stats(mins, 24)
                maxyear = yrs[maxstat[5]][0] if len(maxstat[5]) > 0 else float('nan')
                minyear = yrs[minstat[6]][0] if len(minstat[6]) > 0 else float('nan')
                                            
                sout += ('{:.0f}\t{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}' + \
                         '\t{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}\t' + \
                             '{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}\t{:.0f}\t{:.0f}\n').format(13, mnstat[2], mnstat[0], \
                        mnstat[1], mnstat[3], mnstat[4], \
                        maxstat[2], maxstat[0], maxstat[1], \
                        maxstat[3], maxstat[4], \
                        minstat[2], minstat[0], minstat[1], \
                        minstat[3], minstat[4], maxyear, minyear)

                sout_kantaan += ('{:.0f};{:.0f};{:.0f};{:.0f};{:.0f}' + \
                         ';{:.0f};{:.0f};' + \
                             '{:.0f};{:.0f};{:.0f};{:.0f}\n').format(msid, 13, 
                        mnstat[2] * 10.0, mnstat[0] * 10.0, mnstat[1] * 10.0, \
                        maxstat[0] * 10.0, maxstat[1] * 10.0, \
                        minstat[2] * 10.0, minstat[0] * 10.0, \
                        maxyear, minyear)

                                                                                                                                                                            
            with open(outdir + mareo_names[mar] + '_results_' + hsyst + \
                      ttxt + ('_wmo' if use_wmo else '') + \
                         ystr + '.txt', 'w') as f:
                f.write(sout)
            with open(outdir + mareo_names[mar] + '_results_' + hsyst + \
                      ttxt + ('_wmo' if use_wmo else '') + ystr + '_kantaan.txt', 'w') as f:
                f.write(sout_kantaan)


def pick_distribution_levels(lev, y1 = 1991, y2 = 2020, 
                             didir = 'd:/2021/vertailukausi30/data_from1891/jakaumat/'):
# Pick certain probability levels from the probability distributions
# produced by the Sealevel-program.
# lev: probability levels: either MAX, MIN, or percents

    # Read annual distribution
    alevels = read_distribution(didir + 'Vuosijakauma_' + str(y1) + '-' \
                                + str(y2) + '_mw.txt', 
                                monthly = False)
    afound = (len(alevels) > 0)

    # Note Föglö and Mäntyluoto differ from the usual definition in the file names.
    mnames = {'a': 'Kemi', 'o': 'Oulu', 'b': 'Raahe', 'p': 'Pietarsaari', 
           'v': 'Vaasa', 's': 'Kaskinen', 'm': 'Mäntyluoto', 
           'r': 'Rauma', 't': 'Turku', 'd': 'Föglö', 'h': 'Hanko', 
           'e': 'Helsinki', 'f': 'Hamina', 'c': 'Porvoo'}

    # Tide gauge columns in the annual distribution table.
    mnro = {'a': 1, 'o': 2, 'b': 3, 'p': 4, 'v': 5, 's': 6, 'm': 7,
            'r': 8, 't': 9, 'd': 10, 'h': 11, 'e': 12, 'c': 13, 'f': 14}
        
    for mar in 'aobpvsmrtdhef':
        
        msid = fmisids[mar]

        outt = 'FMISID;Pros'
        for i in range(12):
            outt += ';M{:02d}'.format(i + 1)
        outt += ';All\n'

        # Read tide gauge specific monthly distribution
        levels = read_distribution(didir + mnames[mar] + \
                               '_jakauma_' + str(y1) + '-' \
                                + str(y2) + '_mw.txt', monthly = True)
        mfound = (len(levels) > 0)
            

        # Pick the probability levels        
        for l in lev:
            mdata = levels[l] if mfound else np.nan * np.ones(12)
            adata = alevels[l][mnro[mar] - 1] if afound else np.nan
            outt += '{:.0f};'.format(msid)
            outt += (l if l in ('MAX', 'MIN', 'MHW', 'MLW') \
                     else '{:.1f}'.format(l))
            for m in mdata:
                outt += ';{:.0f}'.format(m * 10.0)
            outt += ';{:.0f}\n'.format(adata * 10.0)

        # Write results to file
        with open('d:/tulos/' + mareo_names[mar] + '_distribution_' \
                  + str(y1) + '_' + str(y2) + '_kantaan.txt', 'w') as f:
            f.write(outt)


############### Visualisation of the results: ###############################

def plot_timeseries(yr, mar, hsyst = 'MW', hourly_from_txt = False,
                    ddir = 'd:/2021/Vertailukausi30/Data/Raw/',
                    ddir2 = 'd:/2021/vertailukausi30/Data_other/Raw_txt/',
                    tdir = 'd:/2021/vertailukausi30/Data_other/Tekstikanta/',
                    outdir = 'd:/tulos/'):
# Plot one-year timeseries of hourly values, monthly statistics and annual 
# statistics, for quality checking and comparison.

    if hourly_from_txt:
        hdata = read_hourlydata(ddir2 + mareo_names[mar] + '_all_' + hsyst + '.txt')
    else:
        hdata = read_cldbdata(ddir + 'sealevels' + \
                                    str(fmisids[mar]) + '.txt',
                                    fmisids[mar])
    ydata = read_cldbdata(ddir + 'annual' + \
                                    str(fmisids[mar]) + '.txt',
                                    fmisids[mar])
    mdata = read_cldbdata(ddir + 'monthly' + \
                                    str(fmisids[mar]) + '.txt',
                                    fmisids[mar])

    adata = read_txtdata(tdir + mareo_names[mar] + '_aariarvot_ref.txt', 
                                 aari = True)
    kdata = read_txtdata(tdir + mareo_names[mar] + '_keskiarvot_ref.txt', 
                                 aari = False)

    mw = read_mw('d:/2021/Vertailukausi30/Data/mwtaulu_1887-2021.txt')

    reflev = N2000mar[mar] if hsyst == 'N2000' else mw[mar][mw['Year'] == yr]

    fig, ax = plt.subplots(figsize = [10, 8])

    kerr = 24 * 3600   # x-akseli vuorokausina
    
    slname = 'Sl' + ('' if hourly_from_txt else hsyst)
    
    t1 = timegm(strptime(str(yr) + '-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
    t2 = timegm(strptime(str(yr + 1) + '-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
    ih = (hdata['Tstamp'] >= t1) & (hdata['Tstamp'] < t2)
    ax.plot((hdata['Tstamp'][ih] - t1)/kerr, 
            hdata[slname][ih], 'k-')

    iy = (ydata['Tstamp'] == t1)
    yts = ydata['Tstamp'][iy]
    ysl = ydata['Sl' + hsyst][iy]
    ax.plot(([(yts - t1)/kerr, (yts - t1)/kerr + 365]), 
            [ysl, ysl], 'r-')

    im = (mdata['Tstamp'] >= t1) & (mdata['Tstamp'] < t2)
    mts = mdata['Tstamp'][im]
    msl = mdata['Sl' + hsyst][im]
    for i in range(0, len(mts)):
        ax.plot([(mts[i] - t1)/kerr, (mts[i] - t1)/kerr + 30], 
                [msl[i], msl[i]], 'b-')

    iyt = (kdata['Yr'] == yr) & (kdata['Mon'] == 13)
    ytst = timegm(strptime(str(yr) + '-01-01T00:00:00', '%Y-%m-%dT%H:%M:%S'))
    yslt = kdata['Mean'][iyt] - reflev
    ax.plot(([(ytst - t1)/kerr, (ytst - t1)/kerr + 365]), 
            [yslt, yslt], 'm--')
    iyt2 = (adata['Yr'] == yr) & (adata['Mon'] == 13)
    yslt2 = adata['Min'][iyt2] - reflev
    yslt3 = adata['Max'][iyt2] - reflev
    ax.plot(([(ytst - t1)/kerr, (ytst - t1)/kerr + 365]), 
            [yslt2, yslt2], 'm--')
    ax.plot(([(ytst - t1)/kerr, (ytst - t1)/kerr + 365]), 
            [yslt3, yslt3], 'm--')

    for mon in np.arange(1,13):
        imt = (kdata['Yr'] == yr) & (kdata['Mon'] == mon)
        mtst = timegm(strptime(str(yr) + '-{:02d}-01T00:00:00'.format(mon), 
                               '%Y-%m-%dT%H:%M:%S'))
        mslt = kdata['Mean'][imt] - reflev
        ax.plot(([(mtst - t1)/kerr, (mtst - t1)/kerr + 30]), 
            [mslt, mslt], 'g--')
        imt2 = (adata['Yr'] == yr) & (adata['Mon'] == mon)
        mslt2 = adata['Min'][imt2] - reflev
        mslt3 = adata['Max'][imt2] - reflev
        ax.plot(([(mtst - t1)/kerr, (mtst - t1)/kerr + 30]), 
            [mslt2, mslt2], 'g--')
        ax.plot(([(mtst - t1)/kerr, (mtst - t1)/kerr + 30]), 
            [mslt3, mslt3], 'g--')

    ax.grid()
    ax.set_title(mareo_names[mar])
    ax.set_ylabel('Sea level (mm, ' + hsyst + ')')
    ax.set_xlabel('Days (' + str(yr) + ')')
    fig.savefig(outdir + mareo_names[mar] + str(yr) + 'timeseries.png', dpi = 300)
    plt.close(fig)
    
    
def plot_statistics(mar, hsyst, y1 = 1991, y2 = 2020, use_wmo = False,
                    datadir = 'd:/2021/Vertailukausi30/Data_from1891/Results/',
                    outdir = 'd:/tulos/'):
    
    ystr = '_' + str(y1) + '_' + str(y2)
    res = read_resultfile(datadir + mareo_names[mar] + '_results_' + \
                          hsyst + ('_wmo' if use_wmo else '') + ystr + '.txt')    

    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.9, 'bottom': 0.15})

    plotsym = {'Mean': 'x-',
               'Highest mean': 'x--',
               'Lowest mean': 'x--',
               'Average max': '.-',
               'Highest max': '.--',
               'Average min': '.-',
               'Lowest min': '.--'}
    
    plotcol = {'Mean': '#303193',
               'Highest mean': '#303193',
               'Lowest mean': '#303193',
               'Average max': '#02B8CE',
               'Highest max': '#02B8CE',
               'Average min': '#3A66E3',
               'Lowest min': '#3A66E3'}

    plotlabel = {'Mean': '$M$',
               'Highest mean': '$M_{max}$',
               'Lowest mean': '$M_{min}$',
               'Average max': 'MHW',
               'Highest max': 'MAX',
               'Average min': 'MLW',
               'Lowest min': 'MIN'}

    plotmarname = {'a': 'Kemi',
                 'o': 'Oulu / Uleåborg',
                 'b': 'Raahe / Brahestad',
                 'p': 'Pietarsaari / Jakobstad',
                 'v': 'Vaasa / Vasa',
                 's': 'Kaskinen / Kaskö',
                 'm': 'Mäntyluoto',
                 'r': 'Rauma / Raumo',
                 't': 'Turku / Åbo',
                 'd': 'Föglö',
                 'h': 'Hanko / Hangö',
                 'e': 'Helsinki / Helsingfors',
                 'f': 'Hamina / Fredrikshamn'}

    im = (res['Month'] < 13)
    xmon = res['Month'][im]

    for ps in plotsym:
        ax.plot(xmon, res[ps][im], plotsym[ps], color = plotcol[ps], label = plotlabel[ps])

    ax.grid()
    ax.legend(ncol = 3, loc = 'upper center')
    ax.set_title(plotmarname[mar] + ', ' + hsyst.upper())
    ax.set_ylabel('cm, ' + hsyst.upper())
    ax.set_xlabel('Kuu / Månad / Month')
    
    fig.savefig(outdir + mareo_names[mar] + '_stats_' + hsyst.upper() + \
                ('_wmo' if use_wmo else '') + ystr + '.png', dpi = 300)
    fig.savefig(outdir + mareo_names[mar] + '_stats_' + hsyst.upper() + \
                ('_wmo' if use_wmo else '') + ystr + '.eps')
    plt.close(fig)


def plot_annualstat(hsyst, 
                    datadir = 'd:/2021/Vertailukausi30/Data/Results/',
                    outdir = 'd:/tulos/'):

    plotsym = {'Mean': 'x-',
               'Highest mean': 'x--',
               'Lowest mean': 'x--',
               'Average max': '.-',
               'Highest max': '.--',
               'Average min': '.-',
               'Lowest min': '.--'}
    
    plotcol = {'Mean': '#303193',
               'Highest mean': '#303193',
               'Lowest mean': '#303193',
               'Average max': '#02B8CE',
               'Highest max': '#02B8CE',
               'Average min': '#3A66E3',
               'Lowest min': '#3A66E3'}

    plotlabel = {'Mean': '$M$',
               'Highest mean': '$M_{max}$',
               'Lowest mean': '$M_{min}$',
               'Average max': 'MHW',
               'Highest max': 'MAX',
               'Average min': 'MLW',
               'Lowest min': 'MIN'}

    plotmarname = {'a': 'Kemi',
                 'o': 'Oulu\nUleåborg',
                 'b': 'Raahe\nBrahestad',
                 'p': 'Pietarsaari\nJakobstad',
                 'v': 'Vaasa\nVasa',
                 's': 'Kaskinen\nKaskö',
                 'm': 'Pori\nBjörneborg',
                 'r': 'Rauma\nRaumo',
                 't': 'Turku\nÅbo',
                 'd': 'Föglö',
                 'h': 'Hanko\nHangö',
                 'e': 'Helsinki\nHelsingfors',
                 'f': 'Hamina\nFredrikshamn'}

    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

    xmar = []
    ydata = {}
    for ps in plotsym:
        ydata[ps] = []
        
    for mar in 'aobpvsmrtdhef':
        res = read_resultfile(datadir + mareo_names[mar] + '_results_' + \
                          hsyst + '.txt')    

        im = (res['Month'] == 13)
        for ps in plotsym:
            ydata[ps].append(res[ps][im])
        xmar.append(plotmarname[mar])

    fig, ax = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.95, 'bottom': 0.2})

    for ps in plotsym:
        ax.plot(xmar, ydata[ps], plotsym[ps], color = plotcol[ps], 
                label = plotlabel[ps])

    plt.xticks(rotation = 90)
    ax.grid()
    ax.legend(ncol = 3)
    ax.set_ylabel('cm, ' + hsyst.upper())
    fig.savefig(outdir + 'All_stats_' + hsyst.upper() + '.png', dpi = 300)
    fig.savefig(outdir + 'All_stats_' + hsyst.upper() + '.eps')
    plt.close(fig)





def tabulate_statistics(use_wmo = False, 
                    datadir = 'd:/2021/Vertailukausi30/Data/Results/',
                    outdir = 'd:/tulos/'):

    ttext_all = '<p>'

    for hsyst in ('Mw', 'N2000'):
        ttext_year = '<table>\n'
        ttext_year += '<tr><td colspan=10>Korkeusjärjestelmä/höjdsystem/height system ' + hsyst.upper() + '</td></tr>\n'
        ttext_year += '<tr><td></td><td colspan=3>Keskiarvot<br>Medelvärden<br>Averages</td>'
        ttext_year += '<td colspan=6>Ääriarvot<br>Extremvärden<br>Extremes</td></tr>\n'
        ttext_year += '<tr><td>Mareografi<br>Mareograf<br>Tide gauge</td>'
        ttext_year += '<td>Alin<br>Lägsta<br>Lowest</td><td>Keskim.<br>I medeltal<br>Average</td>'
        ttext_year += '<td>Ylin<br>Högsta<br>Highest</td>'
        ttext_year += '<td>Keskim. ylin<br>I medeltal högsta<br>Average max</td>'
        ttext_year += '<td>Keskim. alin<br>I medeltal lägsta<br>Average min</td>'
        ttext_year += '<td>Absol. ylin<br>Absol. högsta<br>Absol. max</td>'
        ttext_year += '<td>Vuosi<br>År<br>Year</td>'
        ttext_year += '<td>Absol. alin<br>Absol. lägsta<br>Absol. min</td>'
        ttext_year += '<td>Vuosi<br>År<br>Year</td></tr>'
        
        for mar in 'aobpvsmrtdhef':

            res = read_resultfile(datadir + mareo_names[mar] + '_results_' + \
                              hsyst + ('_wmo' if use_wmo else '') + '.txt')    
        
            ttext = '<table>\n'
            ttext += '<tr><td colspan=10>' + mareo_names_tables[mar] + \
                ', korkeusjärjestelmä/höjdsystem/height system ' + hsyst.upper() + '</td></tr>\n'
            ttext += '<tr><td></td><td colspan=3>Keskiarvot<br>Medelvärden<br>Averages</td>'
            ttext += '<td colspan=6>Ääriarvot<br>Extremvärden<br>Extremes</td></tr>\n'
            ttext += '<tr><td>Kk<br>Månad<br>Month</td>'
            ttext += '<td>Alin<br>Lägsta<br>Lowest</td><td>Keskim.<br>I medeltal<br>Average</td>'
            ttext += '<td>Ylin<br>Högsta<br>Highest</td>'
            ttext += '<td>Keskim. ylin<br>I medeltal högsta<br>Average max</td>'
            ttext += '<td>Keskim. alin<br>I medeltal lägsta<br>Average min</td>'
            ttext += '<td>Absol. ylin<br>Absol. högsta<br>Absol. max</td>'
            ttext += '<td>Vuosi<br>År<br>Year</td>'
            ttext += '<td>Absol. alin<br>Absol. lägsta<br>Absol. min</td>'
            ttext += '<td>Vuosi<br>År<br>Year</td></tr>'
        
            # tsym = ('Lowest mean', 'Mean', 'Highest mean',
            #            'Average max', 'Average min', 'Highest max', 'Maxyear', 
            #            'Lowest mean', 'Minyear')
            tsym = ('Lowest mean', 'Mean', 'Highest mean',
                       'Average max', 'Average min', 'Highest max', 'Year(max)', 
                       'Lowest min', 'Year(min)')
        
            for mon in np.arange(1,14,1):
                
                trow = '<tr><td>' + ('Vuosi<br>År<br>Year' if mon == 13 else '{:02d}'.format(mon)) + '</td>'
                for sym in tsym:
                    trow += '<td>{:.0f}</td>'.format(res[sym][mon - 1])
                trow += '</tr>'
                ttext += trow
                if mon == 13:
                    ttext_year += '<tr><td>' + mareo_names[mar] + '</td>'
                    for sym in tsym:
                        ttext_year += '<td>{:.0f}</td>'.format(res[sym][mon - 1])
                    ttext_year += '</tr>'
            ttext += '</table>'
            
            with open(outdir + mareo_names[mar] + '_table1_' + hsyst.upper() + \
                        ('_wmo' if use_wmo else '') + '.html', 'w') as f:
                f.write(ttext)
                
            ttext_all += ttext + '</p><p>\n'

        ttext_year += '</table>\n'

        with open(outdir + 'Annual_table1_' + hsyst.upper() + \
                    ('_wmo' if use_wmo else '') + '.html', 'w') as f:
            f.write(ttext_year)
        
    ttext_all += '</p>\n'

    with open(outdir + 'All_tables1' + \
              ('_wmo' if use_wmo else '') + '.html', 'w') as f:
        f.write(ttext_all)


def plot_all_stats(hsyst = 'mw', y1 = 1991, y2 = 2020, use_wmo = False):
    
    for mar in 'aobpvsmrtdhef':
        plot_statistics(mar, hsyst, y1 = y1, y2 = y2, use_wmo = use_wmo)




def plot_distributions(mar, y1 = 1991, y2 = 2020,
                    datadir = 'd:/2021/Vertailukausi30/Data/Jakaumat/',
                    outdir = 'd:/tulos/'):

    ystr = '_' + str(y1) + '_' + str(y2)
    res = read_distribution2(datadir + mareo_names[mar] + ystr + '_distribution.txt')    


    plotsym = {'MAX': '+-',
               0.1: '.:',
               1.0: '.--',
               5.0: 'o--',
               25.0: 'x:',
               50.0: 'o-',
               75.0: 'x:',
               95.0: 'o--',
               99.0: '.--',
               99.9: '.:',
               'MIN': '+-'}

    plotcol = {'MAX': '#52C1A1',
               0.1: '#02B8CE',
               1.0: '#3A66E3',
               5.0: '#3A66E3',
               25.0: '#303193',
               50.0: '#303193',
               75.0: '#303193',
               95.0: '#3A66E3',
               99.0: '#3A66E3',
               99.9: '#02B8CE',
               'MIN': '#52C1A1'}

    plotlabel = {'MAX': 'MAX',
               0.1: '0.1 %',
               1.0: '1 %',
               5.0: '5 %',
               25.0: '25 %',
               50.0: '50 %',
               75.0: '75 %',
               95.0: '95 %',
               99.0: '99 %',
               99.9: '99.9 %',
               'MIN': 'MIN'}

    plotmarname = {'a': 'Kemi',
                 'o': 'Oulu / Uleåborg',
                 'b': 'Raahe / Brahestad',
                 'p': 'Pietarsaari / Jakobstad',
                 'v': 'Vaasa / Vasa',
                 's': 'Kaskinen / Kaskö',
                 'm': 'Mäntyluoto',
                 'r': 'Rauma / Raumo',
                 't': 'Turku / Åbo',
                 'd': 'Föglö',
                 'h': 'Hanko / Hangö',
                 'e': 'Helsinki / Helsingfors',
                 'f': 'Hamina / Fredrikshamn'}

    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.9, 'bottom': 0.15})

    xmon = np.arange(1, 13, 1)
    
    for ps in plotsym:
        ax.plot(xmon, res[ps][:-1], plotsym[ps], color = plotcol[ps], 
                label = plotlabel[ps])

    ax.grid()
    ax.legend(ncol = 4, loc = 'upper center')
    ax.set_ylabel('cm, MW')

    ax.set_title(plotmarname[mar])
    ax.set_xlabel('Kuu / Månad / Month')
    fig.savefig(outdir + mareo_names[mar] + '_distribution.png', dpi = 300)
    fig.savefig(outdir + mareo_names[mar] + '_distribution.eps')
    plt.close(fig)


def plot_all_distributions():
    
    for mar in 'aobpvsmrtdhef':
        plot_distributions(mar)


def plot_annualdistr(datadir = 'd:/2021/Vertailukausi30/Data/Jakaumat/',
                    outdir = 'd:/tulos/'):

    plotsym = {'MAX': '+-',
               0.1: '.:',
               1.0: '.--',
               5.0: 'o--',
               25.0: 'x:',
               50.0: 'o-',
               75.0: 'x:',
               95.0: 'o--',
               99.0: '.--',
               99.9: '.:',
               'MIN': '+-'}

    plotcol = {'MAX': '#52C1A1',
               0.1: '#02B8CE',
               1.0: '#3A66E3',
               5.0: '#3A66E3',
               25.0: '#303193',
               50.0: '#303193',
               75.0: '#303193',
               95.0: '#3A66E3',
               99.0: '#3A66E3',
               99.9: '#02B8CE',
               'MIN': '#52C1A1'}

    plotlabel = {'MAX': 'MAX',
               0.1: '0.1 %',
               1.0: '1 %',
               5.0: '5 %',
               25.0: '25 %',
               50.0: '50 %',
               75.0: '75 %',
               95.0: '95 %',
               99.0: '99 %',
               99.9: '99.9 %',
               'MIN': 'MIN'}

    plotmarname = {'a': 'Kemi',
                 'o': 'Oulu\nUleåborg',
                 'b': 'Raahe\nBrahestad',
                 'p': 'Pietarsaari\nJakobstad',
                 'v': 'Vaasa\nVasa',
                 's': 'Kaskinen\nKaskö',
                 'm': 'Pori\nBjörneborg',
                 'r': 'Rauma\nRaumo',
                 't': 'Turku\nÅbo',
                 'd': 'Föglö',
                 'h': 'Hanko\nHangö',
                 'e': 'Helsinki\nHelsingfors',
                 'f': 'Hamina\nFredrikshamn'}


    xmar = []
    ydata = {}
    for ps in plotsym:
        ydata[ps] = []
        
    for mar in 'aobpvsmrtdhef':
        
        res = read_distribution2(datadir + mareo_names[mar] + '_distribution.txt')    

        for ps in plotsym:
            ydata[ps].append(res[ps][12])
        xmar.append(plotmarname[mar])


    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.95, 'bottom': 0.2})

    for ps in plotsym:
        ax.plot(xmar, ydata[ps], plotsym[ps], color = plotcol[ps], label = plotlabel[ps])
    
    plt.xticks(rotation = 90)
    ax.grid()
    ax.legend(ncol = 4)
    ax.set_ylabel('cm, MW')

#    ax.set_title('All tide gauges, annual distribution')
#    ax.set_xlabel('Tide gauge')
    fig.savefig(outdir + 'All_annual_distribution.png', dpi = 300)
    fig.savefig(outdir + 'All_annual_distribution.eps')
    plt.close(fig)

################################################################
### Extra: plot comparisons between different 30-year periods
################################################################

def plot_annualstat_change(hsyst, mars = 'avdhe', 
                    datadir = 'd:/2021/Vertailukausi30/Data_from1891/Results/',
                    outdir = 'd:/tulos/'):

    plotsym = {'Mean': 'x-',
               'Highest mean': 'x--',
               'Lowest mean': 'x--',
               'Average max': '.-',
               'Highest max': '.--',
               'Average min': '.-',
               'Lowest min': '.--'}
    
#    plotcol = {'Mean': '#303193',
#               'Highest mean': '#303193',
#               'Lowest mean': '#303193',
#               'Average max': '#02B8CE',
#               'Highest max': '#02B8CE',
#               'Average min': '#3A66E3',
#               'Lowest min': '#3A66E3'}

    plotcol = {1891: 'black',
               1901: 'blue',
               1911: 'cyan',
               1921: 'green',
               1931: 'lightgreen',
               1941: 'gray',
               1951: 'brown',
               1961: 'yellow',
               1971: 'orange',
               1981: 'red',
               1991: 'magenta'}
    
    # plotlabel = {'Mean': '$M$',
    #            'Highest mean': '$M_{max}$',
    #            'Lowest mean': '$M_{min}$',
    #            'Average max': 'MHW',
    #            'Highest max': 'MAX',
    #            'Average min': 'MLW',
    #            'Lowest min': 'MIN'}

    plotmarname = {'a': 'Kemi',
                 'o': 'Oulu\nUleåborg',
                 'b': 'Raahe\nBrahestad',
                 'p': 'Pietarsaari\nJakobstad',
                 'v': 'Vaasa\nVasa',
                 's': 'Kaskinen\nKaskö',
                 'm': 'Pori\nBjörneborg',
                 'r': 'Rauma\nRaumo',
                 't': 'Turku\nÅbo',
                 'd': 'Föglö',
                 'h': 'Hanko\nHangö',
                 'e': 'Helsinki\nHelsingfors',
                 'f': 'Hamina\nFredrikshamn'}

    plotcolmar = {'a': 'black',
                 'o': 'blue',
                 'b': 'red',
                 'p': 'green',
                 'v': 'magenta',
                 's': 'orange',
                 'm': 'black',
                 'r': 'blue',
                 't': 'red',
                 'd': 'green',
                 'h': 'blue',
                 'e': 'orange',
                 'f': 'brown'}

    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.95, 'bottom': 0.2})

    tseries = {}
    for mar in 'aobpvsmrtdhef':
        tseries[mar] = {}
        for ps in plotsym:
            tseries[mar][ps] = []

    hl = []
    hlabel = []
    for y1 in np.arange(1891, 2001, 10):

        xmar = []
        ydata = {}
        for ps in plotsym:
            ydata[ps] = []
        
        for mar in 'aobpvsmrtdhef':
            res = read_resultfile(datadir + mareo_names[mar] + '_results_' + \
                          hsyst + '_' + str(y1) + '_' + str(y1 + 29) + '.txt')

            im = (res['Month'] == 13)
            for ps in plotsym:
                ydata[ps].append(res[ps][im])
            xmar.append(plotmarname[mar])
            

            for ps in plotsym:
                tseries[mar][ps].append(res[ps][im])
                hh, = ax.plot(xmar, ydata[ps], plotsym[ps], color = plotcol[y1], 
                    label = str(y1))
                if ps == 'Mean':
                    if mar == 'a':
                        hl.append(hh)
                        hlabel.append(str(y1))

    plt.xticks(rotation = 90)
    ax.grid()
    ax.legend(handles = hl, labels = hlabel, ncol = 3)
    ax.set_ylabel('cm, ' + hsyst.upper())
    fig.savefig(outdir + 'All_stats_' + hsyst.upper() + '_change.png', dpi = 300)
#    fig.savefig(outdir + 'All_stats_' + hsyst.upper() + '_change.eps')
    plt.close(fig)

    fig2, ax2 = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.95, 'bottom': 0.2})
    for mar in mars:
        for sym in ['Highest max', 'Average max', 'Mean', 'Average min', 'Lowest min']:
            ax2.plot(np.arange(1891, 2001, 10), tseries[mar][sym], 
                     plotsym[sym],
                     color = plotcolmar[mar],
                     label = mareo_names[mar] + ', ' + sym)

    ax2.grid()
    ax2.legend()
    ax2.set_xlabel('Year')
    ax2.set_ylabel('cm, ' + hsyst)
    fig2.savefig(outdir + 'Timeseries_' + hsyst.upper() + '_change_' + mars + '.png', dpi = 300)
    plt.close(fig2)


def plot_annualdistr_change(mars = 'avdhe', datadir = 'd:/2021/Vertailukausi30/Data_from1891/Jakaumat/',
                    outdir = 'd:/tulos/'):

    plotsym = {'MAX': '+-',
               0.1: '.:',
               1.0: '.--',
               5.0: 'o--',
               25.0: 'x:',
               50.0: 'o-',
               75.0: 'x:',
               95.0: 'o--',
               99.0: '.--',
               99.9: '.:',
               'MIN': '+-'}

    plotsym = {0.1: '.:',
               1.0: '.--',
               5.0: 'o--',
               25.0: 'x:',
               50.0: 'o-',
               75.0: 'x:',
               95.0: 'o--',
               99.0: '.--',
               99.9: '.:'}


    # plotcol = {'MAX': '#52C1A1',
    #            0.1: '#02B8CE',
    #            1.0: '#3A66E3',
    #            5.0: '#3A66E3',
    #            25.0: '#303193',
    #            50.0: '#303193',
    #            75.0: '#303193',
    #            95.0: '#3A66E3',
    #            99.0: '#3A66E3',
    #            99.9: '#02B8CE',
    #            'MIN': '#52C1A1'}

    plotcol = {1891: 'black',
               1901: 'blue',
               1911: 'cyan',
               1921: 'green',
               1931: 'lightgreen',
               1941: 'gray',
               1951: 'brown',
               1961: 'yellow',
               1971: 'orange',
               1981: 'red',
               1991: 'magenta'}
    

    plotcolmar = {'a': 'black',
                 'o': 'blue',
                 'b': 'red',
                 'p': 'green',
                 'v': 'magenta',
                 's': 'orange',
                 'm': 'black',
                 'r': 'blue',
                 't': 'red',
                 'd': 'green',
                 'h': 'blue',
                 'e': 'orange',
                 'f': 'brown'}

    plotlabel = {'MAX': 'MAX',
               0.1: '0.1 %',
               1.0: '1 %',
               5.0: '5 %',
               25.0: '25 %',
               50.0: '50 %',
               75.0: '75 %',
               95.0: '95 %',
               99.0: '99 %',
               99.9: '99.9 %',
               'MIN': 'MIN'}

    plotmarname = {'a': 'Kemi',
                 'o': 'Oulu\nUleåborg',
                 'b': 'Raahe\nBrahestad',
                 'p': 'Pietarsaari\nJakobstad',
                 'v': 'Vaasa\nVasa',
                 's': 'Kaskinen\nKaskö',
                 'm': 'Pori\nBjörneborg',
                 'r': 'Rauma\nRaumo',
                 't': 'Turku\nÅbo',
                 'd': 'Föglö',
                 'h': 'Hanko\nHangö',
                 'e': 'Helsinki\nHelsingfors',
                 'f': 'Hamina\nFredrikshamn'}

    tseries = {}
    for mar in 'aobpvsmrtdhef':
        tseries[mar] = {}
        for ps in plotsym:
            tseries[mar][ps] = []

    xmar = []
    for mar in 'aobpvsmrtdhef':
        xmar.append(plotmarname[mar])

    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.95, 'bottom': 0.2})

    hl = []
    hlabel = []
    for y1 in np.arange(1891, 2001, 10):

        ydata = {}
        for ps in plotsym:
            ydata[ps] = []
        
        for mar in 'aobpvsmrtdhef':

            res = read_distribution2(datadir + mareo_names[mar] + \
                    '_distribution_' + str(y1) + '_' + str(y1 + 29) + \
                        '_kantaan.txt')    
            
            for ps in plotsym:
                tseries[mar][ps].append(res[ps][12] / 10.0)
                ydata[ps].append(res[ps][12] / 10.0)
            
        for ps in plotsym:

            hh, = ax.plot(xmar, ydata[ps], plotsym[ps], color = plotcol[y1], 
                          label = str(y1))
            if ps == 0.5:
                if mar == 'a':
                    hl.append(hh)
                    hlabel.append(str(y1))
    
    plt.xticks(rotation = 90)
    ax.grid()
    ax.legend(handles = hl, labels = hlabel, ncol = 3)
    ax.set_ylabel('cm, MW')

    fig.savefig(outdir + 'All_annual_distribution_change.png', dpi = 300)
#    fig.savefig(outdir + 'All_annual_distribution_change.eps')
    plt.close(fig)
     
    fig2, ax2 = plt.subplots(figsize = [8, 5.5], gridspec_kw = {'left': 0.1, 'right': 0.95,
                                                              'top': 0.95, 'bottom': 0.2})
    for mar in mars:
        for sym in [0.1, 1.0, 5.0, 25.0, 50.0, 75.0, 95.0, 99.0, 99.9]:
            ax2.plot(np.arange(1891, 2001, 10), tseries[mar][sym], 
                     plotsym[sym],
                     color = plotcolmar[mar],
                     label = mareo_names[mar] + ', ' + str(sym))

    ax2.grid()
    ax2.legend()
    ax2.set_xlabel('Year')
    ax2.set_ylabel('cm, MW')
    fig2.savefig(outdir + 'Distribution_timeseries_change_' + mars + '.png', dpi = 300)
    plt.close(fig2)

    