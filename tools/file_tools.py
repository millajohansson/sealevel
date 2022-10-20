# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:07:34 2022

@author: johanssm
"""

# Tools for reading different types of sea level data files.

import numpy as np
from mareos import mareo_names, martunn
from datetime import datetime, timezone, timedelta


def loadfile(fname):
# Load a file into numpy array.
# Exclude commas and asterisks. (Tailored for sea level data files which 
# contain these as quality flags.)
# Lines starting with #, % or any alphabet are considered comments.

    tt = str.maketrans(',*', '  ')    # Muunnetaan pilkut ja tahdet välilyönneiksi
    data = []
    
    with open(fname,'r') as f:
        for row0 in f:
            row = row0.translate(tt).strip()
            if (len(row) > 0):
                if (row[0] not in ('%', '#', '"', "'")) & (not str.isalpha(row[0])):
                    rd = row.split()
                    drow = []
                    for rdi in rd:
                        drow.append(float(rdi))
                    data.append(drow)

    datam = np.array(data)
    return(datam)


def lue_max(mar, hakem = 'd:/data_v/vedenk/aariarvot/', q = 0):

# Luetaan aariarvotiedosto *_aariarvot_ref.txt mareografille mar.
# (Nama ovat siis vanhoja "tekstikannan" aariarvotiedostoja, kaytossa n. vuoteen
# 2020 asti.)
# mar: mareografi (kirjaintunnus esim. Kemi = 'a', tai järjestysnumero esim. Kemi = 1)
# hakem: hakemisto josta luetaan
# jos q>0, "quiet" eli ei tulosta otsikkorivia naytolle


    with open(hakem + '/' + mareo_names[mar] + '_aariarvot_ref.txt', 'r') as f:
        rivit = f.readlines()

    if q <= 0:
        print(rivit[0])

    yyr = []
    ymax = []
    ymin = []
    kkyr = []
    kkmon = []
    kkmax = []
    kkmin = []
    

    for r in range(1, len(rivit), 2):
        rivi1 = np.array([float(x.strip()) for x in rivit[r].replace(',', ' ').split()])
        rivi2 = np.array([float(x.strip()) for x in rivit[r + 1].replace(',', ' ').split()])
        drivi = np.concatenate((rivi1, rivi2))
        drivi[drivi < -990] = float('nan')

        yyr.append(drivi[0])
        ymax.append(drivi[13])
        ymin.append(drivi[27])
        for k in range(0, 12):
            kkyr.append(drivi[0])
            kkmon.append(k + 1)
            kkmax.append(drivi[k + 1])
            kkmin.append(drivi[k + 15])

    yyr = np.array(yyr)
    ymax = np.array(ymax)
    ymin = np.array(ymin)
    kkyr = np.array(kkyr)
    kkmon = np.array(kkmon)
    kkmax = np.array(kkmax)
    kkmin = np.array(kkmin)
    
    return yyr, ymax, ymin, kkyr, kkmon, kkmax, kkmin


def lue_kk3(hakem, marnro, v, kk, pkorj = False, tunn = '', q = False):
# Luetaan vedenkorkeuden kuukausitiedosto (ascii-rekisterista).
# Muokattu Matlab lue_kk3.m:sta.

# jos pkorj=True, tehdaan Pietarsaareen 15 mm korjaus vuoteen 1953 asti.
# tunn: jos ei annettu, oletetaan v (vuoteen 1970) tai i (v:sta 1971 eteenpain).
# jos q = True, ei tulosteta tiedoston otsikkorivia naytolle
# Puuttuvat = NaN

    kuupv = [31,28,31,30,31,30,31,31,30,31,30,31]
    if (v % 4 == 0 and v % 100 != 0) or v % 400 == 0:
        kuupv[1] = 29

    mar = martunn[marnro - 1]
    if v<1977 and marnro == 1:
        mar = 'k'

    if tunn == '':
        if v <= 1970:
            tunn = 'v'
        else:
            tunn = 'i'

    nimi = hakem + '/' + mar + '{:04d}{:02d}'.format(v, kk) + tunn + '.dat'

    values = []
    iflag = []
    pvs = []

    try:
        with open(nimi, 'r') as f:
            data0 = f.readlines()
    except IOError:
        data0 = []

    if any(data0):
        if not q:
            print(data0[0])
        data0 = data0[1:]  # poistetaan 1. rivi (otsikkorivi)
        
    lc=len(data0)
    
    if lc > 0:
        ic = 0
        while ic < lc - 1:
            rflag = []
            rval = []

            row1 = data0[ic].rstrip(' \n')
            row2 = data0[ic + 1].rstrip(' \n')
            if row1[-1] != ',':
                row1 += ' '
            if row2[-1] != ',':
                row2 += ' '
            row= row1 + row2[3:]
            try:
                pv = int(row[0:3])
                if pv < 1 or pv > 31:
                    raise ValueError
            except ValueError:
                print('lue_kk3: virheellinen paivatieto!')
                return [], [], []

            pvs.append(pv)                
            ir = 3
            while ir < len(row) - 4:
                r2 = row[ir:ir+5]
                rflag.append(r2[-1] == ',')
                try:
                    val0 = float(r2[0:4])
                except ValueError:
                    print('lue_kk3: virhe tiedoston rakenteessa!')
                    return [], [], []
                if val0 < -990:
                    rval.append(float('nan'))
                else:
                    rval.append(val0)
                ir += 5
            values.append(rval)
            iflag.append(rflag)
            ic += 2

    pvs = np.array(pvs)
    values = np.array(values)
    iflag = np.array(iflag)

    # Pietarsaaren korjaus -15 mm
    if pkorj and marnro == 4 and v <= 1953:
        values -= 15
        
    return pvs, values, iflag
    

def read_hdata(fname):
    # Luetaan sealevel-ohjelman tuottama tuntidatatiedosto.

    yr = []
    mon = []
    day = []
    hr = []
    sealevel = []
    intp = []
    timesec = []
    
    with open(fname,'r') as f:
        for row0 in f:
            row = row0.strip()
            if (len(row) > 3):
                if (row[0] not in ('%', '#')) & (not str.isalpha(row[0])):
                    rd = [x.strip() for x in row.split()]
                    yr.append(float(rd[0]))
                    mon.append(float(rd[1]))
                    day.append(float(rd[2]))
                    hr.append(float(rd[3]))
                    intp.append(rd[4][-1] == ',')
                    sealevel.append(float(rd[4].strip(',')))
                    timesec.append(datetime(int(rd[0]), int(rd[1]), int(rd[2]), \
                                                int(rd[3]) - 1, \
                        tzinfo = timezone(timedelta(0, 3600))).timestamp())
#                    timesec.append(pvm2sec(int(rd[0]), int(rd[1]), int(rd[2]), \
#                                  int(rd[3]) - 1, 0.0, 0.0, 1.0))                
    return np.array(yr), np.array(mon), np.array(day), \
        np.array(hr), np.array(sealevel), np.array(intp), np.array(timesec)
    

def read_vrkdata(fname):
    # Luetaan sealevel-ohjelman tuottama vrk-datatiedosto.

    yr = []
    mon = []
    day = []
    smean = []
    smax = []
    smin = []
    sdev = []
    intp = []
    intdev = []
    intmax = []
    intmin = []
    
    with open(fname,'r') as f:
        for row0 in f:
            row = row0.strip()
            if (len(row) > 3):
                if (row[0] not in ('%', '#')) & (not str.isalpha(row[0])):
                    rd = [x.strip() for x in row.split()]
                    yr.append(float(rd[0]))
                    mon.append(float(rd[1]))
                    day.append(float(rd[2]))
                    intp.append(rd[3][-1] == ',')
                    smean.append(float(rd[3].strip(',')))
                    if len(rd) == 7:
                        intdev.append(rd[4][-1] == ',')
                        sdev.append(float(rd[4].strip(',')))
                        intmax.append(rd[5][-1] == ',')
                        smax.append(float(rd[5].strip(',')))
                        intmin.append(rd[6][-1] == ',')
                        smin.append(float(rd[6].strip(',')))
                    else:
                        intdev.append(False)
                        sdev.append(float('nan'))
                        intmax.append(rd[4][-1] == ',')
                        smax.append(float(rd[4].strip(',')))
                        intmin.append(rd[5][-1] == ',')
                        smin.append(float(rd[5].strip(',')))
                
    return np.array(yr), np.array(mon), np.array(day), \
        np.array(smean), np.array(sdev), np.array(smax), \
        np.array(smin), np.array(intp), np.array(intdev), np.array(intmax), \
        np.array(intmin)


