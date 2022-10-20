# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 12:18:50 2022

@author: johanssm
"""

# Functions related to mean sea level calculations and definitions

from file_tools import loadfile

def hae_imw(mar, imwdir = 'd:/data/vedenk/keskivedet2021/', 
            imwname = '_imw_yhdistetty15_d3_1899_2020_l15.dat'):
    
    imw = loadfile(imwdir + mar + imwname)
    
    return imw

    
    
def read_mw(fname = 'teormw.dat'):
    # Luetaan teoreettisen keskiveden arvot tiedostosta.
    marn = []
    cs = []
    mws = []
    f = open('d:\\data\\vedenk\\' + fname, 'r')
    for row in f:
        rd = row.split()
        marn.append(rd[0])
        cs.append(rd[1:5])
        mws.append(rd[5:])
    f.close()
    nx = len(mws[0])
    mwdic = {}
    for i in range(nx):
        mm = []
        for j in range(1, 14):
            mm.append(float(mws[j][i]))            
#        d0 = {mws[0][i]: mm}
        mwdic[float(mws[0][i])] = mm
    return(mwdic)

def read_mw2(fname = 'mwtaulu_1887-2023.txt'):
    # Luetaan teoreettisen keskiveden arvot uudenmuotoisesta 
    # (sealevel-ohjelman tekemästä) tiedostosta.
    martunn = 'aobpvsmrtdhefc'

    with open('d:\\data\\vedenk\\' + fname, 'r') as f:
        mwdata = f.readlines()

    mwdic = {}
    for row in mwdata[3:]:
        rd = row.split()
        mwdic[float(rd[0])] = {}
        for i in range(14):
            mwdic[float(rd[0])][martunn[i]] = float(rd[i + 1])

    return(mwdic)

