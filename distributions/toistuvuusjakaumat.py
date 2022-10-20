# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 15:04:11 2017

@author: johanssm
"""

# Calculate probability distributions for sea level maxima/minima.
#
# Partly translated from earlier Matlab scripts.
# 
# Function laske_ja_taulukoi is the easiest-to-use interface for these,
# it calls the other functions to do the actual calculations.


import matplotlib.pyplot as plt
import sys
import numpy as np

sys.path.append('../tools/') 

from file_tools import loadfile, lue_max
from mareos import mareo_names, N2000mar, mareo_names_ver2
#from weibulls import weibfit3, weibfit2, cweibull
from mw_tools import hae_imw, read_mw2
from expfits import cum_expdistr, expfit2

# Note: the Weibull fit functions (weibulls.py) have not been included in my 
# public repository, as they are not originally my code.
# The exponential fits can be used instead.
# If using weibulls.py, comment the dummy functions from below.
def weibfit3(x1, x2, x3, x4):
    return []
def weibfit2(x1, x2, x3, x4):
    return [], [], [], []
def cweibull(x1, x2):
    return []


# Some default directory definitions
# (These can be given as function parameters too)
# Results will be saved here:
resultdir = 'd:/tulos/'

# iMW (mean sea level estimate) default definitions.
# Directory and file name identifier for mean sea level estimate (iMW) files
# (the actual file names being "mar + imwname", where mar = tide gauge identifier,
# e.g. "o_imw_yhdistetty15_d3_1899_2020_l15.dat" for Oulu)
# (These can be given as function parameters too)
imwdir_def = 'd:/data/vedenk/keskivedet2021/'
imwname_def = '_imw_yhdistetty15_d3_1899_2020_l15.dat'

# Directory for monthly/annual maxima/minima obtained from Smartmet Server.
# (To obtain the data files, run datahakuja.py/get_smartmet_statistics)
smdir = 'D:/Data/Vedenk/Aari_keskiarvot/'


def laske_ja_taulukoi(minmax = 0, jakaumassa_n = False, 
                      vanhasovitus = False, n2000vuosi = 2011, avuosi = 1982,
                      lvuosi = 2011, tasot = [0.2, 0.1, 0.05, 0.02, 0.01, 0.004, 0.001],
                      dir_out = resultdir, lahde = 'smartmet', output_mw = False,
                      imwname = imwname_def, imwdir = imwdir_def):

# Lasketaan kaikille mareografeille minimien tai maksimien jakaumat halutuilta vuosilta.
# Interpoloidaan jakaumista vedenkorkeudet halutuille todennäköisyystasoille.
    
    for mar in 'aobpvsmrtdhef': 
        laske_maxjak(mar = mar, mm = minmax, avuosi = avuosi, lvuosi = lvuosi, 
                     n2000vuosi = n2000vuosi, 
                   jakaumassa_n = jakaumassa_n,
                   vanhasovitus = vanhasovitus, lahde = lahde, 
                   imwname = imwname, imwdir = imwdir)
    tloppu = ('_oldn' if jakaumassa_n else '') \
        + ('_oldfit' if vanhasovitus else '') \
        + ('_sm' if lahde == 'smartmet' else '')
    luvut = poimitasot(tasot = tasot, hakem = dir_out, minimit = (minmax == 1),
               avuosi = avuosi, lvuosi = lvuosi, n2000vuosi = n2000vuosi, 
               tdstoloppu = tloppu, output_mw = output_mw, 
               imwname = imwname, imwdir = imwdir)
    
    return luvut


def laske_maxjak(mar = 'e', mm = 0, avuosi = 1991, lvuosi = 2020, 
               aineisto = 'kk', n2000vuosi = 0.0,
               jakaumassa_n = False, vanhasovitus = False,
               dir_out = resultdir, lahde = 'smartmet', imwname = imwname_def, 
               imwdir = imwdir_def):

# Lasketaan maksimien/minimien jakaumat iMW-tasossa.
# Lasketaan kuukausimaksimien tai minimien jakaumat yhdelle mareografille.
# mar: yksi kirjain = mareografi
# Jos mm==0, lasketaan maksimien jakauma.
# Jos mm==1, lasketaan minimien jakauma.
# avuosi-lvuosi: havaintovuodet jotka otetaan mukaan
# aineisto: 'kk' => kk-ääriarvot, 'vrk' => vrk-ääriarvot, 
#       'hr1' => tuntiarvot, 'hr4' => 4 tunnin valiset
# (Vrk-, tunti- ja 4h-arvoja ei ole koodattu python-versioon.)
# n2000vuosi: vuosi jolle lasketaan N2000-tason jakauma (joka muuttuu vuosittain)
# jakaumassa_n: käytetään kumul. jakaumalle 1/n-mallia, ei 1/(n+1)
# (n = vanha tapa, n+1 = oikeampi otosjakauma eli tämä parametri on mukana
# taaksepäin yhteensopivuuden vuoksi)
# vanhasovitus: käytetään 2014 raportin sovitustapaa (sovitus mareografikohtaiseen
# määrään jakauman korkeimpia pisteitä), tämäkin taaksepäin yhteensopivuutta
# dir_out: tänne tulokset
# lahde: data-aineiston lähde: 'txt' = tekstikanta, 'smartmet' = Smartmet Server
# imwname: iMW (tarkennettu keskivesi) -tiedostojen nimirakenne
# imwdir: iMW-tiedostojen hakemisto
#
# Kaytetaan apufunktiota maxjak_uusi, joka laskee jakauman aineistosta.
# Lisaksi kaytetaan apufunktiota jak_ekstrap, joka sovittaa Weibull- tai
# eksponentiaalisen jakauman.
#
# Tallentaa tulokset ascii-tiedostoihin ja piirtää kuvat jakaumista

    # Haetaan keskivedet (N2000-tasossa):
    imw = hae_imw(mar, imwdir, imwname)

    #% Lasketaan havainnoista ylitystodennakoisyysjakauma erillisella
    #% funktiolla
    vkj, ccdf = maxjak_uusi(mar, mm = mm, avuosi = avuosi, lvuosi = lvuosi, 
                            ntaso = imw, aineisto = aineisto, 
                            jakaumassa_n = jakaumassa_n, lahde = lahde)

    #% Muunnetaan mm => cm
    vkj = vkj / 10.0

    if mm == 1:                
        #% (Minimit ekstrapoloidaan "kaantamalla" alitustod.nak.
        #% jakauma ylitystod.nak.jakaumaksi.)
        vkj = -vkj

    #% Kerroin jolla todennakoisyysjakauma kerrotaan jotta saadaan "events/year"
    if aineisto == 'kk':
        tkerr = 12
    elif aineisto == 'vrk':
        tkerr = 365.25
    elif aineisto == 'hr1':
        tkerr = 365.25*24
    elif aineisto == 'hr4':
        tkerr = 365.25*6

    smax = 500 # korkein vedenkorkeusarvo joka ekstrapoloidaan

    #% Maaritellaan kuinka moneen pisteeseen sovitus tehdaan
    if aineisto == 'kk':
        if vanhasovitus:
        #% kaytetaan Katrin sovitustapaa kk-maksimeille
            nsovmar = {'a': 17, 'o': 45, 'b': 40, 'p': 22, 'v': 45, 's': 40, \
                   'm': 30, 'r': 13, 't': 25, 'd': 13, 'h': 40, 'e': 22, 'f': 33}
            nsov = nsovmar[mar]
            iextvk, iextccdf, p1, p2, p3 = jak_ekstrap(vkj, ccdf, smax, n = nsov, tapa = 'exp')
        else:
            # kaytetaan sovitusta 1 events/yr alaspain
            iextvk, iextccdf, p1, p2, p3 = jak_ekstrap(vkj, ccdf, smax, 
                                            raja = [0, 1/tkerr], tapa = 'exp0')
    else:
        # sovitetaan todennakoisyystasosta 2 events/yr alaspain
        iextvk, iextccdf, p1, p2, p3 = jak_ekstrap(vkj, ccdf, smax, 
                                            raja = [0, 2/tkerr], tapa = 'exp0')

    #% Kaannetaan minimijakaumat takaisin
    if mm == 1:
        vkj = -vkj
        iextvk = -iextvk

    #% Muunnetaan todennakoisyys muotoon "events/year":
    ccdf *= tkerr
    iextccdf *= tkerr

    #% Tallennetaan jakauma ja ekstrapolaatio omiin tiedostoihinsa.

    ttunnus = mar + '_' + ('max' if mm == 0 else 'min')
    ttunn2 = str(avuosi) + '_' + str(lvuosi) + aineisto + \
            ('_oldn' if jakaumassa_n else '') \
            + ('_oldfit' if vanhasovitus else '') + ('_sm' if lahde == 'smartmet' else '')
            

    with open(dir_out + '/' + ttunnus + 'jak' + ttunn2 + '.dat', 'w') as fout:
        fout.write('%%Cm\tProb\n')
        for i in range(0, len(vkj)):
            fout.write('{:.4f}\t{:.8g}\n'.format(vkj[i], ccdf[i]))

    with open(dir_out + '/' + ttunnus + 'sov' + ttunn2 + '.dat', 'w') as fout:
        fout.write('%%Cm\tProb\n')
        for i in range(0, len(iextvk)):
            fout.write('{:.4f}\t{:.8g}\n'.format(iextvk[i], iextccdf[i]))


    #% Piirretaan ja tallennetaan kuva jakaumasta
    if n2000vuosi > 0:
        imwv = float((imw[imw[:, 0] == n2000vuosi, 1]) / 10.0)
    else:
        imwv = 0.0
            
    [fig, ax] = plt.subplots()
    ax.semilogy(vkj + imwv, ccdf, 'b+', iextvk + imwv, iextccdf, 'r-')
    ax.grid('on', which = 'both')
    if n2000vuosi == 0.0:
        ax.set_xlabel('Sea level (cm, iMW)')
        ax.set_title(mareo_names[mar])
    else:
        ax.set_xlabel('Sea level (cm, N2000)')
        ax.set_title(mareo_names[mar] + ', distribution in ' + str(n2000vuosi))

    ax.set_ylabel('Exceedance freq.')
    ax.set_ylim([1e-3, 1e1])
    ax.set_xlim([0, 350] if mm == 0 else [-300, 0])
    fig.savefig(dir_out + '/' + ttunnus + ttunn2 + '.png', dpi = 300)
    plt.close()


def maxjak_uusi(mar = 'e', mm = 0, avuosi = 1887, lvuosi = 2016, 
                ntaso = 0.0, aineisto = 'kk', jakaumassa_n = False,
                lahde = 'smartmet'):

# Lasketaan kuukausimaksimien tai minimien jakauma, kayttaen
# havaintoaineistosta mareografilta mar, vuodesta avuosi vuoteen lvuosi asti.
#
# mm: 0 => tarkastellaan maksimeita
#     1 => tarkastellaan minimeita
# ntaso: aikasarja [vuosi, vedenk], jota kaytetaan nollatasona, tai
#       nollatason vakioarvo kaikille vuosille (esim. N60)
# aineisto: 'kk' => kk-maksimit, 'vrk' => vrk-maksimit
#
# Palauttaa:
# vky: vedenkorkeus (cm)
# ccp: complementary cumulative probability (maksimit) tai cumulative probability (minimit)

    if lahde not in ('txt', 'smartmet'):
        return [], []

    #% Luetaan maksimit/minimit tiedostosta:

    #% Kuukausiaariarvot vastaavasta aariarvotiedostosta
    if aineisto == 'kk':
        # Tekstikanta
        if lahde == 'txt':
            yyr, ymax, ymin, kkyr, kkmon, kkmax, kkmin = lue_max(mar)
            yr = kkyr[(kkyr >= avuosi) & (kkyr <= lvuosi)]
            # Valitaan min/max ja muunnetaan korkeusjärjestelmä ref => N2000
            if mm == 1: # minimit
                vk = kkmin[(kkyr >= avuosi) & (kkyr <= lvuosi)] - N2000mar[mar]
            elif mm==0: # maksimit
                vk = kkmax[(kkyr >= avuosi) & (kkyr <= lvuosi)] - N2000mar[mar]
        elif lahde == 'smartmet':
            mondata = loadfile(smdir + \
                               mareo_names_ver2[mar].replace('ä','a').replace('ö','o') + \
                                   '_monthlystats_N2000.txt')
#            yrdata = loadfile(smdir + mareo_names[mar] + '_annualstats_N2000.txt')
            ind = (mondata[:, 0] >= avuosi) & (mondata[:, 0] <= lvuosi)
            vk = mondata[ind, 3 + mm]  # maksimit sarake 3, minimit sarake 4
            yr = mondata[ind, 0]

# Huom. vuorokausi- ja tuntiarvolaskentojen alkuperaiset Matlab-koodit alla,
# ei ole toteutettu Pythonille
#
#% Vuorokausiaariarvot tiedostosta joka on laskettu sealevel-ohjelmalla (vuorokausikeski- ja 
#% aariarvot seka hajonnat N2000-tasossa). Tiedoston otsikkorivien alkuun on lisatty %-merkit 
#% jotta load suostuu lukemaan sen.
#elseif strcmp(aineisto, 'vrk') == 1
#    data=load(['k:/vrkdata/',martunn(mar),'_vrk.txt']);
#    % valitaan haluttu aikajakso
#    data=data(data(:,1)>=avuosi&data(:,1)<=lvuosi,:);
#    if mm==1 % minimit: sarake 7
#        data=data(:,[1,7]);
#    elseif mm==0 % maksimit: sarake 6 (Huom! Jos tiedostossa ei ole hajontaa, sarakkeet muuttuvat!)
#        data=data(:,[1,6]);
#    end;
#
#% Tuntiarvot tai 4 tunnin valiset arvot. Tata ei ole kunnolla testattu!
#% Luetaan sealevel-ohjelmalla laskettu tiedosto jossa on tuntiarvojen aikasarja.
#elseif strcmp(aineisto, 'hr1') == 1 || strcmp(aineisto, 'hr4') == 1
#    data=load(['k:/vrkdata/',martunn(mar),'_tunti.txt']);
#    % valitaan haluttu aikajakso
#    data=data(data(:,1)>=avuosi&data(:,1)<=lvuosi,:);
#    % valitaan 4 tunnin valiset arvot klo 02, 06 jne.
#    if strcmp(aineisto, 'hr4') == 1 
#        data=data(ismember(data(:,4),[2,6,10,14,18,22]),:);
#    end;
#    data=data(:,[1,5]); % sarake 1 = vuosi, sarake 5 = tuntiarvot (muuta aikatietoa ei tarvita)
#end;    
#
#% Nyt muuttujassa data: sarake 1 = vuosi, sarakkeet 2-n = vedenkorkeusarvoja
#

# Jos nollataso on vakioarvo, vahennetaan se
    if isinstance(ntaso, float):
        vk = vk - ntaso

#% Jos nollataso on annettu aikasarjana, vahennetaan vuosittaiset arvot.
    else:
        for ii in range(0, len(yr)):
            iy = np.where(ntaso[:, 0] == yr[ii])[0]
            if len(iy) == 0:
                nt = float('nan')
            else:
                nt = ntaso[iy[0], 1]
            vk[ii] = vk[ii] - nt


    vk = vk[~np.isnan(vk)]

    if mm==1: #% kaannetaan minimit
        vk = -vk

    # Lasketaan jakauma
    vky, ccp, cp = ylijak(vk, samanraja = 0.001, jakaumassa_n = jakaumassa_n)

    #% Kaannetaan minimit takaisin.
    if mm == 1:
        vky = -vky

    return vky, ccp


def jak_ekstrap(vkj, ccdf, smax, raja = [], n = 0, tapa = 'exp0'):

# Sovitetaan todennakoisyysjakaumaan (ccdf) loppupaahan Weibull-jakauma.
# (jos length(raja)==2, sovitetaan todennakoisyysvalille raja(1) - raja(2))
# Jos n>0, tehdaan sovitus n korkeimpaan pisteeseen em. vedenkorkeusrajan
# sijasta.
# smax: korkein x-arvo mihin asti sovitetaan
# tapa: 'exp': eksponenttifunktio Weibull-funktiolla, 
#   'weib3': 3-parametrin Weibull,
#   'exp0': eksponenttifunktio suoralla sovituksella
# t, fparam: sovituksessa kaytettavat parametrit fparam ja t
# (vrt. weibfit2)
#
# Palauttaa:
# extvk, extccdf: ekstrapoloitu jakauma: x-arvot ja todennakoisyydet (ccdf)
# pweib, pcon, wextcon: sovituksen parametrit
# 

    if tapa not in ('exp', 'weib3', 'exp0'):
        print('jak_ekstrap: vain sovitukset exp0, exp ja weib3 toteutettu.')
        return [], [], [], [], []
    
    iok = ~np.isnan(vkj)
    vkj = vkj[iok]
    ccdf = ccdf[iok]

    #% Weibull-jakauman sovitus:

    #% Valitaan pisteet joihin sovitetaan:
    #% Ensimmainen vaihtoehto: sovitetaan n korkeimpaan pisteeseen.
    if n > 0:
        jsov = vkj
        csov = ccdf
        if len(jsov) > n:
            jsov = jsov[-n:]
            csov = csov[-n:]
    
    #% Toinen vaihtoehto: sovitetaan tiettyjen rajojen valille.
    else:
        if len(raja) == 2:
            raja2 = raja[1]
            raja = raja[0]
        else:
            raja2 = float('nan')
            
        if np.isnan(raja2):
            i1 = (ccdf > raja)
        elif np.isnan(raja):
            i1 = (ccdf < raja2)
        else:
            i1 = (ccdf > raja) & (ccdf < raja2)
        jsov = vkj[i1]
        csov = ccdf[i1]

    #% Ei soviteta jos pisteita oli liian vahan.
    if len(jsov) < 2:
        print('jak_ekstrap: jakaumaa ei voi sovittaa alle 2 pisteeseen.')
        return [], [], [], [], []

    smin = jsov[0]    
    svk = np.flipud(np.arange(smax, smin, -5))

    # Tehdaan 3 parametrin Weibull-sovitus (jos parametri tapa == 'weib3').
    if tapa == 'weib3':
        pweib = weibfit3(jsov, 1 - csov, [-1000,1000], 5e-4)
        pcon = []

    # Toinen vaihtoehto: tehdaan eksponentiaalinen sovitus 
    # eli Weibull jonka eksponenttiparametri on 1,
    # tulee kun t=1 ja fparam=1.
    elif tapa == 'exp':
        pweib, pp, p2, con2 = weibfit2(jsov, 1 - csov, t = 1, fparam = 1)
        pcon = [p2, con2]
    # Kolmas vaihtoehto: eksponentiaalinen sovitus suoraan
    # ilman Weibull-funktioiden käyttöä.
    # (Tulee sama lopputulos kuin exp -vaihtoehdolla)
    elif tapa == 'exp0':
        pweib = expfit2(jsov, 1 - csov)
        pcon = []

    # Lasketaan Weibull/eksponentiaalisen funktion arvot halutuille vedenkorkeuksille.
    extvk = svk
    if tapa == 'exp0':
        extccdf = 1 - cum_expdistr(pweib, svk)
    else:
        extccdf = 1 - cweibull(pweib[0:2], svk - pweib[2])
    if pcon != []:
        wextcon = [svk, np.exp(-p2[0] * svk - p2[1]), 
                    np.exp(-(p2[0] + con2[0]) * svk -(p2[1] - con2[1])),
                    np.exp(-(p2[0] - con2[0]) * svk - (p2[1] + con2[1]))]
    else:
        wextcon = []

    return extvk, extccdf, pweib, pcon, wextcon


        
def ylijak(jak, jtn = [], samanraja = 0, jakaumassa_n = False):

# Lasketaan ylitystodennakoisyysjakauma.
# Jakauma lasketaan lukumaarajakaumasta tai vedenkorkeuden havaintosarjasta.
# jak: vedenkorkeudet tms.
# jtn: havaintojen lukumaarat
# Jos jtn == [] => jak tulkitaan havaintosarjaksi (pelkkia vedenkorkeusarvoja
# sellaisenaan).
# samanraja: eroavaisuus, jota lahempana toisiaan olevat arvot
# yhdistetaan samaksi
# jakaumassa_n: laskee "vanhalla" (vuoteen 2016) laskentatavalla jossa 
# kaytettiin jakajaa N
# eli ccdf-jakauman tn-tasot olivat 1/N => N/N, kun uudessa tavassa 
# ne ovat 1/(N+1) => N/(N+1)

# Palauttaa:
#   vks: vedenkorkeusarvot (jakauman x-akseli)
#   ccps: complementary cumulative distribution (ccdf)
#   cps: kumulatiivinen jakauma (cdf), cps = 1 - ccps

    #% Muutetaan havaintosarja lukumaarajakaumaksi:
        # (Tama toimii nopeasti, jos samoja arvoja ei esiinny usein.)
    if len(jtn) == 0:
        vks = np.sort(jak)
        jtn = np.ones(len(jak))
        dvks = np.diff(vks)
        ind = dvks <= samanraja
        while any(ind):
            ind2 = np.concatenate(([False], ind))
            ind1 = np.concatenate((ind, [False]))
            jtn[ind2] = jtn[ind2] + 1
            jtn[ind1] = jtn[ind1] - 1
            vks = vks[jtn > 0]
            jtn = jtn[jtn > 0]
            dvks = np.diff(vks)
            ind = dvks <= samanraja
    else:
        vks = jak
    
    # Poistetaan nollahavainnot ja NaN:t:
    vks = vks[jtn > 0]
    jtn = jtn[jtn > 0]
    jtn = jtn[~np.isnan(vks)]
    vks = vks[~np.isnan(vks)]

    # vks on nyt ascending toisin kuin alkuperaisessa!!

    #% Lasketaan kumulatiivinen jakauma:
    #% (sellaisten havaintojen tod.nak., jotka ovat <=x)

    if jakaumassa_n:
        # Vanha laskentamenetelma ccdf:lle
        lkm = np.sum(jtn)
        ccps = np.flipud(np.cumsum(np.flipud(jtn)))
        ccps = ccps / lkm
        cps = 1.0 - ccps
    else:
        lkm = np.sum(jtn) + 1  #% Tahan lisatty N+1 22.12.2016!!
        cps = np.cumsum(jtn)
        cps = cps / lkm
        # Lasketaan compementary cum. jakauma
        ccps = 1.0 - cps

    return vks, ccps, cps



def poimitasot(tasot = [0.01], avuosi = 1982, lvuosi = 2011, hakem = 'k:/jakaumat/',
               minimit = False, n2000vuosi = 0, aineisto = 'kk', tdstoloppu = '',
               imwname = imwname_def, imwdir = imwdir_def, output_mw = False,
               dir_out = resultdir):

# Taulukko tehdään N2000-järjestelmässä.
# Koska keskivesi muuttuu N2000:n suhteen vuosittain, myös jakaumat muuttuvat.
# Siksi taulukko pätee tietylle vuodelle.
# n2000vuosi: vuosi jolle taulukko tehdään.
# output_mw: jos True, tehdään myös MW:ssa (huom. myös tämä muuttuu hieman vuosittain 
#           koska iMW ja MW eivät ole samat; vain jakauma iMW:ssä ei riipu vuodesta)
    
    tulos = {}
    if minimit:
        mstr = 'min'
    else:
        mstr = 'max'
    for mar in 'aobpvsmrtdhef':
        tulos[mar] = {}            
        msov = loadfile(hakem + mar + '_' + mstr + 'sov' + str(avuosi) + '_' + str(lvuosi) + \
                        aineisto + tdstoloppu + '.dat')
        if msov[0, 1] > msov[-1, 1]:
            msov = np.flipud(msov)
        if n2000vuosi > 0:
            imw = hae_imw(mar, imwdir, imwname)
            imwv = float(imw[imw[:, 0] == n2000vuosi, 1] / 10) # iMW N2000-järjestelmässä cm
            mw = read_mw2()
            mwv = (mw[n2000vuosi][mar] - N2000mar[mar]) / 10.0 # MW N2000-järjestelmässä cm
        else:
            imwv = 0.0
        
        for t in tasot:
            vk = np.interp(np.log10(t), np.log10(msov[:, 1]), msov[:, 0])
            tulos[mar][t] = vk + imwv

    luvut = {}
    for hs in (['N2000', 'MW'] if output_mw else ['N2000',]):
        prstr = 'Vuodelle ' + str(n2000vuosi) + \
            ', aineistosta ' + str(avuosi) + '-' + str(lvuosi) + ', ' + hs + '\n'
        prstr += 'Mar'
        for t in tasot:
            prstr += '\t' + str(t)
        prstr += '\n'
    
        luvut[hs] = {}
        for mar in 'aobpvsmrtdhef':
            luvut[hs][mar] = []
            prstr += '{:11s}'.format(mareo_names[mar])
            for t in tasot:
                vk = tulos[mar][t] - (mwv if hs == 'MW' else 0.0)
                prstr += '\t{:.0f}'.format(vk)
                luvut[hs][mar].append(vk)
            prstr += '\n'

        print(prstr)
        
        with open(dir_out + '/jakauma_vuodelle' + str(n2000vuosi) + '_aineisto' \
                       + str(avuosi) + '-' + str(lvuosi) + '_' + hs + '.txt', 'w') as f:
            f.write(prstr)

    return luvut

