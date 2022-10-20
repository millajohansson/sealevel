# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 14:38:41 2016

@author: johanssm
"""

# Various information related to the mareographs.

# Mareograph character identifiers as a list
martunn = ('a', 'o', 'b', 'p', 'v', 's', 'm', 'r', 't', 'd', 'h', 
           'e', 'f', 'c')

# Mareograph names
mareo_names = {'a': 'Kemi', 
               'o': 'Oulu', 
               'b': 'Raahe', 
               'p': 'Pietarsaari', 
               'v': 'Vaasa', 
               's': 'Kaskinen',
               'm': 'Mantyluoto', 
               'r': 'Rauma', 
               't': 'Turku', 
               'd': 'Degerby', 
               'h': 'Hanko', 
               'e': 'Helsinki',
               'f': 'Hamina', 
               'c': 'Porvoo',
               1: 'Kemi', 
               2: 'Oulu', 
               3: 'Raahe', 
               4: 'Pietarsaari', 
               5: 'Vaasa', 
               6: 'Kaskinen',
               7: 'Mantyluoto', 
               8: 'Rauma', 
               9: 'Turku', 
               10: 'Degerby', 
               11: 'Hanko', 
               12: 'Helsinki',
               13: 'Hamina', 
               14: 'Porvoo'}

# Another version with Föglö and Mäntyluoto
mareo_names_ver2 = {'a': 'Kemi', 
               'o': 'Oulu', 
               'b': 'Raahe', 
               'p': 'Pietarsaari', 
               'v': 'Vaasa', 
               's': 'Kaskinen',
               'm': 'Mäntyluoto', 
               'r': 'Rauma', 
               't': 'Turku', 
               'd': 'Föglö', 
               'h': 'Hanko', 
               'e': 'Helsinki',
               'f': 'Hamina', 
               'c': 'Porvoo',
               1: 'Kemi', 
               2: 'Oulu', 
               3: 'Raahe', 
               4: 'Pietarsaari', 
               5: 'Vaasa', 
               6: 'Kaskinen',
               7: 'Mäntyluoto', 
               8: 'Rauma', 
               9: 'Turku', 
               10: 'Föglö', 
               11: 'Hanko', 
               12: 'Helsinki',
               13: 'Hamina', 
               14: 'Porvoo'}

# Swedish names
mareo_names_sv = {'a': 'Kemi', 
               'o': 'Uleåborg', 
               'b': 'Brahestad', 
               'p': 'Jakobstad', 
               'v': 'Vasa', 
               's': 'Kaskö',
               'm': 'Mäntyluoto', 
               'r': 'Raumo', 
               't': 'Åbo', 
               'd': 'Föglö', 
               'h': 'Hangö', 
               'e': 'Helsingfors',
               'f': 'Fredrikshamn', 
               'c': 'Borgå'}

# Mareografien järjestysnumerot datatiedostoissa ym.
mareo_num = {'a': 1,
          'o': 2,
          'b': 3,
          'p': 4,
          'v': 5,
          's': 6,
          'm': 7,
          'r': 8,
          't': 9,
          'd': 10,
          'h': 11,
          'e': 12,
          'f': 13,
          'c': 14}

# Toiminnan aloitusvuodet
avuodet=(1922, 1922, 1922, 1922, 1922, 1926, 1925, 
         1933, 1922, 1923, 1887, 1904, 1928, 2014)

avuodetmar = {'a': 1922, 'o': 1922, 'b': 1922, 'p': 1922, 'v': 1922, 's': 1926, 
              'm': 1925, 'r': 1933, 't': 1922, 'd': 1923, 'h': 1887, 'e': 1904, 
              'f': 1928, 'c': 2014}

# Datahakemistot Prut:lla
marhake=('Kemi', 'Oulu', 'Raahe', 'Pietarsa', 'Vaasa', 'Kaskinen', 
         'Mantyluo', 'Rauma', 'Turku', 'Degerby', 'Hanko', 
         'Helsinki', 'Hamina', 'Porvoo')

# Mareografien keskinaiset etaisyydet
eta = (82, 61, 141, 87, 82, 89, 52, 87, 110, 147, 116, 130)
et_th = 82  # Turku - Hanko
et_rd = 136 # Rauma - Degerby

# Mareografien koordinaatit
msij = {'a': [65.6733, 24.5153], 
   'o': [65.0402,   25.4183],
   'b': [64.6663,   24.4072],
   'p': [63.7085,   22.6895],
   'v': [63.0817,   21.5712],
   's': [62.3440,   21.2147],
   'm': [61.5943,   21.4667],
   'r': [61.1335,   21.4258],
   't': [60.4282,   22.1005],
   'd': [60.0320,   20.3847],
   'h': [59.8228,   22.9767],
   'e': [60.1533,   24.9597],
   'f': [60.5628,   27.1792]}


# N60-taso referenssitasossa
N60ref = (1675,1687,1667,1642,1667,1695,1715,1727,1804,1795,1861,1894,1886)
N60mar = {'a': 1675, 'o': 1687, 'b': 1667, 'p': 1642, 'v': 1667, 's': 1695, \
            'm': 1715, 'r': 1727, 't': 1804, 'd': 1795, 'h': 1861, 'e': 1894, \
            'c': 1893, 'f': 1886}

# N2000 referenssitasossa
N2000ref = (1264,1288,1242,1202,1233,1270,1331,1386,1512,1522,1609,1642,1674)
N2000mar = {'a': 1264, 'o': 1288, 'b': 1242, 'p': 1202, 'v': 1233, 's': 1270, \
            'm': 1331, 'r': 1386, 't': 1512, 'd': 1522, 'h': 1609, 'e': 1642, \
            'c': 1651, 'f': 1674}

# N43 ja NN referenssitasossa
N43mar = {'a': 1817, 'o': 1815, 'b': 1777, 'p': 1787, 'v': 1794, 's': 1824, \
            'm': 1851, 'r': 1865, 't': 1887, 'd': 1879, 'h': 1910, 'e': 1944, \
            'c': 1938, 'f': 1949}

NNmar = {'a': 2047, 'o': 2034, 'b': 2036, 'p': 2037, 'v': 2031, 's': 2022, \
            'm': 2011, 'r': 2000, 't': 1970, 'd': 1970, 'h': 1937, 'e': 1944, \
            'c': 1928, 'f': 1940}


# Convert FMISIDs to character identifiers
fmisid_to_chr = {134254: 'f', 
           100669: 'c', 
           132310: 'e', 
           134253: 'h', 
           134252: 'd',
           134225: 't', 
           134224: 'r', 
           134266: 'm', 
           134251: 's', 
           134223: 'v',
           134250: 'p', 
           100540: 'b', 
           134248: 'o', 
           100539: 'a'}
