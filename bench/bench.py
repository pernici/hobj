import sys
sys.path.insert(0,'../src')
sys.path.insert(0,'../tests')

from time import time
from hobj import (hobj_str, dup_permanental_minor_poly,
  dup_matching_generating_poly, dup_independence_poly)
from active_nodes import ordered_links
from graphs_gen import sq_mat, dict_fuller, hexagon_p, d_from_m_bip
from test_hobj import band_mat1
from domains import ZZ, QQ

def band_mat(n, k, modulus=11, den=1):
    m = [[0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            if abs(i - j) < k:
                if den == 1:
                    m[i][j] = i*j % modulus + 1
                else:
                    m[i][j] = QQ(i*j, den) + 1
    return m

def test_perm_minor1():
    m = band_mat(12, 12, modulus=7)
    r = dup_permanental_minor_poly(m, ZZ)
    assert r == [176737059423072, 907717809454560, 1016936320591440, 438391944900848, 91329708507784, 10371759217336, 688884839796, 27859313948, 697946924, 10771496, 98743, 489, 1]

def test_perm_minor2():
    m = band_mat(10, 10, den=7)
    #print 'm=', m
    r = dup_permanental_minor_poly(m, QQ)
    #print 'r=', r
    assert r == [QQ(200560062627334080,823543), QQ(856091105852040000,823543), 
        QQ(776434012513248960,823543), QQ(37461268420473600,117649), 
        QQ(5865790024153800,117649), QQ(1384177356648,343), 
        QQ(8737202304,49), QQ(214403400,49), QQ(2860425,49), QQ(2725,7), QQ(1,1)]

def test_matching_generating_poly1():
    d = dict_fuller(70)
    p = dup_matching_generating_poly(d)
    assert p == [52168, 24949770, 2175411410, 76188054720, 1414031839415, 16060741380337, 121977585926655, 657705883833585, 2628626806075280, 8042439217209190, 19311918160511451, 37116989376321905, 58004771574641390, 74650239997165790, 79943103597641110, 71843167484323524, 54553421377655995, 35194282142522435, 19372757557305170, 9127603794188980, 3688877411254177, 1280236743404535, 381580564799125, 97578536447295, 21362502290720, 3989784929340, 632435470975, 84479881155, 9417493125, 864583005, 64182939, 3753855, 166390, 5250, 105, 1]

def test_matching_generating_poly2():
    m, d1, d2 = hexagon_p(6, 8)
    d = d_from_m_bip(m)
    #k1, k2, links, max_active = ordered_links_all(d)
    #print 'DB41 links=', links
    r = dup_matching_generating_poly(d,K=ZZ, val=-1)
    assert r == 3296

def test_independent_poly1():
    d = dict_fuller(70)
    p = dup_independence_poly(d)
    assert p == [58250, 2335845, 41926970, 454409710, 3359594518, 18085177930, 73760007630, 234019294980, 588519920150, 1189774688888, 1955030446640, 2634407161260, 2932125022400, 2711236103195, 2092155247990, 1351710396475, 732697251110, 333467332370, 127354984490, 40728770134, 10865146960, 2403369030, 437000270, 64516615, 7599972, 696535, 47810, 2310, 70, 1]


if __name__ == '__main__':
    for f in [test_perm_minor1, test_perm_minor2,
      test_matching_generating_poly1, test_matching_generating_poly2,
      test_independent_poly1]:
        t0 = time()
        f()
        t1 = time()
        print('%-20s %.2f' %(f.__name__, t1 - t0))
