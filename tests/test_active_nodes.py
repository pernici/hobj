import sys
sys.path.insert(0,'../src')
from active_nodes import (ordered_links, num_active_nodes, ip_ordered_vertices,
    ip_list_objects_from_vlist, ip_get_dn, ip_num_active_elements)

from domains import ZZ
from graphs_gen import dict_fuller

def test_ordered_links():
    d = dict_fuller(60)
    links = ordered_links(d, 0, 1)
    nand = num_active_nodes(d, links)
    assert nand == 10

def test_ip_num_active_elements1():
    d = dict_fuller(60)
    dn = ip_get_dn(d)
    links0 = [0, d[0][0]]
    vlist = ip_ordered_vertices(d, *links0)
    objects = ip_list_objects_from_vlist(d, vlist)
    nu = ip_num_active_elements(objects)
    assert nu == 11


def test_ip_num_active_elements2():
    d = dict_fuller(36)
    dn = ip_get_dn(d)
    links0 = [0, d[0][0]]
    vlist = ip_ordered_vertices(d, *links0)
    objects = ip_list_objects_from_vlist(d, vlist)
    nu = ip_num_active_elements(objects)
    assert nu == 8


if __name__ == '__main__':
    test_ordered_links()
    test_ip_num_active_elements1()
    test_ip_num_active_elements2()
    print('test_active_nodes ok')
