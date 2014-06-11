from collections import defaultdict

def d_from_dir_adj_list(a):
    n = len(a)
    d = defaultdict(list)
    for i, v in enumerate(a):
        for j in v:
            d[i].append(j)
            d[j].append(i)
    return dict(d)

def d_from_links(links):
    d = defaultdict(list)
    for k1, k2 in links:
        d[k1].append(k2)
        d[k2].append(k1)
    return dict(d)

def d_from_m_bip(m):
    d = defaultdict(list)
    n1 = len(m)
    for i in range(n1):
        for j in range(len(m[0])):
            if m[i][j]:
                d[i].append(j+n1)
                d[j+n1].append(i)
    return dict(d)

def is_regular(d, r):
    """
    check that ``d`` is regular with degree ``r``
    """
    for k, v in d.items():
        if len(v) != r:
            return False
        for k2 in v:
            if k not in d[k2]:
                raise ValueError('%s is not in the adjacency list of %s' %(k, k2))
    return True

def complete_d(d):
    """
    transform a directed graph in the corresponding undirected graph
    """
    a = d.items()
    for k, v in a:
        for i in v:
            if k not in d[i]:
                d[i].append(k)
    for k, v in d.items():
        v.sort()


def sq_mat(n1, n2, typ):
    """
    reduced adjacency matrix of a rectangle (n1, n2)

    n1 : number of vertices on side ``(x, 0)``
    n2 : number of vertices on side ``(0, y)``
    typ : type of boundary condition: ``pp`` periodic in both
    directions; ``np`` open boundary conditions
    """
    d1 = {}
    d2 = {}
    c1 = 0
    c2 = 0
    v1 = []
    v2 = []
    if typ == 'pp':
        for i in range(n1):
            for j in range(n2):
                if (i + j ) % 2 == 0:
                    v1.append((i, j))
                    d1[(i,j)] = c1
                    c1 += 1
                else:
                    v2.append((i, j))
                    d2[(i,j)] = c2
                    c2 += 1
        m = [[0]*len(d2) for i in range(len(d1))]
        for i in range(n1):
            for j in range(n2):
                if (i + j)%2 == 0:
                    r1 = d1[(i, j)]
                    m[r1][d2[(i, (j + 1) % n2)]] = 1
                    m[r1][d2[(i, (j - 1) % n2)]] = 1
                    m[r1][d2[((i + 1) % n1, j)]] = 1
                    m[r1][d2[((i - 1) % n1, j)]] = 1
        return m, d1, d2
    elif typ == 'np':
        for i in range(n1):
            for j in range(n2):
                if (i + j )%2 == 0:
                    v1.append((i, j))
                    d1[(i,j)] = c1
                    c1 += 1
                else:
                    v2.append((i, j))
                    d2[(i,j)] = c2
                    c2 += 1
        m = [[0]*len(d2) for i in range(len(d1))]
        for i in range(n1):
            for j in range(n2):
                if (i + j)%2 == 0:
                    r1 = d1[(i, j)]
                    if j < n2 - 1:
                        m[r1][d2[(i, (j + 1))]] = 1
                    if j >= 1:
                        m[r1][d2[(i, j - 1)]] = 1
                    if i < n1 - 1:
                        m[r1][d2[(i + 1, j)]] = 1
                    if i >= 1:
                        m[r1][d2[(i - 1, j)]] = 1
        return m, d1, d2
    else:
        raise NotImplementedError

def sq_d_np(n1, n2):
    def _append(d, r1, r2):
        if r2 not in d[r1]:
            d[r1].append(r2)
        if r1 not in d[r2]:
            d[r2].append(r1)

    dc = {}
    d = defaultdict(list)
    c = 0
    v = []
    for i in range(n1):
        for j in range(n2):
            v.append((i, j))
            dc[(i, j)] = c
            c += 1
    for i in range(n1):
        for j in range(n2):
            r1 = dc[(i, j)]
            if j < n2 - 1:
                r2 = dc[(i, j + 1)]
                _append(d, r1, r2)
            if j >= 1:
                r2 = dc[(i, j - 1)]
                _append(d, r1, r2)
            if i < n1 - 1:
                r2 = dc[(i + 1, j)]
                _append(d, r1, r2)
            if i >= 1:
                r2 = dc[(i - 1, j)]
                _append(d, r1, r2)
    d1 = {}
    for k, v in d.items():
        v.sort()
        d1[k] = v
    return d1

def hexagon_p(n1, n2):
    """
    hexagon lattice in brick wall representation, with periodic b.c.
    """
    assert n1 % 2 == 0
    assert n2 %2 == 0
    m, d1, d2 = sq_mat(n1, n2, 'pp')
    for i in range(n1):
        for j in range(n2):
            if (i + j ) % 2 == 1:
                r2 = d2[(i, j % n2)]
                r1 = d1[(i, (j + 1) % n2)]
                m[r1][r2] = 0
    return m, d1, d2

def triangle_lattice_pp(nx, ny):
    """
    triangle lattice with open boundary conditions
    return m, d, a
    m adjacency matrix
    d dict ((i, j), node_index)
    a  list of links
    """
    c = 0
    d = {}
    for i in range(nx):
        for j in range(ny):
            d[(i, j)] = c
            c += 1
    n = nx*ny
    a = []
    m = [[0]*n for i in range(n)]
    for i in range(nx):
        for j in range(ny):
            r1 = d[(i, j)]
            r2 = d[(i, (j+1)%ny)]
            m[r1][r2] = 1
            m[r2][r1] = 1
            a.append((r1, r2))
            r2 = d[((i-1)%nx, (j+1)%ny)]
            m[r1][r2] = 1
            m[r2][r1] = 1
            a.append((r1, r2))
            r2 = d[((i+1)%nx, j)]
            m[r1][r2] = 1
            m[r2][r1] = 1
            a.append((r1, r2))
    return m, d, a

def kings_sq_mat_np(nx, ny):
    """
    kings on a rectangle with open b.c.
    return m, d, a
    m adjacency matrix
    d dict ((i, j), node_index)
    a  list of links
    """
    c = 0
    d = {}
    for i in range(nx):
        for j in range(ny):
            d[(i, j)] = c
            c += 1
    n = nx*ny
    a = []
    m = [[0]*n for i in range(n)]
    for i in range(nx):
        for j in range(ny):
            r1 = d[(i, j)]
            if j < ny - 1:
                r2 = d[(i, j+1)]
                m[r1][r2] = 1
                m[r2][r1] = 1
                if i < nx - 1:
                    r3 = d[(i+1, j)]
                    r4 = d[(i+1, j+1)]
                    a.append((r1, r2, r3, r4))
            if i < nx - 1:
                r2 = d[(i+1, j)]
                m[r1][r2] = 1
                m[r2][r1] = 1
    return m, d, a

def hard_sq_np(nx, ny):
    """
    hard squares on a rectangle with open b.c.
    return m, d, a
    m adjacency matrix
    d dict ((i, j), node_index)
    a  list of links
    """
    d = {}
    c = 0
    for i in range(2*nx + 1):
        for j in range(2*ny + 1):
            if (i + j) % 2 == 1:
                d[(i, j)] = c
                c += 1
    a = []
    n = len(d)
    m= [[0]*n for i in range(n)]
    for i in range(1, 2*nx):
        for j in range(2*ny - 1):
            if (i + j) % 2 == 0:
                continue
            r1 = d[(i, j)]
            r2 = d[(i + 1, j + 1)]
            r3 = d[(i - 1, j + 1)]
            r4 = d[(i, j + 2)]
            if i % 2 == 1:
                a.append((r1,r2,r3,r4))
            m[r1][r2] = 1
            m[r2][r1] = 1
            m[r1][r3] = 1
            m[r3][r1] = 1
            m[r2][r4] = 1
            m[r4][r2] = 1
            m[r3][r4] = 1
            m[r4][r3] = 1
    return m, d, a

def wanless_example1(k, f):
    c1 = c2 = 0
    d1 = {}
    d2 = {}
    d1['s'] = 0
    c1 += 1
    d2['t'] = 0
    c2 += 1
    for a in range(1, k + 1):
        for b in range(1, k + 1):
            for c in range(1, f + 1):
                d1[(a,b,c)] = c1
                c1 += 1
                d2[(a,b,c)] = c2
                c2 += 1
    m = [[0]*len(d2) for i in range(len(d1))]
    for a in range(1, k + 1):
        for b in range(1, k + 1):
            for b1 in range(1, k + 1):
                for c in range(1, f + 1):
                    if b == 1 and b1 == 1:
                        continue
                    r1 = d1[(a, b, c)]
                    r2 = d2[(a, b1, c)]
                    m[r1][r2] = 1
    for a in range(1, k + 1):
        for c in range(1, f):
            r1 = d1[(a, 1, c)]
            r2 = d2[(a, 1, c+1)]
            m[r1][r2] = 1
    for a in range(1, k + 1):
        r1 = d1['s']
        r2 = d2[(a, 1, 1)]
        m[r1][r2] = 1
    for a in range(1, k + 1):
        r1 = d1[(a, 1, f)]
        r2 = d2['t']
        m[r1][r2] = 1
    return m

def sc_mat_np(n1, n2, n3):
    """
    return reduced adjacency matrix and vertices for slab
    """
    d1 = {}
    d2 = {}
    c1 = 0
    c2 = 0
    # v1 vertices in first set, v2 in second set of bipartite lattice
    v1 = []
    v2 = []
    for i1 in range(n1):
        for i2 in range(n2):
            for i3 in range(n3):
                if (i1 + i2 + i3 )%2 == 0:
                    v1.append((i1, i2, i3))
                    d1[(i1,i2,i3)] = c1
                    c1 += 1
                else:
                    v2.append((i1, i2, i3))
                    d2[(i1,i2,i3)] = c2
                    c2 += 1
    m = [[0]*len(d2) for i in range(len(d1))]
    #print 'm=', m
    for i1 in range(n1):
        for i2 in range(n2):
            for i3 in range(n3):
                if (i1 + i2 + i3)%2 == 0:
                    r1 = d1[(i1, i2, i3)]
                    if i1 < n1 - 1:
                         m[r1][d2[(i1 + 1, i2, i3)]] = 1
                    if i1 >= 1:
                        m[r1][d2[(i1 - 1, i2, i3)]] = 1
                    if i2 < n2 - 1:
                        m[r1][d2[(i1, i2 + 1, i3)]] = 1
                    if i2 >= 1:
                        m[r1][d2[(i1, i2 - 1, i3)]] = 1
                    if i3 < n3 - 1:
                        m[r1][d2[(i1, i2, i3 + 1)]] = 1
                    if i3 >= 1:
                        m[r1][d2[(i1, i2, i3 - 1)]] = 1
    return m, v1, v2, d1, d2


# fullerenes
def dict_fuller(n):
    """
    some fullerene isomers
    ``24,26,28,30,34,36`` used in [1]
    ``40`` used in d40nano_get_links
    ``60`` buckminster fullerene
    ``70`` used in [2]

    References
    ==========

    [1] M.B. Ahmadi, H. Alimorad Daskhezr,
         {\it Commun. Math. Comput. Chem.} {\bf 71}(2014) 355.
    [2] D. Babic and O. Ori, Chem. Phys. Lett. 234, 240 (1995).
    """
    if n == 20:
        d = {0:[1,2,3], 1:[0,4,19], 2:[0,7,18], 3:[0,5,6], 4:[1,5,15],
    5:[3,4,8], 6:[3,7,9], 7:[2,6,10], 8:[5,9,12], 9:[6,8,11], 10:[7,11,17],
    11:[9,10,13], 12:[8,13,15], 13:[11,12,14], 14:[13,16,17],
    15:[4,12,16], 16:[14,15,19],17:[10,14,18],18:[2,17,19], 19:[1,16,18]}

    elif n == 24:
        d = {0:[1,2,3],1:[0,7,22],2:[0,5,6],3:[0,4,21],4:[3,5,8],5:[2,4,10],
    6:[2,7,11],7:[1,6,13],8:[4,9,14],9:[8,10,16],10:[5,9,11],
    11:[6,10,12],12:[11,13,17],13:[7,12,19],14:[8,15,21],15:[14,16,20],
    16:[9,15,17],17:[12,16,18],18:[17,19,20],19:[13,18,22],
    20:[15,18,23],21:[3,14,23],22:[1,19,23],23:[20,21,22]}

    elif n == 26:
        d = {0:[1,2,3],1:[0,4,24], 2:[0,12,25],3:[0,5,6],4:[1,5,7],
        5:[3,4,9],6:[3,11,12],
    7:[4,8,13],8:[7,9,15],9:[5,8,10],10:[9,11,16],11:[6,10,19],
    12:[2,6,20],13:[7,14,21],14:[13,15,17],15:[8,14,16],16:[10,15,18],
    17:[14,18,22],18:[16,17,19],19:[11,18,20],20:[12,19,23],21:[13,22,24],
    22:[17,21,23], 23:[20,22,25],24:[1,21,25],25:[2,23,24]}

    elif n == 28:
        d = {0:[1,2,3],1:[0,19,26],2:[0,5,6],3:[0,4,27],4:[3,5,12],
    5:[2,4,7],6:[2,8,9],7:[5,8,10],8:[6,7,11],9:[6,15,19],10:[7,13,14],
    11:[8,14,15],12:[4,13,20],13:[10,12,16],14:[10,11,17],15:[9,11,18],
    16:[13,17,21],17:[14,16,18],18:[15,17,22],19:[1,9,24],20:[12,23,27],
    21:[16,22,23],22:[18,21,24],23:[20,21,25],24:[19,22,25],25:[23,24,26],
    26:[1,25,27],27:[3,20,26]}

    elif n == 30:
        d = {0:[1,2,3],1:[0,5,6],2:[0,4,29],3:[7,28],4:[5,14],5:[8],
    6:[7,10],7:[25],8:[9,11],9:[10,13],10:[17],11:[12,15],12:[13,19],
    13:[16],14:[15,24],15:[18],16:[17,20],17:[21],18:[19,22],
    19:[20],20:[23],21:[23,25],22:[23,26],23:[],24:[26,29],
    25:[27],26:[27],27:[28],28:[29],29:[]}
        complete_d(d)

    elif n == 34:
        d = {0:[1,2,8],1:[4,5],2:[3,32],3:[4,16],4:[9],5:[6,10],6:[7,15],
    7:[8,28],8:[31],9:[10,11],10:[12],11:[13,17],12:[14,15],
    13:[14,18],14:[20],15:[21],16:[17,26],17:[22],18:[19,23],
    19:[20,24],20:[21],21:[25],22:[23,27],23:[33],24:[25,33],
    25:[28],26:[27,32],27:[29],28:[30],29:[30,33],30:[31],31:[32],
    32:[],33:[]}
        complete_d(d)

    elif n == 36:
        d = {0:[1,2,3],1:[4,32],2:[5,6],3:[7,34],4:[5,12],5:[8],6:[7,10],
        7:[17],8:[9,11],9:[10,14],10:[35],11:[12,13],12:[18],
        13:[14,20],14:[15],15:[22,35],16:[17,23,35],17:[31],
        18:[19,28],19:[20,24],20:[21],21:[22,25],22:[23],23:[26],
        24:[25,27],25:[26],26:[30],27:[28,29],28:[32],29:[30,33],
        30:[31],31:[34],32:[33],33:[34],34:[],35:[]}
        complete_d(d)

    elif n == 40:
        d = {0:[1,2,3], 1:[0,14,17], 2:[0,5,8], 3:[0,4,6], 4:[3,18,21],
    5:[2,10,13], 6:[3,7,8], 7:[6,22,25], 8:[2,6,9], 9:[8,26,29],
    10:[5,11,29], 11:[10,12,34], 12:[11,13,37], 13:[5,12,14],14:[1,13,15],
    15:[14,16,38],16:[15,17,41],17:[1,16,18], 18:[4,17,19], 19:[18,20,42],
    20:[19,21,45], 21:[4,20,22], 22:[7,21,23], 23:[22,24,46],
    24:[23,25,49], 25:[7,24,26],26:[9,25,27], 27:[26,28,30],
    28:[27,29,33], 29:[9,10,28], 30:[27,31,49], 31:[30,32,37],
    32:[19,30,31], 33:[32,35,38], 34:[30,31,37], 35:[28,33,37],
    36:[24,27,38], 37:[11,34,35], 38:[33,36,39], 39:[20,23,38]}


    elif n == 60:
        d = {0:[1,2,3], 1:[0,4,58], 2:[0,8,59], 3:[0,5,6], 4:[1,9,32], 5:[3,9,10],
    6:[3,7,11], 7:[6,8,14], 8:[2,7,16], 9:[4,5,17], 10:[5,11,19], 11:[6,10,12],
    12:[11,13,21], 13:[12,14,23], 14:[7,13,15], 15:[14,16,25], 16:[8,15,44],
    17:[9,18,33], 18:[17,19,26], 19:[10,18,20], 20:[19,21,27], 21:[12,20,22],
    22:[21,23,29],23:[13,22,24],24:[23,25,31], 25:[15,24,42], 26:[18,27,35],
    27:[20,26,28],28:[27,29,37], 29:[22,28,30], 30:[29,31,38], 31:[24,30,41],
    32:[4,33,48], 33:[17,32,34], 34:[33,35,45], 35:[26,34,36], 36:[35,37,46],
    37:[28,36,38], 38:[30,37,39], 39:[38,40,47], 40:[39,41,52], 41:[31,40,42],
    42:[25,41,43], 43:[42,44,53], 44:[16,43,57], 45:[34,46,49], 46:[36,45,47],
    47:[39,46,51], 48:[32,49,54], 49:[45,48,50], 50:[49,51,55], 51:[47,50,52],
    52:[40,51,53], 53:[43,52,56], 54:[48,55,58], 55:[50,54,56], 56:[53,55,57],
    57:[44,56,59], 58:[1,54,59], 59:[2,57,58]}

    elif n == 70:
        a70 = [[1,12,19],[2,3],[4,6],[5,11],[5,7],[10],[13,14],[8,15],[9,22], [10,24],
    [16],[17,18],[13,67],[69],[15,27],[20],[17,26],[35],[19,42],[68], [21,28],
    [22,29],[23],[24,31],[25],[26,33],[34],[36,69],[36,37],[30,37],[31,38],
    [32],[33,39],[40],
    [40,41],[41,42],[48],[43],[39,44],[45],[46],[47],[56],[44,49],[50],[46,51],
    [52],[52,55],[49,58],[53],[51,53],[54],[54],[59],[60],[56,62],[65],
    [58,66,69],[61],[60,61],[62],[63],[64],[64,66],[65],[68],[67],[68],[],[]]
        d = d_from_dir_adj_list(a70)

    else:
        raise ValueError('case not existing')

    assert is_regular(d, 3)
    return d

def d40nano_get_links(n, typ):
    m = 20
    if typ == 1:
        # links for north cap
        links = [(0,2),(2,8),(6,8),(3,6),(0,3),(0,1),(3,4),(6,7),(8,9),(2,5),
                (4,18),(17,18),(1,17),(1,14),(16,17),(18,19),(4,21),
                (5,13),(12,13),(13,14),(5,10),
                (20,21),(19,20),(21,22),(7,22),(7,25),
                (22,23),(23,24),(24,25),(25,26),(9,26),(9,29),(26,27),(10,29),
                (28,29),(27,28),
                (10,11),(11,12),(14,15),(15,16)]
        # links for the strips
        for i in range(n):
            a = [(11 + m*i, 14 + m*(i+1)), (12 + m*i, 17 + m*(i+1)),
        (15 + m*i, 18 + m*(i+1)), (17 + m*(i+1), 18 + m*(i+1)),
        (16 + m*(i+1), 17 + m*(i+1)),
        (19 + m*i, 22 + m*(i+1)), (16 + m*i, 21 + m*(i+1)),
        (21 + m*(i+1), 22 + m*(i+1)), (20 + m*(i+1), 21 + m*(i+1)),
        (18 + m*(i+1), 19 + m*(i+1)), (19 + m*(i+1), 20 + m*(i+1)),
        (22 + m*(i+1), 23 + m*(i+1)), (20 + m*i, 25 + m*(i+1)),
        (23 + m*i, 26 + m*(i+1)), (25 + m*(i+1), 26 + m*(i+1)),
        (24 + m*(i+1), 25 + m*(i+1)), (23 + m*(i+1), 24 + m*(i+1)),
        (26 + m*(i+1), 27 + m*(i+1)),
        (24 + m*i, 29 + m*(i+1)),
        (27 + m*i, 10 + m*(i+1)), (10 + m*(i+1), 29 + m*(i+1)),
        (28 + m*(i+1), 29 + m*(i+1)), (27 + m*(i+1), 28 + m*(i+1)),
        (28 + m*i, 13 + m*(i+1)), (13 + m*(i+1), 14 + m*(i+1)),
        (10 + m*(i+1), 11 + m*(i+1)), (11 + m*(i+1), 12 + m*(i+1)),
        (12 + m*(i+1), 13 + m*(i+1)), (14 + m*(i+1), 15 + m*(i+1)),
        (15 + m*(i+1), 16 + m*(i+1))]
            links.extend(a)
        N = m*n
        a = [(11,37),(28,35),(35,37),(33,35),(34,37),(27,36),(24,36),
            (36,38),(33,38),(38,39),(23,39),(20,39),(32,33),(19,32),
            (31,32),(16,31),(31,34), (30,34), (12,30), (15,30)]
        a = [(i + N, j + N) for i, j in a]
        links.extend(a)

    elif typ == 2:
        # links for north cap
        links = [(0,2),(2,8),(8,9),(6,9),(3,6),(0,3),(0,1),(3,4),(6,7),(2,5),
                (4,18),(17,18),(1,17),(1,14),(16,17),(18,19),(4,21),
                (5,13),(12,13),(13,14),(5,10),
                (20,21),(19,20),(21,22),(7,22),(7,25),
                (22,23),(23,24),(24,25),(25,26),(9,26),(8,29),(26,27),(10,29),
                (28,29),(27,28),
                (10,11),(11,12),(14,15),(15,16)]
        # links for the strips
        for i in range(n):
            a = [(11 + m*i, 14 + m*(i+1)), (12 + m*i, 17 + m*(i+1)),
        (15 + m*i, 18 + m*(i+1)), (17 + m*(i+1), 18 + m*(i+1)),
        (16 + m*(i+1), 17 + m*(i+1)),
        (19 + m*i, 22 + m*(i+1)), (16 + m*i, 21 + m*(i+1)),
        (21 + m*(i+1), 22 + m*(i+1)), (20 + m*(i+1), 21 + m*(i+1)),
        (18 + m*(i+1), 19 + m*(i+1)), (19 + m*(i+1), 20 + m*(i+1)),
        (22 + m*(i+1), 23 + m*(i+1)), (20 + m*i, 25 + m*(i+1)),
        (23 + m*i, 26 + m*(i+1)), (25 + m*(i+1), 26 + m*(i+1)),
        (24 + m*(i+1), 25 + m*(i+1)), (23 + m*(i+1), 24 + m*(i+1)),
        (26 + m*(i+1), 27 + m*(i+1)),
        (24 + m*i, 29 + m*(i+1)),
        (27 + m*i, 10 + m*(i+1)), (10 + m*(i+1), 29 + m*(i+1)),
        (28 + m*(i+1), 29 + m*(i+1)), (27 + m*(i+1), 28 + m*(i+1)),
        (28 + m*i, 13 + m*(i+1)), (13 + m*(i+1), 14 + m*(i+1)),
        (10 + m*(i+1), 11 + m*(i+1)), (11 + m*(i+1), 12 + m*(i+1)),
        (12 + m*(i+1), 13 + m*(i+1)), (14 + m*(i+1), 15 + m*(i+1)),
        (15 + m*(i+1), 16 + m*(i+1))]
            links.extend(a)
        N = m*n
        a = [(11,37),(28,35),(35,37),(33,35),(34,37),(27,36),(24,36),
            (36,38),(33,38),(38,39),(23,39),(20,39),(32,33),(19,32),
            (31,32),(16,31),(31,34), (30,34), (12,30), (15,30)]
        a = [(i + N, j + N) for i, j in a]
        links.extend(a)
    else:
        raise ValueError('case not existing')

    return links

def line_graph(d):
    """
    return the dict for the line graph of the graph with dict ``d``
    """
    d1 = {}
    edges = []
    c = 0
    dn = {}
    for k1, v in d.items():
        for k2 in v:
            if k2 < k1:
                continue
            t = (k1, k2)
            edges.append(t)
            dn[t] = c
            c += 1
    v = [0]*c
    for k, cx in dn.items():
        v[cx] = k
    for i in range(c):
        w = []
        a = v[i]
        for j in range(c):
            if j == i:
                continue
            b = v[j]
            if a[0] in b or a[1] in b:
                w.append(j)
        d1[i] = w
    return d1
