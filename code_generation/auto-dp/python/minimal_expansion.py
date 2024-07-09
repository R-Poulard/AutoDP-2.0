
parenthesis_systems = ['()','[]','{}','<>','Aa','Bb','Cc']

class MinimalExpansion:
    """
    A class simply meant to contain all info related to
    the minimal expansion of a fatgraph.
    """
    def __init__(self):
        self.dbn = "" # underlying dbn, one pair of parenthesis per helix.
        self.vertices = set([]) # 10 vertices in each helix (5-5)
        self.edges = set([]) # backbone + base pairs
        self.adj = {} # same but different representation
        self.helices = set([]) # helix coordinates.

    def connect(self, i, j):

        # store in ascending order
        if i > j:
            i, j = j, i

        # edge list (set) representation
        self.edges.add((i,j))

        # adj dict represenation
        try:
            self.adj[i].append(j)
        except KeyError:
            self.adj[i] = [j]

        try:
            self.adj[j].append(i)
        except KeyError:
            self.adj[j] = [i]


    def extract_helices(self, inter_helix_gap=True):
        """
        Assumes a single parenthesis pair per
        helix. 
        """
        # linking position in dbn to helix position in minimal expansion:
        str2helix_pos = {}
        anchor = 1
        for k in range(len(self.dbn)):
            str2helix_pos[k] = anchor
            anchor += 4
            # where inter helix gap matters
            if inter_helix_gap:
                anchor += 1

        # initializing stacks
        stacks = {}
        for s in parenthesis_systems:
            stacks[s] = []

        for k, c in enumerate(self.dbn):
            for s in parenthesis_systems:
                # opening parenthesis: add to stack
                if c==s[0]:
                    stacks[s].append(str2helix_pos[k])
                # closing: take corresponding opening and make helix
                if c==s[1]:
                    helix_start = stacks[s].pop()
                    helix_end = str2helix_pos[k]
                    self.helices.add((helix_start,helix_start+4,helix_end,helix_end+4))

    def from_str(self, dbn_str, inter_helix_gap=True, overarching=True):

        # setting dbn
        self.dbn = dbn_str

        print("dbn", dbn_str)

        # precomputation: helices
        self.extract_helices(inter_helix_gap=inter_helix_gap)

        print("helices", self.helices)

        for i,ip,jp,j in self.helices:
            # first half
            for vertex in range(i,ip+1,1):
                self.vertices.add(vertex)

            # second half
            for vertex in range(jp, j+1, 1):
                self.vertices.add(vertex)

            # base pairs
            for k in range(5):
                self.connect(i+k,j-k)

        # backbone
        lwb = min([min(helix) for helix in self.helices])
        upb = max([max(helix) for helix in self.helices])
    
        for pos in range(lwb, upb, 1):
            self.connect(pos, pos+1)

        # overarching edge
        if overarching:
            self.connect(lwb, upb)

    def dump_to_gr(self, gr_fname):
        f = open(gr_fname, 'w')

        f.write('p tw '+str(len(self.vertices))+' '+str(len(self.edges))+'\n')
        for u,v in self.edges:
            print(u,v,file=f)
