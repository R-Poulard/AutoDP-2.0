
class TreeOfBags:

    def __init__(self):

        # empty
        self.bag_content = {}
        self.bag_adj = {}

        # silly first values:
        self.root = '-1'
        self.bag_content['-1'] = []

    def dfs_edge_iterator(self, no_minus_one=False):

        queue = [('-1', self.root)]
    
        while len(queue) > 0:
            u,v = queue.pop()
    
            # dummy edge choice
            if not no_minus_one or u!='-1':
                # yield
                yield u,v
        
            for w in self.bag_adj[v]:
                if w!=u:
                    queue.append((v,w))

    def dfs_bag_iterator(self):
        
        queue = [('-1', self.root)]
    
        while len(queue) > 0:
            u,v = queue.pop()
    
            # yield
            yield v
        
            for w in self.bag_adj[v]:
                if w!=u:
                    queue.append((v,w))
    
    def read_from_file(self, td_file, default_root='1'):

        # setting root:
        first_line = open(td_file).readlines()[0]
        if first_line.split(' ')[0]=='root':
            self.root = first_line.split(' ')[1].rstrip('\n')
        else:
            self.root = default_root
        
        # extracting bags and tree from td file
        started_bs = False

        for line in open(td_file).readlines():
            if line[0]=='b':
                started_bs = True
                self.bag_content[line.split(' ')[1].rstrip('\n')] = [vertex.rstrip('\n') for vertex in line.split(' ')[2:]]

            else:
                if started_bs:
                    i = line.split(' ')[0]
                    j = line.split(' ')[1].rstrip('\n')
                    try:
                        self.bag_adj[i].append(j)
                    except KeyError:
                        self.bag_adj[i] = [j]
                    try:
                        self.bag_adj[j].append(i)
                    except KeyError:
                        self.bag_adj[j] = [i]
