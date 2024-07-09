from autodp.minimal_expansion import MinimalExpansion
from autodp.tree_of_bags import TreeOfBags 
from enum import Enum


class TreeDecomposition(TreeOfBags):

    def __init__(self):

        # initializing content dict, adj, dict
        super().__init__()

        self.graph = MinimalExpansion()

        # sub-tree contents: compute them once, and then use them a lot
        self.vertices_above = {} # dict: TD edge --> set of vertices above it.
        self.vertices_below = {} # dict: TD edge --> set of vertices below it

    def set_graph(self, min_exp):
        self.graph = min_exp

    def pick_root(self):
        m = min(self.graph.vertices)
        M = max(self.graph.vertices)

        print("content ",self.bag_content)
        for key, val in self.bag_content.items():
            print(m, M, val)
            if str(m) in val and str(M) in val:
                self.root = key
                break

    def fill_vertices_above(self):
        """
        going to define and use a recursive function.
        """

        def fill_above(t,u,v,acc):
            """
            recursive function, with accumulator, to compute the 
            set of vertices above a certain edge.
            the edge in question is u-->v (root to leaf direction)
            the set of vertices 'above' is the union of bag contents
            on the u-side when cutting the edge u,v in the tree.
            
            
            """
    
            # initialization
            self.vertices_above[(u,v)] = set(self.bag_content[u])
            self.vertices_above[(u,v)] = self.vertices_above[(u,v)].union(acc)

            # include the vertices below the cousins of u.
            for w in self.bag_adj[u]:
                if w!=v and w!=t:
                    self.vertices_above[(u,v)] = self.vertices_above[(u,v)].union(self.vertices_below[(u,w)])
        
            # recursive calls
            for w in self.bag_adj[v]:
                if w!=u:
                    fill_above(u,v,w,self.vertices_above[(u,v)])


        # actual use of the recursive function. '-1' is a dummy value that is never a bag label.
        self.vertices_above[('-1',self.root)] = set([])
        for u in self.bag_adj[self.root]:
            fill_above('-1',self.root, u, set([]))

    def fill_vertices_below(self):

        def fill_below(u,v):
            """
            recursive DP function to fill dictionary.
            """
            # if already computed just return
            if (u,v) in self.vertices_below.keys():
                return self.vertices_below[(u,v)]
    
            # start with the most obvious part.
            self.vertices_below[(u,v)] = set(self.bag_content[v])    

            # recursive calls to add what is below v.
            for w in self.bag_adj[v]:
                if w!=u:
                    self.vertices_below[(u,v)] = self.vertices_below[(u,v)].union(fill_below(v,w))
    
            return self.vertices_below[(u,v)]

        # main call.
        fill_below('-1', self.root)

    def replace(self,u,v,lo,hi,val,to_keep_safe):
        """
        (u,v) must be an edge of the tree
        decomposition.
        starting from u and in the direction
        of v, all indices/vertices
        contained in the interval [lo,hi]
        are replaced with val
        """
    
        queue = [(u,v)]
    
        while len(queue) > 0:
            u,v = queue.pop()
            assert(u in self.bag_adj[v])
            new_content = []
            for vertex in self.bag_content[v]:
                if vertex[0]=='H':
                    new_content.append(vertex)
                    continue
                if int(vertex) >= lo and int(vertex) <= hi and vertex not in to_keep_safe:
                    new_content.append(str(val))
                    print("replacing",vertex,"with",val,"in bag",v)
                else:
                    new_content.append(vertex)
    
            self.bag_content[v] = new_content
    
            for w in self.bag_adj[v]:
                if w!=u:
                    queue.append((v,w))

    def diag_canonicize(self,u,v,i,j,ip,jp,helixname):
        """
        in direction of u, replaces all helix occurences
        with i,j. same in direction of v with ip,jp
        puts intermediary bags in between for helix.
        """
        inter = set(self.bag_content[u]).intersection(set(self.bag_content[v]))

        to_keep_safe = inter.intersection(set([str(c) for c in [i,j,ip,jp]]))

        ## replacements 
        self.replace(v,u,i,ip,i,to_keep_safe)
        self.replace(v,u,jp,j,j,to_keep_safe)
        self.replace(u,v,i,ip,ip,to_keep_safe)
        self.replace(u,v,jp,j,jp,to_keep_safe)

        ## insertions of new bags
    
        # disconnecting
        self.bag_adj[u] = [b for b in self.bag_adj[u] if b != v]
        self.bag_adj[v] = [b for b in self.bag_adj[v] if b != u]
    
        # constant part, inter without helix, except extremities
        const = [vertex for vertex in inter if int(vertex) <= i or int(vertex) >= j or (int(vertex) >= ip and int(vertex) <= jp)]

        # let it be clear:
        self.bag_content[u] = list(const) + [str(i),str(j)]
        self.bag_content[v] = list(const) + [str(ip),str(jp)]

        # new intermediary bags
        cur_i = i
        cur_j = j
        cur_bag = u
        cnt = 0
    
        i_turn = True
    
        while cur_i < ip or cur_j > jp:
            prev_bag = cur_bag
            cur_bag = helixname +'_'+ str(cnt)
            if i_turn:
                self.bag_content[cur_bag] = list(const) + [str(cur_i), str(cur_j), str(cur_i+1)]
                self.bag_adj[prev_bag].append(cur_bag)
                self.bag_adj[cur_bag] = [prev_bag]
                cur_i += 1
                i_turn = False
    
            else:
            # turn of j to move
                self.bag_content[cur_bag] = list(const) + [str(c) for c in [cur_i, cur_j, cur_j-1]]
                self.bag_adj[prev_bag].append(cur_bag)
                self.bag_adj[cur_bag] = [prev_bag]
                cur_j -= 1
                i_turn = True
            cnt += 1
    
        # connecting to last bag (v)
        self.bag_adj[cur_bag].append(v)
        self.bag_adj[v].append(cur_bag)

    def add_clique_case_subtree(self,subtree_root, 
                                helixname,
                                i,ip,jp,j):
        """
        Args:
            - subtree_root: the index of the bag that was
            seen to contain the four extremities of the 
            helix. we will attach to it a bag containing
            just the four extremities of the helix,
            where will be rooted a sub-tree decomposition
            containing the helix
    
            - clique_case_index: The index for the
            bag containing only the 4 extremities. I.e.
            the separator. This index is not just a number:
            it contains the name of the helix as well.
    
            - abcd: the separator, as a list.
    
            `- adj: the adjacency of the tree. will be modified
            to include the helix.
    
            - index2bag: dictrionary from indices to bag content
    
            - graph_adj: graph adjacency dictionnary: used to detect
            when in bulges (vertices have degree=2 in them)
        """
        
        cur_i = i
        cur_j = j
        cnt = 0
        prev_bag = subtree_root
        cur_bag = helixname+'_'+str(cnt)
        iturn = True
        const = [str(ip),str(jp)]
    
        while cur_i < ip or cur_j > jp:
            if iturn:
                # build
                bag_content = [str(cur_i),str(cur_j),str(cur_i+1)] + const
                self.bag_content[cur_bag] = bag_content
                self.bag_adj[prev_bag].append(cur_bag)
                self.bag_adj[cur_bag] = [prev_bag]
    
                # update
                cur_i += 1
                prev_bag = cur_bag
                cnt += 1
                cur_bag = helixname+'_'+str(cnt)
                iturn = False
            else:
                #jturn, build
                bag_content = [str(cur_i),str(cur_j),str(cur_j-1)] + const
                self.bag_content[cur_bag] = bag_content
                self.bag_adj[prev_bag].append(cur_bag)
                self.bag_adj[cur_bag] = [prev_bag]
    
                # update
                cur_j -= 1
                prev_bag = cur_bag
                cnt += 1
                cur_bag = helixname+'_'+str(cnt)
                iturn = True

    def take_clique_minor(self, i, ip, jp, j, mid_point, ksupm):

        v = self.bag_adj[self.root][0]

        if ksupm:
            # to make substitution in whole tree: pick an edge and 
            # go in both directions from it

            # one direction
            self.replace(self.root,v,i,i+mid_point,i,[])
            self.replace(self.root,v,i+mid_point+1,ip,ip,[])
            self.replace(self.root,v,jp,j-mid_point,jp,[])
            self.replace(self.root,v,j-mid_point+1,j,j,[])
            # the other
            self.replace(v,self.root,i,i+mid_point,i,[])
            self.replace(v,self.root,i+mid_point+1,ip,ip,[])
            self.replace(v,self.root,jp,j-mid_point,jp,[])
            self.replace(v,self.root,j-mid_point+1,j,j,[])
        else: 
            self.replace(self.root,v,i,i+mid_point-1,i,[])
            self.replace(self.root,v,i+mid_point,ip,ip,[])
            self.replace(self.root,v,jp,j-mid_point-1,jp,[])
            self.replace(self.root,v,j-mid_point,j,j,[])
            self.replace(v,self.root,i,i+mid_point-1,i,[])
            self.replace(v,self.root,i+mid_point,ip,ip,[])
            self.replace(v,self.root,jp,j-mid_point-1,jp,[])
            self.replace(v,self.root,j-mid_point,j,j,[])

    def contract_identical_bags(self):

        def substitute(vertex, v, u):
            if vertex==v:
                return u
            return vertex

        smth_contracted = True
        
        while smth_contracted:
            print('------------')
            smth_contracted = False
            queue = [('1',ngbh) for ngbh in self.bag_adj['1']]
        
            while len(queue) > 0:
                u, v = queue.pop()
                if set(self.bag_content[u])==set(self.bag_content[v]) or set(self.bag_content[v]).issubset(set(self.bag_content[u])) or set(self.bag_content[u]).issubset(self.bag_content[v]):
                    if (u[0]!='H') and (v[0]=='H'):
                    # no hybrid contractions for clarity
                        for w in self.bag_adj[v]:
                            if w!=u:
                                queue.append((v,w))
                        continue
        
                    print("must contract !",v,"into",u,"(",self.bag_content[u],"<-",self.bag_content[v],")")
                    smth_contracted = True
                    self.bag_adj[u] += self.bag_adj[v]
                    self.bag_adj[u].remove(u)
                    self.bag_adj[u].remove(v)
                    if set(self.bag_content[u]).issubset(set(self.bag_content[v])):
                        self.bag_content[u] = self.bag_content[v]
                    self.bag_content.pop(v)
                    for w in self.bag_adj[v]:
                        self.bag_adj[w] = [substitute(vertex, v,u) for vertex in self.bag_adj[w]]
                        if w!=u:
                            queue.append((u,w))
                    self.bag_adj.pop(v)
                    continue
        
                for w in self.bag_adj[v]:
                    if w!=u:
                        queue.append((v,w))


    def assert_vertex_representation(self, i):
        for u in self.dfs_bag_iterator():
            if i in self.bag_content[u]:
                return True
        return False

    def assert_vertex_connectivity(self, i):
    
        def recursive_connectivity_assert(prev, u, seen_above, currently_active):
            if i in self.bag_content[u]:
                if seen_above:
                    if not currently_active:
                        return False
                else:
                    seen_above = True
                    currently_active = True
            else:
                currently_active = False

            ans = True

            for v in self.bag_adj[u]:
                if v!=prev:
                    ans = ans and recursive_connectivity_assert(u,v, seen_above, currently_active)

            return ans
        
        return recursive_connectivity_assert('-1',self.root, False, False)

    def assert_edge_representation(self, i,j):

        queue = [('-1', self.root)]
        for u in self.dfs_bag_iterator():
            if i in self.bag_content[u] and j in self.bag_content[u]:
                return True

        return False
