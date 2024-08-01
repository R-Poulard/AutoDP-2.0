from utils import read_gr_file, add_element_to_adj
from autodp.minimal_expansion import MinimalExpansion
from autodp.tree_decomposition import TreeDecomposition
import sys

sys.setrecursionlimit(10000)

# retrieve graph
inter_helix_gap=snakemake.config["inter_helix_gap"]
underlying_graph = MinimalExpansion()
dbn = open(snakemake.input.dbn).readlines()[0].rstrip('\n')
underlying_graph.from_str(dbn, inter_helix_gap=inter_helix_gap)

print(underlying_graph.helices)

# tree decomposition object init
tree_dec = TreeDecomposition()
tree_dec.set_graph(underlying_graph)

tree_dec.read_from_file(snakemake.input.tdname)

# preprocessing for faster computations later on.
tree_dec.fill_vertices_below()
tree_dec.fill_vertices_above()

vert_below = tree_dec.vertices_below
vert_above = tree_dec.vertices_above

width = max([len(val) for key, val in tree_dec.bag_content.items()])-1

# helix processing
for hline in open(snakemake.input.helix).readlines():

    processed = False

    extremities = hline.split('(')[1].split(')')[0].split(',')
    helixname = hline.split(' ')[0]
    print("helix ", hline)

    i = int(extremities[0]) 
    ip = int(extremities[1])
    jp = int(extremities[2])
    j = int(extremities[3])

    if width <= 3:
    # only diag case

        for k in range(ip-i):
            if processed:
                break
            for m in range(k+2,ip+1-i,1):
                if processed:
                    break
                # looking for a separator separating (i+k,j-k) and (i+m,j-m)

                # iterating over all edges of the tree decomposition
                # in a dfs way
                print(i+k,j-k,"sep",i+m, j-m," ?")
                for u,v in tree_dec.dfs_edge_iterator():
                    if u[0]=='H' or v[0]=='H':
                        continue 

                    above = vert_above[(u,v)]
                    below = vert_below[(u,v)]

                    if str(i+k) in above-below and str(j-k) in above-below and str(i+m) in below-above and str(j-m) in below-above:
                        print("-->yes does sep !")
                        inter = set(tree_dec.bag_content[u]).intersection(set(tree_dec.bag_content[v]))
#                        
                        tree_dec.diag_canonicize(u,v,i,j,ip,jp,helixname)
                        processed = True
                        break
                    if str(i+k) in below-above and str(j-k) in below-above and str(i+m) in above-below and str(j-m) in above-below:
                        print("-->yes does sep !")

                        inter = set(tree_dec.bag_content[u]).intersection(set(tree_dec.bag_content[v]))

                        tree_dec.diag_canonicize(v,u,i,j,ip,jp,helixname)
                        processed = True
                        break

        # no need to check other cases. will do same for them, and check completeness at the end
        if processed:
            continue

    else:
    # else width >=4

        # detecting clique case = detecting hop edge
        print("i,ip,jp,j", i, ip, jp, j)
        found_hop = False
        for k in range(ip-i+1):
            if found_hop:
                break
            for m in range(ip-i+1):
                if found_hop:
                    break
                if abs(m-k) > 1:
                    # go over bags to see if one represents edge (i+k, j-m):
                    print("looking for edge ", i+k,j-m)
                    for u in tree_dec.dfs_bag_iterator():
                        if str(i+k) in tree_dec.bag_content[u] and str(j-m) in tree_dec.bag_content[u]:
                            found_hop = True
                            print("found in ", u, tree_dec.bag_content[u])
                            # point to build G_clique minor
                            if k < m:
                                mid_point = k+1
                                ksupm = False
                            else:
                                mid_point = m+1
                                ksupm = True
                            break
        
        # then clique case
        if found_hop:
            print("found clique case with mid point", mid_point, ksupm)
            for u in tree_dec.dfs_bag_iterator():
                print(u, tree_dec.bag_content[u])
            tree_dec.take_clique_minor(i,ip,jp,j,mid_point,ksupm)
            # find bag with i,ip,jp,j
            print("after minor")
            for u in tree_dec.dfs_bag_iterator():
                print(u, tree_dec.bag_content[u])

            for u in tree_dec.dfs_bag_iterator():
                print([str(vert) for vert in [i,ip,jp,j]], set(tree_dec.bag_content[u]))
                if set([str(vert) for vert in [i,ip,jp,j]]).issubset(set(tree_dec.bag_content[u])):
                    tree_dec.add_clique_case_subtree(u, helixname, i,ip,jp,j)
                    print("setting processed to true")
                    processed = True
                    break
            
        else:
            print("did not find clique case, resorting to diag case.")
        # diag case. have to find ij/ipjp separator
            
            for u, v in tree_dec.dfs_edge_iterator():

                if u[0]=='H' or v[0]=='H':
                    continue 

                print(tree_dec.vertices_above.keys())
                above = tree_dec.vertices_above[(u,v)]
                below = tree_dec.vertices_below[(u,v)]

                if str(i) in above-below and str(j) in above-below and str(ip) in below-above and str(jp) in below-above:
                    tree_dec.diag_canonicize(u,v,i,j,ip,jp,helixname)
                    processed = True
                    break
                if str(i) in below-above and str(j) in below-above and str(ip) in above-below and str(jp) in above-below:
                    tree_dec.diag_canonicize(v,u,i,j,ip,jp,helixname)
                    processed = True
                    break

    # assert that helix has been processed, i.e. assert completeness of disjunction
    try:
        assert(processed) 
    except AssertionError:
        raise AssertionError

# removing adjacent redundancies
tree_dec.contract_identical_bags()

####################################
####################################
####################################


# writing output (finally)
f = open(snakemake.output[0],'w')

# header line
tree_dec.pick_root()
print("root_picked", tree_dec.root)
print('root', tree_dec.root, file=f)
print('s td '+str(len(tree_dec.bag_content.keys()))+' ', end="", file=f)
# largest bag size
width = max([len(val) for _, val in tree_dec.bag_content.items()])
print(str(width)+' NVERTICES', file=f)


for key, val in tree_dec.bag_content.items():
    print('b',key,end="",file=f)
    for vertex in sorted(list(set(val))):
        print('',vertex, end="",file=f)
    print("", file=f)

for u,v in tree_dec.dfs_edge_iterator():
    print(u,v,file=f)


####################################
####################################
####################################

#### CHECKS ON NEWLY COMPUTED TREE DEC ####
print(tree_dec.bag_content)
max_vertex = max([max([int(i) for i in val]) for key ,val in tree_dec.bag_content.items() if key!='-1'])
tree_dec.assert_edge_representation('1', str(max_vertex))


for u in underlying_graph.vertices:
    # represented
    try:
        assert(tree_dec.assert_vertex_representation(str(u)))
    except AssertionError:
        print(u, " not represened")
        raise AssertionError
    # in a connected way
    try:
        assert(tree_dec.assert_vertex_connectivity(str(u)))
    except AssertionError:
        print("subtree", u, " not connected")
        raise AssertionError

for u,v in underlying_graph.edges:
    try:
        assert(tree_dec.assert_edge_representation(str(u),str(v)))
    except AssertionError:
        print(u,v," is not represented")
        raise AssertionError
