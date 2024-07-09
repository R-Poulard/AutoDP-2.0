def add_element_to_adj(adj, vertex, new_neighbor):
    """
    simple try catch, just annoying to write all the time
    """
    try:
        adj[vertex].append(new_neighbor)
    except KeyError:
        adj[vertex] = [new_neighbor]

    return adj

def read_gr_file(lines):
    """
    Reads lines of a .gr file and 
    returns an adjacency graph
    """
    graph_adj = {}
    for line in lines[1:]:
        u = line.split(' ')[0]
        v = line.split(' ')[1].rstrip('\n')

        try:
            graph_adj[u].append(v)
        except KeyError:
            graph_adj[u] = [v]

        try:
            graph_adj[v].append(u)
        except KeyError:
            graph_adj[v] = [u]

    return graph_adj

def read_td_lines(lines):

    # extracting bags and tree from td file
    index2bag = {}
    adj = {}

    started_bs = False

    for line in lines:
        if line[0]=='b':
            started_bs = True
            index2bag[line.split(' ')[1].rstrip('\n')] = [vertex.rstrip('\n') for vertex in line.split(' ')[2:-1]]

        else:
            if started_bs:
                i = line.split(' ')[0]
                j = line.split(' ')[1].rstrip('\n')
                try:
                    adj[i].append(j)
                except KeyError:
                    adj[i] = [j]
                try:
                    adj[j].append(i)
                except KeyError:
                    adj[j] = [i]

    return adj, index2bag

### TREE DEC PROCESSING FUNCTIONS


def full_compact_diag(label, extremities, adj, index2bag, extremities_label, root):

    def replace_extremities(vert, extremities, indices_to_keep):
        print("extremities",extremities)
        if vert not in extremities or vert in indices_to_keep:
            return vert
        else:

            subs = {extremities[0]:extremities[2],
                    extremities[1]:extremities[3],
                    extremities[2]:extremities[0],
                    extremities[3]:extremities[1]}

            return subs[vert]

    below_helix = False
    queue = [('-1',root)]
    while len(queue) > 0:
        prev,u = queue.pop()

#        if below_helix:
#            index2bag[u] = list(set([replace_extremities(vert, extremities, indices_to_keep) for vert in index2bag[u]])) 
#            if u[0]=='H':
#                cluster_extremities[u.split('_')[0]] = list([replace_extremities(vert, extremities, indices_to_keep) for vert in cluster_extremities[u.split('_')[0]]])

        if u.split('_')[0]==label and prev.split('_')[0]!=label:
            if not set(extremities).issubset(set(index2bag[prev])):
            # diag case
                below_helix = True
                indices_to_keep = set(extremities).intersection(set(index2bag[u]))
                
                for e in set(extremities) - indices_to_keep:
                    extremities_label[e] = extremities_label[replace_extremities(e, extremities, indices_to_keep)]

                keep_going = True
                while keep_going:
                    keep_going = False
                    for v in adj[u]:
                        if v.split('_')[0]==label:
                            for w in adj[v]:
                                if w!=u:
                                    print("diag contracting ", index2bag[v], "into", index2bag[u])
                                    adj[u] = [bag for bag in adj[u] if bag!=v]+[w]
                                    adj[w] = [bag for bag in adj[w] if bag!=v]+[u]
                                    keep_going = True

        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))

    return adj, index2bag, extremities_label


def full_compact_clique(label, extremities, adj, index2bag,root):
    queue = [('-1',root)]
    while len(queue) > 0:
        prev,u = queue.pop()
        if u.split('_')[0]==label:
            if set(extremities).issubset(set(index2bag[prev])):
            # clique case
                index2bag[u] = extremities
                adj[u] = []
        
        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))

    return adj, index2bag

def filter_extremities(all_extremities, adj, index2bag, root):
    queue = [('-1',root)]
    print(index2bag.keys())
    while len(queue) > 0:
        prev,u = queue.pop()

        index2bag[u] = [vert for vert in index2bag[u] if vert in all_extremities]

        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))

    return index2bag

def subgraph_info_extraction(label,  adj, root):
    
    subgraph_string = ""
    queue = [('-1',root)]
    content = []
    while len(queue) > 0:
        prev,u = queue.pop()

        # final filtering: extremities only
        print(label, u[:2], u) 
        if u.split('_')[0]==label:
            content.append(u)
            if len(subgraph_string) > 0:
                subgraph_string += ' -> '
                subgraph_string += u
            else:
                root = prev
                subgraph_string += u

        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))

    subgraph_string += ';'

    return subgraph_string, content, root

def equations_prep_work(tdname, helix, root):
    # extracting bags and tree from td file
    adj, index2bag = read_td_lines(open(tdname).readlines())

    # extracting extremities
    cluster_extremities = {}
    for helixline in open(helix).readlines():
        label = helixline.split(' ')[0]
        extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]
        cluster_extremities[label] = extremities

    all_extremities = set([])

    for key, val in cluster_extremities.items():
        all_extremities = all_extremities.union(set(val))

    # extremities label dictionary
    extremities_label = {}

    for k, e in enumerate(sorted(list(all_extremities), key=lambda x:int(x))):
        extremities_label[e] = chr(ord('a')+k) 

    # for all helices, contraction to one bag (the highest in the tree)
    for helixline in open(helix).readlines():
        label = helixline.split(' ')[0]
        extremities = cluster_extremities[label]

        adj, index2bag, extremities_label = full_compact_diag(label, extremities, adj, index2bag, extremities_label, root)
        adj, index2bag = full_compact_clique(label, extremities, adj, index2bag, root)

    index2bag = filter_extremities(all_extremities, adj, index2bag, root)
    # final filtering: extremities only

    # for all helices: building the colored box around it
    subgraph = {}
    cluster_content = {}
    cluster_root = {}

    # which cluster filling
    which_cluster = {}
    for u in index2bag.keys():
        if u[0]=='H':
            which_cluster[u] = u.split('_')[0]

    # subgraph constructions
    for helixline in open(helix).readlines():

        label = helixline.split(' ')[0]
        subgraph_string, content, croot = subgraph_info_extraction(label, adj, root)

        if subgraph_string!=';':
            subgraph[label] = subgraph_string
        cluster_content[label] = content
        cluster_root[label] = croot


    # in the following working with bags of letters
    index2letters = {}
    for key, val in index2bag.items():
        index2letters[key] = []
        for vert in val:
            if vert in all_extremities:
                index2letters[key].append(extremities_label[vert])
            else:
                index2letters[key].append(vert)

    return adj, all_extremities, cluster_root, cluster_content, cluster_extremities, extremities_label, index2bag, index2letters, subgraph, which_cluster
