from utils import read_td_lines


class DPTable:

    def __init__(self, name, indices):
        self.name = name
        self.indices = indices
        self.child_tables = []
         
    def add_child_table(table, substitution_matrix={}):
        self.child_tables.append((table, substitution_matrix))
    

def compute_dp_table(td_lines, helix_lines, f):


    # extracting bags and tree from td file
    root = td_lines[0].split(' ')[1].rstrip('\n')
    adj, index2bag = read_td_lines(td_lines)

    # extracting helix extemities information 
    cluster_extremities = {}
    for helixline in helix_lines:
        label = helixline.split(' ')[0]
        extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]
        cluster_extremities[label] = extremities

    # putting all extremities in one set (practical)
    all_extremities = set([])
    for key, val in cluster_extremities.items():
        all_extremities = all_extremities.union(set(val))

    # contraction: 
    for helixline in helix_lines:
        label = helixline.split(' ')[0]
        extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]

        queue = [('-1',root)]
        while len(queue) > 0:
            prev,u = queue.pop()

            if u.split('_')[0]==label:
                if not set(extremities).issubset(set(index2bag[prev])):
                # diag case
                    keep_going = True
                    while keep_going:
                        keep_going = False
                        for v in adj[u]:
                            if v.split('_')[0]==label:
                                for w in adj[v]:
                                    if w!=u and w.split('_')[0] == label:
                                        # grand-child is still in helix, contracting v into u
                                        adj[u] = [bag for bag in adj[u] if bag!=v]+[w]
                                        adj[w] = [bag for bag in adj[w] if bag!=v]+[u]
                                        keep_going = True

            for v in adj[u]:
                if v!=prev:
                    queue.append((u,v))

    cluster_content = {}
    which_cluster = {}
    for helixline in helix_lines:
        label = helixline.split(' ')[0]
       
        cluster_content[label] = []
        extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]

        queue = [('-1',root)]
        while len(queue) > 0:
            prev,u = queue.pop()
            print(label, u[:2], u) 
            if u.split('_')[0]==label:
                cluster_content[label].append(u)
                which_cluster[u] = label

                if set(extremities).issubset(set(index2bag[prev])):
                # clique case
                    index2bag[u] = extremities
                    adj[u] = []
                else:
                # diag case
                    index2bag[u] = [vert for vert in index2bag[u] if vert in all_extremities]

            for v in adj[u]:
                if v!=prev:
                    queue.append((u,v))

    const_part = {}

    for key, val in cluster_content.items():
        if len(val) > 0:
            const_part[key] = set.intersection(*[set(index2bag[u]) for u in val])
        else:
            const_part[key] = set([])

    # adj, index2bag, const_part, which_cluster, cluster_root, cluster_extremities, all_extremities have been
    # constructed

    cnt = 0

    queue = [('-1',root)]

    index2bag['-1'] = []

    def num_to_letters(cnt):
        if cnt < 26:
            return chr(ord('a') + cnt).upper() #'\\colorbox{c'+str(cnt)+'}{'+chr(ord('a') + cnt).upper()+'}'
        else:
            return chr(ord('a') + int(cnt/26)).upper()+chr(ord('a') + int(cnt%26)).upper()

    cnt = 0
    cnt_col = 0
    bag_letter = {}

    queue = [('-1',root)]

    while len(queue) > 0:
        prev,u = queue.pop()
        if not u[:1]=='H':
            bag_letter[u] = num_to_letters(cnt)
            cnt += 1
    #        cnt_col += 1
        else:
            if set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
            # clique
                bag_letter[u] = "\\colorbox{c"+u.split('_')[0][1:]+"}{$C_{\\boxtimes}$}"
                cnt_col += 1
            else:
            # diag
                if prev.split('_')[0]!=u.split('_')[0]:
                    bag_letter[u] = "\\colorbox{c"+u.split('_')[0][1:]+"}{$"+num_to_letters(cnt)+"$}"
                    cnt += 1
                    cnt_col += 1

        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))

    ext_to_letter = {}
    for k, e in enumerate(sorted(list(all_extremities),key=lambda x: int(x))):
        ext_to_letter[e] = chr(ord('a')+k)

    def partner(e, sorted_ext):
        if e==sorted_ext[0]:
            return sorted_ext[3]
        if e==sorted_ext[1]:
            return sorted_ext[2]
        if e==sorted_ext[3]:
            return sorted_ext[0]
        if e==sorted_ext[2]:
            return sorted_ext[1]

    def subs(e, sorted_exs):
        if e==sorted_exs[0]:
            return sorted_exs[1]
        if e==sorted_exs[1]:
            return sorted_exs[0]
        if e==sorted_exs[2]:
            return sorted_exs[3]
        if e==sorted_exs[3]:
            return sorted_exs[2]

    def increment(e, sorted_exs):
        if e==sorted_exs[0]:
            return '+1'
        if e==sorted_exs[1]:
            return '-1'
        if e==sorted_exs[2]:
            return '+1'
        if e==sorted_exs[3]:
            return '-1'

    queue = [('-1',root)]
    while len(queue) > 0:
        prev,u = queue.pop()
        if u[:1]=='H' and not set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
        # diag case
            if u.split('_')[0]!=prev.split('_')[0]:

                local_index_label = {}

                absent_ex = (set(cluster_extremities[which_cluster[u]])-set(index2bag[prev])).pop() 
                sorted_exs = sorted(cluster_extremities[which_cluster[u]], key=lambda x: int(x))
                
                indices = sorted(list(set(cluster_extremities[which_cluster[u]]) - set([absent_ex, partner(absent_ex, sorted_exs)])),key=lambda x: int(x))
                increments = [increment(e, sorted_exs) for e in indices]
                letters = [ext_to_letter[e] for e in indices]

                for i in indices:
                    local_index_label[i] = ext_to_letter[i]
                    local_index_label[subs(i, sorted_exs)] = ext_to_letter[i]

                const = []
                for e in sorted(const_part[which_cluster[u]], key=lambda x : int(x)):
                    if e in indices:
                        const.append(ext_to_letter[e]+"'")
                        local_index_label[e] = ext_to_letter[e]+"'"
                    else:
                        const.append(ext_to_letter[e])
                        local_index_label[e] = ext_to_letter[e]

                second_helix_bag = [v for v in adj[u] if v!=prev].pop()
                if len(adj[second_helix_bag]) >= 2:
                    normal_childs = []
                    normal_child_letters = {}
                    diag_childs = []
                    diag_child_letters = {}
                    diag_child_const_parts = {}
                    clique_childs = []
                    clique_childs_letters = {}
                    for w in adj[second_helix_bag]:
                        if w!=u:

                            child_table = [w for w in adj[second_helix_bag] if w!=u].pop()
                            child_indices = sorted(index2bag[child_table], key=lambda x:int(x))
                            child_indices = [i for i in child_indices if i in index2bag[second_helix_bag]]


                            if child_table[0]=='H' and not set(cluster_extremities[which_cluster[child_table]]).issubset(set(index2bag[second_helix_bag])):
                            # diag case below

                                diag_childs.append(child_table)

                                absent_ex = (set(cluster_extremities[which_cluster[child_table]])-set(index2bag[second_helix_bag])).pop() 
                                sorted_exs = sorted(cluster_extremities[which_cluster[child_table]], key=lambda x: int(x))
                                child_indices = sorted(list(set(cluster_extremities[which_cluster[child_table]]) - set([absent_ex, partner(absent_ex, sorted_exs)])),key=lambda x: int(x))
                                child_letters = [local_index_label[i] for i in child_indices]
                    
                                diag_child_letters[child_table] = child_letters

                                child_const = [local_index_label[i] for i in sorted(const_part[which_cluster[child_table]],key=lambda x : int(x))]
                    
                                diag_child_const_parts[child_table] = child_const

                            elif child_table[0]=='H' and set(cluster_extremities[which_cluster[child_table]]).issubset(set(index2bag[second_helix_bag])):
                                clique_childs.append(child_table)
                                child_letters = [local_index_label[i] for i in child_indices]
                                clique_childs_letters[child_table] = child_letters
                            else:
                            # normal table below
                                normal_childs.append(child_table)
                                child_letters = [local_index_label[i] for i in child_indices]
                                normal_child_letters[child_table] = child_letters


                if increments[0]=='-1':
                    print(bag_letter[u]+"'"+'[',letters[0]+increments[0],',',letters[1]+'|'+",".join(const),'], &\\text{if }',letters[0]+increments[0],'\\notin\{',letters[1],",",",".join(const),'\}','\\\\', file=f, end="")
                if increments[1]=='-1':
                    print(bag_letter[u]+"'"+'[',letters[0],',',letters[1]+increments[1]+'|'+",".join(const),'], &\\text{if }',letters[1]+increments[1],',\\notin\{',letters[0],",",",".join(const),'\}','\\\\', file=f, end="")
                print(bag_letter[u]+'[',letters[0]+increments[0],',',letters[1]+increments[1]+'|'+",".join(const),']+\\Delta G('+letters[0]+','+letters[1]+') &\\text{if }','\{',letters[0]+increments[0],',',letters[1]+increments[1],'\}\\cap','\{',",".join(const),'\}=\\emptyset', file=f,end="")

                print('\\end{cases}',file=f, end="")

                print('$$',file=f)
                


                print('$$', bag_letter[u], '\\left[', file=f, end="")
                print(",".join(letters),"|",",".join(const),'\\right]',file=f, end="")
                print(' =  \\min\\begin{cases}', file=f, end="")
                if increments[0]=='+1':
                    print(bag_letter[u]+'[',letters[0]+increments[0],',',letters[1]+'|'+",".join(const),'], &\\text{if }',letters[0]+increments[0],'\\notin\{',letters[1],",",",".join(const),'\}','\\\\', file=f, end="")
                if increments[1]=='+1':
                    print(bag_letter[u]+'[',letters[0],',',letters[1]+increments[1]+'|'+",".join(const),'], &\\text{if }',letters[1]+increments[1],',\\notin\{',letters[0],",",",".join(const),'\}','\\\\', file=f, end="")
                if increments[0]=='-1':
                    print(bag_letter[u]+"'"+'[',letters[0]+increments[0],',',letters[1]+'|'+",".join(const),'], &\\text{if }',letters[0]+increments[0],'\\notin\{',letters[1],",",",".join(const),'\}','\\\\', file=f, end="")
                if increments[1]=='-1':
                    print(bag_letter[u]+"'"+'[',letters[0],',',letters[1]+increments[1]+'|'+",".join(const),'], &\\text{if }',letters[1]+increments[1],',\\notin\{',letters[0],",",",".join(const),'\}','\\\\', file=f, end="")
                print(bag_letter[u]+'[',letters[0]+increments[0],',',letters[1]+increments[1]+'|'+",".join(const),']+\\Delta G('+letters[0]+','+letters[1]+') &\\text{if }','\{',letters[0]+increments[0],',',letters[1]+increments[1],'\}\\cap','\{',",".join(const),'\}=\\emptyset', file=f,end="")


                if len(adj[second_helix_bag]) >= 2:
                    print(',\\\\',file=f,end="")

                    sub_terms = []
                    for child_table in normal_childs:
                        child_letters = normal_child_letters[child_table]
                        for k, e in enumerate(child_letters):
                            if e==letters[1]:
                                child_letters[k]=e+'+1'
                        term = bag_letter[child_table]+"'["+",".join(child_letters)+']'
                        sub_terms.append(term)
                    for child_table in diag_childs:
                        child_letters = diag_child_letters[child_table]
                        for k, e in enumerate(child_letters):
                            if e==letters[1]:
                                child_letters[k]=e+'+1'
                        child_const = diag_child_const_parts[child_table]
                        term = bag_letter[child_table]+"'["+",".join(child_letters)+'|'+",".join(child_const)+']'
                        sub_terms.append(term)
                    for child_table in clique_childs:
                        child_letters = clique_childs_letters[child_table]
                        for k, e in enumerate(child_letters):
                            if e==letters[1]:
                                child_letters[k]=e+'+1'
                        child_letters[1] =  child_letters[1]+'-1'
                        child_letters[3] =  child_letters[3]+'-1'
                        term = bag_letter[child_table]+"'["+",".join(child_letters)+']'
                        sub_terms.append(term)

                    print('+'.join(sub_terms),end="",file=f)

                print('\\end{cases}',file=f, end="")

                print('$$',file=f)

        elif u[:1]=='H' and set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
        # clique case
            pass
        else: 

            indices = set(index2bag[u]).intersection(set(index2bag[prev]))
            new_vars = set(index2bag[u])-set(index2bag[prev])
    #        new_vars -= set([first_anchor])
    #        new_vars -= set([last_anchor])


            indices = [ext_to_letter[e] for e in sorted(list(indices),key=lambda x: int(x))]
            new_vars = [ext_to_letter[e] for e in sorted(list(new_vars),key=lambda x: int(x))]

            if len(indices)==0:
                print("$$",bag_letter[u], end=" ", file=f)
            else:
                print("$$",bag_letter[u]+'\\left[',",".join(indices),'\\right]',file=f,end = " ")
            
            if len(set(adj[u])-set([prev])) > 0:
                print("=\\min_{",",".join(new_vars),"}","\\left(",file=f , end=" ")
            
                terms = []

                for v in adj[u]:
                    if v!=prev:
                        if v[:1]=='H' and not set(cluster_extremities[which_cluster[v]]).issubset(set(index2bag[u])):
                            # diag bag
                            absent_ex = (set(cluster_extremities[which_cluster[v]])-set(index2bag[u])).pop() 
                            sorted_exs = sorted(cluster_extremities[which_cluster[v]], key=lambda x: int(x))
                
                            indices = sorted(list(set(cluster_extremities[which_cluster[v]]) - set([absent_ex, partner(absent_ex, sorted_exs)])),key=lambda x: int(x))
                            letters = [ext_to_letter[e] for e in indices]


                            const = []
                            for e in sorted(const_part[which_cluster[v]], key=lambda x : int(x)):
    #                            if e in indices:
    #                                const.append(ext_to_letter[e]+"'")
    #                            else:
                                const.append(ext_to_letter[e])

                            terms.append(bag_letter[v]+'\\left['+",".join(letters)+"|"+",".join(const)+'\\right]')

                        elif v[:1]=='H' and set(cluster_extremities[which_cluster[v]]).issubset(set(index2bag[u])):
                            indices_v = set(index2bag[u]).intersection(set(index2bag[v]))
                            indices_v = [ext_to_letter[e] for e in sorted(list(indices_v),key=lambda x: int(x))]
                            indices_v[1] = indices_v[1]+'-1'
                            indices_v[3] = indices_v[3]+'-1'
                            terms.append(bag_letter[v]+'\\left['+",".join(indices_v)+'\\right]')

                        else:
                            #normal bag
                            indices_v = set(index2bag[u]).intersection(set(index2bag[v]))
                            indices_v = [ext_to_letter[e] for e in sorted(list(indices_v),key=lambda x: int(x))]
                            terms.append(bag_letter[v]+'\\left['+",".join(indices_v)+'\\right]')
                print("+".join(terms), file=f, end="")


            cnt += 1

        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))
