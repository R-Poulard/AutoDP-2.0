import sys
import copy
from utils import read_td_lines
from colors import *

colors = hex_Set3

root = open(snakemake.input.tdname).readlines()[0].split(' ')[1].rstrip('\n')
print(root)

# extracting bags and tree from td file
adj, index2bag = read_td_lines(open(snakemake.input.tdname).readlines())


subgraph = {}
cluster_content = {}
which_cluster = {}
cluster_root = {}
cluster_extremities = {}
for helixline in open(snakemake.input.helix).readlines():
    label = helixline.split(' ')[0]
    extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]
    cluster_extremities[label] = extremities

all_extremities = set([])

for key, val in cluster_extremities.items():
    all_extremities = all_extremities.union(set(val))

for helixline in open(snakemake.input.helix).readlines():
    label = helixline.split(' ')[0]
    extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]

    queue = [('-1', root)]
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

        print(adj.keys())
        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))

for helixline in open(snakemake.input.helix).readlines():
    label = helixline.split(' ')[0]
   
    cluster_content[label] = []
    extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]
    cluster_extremities[label] = extremities

    subgraph[label] = ""
    queue = [('-1', root)]
    while len(queue) > 0:
        prev,u = queue.pop()
        print(label, u[:2], u) 
        if u.split('_')[0]==label:
            cluster_content[label].append(u)
            which_cluster[u] = label
            if len(subgraph[label]) > 0:
                subgraph[label] += ' -> '
                subgraph[label] += u
            else:
                cluster_root[label] = prev
                subgraph[label] += u

            if set(extremities).issubset(set(index2bag[prev])):
            # clique case
                index2bag[u] = extremities
                adj[u] = []
            else:
            # diag case
                index2bag[u] = [vert for vert in index2bag[u] if vert in all_extremities]
#                new_prev = copy.copy(prev)
#                new_u = copy.copy(u)
#                keep_going = True
#                while keep_going:
#                    keep_going = False
#                    for v in adj[new_u]:
#                        if v!=new_prev:
#                            if v.split('_')[0]==label:
#                                new_prev = copy.copy(new_u)
#                                new_u = copy.copy(v)
#                                keep_going = True
#                print("u,new_u,new_prev",u,new_prev,new_u)
##                input()
#                adj[u] = [prev, new_u]
#                adj[new_u] = [vert for vert in adj[new_u] if vert!=new_prev]+[u] 
#                prev = copy.copy(u)
#                u = copy.copy(new_u)

        for v in adj[u]:
            if v!=prev:
                queue.append((u,v))

    subgraph[label] += ';'
    if subgraph[label]==';':
        subgraph.pop(label)



const_part = {}

for key, val in cluster_content.items():
    print("cluster content ",key, val)
    if len(val) > 0:
        const_part[key] = set.intersection(*[set(index2bag[u]) for u in val])
    else:
        const_part[key] = set([])

cnt = 0
f = open(snakemake.output[0],'w')

f.write('digraph G {\n')
print("    node [shape=box];",file=f)
for key, val in subgraph.items():
    color = colors[int(key[1:])]
    if set(cluster_extremities[key]).issubset(set(index2bag[cluster_root[key]])) :
        case = ' (clique)'
    else:
        case = ' (diag)'
    print("val", val)
    print("    subgraph cluster"+str(cnt)+' {',file=f)
    print("        node [style=filled,fillcolor=white];",file=f)
    print('        labeljust="l";',file=f)
    print("        style=filled;",file=f)
    print('        color="'+color+'";',file=f)
    print("        "+val,file=f)
    print(key,set(cluster_extremities[key]))
    print(set(index2bag[cluster_root[key]]))
    ext = " ("
    for c in sorted(cluster_extremities[key],key=lambda x:int(x)):
        ext += c +"-"
    ext = ext[:-1]
    ext += ')'
    print('        label="'+key+ext+case+'";',file=f)
    print("    }",file=f)
    cnt += 1

queue = [('-1',root)]


def boldified(c, all_extremities):
    if c in all_extremities:
        return "<b>"+c+"</b>"
    else:
        return c
index2bag['-1'] = []

def num_to_letters(cnt):
    if cnt < 26:
        return chr(ord('a') + cnt).upper()
    else:
        return chr(ord('a') + int(cnt/26)).upper()+chr(ord('a') + int(cnt%26)).upper()

cnt = 0

def partner(e, sorted_ext):
    if e==sorted_ext[0]:
        return sorted_ext[3]
    if e==sorted_ext[1]:
        return sorted_ext[2]
    if e==sorted_ext[3]:
        return sorted_ext[0]
    if e==sorted_ext[2]:
        return sorted_ext[1]

while len(queue) > 0:
    prev,u = queue.pop()
    label = "<{"
    if u[:1]=='H' and not set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
    # diag case
        absent_ex = (set(cluster_extremities[which_cluster[u]]) - set(index2bag[prev])).pop()
        sorted_exs = sorted(cluster_extremities[which_cluster[u]], key=lambda x: int(x))
        if prev.split('_')[0]==u.split('_')[0]:
            for c in sorted(set([absent_ex, partner(absent_ex, sorted_exs)]),key=lambda x: int(x)):
                label += " "+boldified(c, all_extremities)
        else:
            label += '  <FONT COLOR="RED">'+num_to_letters(cnt)+'</FONT>'
            cnt += 1
            for c in set(cluster_extremities[which_cluster[u]])-set([absent_ex, partner(absent_ex, sorted_exs)]):
                label += " "+boldified(c, all_extremities)
        label += "| "
        for c in const_part[which_cluster[u]]:
            label += " "+boldified(c, all_extremities)

    elif u[:1]=='H' and set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
    # clique case
        for c in index2bag[u]:
            if c in set(index2bag[u]).intersection(set(index2bag[prev])):
                label += " "+boldified(c, all_extremities)
            else:
                label += '  <FONT COLOR="DARKGREEN">'+c+'</FONT>'

    else: 
        label += '  <FONT COLOR="RED">'+num_to_letters(cnt)+'</FONT>'
        cnt += 1

        for c in index2bag[u]:
            if c in set(index2bag[u]).intersection(set(index2bag[prev])):
                label += " "+boldified(c, all_extremities)
            else:
                label += '  <FONT COLOR="DARKGREEN">'+c+'</FONT>'
    label += "}>"
    print("    ",u,'[shape=record,label=',label,'];', file=f)
    if u.split('_')[0]!=prev.split('_')[0]:
        print("    ",prev," -> ",u+';',file=f)

    for v in adj[u]:
        if v!=prev:
            queue.append((u,v))

print('}',file=f)
