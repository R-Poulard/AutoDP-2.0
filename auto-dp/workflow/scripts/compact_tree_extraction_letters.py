import sys
import copy
from colors import *

colors = hex_Set3
from utils import read_td_lines, full_compact_diag, full_compact_clique, subgraph_info_extraction, filter_extremities, equations_prep_work


# constant part of each helix
root = open(snakemake.input.tdname).readlines()[0].split(' ')[1].rstrip('\n')
print(root)

adj, all_extremities, cluster_root, cluster_content, cluster_extremities, extremities_label, index2bag, index2letters, subgraph, which_cluster = equations_prep_work(snakemake.input.tdname, snakemake.input.helix, root) 

# STARTING TO WRITE OUTPUT DOT FILE
f = open(snakemake.output[0],'w')

f.write('digraph G {\n')
print("    node [shape=box];",file=f)

cnt = 0
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
    ext = " ("
    for c in sorted(cluster_extremities[key],key=lambda x:int(x)):
        ext += c +"-"
    ext = ext[:-1]
    ext += ')'
    print('        label="'+key+ext+case+'";',file=f)
    print("    }",file=f)
    cnt += 1


def rank(c,all_extremities):
    for k, val in enumerate(sorted(list(all_extremities), key=lambda x: int(x))):
        if val==c:
            return k

def boldified(c, all_extremities, inc=None):
    if c in all_extremities:
        if inc:
            return "<b>"+chr(ord('a')+rank(c, all_extremities))+inc+"</b>"
        else:
            return "<b>"+chr(ord('a')+rank(c, all_extremities))+"</b>"
    else:
        return "<b>"+extremities_label[c]+"</b>"


def num_to_letters(cnt):
    if cnt < 26:
        return chr(ord('a') + cnt).upper()
    else:
        return chr(ord('a') + int(cnt/26)).upper()+chr(ord('a') + int(cnt%26)).upper()

cnt = 0

bag_letter = {}

queue = [('-1','1')]

while len(queue) > 0:
    prev,u = queue.pop()
    if not u[:1]=='H':
        bag_letter[u] = '<FONT COLOR="RED"><b>'+num_to_letters(cnt)+'</b></FONT>'
        cnt += 1
    else:
        if set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
        # clique
            bag_letter[u] = '<FONT COLOR="RED"><b>'+'CLIQUE'+'</b></FONT>'
        else:
        # diag
            bag_letter[u] = '<FONT COLOR="RED"><b>'+'DIAG'+'</b></FONT>'

    for v in adj[u]:
        if v!=prev:
            queue.append((u,v))

queue = [('-1','1')]
index2bag['-1'] = []
index2letters['-1'] = []

while len(queue) > 0:
    prev,u = queue.pop()
    label = "<{"
    if u[:1]=='H' and not set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
    # diag case
        label += '  <FONT COLOR="RED"><b>'+num_to_letters(cnt)+'</b></FONT>'
        cnt += 1
        for c in set(index2bag[u]):
            if c in cluster_extremities[which_cluster[u]]:
                label += " <b>"+extremities_label[c]+" </b>"
        label += "| "
        for c in set(index2bag[u]):
            if c not in cluster_extremities[which_cluster[u]]:
                label += " <b>"+extremities_label[c]+" </b>"

    elif u[:1]=='H' and set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
        for c in sorted(index2letters[u]):
            label += " <b>"+c+" </b>"

    else: 
    # other bags
        label += '  <FONT COLOR="RED"><b>'+num_to_letters(cnt)+'</b></FONT>'
        cnt += 1

        for c in index2letters[u]:
            if c in set(index2letters[u]).intersection(set(index2letters[prev])):
                label += " <b>"+c+" </b>"
            else:
                label += ' <b><FONT COLOR="DARKGREEN">'+c+" </FONT></b>"
    label += "}>"
    print("    ",u,'[shape=record,label=',label,'];', file=f)
    if u.split('_')[0]!=prev.split('_')[0]:
        print("    ",prev," -> ",u+';',file=f)

    for v in adj[u]:
        if v!=prev:
            queue.append((u,v))

label = "<{"
for c in sorted(list(all_extremities),key=lambda x:int(x)):
    label += "<b>"+extremities_label[c]+"</b>"
    label += '=&#956;('+c+')'+' | '
label += "}>"

print("    labels",'[shape=record,label=',label,'];', file=f)

print('}',file=f)
