import sys
import os
import copy
from colors import *

colors = hex_Set3
from utils import read_td_lines, full_compact_diag, full_compact_clique, subgraph_info_extraction, filter_extremities, equations_prep_work


# constant part of each helix

adj, all_extremities, cluster_root, cluster_content, cluster_extremities, extremities_label, index2bag, index2letters, subgraph, which_cluster = equations_prep_work(snakemake.input.tdname, snakemake.input.helix) 
const_part = {}

for key, val in cluster_content.items():
    print("cluster content ",key, val)
    if len(val) > 0:
        const_part[key] = set.intersection(*[set(index2bag[u]) for u in val])
    else:
        const_part[key] = set([])

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

def boldified(c, extremities_label, inc=None):
    if inc:
        return "<b>"+extremities_label[c]+inc+"</b>"
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
        bag_letter[u] = num_to_letters(cnt)
        cnt += 1
    else:
        if set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
        # clique
            bag_letter[u] = 'CLIQUE'
        else:
        # diag
            bag_letter[u] = num_to_letters(cnt)
            cnt += 1

    for v in adj[u]:
        if v!=prev:
            queue.append((u,v))


queue = [('-1','1')]
index2bag['-1'] = []
index2letters['-1'] = []

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


while len(queue) > 0:
    prev,u = queue.pop()
    label = "<{"
    label += "}>"
    image_name = snakemake.wildcards.family+'_'+u+'.jpg'

    latex_string = ""

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

            latex_string += '$\displaystyle '+bag_letter[u]+'\\left['
            latex_string += ",".join(letters)+"|"+",".join(const)+'\\right]'
            latex_string += '\\min\\left('
            latex_string += bag_letter[u]+'['+letters[0]+increments[0]+','+letters[1]+'|'+",".join(const)+'],'
            latex_string += bag_letter[u]+'['+letters[0]+','+letters[1]+increments[1]+'|'+",".join(const)+'],'
            latex_string += bag_letter[u]+'['+letters[0]+increments[0]+','+letters[1]+increments[1]+'|'+",".join(const)+']+bp('+letters[0]+','+letters[1]+')'

            if len(adj[u]) >= 2:
                latex_string += ','
                child_table = [v for v in adj[u] if v!=prev].pop()
                child_indices = sorted(index2bag[child_table], key=lambda x:int(x))
                child_indices = [i for i in child_indices if i in local_index_label.keys()]

                if child_table[0]=='H' and not set(cluster_extremities[which_cluster[child_table]]).issubset(set(index2bag[u])):
                # diag case below
                    absent_ex = (set(cluster_extremities[which_cluster[child_table]])-set(index2bag[u])).pop() 
                    sorted_exs = sorted(cluster_extremities[which_cluster[child_table]], key=lambda x: int(x))
                    child_indices = sorted(list(set(cluster_extremities[which_cluster[child_table]]) - set([absent_ex, partner(absent_ex, sorted_exs)])),key=lambda x: int(x))
                    child_letters = [local_index_label[i] for i in child_indices]
                    child_const = [local_index_label[i] for i in sorted(const_part[which_cluster[child_table]],key=lambda x : int(x))]
        
                    latex_string += bag_letter[child_table]+'['+",".join(child_letters)+'|'+",".join(child_const)+']'

                else:
                # normal table below
                    child_letters = [local_index_label[i] for i in child_indices]
                    latex_string += bag_letter[child_table]+'['+",".join(child_letters)+']'

            latex_string += '\\right) $'

    elif u[:1]=='H' and set(cluster_extremities[which_cluster[u]]).issubset(set(index2bag[prev])):
    # clique case
        indices = set(index2bag[u]).intersection(set(index2bag[prev]))
        indices = [ext_to_letter[e] for e in sorted(list(indices),key=lambda x: int(x))]
        latex_string += "$\displaystyle "+bag_letter[u]+'\\left['+",".join(indices)+'\\right] $'
    else: 

        indices = set(index2bag[u]).intersection(set(index2bag[prev]))
        new_vars = set(index2bag[u])-set(index2bag[prev])

        indices = [ext_to_letter[e] for e in sorted(list(indices),key=lambda x: int(x))]
        new_vars = [ext_to_letter[e] for e in sorted(list(new_vars),key=lambda x: int(x))]

        if len(indices)==0:
            latex_string += "$\displaystyle "+bag_letter[u]
        else:
            latex_string += "$\displaystyle "+bag_letter[u]+'\\left['+",".join(indices)+'\\right]'
        
        if len(adj[u]) > 0:
            latex_string += "=\\min_{"+",".join(new_vars)+"}"+"\\left("
        
            terms = []

            for v in adj[u]:
                if v!=prev:
                    indices_v = set(index2bag[u]).intersection(set(index2bag[v]))
                    indices_v = [ext_to_letter[e] for e in sorted(list(indices_v),key=lambda x: int(x))]
                    terms.append(bag_letter[v]+'\\left['+",".join(indices_v)+'\\right]')
            latex_string += "+".join(terms)


        latex_string += "\\right) $"

    with open('miscellani_latex/'+image_name.split('.')[0]+'.tex','w') as latex_file:
        print("\\documentclass{standalone}",file=latex_file)
        print("\\begin{document}",file=latex_file)
        print(latex_string,file=latex_file)
        print("\\end{document}",file=latex_file)
    os.system('latex -output-directory=miscellani_latex/ miscellani_latex/'+image_name.split('.')[0]+'.tex')
    os.system('dvips -o miscellani_latex/'+image_name+' miscellani_latex/'+image_name.split('.')[0]+'.dvi')
    os.system('gs -sDEVICE=jpeg -dJPEGQ=100 -dNOPAUSE -dBATCH -dSAFER -r900 -sOutputFile=miscellani_latex/'+image_name+' miscellani_latex/'+image_name.split('.')[0]+'.ps')

    print("    ",u,'[image="miscellani_latex/'+image_name+'"];', file=f)
    if u.split('_')[0]!=prev.split('_')[0]:
        print("    ",prev," -> ",u+';',file=f)

    for v in adj[u]:
        if v!=prev:
            queue.append((u,v))

label = "<{"
for c in sorted(list(all_extremities),key=lambda x:int(x)):
    label += boldified(c, extremities_label)
    label += '=&#956;('+c+')'+' | '
label += "}>"

print("    labels",'[shape=record,label=',label,'];', file=f)

print('}',file=f)
