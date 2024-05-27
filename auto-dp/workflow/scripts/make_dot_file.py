import sys

Paired = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928']
Accent = ['#7FC97F', '#BEAED4', '#FDC086', '#FFFF99', '#386CB0', '#F0027F', '#BF5B16', '#666666']
Dark2 = ['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666']
Set1 = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999']
Set2 = ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3']
Set3 = ['#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F']
tab10 = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']
tab20 = ['#1F77B4', '#AEC7E8', '#FF7F0E', '#FFBB78', '#2CA02C', '#98DF8A', '#D62728', '#FF9896', '#9467BD', '#C5B0D5', '#8C564B', '#C49C94', '#E377C2', '#F7B6D2', '#7F7F7F', '#C7C7C7', '#BCBD22', '#DBDB8D', '#17BECF', '#9EDAE5']
tab20b = ['#393B79', '#5254A3', '#6B6ECF', '#9C9EDE', '#637939', '#8CA252', '#B5CF6B', '#CEDB9C', '#8C6D31', '#BD9E39', '#E7BA52', '#E7CB94', '#843C39', '#AD494A', '#D6616B', '#E7969C', '#7B4173', '#A55194', '#CE6DBD', '#DE9ED6']
tab20c = ['#3182BD', '#6BAED6', '#9ECAE1', '#C6DBEF', '#E6550D', '#FD8D3C', '#FDAE6B', '#FDD0A2', '#31A354', '#74C476', '#A1D99B', '#C7E9C0', '#756BB1', '#9E9AC8', '#BCBDDC', '#DADAEB', '#636363', '#969696', '#BDBDBD', '#D9D9D9']

colors = Set3

# extracting bags and tree from td file
index2bag = {}
adj = {}

started_bs = False

root = open(snakemake.input.tdname).readlines()[0].split(' ')[1].rstrip('\n')

for line in open(snakemake.input.tdname).readlines():
    if line[0]=='b':
        started_bs = True
        index2bag[line.split(' ')[1].rstrip('\n')] = [vertex.rstrip('\n') for vertex in line.split(' ')[2:]]

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

root = list(adj.keys())[0]

f = open(snakemake.output[0],'w')

f.write('digraph G {\n')
print("    node [shape=box];",file=f)
subgraph = {}
cluster_content = {}
which_cluster = {}
cluster_root = {}
cluster_extremities = {}
for helixline in open(snakemake.input.helix).readlines():
    label = helixline.split(' ')[0]
   
    cluster_content[label] = []
    extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]
    cluster_extremities[label] = extremities

    subgraph[label] = ""
    queue = [('-1',root)]
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
def color(u):
    return colors[int(u.split('_')[0][1:])%len(colors)]
for key, val in subgraph.items():
    print("    subgraph cluster"+str(cnt)+' {',file=f)
    print("        node [style=filled,fillcolor=white];",file=f)
    print('        labeljust="l";',file=f)
    print("        style=filled;",file=f)
    print('        color="'+str(colors[int(key[1:])%len(colors)])+'";',file=f)
    print("        "+val,file=f)
    print(key,set(cluster_extremities[key]))
    print(set(index2bag[cluster_root[key]]))
    if set(cluster_extremities[key]).issubset(set(index2bag[cluster_root[key]])) :
        case = ' (clique)'
    else:
        case = ' (diag)'
    print(case)
    ext = " ("
    for c in sorted(cluster_extremities[key],key=lambda x:int(x)):
        ext += c +"-"
    ext = ext[:-1]
    ext += ')'
    print('        label="'+key+ext+case+'";',file=f)
    print("    }",file=f)
    cnt += 1

queue = [('-1',root)]


def color(u):
    return colors[int(u.split('_')[0][1:])%len(colors)]

all_extremities = set([])

def boldified(c, all_extremities):
    if c in all_extremities:
        return "<b>"+c+"</b>"
    else:
        return c

for key, val in cluster_extremities.items():
    all_extremities = all_extremities.union(set(val))

while len(queue) > 0:
    prev,u = queue.pop()
    label = "<{"
    if u[:1]=='H':
        print(index2bag[u])
        print(const_part[which_cluster[u]])
        label+=''
        for c in set(index2bag[u])-const_part[which_cluster[u]]:
            label += " "+boldified(c, all_extremities)
        label += "| "
        for c in const_part[which_cluster[u]]:
            label += " "+boldified(c, all_extremities)

    else:
        print(u,index2bag[u])
        for c in index2bag[u]:
            label += " "+boldified(c, all_extremities)
#    if u[0]=='H':
#        print(u,color(u))
#        print("    ",u,'[color="'+color(u)+'",style=filled,label="',label,'"];', file=f)
#    else:
    label += "}>"
    print("    ",u,'[shape=record,label=',label,'];', file=f)
    if u.split('_')[0]!=prev.split('_')[0] and prev!='-1':
        print("    ",prev," -> ",u+';',file=f)

    for v in adj[u]:
        if v!=prev:
            queue.append((u,v))


print('}',file=f)
