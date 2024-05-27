from colors import *
def find_helices(adj):

    # basal indices (largest overarching arc)
    basal_indices = []

    for i in range(1, max(adj.keys())+1,1):
        for j in range(i+1, max(adj.keys())+1,1):
            if j in adj[i] and abs(j-i) > 1:
                try:
                    if j+1 in adj[i-1]:
                        continue
                    else:
                        basal_indices.append((i,j))
                except KeyError:
                    basal_indices.append((i,j))


    print("basal indices", basal_indices)
    
    helices = []
    # complete
    for u, v in basal_indices:
        w = u 
        x = v
        while x-1 in adj[w+1] and w+1 < x-1:
            w+=1
            x-=1

        helices.append((u,v,w,x))

    return helices

f = open(snakemake.output[0],'w')



dbn = open(snakemake.input.dbn).readlines()[0]

symbols = ['()','[]','{}','<>','Aa','Bb','Cc']

stacks = {}

for s in symbols:
    stacks[s] = []

edges = set([])
vertices =set([])

k = 0
int_label = {}
extremities_pos = [0]

dbn_list = list(dbn)
adj = {}

pos = 0
while len(dbn_list) > 0:
    print(dbn_list)
    int_label[k+1] = pos+1
    c = dbn_list.pop(0)
    vertices.add(k+1)

    if len(dbn_list) > 0:
        vertices.add((k+2))
        edges.add((k+1,k+2))
        try:
            adj[k+1].append(k+2)
        except KeyError:
            adj[k+1] = [k+2]
        try:
            adj[k+2].append(k+1)
        except KeyError:
            adj[k+2] = [k+1]
     
    for s in symbols:
        if c==s[0]:
            stacks[s].append(k+1)
        if c==s[1]:
            l = stacks[s].pop()
            edges.add((l, k+1))
            try:
                adj[k+1].append(l)
            except KeyError:
                adj[k+1] = [l]
            try:
                adj[l].append(k+1)
            except KeyError:
                adj[l] = [k+1]

    if len(dbn_list) > 0 and c==dbn_list[0]:
        pos += 1
    else:
        extremities_pos.append(k)
    k+= 1


helices = find_helices(adj)

pos_to_helix = {}
helix_to_pos = {}

all_pos = []
for k, helix in enumerate(helices):
    label = 'H'+str(k)
    poss = []
    for i in sorted(list(helix)):
        i = int(i)
        pos_to_helix[i] = label
        poss.append(i)
    helix_to_pos[label] = poss
    all_pos += poss

print('\\documentclass{standalone}',file=f)
print('\\usepackage{tikz}',file=f)
print('\\begin{document}',file=f)
print('\\begin{tikzpicture}',file=f)

print(pos_to_helix)
for k, pos in enumerate(sorted(all_pos)):
    label = pos_to_helix[pos]
    if pos==min(helix_to_pos[label]):
        m = min(helix_to_pos[label])
        M = max(helix_to_pos[label])
        kp = -1
        for kpp, val  in enumerate(sorted(all_pos)):
            if all_pos[kpp]==M:
                kp = kpp
                break

        print('\\filldraw[fill=black!20!white] (2*', file=f, end="")
        print(str(k)+',0) arc (180:0:'+str(kp-k+0.5)+') -- ',file=f,end="")
        print('(2*'+str(kp)+'-1,0)', file=f, end="")
        print(' arc (0:180:'+str(kp-k-0.5)+') -- cycle;', file=f)

print('\\end{tikzpicture}',file=f)
print('\\end{document}',file=f)
