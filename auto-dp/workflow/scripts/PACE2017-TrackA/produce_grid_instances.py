import os
import sys

lateral_size = int(sys.argv[-1])

outfile = open(str(lateral_size)+'x'+str(lateral_size)+'.gr','w')

vertices = set([])
edges = set([])

index = 1

kl_to_ind = {}
ind_to_kl = {}

for k in range(lateral_size):
    for l in range(lateral_size):
        vertices.add(index)
        kl_to_ind[(k,l)] = index
        ind_to_kl[index] = (k,l)
        index += 1

for k in range(lateral_size):
    for l in range(lateral_size):
        if k < lateral_size-1:
            edges.add((kl_to_ind[(k,l)],kl_to_ind[(k+1,l)]))
        if l < lateral_size-1:
            edges.add((kl_to_ind[(k,l)],kl_to_ind[(k,l+1)]))
        if k > 0:
            edges.add((kl_to_ind[(k,l)],kl_to_ind[(k-1,l)]))
        if l > 0:
            edges.add((kl_to_ind[(k,l)],kl_to_ind[(k,l-1)]))

outfile.write('p tw '+str(len(vertices))+' '+str(len(edges))+'\n')
for u,v in edges:
    outfile.write(str(u)+' '+str(v)+'\n')

outfile.close()
