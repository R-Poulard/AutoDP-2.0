from colors import *
from autodp.minimal_expansion import MinimalExpansion

f = open(snakemake.output[0],'w')

colors = Set3

print('\\documentclass{standalone}', file=f)
print('\\usepackage[x11names]{xcolor}', file=f)
print('\\usepackage{tikz}', file=f)
print('\\usepackage{relsize}', file=f)
print('\\begin{document}', file=f)
for k, c in enumerate(colors):
    print('\\definecolor{c'+str(k)+'}{rgb}{'+str(c)[1:-1]+'}',file=f)

print('\\begin{tikzpicture}', file=f)
inter_helix_gap=snakemake.config["inter_helix_gap"]
graph = MinimalExpansion()
graph.from_str(open(snakemake.input.dbn).readlines()[0].rstrip('\n'), inter_helix_gap=inter_helix_gap)
graph.extract_helices(inter_helix_gap=inter_helix_gap)


for k, e in enumerate(sorted(list(set.union(*[set(extremities) for extremities in graph.helices])))):
    print("\\fill ("+str(e-1)+",0) circle (0.) node[below=1cm] {\\relsize{+5}{\\textbf{"+chr(ord('a')+k)+"}}};", file=f)

for e in sorted(list(set.union(*[set(extremities) for extremities in graph.helices]))):
    print("\\node at ("+str(e-1)+",-2.5) {("+str(e)+")};", file=f)



for k in range(len(graph.vertices)-2):
    print("\\draw[thick] ("+str(k)+",0) -- ("+str(k+1)+",0);", file=f )

for k,l in graph.edges:
    if abs(k-l) <=1:
        continue
    print("\\draw[very thick] ("+str(l-1)+",0) arc \
                (0:180:"+str(abs(0.5*(k-l)))+");",file=f)

for k, extremities in enumerate(graph.helices):
    label = 'H'+str(k)

    i = int(extremities[0])-1 
    ip = int(extremities[1])-1
    jp = int(extremities[2])-1
    j = int(extremities[3])-1

    print('\\filldraw[fill=c'+str(k%len(colors))+',rounded corners] ', file=f, end="")
    print(' ('+str(i-0.2)+',-0.5) rectangle ('+str(ip+0.2)+',+0.5);', file=f)
    print('\\filldraw[fill=c'+str(k%len(colors))+',rounded corners] ', file=f, end="")
    print(' ('+str(jp-0.2)+',-0.5) rectangle ('+str(j+0.2)+',+0.5);', file=f)


print('\\end{tikzpicture}', file=f)
print('\\end{document}', file=f)
