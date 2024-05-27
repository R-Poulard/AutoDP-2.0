f = open(snakemake.output,'w')

print('\\documentclass{article}',file=f)
print('\\usepackage{graphics}',file=f)
print('\\begin{document}',file=f)


for k in range(len(snakemake.input.tdimages)):


print('\\end{document}',file=f)
