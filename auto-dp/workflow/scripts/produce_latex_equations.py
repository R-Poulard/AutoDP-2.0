import sys
from autodp.equation_tree import TreeOfEquations, BagType
import copy
from utils import read_td_lines
from colors import *

colors = Set3

# extracting bags and tree from td file
tree_dec = TreeOfEquations()
print("init bag_content",tree_dec.bag_content)
tree_dec.read_from_file(snakemake.input.tdname)
print("root",tree_dec.root)

print("after reading bag_content",tree_dec.bag_content)
# helices: sets a lot of useful variables
tree_dec.set_helices(open(snakemake.input.helix).readlines())

# extracting anchors
comp_key = lambda x : int(x)
print("after set helices bag_content",tree_dec.bag_content)
first_anchor = min([min([vertex for vertex in val], key=comp_key) for key, val in tree_dec.bag_content.items() if key!='-1'], key=comp_key)
last_anchor = max([max([vertex for vertex in val], key=comp_key) for key, val in tree_dec.bag_content.items() if key!='-1'], key=comp_key)

# contraction: 
tree_dec.contract_to_skeleton()

# filter to anchor vertices only
tree_dec.filter_anchors()

# tree_dec.bag_adj, tree_dec.bag_content, tree_dec.const_part, tree_dec.which_helix, cluster_root, tree_dec.cluster_extremities, all_extremities have been
# constructed

# LATEX_PREAMBULE BEGIN
f = open(snakemake.output[0],'w')

print('\\documentclass{article}', file=f)
print('\\usepackage[x11names]{xcolor}', file=f)
print('\\usepackage{amsmath}', file=f)
print('\\usepackage{graphicx}', file=f)
print('\\usepackage{amssymb}', file=f)
print('\\begin{document}', file=f)
print('\\textbf{fatgraph name: '+snakemake.wildcards.family+'}', file=f)

for k, c in enumerate(colors):
    print('\\definecolor{c'+str(k)+'}{rgb}{'+str(c)[1:-1]+'}',file=f)

print('\\begin{center}', file=f)
print('\\begin{figure}[h]', file=f)
print('\\includegraphics[width=\\textwidth]{'+snakemake.input.colored_dbn+'}', file=f)
print('\\end{figure}', file=f)
print('\\end{center}', file=f)

tree_dec.set_ext_to_letter()
print('ext to letter', tree_dec.ext_to_letter)
print('first and last anchors, already given: $',tree_dec.ext_to_letter[first_anchor]
                                                ,','
                                                ,tree_dec.ext_to_letter[last_anchor]
                                                ,'$'
                                                ,file=f)
# LATEX_PREAMBULE END

# giving a name and a color to each table
tree_dec.set_dp_tables()

for prev,u in tree_dec.dfs_edge_iterator(no_minus_one=False):
    print("",prev, u)
    print(u,tree_dec.equations[u].latex_print({},tree_dec.ext_to_letter))
    if tree_dec.bag_type[u]==BagType.DIAG_FIRST:

        # precomputed by set dp table
        equation = tree_dec.equations[u]

        # const variables
        const = []
        for e in equation.constant_indices:
            if e in equation.variable_indices:
                const.append(tree_dec.ext_to_letter[e]+"'")
            else:
                const.append(tree_dec.ext_to_letter[e])

        # variable indices
        variables = [tree_dec.ext_to_letter[e] for e in equation.variable_indices]

        print('$$', equation.latex_name+"'", '\\left[', file=f, end="")
        print(",".join(variables),"\mid ",",".join(const),'\\right]',file=f, end="")
        print(' =  \\min\\begin{cases}', file=f, end="")

        if not equation.inward:
            print(equation.latex_name+"'"+'[',
                  variables[0]+equation.increments[0],
                  ',',
                  variables[1]+'\mid '+",".join(const),
                  '], &\\text{if }',
                  variables[0]+equation.increments[0],
                  '\\notin\{',variables[1],
                  ",",
                  ",".join(const),
                  '\}',
                  '\\\\', 
                  file=f, end="")
        else:
            print(equation.latex_name+"'"+'[',
                  variables[0],
                  ',',
                  variables[1]+equation.increments[1]+'\mid '+",".join(const),
                  '], &\\text{if }',
                  variables[1]+equation.increments[1],
                  ',\\notin\{',
                  variables[0],
                  ",",
                  ",".join(const),
                  '\}',
                  '\\\\', 
                  file=f, end="")

            print(equation.latex_name+'[',
                  variables[0]+equation.increments[0],
                  ',',
                  variables[1]+equation.increments[1]+'\mid '+",".join(const),
                  ']+\\Delta G('+variables[0]+','+variables[1]+') &\\text{if }',
                  '\{',variables[0]+equation.increments[0],
                  ',',
                  variables[1]+equation.increments[1],
                  '\}\\cap',
                  '\{',
                  ",".join(const),'\}=\\emptyset', file=f,end="")

        print('\\end{cases}',file=f, end="")
        print('$$',file=f)
        
        print('$$', equation.latex_name, '\\left[', file=f, end="")
        print(",".join(variables),"\mid ",",".join(const),'\\right]',file=f, end="")
        print(' =  \\min\\begin{cases}', file=f, end="")

        if equation.inward:

            print(equation.latex_name+'[',
                  variables[0]+equation.increments[0],
                  ',',variables[1]+'\mid '+",".join(const),
                  '], &\\text{if }',
                  variables[0]+equation.increments[0],
                  '\\notin\{',variables[1],
                  ",",
                  ",".join(const),
                  '\}','\\\\', file=f, end="")

        else:
            print(equation.latex_name+'[',
                  variables[0],
                  ',',variables[1]+equation.increments[1]+'\mid '+",".join(const),
                  '], &\\text{if }',
                  variables[1]+equation.increments[1],
                  ',\\notin\{',
                  variables[0],
                  ",",
                  ",".join(const),'\}','\\\\', file=f, end="")

        if not equation.inward:
            print(equation.latex_name+"'"+'[',
                  variables[0]+equation.increments[0],
                  ',',
                  variables[1]+'\mid '+",".join(const),
                  '], &\\text{if }',
                  variables[0]+equation.increments[0],
                  '\\notin\{',variables[1],
                  ",",",".join(const),'\}','\\\\', file=f, end="")

        else:
            print(equation.latex_name+"'"+'[',
                  variables[0],
                  ',',
                  variables[1]+equation.increments[1]+'\mid '+",".join(const),
                  '], &\\text{if }',variables[1]+equation.increments[1],
                  ',\\notin\{',variables[0],
                  ",",
                  ",".join(const),'\}','\\\\', file=f, end="")

        print(equation.latex_name+'[',
              variables[0]+equation.increments[0],
              ',',
              variables[1]+equation.increments[1]+'\mid '+",".join(const),
              ']+\\Delta G('+variables[0]+','+variables[1]+') &\\text{if }',
              '\{',
              variables[0]+equation.increments[0],
              ',',
              variables[1]+equation.increments[1],
              '\}\\cap',
              '\{',
              ",".join(const),'\}=\\emptyset', file=f,end="")

        
        if len(tree_dec.bag_adj[equation.second_bag]) >= 2:
            sub_terms = []
            for sub_eq in tree_dec.equations[equation.second_bag].subterms:
                letter_table = {}
                for e in equation.absent_indices:
                    letter_table[e] = tree_dec.ext_to_letter[equation.subs_table[e]]
                for e in equation.variable_indices:
                    letter_table[e] = tree_dec.ext_to_letter[e]+"'"
                sub_terms.append(sub_eq.latex_print(letter_table, tree_dec.ext_to_letter))

            print(',\\\\',file=f,end="")
            print('+'.join(sub_terms),end="",file=f)

        print('\\end{cases}',file=f, end="")

        print('$$',file=f)

    elif tree_dec.bag_type[u]==BagType.CLIQUE:
    # clique case
        pass
    elif tree_dec.bag_type[u]==BagType.TRANSITIONAL:

        equation = tree_dec.equations[u]

        indices = [tree_dec.ext_to_letter[e] for e in sorted(list(equation.indices),key=lambda x: int(x))]
        new_vars = [tree_dec.ext_to_letter[e] for e in sorted(list(equation.marginalization),key=lambda x: int(x)) if e not in [first_anchor, last_anchor]]

        if len(indices)==0:
            print("$$",equation.latex_name, end=" ", file=f)
        else:
            print("$$",equation.latex_name+'\\left[',",".join(indices),'\\right]',file=f,end = " ")
        
        if len(equation.subterms) > 0:
            print("=\\min_{",",".join(new_vars),"}","\\left(",file=f , end=" ")
            terms = [eq.latex_print({}, tree_dec.ext_to_letter) for eq in equation.subterms]
            print("+".join(terms), file=f, end="")


        print("\\right)", file=f, end=" ")

        print("$$", file =f)

#print('$$',end=" ",file=f)
#print("C_\\boxtimes'[i,i',j',j]=\\begin{cases} ",file=f,end=" ")
#print("C_\\boxtimes'[i,i',j',j-1] \\\\",file=f)
#print("C_\\boxtimes[i+1,i',j',j-1] + \\Delta G(i,j) \\\\",file=f)
#print("\\Delta G(i,j) \\\\",file=f)
#print('\\end{cases}',file=f)
#print('$$',end=" ",file=f)
#
#print('$$',end=" ",file=f)
#print('C_\\boxtimes=\\begin{cases}',file=f,end=" ")
#
#print('end{cases}',file=f)
#print('$$',end=" ",file=f)


print('\\end{document}',file=f)
