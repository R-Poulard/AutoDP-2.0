from autodp.equation_tree import TreeOfEquations, BagType

# extracting bags and tree from td file
tree_dec = TreeOfEquations()
tree_dec.read_from_file(snakemake.input.tdname)

# helices: sets a lot of useful variables
tree_dec.set_helices(open(snakemake.input.helix).readlines())

# contraction: 
tree_dec.contract_to_skeleton()

# filter to anchor vertices only
tree_dec.filter_anchors()

tree_dec.set_ext_to_letter()
tree_dec.set_dp_tables()

clique_case = False
filtered_bag_list = []
for prev, u in tree_dec.dfs_edge_iterator():
    print("iter ", prev, u)
    if tree_dec.bag_type[u]==BagType.CLIQUE:
        clique_case = True
        continue
    if tree_dec.bag_type[u]==BagType.DIAG_SECOND:
        continue
    if tree_dec.bag_type[u]!=BagType.DIAG_FIRST:
        if len(tree_dec.equations[u].indices) == 0:
            continue
    filtered_bag_list.append(u)

############ start file #############
############# BEGIN INCLUDES ########
f = open('results/c_code/C5_folding.c','w')

print('#include <stdio.h>', file = f)
print('#include <limits.h>', file = f)
print('#include <stdlib.h>', file = f)
print('#include <stdbool.h>', file = f)
print("", file=f)

############# END INCLUDES ##########
#####################################

########################################
######## FUNCTION SIGNATURES BEGIN #####
print("int min(int a, int b) { if (a<b) {return a;} else {return b;}};", file=f)
print("", file=f)

# some declarations
for u in filtered_bag_list:
    print(tree_dec.equations[u].c_index_fun_signature(tree_dec.ext_to_letter), 
          file=f)
if clique_case:
    print("int index_CLIQUE(int i, int j, int k, int l);", file=f)
print("", file=f)
# INIT
for u in filtered_bag_list:
    print('void init_fill_'+tree_dec.equations[u].main_name+'();', file=f)
if clique_case:
    print("void init_fill_CLIQUE();", file=f)
print("", file=f)

root_child = tree_dec.equations[tree_dec.bag_adj['-1'][0]]
print('int compute_'+root_child.main_name+"();",file=f)
for u in filtered_bag_list:
    indices_raw = ['int '+tree_dec.ext_to_letter[v] for v in tree_dec.equations[u].sorted_indices() if v in tree_dec.ext_to_letter]
    indices = []
    for k, i in enumerate(indices_raw):
        if i in indices_raw[:k]:
            indices.append(i+'2')
        else:
            indices.append(i)
    print('int compute_'+tree_dec.equations[u].main_name+'('+','.join(indices)+');', file=f)
    if tree_dec.bag_type[u]==BagType.DIAG_FIRST:
        print('int compute_'+tree_dec.equations[u].main_name+'2('+','.join(indices)+');', file=f)



print("int fold();\n", file=f)
######## FUNCTION SIGNATURES END #####
########################################

########################################################
######### GLOBAL VARIABLES DECLARATION BEGIN ########### 
print("int n;\n", file=f)
print("char * line = NULL;\n", file=f)
for u in filtered_bag_list:
    print("double * "+tree_dec.equations[u].main_name+';', file=f)
    if tree_dec.bag_type[u]==BagType.DIAG_FIRST:
        print("double * "+tree_dec.equations[u].main_name+'2;', file=f)
        
if clique_case:
    print("double * CLIQUE;", file=f)
    print("double * CLIQUE2;", file=f)
print("", file=f)
######### GLOBAL VARIABLES DECLARATION END ############# 
########################################################


########################################################
######### BEGIN MAIN ###################################

print("int main(int argc, char ** argv) {", file=f)
print("    char * line = NULL;", file=f)
print("    size_t len = 0;", file=f)
print('    FILE * fp = fopen(argv[1], "r");', file=f)
print("    if (fp == NULL)", file=f)
print("        exit(EXIT_FAILURE);", file=f)
print("    getline(&line, &len, fp); ", file=f)
print("    n = (int) len;")
print('    printf("%s", line);', file=f)

indent = "    "

# MALLOC
for u in filtered_bag_list:
    print(u, "equation", tree_dec.equations[u].main_name )
    tree_dec.equations[u].terminal_print(tree_dec.ext_to_letter)
    print(indent+tree_dec.equations[u].c_allocation_print(), file=f)

if clique_case:
    total_size = '*'.join(['n' for _ in range(4)])
    print("".join([indent,'double * CLIQUE = malloc(',
                   str(total_size),'*sizeof(double));']),file=f)

# INIT
for u in filtered_bag_list:
    print(indent+'init_fill_'+tree_dec.equations[u].main_name+'();', file=f)
print(indent+'init_fill_CLIQUE();', file=f)

# ACTUAL CONPUTATION
print(indent+"int score = fold();",file=f)
print(indent+"char * structure = NULL;",file=f)
print(indent+"structure = backtrace();",file=f)

# free
for u in filtered_bag_list:
    print(indent+tree_dec.equations[u].c_free_print(), file=f)

if clique_case:
    print(indent+"free(CLIQUE);", file=f)
    total_size = '*'.join(['n' for _ in range(4)])

print("}", file=f)
################ END MAIN ##############################
########################################################


for line in open('resources/bp_score.c').readlines():
    print(line, file=f, end="")

print("int fold() {",file=f)
print(indent+'compute_'+root_child.main_name+"();",file=f)
print(indent+"}\n",file=f)

if clique_case:
    for line in open('resources/clique_case_c.c').readlines():
        print(line, file=f, end="")

for u in filtered_bag_list:
    print(tree_dec.equations[u].c_array_init_fill(tree_dec.ext_to_letter), file=f)
    print(tree_dec.equations[u].c_index_function(tree_dec.ext_to_letter), file=f)

if clique_case:
    print('void init_fill_CLIQUE() {', file=f)
    print(indent+'for (int i=0; i < n;i++) {', file=f)
    print(indent+'    for (int j=i; j < n;j++) {', file=f)
    print(indent+'        for (int k=j; k < n;k++) {', file=f)
    print(indent+'            for (int l=k; l < n;l++) {', file=f)
    print(indent+'                CLIQUE[i,j,k,l] = INT_MIN;', file=f)
    print(indent+'            }', file=f)
    print(indent+'        }', file=f)
    print(indent+'    }', file=f)
    print(indent+'}',file=f)
    print('}',file=f)
    print('int index_CLIQUE(int i, int j, int k, int l) {', file=f)
    print('    return n*n*n*i+n*n*j+n*k+l;', file=f)
    print('}',file=f)

for prev, u in tree_dec.dfs_edge_iterator():
    print("",file=f)
    equation = tree_dec.equations[u]
    print("iter ", prev, u, equation.main_name, tree_dec.bag_type[u])

    if tree_dec.bag_type[u]==BagType.DIAG_FIRST:
        equation = tree_dec.equations[u]
        
        # const variables
        const = []
        for e in equation.constant_indices:
            if e in equation.variable_indices:
                const.append(tree_dec.ext_to_letter[e]+"2")
            else:
                const.append(tree_dec.ext_to_letter[e])
        
        # variable indices
        variables = [tree_dec.ext_to_letter[e] for e in equation.variable_indices]

        # booleans for \neq sep
        terms = []
        for c in const:
            terms.append('('+variables[1]+equation.increments[1]+'!='+c+')')
        eq_some_const1 = " && ".join(terms)  
        print("eq some const terms", terms)
        terms = []
        for c in const:
            terms.append('('+variables[1]+equation.increments[1]+'!='+c+')')
        eq_some_const2 = " && ".join(terms)  

        # children sum
        CHILDREN_SUM = 'INT_MAX // no children'
        if len(tree_dec.bag_adj[equation.second_bag]) >= 2:
            sub_terms = []
            for sub_eq in tree_dec.equations[equation.second_bag].subterms:
                letter_table = {}
                for e in equation.absent_indices:
                    letter_table[e] = equation.subs_table[e]
                for e in equation.variable_indices:
                    letter_table[e] = tree_dec.ext_to_letter[e]+"2"
                sub_terms.append(sub_eq.c_code_print(letter_table, tree_dec.ext_to_letter))

            CHILDREN_SUM = '+'.join(sub_terms)
        

        for line in open('resources/template_diag.c').readlines():
            line = line.replace('MAINNAME',equation.main_name)
            line = line.replace('V1_INC1_CONST_COMP', eq_some_const1)
            line = line.replace('V2_INC2_CONST_COMP', eq_some_const2)
            line = line.replace('V1',variables[0])
            line = line.replace('V2',variables[1])
            line = line.replace('INC1', equation.increments[0])
            line = line.replace('INC2', equation.increments[1])
            line = line.replace('CONST_INT', ",".join(['int '+v for v in const]))
            line = line.replace('CONST', ",".join(const))
            line = line.replace('CHILDREN_SUM', CHILDREN_SUM)
            print(line, file=f,end="")

    elif tree_dec.bag_type[u]==BagType.TRANSITIONAL:

        # iteration over new variables
        marginalization = ""
        marginalization_indices = sorted(list(equation.marginalization))
        introduced_marginalized = []
        print("marginalization", marginalization_indices)
        for k, new_variable in enumerate(marginalization_indices):
            # figuring out bounds
            minimum = '0'
            maximum = 'n'
            for constraint in equation.sorted_indices()+introduced_marginalized:
                if constraint < new_variable:
                    minimum = tree_dec.ext_to_letter[constraint]
                if constraint > new_variable:
                    maximum = tree_dec.ext_to_letter[constraint]
            new_letter = tree_dec.ext_to_letter[new_variable]
            marginalization += ''.join([indent*(k+1),'for (int ',
                           new_letter,'=',
                           minimum,
                           ';',
                           new_letter,'<',maximum,';',new_letter,'++) {\n'])
            introduced_marginalized.append(new_variable)

        print("marginalization",repr(marginalization))
        marginalization = marginalization.rstrip('\n')
        indent_marginalization = len(marginalization_indices)*indent

        # closing the for loop
        for_loop_closing = ""
        for k in range(len(marginalization_indices)-1,-1,-1):
            for_loop_closing += indent*(k+1)+'}\n'

        # children sum
        CHILDREN_SUM = 'INT_MAX // no children'
        if len(tree_dec.bag_adj[u]) >= 2:
            sub_terms = []
            for sub_eq in equation.subterms:
                sub_terms.append(sub_eq.c_code_print({}, 
                                                     tree_dec.ext_to_letter))
            CHILDREN_SUM = '+'.join(sub_terms)

        # special case: root
        skip_index_parts = False
        if len(equation.indices)==0:
            skip_index_parts = True
        index_part = False

        for line in open('resources/template_transitional.c').readlines():
            if line.find('index_start') >= 0:
                index_part = True
            if line.find('index_end') >= 0:
                index_part = False
            if index_part and skip_index_parts:
                continue
            if line.find('//') >= 0:
                continue

            line = line.replace('MAINNAME',equation.main_name)
            line = line.replace('INT_INDICES', 
                                equation.c_int_indices(tree_dec.ext_to_letter))
            line = line.replace('INDICES',
        ','.join([tree_dec.ext_to_letter[i] for i in equation.sorted_indices()]))
            line = line.replace('FOR_LOOP_NEW_VARIABLES_OPEN', marginalization)
            line = line.replace('CHILDREN_SUM', CHILDREN_SUM)
            line = line.replace('FOR_LOOP_NEW_VARIABLES_CLOSE',for_loop_closing)
            line = line.replace('INDENT', indent_marginalization)
            print(line, file=f, end="")

