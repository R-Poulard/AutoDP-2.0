from autodp.equation_tree import TreeOfEquations, BagType

PSEUDO="M"
print("PSEUDO=",PSEUDO)
# extracting bags and tree from td file
tree_dec = TreeOfEquations()
tree_dec.read_from_file('results/processed_td_files/processed_'+PSEUDO+'.td')

# helices: sets a lot of useful variables
tree_dec.set_helices(open('results/helix_annotations/'+PSEUDO+'.helix').readlines())
# contraction: 
tree_dec.contract_to_skeleton()
# filter to anchor vertices only
tree_dec.filter_anchors()
if tree_dec.root[0]=="H":
    new_root='0'
    tmp=tree_dec.root
    tree_dec.bag_adj[new_root]=['-1',tree_dec.root]
    tree_dec.bag_adj[tree_dec.root].remove('-1')
    tree_dec.root=new_root
    bg=[int(x) for x in tree_dec.bag_content[tmp]]
    tree_dec.bag_content[new_root]=tree_dec.bag_content[tmp]

tree_dec.set_ext_to_letter()
tree_dec.set_dp_tables()
sorted_keys = sorted(tree_dec.ext_to_letter.keys(), key=int)

sorted_letters=[tree_dec.ext_to_letter[key] for key in sorted_keys]

reverse_ext_to_letter={tree_dec.ext_to_letter[key]:key for key in tree_dec.ext_to_letter.keys()}
reverse_equation={tree_dec.equations[key]:key for key in tree_dec.bag_type.keys()}
children_depth={i[1]:0 for i in tree_dec.equations.items()}
print(sorted_keys)
print(sorted_letters)

def depth(a,prev):
    if len(tree_dec.bag_adj[a])<=1:
        print(a,"ended")
        return 1
    print(a)
    size=0
    for i in tree_dec.bag_adj[a]:
        if i!=a and i!=prev:
            size+=depth(i,a)
    children_depth[tree_dec.equations[a]]=size
    if tree_dec.bag_type[a]==BagType.DIAG_FIRST:
        return size
    else:
        return size+1
depth(tree_dec.root,'-1')

header=""#will be taken from a template
declarations=["//declarations"]
main="//Here would be main\n"#main should be a constant as you don't have variable allocation (only for the hash table)
functions=[]
backtrace=[]
#principle, declarations and function will be apppended during the tree parcours
#then joined by space and put in the file at once
clique=False
bag_size=[0] #length of the max tuple put as key into the hashing table

def print_children(lst,x,indent_marginalization):
    if len(lst) == 1:
        return indent_marginalization+"  int tmp"+str(x)+"= "+lst[0]+";\n"+indent_marginalization+"  if(tmp"+str(x)+"==INT_MAX){continue;}\n"
    else:
        return indent_marginalization+"  int tmp"+str(x)+"= "+lst[-1]+";\n"+indent_marginalization+"  if(tmp"+str(x)+"==INT_MAX){continue;}\n"+print_children(lst[:-1],x+1,indent_marginalization)

def print_children_diag(lst,x,indent_marginalization):
    if len(lst) == 1:
        return indent_marginalization+"  int tmp"+str(x)+"= "+lst[0]+";\n"+indent_marginalization+"  if(tmp"+str(x)+"==INT_MAX){return min_value;}\n"
    else:
        return indent_marginalization+"  int tmp"+str(x)+"= "+lst[-1]+";\n"+indent_marginalization+"  if(tmp"+str(x)+"==INT_MAX){return min_value;}\n"+print_children(lst[:-1],x+1,indent_marginalization)

def print_nested(x):
    if x== 0:
        return "tmp"+str(x)
    else:
        return f"add({print_nested(x-1)},tmp{x})"
                    
def treat_bag(prev,node):  
    child_info=[]

    #print("(",prev,")",node,tree_dec.bag_type[node],tree_dec.bag_content[node])

    for w in tree_dec.bag_adj[node]:
        ##print(tree_dec.bag_adj[node])
        if w!=prev and w!=node:
            treat_bag(node,w)
    
    if tree_dec.bag_type[node]==BagType.CLIQUE:
        ##print("clique changed here")
        global clique
        clique=True #no need to do more as the clique implementation will be added to the lists in the beginning once the tree done
        if bag_size[0]<5:
            #print("chgmt in clique")
            bag_size[0]=5
        pass
    elif tree_dec.bag_type[node]==BagType.TRANSITIONAL:
         # iteration over new variables
        equation = tree_dec.equations[node]
        indent="    "
        marginalization = ""
        marginalization_indices = sorted(list(equation.marginalization),key=int)
        introduced_marginalized = []
        #print("marginalization", marginalization_indices)
        for k, new_variable in enumerate(marginalization_indices):
            # figuring out bounds
            minimum = 'START'
            minimum_val=-1
            maximum = 'MAX'
            maximum_val=9999999
            #print("CONSTRAINT",equation.sorted_indices())
            #print("margi",introduced_marginalized)
            for constraint in equation.sorted_indices()+introduced_marginalized:
                if int(constraint) < int(new_variable) and int(constraint)>int(minimum_val):
                    #print("loooo",tree_dec.ext_to_letter[constraint],add_on[tree_dec.ext_to_letter[constraint]])
                    minimum = tree_dec.ext_to_letter[constraint]#+add_on[tree_dec.ext_to_letter[constraint]]
                    minimum_val = int(constraint)
                if int(constraint) > int(new_variable) and int(constraint)<int(maximum_val):
                    maximum = tree_dec.ext_to_letter[constraint]
                    maximum_val=int(constraint)
            new_letter = tree_dec.ext_to_letter[new_variable]
            add_before=""
            add_after=""
            
            idx_maximum=sorted_letters.index(maximum)
            idx_new=sorted_letters.index(new_letter)

            if minimum!='START':
                cmp=0
                idx_minium=sorted_letters.index(minimum)
                for i in range(idx_new-1,-1,-1):
                    if i%2==0:
                        cmp+=1
                    if i==idx_minium:
                        break
                if cmp!=0:
                    add_before="+"+str(cmp)
            
            if idx_new%2==1 and idx_new==idx_maximum-1:
                add_after="+1"
            else:
                cmp=0
                for i in range(idx_new+idx_new%2,idx_maximum):
                    #print("I=",i)
                    if i==idx_new+idx_new%2:
                        continue
                    if i==idx_maximum:
                        break
                    if i%2==0:
                        #print("CMP+1")
                        cmp+=1
                if cmp!=0:
                    add_after="-"+str(cmp)
            #print("result=",add_after)            
            marginalization += ''.join([indent*(k+1),'for (int ',
                           new_letter,'=',
                           minimum,add_before,
                           ';',
                           new_letter,'<',maximum,add_after,';',new_letter,'++) {\n'])
            introduced_marginalized.append(new_variable)

        ##print("marginalization",repr(marginalization))
        marginalization = marginalization.rstrip('\n')
        indent_marginalization = len(marginalization_indices)*indent

        # closing the for loop
        for_loop_closing = ""
        for k in range(len(marginalization_indices)-1,-1,-1):
            for_loop_closing += indent*(k+1)+'}\n'

        # children sum
        CHILDREN_SUM = ''
        CHILDREN_BACKTRACE=''
        if len(tree_dec.bag_adj[node]) >= 2:
            sub_terms = []
            bag_types=[]
            for sub_eq in sorted(equation.subterms,key= lambda x:children_depth[x]):
                sub_terms.append(sub_eq.c_code_print({},tree_dec.ext_to_letter).replace("(","( hashTable,"))
                bag_types.append(tree_dec.bag_type[reverse_equation[sub_eq]])

            for idv in range(0,len(sub_terms)):
                i=sub_terms[len(sub_terms)-1-idv]
                if bag_types[len(sub_terms)-1-idv]==BagType.TRANSITIONAL:
                    delim=i.find('(')
                    delim2=i.find('_')
                    CHILDREN_BACKTRACE+=indent_marginalization+"    backtrace_"+i[delim2+1:delim+1]+"hashTable,tmp"+str(idv)+","+i[delim+1:].replace("hashTable,","")+";\n"
                else:
                    delim=i.find('(')
                    delim2=i.find('_')
                    CHILDREN_BACKTRACE+=indent_marginalization+"    backtrace_"+i[delim2+1:delim+1]+"hashTable,tmp"+str(idv)+","+i[delim+1:].replace("hashTable,","")+";\n"+indent_marginalization+"    bracket+=1;\n"     
            CHILDREN_SUM=print_children(sub_terms,0,indent_marginalization)
            CHILDREN_MAX=print_nested(len(sub_terms)-1)

        # special case: root
        skip_index_parts = False
        if len(equation.indices)==0:
            skip_index_parts = True
        index_part = False
        tmp=[]
        first=True
        for line in open('resources/memo/template_transitional.c').readlines():
            
            if line.find('index_start') >= 0:
                index_part = True
            if line.find('index_end') >= 0:
                index_part = False
            if index_part and skip_index_parts:
                continue
            if line.find('//') >= 0:
                continue

            line = line.replace('MAINNAME',equation.main_name)
            idx=equation.c_int_indices(tree_dec.ext_to_letter)
            line = line.replace('INT_INDICES',idx)
            line = line.replace('INDICES',
        ','.join([tree_dec.ext_to_letter[i] for i in equation.sorted_indices() if i in tree_dec.ext_to_letter]))
            line = line.replace('FOR_LOOP_NEW_VARIABLES_OPEN', marginalization)
            line = line.replace('CHILDREN_SUM', CHILDREN_SUM)
            line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
            line = line.replace('FOR_LOOP_NEW_VARIABLES_CLOSE',for_loop_closing)
            line = line.replace('INDENT', indent_marginalization)
            line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
            line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
            line = line.replace('SIZE_ARRAY',str(len([tree_dec.ext_to_letter[i] for i in equation.sorted_indices() if i in tree_dec.ext_to_letter])+1))
            if first:
                if skip_index_parts:
                    line=line.replace("hashTable,","hashTable")
                first=False
                declarations.append(line.replace("{",";"))
            tmp.append(line)
        functions.append("".join(tmp))
        tmp=[]
        first=True
        for line in open('resources/memo/backtrace/template_transitional.c').readlines():
            
            if line.find('index_start') >= 0:
                index_part = True
            if line.find('index_end') >= 0:
                index_part = False
            if index_part and skip_index_parts:
                continue
            if line.find('//') >= 0:
                continue

            line = line.replace('MAINNAME',equation.main_name)
            idx=equation.c_int_indices(tree_dec.ext_to_letter)
            line = line.replace('INT_INDICES',idx)
            line = line.replace('INDICES',
        ','.join([tree_dec.ext_to_letter[i] for i in equation.sorted_indices() if i in tree_dec.ext_to_letter]))
            line = line.replace('FOR_LOOP_NEW_VARIABLES_OPEN', marginalization)
            line = line.replace('CHILDREN_SUM', CHILDREN_SUM)
            line = line.replace('CHILDREN_BACKTRACE', CHILDREN_BACKTRACE)
            line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
            line = line.replace('FOR_LOOP_NEW_VARIABLES_CLOSE',for_loop_closing)
            line = line.replace('INDENT', indent_marginalization)
            line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
            line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
            line = line.replace('SIZE_ARRAY',str(len([tree_dec.ext_to_letter[i] for i in equation.sorted_indices() if i in tree_dec.ext_to_letter])+1))
            if first:
                first=False
                declarations.append(line.replace("{",";"))
            tmp.append(line)
            bg=len([tree_dec.ext_to_letter[i] for i in equation.sorted_indices() if i in tree_dec.ext_to_letter])
            if bg+1>bag_size[0]:
                #print("ICIIIII CHG BAGSIZE transitional",bg+1)
                bag_size[0]=bg+1
        backtrace.append("".join(tmp))
    elif tree_dec.bag_type[node]==BagType.DIAG_FIRST:
        #do diagonal things 
        equation = tree_dec.equations[node]
        if equation.inward is True:
            # const variables
            const = []
            for e in equation.constant_indices:
                if e in equation.variable_indices:
                    const.append(tree_dec.ext_to_letter[e]+"2")
                else:
                    const.append(tree_dec.ext_to_letter[e])
            
            # variable indices
            variables = [tree_dec.ext_to_letter[e] for e in sorted(equation.variable_indices,key=int)]

            # booleans for \neq sep
            terms = []

            for c in sorted(const):
                if int(reverse_ext_to_letter[variables[0]])<int(reverse_ext_to_letter[c]):
                    terms.append('('+variables[0]+'<'+c+')')
                else:
                    terms.append('('+variables[0]+'>'+c+'-1)')
                
            eq_some_const2 = " && ".join(terms)  

            terms = []


            for c in sorted(const):
                if int(reverse_ext_to_letter[variables[1]])<int(reverse_ext_to_letter[c]):
                    terms.append('('+variables[1]+'<'+c+')')
                else:
                    terms.append('('+variables[1]+'>'+c+'-1)')
                
            eq_some_const1 = " && ".join(terms)  

            # children sum
            CHILDREN_SUM = ''
            CHILDREN_MAX = '0'
            CHILDREN_BACKTRACE=''
            #print("equation_subtable",equation.subs_table)
            if len(tree_dec.bag_adj[equation.second_bag]) >= 2:
                sub_terms = []
                bag_types = []
                for sub_eq in sorted(tree_dec.equations[equation.second_bag].subterms,key= lambda x:children_depth[x]):
                    letter_table = {}
                    
                    for e in equation.absent_indices:
                        #print("abscent",e)
                        letter_table[e] = tree_dec.ext_to_letter[equation.subs_table[e]]
                    for e in equation.variable_indices:
                        #print("variable",e)
                        letter_table[e] = tree_dec.ext_to_letter[e]
                    #print("iciiiiiiii",letter_table,tree_dec.ext_to_letter)
                    sub_terms.append(sub_eq.c_code_print(letter_table, tree_dec.ext_to_letter).replace("(","( hashTable,").replace(","+variables[1],","+variables[1]+"+1"))
                    bag_types.append(tree_dec.bag_type[reverse_equation[sub_eq]])

                for idv in range(0,len(sub_terms)):
                    i=sub_terms[len(sub_terms)-1-idv]
                    if bag_types[len(sub_terms)-1-idv]==BagType.TRANSITIONAL:
                        delim=i.find('(')
                        delim2=i.find('_')
                        CHILDREN_BACKTRACE+="    backtrace_"+i[delim2+1:delim+1]+"hashTable,tmp"+str(idv)+","+i[delim+1:].replace("hashTable,","")+";\n"
                    else:
                        delim=i.find('(')
                        delim2=i.find('_')
                        CHILDREN_BACKTRACE+="    backtrace_"+i[delim2+1:delim+1]+"hashTable,tmp"+str(idv)+","+i[delim+1:].replace("hashTable,","")+";\n"+indent_marginalization+"    bracket+=1;\n"     
                CHILDREN_SUM=print_children_diag(sub_terms,0,"  ")
                CHILDREN_MAX=print_nested(len(sub_terms)-1)
            #print("heeeeeeerrrrreeeee",CHILDREN_BACKTRACE)
            first=False
            #print(equation.main_name,equation.sorted_indices())
            #print(equation.main_name,equation.sorted_indices())
            if bag_size[0]<len(const)+3:
                #print("ICIIIII CHG BAGSIZE",len(const)+3,const)
                bag_size[0]=len(const)+3
            with open('resources/memo/template_diag.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1_INC1_CONST_COMP', eq_some_const1)
                line = line.replace('V2_INC2_CONST_COMP', eq_some_const2)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('INC1', equation.increments[0])
                line = line.replace('INC2', equation.increments[1])
                line = line.replace('CONST_INT', ",".join(['int '+v for v in sorted(const)]))
                line = line.replace('CONST', ",".join(sorted(const)))
                line = line.replace('CHILDREN_SUM', CHILDREN_SUM)
                line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
                line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
                line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
                line = line.replace('SIZE_ARRAY',str(len(const)+3))
                
                functions.append(line)
                lp=line.split('\n')[0].replace("{",";\n")
                declarations.append(lp.replace("(","2("))
                declarations.append(lp)
            with open('resources/memo/backtrace/template_diag.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1_INC1_CONST_COMP', eq_some_const1)
                line = line.replace('V2_INC2_CONST_COMP', eq_some_const2)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('INC1', equation.increments[0])
                line = line.replace('INC2', equation.increments[1])
                line = line.replace('CONST_INT', ",".join(['int '+v for v in sorted(const)]))
                line = line.replace('CONST', ",".join(sorted(const)))
                line = line.replace('CHILDREN_SUM', CHILDREN_SUM.replace("min_value",""))
                line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
                line = line.replace('CHILDREN_BACKTRACE', CHILDREN_BACKTRACE)
                line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
                line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
                line = line.replace('SIZE_ARRAY',str(len(const)+3))
                
                backtrace.append(line)
                lp=line.split('\n')[0].replace("{",";\n")
                declarations.append(lp.replace("(","2("))
                declarations.append(lp)
            pass
        else:
            print('ici')
            # const variables
            const = []
            for e in equation.constant_indices:
                if e in equation.variable_indices:
                    const.append(tree_dec.ext_to_letter[e]+"2")
                else:
                    const.append(tree_dec.ext_to_letter[e])
            
            # variable indices
            variables = [tree_dec.ext_to_letter[e] for e in sorted(equation.variable_indices,key=int)]

            # booleans for \neq sep
            terms = []

            for c in sorted(const):
                if int(reverse_ext_to_letter[variables[0]])<int(reverse_ext_to_letter[c]):
                    terms.append('('+variables[0]+'<'+c+'+1)')
                else:
                    terms.append('('+variables[0]+'>'+c+')')
                
            eq_some_const2 = " && ".join(terms)  

            terms = []


            for c in sorted(const):
                if int(reverse_ext_to_letter[variables[1]])<int(reverse_ext_to_letter[c]):
                    terms.append('('+variables[1]+'<'+c+')')
                else:
                    terms.append('('+variables[1]+'>'+c+'-1)')
                
            eq_some_const1 = " && ".join(terms)  

            # children sum
            CHILDREN_SUM = ''
            CHILDREN_MAX = '0'
            CHILDREN_BACKTRACE=''
            #print("equation_subtable",equation.subs_table)
            if len(tree_dec.bag_adj[equation.second_bag]) >= 2:
                sub_terms = []
                bag_types = []
                for sub_eq in sorted(tree_dec.equations[equation.second_bag].subterms,key= lambda x:children_depth[x]):
                    letter_table = {}
                    
                    for e in equation.absent_indices:
                        #print("abscent",e)
                        letter_table[e] = tree_dec.ext_to_letter[equation.subs_table[e]]
                    for e in equation.variable_indices:
                        #print("variable",e)
                        letter_table[e] = tree_dec.ext_to_letter[e]
                    #print("iciiiiiiii",letter_table,tree_dec.ext_to_letter)
                    sub_terms.append(sub_eq.c_code_print(letter_table, tree_dec.ext_to_letter).replace("(","( hashTable,"))
                    bag_types.append(tree_dec.bag_type[reverse_equation[sub_eq]])

                for idv in range(0,len(sub_terms)):
                    i=sub_terms[len(sub_terms)-1-idv]
                    if bag_types[len(sub_terms)-1-idv]==BagType.TRANSITIONAL:
                        delim=i.find('(')
                        delim2=i.find('_')
                        CHILDREN_BACKTRACE+="    backtrace_"+i[delim2+1:delim+1]+"hashTable,tmp"+str(idv)+","+i[delim+1:].replace("hashTable,","")+";\n"
                    else:
                        delim=i.find('(')
                        delim2=i.find('_')
                        CHILDREN_BACKTRACE+="    backtrace_"+i[delim2+1:delim+1]+"hashTable,tmp"+str(idv)+","+i[delim+1:].replace("hashTable,","")+";\n"+indent_marginalization+"    bracket+=1;\n"     
                CHILDREN_SUM=print_children_diag(sub_terms,0,"  ")
                CHILDREN_MAX=print_nested(len(sub_terms)-1)
            #print("heeeeeeerrrrreeeee",CHILDREN_BACKTRACE)
            first=False
            #print(equation.main_name,equation.sorted_indices())
            #print(equation.main_name,equation.sorted_indices())
            if bag_size[0]<len(const)+3:
                #print("ICIIIII CHG BAGSIZE",len(const)+3,const)
                bag_size[0]=len(const)+3
            with open('resources/memo/template_diag_reverse.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1_INC1_CONST_COMP', eq_some_const1)
                line = line.replace('V2_INC2_CONST_COMP', eq_some_const2)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('INC1', equation.increments[0])
                line = line.replace('INC2', equation.increments[1])
                line = line.replace('CONST_INT', ",".join(['int '+v for v in sorted(const)]))
                line = line.replace('CONST', ",".join(sorted(const)))
                line = line.replace('CHILDREN_SUM', CHILDREN_SUM)
                line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
                line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
                line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
                line = line.replace('SIZE_ARRAY',str(len(const)+3))
                
                functions.append(line)
                lp=line.split('\n')[0].replace("{",";\n")
                declarations.append(lp.replace("(","2("))
                declarations.append(lp)
            with open('resources/memo/backtrace/template_diag_reverse.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1_INC1_CONST_COMP', eq_some_const1)
                line = line.replace('V2_INC2_CONST_COMP', eq_some_const2)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('INC1', equation.increments[0])
                line = line.replace('INC2', equation.increments[1])
                line = line.replace('CONST_INT', ",".join(['int '+v for v in sorted(const)]))
                line = line.replace('CONST', ",".join(sorted(const)))
                line = line.replace('CHILDREN_SUM', CHILDREN_SUM.replace("min_value",""))
                line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
                line = line.replace('CHILDREN_BACKTRACE', CHILDREN_BACKTRACE)
                line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
                line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
                line = line.replace('SIZE_ARRAY',str(len(const)+3))
                
                backtrace.append(line)
                lp=line.split('\n')[0].replace("{",";\n")
                declarations.append(lp.replace("(","2("))
                declarations.append(lp)
            pass
    return node

last=treat_bag('-1',tree_dec.root)


if clique:
    with open('resources/memo/clique_case.c','r') as f:
        functions.insert(0,f.read())
    declarations.insert(1,"int compute_CLIQUE(HashTable *hashTable, int i, int j, int k, int l);\n\nint compute_CLIQUE2(HashTable *hashTable, int i, int j, int k, int l);\n")
    declarations.insert(2,"void backtrace_CLIQUE(HashTable *hashTable, int score, int i,int j, int k, int l);\n\nvoid backtrace_CLIQUE2(HashTable *hashTable,int score,int i,int j,int k,int l);\n")
    with open('resources/memo/backtrace/clique_case.c','r') as f:
        backtrace.insert(0,f.read())

with open('resources/memo/header_template.c','r') as f:
    header=f.read()
#print("bag_size=",bag_size[0])
header=header.replace("##MODIFY_SIZE_TUPLE##",str(bag_size[0]))
#print("ici chg",bag_size[0])
with open('resources/memo/main_template.c','r') as f:
    main=f.read()
main=main.replace('ROOT',tree_dec.equations[tree_dec.root].main_name)
main=main.replace('END',sorted_letters[-1])
print("iciiii",sorted_letters,str(int(len(sorted_letters)/2)-1))
main=main.replace('ANCH',str(int(len(sorted_letters)/2)-1))

with open('../'+PSEUDO+'_better.c',"w") as f:
    f.write(header+"\n")
    f.write("\n".join(declarations)+"\n")
    f.write(main+"\n")
    f.write("\n".join(functions))
    f.write("\n".join(backtrace))
#print("ext_to_letter",tree_dec.ext_to_letter)
#print("all_extremities",tree_dec.all_extremities)

#print("done")
exit()
