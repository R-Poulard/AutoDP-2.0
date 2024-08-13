from autodp.equation_tree import TreeOfEquations, BagType
import math

import os,sys

#This is the start of the Boltzmann repartition extensions, not working yet

os.chdir("../auto-dp/")
if len(sys.argv)==1:
    
    PSEUDO="H"
    DIRECTORY="../test_exp/"
else:

    PSEUDO=sys.argv[1]  
    DIRECTORY="../Turner/"+sys.argv[2]

tree_dec = TreeOfEquations()

tree_dec.read_from_file('results/processed_td_files/processed_'+PSEUDO+'.td')

# helices: sets a lot of useful variables
tree_dec.set_helices(open('results/helix_annotations/'+PSEUDO+'.helix').readlines())

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

def depth(a,prev):
    if len(tree_dec.bag_adj[a])<=1:
        #print(a,"ended")
        return 1
    #print(a)
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
stack_cases=[]
#principle, declarations and function will be apppended during the tree parcours
#then joined by space and put in the file at once
clique=False
bag_size=[0] #length of the max tuple put as key into the hashing table

def print_children(lst,x,indent_marginalization):
    if len(lst) == 1:
        return indent_marginalization+"  double tmp"+str(x)+";\n  "+indent_marginalization+"if("+lst[0].replace("DATA","&tmp"+str(x))+"==0){possible=0;}\n"+indent_marginalization+"  else if(tmp"+str(x)+"==0.0){continue;}\n"
    else:
        return indent_marginalization+"  double tmp"+str(x)+";\n  "+indent_marginalization+"if("+lst[-1].replace("DATA","&tmp"+str(x))+"==0){possible=0;}\n"+indent_marginalization+"  else if(tmp"+str(x)+"==0.0){continue;}\n"+print_children(lst[:-1],x+1,indent_marginalization)

def print_children_backtrace(lst,x,indent_marginalization):
    if len(lst) == 1:
        return indent_marginalization+"  double tmp"+str(x)+";\n  "+indent_marginalization+lst[0].replace("DATA","&tmp"+str(x))+";\n"+indent_marginalization+"  if(tmp"+str(x)+"==0.0){continue;}\n"
    else:
        return indent_marginalization+"  double tmp"+str(x)+";\n  "+indent_marginalization+lst[-1].replace("DATA","&tmp"+str(x))+";\n"+indent_marginalization+"  if(tmp"+str(x)+"==0.0){continue;}\n"+print_children_backtrace(lst[:-1],x+1,indent_marginalization)

def print_children_diag(lst,x,indent_marginalization):
    if len(lst) == 1:
        return "double tmp"+str(x)+";\n   "+indent_marginalization+"if("+lst[0].replace("DATA","&tmp"+str(x))+"==0){possible=0;}\n"+indent_marginalization+"    else if(tmp"+str(x)+"==0.0){ goto loop;}\n"
    else:
        return "  double tmp"+str(x)+";\n   "+indent_marginalization+"if("+lst[-1].replace("DATA","&tmp"+str(x))+"==0){possible=0;}\n"+indent_marginalization+"    else if(tmp"+str(x)+"==0.0){ goto loop;}\n"+print_children_diag(lst[:-1],x+1,indent_marginalization)

def print_children_diag_backtrace(lst,x,indent_marginalization):
    if len(lst) == 1:
        return "double tmp"+str(x)+";\n   "+indent_marginalization+lst[0].replace("first,last,DATA","NULL,NULL,&tmp"+str(x))+";\n"+indent_marginalization+"    if(tmp"+str(x)+"==0.0){ goto loop;}\n"
    else:
        return "  double tmp"+str(x)+";\n   "+indent_marginalization+lst[-1].replace("first,last,DATA","NULL,NULL,&tmp"+str(x))+";\n"+indent_marginalization+"    if(tmp"+str(x)+"==0.0){ goto loop;}\n"+print_children_diag_backtrace(lst[:-1],x+1,indent_marginalization)

def print_nested(x):
    if x== 0:
        return "tmp"+str(x)
    else:
        return f"{print_nested(x-1)}*tmp{x}"
    
def print_nested_mfe(x):
    if x==-1:
        return "1.0"
    if x== 0:
        return "mfe0"
    else:
        return f"add({print_nested_mfe(x-1)},mfe{x})"
                       
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
        if clique is False:
            stack_cases.insert(0,tree_dec.equations[node].c_code_print_stack({},tree_dec.ext_to_letter))
        clique=True #no need to do more as the clique implementation will be added to the lists in the beginning once the tree done
        if bag_size[0]<5:
            #print("chgmt in clique")
            bag_size[0]=5
        pass
    elif tree_dec.bag_type[node]==BagType.TRANSITIONAL:
         # iteration over new variables
        
        equation = tree_dec.equations[node]
        stack_cases.append(equation.c_code_print_stack({},tree_dec.ext_to_letter))
        indent="    "
        marginalization = ""
        marginalization_indices = sorted(list(equation.marginalization),key=int)
        introduced_marginalized = []
        #print("marginalization", marginalization_indices)
        MFE_SUM=[]
        MFE_BACKTRACE=[]
        mfe_itg=0
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
            if ord(new_letter)==ord(minimum)+1 and sorted_letters.index(new_letter)%2==0:
                MFE_SUM.append("      double mfe"+str(mfe_itg)+" = MFEFree(pk,"+minimum+","+new_letter+"-1);\n")                
                MFE_BACKTRACE.append("        backtrace_MFEFree(pk,bt,mfe"+str(mfe_itg)+","+minimum+","+new_letter+"-1);\n")
                mfe_itg+=1
            elif ord(new_letter)==ord(maximum)-1 and sorted_letters.index(new_letter)%2==1:
                MFE_SUM.append("      double mfe"+str(mfe_itg)+" = MFEFree(pk,"+new_letter+","+maximum+"-1);\n")
                MFE_BACKTRACE.append("        backtrace_MFEFree(pk,bt,mfe"+str(mfe_itg)+","+new_letter+","+maximum+"-1);\n")
                mfe_itg+=1
        ##print("marginalization",repr(marginalization))
        marginalization = marginalization.rstrip('\n')
        indent_marginalization = len(marginalization_indices)*indent
        MFE_SUM=indent_marginalization.join(MFE_SUM)
        MFE_BACKTRACE=indent_marginalization.join(MFE_BACKTRACE)
        MFE_MAX=print_nested_mfe(mfe_itg-1)

        # closing the for loop
        for_loop_closing = ""
        for k in range(len(marginalization_indices)-1,-1,-1):
            for_loop_closing += indent*(k+1)+'}\n'

        # children sum
        CHILDREN_SUM = ''
        CHILDREN_BACKTRACE=''
        CONDITIONS=[]
        if len(tree_dec.bag_adj[node]) >= 2:
            sub_terms = []
            bag_types=[]
            for sub_eq in sorted(equation.subterms,key= lambda x:children_depth[x]):
                #print("iciiiiiiiiiiiiiii",reverse_equation[sub_eq],sub_eq,tree_dec.bag_type[reverse_equation[sub_eq]])
                if tree_dec.bag_type[reverse_equation[sub_eq]]==BagType.TRANSITIONAL:
                    sub_terms.append(sub_eq.c_code_print_get({},tree_dec.ext_to_letter).replace("(","( pk,first,last,DATA,"))
                else:
                    sub_terms.append(sub_eq.c_code_print_get({},tree_dec.ext_to_letter).replace("(","( pk,first,last,DATA,"))     
                    if tree_dec.bag_type[reverse_equation[sub_eq]]==BagType.DIAG_SECOND:
                        if sub_eq.inward==False:
                            CONDITIONS.append("!evaluate(pk,"+tree_dec.ext_to_letter[sub_eq.sorted_indices()[0]]+"-1,"+tree_dec.ext_to_letter[sub_eq.sorted_indices()[1]]+")")
                        else:
                            CONDITIONS.append("!evaluate(pk,"+tree_dec.ext_to_letter[sub_eq.sorted_indices()[0]]+","+tree_dec.ext_to_letter[sub_eq.sorted_indices()[1]]+"-1)")
                    else:
                        CONDITIONS.append("!evaluate(pk,"+tree_dec.ext_to_letter[sub_eq.sorted_indices()[0]]+","+tree_dec.ext_to_letter[sub_eq.sorted_indices()[3]]+"-1)")
                        CONDITIONS.append("!evaluate(pk,"+tree_dec.ext_to_letter[sub_eq.sorted_indices()[1]]+"-1,"+tree_dec.ext_to_letter[sub_eq.sorted_indices()[2]]+")")
                   
                bag_types.append(reverse_equation[sub_eq])
            if len(CONDITIONS)!=0:
                CONDITIONS="||".join(CONDITIONS)
                CONDITIONS=indent_marginalization+"  if("+CONDITIONS+"){continue;}"
            else:
                CONDITIONS=""
           #print("condition",CONDITIONS)
            
            for idv in range(0,len(sub_terms)):
                i=sub_terms[len(sub_terms)-1-idv]
                if tree_dec.bag_type[bag_types[len(sub_terms)-1-idv]]==BagType.TRANSITIONAL:
                    delim=i.find('(')
                    delim2=i.find('_')
                    CHILDREN_BACKTRACE+=indent_marginalization+"    int res"+str(idv)+"="+i.replace("get","backtrace").replace("first,last,DATA","bt,tmp"+str(idv))+";\n"
                else:
                   #print("nooo",tree_dec.bag_type[bag_types[len(sub_terms)-1-idv]],len(tree_dec.bag_adj[bag_types[len(sub_terms)-1-idv]]))
                    if tree_dec.bag_type[bag_types[len(sub_terms)-1-idv]]==BagType.DIAG_SECOND and len(tree_dec.bag_adj[bag_types[len(sub_terms)-1-idv]])>=2:
                        delim=i.find('(')
                        delim2=i.find('_')
                        CHILDREN_BACKTRACE+=indent_marginalization+"    int res"+str(idv)+"="+i.replace("get","backtrace").replace("first,last,DATA","bt,tmp"+str(idv))+";\n"
                    else:
                        delim=i.find('(')
                        delim2=i.find('_')
                        CHILDREN_BACKTRACE+=indent_marginalization+"    int res"+str(idv)+"="+i.replace("get","backtrace").replace("first,last,DATA","bt,tmp"+str(idv))+";\n"+indent_marginalization+"    bt->bracket+=1;\n"     
            CHILDREN_SUM=print_children(sub_terms,0,indent_marginalization)
            CHILDREN_SUM_BT=print_children_backtrace(sub_terms,0,indent_marginalization)
            CHILDREN_MAX=print_nested(len(sub_terms)-1)

        # special case: root
        skip_index_parts = False
        if len(equation.indices)==0:
            skip_index_parts = True
        index_part = False
        tmp=[]
        first=True
        for line in open('resources/stack_Turner_exp/template_transitional.c').readlines():
            
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
            line = line.replace('MFE_SUM', MFE_SUM)
            line = line.replace('MFE_MAX', MFE_MAX)
            line = line.replace('CONDITIONS', CONDITIONS)
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
                declarations.append(line.replace("{",";").replace('0','1').replace("get_","compute_").replace(",double * value",""))
                
            tmp.append(line)
        functions.append("".join(tmp))
        tmp=[]
        first=True
        for line in open('resources/stack_Turner_exp/backtrace/template_transitional.c').readlines():
            
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
            line = line.replace('CHILDREN_SUM', CHILDREN_SUM_BT.replace("first,last","NULL,NULL"))
            line = line.replace('CHILDREN_BACKTRACE', CHILDREN_BACKTRACE)
            line = line.replace('RETURN_CHECK', " && ".join(["res"+str(i) for i in range(0,len(sub_terms))]))
            line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
            line = line.replace('MFE_SUM', MFE_SUM)
            line = line.replace('MFE_MAX', MFE_MAX)
            line = line.replace('MFE_BACKTRACE', MFE_BACKTRACE)  
            line = line.replace('CONDITIONS', CONDITIONS)         
            line = line.replace('FOR_LOOP_NEW_VARIABLES_CLOSE',for_loop_closing)
            line = line.replace('INDENT', indent_marginalization)
            line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
            line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
            line = line.replace('SIZE_ARRAY',str(len([tree_dec.ext_to_letter[i] for i in equation.sorted_indices() if i in tree_dec.ext_to_letter])+1))
            if first:
                first=False
                declarations.append(line.replace("{",";"))
               #print(declarations)
            tmp.append(line)
            bg=len([tree_dec.ext_to_letter[i] for i in equation.sorted_indices() if i in tree_dec.ext_to_letter])
            if bg+1>bag_size[0]:
                #print("ICIIIII CHG BAGSIZE transitional",bg+1)
                bag_size[0]=bg+1
        backtrace.append("".join(tmp))
    elif tree_dec.bag_type[node]==BagType.DIAG_FIRST:
        #do diagonal things 
       
        equation = tree_dec.equations[node]
        stack_cases.append(equation.c_code_print_stack({},tree_dec.ext_to_letter))
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

            terms = []

            BOUNDI="INF"
            BOUNDI_value=99999999
            BOUNDJ_value=-1
            BOUNDJ="INT_MIN"
            nb_anch=0
            for c in reverse_ext_to_letter.keys():
                if int(reverse_ext_to_letter[variables[1]])>int(reverse_ext_to_letter[c]) and int(reverse_ext_to_letter[variables[0]])<int(reverse_ext_to_letter[c]):
                    nb_anch+=1
                    if c not in sorted(const):
                        continue
                    if int(reverse_ext_to_letter[c])>BOUNDJ_value:
                        BOUNDJ_value=int(reverse_ext_to_letter[c])
                        BOUNDJ=c
                    if int(reverse_ext_to_letter[c])<BOUNDI_value:
                        BOUNDI_value=int(reverse_ext_to_letter[c])
                        BOUNDI=c
           #print("res ",BOUNDI,BOUNDJ)
            nb_anch=int((nb_anch-2)/2+1)
           #print(nb_anch)
            # children sum
            add_before=""
            add_after=""
            ADJJ=''
            if BOUNDJ!='INT_MIN':
                cmp=0
                idx_boundJ=sorted_letters.index(BOUNDJ)
                idx_V2=sorted_letters.index(variables[1])
                #print("ici",BOUNDJ,variables[1],idx_V2,idx_boundJ)
                for i in range(idx_V2-1,-1,-1):
                    if i==idx_boundJ:      
                        break
                    if i%2==1:
                        cmp+=1
                if cmp!=0:
                    add_before="+"+str(cmp)
                    ADJJ=add_before
            ADJI=''
            if BOUNDI!='INF':
                cmp=0
                idx_boundI=sorted_letters.index(BOUNDI)
                idx_V1=sorted_letters.index(variables[0])
               #print("ici",BOUNDI,variables[0],idx_V1,idx_boundI)
                for i in range(idx_V1+idx_V1%2,idx_boundI):
                   #print("I=",i)
                    if i==idx_V1+idx_V1%2:
                        continue
                    if i%2==1:
                        cmp+=1
                    if i==idx_boundI:
                        break
                    
                if cmp!=0:
                    add_after="-"+str(cmp)
                    ADJI=add_after
            
           #print("res ",BOUNDI,ADJI,BOUNDJ,ADJJ)
            CONST_SUM = ''
            CONST_MAX = '0'
            CONST_BACKTRACE = ""
            nbres=0
            if BOUNDI!='INF' and ord(BOUNDI)==ord(variables[0])+2:
                CONST_SUM="double mfe1=MFEFree(pk,"+variables[0]+","+BOUNDI+"-1);"
                CONST_MAX="mfe1"
                CONST_BACKTRACE="int res"+str(nbres)+"=backtrace_MFEFree(pk,bt,mfe1,"+variables[0]+","+BOUNDI+"-1);\n"
                nbres+=1
            if BOUNDJ!='INT_MIN' and ord(BOUNDJ)==ord(variables[1])-2:
                CONST_SUM+="\n    double mfe2=MFEFree(pk,"+BOUNDJ+","+variables[1]+");"
                CONST_MAX="add("+CONST_MAX+",mfe2)"
                CONST_BACKTRACE+="        int res"+str(nbres)+"=backtrace_MFEFree(pk,bt,mfe2,"+BOUNDJ+","+variables[1]+");\n"
                nbres+=1
            CHILDREN_SUM = ''
            CHILDREN_SUM_BT = ''
            CHILDREN_MAX = '1.0'
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
                    sub_terms.append(sub_eq.c_code_print_get(letter_table, tree_dec.ext_to_letter).replace(","+variables[1],","+variables[1]+"+1").replace("(","( pk,first,last,DATA,"))
                    bag_types.append(tree_dec.bag_type[reverse_equation[sub_eq]])

                for idv in range(0,len(sub_terms)):
                    i=sub_terms[len(sub_terms)-1-idv]
                    delim=i.find('(')
                    delim2=i.find('_')
                    if len(CHILDREN_BACKTRACE)==0:
                        CHILDREN_BACKTRACE+="       bt->bracket+=1;\n"+"int res"+str(nbres)+"="+i.replace("get","backtrace").replace("first,last,DATA","bt,tmp"+str(idv))+";\n"     
                    else:
                        CHILDREN_BACKTRACE+="      backtrace_"+"int res"+str(nbres)+"="+i.replace("get","backtrace").replace("first,last,DATA","bt,tmp"+str(idv))+";\n"        
                    nbres+=1
                CHILDREN_SUM=print_children_diag(sub_terms,0," ")
                CHILDREN_SUM_BT=print_children_diag_backtrace(sub_terms,0," ")
                CHILDREN_MAX=print_nested(len(sub_terms)-1)
            #print("heeeeeeerrrrreeeee",CHILDREN_BACKTRACE)
            first=False
            #print(equation.main_name,equation.sorted_indices())
            #print(equation.main_name,equation.sorted_indices())
            if bag_size[0]<len(const)+4:
                #print("ICIIIII CHG BAGSIZE",len(const)+3,const)
                bag_size[0]=len(const)+4
            with open('resources/stack_Turner_exp/template_diag.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('BOUNDI',BOUNDI)
                line = line.replace('BOUNDJ',BOUNDJ)
                line = line.replace('CONST_SUM',CONST_SUM)
                line = line.replace('CONST_MAX',CONST_MAX)
                line = line.replace('ADJJ',ADJJ)
                line = line.replace('ADJI',ADJI)
                line = line.replace('NB_HELIX',str(nb_anch))
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
                declarations.append(lp)
                declarations.append(lp.replace('0','1'))
                declarations.append(lp.replace('0','1').replace("get_","compute_").replace(",double *value",""))
                
            with open('resources/stack_Turner_exp/backtrace/template_diag.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('BOUNDI',BOUNDI)
                line = line.replace('BOUNDJ',BOUNDJ)
                line = line.replace('CONST_SUM',CONST_SUM)
                line = line.replace('CONST_MAX',CONST_MAX)
                line = line.replace('CONST_BACKTRACE',CONST_BACKTRACE)
                line = line.replace('RETURN_CHECK', " && ".join(["res"+str(i) for i in range(0,nbres)]))
                line = line.replace('ADJJ',ADJJ)
                line = line.replace('ADJI',ADJI)
                line = line.replace('NB_HELIX',str(nb_anch))
                line = line.replace('INC1', equation.increments[0])
                line = line.replace('INC2', equation.increments[1])
                line = line.replace('CONST_INT', ",".join(['int '+v for v in sorted(const)]))
                line = line.replace('CONST', ",".join(sorted(const)))
                line = line.replace('CHILDREN_SUM', CHILDREN_SUM_BT.replace("min_value",""))
                line = line.replace('CHILDREN_MAX', CHILDREN_MAX)
                line = line.replace('CHILDREN_BACKTRACE', CHILDREN_BACKTRACE)
                line = line.replace('INT_VALUE_NAME',str(ord(equation.main_name)))
                line = line.replace('INT_VALUE_SECOND_NAME',str(-1*ord(equation.main_name)))
                line = line.replace('SIZE_ARRAY',str(len(const)+3))
                
                backtrace.append(line)
                lp=line.split('\n')[0].replace("{",";\n")
                declarations.append(lp)
                declarations.append(lp.replace('0','1'))
                
            pass
        else:

            const = []
            for e in equation.constant_indices:
                if e in equation.variable_indices:
                    const.append(tree_dec.ext_to_letter[e]+"2")
                else:
                    const.append(tree_dec.ext_to_letter[e])
            
            # variable indices
            variables = [tree_dec.ext_to_letter[e] for e in sorted(equation.variable_indices,key=int)]

            terms = []

            BOUNDI="0"
            BOUNDI_value=-1
            BOUNDJ_value=99999
            BOUNDJ="MAX-1"
            nb_anch=0
            for c in reverse_ext_to_letter.keys():
                if int(reverse_ext_to_letter[variables[1]])<int(reverse_ext_to_letter[c]) and BOUNDJ_value>int(reverse_ext_to_letter[c]):
                    if c not in sorted(const):
                        continue
                    
                    BOUNDJ_value=int(reverse_ext_to_letter[c])
                    BOUNDJ=c
                elif int(reverse_ext_to_letter[variables[0]])>int(reverse_ext_to_letter[c]) and BOUNDI_value<int(reverse_ext_to_letter[c]):
                    if c not in sorted(const):
                        continue
                    BOUNDI_value=int(reverse_ext_to_letter[c])
                    BOUNDI=c
           #print("res ",BOUNDI,BOUNDJ)
           #print(nb_anch)
            # children sum
            add_before=""
            add_after=""
            ADJI=''
            if BOUNDI!='0':
                cmp=0
                idx_boundI=sorted_letters.index(BOUNDI)
                idx_V1=sorted_letters.index(variables[0])
                #print("ici",BOUNDJ,variables[1],idx_V2,idx_boundJ)
                for i in range(idx_V1-1,-1,-1):
                    if i==idx_boundI:      
                        break
                    if i%2==1:
                        cmp+=1
                if cmp!=0:
                    add_before="+"+str(cmp)
                    ADJI=add_before
            ADJJ=''
            if BOUNDJ!='MAX-1':
                cmp=0
                idx_boundJ=sorted_letters.index(BOUNDJ)
                idx_V2=sorted_letters.index(variables[1])
               #print("ici",BOUNDJ,variables[1],idx_V2,idx_boundJ)
                for i in range(idx_V2+idx_V2%2,idx_boundJ):
                   #print("I=",i)
                    if i==idx_V2+idx_V2%2:
                        continue
                    if i%2==1:
                        cmp+=1
                    if i==idx_boundJ:
                        break
                    
                if cmp!=0:
                    add_after="-"+str(cmp)
                    ADJJ=add_after
            
           #print("res ",BOUNDI,ADJI,BOUNDJ,ADJJ)
            CONST_SUM = ''
            CONST_MAX = '0'
            CONST_BACKTRACE = ""
            #if BOUNDI!='MAX-1' and ord(BOUNDI)==ord(variables[0])+2:
            CONST_SUM="double mfe1=MFEFree(pk,"+variables[1]+","+BOUNDJ+");"
            CONST_MAX="mfe1"
            CONST_BACKTRACE="backtrace_MFEFree(pk,bt,mfe1,"+variables[1]+","+BOUNDJ+");\n"
            #if BOUNDJ!='0' and ord(BOUNDI)==ord(variables[0])-2:
            CONST_SUM+="\n    double mfe2=MFEFree(pk,"+BOUNDI+","+variables[0]+");"
            CONST_MAX="add("+CONST_MAX+",mfe2)"
            CONST_BACKTRACE+="        backtrace_MFEFree(pk,bt,mfe2,"+BOUNDI+","+variables[0]+");\n"
            CHILDREN_SUM = ''
            CHILDREN_MAX = '1.0'
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
                    sub_terms.append(sub_eq.c_code_print(letter_table, tree_dec.ext_to_letter).replace("(","( pk,").replace(","+variables[0],","+variables[0]+"-1"))
                    bag_types.append(tree_dec.bag_type[reverse_equation[sub_eq]])

                for idv in range(0,len(sub_terms)):
                    i=sub_terms[len(sub_terms)-1-idv]
                    i=sub_terms[len(sub_terms)-1-idv]
                    delim=i.find('(')
                    delim2=i.find('_')
                    if len(CHILDREN_BACKTRACE)==0:
                        CHILDREN_BACKTRACE+="       bt->bracket+=1;\n       backtrace_"+i[delim2+1:delim+1]+"pk,bt,tmp"+str(idv)+","+i[delim+1:].replace("pk,","")+";\n"     
                    else:
                        CHILDREN_BACKTRACE+="      backtrace_"+i[delim2+1:delim+1]+"pk,bt,tmp"+str(idv)+","+i[delim+1:].replace("pk,","")+";\n"         
                CHILDREN_SUM=print_children_diag(sub_terms,0,"      ")
                CHILDREN_MAX=print_nested(len(sub_terms)-1)
            #print("heeeeeeerrrrreeeee",CHILDREN_BACKTRACE)
            first=False
            #print(equation.main_name,equation.sorted_indices())
            #print(equation.main_name,equation.sorted_indices())
            if bag_size[0]<len(const)+4:
                #print("ICIIIII CHG BAGSIZE",len(const)+3,const)
                bag_size[0]=len(const)+4
            with open('resources/stack_Turner_exp/template_diag_reverse.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('BOUNDI',BOUNDI)
                line = line.replace('BOUNDJ',BOUNDJ)
                line = line.replace('CONST_SUM',CONST_SUM)
                line = line.replace('CONST_MAX',CONST_MAX)
                line = line.replace('ADJJ',ADJJ)
                line = line.replace('ADJI',ADJI)
                line = line.replace('NB_HELIX',str(nb_anch))
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
                declarations.append(lp)
                declarations.append(lp.replace("0(","1("))
                declarations.append(lp.replace("0(","2("))

            with open('resources/stack_Turner_exp/backtrace/template_diag_reverse.c','r') as f:
                line = f.read()
                line = line.replace('MAINNAME',equation.main_name)
                line = line.replace('V1',variables[0])
                line = line.replace('V2',variables[1])
                line = line.replace('BOUNDI',BOUNDI)
                line = line.replace('BOUNDJ',BOUNDJ)
                line = line.replace('CONST_SUM',CONST_SUM)
                line = line.replace('CONST_MAX',CONST_MAX)
                line = line.replace('CONST_BACKTRACE',CONST_BACKTRACE)
                line = line.replace('ADJJ',ADJJ)
                line = line.replace('ADJI',ADJI)
                line = line.replace('NB_HELIX',str(nb_anch))
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
                declarations.append(lp)
                declarations.append(lp.replace("0(","1("))
                declarations.append(lp.replace("0(","2("))

            pass
    return node

last=treat_bag('-1',tree_dec.root)


if clique:

    with open('resources/stack_Turner_exp/clique_case.c','r') as f:
        functions.insert(0,f.read())
    declarations.insert(1,"""PRIVATE int get_CLIQUE1(pk_compound *pk, Node **first, Node **last, double *value, int i, int i2, int j2, int j);

PRIVATE int get_CLIQUE0(pk_compound *pk, Node **first, Node **last, double *value, int i, int i2, int j2, int j);

PRIVATE int compute_CLIQUE1(pk_compound *pk,Node **first,Node** last, int i, int j, int k, int l);

PRIVATE int backtrace_CLIQUE0(pk_compound *pk,bt_struct* bt, double score, int i,int j, int k, int l);

PRIVATE int backtrace_CLIQUE1(pk_compound *pk,bt_struct* bt, double score, int i,int j, int k, int l);""")
    with open('resources/stack_Turner_exp/backtrace/clique_case.c','r') as f:
        backtrace.insert(0,f.read())

with open('resources/stack_Turner_exp/header_template.c','r') as f:
    header=f.read()
#print("bag_size=",bag_size[0])
header=header.replace("##MODIFY_SIZE_TUPLE##",str(bag_size[0]))
#print("ici chg",bag_size[0])
with open('resources/stack_Turner_exp/main_template.c','r') as f:
    main=f.read()
header=header.replace('ROOT',tree_dec.equations[tree_dec.root].main_name)
header=header.replace('PSEUDO',PSEUDO)
header=header.replace('END',sorted_letters[-1])
header=header.replace('ALL_ANC',",".join(sorted_letters))
header=header.replace("VARIABLE_PART","\n".join(stack_cases))
print("iciiii",sorted_letters,str(int(len(sorted_letters)/2)-1))
header=header.replace('ANCH',str(int(len(sorted_letters)/2)-1))
header=header.replace('##DECLARATIONS##',"\n".join(declarations))
header=header.replace('##NAME##',PSEUDO)

#f = open(snakemake.output[0],'w')
f = open(DIRECTORY+'/'+PSEUDO+'.c',"w")

print(header+"\n",file=f)
    #f.write("\n".join(declarations)+"\n")
print(main+"\n",file=f)
print("\n".join(functions),file=f)
print("\n".join(backtrace),file=f)
f.close()
"""
with open(DIRECTORY+'/'+PSEUDO+'.c',"w") as f:
    f.write(header+"\n")
    #f.write("\n".join(declarations)+"\n")
    f.write(main+"\n")
    f.write("\n".join(functions))
    f.write("\n".join(backtrace))
"""
#print("ext_to_letter",tree_dec.ext_to_letter)
#print("all_extremities",tree_dec.all_extremities)

#print("done")

f = open(DIRECTORY+"/"+PSEUDO+"_folding.h",'w')
print("#include <stdio.h>\n#include <stdlib.h>\n#include <stdbool.h>\n#include <math.h>\n#include <limits.h>\n#include <time.h>\n#include <string.h>\n#include <ViennaRNA/fold_compound.h>\n#include <ViennaRNA/utils/basic.h>\n#include <ViennaRNA/utils/strings.h>\n#include <ViennaRNA/mfe.h>\n#include <ViennaRNA/loop_energies.h>\n#include \"aux.h\"\n\nint setup_"+PSEUDO+"(vrna_fold_compound_t *fc,const char* seq, const char* con,bt_struct* BT);\n",file=f)
f.close()

if True or not os.path.exists(DIRECTORY+"aux.c"):
    f = open(DIRECTORY+"aux.c",'w')
    with open('resources/stack_Turner/aux.c','r') as f2:
        print("".join(f2.readlines()),file=f)
if True or not os.path.exists(DIRECTORY+"aux.h"):
    f = open(DIRECTORY+"aux.h",'w')
    with open('resources/stack_Turner/aux.h','r') as f2:
        print("".join(f2.readlines()),file=f)
exit()
