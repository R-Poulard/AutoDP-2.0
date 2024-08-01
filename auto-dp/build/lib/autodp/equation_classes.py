def ext_to_letter(k):
    """
    Simply returns the k-th letter of the alphabet.
    """
    return chr(ord('a')+int(k))

class CommonEquationFeatures():
    """
    Basic equation class from which other equation object
    classes will inherit. 

    Its main purpose is to require the implementation 
    of some methods (for now, just latex_print, but
    in the future, c_code_print maybe ?)
    """

    # class specific methods
    def terminal_print(self, ext_to_letter):
        # request implementation of terminal print
        raise NotImplementedError('terminal_print must be implemented for class '+self.__class__.__name__)
        
    def latex_print(self, letter_table, ext_to_letter):
        # request implementation of latex_print
        raise NotImplementedError('latex_print must be implemented for class '+self.__class__.__name__)
    
    def c_allocation_print(self):
        # request implementation of c_allocation_print
        raise NotImplementedError('c_allocation_print must be implemented for class '+self.__class__.__name__)

    def c_free_print(self):
        # request implementation of c_free_print
        raise NotImplementedError('c_free_print must be implemented for class '+self.__class__.__name__)
    
    def c_code_print(self, letter_table, ext_to_letter):
        # request implementation of c_code
        raise NotImplementedError('c_code_print must be implemented for class '+self.__class__.__name__)

    def sorted_indices(self, letter_table, ext_to_letter):
        # request implementation of sorted_indices
        raise NotImplementedError('sorted_indices must be implemented for class '+self.__class__.__name__)

    # common methods
    
    def c_index_function(self, ext_to_letter):
        indices_raw = [ext_to_letter[v] for v in self.sorted_indices()]
        indices = []
        for k, i in enumerate(indices_raw):
            if i in indices_raw[:k]:
                indices.append(i+'2')
            else:
                indices.append(i)
        indent = "    "
        res = "".join(["int index_",self.main_name])
        res += "("+",".join(['int '+i for i in indices])+')  {\n'
        res += indent+'return '
        terms = []
        for k, index in enumerate(indices):
            terms.append('*'.join(['n' for _ in range(len(indices)-1-k)]+[index]))
        res += '+'.join(terms)+';\n'
        res += '}\n'
        return res

    def c_index_fun_signature(self, ext_to_letter):
        indices_raw = [ext_to_letter[v] for v in self.sorted_indices()]
        indices = []
        for k, i in enumerate(indices_raw):
            if i in indices_raw[:k]:
                indices.append(i+'2')
            else:
                indices.append(i)
        indent = "    "
        res = "".join(["int index_",self.main_name])
        res += "("+",".join(['int '+i for i in indices])+');'
        return res

    def c_int_indices(self, ext_to_letter):
        indices = self.sorted_indices()
        indices = [ext_to_letter[i] for i in indices]
        return ','.join(['int '+i for i in indices])

    
    def c_array_init_fill(self, ext_to_letter):
        indent = "    "
        res = "".join(["void init_fill_",self.main_name,"() {\n"])
        indices = self.sorted_indices()
        indices = [ext_to_letter[i] for i in indices]
        for k,index in enumerate(indices):
            res += (k+1)*indent
            if k==0:
                res += "".join(["for (int ",index,"=0;",index,"<n;",index,"++) {\n"])
            else:
                res += "".join(["for (int ",index,"=",indices[k-1],";",index,"<n;",index,"++) {\n"])
        res+= (len(indices)+1)*indent+self.main_name+'[index_'+self.main_name+'('+",".join(indices)+')] = INT_MIN;\n' 
        for k in range(len(indices),-1,-1):
            res += k*indent+'}\n'

        return res


class TransitionalEquation(CommonEquationFeatures):
    """
    Class for representing 
    DP equations linked to ``transitional''
    bags (i.e. bags that are not part of a clique
    representation)

    Its fields contain all of the information 
    needed to reconstruct the equation. They
    are initialized to dummy values, and
    will be updated based on what read 
    in canonical tree decompositions.

    As any class inheriting from CommonEquationFeatures,
    it needs to implement latex_print
    """

    def __init__(self):

        # type of equation
        self.type = "TRANSITIONAL"
        
        # name of the main table of the equation
        self.main_name = ""
        self.latex_name = ""

        # table indices
        self.indices = set([])

        # new variable (marginalization)
        self.marginalization = set([])

        # sub-terms: list of other equations
        self.subterms = []

    def terminal_print(self, ext_to_letter):
        print(self.main_name, self.latex_name, [ext_to_letter[i] for i in self.indices], self.subterms)
        print(self.main_name, self.latex_name, self.indices, self.subterms)

    def c_allocation_print(self):
        total_size = '*'.join(['n' for _ in range(len(self.indices))])
        res = "double * "+self.main_name+" = malloc("+total_size+"*sizeof(double));"
        return res

    def c_free_print(self):
        return "".join(["free(",self.main_name,");"])

    def sorted_indices(self):
        return sorted(list(self.indices),key=lambda x: int(x))

    def c_code_print(self, letter_table, ext_to_letter):
        for e in self.indices:
            if e not in letter_table.keys():
                letter_table[e] = ext_to_letter[e]
        return "compute_"+self.main_name+'('+','.join([letter_table[e] for e in self.indices])+')'

    def latex_print(self, letter_table, ext_to_letter):

        for e in self.indices:
            if e not in letter_table.keys():
                letter_table[e] = ext_to_letter[e]

        res = ""
        res += self.latex_name+'['
        res += ','.join([letter_table[e] for e in self.indices])
        res += ']'
        return res

class CliqueCaseHelix(CommonEquationFeatures):

    def __init__(self):

        # tupe
        self.type = "CLIQUE"

        # table indices
        self.indices = []
        
        # name of the main table of the equation
        self.main_name = ""
        self.latex_name = ""

    def terminal_print(self, ext_to_letter):
        print(self.main_name, self.latex_name, [ext_to_letter[i] for i in self.indices])
    
    def c_allocation_print(self):
        total_size = '*'.join(['n' for _ in range(len(self.indices))])
        res = "double * "+self.main_name+" = malloc("+total_size+"*sizeof(double));"
        return res
        
    def c_free_print(self):
        return "".join(["free(",self.main_name,");"])
    
    def sorted_indices(self):
        return sorted(list(self.indices))

    def c_code_print(self, letter_table, ext_to_letter):
        for e in self.indices:
            if e not in letter_table.keys():
                letter_table[e] = ext_to_letter[e]
        mask=["","-1","","-1"]
        args=[letter_table[e]+i for i in mask for e in self.indices]
        for i in len(mask):
            args[i]=args[i]+mask[i]
        return "compute_CLIQUE2("+','.join(args)+')'
    
    def latex_print(self, letter_table, ext_to_letter):


        for e in self.indices:
            if e not in letter_table.keys():
                letter_table[e] = ext_to_letter[e]

        res = ""
        res += self.latex_name+'['
        res += ','.join([letter_table[e] for e in self.indices])
        res += ']'
        return res

class DiagCaseHelix(CommonEquationFeatures):

    def __init__(self):

        # type
        self.type = "DIAG"

        # name of the main table of the equation
        self.main_name = ""
        self.latex_name = ""

        # the ones before the \mid
        self.variable_indices = []

        # absent indices
        self.absent_indices = []

        # to connect the two
        self.corresponding_variable = {}

        # the ones after the \mid (the constraints, constant)
        self.constant_indices = []

        # whether the variable indices get a +1 (external bp) or -1 (internal)
        self.increments = []

        # the recursion terms - "children" equations
        self.subterms = []

        # boolean for in which sense it goes
        self.inward = True

        # substitution matrix for subterms
        self.subs_table = {}

        # diag helix canonical representation: two bags
        self.first_bag = '-1'
        self.second_bag = '-1'
        
    def terminal_print(self, ext_to_letter):
        print(self.main_name, self.latex_name, 
              [ext_to_letter[i] for i in self.variable_indices],
              [ext_to_letter[i] for i in self.constant_indices])
    
    def c_allocation_print(self):
        total_size = '*'.join(['n' for _ in range(len(self.variable_indices)+len(self.constant_indices))])
        res = "double * "+self.main_name+" = malloc("+total_size+"*sizeof(double));\n"
        res += "double * "+self.main_name+"2 = malloc("+total_size+"*sizeof(double));"
        return res
    
    def sorted_indices(self):
        return sorted(list(self.variable_indices), key=lambda x : int(x)) + sorted(list(self.constant_indices),key=lambda x: int(x))
    
    def c_free_print(self):
        return "".join(["free(",self.main_name,");"])
    
    def c_code_print(self, letter_table, ext_to_letter):
        for e in self.sorted_indices():
            if e not in letter_table.keys():
                letter_table[e] = ext_to_letter[e]
        return "compute_"+self.main_name+'('+','.join([letter_table[e] for e in self.sorted_indices()])+')'

    def latex_print(self, letter_table, ext_to_letter):

        for e in list(self.variable_indices)+list(self.constant_indices):
            if e not in letter_table.keys():
                letter_table[e] = ext_to_letter[e]

        res = self.latex_name+'['
        res += ','.join([letter_table[e] for e in self.variable_indices])
        res += '\mid '
        res += ','.join([letter_table[e] for e in self.constant_indices])
        res += ']'

        return res
