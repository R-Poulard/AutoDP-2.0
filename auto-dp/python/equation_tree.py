from autodp.tree_of_bags import TreeOfBags
from autodp.equation_classes import DiagCaseHelix, CliqueCaseHelix, TransitionalEquation
from enum import Enum

class BagType(Enum):
    TRANSITIONAL = 1
    DIAG_FIRST = 2
    DIAG_SECOND = 3
    CLIQUE = 4

class TreeOfEquations(TreeOfBags):

    def __init__(self):
        """
        For clarity, all attributes specific to TreeOfEquations
        (other attributes are inherited from TreeOfBags) are set here to
        empty/dummy values. No other attributes will be defined in the methods.
        """
        super().__init__()

        # will store for each helix label the set of bags representing the helix
        self.representation_bags = {}

        # for a bag, tells you for which helix it is part of the representation. 
        self.which_helix = {}

        # helices
        self.helices = []

        # all extremities: the set of all anchors, for practicality
        self.all_extremities = set([])

        # helix extremities: gives you the set of extremities associated to an helix
        self.helix_extremities = {}

        # const part: the anchors present in all bags representing an helix (for diag case mostly)
        self.const_part = {}

        # bag type
        self.bag_type = {}

        # equations
        self.equations = {}

        # extremity (position/integer) to variable name (latin letter)
        self.ext_to_letter = {}

    def set_ext_to_letter(self):
        for k, e in enumerate(sorted(list(self.all_extremities),key=lambda x: int(x))):
            self.ext_to_letter[e] = chr(ord('a')+k)

    def set_dp_tables(self):

        cnt = 0 # alphabet/variable increment variable
        cnt_col = 0 # color increment variable (to match colors of other drawings)

        self.equations['-1'] = TransitionalEquation()

        for prev, u in self.dfs_edge_iterator():
            if not u[:1]=='H':
                # bag type
                self.bag_type[u] = BagType.TRANSITIONAL

                # equation object init
                equation = TransitionalEquation()
                equation.main_name = self.num_to_letters(cnt)
                equation.latex_name = self.num_to_letters(cnt)

                # indices and marginalizatin
                equation.indices = set(self.bag_content[u]).intersection(set(self.bag_content[prev]))
                equation.marginalization = set(self.bag_content[u]) - set(self.bag_content[prev])

                # connecting
                self.equations[u] = equation
                self.equations[prev].subterms.append(equation)

                # incrementing name but not color
                cnt += 1
            else:
                if set(self.helix_extremities[self.which_helix[u]]).issubset(set(self.bag_content[prev])):
                # clique
                    # set bag type
                    self.bag_type[u] = BagType.CLIQUE

                    # equation init
                    equation = CliqueCaseHelix()
                    equation.main_name = "CLIQUE"
                    equation.latex_name = "\\colorbox{c"+u.split('_')[0][1:]+"}{$C_{\\boxtimes}$}"

                    # indices
                    self.equations[u] = equation
                    equation.indices = self.helix_extremities[self.which_helix[u]]

                    # connecting
                    self.equations[prev].subterms.append(equation)

                    # updating color but not name
                    cnt_col += 1

                else:
                # diag: only setting the letter for the first bag out of the two representating of the helix.
                    if prev.split('_')[0]!=u.split('_')[0]:
                        # bag type
                        self.bag_type[u] = BagType.DIAG_FIRST

                        # equation init
                        equation = DiagCaseHelix()
                        equation.first_bag = u
                        equation.main_name = self.num_to_letters(cnt)
                        equation.latex_name = "\\colorbox{c"+u.split('_')[0][1:]+"}{$"+self.num_to_letters(cnt)+"$}"

                        # indices
                        absent_ex = (set(self.helix_extremities[self.which_helix[u]]) - set(self.bag_content[u])).pop()
                        sorted_exs = sorted(list(self.helix_extremities[self.which_helix[u]]), key=lambda x : int(x))
                        equation.variable_indices = set(self.helix_extremities[self.which_helix[u]]) - set([absent_ex, partner(absent_ex, sorted_exs)])
                        equation.constant_indices = self.const_part[self.which_helix[u]]
                        equation.absent_indices = set([absent_ex, partner(absent_ex, sorted_exs)])
                        
                        # increments
                        equation.increments = [increment(e, sorted_exs) for e in sorted(list(equation.variable_indices))]
                        equation.inward = not absent_ex in [min(self.helix_extremities[self.which_helix[u]]),max(self.helix_extremities[self.which_helix[u]])]

                        # subs table
                        for e in equation.absent_indices:
                            equation.subs_table[e] = subs(e, sorted_exs)

                        # connecting
                        self.equations[u] = equation
                        self.equations[prev].subterms.append(equation)

                        cnt += 1
                        cnt_col += 1
                    else:
                        # bag type
                        self.bag_type[u] = BagType.DIAG_SECOND
                        self.equations[prev].second_bag = u

                        self.equations[u] = self.equations[prev]

    def set_helices(self, helices):
        """
        sets a lot of fields: helices, helix_extremities, all_extremities, representation_bags
        and which_helix
        
        Input:
            - helices: lines of helix annotation files. typically open(snakemake.input.helix).readlines()
        """
        self.helices = helices

        for helixline in self.helices:

            # helix info extraction
            label = helixline.split(' ')[0]
            extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]
            
            # set helix_extremities
            self.helix_extremities[label] = extremities
            self.all_extremities = self.all_extremities.union(extremities)
           
            # set representation_bags and which helix
            self.representation_bags[label] = []
            print(label)
            for u in self.dfs_bag_iterator():
                print(u)
                if u.split('_')[0]==label:
                    self.representation_bags[label].append(u)
                    self.which_helix[u] = label
            # set const part
            print(self.representation_bags[label])
            self.const_part[label] =set.intersection(*[set(self.bag_content[u]) for u in self.representation_bags[label]])

    @staticmethod
    def num_to_letters(cnt):
        if cnt < 26:
            return chr(ord('a') + cnt).upper() #'\\colorbox{c'+str(cnt)+'}{'+chr(ord('a') + cnt).upper()+'}'
        else:
            return chr(ord('a') + int(cnt/26)).upper()+chr(ord('a') + int(cnt%26)).upper()
         

    def contract_to_skeleton(self):
        """
        Contracts each helix representation to 2 bags in the diag
        case and one bag in the clique case.
        """

        for helixline in self.helices:
            label = helixline.split(' ')[0]
            extremities = [c.replace(' ','') for c in helixline.split('(')[1].split(')')[0].split(',')]

            queue = [('-1',self.root)]
            while len(queue) > 0:
                prev,u = queue.pop()

                if u.split('_')[0]==label:
                    if not set(extremities).issubset(set(self.bag_content[prev])):
                    # diag case
                        keep_going = True
                        while keep_going:
                            keep_going = False
                            for v in self.bag_adj[u]:
                                if v.split('_')[0]==label:
                                    for w in self.bag_adj[v]:
                                        if w!=u and w.split('_')[0] == label:
                                            # grand-child is still in helix, contracting v into u
                                            self.bag_adj[u] = [bag for bag in self.bag_adj[u] if bag!=v]+[w]
                                            self.bag_adj[w] = [bag for bag in self.bag_adj[w] if bag!=v]+[u]
                                            keep_going = True

                    else: 
                    # clique case    
                        self.bag_adj[u] = []



                for v in self.bag_adj[u]:
                    if v!=prev:
                        queue.append((u,v))

    def filter_anchors(self):
        """
        Processes the bags so that they contain anchors (helix extremities only)
        """
        for u in self.dfs_bag_iterator():
            self.bag_content[u] = [vert for vert in self.bag_content[u] if vert in self.all_extremities]
        

def increment(e, sorted_exs):
    """
    depending on wheter internal
    or external extremal base pair,
    the increment in the dp equations
    that update variables are not the
    same.
    """
    if e==sorted_exs[0]:
        return '+1'
    if e==sorted_exs[1]:
        return '-1'
    if e==sorted_exs[2]:
        return '+1'
    if e==sorted_exs[3]:
        return '-1'

def partner(e, sorted_ext):
    """
    in the standard (i,ip,jp,j) representation
    of helices that we use, 
    with i < ip < jp < j, i.e.
    i,j the outer arc and ip,jp the inner arc,
    i and j are partners, 
    and ip, jp as well.
    """
    if e==sorted_ext[0]:
        return sorted_ext[3]
    if e==sorted_ext[1]:
        return sorted_ext[2]
    if e==sorted_ext[3]:
        return sorted_ext[0]
    if e==sorted_ext[2]:
        return sorted_ext[1]

def subs(e, sorted_exs):
    """
    In the diag case, i is incrementally
    transformed into ip, or vice-versa.
    same for j and jp.
    subs: i -> ip, ip -> o, jp ->j, j->jp
    """
    if e==sorted_exs[0]:
        return sorted_exs[1]
    if e==sorted_exs[1]:
        return sorted_exs[0]
    if e==sorted_exs[2]:
        return sorted_exs[3]
    if e==sorted_exs[3]:
        return sorted_exs[2]
