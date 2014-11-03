#-------------------------------------------------------------------------------------
#  Copyright 2014 Michael Peeri
#
#  This file is part of hmmdsl.
#  hmmdsl is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  hmmdsl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with hmmdsl.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------------
import hmmdsl_py


class _NodeNamer:
    def __init__(self, base="Node"):
        self.base = base
        self.id = 0
    def __call__(self):
        ret = "{}{}".format(self.base, self.id)
        self.id += 1
        return ret
_node_namer = _NodeNamer()

class Node:
    def __init__(self, descendents, model):
        self.model = model
        self.descendents = descendents
        self.props = {}
        self.name = _node_namer()
        self.model_up_to_date = False

    def __getitem__(self, key):
        return self._props[key]

    def __repr__(self):
        return "Node(descendents=[{}], props={}, model_up_to_date={})".format(",".join(map(lambda x: repr(x), self.elems())), self.props, self.model_up_to_date)

    def elems(self):
        # TODO: add descendents recursively
        elems = []
        for node in self.descendents:
            if isinstance(node,int):
                elems.append(node)
            else:
                if isinstance(node,Node):
                    elems.extend(node.elems())
                else:
                    raise "Unsupport type in Node.descendents"
            
        return elems

def write_node_edges_to_dot(node):
    ret = ""
    src = node.name
    for tgt in node.descendents:
        if isinstance(tgt, int):
            ret += "{} -> Seq{};\n".format(src, tgt)
        else:
            if isinstance(node,Node):
                ret += "{} -> {};\n".format(src, tgt.name)
                ret += write_node_edges_to_dot(tgt)
    return ret

def write_node_descriptions_to_dot(node):
    ret = ""
    ret += "{0} [label=\"{0}\" size={1}];\n".format(node.name, len(node.elems()))

    for tgt in node.descendents:
        if isinstance(tgt, int):
            ret += "Seq{} [size=1];\n".format(tgt)
        else:
            if isinstance(node,Node):
                ret += write_node_descriptions_to_dot(tgt)
    return ret

def write_nodes_to_dot(nodes):
    ret = "digraph G {\n"
    for n in nodes:
        ret += write_node_descriptions_to_dot(n)
    for n in nodes:
        ret += write_node_edges_to_dot(n)

    ret += "}\n"

    return ret
    
            

def make_square_matrix(N, val=0.0):
    rows = []
    for i in range(N):
        newrow = [val] * N
        rows.append(newrow)
    return rows
        
    
class UnionFind:
    # Union-Find for Kruskal impl. (see Kleinberg&Tardos section 4.6 p. 151)
    # TODO: Implement this efficiently, as K&T show
    def __init__(self, S):
        self._sets = []
        for s in S:
            self._sets.append( set([s]) )

    def find(self, u):
        for s in self._sets:
            if u in s:
                return s
        raise "%s not in any set in the union!".format(str(u))

    def union_old(self, A, B):
        N = len(self._sets)
        for i in range(N):
            if self._sets[i] == A:
                self._sets[i] = self._sets[i].union(B)
                break
        self._sets.remove(B)

    def union(self, a, b):
        N = len(self._sets)
        for i in range(N):
            if a in self._sets[i]:
                break
        assert(a in self._sets[i])

        for j in range(N):
            if b in self._sets[j]:
                break
        assert(b in self._sets[j])

        if i==j: return
        self._sets[i] = self._sets[i].union( self._sets[j])
        del self._sets[j]

    def __iter__(self):
        return iter(self._sets)

    def __repr__(self):
        return "UnionFind({})".format(self._sets)

    
class Algo:
    def __init__(self, optimizer, prototype_maker, score, relaxer, fasta_file, threshold=0.0):
        self._optimizer = optimizer
        self._prototype_maker = prototype_maker
        self._score = score
        self._relaxer = relaxer
        self.fasta_file = fasta_file
        self.threshold = threshold
        
        self.nodes = []
        N = hmmdsl_py.count_fasta(fasta_file)
        for n in range(N):
            newnode = Node([n], self._prototype_maker())
            self.nodes.append( newnode )
        #self.sets = UnionFind(range(N))
            
    def __repr__(self):
        ret = "Algo("
        ret += ",".join(map(lambda x: repr(x), self.nodes))
        ret += ")"
        return ret

    def iter(self, relax=False):
        # Update models for all nodes
        for node in self.nodes:
            if( not node.model_up_to_date ):
                # Optimize the model for the sequences in the cluster
                self._optimizer( node.model, node.elems() ) # TODO - set initial model
                if relax:
                    # Relax the model
                    self._relaxer( node.model )
                # Mark the model as up-to-date
                node.model_up_to_date = True

        # Calculate NxN similarity matrix
        N = len(self.nodes)
        scores = make_square_matrix(N)
        
        total_scores = N**2
        scores_so_far = 0
        for i in range(N):
            mi = self.nodes[i].model
            sigma_score = 0.0
            for j in range(N):
                scores[i][j] = self._score( mi, self.nodes[j].elems() ) # score model i against the members of each node
                if scores_so_far % 10 == 9:
                    print("(calculating scores, {}% done)".format(float(scores_so_far)/float(total_scores)*100))
                scores_so_far += 1
                sigma_score += scores[i][j]
            print("Node {}: avg. score = {}".format(i, sigma_score/float(N)) )


        # Find clusters to merge
        # Note: similarity is not guaranteed to be symmetric, so we need to iterate over all Nx(N-1) pairs

        sets = UnionFind(range(N))

        for i in range(N):
            for j in range(N):
                if i==j: continue
                if scores[i][j] < self.threshold: continue
                # Mark i and j for merging
                sets.union(i, j)

        newnodes = []
        delnodes = []

                
        #print("Nodes: ", self.nodes)
        #print("Sets: ", sets)
        # Update the clusters according to the disjoint sets
        for s in sets:
            if len(s) > 1:
                members = []
                members.extend(map( lambda x: self.nodes[x], s ))
                newnode = Node( members, self._prototype_maker())
                newnodes.append(newnode)
                delnodes.extend(s)

        delnodes.sort()
        delnodes.reverse()        
        for i in delnodes:
            del self.nodes[i]

        self.nodes.extend(newnodes)
        return len(newnodes)>0
        
        

        
