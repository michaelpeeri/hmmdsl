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
from math import exp, sqrt
import hmmdsl_py
import pickle

def WriteGraphToDOT( graph, filename ):
    f = open( filename, 'w')
    f.write('digraph hmm {\n')
    f.write('//plot using: dot -Tpng {0} >{0}.png\n'.format( filename ) )
    f.write('\trankdir=LR;\n')

    for s in graph.keys():
        targets = graph[s]
        p = 1. / len(targets)
        for t in targets:
            f.write( '\tS{0} -> S{1} [penwidth={2}, label=\"{3:1.3f}\"];\n'.format( s, t, p*10.0, p ) );

    f.write('}')

def AddAlphabet(model, alphabet):
    for symbol in alphabet:
        model.AddAlphabetSymbol( symbol )

class StorableModel:
    def __init__(self, model):
        self.symbols = []
        for c in range(model.num_symbols()):
            self.symbols.append( model.get_symbol(c) )

        self.states = []
        for s in range(model.num_states()):
            self.states.append( model.GetStateName(s) )

        # Initialize 'a'
        self.a = []
        for k in range(model.num_states()):
            row = []
            for l in range(model.num_states()):
                row.append( model.a(k, l) )
            self.a.append(row)

        # Initialize 'e'
        self.e = []
        for k in range(model.num_states()):
            row = []
            for s in range(model.num_symbols()):
                row.append( model.e(k, model.get_symbol(s) ) )
            self.e.append(row)

        if isinstance(model, hmmdsl_py.Model):
            # Initialize 'eta'
            self.eta = []
            for k in range(model.num_states()):
                self.eta.append( model.GetEta(k) )
                
            # Initialize 'mu'
            self.mu = []
            for k in range(model.num_states()):
                self.mu.append( model.GetMu(k) )

        self.model_type = type(model)

    # Return an equivalent model to the pickled one
    def convert(self):
        try:
            model_type = self.model_type
        except AttributeError as err:
            self.model_type = hmmdsl_py.Model

        if self.model_type == hmmdsl_py.Model:
            model = hmmdsl_py.Model()
        else:
            model = hmmdsl_py.HMMModel()

        AddAlphabet( model, self.symbols );
        for s in range( 2, len(self.states ) ):
            model.AddState( s, self.states[s] )

        for k in range(len(self.states)):
            #if( k == model.GetTerminalState() ): continue
            for l in range(len(self.states)):
                #if( l == model.GetInitialState() ): continue
                model.SetTransitionLogspace(k, l, self.a[k][l] )

            if isinstance(model, hmmdsl_py.Model):
                model.SetEta( k, self.eta[k] )
                model.SetMu(  k, self.mu[k] )

        for k in range(len(self.states)):
            if( model.is_silent(k) ): continue
            for s in range(len(self.symbols)):
                model.SetEmissionProbabilityLogspace( k, self.symbols[s], self.e[k][s] )
            
        return model
        

def store(model, path):
    with open(path, "wb") as f:
        pickle.dump( StorableModel(model), f, -1 )

def load(path):
    with open(path, "rb") as f:
        sm = pickle.load( f )
    return sm.convert()
    

def GetStateLabel( model, state ):
    label = '<<table border="0" cellborder="0" cellpadding="0" cellspacing="0">'

    if isinstance(model, hmmdsl_py.Model):
        mu = model.GetMu(state)
        eta = model.GetEta(state)
        mean = mu/eta
        variance = mu/(eta**2)
        stddev = sqrt(variance)
        label += '<tr><td colspan="2">&#x3C4;={:.3g}&#x0B1;{:.2g}</td></tr>'.format(mean, 2*stddev)
        label += '<tr><td colspan="2">&#956;={:.3g}</td></tr>'.format(mu)
        label += '<tr><td colspan="2">&#x3B7;={:.3g}</td></tr>'.format(eta)

    entropy = model.emissions_entropy(state)
    for c in range(model.num_symbols()):
        symbol = model.get_symbol(c);
        e = exp(model.e( state, symbol ))
        if( e >= 0.01 ):
            label += '<tr><td width="{:.3g}" height="5" fixedsize="true" bgcolor="black" align="left"></td><td width="20" height="10">{}</td></tr>'.format( e*80, symbol )
    label += '<tr><td colspan="2">H(e)={:.2g}b</td></tr>'.format(entropy)
    label += '</table>>'
        
    return label

def WriteDOT( model, filename ):
    f = open( filename, 'w')
    f.write('//plot using: dot -Tpng {0} >{0}.png\n'.format( filename ) )
    f.write('digraph hmm {\n')
    f.write('\trankdir=LR;\n')


    # Write the state definitions
    # Note: The initial state is included in these states
    for s in range(model.num_states()):
        if( not model.is_silent(s) ):
            f.write( '\tS{0} [shape=Mrecord label={1}];\n'.format( s, GetStateLabel( model, s) ) );
        else:
            f.write( '\tS{0} [shape=Mrecord label={1}];\n'.format( s, s) ); # Don't write emission probs for silent states

    # Write the transitions
    for s in range(model.num_states()):
        for t in range(model.num_states()):
            p = exp(model.a(s, t));

            if( p >= 5.e-3 ):      # Only include transitions with probability >= 0.005 (this should exclude pseudocounts)
                f.write( '\tS{0} -> S{1} [penwidth={2}, label=\"{3:1.3f}\"];\n'.format( s, t, p*10.0, p ) );

    f.write('}')

