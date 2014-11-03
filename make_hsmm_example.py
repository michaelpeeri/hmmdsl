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
import sys
import getopt
import collections # for deque
from math import exp
import hmmdsl_py
import dotio
import training

# Source: Wikipedia
#background_emissions = {'A': 0.078,'C': 0.019,'D': 0.053,'E': 0.063,'F': 0.039,'G': 0.072,'H': 0.023,'I': 0.053,'K': 0.059,'L': 0.091,'M': 0.023,'N': 0.043,'P': 0.052,'Q': 0.042,'R': 0.051,'S': 0.068,'T': 0.059,'V': 0.066,'W': 0.014,'Y': 0.032 }

# Source: python3 measure_background_distribution.py on all S. cer. proteins
background_emissions = {'A': 0.054771, 'P': 0.043866, 'F': 0.045158, 'I': 0.065673, 'C': 0.013187, 'H': 0.021793, 'Y': 0.033828, 'W': 0.010448, 'Q': 0.039149, 'D': 0.057591, 'S': 0.090586, 'M': 0.021020, 'E': 0.064357, 'R': 0.044566, 'G': 0.049523, 'N': 0.061082, 'T': 0.059083, 'K': 0.072734, 'L': 0.095785, 'V': 0.055800}


def SetBackgroundEmissionProbs1(model, state):
    for k, v in background_emissions.items():
        model.SetEmissionProbability( state, k, v )

def AddAlphabet(model, alphabet):
    verify = {}
    for symbol in alphabet:
        model.AddAlphabetSymbol( symbol )
        verify[symbol] = 1
    assert(len(alphabet) == len(verify.keys()))

def AddStates(model, states, captions):
    assert(len(states)==len(captions))

    for i in range(len(states)):
        model.AddState( states[i], captions[i] )

def SetBackgroundEmissionProbs(model, states):    
    for s in states:
        SetBackgroundEmissionProbs1( model, s );

def SetNeutralTransitions( model, graph ):
    for s in graph.keys():
        targets = graph[s]
        p = 1. / len(targets)
        for t in targets:
            model.SetTransition( s, t, p )

def SetEmissionsClass( model, *group ):
    for i in range(len(group)-1):
        model.MergeEmissionsClass( group[i], group[i+1] )


# Add constrained state, with background probability for all other states
# Make sure the resulting probabilities are normalized
# Note: dotio prints emissions > 0.005, so background should probably be lower than that
def SetConstrainedEmissions( model, state, constraints, background, alphabet ):
    e = {}
    sum = 0.0

    # Set the probability for each symbol (either as set in the constraint, or at background level)
    for k in alphabet:
        p = 0.0
        if( k in constraints):
            p = constraints[k]
        else:
            p = background(k)
        assert((p>=0.0) & (p<=1.0))
        sum += p
        e[k] = p

    verify = 0.0
    factor = 1.0/sum
    # Send the normalized probabilities to the model
    for k, v in e.items():
        verify += v * factor
        model.SetEmissionProbability( state, k, v * factor )
    assert(abs(verify-1.0)<1e-10)


def InitializePrototype(model, backgroundModel = False):	
    # Define alphabet
    alphabet = ['E', 'D', 'R', 'H', 'K', 'I', 'L', 'V', 'S', 'T', 'Q', 'W', 'G', 'P', 'F', 'Y', 'C', 'M', 'A', 'N']
    AddAlphabet( model, alphabet )
    # Define node identifiers
    states = range(2,31+2)
    I = 0
    T = 1
    (                 Sa0, Sa1, Sa2, Sa3, Sa4, Sa5, Sa6, Sa7, Sa8, Sa9, Sa10, Sa11, Sa12, Sa13, Sa14, K1,  Sb0, Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14   ) = states
    captions = [      "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",  "-",  "-",  "-",  "-",  "K", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",  "-",  "-",  "-",  "-"  ]

    if( backgroundModel ):
        captions[14] = "-"

    AddStates( model, states, captions );

    model.SetDefaultProbabilities();

    graph = { I:     [Sa0, Sa1, Sa2, Sa3, Sa4, Sa5, Sa6, Sa7, Sa8, Sa9, Sa10, Sa11, Sa12, Sa13, Sa14, K1],
              Sa0:   [Sa1],
              Sa1:   [Sa2],
              Sa2:   [Sa3],
              Sa3:   [Sa4],
              Sa4:   [Sa5],
              Sa5:   [Sa6],
              Sa6:   [Sa7],
              Sa7:   [Sa8],
              Sa8:   [Sa9],
              Sa9:   [Sa10],
              Sa10:  [Sa11],
              Sa11:  [Sa12],
              Sa12:  [Sa13],
              Sa13:  [Sa14],
              Sa14:  [K1],
              K1:    [Sb0, Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14, T],
              Sb0:   [Sb1],
              Sb1:   [Sb2],
              Sb2:   [Sb3],
              Sb3:   [Sb4],
              Sb4:   [Sb5],
              Sb5:   [Sb6],
              Sb6:   [Sb7],
              Sb7:   [Sb8],
              Sb8:   [Sb9],
              Sb9:   [Sb10],
              Sb10:  [Sb11],
              Sb11:  [Sb12],
              Sb12:  [Sb13],
              Sb13:  [Sb14],
              Sb14:  [T ] }
    SetNeutralTransitions( model, graph )

    SetBackgroundEmissionProbs( model,
                                [Sa0, Sa1, Sa2, Sa3, Sa4, Sa5, Sa6, Sa7, Sa8, Sa9, Sa10, Sa11, Sa12, Sa13, Sa14] )
    SetBackgroundEmissionProbs( model,
                                [Sb0, Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14] )


    background = lambda x: 0.004
    default = lambda x: background_emissions[x]*0.2


    if( not backgroundModel ):
        SetConstrainedEmissions( model, K1, { 'R': 0.35, 'K': 0.65  }, default, alphabet )  # Arbitrary distribution for positively-charged residues
    else:
        SetBackgroundEmissionProbs( model,
                                    [K1] )
        
    #SetConstrainedEmissions( model, N1, { 'D': 0.4, 'E': 0.6 }, background, alphabet ) # Arbitrary distribution for negatively-charged residues


    #for k in [S0,  S2 ]:
    #    SetConstrainedEmissions( model, k, { 'S': 0.0001, 'T': 0.0001 }, default, alphabet )
    #
    #for k in [S1]:
    #    SetConstrainedEmissions( model, k, { 'S': 0.54, 'T': 0.46 }, background, alphabet )

    (Learned, Fixed) = (0,1)
    (Emissions, Transitions, Durations) = (0, 1, 2)

    for k in [Sa0, Sa1, Sa2, Sa3, Sa4, Sa5, Sa6, Sa7, Sa8, Sa9, Sa10, Sa11, Sa12, Sa13, Sa14, Sb0, Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14 ]:
        #model.SetMu (k, 1.0)
        #model.SetEta(k, 0.069)
        model.SetMu (k, 1.0)
        model.SetEta(k, 0.01)

        model.SetStateLearningMode( k, Emissions,   Fixed )
        model.SetStateLearningMode( k, Transitions, Fixed )
        model.SetStateLearningMode( k, Durations,   Learned )
        

    for k in [K1]:
        model.SetMu (k, 12.0)
        model.SetEta(k, 1.5)

        model.SetStateLearningMode( k, Emissions,   Fixed )
        model.SetStateLearningMode( k, Transitions, Fixed )
        model.SetStateLearningMode( k, Durations,   Learned )

    if( not backgroundModel):
        model.debug_print()


# Simple nullary function returning RNG-generated numbers, seeded using /dev/random
class _rng:
    def __init__(self):
        from random import random, Random, SystemRandom
        init = SystemRandom()
        self.rng = Random(init.random())
    def __call__(self):
        return self.rng.random()

#a = UtilityMaximiser(3)
#a.SetUtility(0, 5)
#a.SetUtility(1, 2)
#a.SetUtility(2, 3)

def main_test(fasta_file, repeats, args = None):

    perturb_amount = 0.8
    rnd = _rng()

    for i in repeats:


        model = hmmdsl_py.Model()
        InitializePrototype(model)

        #perturb = hmmdsl_py.Perturb(model)
        #perturb.perturb_emissions(perturb_amount)

        #for k in range(2,14+2):
        #    model.SetMu (k, rnd()*3.0 + 1.0 )

        dotio.WriteDOT(model, "hsmm_local_prototype_04.dot" )
        dotio.store(model, "hsmm_local_prototype_04.pickle" )


        backgroundmodel = hmmdsl_py.Model()
        InitializePrototype(backgroundmodel, True)

        dotio.WriteDOT(backgroundmodel, "hsmm_local_prototype_04_background.dot" )
        dotio.store(backgroundmodel, "hsmm_local_prototype_04_background.pickle" )

        return

        #model.split_state( 4, 0.7 )
        
        #dotio.WriteDOT(model, "model_split.dot" )
        #dotio.store(model, "model_split.pickle" )

        training.train(fasta_file, model, [])
        
        #model.debug_print()
        dotio.WriteDOT(model, "hsmm14_perturb{}_vs_effectros.faa.dot".format(i) )
        
        dotio.store(model, "hsmm14_perturb{}_vs_effectors.faa.pickle".format(i) )


class UsageError(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError as err:
            raise UsageError(err)

        if( len(argv) < 2 ):
            raise UsageError("No file specified");

        if( len(argv) < 4 ):
            raise UsageError("No range specified");

        i = int(argv[2])
        j = int(argv[3])
        repeats = range(i,j)

        main_test(argv[1], repeats, argv)

    except UsageError as err:
        print( err.msg )
        return 2;



if __name__ == "__main__":
    sys.exit(main())

