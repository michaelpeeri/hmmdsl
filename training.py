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


class CircularCounter:
    def __init__(self, size):
        self._items = collections.deque([], size)
        self._counts = {}
        self._size = size

    def push(self, item):
        if( len(self._items) == self._items.maxlen ):
            removed_item = self._items.pop()
            self._counts[removed_item] -= 1
        self._items.appendleft(item)
        self._update(item, 1)

    def _update(self, item, increment):
        if item in self._counts:
            self._counts[item] += increment
        else:
            self._counts[item] = 1

    def counts(self):
        return self._counts

    def size(self):
        return self._size

    def __iter__(self):
        return self._items.__iter__()

    def __len__(self):
        return len(self._items)

def normalize(vector):
    factor = max(vector)
    return list(map(lambda x: x/factor, vector))


class UtilityMaximiser:
    def __init__(self, num_options):
        self._history = CircularCounter(25) # Store the last N options
        self._utility = [0] * num_options # Store the last utility for each option
        self._num_options = num_options

    def SetUtility(self, option, utility):
        self._utility[option] = utility

    def _CalcOverdueOption(self, duration=6):
        # If we're only starting out, no option can be overdue
        if( len(self._history) < duration ):
            return None

        # Are there any elements that haven't been returned recently?
        not_returned_recently = set(range(self._num_options)).difference(set(self._history))
        if( len(not_returned_recently) > 0 ):
            chosen = not_returned_recently.pop() # return an arbitrary item
            
        else:
            # Choose the oldest element
            hist = list(self._history)
            last_return = list(map( lambda x: hist.index(x), range(self._num_options)))
            if( max(last_return) < duration ):
                return None

            chosen = last_return.index(max(last_return))

        return chosen

    def GetOverdueOption(self, duration=6):
        chosen = self._CalcOverdueOption(duration)
        if( chosen != None ):
            self._history.push(chosen)
        return chosen

    def _CalcNextOption(self):
        N = self._history.size()
        counts = self._history.counts()

        current = [0] * len(self._utility)

        for k, v in counts.items():
            current[k] = v

        actual = list(map(lambda x: float(x)/N, current ))
        wanted = normalize(self._utility)

        difference = list(map(lambda x,y:x-y, wanted, actual))
        chosen = difference.index(max(difference))
        return chosen

    def GetNextOption(self):
        # Try overdue options first
        chosen = self._CalcOverdueOption()
        if( chosen == None ):
            # Otherwise, return the option with the greatest utility
            print(" (utility)")
            chosen = self._CalcNextOption()
        else:
            print(" (overdue)")


        self._history.push(chosen)
        return chosen
    
def train(fasta_file, model, elems=[], radius=5e-4):

    print("Starting training heuristic")
    print("Training-set: {}".format(fasta_file))
    if( len(elems)==0 ):
        print("Training subset: Full")
    else:
        print("Training subset: {} (N={})".format(elems, len(elems)))

    print("Model class: {}".format(type(model)))
    print("Convergence radius: {}".format(radius))

    #InitializeTestModel(model)
    #model.debug_print()
    em = None
    if isinstance(model, hmmdsl_py.Model):
        # HSMM
        em = hmmdsl_py.EM(model, fasta_file, elems)
        (train_e, train_scale, train_shape) = range(3)
        mode_desc = ["e", "scale", "shape"]
        f1 = lambda em: em.reestimate_b()
        f2 = lambda em: em.reestimate_scale()
        f3 = lambda em: em.reestimate_shape()
        train_func = [f1,f2,f3]

        id = UtilityMaximiser(3)
        id.SetUtility(train_e, 0.10)
        id.SetUtility(train_scale, 0.02)
        id.SetUtility(train_shape, 0.01)
    else:
        # HMM
        em = hmmdsl_py.HMMEM(model, fasta_file, elems)

        (train_e, train_a) = range(2)
        #(train_e, ) = range(1)
        mode_desc = ["e", "a"]
        #mode_desc = ["e"]
        f1 = lambda em: em.reestimate_b()
        f2 = lambda em: em.reestimate_a()
        train_func = [f1,f2]
        #train_func = [f1]

        id = UtilityMaximiser(2)
        #id = UtilityMaximiser(1)
        #id.SetUtility(train_e, 0.10)
        id.SetUtility(train_a, 0.02)


    #dotio.WriteDOT(model, "model_initial.dot" )

    from os import getpid
    #log = open("trainlog{}.dat".format(getpid()), "w")

    prev = 0.0
    iter = 0
    train_mode = train_e
    stable_count = 0
    phist = collections.deque([], 8)
    overdue_only_mode = False
    
    
    while(True):
        print("Iter {}{}".format(iter, mode_desc[train_mode]))
        
        p =  train_func[train_mode](em)
        #log.write("{} {}\n".format(iter, p))
        #log.flush()

        if( iter > 0 ):
            delta = (prev - p)/prev
            print( "     {} ({}%) {} {}".format(p-prev, delta*100, prev, p) )
            id.SetUtility( train_mode, max(delta*100, 0.00001))

            if( delta*100 < radius ):
                stable_count += 1
            else:
                stable_count = 0

            if( ( stable_count > 9) or ( min(phist) > p) ):
                print("Convergence criteria met; will exit if no training params are overdue; delta={}%, stable_count={}, min(phist)>p={}".format(delta*100, stable_count, min(phist)>p))
                overdue_only_mode = True
            else:
                overdue_only_mode = False


        phist.appendleft(p)

        # Setup next iter
        if( (iter%25 == 0) and (iter>0) ):
            #model.debug_print()
            #dotio.WriteDOT(model, "model_iter{}.dot".format(iter) )
            pass

        if( overdue_only_mode ):
            train_mode = id.GetOverdueOption()
            if( train_mode == None ):
                break
        else:
            train_mode = id.GetNextOption()

        prev = p
        iter += 1
         
    return 0.0
