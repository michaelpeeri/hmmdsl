-- -------------------------------------------------------------------------------------
--   Copyright 2014 Michael Peeri
-- 
--   This file is part of hmmdsl.
--   hmmdsl is free software: you can redistribute it and/or modify
--   it under the terms of the GNU General Public License as published by
--   the Free Software Foundation, either version 3 of the License, or
--   (at your option) any later version.
-- 
--   hmmdsl is distributed in the hope that it will be useful,
--   but WITHOUT ANY WARRANTY; without even the implied warranty of
--   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
--   GNU General Public License for more details.
-- 
--   You should have received a copy of the GNU General Public License
--   along with hmmdsl.  If not, see <http://www.gnu.org/licenses/>.
-- -------------------------------------------------------------------------------------
import Data.Matrix (Matrix, fromList, fromLists, nrows, ncols, (!))
import Math.Gamma (gamma)
import Math.Gamma.Incomplete (pHypGeom)
import Data.List (intersect)
import qualified Data.Map (fromList, findWithDefault)


-- ===============================================================
-- Discretized gamma distribution
-- ===============================================================
q_gamma :: Int -> Double -> Double -> Double
q_gamma d eta mu = discretize d eta mu
    where
    gamma_dist::Double->Double->Double->Double
    gamma_dist _d eta mu = (mu**eta)*(_d**(eta-1.0))*exp(-mu*_d)/gamma(eta)

    gamma_cdf::Double->Double->Double->Double
    gamma_cdf x eta mu = (pHypGeom eta (x/mu)) / gamma(eta)

    discretize::Int->Double->Double->Double
    discretize  0 _   _  = error "q_gamma called with d=0!"
    discretize  1 eta mu = (gamma_cdf 1.5 eta mu) - 
                           0.0   -- (gamma_cdf 0 eta mu) == 0
    discretize _d eta mu = (gamma_cdf ((fromIntegral _d)+0.5) eta mu) - 
                           (gamma_cdf ((fromIntegral _d)-0.5) eta mu)
_D = 100

hsmm_algo :: [Char] -> [Char] -> Matrix Double -> Matrix Double -> Matrix Double -> ((Int -> Int -> Double), (Int -> Int -> Double), (Int -> Int -> Double), (Int -> Int -> Double), Double, IO ())
hsmm_algo alphabet seq a e d' = (forward, forward_begin, backward, backward_begin, pseq, debug_print)
    where

    init = 1
    term = 2

    -- ===============================================================
    --  forward [HSMM]
    --  ref: Rabiner eq. (68)
    -- ===============================================================
    forward :: Int -> Int -> Double
    forward t 1 = if (t==1) then 1.0 else 0.0
    forward t 2 = 0.0
    forward t i = sum [
            (fb' ! (t-d, i)) *
            (p i d) *
            (product [e!(i, (sym (seq!!s))) | s <- (t_to_seq (intersect [(t-d+1)..t] [2..(_T+1)] )) ] )
            | d <- [1..(min _D (t-1))]
            ]

    -- ===============================================================
    --  forward_begin [HSMM]
    --  ref: Rabiner eq. (71)
    -- ===============================================================
    forward_begin :: Int -> Int -> Double
    forward_begin t j | (t<(_T+1)) && (j==2) = 0.0
                      | otherwise = sum [(f' ! (t, i)) * (a ! (i, j)) | i <- [1.._N]]

    
    -- ===============================================================
    --  backward [HSMM]
    --  ref: Rabiner eq. (76)
    -- ===============================================================
    backward :: Int -> Int -> Double
    backward t i | (t==_T+2) && (i/=2) = 0.0
                 | (t==_T+2) && (i==2) = 1.0
                 | (t>1)     && (i==1) = 0.0
                 | otherwise = sum [ (a ! (i, j)) * (bb' ! (t, j))    | j <- [1.._N] ]

    -- convert from 't' coordinates (2..T+1) to 'sequence' coordinates (0..T-1)
    t_to_seq :: [Int] -> [Int]
    t_to_seq v = (map (\x -> x-2) v )

    -- ===============================================================
    --  backward_begin [HSMM]
    --  ref: Rabiner eq. (77)
    -- ===============================================================
    backward_begin :: Int -> Int -> Double
    backward_begin t 2 = if (t==(_T+1)) then 1.0 else 0.0
    backward_begin t 1 = 0.0
    backward_begin t i = sum [
                   (b' ! (t+d, i)) *
                   (p i d) *
                   (product [e!(i, (sym (seq!!s))) |  s <- (t_to_seq (intersect [(t+1)..(t+d)] [2..(_T+1)] ) )])
                     | d <- [1..(min _D (_T-t+1))]
                     ]

    -- ===============================================================
    --  P(O|model)
    -- ===============================================================
    pseq :: Double
    --pseq = sum [(f' ! (_T+1, i)) * (a ! (i, j)) | i <- [1.._N]]
    pseq = fb' ! (_T+2, 2)


    debug_print :: IO ()
    debug_print = do
        let t = 5
        let d = 4
        let i = 3
        --putStrLn("p(d_3)")
        --putStrLn(show [p 3 d | d <- [1..100]])
        putStrLn("forward:")
        putStrLn((show f'))
        putStrLn("forward_begin:")
        putStrLn((show fb'))
        putStrLn("backward:")
        putStrLn((show b'))
        putStrLn("backward_begin:")
        putStrLn((show bb'))
        --putStrLn((show (intersect [(t+1-2)..(t+d-2)] [2..(_T+1)] ) ) )

        --putStrLn((show [( t_to_seq (intersect [(t+1)..(t+dd)] [2..(_T+1)] )) | dd <- [1..20] ] ))

        --putStrLn(show _T)
        --putStrLn(show [ [1..(min _D (_T+2-tt+1))] | tt <- [1..(_N+2)]  ] )

        --putStrLn("------");
        --putStrLn(show [e!(3, (sym (seq!!s))) |  s <- (t_to_seq (intersect [(t+1)..(t+d)] [2..(_T+1)] ) )])
        --putStrLn(show (product [e!(i, (sym (seq!!s))) |  s <- (t_to_seq (intersect [(t+1)..(t+d)] [2..(_T+1)] ) )]))

        --putStrLn(show [[1..(min _D (_T-tt+1))] | tt <- [1..8] ])

    -- ===============================================================
    --  calculate p(i,d) -- gamma durations distribution
    -- ===============================================================
    p :: Int -> Int -> Double
    p i d = q_gamma d (d'!(i,1)) (d'!(i,2))

    _N = nrows a
    _T = length seq

    -- map seq observations to positions in the 'e' matrix
    -- note: parameters must appear in the same order in 'alphabet' and 'e' params
    _convert = Data.Map.fromList (zip alphabet [1..])
    sym :: Char -> Int
    sym c = Data.Map.findWithDefault (-1) c _convert

    -- Use the lazy evaluation of list comprehensions to memoize recurrence values
    -- TODO: add ref for this technique
    f'  = fromList (_T+2) _N [forward        t i| t <- [1..(_T+2)], i <- [1.._N]]
    fb' = fromList (_T+2) _N [forward_begin  t i| t <- [1..(_T+2)], i <- [1.._N]]
    b'  = fromList (_T+2) _N [backward       t i| t <- [1..(_T+2)], i <- [1.._N]]
    bb' = fromList (_T+2) _N [backward_begin t i| t <- [1..(_T+2)], i <- [1.._N]]



main = do
     let alphabet = "abcd"
     let seq = "aaaaaccc"

     let _a1 = fromLists [[0.0::Double, 0.0, 0.5, 0.5], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 0.0]]
     let _b1 = fromLists [[0.0::Double, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.6, 0.2, 0.1, 0.1], [0.05, 0.05, 0.85, 0.05] ]
     let _d1 = fromLists [[0.0::Double, 0.0], [0.0, 0.0], [2.0, 12.0], [2.0, 12.0]]


     --putStrLn ("a:")
     --putStrLn (show _a1)
     --putStrLn ("e:")
     --putStrLn (show _b1)
     --putStrLn ("d:")
     --putStrLn (show _d1)


     --putStrLn ("test:")
     --putStrLn (show (_d1!(3,2)))

     let (f, fb, b, bb, pseq, debug_print) = hsmm_algo alphabet seq _a1 _b1 _d1

     putStrLn ("P(model|O): " ++ (show pseq))

     debug_print
