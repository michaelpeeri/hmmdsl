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

-- -------------------------------------------------------------------------------------
-- Some Hidden Semi-Markov Model functions
-- Read HSMM models from a file, and compute some HSMM-related matrices for them.
-- This was mainly written as a reference to test the C++ impl. of hmmdsl against.
--
-- Haskell requirements (see imports):
--
-- Data.Matrix (cabal install matrix)
-- Math.Gamma (cabal install gamma)
-- Regex.TDFA (cabal install regex-tdfa)
-- -------------------------------------------------------------------------------------
import Data.Matrix (Matrix, fromList, fromLists, nrows, ncols, (!), toLists)
import Math.Gamma (gamma)
import Math.Gamma.Incomplete (pHypGeom)
import Data.List (intersect, intersperse)
import qualified Data.Map (fromList, findWithDefault)
import System.Environment (getArgs)
import System.IO (openFile, hClose, IOMode(..), Handle, hGetLine, hPutStrLn)
import Text.Regex.TDFA ((=~))


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

hsmm_algo :: [Char] -> [Char] -> Matrix Double -> Matrix Double -> Matrix Double -> ((Int -> Int -> Double), (Int -> Int -> Double), (Int -> Int -> Double), (Int -> Int -> Double), (Int -> Int -> Double), (Int -> Int -> Double), Double, IO (), (String -> IO ()))
hsmm_algo alphabet seq a e d' = (forward, forward_begin, backward, backward_begin, gamma_ti, xi, pseq, debug_print, save)
    where

    init = 1
    term = 2

    -- convert from 't' coordinates (2..T+1) to 'sequence' coordinates (0..T-1)
    t_to_seq :: [Int] -> [Int]
    t_to_seq v = (map (\x -> x-2) v )

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
    --  gamma_ti [HSMM]
    --  ref: Rabiner eq. (80)
    --  TODO: optimize summing
    -- ===============================================================
    gamma_ti :: Int -> Int -> Double
    gamma_ti t i = sum [
                   ((fb' ! (tau, i)) * (bb' ! (tau, i))) - 
                   ((f'  ! (tau, i)) * (b'  ! (tau, i)))
                   | tau <- [1..t]
                   ]

    -- ===============================================================
    --  xi [HSMM]
    --  ref: Rabiner eq. (80)
    -- ===============================================================
    xi :: Int -> Int -> Double
    xi i j = sum [
             (f' ! (t, i) ) *
             (a ! (i, j) ) *
             (bb' ! (t, j) )
             | t <- [2..(_T+1)]
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
        putStrLn("gamma:")
        putStrLn((show g'))
        putStrLn("xi:")
        putStrLn((show xi'))

    save :: String -> IO ()
    save path = do
        handle <- openFile path WriteMode
        write_matrix handle f'  "forward"
        write_matrix handle fb' "forward_begin"
        write_matrix handle b'  "backward"
        write_matrix handle bb' "backward_begin"
        write_matrix handle g'  "gamma"
        write_matrix handle xi' "xi"
        hClose handle

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
    g'  = fromList (_T+2) _N [gamma_ti t i| t <- [1..(_T+2)], i <- [1.._N]]
    xi' = fromList _N _N [xi i j| i <- [ 1.._N], j <- [1.._N]]


read_string :: Handle -> String -> IO String
read_string handle expectedHeading = do
    heading <- hGetLine handle
    value <- hGetLine handle
    if( heading == expectedHeading ) then return value else error "Unexpected content!"

read_matrix :: Handle -> String -> (Double->Double) -> IO (Matrix Double)
read_matrix handle expectedHeading conv = do
    -- Read the heading line (e.g. [somematrix:15 20]
    headingLine <- hGetLine handle
    let heading = (headingLine =~ "[[]([a-z]+)[:]([0-9]+) ([0-9]+)[]]" :: [[String]])!!0

    if( heading!!1 /= expectedHeading ) then error "Unexpected heading" else return () -- TODO FIX THIS
    let _nrows = (read (heading!!2))::Int
    let _ncols = (read (heading!!3))::Int

    -- Read nrows lines 
    lines <- (sequence (replicate _nrows (hGetLine handle)))

    -- Parse the lines
    let parseVal = \x -> conv (read (x!!0)::Double)  -- Convert all values to Double
    let parseLine = \line -> map parseVal (line =~ "[-]?[0-9.]+([eE][+-]?[0-9]+)?"  :: [[String]])
    let rows = map parseLine lines

    -- Return the result as a Matrix
    let mtx = fromLists rows
    if( (nrows mtx) /= _nrows ) then error "Unexpected number of rows" else return ()
    if( (ncols mtx) /= _ncols ) then error "Unexpected number of cols" else return ()
    return mtx


read_data_file :: String -> [String] -> IO (String, String, [Matrix Double])
read_data_file path expectedMatrices = do
    handle <- openFile path ReadMode
    alphabet <- read_string handle "[alphabet]"
    seq <- read_string handle "[sequence]"
    a' <- read_matrix handle "a" exp
    e' <- read_matrix handle "e" exp
    d' <- read_matrix handle "d" (\x -> x)
    hClose handle
    return (alphabet, seq, [a', e', d'])

write_matrix :: Handle -> Matrix Double -> String -> IO ()
write_matrix handle mtx heading = do
    hPutStrLn handle ( "[" ++ heading ++ ":" ++ (show (nrows mtx)) ++ " " ++ (show (ncols mtx)) ++  "]"  )

    let writeRow = \row -> concat (intersperse " " (map show (map log row)))
    let x = map writeRow (toLists mtx)
    mapM (hPutStrLn handle) x
    return ()


main = do

    args <- getArgs
    let dataFile = args!!0

    (alphabet, seq, [_a1, _b1, _d1]) <- read_data_file dataFile ["a", "e", "d"]

    putStrLn(show _a1)
    putStrLn(show _b1)
    putStrLn(show _d1)

    let (f, fb, b, bb, g, xi, pseq, debug_print, save) = hsmm_algo alphabet seq _a1 _b1 _d1

    putStrLn ("P(model|O): " ++ (show pseq))

    debug_print

    save "results.dat"
