{-# LANGUAGE FlexibleContexts #-}

module GraphFRI where

import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Data (toList)
import Numeric.LinearAlgebra.Devel (Matrix(..), hmatMul)
import System.Random

-- Graph Fourier Transform (GFT) matrix computation
graphFourierTransform :: Matrix Double -> Matrix Double
graphFourierTransform circulantGraph =
    let laplacian = (diag $ sumElements circulantGraph) - circulantGraph
        (eigvals, eigvecs) = eigSH' laplacian (size laplacian - 1)
        gft = takeColumns (size circulantGraph - 1) eigvecs
    in gft

-- Graph-FRI sampling framework
graphFRISampling :: Matrix Double -> Vector Double -> Int -> Vector Double
graphFRISampling circulantGraph signal k =
    let gft = graphFourierTransform circulantGraph
        signalGFT = gft #> signal
        (_, idxs) = unzip . take k . reverse . sort . zip (toList signalGFT) $ [0..]
        sparseSignalGFT = accum (konst 0 (size signalGFT)) (+) (zip idxs (toList signalGFT))
        reconstructedSignal = pinv gft #> sparseSignalGFT
    in reconstructedSignal

-- Generate a random circulant graph with given size
generateCirculantGraph :: Int -> StdGen -> (Matrix Double, StdGen)
generateCirculantGraph n gen =
    let weights = take n (randoms gen :: [Double])
        circulantGraph = diag weights + (tr . shift 1 . fromList $ take (n-1) (cycle [-1])) + (shift (-1) . fromList $ take (n-1) (cycle [-1]))
    in (circulantGraph, snd $ split gen)

-- Generate a random signal with given size
generateSignal :: Int -> StdGen -> (Vector Double, StdGen)
generateSignal n gen = let (vals, newGen) = splitAt n (randoms gen :: [Double]) in (fromList vals, newGen)

-- Test the Graph-FRI framework
testGraphFRI :: Int -> Int -> Int -> Int -> IO ()
testGraphFRI numGraphs n k seed = do
    gen <- newStdGen
    let testOneGraph _ 0 _ _ acc _ _ = return acc
        testOneGraph gLeft gRight n' k' acc g = do
            let (circulantGraph, g') = generateCirculantGraph n' g
                (signal, g'') = generateSignal n' g'
                reconstructedSignal = graphFRISampling circulantGraph signal k'
                reconstructionError = norm_2 (signal - reconstructedSignal) / fromIntegral n'
            testOneGraph (gLeft - 1) gRight n' k' (acc + reconstructionError) g''
        testGraph _ 0 _ _ acc = return acc
        testGraph gLeft gRight n' k' acc = do
            let (totalError, g') = testOneGraph numGraphs numGraphs n' k' 0 gen
            testGraph (gLeft - 1) gRight n' k' (acc + totalError / fromIntegral numGraphs)
    finalError <- testGraph numGraphs numGraphs n k 0
    putStrLn $ "Average reconstruction error over " ++ show numGraphs ++ " graphs with " ++ show n ++ " vertices and sparsity " ++ show k ++ " is " ++ show finalError

main :: IO ()
main = testGraphFRI 10 20 5 123456789
