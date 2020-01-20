# FasPI

## Overview
FastPI (Fast PseudoInverse) is a novel method computing pseudoinverse for sparse feature matrices used in real-world optimization problems. 
Based on the observation that many real-world feature matrices are sparse and highly skewed, FastPI reorders and divides the feature matrix and incrementally computes low-rank SVD from the divided componenets. 
The details of FastPi are described in the following paper:
* Fast and Accurate Pseudo inverse for Real-world Sparse Matrices  
  29th International Joint Conference on Artificial Intelligence (IJCAI 2020) submitted
  
## Usage

You first need to compile the following c++ code, called `ComputeConnComp.cpp`, using mex in MATLAB. 
Type the following command in MATLAB:
```
>> mexCompile
```

This will generate a compiled file according to your system and OS. 
If you encounter an error about mex, please check your mex setup first. 

You can check a demo about FastPI by typing the following command:

```
>> demo
```

The main code of FastPI is implemented in `FastPI.m`. 
The input and output of the function are as follows:
* Input
  - A: feature matrix (m x n)
  - alpha: target rank ratio (0 < alpha <= 1.0)
* Output
  - [V, pinvS, UT]: results for the pseudoinverse based on SVD such that pinvA = V * pinvS * UT
  - rank: target rank according to alpha

The example of the usage of `FastPI` is as follows:
```
alpha = 0.1;
[V, pinvS, UT, rank] = FastPI(A, alpha);
```

If you need the explicit pseudoinverse `pinvA`, then it is computed as the followings:
```
pinvA = (V * pinvS) * UT
```

However, this requires huge memory space if the given matrix is large. 
The typical usage is to use them as an implicit pseudoinverse. 
For example, suppose we need to perform a matrix-vector multiplication, i.e., `pinvA * x` where `x` is a vector. 
Then, `pinvA * x = (V * (pinvS * (UT * x)))` is much more efficient than the explicit version. 
