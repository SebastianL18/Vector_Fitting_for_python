# Fast Relaxed Vector Fitting for Python

vectfit3.py is the implementation of the Fast Relaxed Vector Fitting algortihm on python. This algorithm was originaly developed by B. Gustavsen [1] and the original code was written in Matlab eviroment, which is aviable in the SINTEF's web page [4]. 

The pourpose of this algorithm is to compute a rational approximation from tabuled data in the frequency domain for scalar or vectorized problems, the resulting model can be expressed in either state-space form or pole-residue form. The idea behind vector fitting comes from the following equation:
![vf_eq1](https://github.com/SebastianL18/vectfit3_for_python/assets/144065519/cb6f269c-a06e-47f0-b66a-45ddc6dddd9c)
where $F(s)$ are the frequency domain samples to be fitted, $s$ is an array with the complex frequency points of evaluation for $F(s)$, $p_k$ are a set of $n$ single poles to fit and $R_k$ are a set of $n$ residues matrixes. As $p_k$ are scalars and $R_k$ are matrixes, the approximation is formulated by using common poles for all elements in $F(s)$ and   Terms $d$ and $e$ are constants and proportional parameters 




This module is created as a complement for the LineParameters program, where propagation and 
characteristic admitance matrixes of frequency dependent transmission lines and cables must be fitted prior to build 
a time domain model for electromagnetic transients simulations.

    - Original Matlab code autor: Bjorn Gustavsen (08/2008)
    - Transcripted and adapted by: Sebastian Loaiza (03/2024)

 * References:
 
    [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency       
        domain responses by Vector Fitting", IEEE Trans. Power Delivery,        
        vol. 14, no. 3, pp. 1052-1061, July 1999.
        
    [2] B. Gustavsen, "Improving the pole relocating properties of vector
        fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,
        July 2006.
        
    [3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,
        "Macromodeling of Multiport Systems Using a Fast Implementation of
        the Vector Fitting Method", IEEE Microwave and Wireless Components 
        Letters, vol. 18, no. 6, pp. 383-385, June 2008.

    [4] B. Gustavsen, "User's Guide for vectfit3.m (Fast, Relaxed Vector 
        fitting)", SINTEF Energy Research, N-7465 Trondheim, Norway, 2008. Aviable 
        online: https://www.sintef.no/en/software/vector-fitting/downloads/#menu
        accesed on: 2/2/2024
 
 * Changes:
    
    - All options for vectfit3 configuration are defined as boolean variables, except asymp which has 3 posible states
    - The new option "lowert_mat" for vectfit3 configuration is added. This indicates when F(s) samples belong to a 
      lower triangular matrix function, that reduces the number of elements to fit for a symmetric matrix function
    - A new method to sort the poles computed during the identification process is implemented
