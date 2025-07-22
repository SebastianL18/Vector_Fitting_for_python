# Fast Relaxed Vector Fitting for Python

_vectfit3.py_ is the implementation of the Fast Relaxed Vector Fitting algorithm in python. This algorithm was originally developed by B. Gustavsen and A. Semlyem [^1] and the original code was written in MATLAB environment, which is available at the [SINTEF's web page](https://www.sintef.no/en/software/vector-fitting/downloads/#menu). 

The purpose of this algorithm is to compute a rational approximation from tabulated data in the frequency domain for scalar or vector problems. The resulting model can be expressed in either state-space form or pole-residue form. 

This module was created for anyone who needs to use vector fitting algorithm in a python project, however **embedding this module in any commercial software is strictly prohibited**, and if it is used **in scientific work, publications[^1],[^2] and[^3] must be referenced**. Additionally, _vectfit3.py_ is a complement for the _LineParameters program_, a personal and academic project in development, where propagation and characteristic admittance matrixes of frequency dependent transmission lines and cables must be fitted prior to build a time domain model for electromagnetic transient simulations. Such as in the previous case, must of vector fitting applications involves the synthesis of rational models for power systems elements, electronic components, or equivalent circuits. 

## Vector Fitting Algorithm

The idea behind vector fitting comes from the concept that any frequency domain function can be approximated as the ratio of two polynomials and that ratio can be rewritten as a sum of partial fractions. Hence, the fitting of a collection of frequency domain samples $F(s)$ can be expressed as:
![vf_eq1](https://github.com/SebastianL18/vectfit3_for_python/assets/144065519/cb6f269c-a06e-47f0-b66a-45ddc6dddd9c)
where $s$ is an array with the complex frequency points of evaluation for $F(s)$, $p_k$ are a set of $n$ single poles and $R_k$ are a set of $n$ residues matrixes. As terms $p_k$ are scalars and $R_k$ are matrixes, the approximation is formulated by using common poles for all elements to fit and individual residues, therefore if $F(s)$ is scalar, $R_k$ results scalar as well. Terms $d$ and $e$ are constants and proportional parameters and are optional. 

The fitting process is divided in two stages, the poles identification, and the residues identification. The poles identification process begins with a set of starting poles and the evaluation of (1). Then, an overdetermined linear equation system is formulated. This system is solved as a least-squares problem and its solution gives the residues of an auxiliar function called $sigma$. Finally, the new poles are computed by finding the associated zeros of $sigma$. For the residues identification those new poles are used to build another overdetermined equations systems to be solved as a least-squares problem, in order to obtain the residues of $F(s)$. For more details see[^1]

## How to use _vectfit3.py_
### _vectfit3.py_ module
This module **requires numpy, scipy and matplotlib modules installed**. To use this module, _vectfit3.py_ code can be cloned or directly downloaded from this repository. Afterward, from the main code the vector fitting function can be imported as follows:

    from vectfit3 import vectfit

### _vectfit_ configuration
The vector fitting function has many options that can be modified via _opts_ dictionary. _opts_ dictionary contains the following keys, and each one of them is related to a modifiable option as described:

* **"symm_mat"**: _type=bool_. Indicates when $F(s)$ samples belong to a symmetric matrix function.
* **"RMO_data"**: _type=bool_. Indicates when matrix function elements are arranged in RMO into F(s)
* **"realx"** _type=bool_. Enables vector fitting with relaxed nontriviality.
* **"stable"**: _type=bool_. Enables stable poles enforcement.
* **"asymp"**: _type=int_, value=1, 2 or 3. Produces [D=0; E=0], [D!=0; E=0] or [D!=0; E!=0] respectively[^4].
* **"skip_pole"**: _type=bool_. Disables pole identification omission.
* **"skip_res"**: _type=bool_.: Disables residue identification omission.
* **"cmplx_ss"**: _type=bool_. Produces complex or real state space model[^4]
* **"spy1"**: _type=bool_. Excludes first stage of vector fitting plot in the results.
* **"spy2"**: _type=bool_. Enables magnitude plot for fitting of $f(s)$ in the results.
* **"logx"**: _type=bool_. Enables logarithmic axis for x in the graphs.
* **"logy"**: _type=bool_. Enables logarithmic axis for y in the graphs.
* **"errplot"**: _type=bool_. Includes relative error graphs in the results.
* **"phaseplot"**: _type=bool_. Excludes plot of phase angle.
* **"legend"**: _type=bool_. Includes legends in plots.

Default options of vector fitting are already defined into vectfit3.py module. There, _opts_ dictionary contains the following configuration:

    opts={
        "symm_mat"   : False, # F(s) samples belong to a full matrix
        "relax"      : True,  # Use vector fitting with relaxed non triviality
        "RMO_data"   : True,  # Matrix elements are organized in RMO into F(s) (asymmetric problems)
        "stable"     : True,  # Enforce stable poles
        "asymp"      : 2,     # Include only D in fitting (not E).
        "skip_pole"  : False, # Do NOT skip pole identification
        "skip_res"   : False, # Do NOT skip residues identification (C,D,E).
        "cmplx_ss"   : True,  # Create complex state space model
        "spy1"       : False, # No plotting for first stage of vector fitting
        "spy2"       : True,  # Create magnitude plot for fitting of f(s)
        "logx"       : True,  # Use logarithmic axis for x
        "logy"       : True,  # Use logarithmic axis for y
        "errplot"    : True,  # Include deviation in magnitude plot
        "phaseplot"  : False, # Exclude plot of phase angle
        "legend"     : True,  # Do include legends in plots
        }

To change any option, _opts_ dictionary can be imported from the module and edited as shown below:

    from vectfit3 import opts
    opts["asymp"]=3        # Modified to include D and E in fitting
    opts["logx"]=False     # Modified to use linear axis for x
    opts["spy2"]=False     # Modified to omit graphs generation into the iterative application of vectfit
    opts["phaseplot"]=True # Modified to include the phase angle graph in the results
    opts["skip_res"]=True  # Modified to skip residue computation during the iterative execution of vector fitting

### _vectfit_ function arguments and output data
The following code line shows a conventional call to _vectfit_ function:

    (SER,poles,rmserr,fit)=vectfit(f,s,poles,weights,opts)

The input arguments are:

* _f_ = $F(s)$ samples of the frequency domain function "1Dim or NDdim" to be fitted. _numpy.array_ of dimensions [Nc x N]. Nc is the number of elements and N is the number of samples
* _s_ = Complex frequency points of evaluation. _numpy.array_ of length [N]
* _poles_ = Initial search poles. _numpy.array_ if length [n]. n is the order of approximation
* _weights_ = priority of each frequency sample in the process. _numpy.array_ of length [N] or shape [Nc x N] according to the weighting strategy used[^4]
* _opts_ = (optional) dictionary with the modifiers to change the default configuration of the algorithm

The function returns a tuple with the following items:

* _SER_ = Dictionary that includes the Space-State model built form the approximation of $F(s)$
* _poles_ = A numpy.array of length [n] with the new poles of $F(s)$
* _rmserr_ = The root mean squared error achieved in the fitting process
* _fit_ = Evaluation of the fitted function that approximates $F(s)$. It has the same shape of _f_

In state-space form the fitted function is expressed as:
![vf_eq2](https://github.com/SebastianL18/vectfit3_for_python/assets/144065519/ecb5c0ca-6dc2-4ebc-ad54-cf75f0573a09)
with A as a [n x n] diagonal and complex matrix containing the poles of $F(s)$, B as the input matrix filled with ones and shaped [n x 1], C as the output matrix of shape [Nc x n], and D and E as 
[Nc x 1] matrixes containing the constants and proportional terms respectively. _vectfit_ internally creates A, B, C, D and E matrixes and puts them in _SER_. if _opts["cmplx_ss"]=False_ then a real only state space model is returned. This modification produces a real A matrix which is block diagonal and contains real poles along its diagonal and complex poles in blocks of shape [2 x 2] with:

    #[[ real(cpole), -imag(cpole) ], [ imag(cpole),  real(cpole) ]]

B is now filled with 1 for real poles and 2 followed by 0 for each complex conjugated pair. values in C are separated in real and imaginary parts and are organized in contiguous pairs into C matrix.

### Iterative implementation of _vectfit_
In order to achieve the best possible fitting _vectfit_ should be used in an iterative manner, because computed poles are always a good set of initial poles for vector fitting algorithm. The next lines of code illustrates how to use _vectfit_ recursively. Note that rms error also may be compared with a tolerance value to set a stop condition into the loop.

    opts["skip_res"]=True  # Modified to skip residue computation during the iterative execution of vector fitting
    opts["spy2"]=False     # Modified to omit graphs generation into the iterative application of vectfit

    Niter=3 #number of iterations
    for itr in range(Niter):
        if itr==Niter-1:
            opts["spy2"]=True       #enabling graphs for the results
            opts["skip_res"]=False  #enabling residue computation in the final iteration
        (SER,poles,rmserr,fit)=vectfit(F,s,poles,weights,opts)
        print("     ...",itr+1," iterations applied")

As shown in the previous example code, skipping the residue computation during the iterative implementation of _vectfit_ is a good practice, because residues of $F(s)$ are obtained only from the best poles computed, and computing time is saved.

## Test cases for _vectfit3.py_

This repository also includes a sample code "_vectfit_testing.py_", with some test cases for _vectfit3.py_ module. Some examples require **pandas** module and external .csv files with input data to run successfully, those data files are also included in this repository. Example cases in_vectfit_testing.py_ are:

 * Fitting of a scalar and artificial frequency domain function
 * 18th order fitting of a two-dimensional frequency response F(s)
 * Fitting of frequency domain samples imported from the measured frequency response of a transformer. **TRANSF_DATA.csv** is needed
 * Element wise approximation of an admitance matrix computed from an equivalent power distribution network. **SYSADMITANCE_DATA.csv** is needed. For more details see[^4]
 * Approximation of a 3x3 propagation matrix of an aerial transmission line. This matrix function corresponds to a single porpagation mode and time delay is already extrated. Requires **MODEH_DATA.csv**

In the following figures vector fitting results for some test cases are shown:

 * Plots obtained from test case 2: Fitting of a two-dimensional function with 18 poles.

![Example2Results](https://github.com/user-attachments/assets/d47ddc39-7102-47b5-a35f-1ee374b61d6f)

 * Plots obtained from test case 4: Element wise approximation of an admittance matrix with 21 elements and using 50 poles.

![Example4Results](https://github.com/user-attachments/assets/a365b6b4-80ab-4981-887e-47683172d549)

* Plots obtained from test case 5: Approximation of a 3x3 propagation matrix of an aerial transmission line.

![Example5Results](https://github.com/user-attachments/assets/972c5d79-3f57-4c17-b98f-1a419aa168c0)

## Changes in this version of the algorithm

 * Differences between _vectfit3.py_ and the original MATLAB implementation are listed below:
    
    - All options for _vectfit3_ configuration are defined as boolean variables, except asymp which has 3 possible states.
    - The new option "symm_mat" for _vectfit3_ configuration is added. This indicates when $F(s)$ samples belong to a symmetric matrix function in CMO. This is a common practice because it reduces the number of elements to fit. See test 4
    - A new option, "RMO_data" is added for vecfit3 configuration. This shows that elements of an asymmetric matrix function are saved in Row Major Order into F(s), set as false for CMO. See test 5
    - A new method to sort the poles computed during the identification process is implemented.
    - Real and complex data is meticulously treated and differentiated through the entire process.
    - General code organization.
    - A new method to compute error plots is implemented. The new error is $log_{10}(error_{relative})$
    - Error graphs are now plotted outside the magnitude axis as a subplot in the same figure.
    - Now "cmplx_ss", "symm_mat" and "RMO_data" flags are also members of SER "dictionary". This is helpful when using the space-state model in external post-processing routines.
    - ss2pr() function is replaced by to buildRES(). The new function only computes residue matrices because vectfit already returns the poles. 
    - The tri2full() function has been renamed to flat2full(). Additionally, it has been modified to handle asymmetric matrix problems.

## Contributions and development

To contribute, give suggestions or report any bug please contact me:
 * Sebastian Loaiza Elejalde
    - _Dsc Student in Power Systems_
    - [_CINVESTAV - GDL_](https://unidad.gdl.cinvestav.mx/)
    - sebloel18@gmail.com 📬
    - sebastian.loaiza@cinvestav.mx 📬

## References
 
This module can be cited through the .cff file attached to the repository.

[^1]: B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.
        
[^2]: B. Gustavsen, "Improving the pole relocating properties of vector fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.
        
[^3]: D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, "Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method", IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008.

[^4]: B. Gustavsen, "User's Guide for vectfit3.m (Fast, Relaxed Vector fitting)", SINTEF Energy Research, N-7465 Trondheim, Norway, 2008. Aviable online: https://www.sintef.no/en/software/vector-fitting/downloads/#menu accesed on: 2/2/2024