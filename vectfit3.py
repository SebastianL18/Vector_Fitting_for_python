# vectfit3.py module.

"""
*** FastRelaxed Vector Fitting for Python v1.1***

vectfit3.py is the implementation of Fast Relaxed Vector Fitting algortihm on python. The original code was written 
in Matlab eviroment. The pourpose of this algorithm is to compute a rational approximation from tabuled data in the 
frequency domain for scalar or vectorized problems. The resulting model can be expressed in either state-space form 
or pole-residue form.

    - Original Matlab code autor: Bjorn Gustavsen (08/2008)
    - Transcripted and adapted by: Sebastian Loaiza (03/2024)
    - Last revision by: Sebastian Loaiza (06/2024)

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
    - A new option, "lowert_mat" is added for vectfit3 configuration. This indicates when F(s) samples belong to a 
      lower triangular matrix function, that reduces the number of elements to fit for a symmetric matrix function.
    - A new option, "RMO_data" is added for vecfit3 configuration. This shows that matrix function elements are saved in 
      Row Major Order into F(s).
    - New options mentioned before and "cmplx_ss" flag are also included in SER.
    - A new method to sort the poles computed during the identification process is implemented.
    - tri2full() function is renamed to flat2full(). It is modified to consider asymmetric matrix problems as well 
      depending on the status of "lower_mat" and "RMO_data" flags.
    - ss2pr() function is replaced by to buildRES(). The new function just compute residues matrixes becouse vectfit returns 
      the poles already.
"""

# Scientific computing modules:
import numpy as np
from scipy.linalg import qr, eigvals, lstsq
from scipy.constants import pi

### ----------------------------------------------------------- Data structures --------------------------------------------------------------  ###

#opts{}: vectfit modifiers dictionary.
# Dictionary which contains the default settings fot vectfit. 
# Any key can be modified to change som options included in vectfit3
opts={
    "symm_mat"   : False, # Indicates when F(s) samples belong to a lower triangular matrix (symmetric problems)
    "RMO_data"   : True,  # Matrix elements are organized in RMO into F(s) (asymmetric problems)
    "relax"      : True,  # Use vector fitting with relaxed non triviality
    "stable"     : True,  # Enforce stable poles
    "asymp"      : 2,     # Include only D in fitting (not E). See [4]
    "skip_pole"  : False, # Do NOT skip pole identification
    "skip_res"   : False, # Do NOT skip residues identification (C,D,E). See [4]
    "cmplx_ss"   : True,  # Create complex state space model
    "spy1"       : False, # No plotting for first stage of vector fitting
    "spy2"       : True,  # Create magnitude plot for fitting of f(s)
    "logx"       : True,  # Use logarithmic axis for x
    "logy"       : True,  # Use logarithmic axis for y
    "errplot"    : True,  # Include deviation in magnitude plot
    "phaseplot"  : False, # Exclude plot of phase angle
    "legend"     : True,  # Do include legends in plots
    }

### ----------------------------------------------------------------- Functions --------------------------------------------------------------- ###

# vectfit() subroutine.
def opts_errorCheck(opts):
    """Fucntion to check any configuration error into opts dictionary.
        *Returns True if an error is found otherwise false"""
    if opts["asymp"]!=1 and opts["asymp"]!=2 and opts["asymp"]!=3:
        print("vectfit3::ERROR::Ilegal value for [asymp] option. It must be 1, 2 or 0")
        return True
    opitems=list(opts.items())
    asymp=("asymp",opts["asymp"])
    opitems.remove(asymp) #asymp already checked
    for option in opitems:
        if not(isinstance(option[1],bool)):
            print("vectfit3::ERROR::Ilegal value for ["+option[0]+"] option. It must be boolean: True or False")
            return True
    return False

# vectfit() subroutine.
def dim_errorCheck(F,s,poles,weights):
    """Function to check dimentions compatibility among vectfit's arguments
        *Returns True if an error is found otherwise false"""
    if len(s.shape)>1 or len(poles.shape)>1:
        print("vectfit3::ERROR::Arguments s and poles must be one-dimentional arrays")
        return True
    if len(F.shape)>1: 
        # Vectorized problem: Nc>1
        if s.size!=F.shape[1]:
            print("vectfit3::ERROR::The number of frequency samples in s and F does not coincide!")
            return True
    else:
        # Scalar problem: Nc=1
        if s.size!=F.size:
            print("vectfit3::ERROR::The number of frequency samples in s and F does not coincide!")
            return True
    if len(weights.shape)>1: 
        # Individual weighting configuration: weights [Nc x N]
        if s.size!=weights.shape[1]:
            print("vectfit3::ERROR::The number of frequency samples in s and elements in weights does not coincide!")
            return True
        if weights.shape[0]!=F.shape[0]:
            print("vectfit3::ERROR::The number of elements in F does not coincide to their weiths for individual weighting")
            return True
    else:
        # Common wighting configuration: weights [1 x N]
        if s.size!=weights.size:
            print("vectfit3::ERROR::The number of frequency samples in s and elements in weights does not coincide!")
            return True
    return False

# vectfit() subroutine.
def identifyPoles(LAMBD):
    """Function to classify poles storaged into the LAMBD diagonal matrix.
       Argument.
        - LAMBD: Diagonal matrix containing the poles of the aproximated model
        *Returns cindex array which marks the poles with the following flags:
         - 0 for real poles
         - 1 and 2 for complex conjugated pairs
    """
    n=LAMBD.shape[0]
    cindex=np.zeros(n, dtype=np.uint16)
    for m in range(n):
        if np.imag(LAMBD[m,m])!=0:
            if m==0:
                cindex[m]=1
            else:
                if cindex[m-1]==0 or cindex[m-1]==2:
                    cindex[m]=1
                    cindex[m+1]=2
                else:
                    cindex[m]=2
    return cindex

# vectfit() subroutine.
def sortPoles(poles):
    """Function to sort the poles obtained in the poles identification process.
       the poles array is sorted by magnitud in asending order, real poles go first and
       for each pair of complex conjugated poles, the one with positive imaginary is placed first
        *Returns a complex array with the required poles arrangement"""      
    poltype=np.abs(poles.imag)<1e-10 # Real poles are marked as True, complex as false
    # *!{it is necesary to use a tolerance instead of np.isreal() due to the numeric error in eigenvalues computation}
    # Real poles extraction
    realpoles=poltype*poles.real                    
    realpoles=np.take_along_axis(realpoles, np.nonzero(realpoles)[0], axis=None) 
    # Complex poles extraction
    complexpoles=np.logical_not(poltype)*poles 
    complexpoles=np.take_along_axis(complexpoles, np.nonzero(complexpoles)[0], axis=None) 
    # Sorting real and complex poles separately: In the following procedure the function sorted() with the list of tuples (abs,pol) 
    #orders by considering first abs, if reapeted magnitudes are found complex poles are ordered considering their imaginary parts.
    realpoles=-1*np.sort(np.abs(realpoles))
    complexpoles=np.array([p for _,p in sorted(zip(np.abs(complexpoles),complexpoles))])
    # Organization of the complex conjugated pairs
    cpolesorder=np.array([], dtype=bool) 
    #desired order [True(+),False(-),...]
    for _ in range(int(complexpoles.size/2)):
        cpolesorder=np.append(cpolesorder,np.array([True,False]))
    polswap=np.logical_xor(complexpoles.imag>0,cpolesorder) #marks the poles that must be swaped as True
    complexpoles=complexpoles-2j*polswap*complexpoles.imag #change of sign for the imaginary parts of the marked poles
    # Organized poles
    return np.append(realpoles,complexpoles)

# vectfit() subroutine.
def vectfitPlot(F,fit,s,opts,initialState=False):
    """Function to plot vector fitting results. 
       The number of graphs that are displayed varies in function of opts configuration
       
       Arguments.
       
         - F: F(s) frequency samples of the original function
         - fit: fit(s) fited function that aproximates F(s)
         - s: Complex frequency points of evaluation for F and fit
         - opts: Dictionary with the modifiers to change the default configuration of Vector Fitting
       
       Results.
        
         - Displays the magnitude plot of F and fit samples with or without error
         - May display the phase angle plot of F and fit samples if enabled
    """
    # Importing plots and graphs module:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    #matplotlib configuration:
    plt.rcParams["font.family"]="Times New Roman" # font type
    mpl.rcParams["font.size"]=9                   # font size in points
    mpl.rcParams["figure.dpi"]=100                # canvas figure resolution
    mpl.rcParams["savefig.dpi"]=300               # saved figure resolution
    
    freq=np.real(s/(2*pi*1j))
    # LogLog plots: Graphs with logarithmic x and y axis
    if opts["errplot"]:
        fig1,ax1=plt.subplots(2,1)
        l1=ax1[0].plot(freq,np.abs(F.T),color='c',linewidth=1.2)
        l2=ax1[0].plot(freq,np.abs(fit.T),color='k',linewidth=1.3,linestyle="dashed")
        ax1[0].set(xlabel="Frequency (Hz)", ylabel="Magnitude",title="Vector Fitting results")
        if opts["legend"]:
            l2[0].set(label="Fitted function")
            l1[0].set(label="F(s) samples")
            ax1[0].legend()
        ax1[0].grid(True)
        #logarithmic error computation
        logError=np.zeros(F.shape,dtype=np.float64)
        for i in range(F.shape[0]):
            logError[i]=np.log10(np.abs(F[i]-fit[i])/np.max(np.abs(F[i])))
        ax1[1].plot(freq,logError.T,color='r',linewidth=1.2)
        ax1[1].set(xlabel="Frequency (Hz)", ylabel="log10(relative error)")
        ax1[1].grid(True)
    else:
        if initialState:
            flabel="Sigma function"
            tlabel="Vector Fitting initial state"
        else:
            flabel="Fitted function"
            tlabel="Vector Fitting results"
        fig1,ax1=plt.subplots()
        l1=ax1.plot(freq,np.abs(F.T),color='c',linewidth=1.2,label="F(s) samples")
        l2=ax1.plot(freq,np.abs(fit.T),color='k',linewidth=1.3,linestyle="dashed")
        ax1.set(xlabel="Frequency (Hz)", ylabel="Magnitude",title=tlabel)
        if opts["legend"]:
            l2[0].set(label=flabel)
            l1[0].set(label="F(s) samples")
            ax1.legend()
        ax1.grid(True)
    if opts["phaseplot"]:
        #fisrt function angles are computed in degrees, and then are unwrapped 
        F_angle=np.unwrap(np.angle(F,deg=True),period=360)
        fit_angle=np.unwrap(np.angle(fit,deg=True),period=360)
        fig2,ax2=plt.subplots()
        l3=ax2.plot(freq,F_angle.T,color='c',linewidth=1.2)
        l4=ax2.plot(freq,fit_angle.T,color='k',linewidth=1.3,linestyle="dashed")
        ax2.set(xlabel="Frequency (Hz)", ylabel="Phase angle (deg)",title="Vector Fitting results")
        if opts["legend"]:
            l3[0].set(label="F(s) samples")
            l4[0].set(label="Fitted function")
            ax2.legend()
        ax2.grid(True)
    if opts["logx"] and opts["logy"]:
        #full logarithmic graphs. Logarithmic x and y axis
        if opts["errplot"]:
            # Magnitude plot:
            ax1[0].set_xscale("log")
            ax1[0].set_yscale("log")
            # Error plot
            ax1[1].set_xscale("log")
        else:
            # Magnitude plot:
            ax1.set_xscale("log")
            ax1.set_yscale("log") 
        if opts["phaseplot"]:
            # Phase plot:
            ax2.set_xscale("log")
    elif opts["logx"]:
        #semilogarithmic graphs. Logarithmic x and linear y axis
        if opts["errplot"]:
            ax1[0].set_xscale("log")
            ax1[1].set_xscale("log")
        else:
            ax1.set_xscale("log")
        if opts["phaseplot"]:
            ax2.set_xscale("log")
    elif opts["logy"]:
        #semilogarithmic graphs. Linear x and logarithmic y axis
        if opts["errplot"]:
            ax1[0].set_yscale("log")
        else:
            ax1.set_yscale("log")
    #else default linear axis configuration of plot() is used
    plt.show()
      
# vectfit3() subroutine.
def buildSER(Ac,Bc,Cc,Dr,Er,complex_format,symmetric_data,RMO_matrixData):
    """Function to build the state-space model of the fitted function such as:
        F(s) = C * (sI-A)^-1 * B + D + sE
    
        Arguments. c subindex indicates complex system
        
         - Ac: Diagonal matrix containing the poles of F(s). Shape [n x n] for n the number fo poles
         - Bc: Input matrix for complex form. Must be a vector of ones. Shape [n x 1]
         - Cc: Residues matrix. Complex and shaped [Nc x n] with Nc as the number of elements in F(s)
         - D: Constant terms matrix. Real matrix of shape [Nc x 1]
         - E: Proportional terms matrix. Real matrix of shape [Nc x 1]
         - complex_format: Boolean option from opts{} which indicates if a complex system is needed
         - symmetric_data: Boolean option from opts{} which indicates if data come from a symmetric matrix problem
         - RMO_matrixData: Boolean option from opts{} which indicates if matrix elements are ordered in RMO 
        
        Results.
        
         - SER: Dictionary that storages A,B,C,D,E system matrixes adapted as indicated by cmplx_ss and symm_mat
    """
    SER=dict(A=Ac,B=Bc,C=Cc,D=Dr,E=Er, cmplx_ss=True, symm_mat=symmetric_data, RMO_data=RMO_matrixData) #complex state-space system
    if not(complex_format):
        # Real state-space system is required so the matrixes are modified
        Ar=np.real(Ac) #importing real poles to Ar
        Br=Bc          #real poles remains realated with ones in B
        Cr=np.real(Cc) #importing real values to Cr
        cindex=identifyPoles(Ac) #poles identification
        # Reorganization of A as a real and block diagonal matrix
        k=0
        for m in range(Ac.shape[0]):
            if cindex[m]==1:
                #complex pole to modify:
                cpole=Ac[k,k]
                re_cpole=np.real(cpole) #real part
                im_cpole=np.imag(cpole) #imaginary part
                # Real system adaptation
                Ar[k:k+2,k:k+2]=np.array([[re_cpole,im_cpole],[-im_cpole,re_cpole]]) #complex pair block in Ar
                # Input matrix modification to consider the complex pair block in Ar
                Br[k,0]=2
                Br[k+1,0]=0
                # Real and imaginary part separation for data in output matrix
                Cr[:,k]=np.real(Cc[:,k])
                Cr[:,k+1]=np.imag(Cc[:,k])
            k+=1
        # Updating SER with new real matrixes:
        SER["A"]=Ar
        SER["B"]=Bc
        SER["C"]=Cr
        SER["cmplx_ss"]=False
    return SER

# build_RES() subroutine.
def flat2full(SER):
    """Function to transform data arrays in SER from its default representation (flattened and element-wise) to a full matrix representation.
        This function works with symmetric or asymmetric matrix functions. If "lower_mat" flag is true, this function considers that the matrix 
        problem need to be decompressed from lower triangular representation. 
                
        *vectfit() function works with flattened element-wise data, therefore if data comes from a matrix problem, this function returns
        a new SER with modified shapes:
            - A shape [n x n] => [Ny*n x Ny*n] for Ny as the original dimentions of the fitted matrix and n as the order of aproximation
            - B shape [n x 1] => [Ny*n x Ny]
            - C shape [Nc x n] => [Ny x Ny*n] for Nc as the number fitted elements into the matrix (organized element-wise)
            - D shape [Nc x 1] => [Ny x Ny]
            - E shape [Nc x 1] => [Ny x Ny]

        Important!: lower trinagular compression of symmetric matrix data need to be applied prior vector fitting application
    """
    A=SER["A"]
    B=SER["B"]
    C=SER["C"]
    D=SER["D"]
    E=SER["E"]
    # Unzip process parameters:
    n=B.shape[0]    # Order of approximation
    Nc=C.shape[0]   # Number of different elements in the flattened problem
    if SER["symm_mat"]:
        # Case for a symmetric problem, data must be decompressed and organized in symmetric format
        sum=0           # Sum counter
        Ny=0            # Shape of the full unzipped matrix function
        while sum<Nc:
            Ny+=1
            sum+=Ny
        # Full versions of the space-state system
        Af=np.zeros((Ny*n,Ny*n), dtype=A.dtype)
        Bf=np.zeros((Ny*n,Ny), dtype=np.float64)
        Cf=np.zeros((Ny,Ny*n), dtype=C.dtype)
        Df=np.zeros((Ny,Ny), dtype=np.float64)
        Ef=np.zeros((Ny,Ny), dtype=np.float64)
        #filling data in full matrixes
        coli=0; k=0 #column and element indicators
        for ind in range(0,Ny*n,n):
            endf=ind+n
            Af[ind:endf,ind:endf]=A
            Bf[ind:endf,coli]=np.ravel(B)
            for rowi in range(coli,Ny):
                Cf[rowi,(coli)*n:(coli+1)*n]=C[k,:]
                Cf[coli,(rowi)*n:(rowi+1)*n]=C[k,:]
                k+=1
            coli+=1
        # Generating an upper trinagular matrix from arrays D and E
        Df[np.triu_indices(Ny)]=D       
        Ef[np.triu_indices(Ny)]=E
        # Converting Df and Ef matrixes to full and symmetric matrixes of shape [Ny x Ny] 
        Df+=Df.T-np.diag(np.diag(Df))   
        Ef+=Ef.T-np.diag(np.diag(Ef))
        SER["symm_mat"]=False
    else:
        # Case for an asymmetric problem. Flattened data must be in CMO or in RMO
        Ny=int(np.sqrt(Nc))
        # Full versions of the space-state system
        Af=np.zeros((Ny*n,Ny*n), dtype=A.dtype)
        Bf=np.zeros((Ny*n,Ny), dtype=np.float64)
        Cf=np.zeros((Ny,Ny*n), dtype=C.dtype)
        Df=np.zeros((Ny,Ny), dtype=np.float64)
        Ef=np.zeros((Ny,Ny), dtype=np.float64)
        # Saving new full matrixes in SER:
        #filling data in full matrixes
        coli=0; k=0 #column and element indicators
        for ind in range(0,Ny*n,n):
            endf=ind+n
            Af[ind:endf,ind:endf]=A
            Bf[ind:endf,coli]=np.ravel(B)
            if SER["RMO_data"]: #data is organized in Row Major Order (default)
                Df[:,coli]=D[coli*Ny:(coli+1)*Ny]
                Ef[:,coli]=E[coli*Ny:(coli+1)*Ny]
                for rowi in range(0,Ny):
                    Cf[rowi,(coli)*n:(coli+1)*n]=C[k,:]
                    k+=1
            else: #data is organized in Column Major Order (transpose filling to iterate colums faster)
                Df[coli,:]=D[coli*Ny:(coli+1)*Ny]
                Ef[coli,:]=E[coli*Ny:(coli+1)*Ny]
                for rowi in range(0,Ny):
                    Cf[coli,(rowi)*n:(rowi+1)*n]=C[k,:]
                    k+=1
            coli+=1
    SER["A"]=Af
    SER["B"]=Bf
    SER["C"]=Cf
    SER["D"]=Df
    SER["E"]=Ef       
    return SER

def buildRES(SERC, SERB):
    """Function to generate residues matrixes of the fitted function computed by vectfit().
        Arguments:
        
        -SERC: Residue matrix C in SER. Generated by vectfit() and modified by flat2full()
        -SERB: Input matrix B in SER. Generated by vectfit() and modified by flat2full()
        
        Important! Arguments must be in full matrix representation instead of flattened and element-wise 
        representation (vectfit default output). Therefore flat2full(SER) must be applyed previously for matrix problems
        
        *To save computational time it is recommended to generate a SER in complex format
        
        *Returns residues matrixes stacked in Res:
            - Res: a 3D array of shape [Ny x Ny x n]. where Ny is the matrix function size and n the aproximation order
    """
    Ny=SERC.shape[0]        #dimentions of the matrix function [Ny x Ny]
    n=int(SERC.shape[1]/Ny) #order of aproximation
    C=SERC                  #residues of Fit(s)
    B=SERB                  #input values (for real SER it distinguish real from complex poles)
    # Data needs to be in complex format. Transformatio is carried out if needed:
    if C.dtype==np.float64: 
        # Real to complex transformation
        Cc=np.zeros(C.shape, dtype=np.complex128)   #complex version of C
        ones=np.ones(n, dtype=np.float64)           #array of ones to compare to values in B
        #realIndex=np.arange(n)[ones==B[0:n,0]]      #real residues' positions
        cmplxIndex=np.arange(n)[ones*2==B[0:n,0]]   #complex residues' positions
        resk=np.zeros(n, dtype=np.complex128)       #transformed residues of each element
        for row in range(Ny):
            for col in range(Ny):
                resk[:]=C[row,col*n:(col+1)*n] #real residues are taken by default
                #then conjugated pairs need to be computed from consecutive real and imaginary parts in real C
                cmplxres=np.take_along_axis(resk, cmplxIndex, axis=None)+1j*np.take_along_axis(resk, cmplxIndex+1, axis=None)
                np.put(resk, cmplxIndex, cmplxres)
                np.put(resk, cmplxIndex+1, np.conj(cmplxres))
                Cc[row,col*n:(col+1)*n]=resk
        C=Cc  
    Res=np.zeros((Ny,Ny,n), dtype=np.complex128) #Residues matrixes
    for row in range(Ny):
        for col in range(Ny):
            Res[row,col,:]=C[row,col*n:(col+1)*n]
    return Res

# * ----------------------------------------------------------  main vectfit3 function ---------------------------------------------------------- *

def vectfit(F,s,poles,weights,opts=opts):
    """ 
    vectfit(): Function to compute a rational aproximation in the frequency domain with the 
    Fast Relaxed Vector Fitting algorithm. Should be used recursively to achieve the best fit.
    
    Arguments. must of the arrays should be dtype np.complex128.
    
        - F: F(s) frequency domain function "1Dim or NDdim" to be fitted of dimentions [Nc x N].
            Nc: number of elements for vector case, otherwise 1.
            N:  number of frequency samples
        - s: array of frequency points of dimentions [1 x N] in [rad/s]
        - poles: array of initial search poles [1 x n] for n as the aproximation order
        - weights: priority of each frequency sample in the process [1 or Nc x N].
            No wighting desired: weight=ones[1 x N]
            Common weighting: weight=array[1 x N]
            Individial elementwise weighting: weight=array[Nc x N]
        - opts: dictionary with the modifiers to change the default configuration of the algorithm
            To be used the opts dictionary should be imported and modified from main module taking
            into account the following value options:
            
                "lowert_mat" <- (bool): Indicates when F(s) samples belong to a lower triangular matrix (symmetric problem)
                "RMO_data"   <- (bool): Indicates when matrix function elements are arranged in RMO into F(s)
                "realx"      <- (bool): Enable vector fitting with relaxed nontriviality
                "stable"     <- (bool): Enable stable poles enforcement
                "asymp"      <- (1, 2 or 3): Produce [D=0; E=0], [D!=0; E=0] or [D!=0; E!=0] respectively [4]
                "skip_pole"  <- (bool): Disable pole identification omission
                "skip_res"   <- (bool): Disable residue identification omission
                "cmplx_ss"   <- (bool): Produce complex or real state space model [4]
                "spy1"       <- (bool): Exclude first stage of vector fitting plot in the results
                "spy2"       <- (bool): Enable magnitude plot for fitting of f(s) in the results
                "logx"       <- (bool): Enable logarithmic axis for x in the graphs
                "logy"       <- (bool): Enable logarithmic axis for y in the graphs
                "errplot"    <- (bool): Include logarithmic error graph in the results
                "phaseplot"  <- (bool): Exclude plot of phase angle
                "legend"     <- (bool): Include legends in plots
    
    Output variables. linked as a tuple
    
        - SER: dictionary containing the Space-State model that aproximates F(s) as:
            F(s) = C * (sI-A)^-1 * B + D + sE
        - poles: array of final poles [1 x n] for n as the aproximation order
        - rmserr: root mean squared error achieved 
        - fit: evaluation of the fitted function that aproximates F(s)
    """
    # Entry errors cheking
    if opts_errorCheck(opts) or dim_errorCheck(F,s,poles,weights):
        print("vecfit() not lunched due to and entry error!")
        return False
    # If no error is found the algorithm computation proceeds
    TOLlow=1e-18; TOLhigh=1e18; # vectfit3 tolerances
    # Initial poles configuration
    if s[0]==0:
        if poles[0]==0 and poles[1]!=0:
            poles[0]=-1
        elif poles[0]!=0 and poles[1]==0:
            poles[1]=-1
        elif poles[0]==0 and poles[1]==0:
            poles[0]=-1+10j
            poles[1]=-1-10j
    N=s.size                                 # Number of samples in the frequency domain
    n=poles.size                             # Order of aproximation
    if len(F.shape)>1:
        Nc=F.shape[0]                        # Number of elements to fit
    else:
        Nc=1
        # Reshaping arrays to ensure a good vectorized computation
        F=np.reshape(F,(1,N))
        
    # Problem arrays declaration
    LAMBD=np.diag(poles)                      # Diagonal matriz with searching poles  
    B=np.ones((n,1), np.float64)              # Column vector of ones
    SERA=poles                                # A Matrix in the space state model
    SERB=np.ones((n,1), np.float64)           # B Matrix in the space state model
    SERC=np.zeros((Nc,n),dtype=np.complex128) # C Matrix in the space state model
    SERD=np.zeros(Nc,dtype=np.float64)        # D Matrix in the space state model
    SERE=np.zeros(Nc,dtype=np.float64)        # E Matrix in the space state model
    fit=np.zeros((Nc,N),dtype=np.complex128)  # Array to store fitted values
    rmserr=-1                                 # Root mean squared error
    # Weighting type identification
    if len(weights.shape)==1:
        commonWeighting=True
        # Reshaping arrays to ensure a good vectorized computation
        #weights=np.reshape(weights,(1,N))
    else:
        commonWeighting=False
    # Space state selected structure
    if opts["asymp"]==1: 
        offs=0 # for [D=0; E=0]
    elif opts["asymp"]==2:
        offs=1 # [D!=0; E=0]
    else:
        offs=2 # [D!=0; E!=0]
        
    # *-- POLES IDENTIFICATION PROCESS
        
    if not(opts["skip_pole"]):
        Escale=np.zeros(Nc+1)
        # Initial complex poles identification: cindex marks 0 for reals, 1 and 2 for complex conjugated pairs 
        cindex=identifyPoles(LAMBD)       
        # Building System Matrix
        Dk=np.zeros((N,n),dtype=np.complex128)
        for m in range(n):
            if cindex[m]==0:
                # real pole defined
                Dk[:,m]=1/(s-LAMBD[m,m])
            elif cindex[m]==1:
                #complex pole first part defined
                Dk[:,m]=1/(s-LAMBD[m,m])+1/(s-np.conj(LAMBD[m,m]))
                Dk[:,m+1]=1j/(s-LAMBD[m,m])-1j/(s-np.conj(LAMBD[m,m]))
        # Depending on the D and E option selected
        if offs==0 or offs==1: 
            # for E=0
            Dk=np.hstack((Dk,np.ones((N,1)))) #note that Dk.shape changed to (N,n+1)
        else:
            # for E!=0
            Dk=np.hstack((Dk,np.ones((N,1))))
            Dk=np.hstack((Dk,np.reshape(s,(N,1)))) #note that Dk.shape changed to (N,n+2)
        # Scaling for last row of LS-Problem
        scale=0
        for m in range(Nc):
            if commonWeighting:
                scale+=np.linalg.norm(weights*F[m,:])**2
            else:
                scale+=np.linalg.norm(weights[m,:]*F[m,:])**2
        scale=np.sqrt(scale)/N
        # Applying relaxed version of the algorithm
        if opts["relax"]:
            AA=np.zeros((Nc*(n+1),n+1),dtype=np.float64)
            bb=np.zeros(Nc*(n+1),dtype=np.float64)
            Escale=np.zeros(n+1,dtype=np.float64)
            offset=n+offs
            for k in range(Nc):
                Ac=np.zeros((N,offset+n+1),dtype=np.complex128)
                if commonWeighting:
                    weig=weights
                else:
                    weig=weights[k,:]
                for m in range(offset): #left block
                    Ac[:,m]=weig*Dk[:,m]
                for m in range(n+1): #right block
                    Ac[:,offset+m]=-weig*Dk[:,m]*F[k,:]
                # Partitioned problem in real and imaginary part, a real array is obtained:
                A=np.vstack((Ac.real,Ac.imag))
                # Integral criterion for sigma
                if k==Nc-1:
                    A=np.vstack((A,np.zeros((1,A.shape[1])))) # adding extra row    
                    for m in range(n+1):
                        A[2*N,offset+m]=np.real(scale*np.sum(Dk[:,m]))
                # Obtaining QR transformation of A
                (Q,R)=qr(A,mode="economic") #routine imported from scipy module
                R22=R[offset:offset+n+1,offset:offset+n+1]
                AA[k*(n+1):(k+1)*(n+1),:]=R22
                if k==Nc-1:
                    bb[k*(n+1):(k+1)*(n+1)]=Q[-1,offset:]*N*scale
            for m in range(AA.shape[1]):
                Escale[m]=1/np.linalg.norm(AA[:,m])
                AA[:,m]=Escale[m]*AA[:,m]
            # LS solution routine imported from scipy module: 
            #   - check_finite option is disabled to improve performance. NaN should not appear into the arrays
            #   - gelsy lapack driver is chosen because is slightly faster than default ("gelsd")
            # The routine returns a tuple (solution,residues,rank,svalues), however just the solution is taken: 
            x=lstsq(AA,bb, check_finite=False, lapack_driver="gelsy")[0] # solution <- result[0]
            x=x*Escale
        # ...end of opts["relax"]=True segment
        # Case for non relaxion version of the algorithm.
        # Also may be needed when D of sigma results extremely small and large. I needs to be solved again.  
        if not(opts["relax"]) or np.abs(x[-1])<TOLlow or np.abs(x[-1])>TOLhigh:
            AA=np.zeros((Nc*n,n),dtype=np.float64)
            bb=np.zeros(Nc*n,dtype=np.float64)
            if opts["relax"]:
                Dnew=1
            else:
                if x[-1]==0:
                    Dnew=1
                elif np.abs(x[-1])<TOLlow:
                    Dnew=np.sign(x[-1])*TOLlow
                elif np.abs(x[-1])>TOLhigh:
                    Dnew=np.sign(x[-1])*TOLhigh
            offset=n+offs
            for k in range(Nc):
                Ac=np.zeros((N,offset+n),dtype=np.complex128)
                Escale=np.zeros(n,dtype=np.float64)
                if commonWeighting:
                    weig=weights
                else:
                    weig=weights[k,:]
                for m in range(offset): #left block
                    Ac[:,m]=weig*Dk[:,m]
                for m in range(n): #right block
                    Ac[:,offset+m]=-weig*Dk[:,m]*F[k,:]
                bc=Dnew*weig*F[k,:]
                # Partitioned problem in real and imaginary part, real arrays are obtained:
                A=np.vstack((Ac.real,Ac.imag))
                b=np.append(bc.real,bc.imag)
                # Obtaining QR transformation of A
                (Q,R)=qr(A,mode="economic")
                R22=R[offset:offset+n,offset:offset+n]
                AA[k*n:(k+1)*n,:]=R22
                # The following arrays must be reshaped in order to perform the matrix product
                bb[k*n:(k+1)*n]=np.ravel(np.transpose(Q[:,offset:offset+n])@np.reshape(b,(b.size,1)))
            for m in range(AA.shape[1]):
                Escale[m]=1/np.linalg.norm(AA[:,m])
                AA[:,m]=Escale[m]*AA[:,m]
            # LS solution routine imported from scipy module: 
            #   - check_finite option is disabled to improve performance. NaN should not appear into the arrays
            #   - gelsy lapack driver is chosen because is slightly faster than default ("gelsd")
            # The routine returns a tuple (solution,residues,rank,svalues), however just the solution is taken: 
            x=lstsq(AA,bb, check_finite=False, lapack_driver="gelsy")[0] # solution <- result[0]
            x=x*Escale
            x=np.append(x,Dnew)
        # ...end of opts["relax"]=False or out of tolerance segment
        # Changing back to make C complex:    
        C=np.zeros(x.size-1,dtype=np.complex128)    
        C[0:C.size]=x[0:-1]
        D=x[-1]
        for m in range(n):
            if cindex[m]==1:
                r1=C[m] #real part of the complex number
                r2=C[m+1] #imaginary part of the complex number
                #building the conjugated complex pair
                C[m]=r1+1j*r2
                C[m+1]=r1-1j*r2
        # Graphs of initial stage of vector fitting process
        if opts["spy1"]:
            print("\nvectfit3::spy1_Enabled::Building and showing graph for initial fitting state...")
            # First fitting state evaluation: 
            # As mentioned in [1] pole identification begins with the aproximation of sigma(s). sigma(s) is an unknown function whose  
            #approximation has the same poles of F(s). Furtheremore it is formulated that (sigma*f)_fit(s) = sigma_fit*f(s), whence can be
            #demostrated that zeros of sigma are a better set of poles to fit f(s), therefore sigma(s) is called the initial fitting state.
            Dk=np.zeros((N,n),dtype=np.complex128)
            for m in range(n):
                Dk[:,m]=1/(s-LAMBD[m,m])
            sigma=D+Dk@C
            #setting temporal options to plot initial vector fitting state
            opts_temp=opts
            opts_temp["errplot"]=False
            opts_temp["phaseplot"]=False
            vectfitPlot(F,sigma,s,opts_temp,True)
        # Calculating the zeros for sigma
        m=0 #auxiliar counter
        for k in range(n):
            if m<n:
                if np.abs(LAMBD[m,m])>np.abs(np.real(LAMBD[m,m])):
                    # Complex poles clasification: diagonal blocks of the form: 
                    # [ real, -imag ]
                    # [ imag,  real ]
                    LAMBD[m+1,m]=-np.imag(LAMBD[m,m])
                    LAMBD[m,m+1]=np.imag(LAMBD[m,m])
                    LAMBD[m,m]=np.real(LAMBD[m,m])
                    LAMBD[m+1,m+1]=LAMBD[m,m]
                    # B vector modification for  complex poles
                    B[m,0]=2
                    B[m+1,0]=0
                    ccval=C[m]
                    # Complex values clasification:
                    C[m]=np.real(ccval)
                    C[m+1]=np.imag(ccval)
                    m+=1
            m+=1
        #to perform the vectorized B*C and obtain a matrix C need to be reshaped
        C=np.reshape(C,(1,C.size))
        #also just the real parts of LAMBD and C are taken in order to obtain a real[float64] matrix for ZER
        ZER=LAMBD.real-B@C.real/D
        # Computation of ZER eigenvalues
        poles=eigvals(ZER) #routine imported from scipy module
        # Unstabla values identification
        unstable=poles.real>0 #generates a logical array
        if opts["stable"]:
            if np.any(unstable):
                #the product of roetter and unstable extracts the unstable poles
                extracted=poles*unstable
                poles=poles-2*extracted.real
        poles=sortPoles(poles)
        SERA=poles
    #...end of poles identification process
        
    # --* RESIDUES IDENTIFICATION PROCESS
    
    if not(opts["skip_res"]):
        # Now SER for f is calculated by using modified zeros of sigma as new poles:
        LAMBD=np.diag(poles)
        # Initial complex poles identification: cindex marks 0 for reals, 1 and 2 for complex conjugated pairs 
        cindex=identifyPoles(LAMBD)
        # Building System Matrix
        Dk=np.zeros((N,n),dtype=np.complex128)
        for m in range(n):
            if cindex[m]==0:
                # real pole defined
                Dk[:,m]=1/(s-LAMBD[m,m])
            elif cindex[m]==1:
                #complex pole first part defined
                Dk[:,m]=1/(s-LAMBD[m,m])+1/(s-np.conj(LAMBD[m,m]))
                Dk[:,m+1]=1j/(s-LAMBD[m,m])-1j/(s-np.conj(LAMBD[m,m]))
        if commonWeighting: #case for common wighting
            C=np.zeros((Nc,n),dtype=np.complex128)
            for m in range(n): #same weight for all frequency samples
                Dk[:,m]=weights*Dk[:,m]
            # The SER for the new fitting is calculated by using the calculated zeros as new poles
            if opts["asymp"]==1:
                A=np.zeros((2*N,n),dtype=np.float64)
                A[0:N,:]=Dk.real
                A[N:2*N,:]=Dk.imag
            elif opts["asymp"]==2:
                A=np.zeros((2*N,n+1),dtype=np.float64)
                A[0:N,0:n]=Dk.real
                A[N:2*N,0:n]=Dk.imag
                A[0:N,n]=weights
            else:
                A=np.zeros((2*N,n+2),dtype=np.float64)
                A[0:N,0:n]=Dk.real
                A[N:2*N,0:n]=Dk.imag
                A[0:N,n]=weights
                A[N:2*N,n+1]=np.imag(weights*s)
            BB=np.zeros((2*N,Nc),dtype=np.float64)
            BBc=np.zeros((N,Nc),dtype=np.complex128)
            for m in range(Nc):
                BBc[:,m]=weights*F[m,:] #complex values for BB
            BB[0:N,:]=BBc.real
            BB[N:2*N,:]=BBc.imag
            Escale=np.zeros(A.shape[1],dtype=np.float64)
            for m in range(A.shape[1]):
                Escale[m]=np.linalg.norm(A[:,m])
                A[:,m]=A[:,m]/Escale[m]
            # LS solution routine imported from scipy module: 
            #   - check_finite option is disabled to improve performance. NaN should not appear into the arrays
            #   - gelsy lapack driver is chosen because is slightly faster than default ("gelsd")
            # The routine returns a tuple (solution,residues,rank,svalues), however just the solution is taken:   
            x=lstsq(A,BB, check_finite=False, lapack_driver="gelsy")[0] # solution <- result[0]
            for m in range(Nc):
                x[:,m]=x[:,m]/Escale
            x=np.transpose(x)
            C[:,0:n]=x[:,0:n]
            if opts["asymp"]==2:
                SERD=x[:,n]
            elif opts["asymp"]==3:
                SERD=x[:,n]
                SERE=x[:,n+1]
        else: #no common wighting used    
            C=np.zeros((Nc,n),dtype=np.complex128)
            for k in range(Nc):
                # The SER for the new fitting is calculated by using the calculated zeros as new poles
                if opts["asymp"]==1:
                    A=np.zeros((2*N,n),dtype=np.float64)
                    A[0:N,:]=Dk.real
                    A[N:2*N,:]=Dk.imag
                elif opts["asymp"]==2:
                    A=np.zeros((2*N,n+1),dtype=np.float64)
                    A[0:N,0:n]=Dk.real
                    A[N:2*N,0:n]=Dk.imag
                    A[0:N,n]=1
                else: #for asymp==3
                    A=np.zeros((2*N,n+2),dtype=np.float64)
                    A[0:N,0:n]=Dk.real
                    A[N:2*N,0:n]=Dk.imag
                    A[0:N,n]=1
                    A[N:2*N,n+1]=np.imag(s)
                for m in range(A.shape[1]):
                    A[0:N,m]=weights[k,:]*A[0:N,m]
                    A[N:2*N,m]=weights[k,:]*A[N:2*N,m]
                BB=np.zeros(2*N,dtype=np.float64)
                BBc=np.zeros(N,dtype=np.complex128) #complex values for BB
                BBc=weights[k,:]*F[k,:]
                BB[0:N]=BBc.real
                BB[N:2*N]=BBc.imag
                Escale=np.zeros(A.shape[1],dtype=np.float64)
                for m in range(A.shape[1]):
                    Escale[m]=np.linalg.norm(A[:,m])
                    A[:,m]=A[:,m]/Escale[m]
                # LS solution routine imported from scipy module: 
                #   - check_finite option is disabled to improve performance. NaN should not appear into the arrays
                #   - gelsy lapack driver is chosen because is slightly faster than default ("gelsd")
                # The routine returns a tuple (solution,residues,rank,svalues), however just the solution is taken: 
                x=lstsq(A,BB, check_finite=False, lapack_driver="gelsy")[0] # solution <- result[0]
                x/=Escale
                C[k,:]=x[0:n]
                if opts["asymp"]==2:
                    SERD[k]=x[n]
                elif opts["asymp"]==3:
                    SERD[k]=x[n]
                    SERE[k]=x[n+1]
        #...end of wighting options for residue computation
        # Changing back to make C complex:
        for m in range(n):
            if cindex[m]==1:
                for k in range(Nc):
                    r1=C[k,m] #real part of the complex number
                    r2=C[k,m+1] #imaginary part of the complex number
                    #bulding the conjugated complex pair:
                    C[k,m]=r1+1j*r2 
                    C[k,m+1]=r1-1j*r2
        # New fitting evaluation:
        SERA=np.diag(LAMBD)   
        SERC=C
        for m in range(n):
            Dk[:,m]=1/(s-SERA[m])
        for k in range(Nc):
            fit[k,:]=Dk@SERC[k,:]
            if opts["asymp"]==2:
                fit[k,:]=fit[k,:]+SERD[k]
            elif opts["asymp"]==3:
                fit[k,:]=fit[k,:]+SERD[k]+s*SERE[k]
        # Root mean squared error computation:
        diff=fit-F #diferrences between samples and the fitted function
        rmserr=np.sqrt(np.sum(np.sum(np.abs(diff**2))))/np.sqrt(Nc*N)
        
        # - Graphs generation for vector fitting results:
        if opts["spy2"]:
            print("\nvectfit3::spy2_Enabled::Building and showing graphs for the results...")
            vectfitPlot(F,fit,s,opts)
    #...end of residue identification process
    
    # - Building the state-space model
    A=np.diag(SERA)
    B=SERB
    C=SERC
    D=SERD
    E=SERE
    SER=buildSER(A,B,C,D,E,opts["cmplx_ss"],opts["symm_mat"], opts["RMO_data"])
    
    # Vector fitting process finished.
    return (SER,poles,rmserr,fit)