### Sample code to test Vector Fitting algorthm implemented in Python by Sebastian Loaiza ###

from vectfit3 import vecfit        # Vector Fitting algorithm imported from its external module vectfit3.py
import numpy as np
import pandas as pd
from scipy.constants import pi

# Example to be applied:
# 1. Scalar and artificial frequency domain function f(s)
# 2. 18th order frequency response F(s) of two dimentions
# 3. Escalar measured response of a transformer 
# 4. Elementwise aproximation of a 6x6 admitance matrix Y(s)

test=4 # <- Desired test

if test==1:
    print("Test 1: Scalar and artificial frequency domain function f(s)") # -------------------------------------------------------------------- #

    N=101                                              # Number of frequency samples
    s=2j*pi*np.logspace(0,4,N,dtype=np.complex128)     # Samples in the frequency domain
    f=np.zeros(N,dtype=np.complex128)                  # Function samples in the frequency domain
    weights=np.ones(N,dtype=np.float64)                # Common and unitary weiths for all frequency samples
    
    # Generating test function values
    n=0       
    for sn in s:
        f[n]=2/(sn+5)+(30+40j)/(sn-(-100+500j))+(30-40j)/(sn-(-100-500j))+0.5
        n+=1
    
    print("\nFrequency domain samples of f(s) = \n",f)
    print("f(s) shape = ",f.shape, "\ndata type in f(s) = ",type(f),type(f[0]))

    n=3                                                # Order of aproximation
    poles=-2*pi*np.logspace(0,4,n,dtype=np.complex128) # Initial searching poles
    
    print("\nInitial searching poles:\n",poles)

    # vector fitting configuration
    from vectfit3 import opts
    opts["asymp"]=3        # Modified to include D and E in fitting
    opts["phaseplot"]=True # Modified to include the phase angle graph
    
    # Remaining options by default

    print("\n * Applying vector fitting...")
    (SER,poles,rmserr,fit)=vecfit(f,s,poles,weights,opts)
    print(" v/ Fitting process completed. Aproximation error achieved = ",rmserr)
    print("\nFinal poles computed:\n",poles)
    
elif test==2:
    print("Test 2: 18th order frequency response F(s) of two dimentions") # -------------------------------------------------------------------- #
    
    w=2*pi*np.linspace(1,1e5,100,dtype=np.complex128) # Frequency points of evaluation
    s=1j*w                                            # Complex frequency points of evaluation
    N=s.size                                          # Samples in the frequency domain
    weights=np.ones(N,dtype=np.float64)               # Common and unitary weiths for all frequency samples

    # Test function generation from its SER representation:
    #Sample poles for the elements in F(s)
    p=np.array([-4500,-41000,-100+5e3j,-100-5e3j,-120+15e3j,-120-15e3j,-3e3+35e3j,-3e3-35e3j],dtype=np.complex128)
    p=np.append(p,np.array([-200+45e3j,-200-45e3j,-1500+45e3j,-1500-45e3j],dtype=np.complex128))
    p=np.append(p,np.array([-5e2+70e3j,-5e2-70e3j,-1e3+73e3j,-1e3-73e3j,-2e3+90e3j,-2e3-90e3j],dtype=np.complex128))
    p=2*pi*p
    N1=10
    p1=p[:N1]   #First element poles
    p2=p[N1-2:] #Second element poles
    #Sample residues for the elements in F(s)
    r=np.array([-3000,-83000,-5+7e3j,-5-7e3j,-20+18e3j,-20-18e3j,6e3+45e3j,6e3-45e3j],dtype=np.complex128)
    r=np.append(r,np.array([40+60e3j,40-60e3j,90+10e3j,90-10e3j],dtype=np.complex128))
    r=np.append(r,np.array([5e4+80e3j,5e4-80e3j,1e3+45e3j,1e3-45e3j,-5e3+92e3j,-5e3-92e3j],dtype=np.complex128))
    r=2*pi*r
    r1=r[:N1]   #First element residues
    r2=r[N1-2:] #Second element residues
    # SER constant values:
    D=0.2
    E=2e-5
    # Generating test function values
    F=np.zeros((2,N),dtype=np.complex128)
    n=0
    for sn in s:
        for k in range(N1):
            F[0,n]+=r1[k]/(sn-p1[k]) #evaluation for the first element
            F[1,n]+=r2[k]/(sn-p2[k]) #evaluation for the second element
        F[0,n]+=sn*E
        F[1,n]+=sn*3*E
        n+=1
    F[0,:]=F[0,:]+3*D
    
    print("\nFrequency domain samples of f(s) = \n",F)
    print("f(s) shape = ",F.shape, "\ndata type in f(s) = ",type(F),type(F[0][0]))
    
    # Initial complex poles:
    n=18 # Order of aproximation
    # Initial
    Bet=np.linspace(w[0],w[N-1],int(n/2))
    poles=np.zeros(n,dtype=np.complex128)
    #setting poles as complex conjugated pairs
    for k in range(int(n/2)):
        alf=-Bet[k]*1e-2
        poles[2*k]=alf-1j*Bet[k]
        poles[2*k+1]=alf+1j*Bet[k]
    
    print("\nInitial searching poles:\n",poles)
    
    # vector fitting configuration
    from vectfit3 import opts
    opts["asymp"]=3        # Modified to include D and E in fitting
    opts["skip_res"]=True  # Modified to skip residue computation in the first iteration
    opts["cmplx_ss"]=False # Modified to build a real-only state space model
    opts["logx"]=False     # Modified to use linear axis for x
    opts["spy2"]=False     # Modified to omit graphs generation into the iterative application of vectfit
    opts["phaseplot"]=True # Modified to include the phase angle graph in the results
    # Remaining options by default
    
    print("\n * Applying 3 iterations of vector fitting...")
    Niter=3
    for itr in range(Niter):
        if itr==Niter-1:
            opts["skip_res"]=False #residue computation in final iteration
            opts["spy2"]=True      #enabling graphs for the results
        (SER,poles,rmserr,fit)=vecfit(F,s,poles,weights,opts)
        print("     ...",itr+1," iterations applied")
    print(" v/ Fitting process completed. Aproximation error achieved = ",rmserr)
    print("\nFinal poles computed:\n",poles)

elif test==3:
    print("Test 3: Escalar measured response of a transformer") # ------------------------------------------------------------------------------ #
    # Importing measured data from a .csv file (TRANSF_DATA.csv)
    csv_path=r"C:\Users\Sebastian\Documents\GitHub\Vector_Fitting_for_python\TRANSF_DATA.csv" #local path! Update with your own path
    Mdata=pd.read_csv(csv_path)                         # Measured data
    Mdata=Mdata.to_numpy()                              # Data transfered into a numpy array
    f=Mdata[:160,0]*np.exp(1j*Mdata[:160,1]*pi/180)     # f(s) samples computation. *Just 160 samples are taken
    
    print("\nFrequency domain samples of f(s) = \n",f)
    print("f(s) shape = ",f.shape, "\ndata type in f(s) = ",type(f),type(f[0]))
    
    N=f.size                                           # Number of frequency samples
    w=2*pi*np.linspace(0,10e6,401)[1:161]              # Frequency evaluation for the range of interest
    s=1j*w                                             # Samples in the frequency domain
    
    print("\nFor the first attempt 6 complex starting poles and no weighting are used")
    
    weights=np.ones(N,dtype=np.float64)                # Common and unitary weiths for all frequency samples
    n=6 # Order of approximation
    # Starting poles generation:
    Bet=np.linspace(w[0],w[N-1],int(n/2))
    poles=np.zeros(n,dtype=np.complex128)
    #setting poles as complex conjugated pairs
    for k in range(int(n/2)):
        alf=-Bet[k]*1e-2
        poles[2*k]=alf-1j*Bet[k]
        poles[2*k+1]=alf+1j*Bet[k]
    
    print("\nInitial searching poles:\n",poles)

    # vector fitting configuration
    from vectfit3 import opts
    opts["asymp"]=3        # Modified to include D and E in fitting
    opts["logx"]=False     # Modified to use linear axis for x
    opts["spy2"]=False     # Modified to omit graphs generation into the iterative application of vectfit
    opts["phaseplot"]=True # Modified to include the phase angle graph in the results
    # Remaining options by default

    print("\n * Applying 5 iterations of vector fitting...")
    Niter=5
    for itr in range(Niter):
        if itr==Niter-1:
            opts["spy2"]=True      #enabling graphs for the results
        (SER,poles,rmserr,fit)=vecfit(f,s,poles,weights,opts)
        print("     ...",itr+1," iterations applied")
    print(" v/ Fitting process completed. Aproximation error achieved = ",rmserr)
    print("\nFinal poles computed:\n",poles)

    print("\nFor better results second attempt is carried out with 30 complex starting poles and inverse weighting")

    weights=1/np.abs(f) # Inverse weighting for frequency samples
    n=30 # Order of approximation
    # Starting poles generation:
    Bet=np.linspace(w[0],w[N-1],int(n/2))
    poles=np.zeros(n,dtype=np.complex128)
    #setting poles as complex conjugated pairs
    for k in range(int(n/2)):
        alf=-Bet[k]*1e-2
        poles[2*k]=alf-1j*Bet[k]
        poles[2*k+1]=alf+1j*Bet[k]
    
    print("\nInitial searching poles:\n",poles)

    opts["spy2"]=False     # Modified to omit graphs generation into the iterative application of vectfit
    # Remaining options as configured before

    print("\n * Applying 5 iterations of vector fitting...")
    Niter=5
    for itr in range(Niter):
        if itr==Niter-1:
            opts["spy2"]=True      #enabling graphs for the results
        (SER,poles,rmserr,fit)=vecfit(f,s,poles,weights,opts)
        print("     ...",itr+1," iterations applied")
    print(" v/ Fitting process completed. Aproximation error achieved = ",rmserr)
    print("\nFinal poles computed:\n",poles)
    
elif test==4:
    print("Test 4: Elementwise aproximation of a 6x6 admitance matrix") # --------------------------------------------------------------------------------- #
    # Importing measured data from a .csv file (SYSADMITANCE_DATA.csv)
    csv_path=r"C:\Users\Sebastian\Documents\GitHub\Vector_Fitting_for_python\SYSADMITANCE_DATA.csv" #local path! Update with your own path
    Mdata=pd.read_csv(csv_path)                         # Measured data
    Mdata=np.ravel(Mdata.to_numpy())                    # Data transfered into a numpy array
    N=int(Mdata[0])                                     # Number of frequency samples
    Ysys=np.zeros((6,6,N),dtype=np.complex128)          # Samples of matrix Y(s) imported from the file
    s=np.zeros(N,dtype=np.complex128)                   # Complex frequency points of evaluation for Y(s)
    
    # Y(s) and s are compressed in row-major order into Mdata array. Samples are organized in blocks of 73 numbers, inside each block the first
    #number matchs the complex frequency samples and the remaining ones corresponds to each element in the admitance matrix Y(s). first comes 
    #the real part and then the imaginary part.
    k=0 #sample index
    for i in range(1,Mdata.size,73):
        s[k]=1j*Mdata[i] #complex frequency sample
        for row in range(6):
            ind=row*12+i
            Ysys[row,:,k]=Mdata[ind+1:ind+12:2]+1j*Mdata[ind+2:ind+13:2] #entire row is append to Y(s)
        k+=1
    # Stacking Y(s) data as elements of a frequency domain function F(s). 
    # Due to Y(s) is symmetric only the members into the lower trinagular submatrix Y(s) are taken    
    F=np.zeros((21,N), dtype=np.complex128)
    k=0 #element index
    for col in range(6):
        for row in range(col,6):
            F[k,:]=Ysys[row,col,:] #all samples are in z axis
            k+=1
            
    print("\nFrequency domain samples of f(s) = \n",F)
    print("f(s) shape = ",F.shape, "\ndata type in f(s) = ",type(F),type(F[0,0]))
    
    weights=1/np.sqrt(np.abs(F)) # Weighting with inverse of the square root of the magnitude of F(s)
    n=50                         # Order of approximation
    w=s.imag                     # Angular frequency samples in rad/s
    # Starting poles generation:
    Bet=np.linspace(w[0],w[N-1],int(n/2))
    poles=np.zeros(n,dtype=np.complex128)
    #setting poles as complex conjugated pairs
    for k in range(int(n/2)):
        alf=-Bet[k]/100
        poles[2*k]=alf-1j*Bet[k]
        poles[2*k+1]=alf+1j*Bet[k]
    
    print("A better set of initial poles are obtained by fitting the weighted column sum of the first column of Y(s)\n * Applying 5 iterations of vector fitting...")
    # These weights are common for all column elements
    g=np.zeros(N, dtype=np.complex128)
    for k in range(6):
       g=g+F[k,:]/np.linalg.norm(F[k,:])
    weights_g=1/np.abs(g)
  
    # vector fitting configuration
    from vectfit3 import opts
    opts["asymp"]=3        # Modified to include D and E in fitting
    opts["logx"]=False     # Modified to use linear axis for x
    opts["spy2"]=False     # Modified to omit graphs generation into the iterative application of vectfit
    opts["phaseplot"]=True # Modified to include the phase angle graph in the results
    opts["skip_res"]=True  # Modified to skip residue computation during the iterative execution of vector fitting
    # Remaining options by default
    
    print("\n * Applying 5 iterations of vector fitting...")
    Niter=5 #number of iterations
    for itr in range(Niter):
        if itr==Niter-1:
            #opts["spy2"]=True       #enabling graphs for the results
            opts["skip_res"]=False  #enabling residue computation in the final iteration
        (SER,poles,rmserr,fit)=vecfit(g,s,poles,weights_g,opts)
        print("     ...",itr+1," iterations applied")
    print(" v/ Fitting process completed. Aproximation error achieved = ",rmserr)
    print("\nInitial poles computed from weighted column sum of Y(s):\n",poles)
    
    # Final fitting with new initial poles set and all elementos of Y(s) in F(s)
    opts["skip_res"]=True  # Modified to skip residue computation during the iterative execution of vector fitting
    opts["spy2"]=False     # Modified to omit graphs generation into the iterative application of vectfit

    Niter=3 #number of iterations
    for itr in range(Niter):
        if itr==Niter-1:
            opts["spy2"]=True       #enabling graphs for the results
            opts["skip_res"]=False  #enabling residue computation in the final iteration
        (SER,poles,rmserr,fit)=vecfit(F,s,poles,weights,opts)
        print("     ...",itr+1," iterations applied")
    print(" v/ Fitting process completed. Aproximation error achieved = ",rmserr)
    print("\nFinal poles computed:\n",poles)