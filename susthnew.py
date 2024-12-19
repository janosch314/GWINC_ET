# This scripts translates into python the Advanced Virgo suspension thermal noise script
# and integrates it into the Einstein Telescope noise-budget. The calculations are given
# in the Virgo Technical Document 'VIR-015A-09' (PPP effect).


import os
import numpy as np
from numpy import pi, sqrt, arctan, sin, cos, roots, size, real, imag, sqrt, exp
import const


# exterior definition of propagated functions for simple calculations
# fibre tension, cross-sectional area and momenta
def tension_per_fibre(mass, fiberNb, g):
    tension = mass*g/fiberNb
    return tension
def cross_section_area(r):
    S = np.pi*(r**2)
    return S
def cross_section_momentum(r):
    I = (np.pi/4)*r**4
    return I

def STNpy(f,ifo):
    WireMat0=ifo.Suspension.Stage[0].WireMaterial
    WireMat1=ifo.Suspension.Stage[1].WireMaterial
    WireMat2=ifo.Suspension.Stage[2].WireMaterial


    g=const.g
    kb=const.kB
    Temp=ifo.Suspension.Temp

# Vertical to horizontal coupling
    thetaVH=ifo.Infrastructure.Length/ const.R_earth
    w=2*pi*f


#____Mirror pendulum stage
#____Fused Silica thermo-mechanical properties
    alpha=ifo.Suspension[WireMat0].Alpha
    beta=ifo.Suspension[WireMat0].dlnEdT
    C=ifo.Suspension[WireMat0].C
    K=ifo.Suspension[WireMat0].K
    rho=ifo.Suspension[WireMat0].Rho
    phi0=ifo.Suspension[WireMat0].Phi
    ds=ifo.Suspension[WireMat0].Dissdepth
    E0=ifo.Suspension[WireMat0].Y

    Optcyl=ifo.Suspension.Optcyl #Switch cylindrical/optimized fibres

    m=ifo.Suspension.Stage[0].Mass
    N=ifo.Suspension.Stage[0].NWires


    T=tension_per_fibre(m,N,g)

    if Optcyl==0:
        L=ifo.Suspension.Stage[0].Length
        d=ifo.Suspension.Stage[0].WireDiameter
        r=d/2
        II=cross_section_momentum(r)
    else:
        L1=ifo.Suspension.Stage[0].Length1
        L2=ifo.Suspension.Stage[0].Length2
        L3=ifo.Suspension.Stage[0].Length3
        d1=ifo.Suspension.Stage[0].WireDiameter1
        d2=ifo.Suspension.Stage[0].WireDiameter2
        d3=ifo.Suspension.Stage[0].WireDiameter3
        r1=d1/2
        r2=d2/2
        r3=d3/2
        L=L1+L2+L3

    delta=1e-3
    def surface_loss_angle(radius, dissipation_depth):
        phiST = phi0*(1+(4/radius)*dissipation_depth)
        #phi0=material intrinsic loss angle
        return phiST
    def dissipation_delta(YM, alpha, beta, mass, fiberNb, radius, Temp, rho, C ):
        T = tension_per_fibre(mass, fiberNb, g)
        S = cross_section_area(radius)
        Delta = YM*(alpha - beta*T/(S*YM))**2*Temp/(rho*C)
        return Delta

    z=np.zeros([len(f),3000],dtype = "complex")
    x=np.zeros([len(f),3000],dtype = "complex")
    TotEnergy=np.zeros([len(f)],dtype = "complex")
    Energy=np.zeros([len(f),3000],dtype = "complex")
    phiST=np.zeros([len(f)],dtype = "complex")
    phithT=np.zeros([len(f)],dtype = "complex")

    for i in range(len(w)):
        if Optcyl==0: # Cylindrical fibres
    #Surface contribution for circular cross section fibres (horizontal
    #motion)
            mu=4/r
            phiST[i]=surface_loss_angle(r,ds)
            S=cross_section_area(r)
    #Thermoelastic contribution
            Delta=dissipation_delta(E0, alpha, beta, m, N, r, Temp, rho, C)
            tau=0.0737*rho*C*(2*r)**2/K
            phithT[i]=Delta*w[i]*tau/(1+w[i]**2*tau**2)
            S0=S
        else:  # Optimized fibres
    # Optimized fibres routine is called
            [ZZ,XX,dummy,ENER]=Optimisedfibres(f[i],0,ifo)
            z[i,:]=np.concatenate(ZZ)
            x[i,:]=np.concatenate(XX)
            Energy[i,:]=np.concatenate(ENER)
            TotEnergy[i]=np.sum(Energy[i,:])
            phiST[i]=0
            phithT[i]=0

            for j in range(len(z[1,:])):
                if z[i,j]<L1:
                    r=r1
                    #break
                elif z[i,j]<L1+L2:
                    r=r2
                    #break
                else:
                    r=r3
                    #break
                mu=4/r
                phiS=phi0*(1+mu*ds)
                S=cross_section_area(r)
                Delta=dissipation_delta(E0, alpha, beta, m, N, r, Temp, rho, C)
                tau=0.0737*rho*C*(2*r)**2/K
            # Surface contribution and thermoelastic contribution are weighted by
            # the bending energy
                phiST[i]=phiST[i]+phiS*Energy[i,j]
                phithT[i]=phithT[i]+Delta*w[i]*tau/(1+w[i]**2*tau**2)*Energy[i,j]

            phiST[i]=phiST[i]/TotEnergy[i];
            phithT[i]=phithT[i]/TotEnergy[i];
            S0=pi*(r1**2*L1+r2**2*L2+r3**2*L3)/L;
        phiT=phiST+phithT

#--------------------- Marionette ---------------------

# Marionette and RM stages parameter are loaded.
    m1=ifo.Suspension.Stage[2].Mass
    LL1=ifo.Suspension.Stage[2].Length
    dd1=ifo.Suspension.Stage[2].WireDiameter
    N1=ifo.Suspension.Stage[2].NWires
    rho1=ifo.Suspension[WireMat2].Rho
    C1=ifo.Suspension[WireMat2].C
    K1=ifo.Suspension[WireMat2].K
    alpha1=ifo.Suspension[WireMat2].Alpha
    beta1=ifo.Suspension[WireMat2].dlnEdT
    phi1=ifo.Suspension[WireMat2].Phi
    E1=ifo.Suspension[WireMat2].Y

    T1 = tension_per_fibre(m1+m,N1,g)
    S1=cross_section_area(dd1/2)

# Thermoelastic dissipation contribution
    phith1=np.zeros([len(w)],dtype = "complex")
    Delta1 = dissipation_delta(E1, alpha1, beta1, m1+m, N1, dd1/2, Temp, rho1, C1)
    tau1=0.0737*rho1*C1*(dd1)**2/K1
    for i in range(len(w)):
        phith1[i]=Delta1*w[i]*tau1/(1+w[i]**2*tau1**2)

    phiT1=phith1+phi1

    M1=m+m1
    I1=pi/4*(dd1/2)**4


# Dilution factors for marionette and RM stages
    def dillution_factor(YM, radius, N, mass, sus_length):
        I = cross_section_momentum(radius)
        dilution_factor = np.sqrt(YM*I*N/mass/g/(sus_length**2))
        return dilution_factor

    dil1 = dillution_factor(E1,dd1/2,N1,M1,LL1)

# Gravitational reaction constants
    def gravitational_reaction_factor(YM, radius, N, mass, sus_length):
        grf = mass*g/sus_length*(1+dillution_factor(YM, radius, N, mass, sus_length))
        return grf


    k0p1 = gravitational_reaction_factor(E1, dd1/2, N1, M1, LL1)

    w01=sqrt(k0p1/M1)

# Horizontal viscous Q's
    Q=ifo.Suspension.Stage[0].Qvh
    Q1=ifo.Suspension.Stage[2].Qvh


# Total loss angle internal+viscous
    phi1tot=dil1*phiT1+w/(Q1*w01)
    kp1=k0p1*(1+1j*phi1tot)

# Reaction constants for mirror stage fused silica fibres VIR-0.15A-09
# sec. 2
    Kpend=np.zeros([len(w)],dtype = "complex")

    for i in range(len(w)):
        if Optcyl==0: #Cylindrical fibres
            E=E0*(1+1j*phiT[i])
            Lambda=sqrt(1/(2*E*II)*(T+T*sqrt(1+4*E*II*w[i]**2*rho*S/T**2)))
            p=sqrt(1/(2*E*II)*(-T+T*sqrt(1+4*E*II*w[i]**2*rho*S/T**2)))
            Kpend[i]=4*E*II*Lambda*p*(Lambda**3*cos(p*L)+Lambda**2*sin(p*L)*p+p**3*sin(p*L)+Lambda*cos(p*L)*p**2)/(-2*Lambda*cos(p*L)*p+Lambda**2*sin(p*L)-sin(p*L)*p**2)
        else: # Optimized fibres
            [ZZ,XX,Fend,ENER]=Optimisedfibres(f[i],phiT[i],ifo)
            z[i,:]=np.concatenate(ZZ)
            x[i,:]=np.concatenate(XX)
            #[z[i,:], x[i,:], Fend]=An_suspn_3(f[i],phiT[i],ifo)Xh
            Kpend[i]=4*Fend/delta

    w0=sqrt(real(Kpend[0])/m);
    Kpend=real(Kpend)+1j*(imag(Kpend)+w/(w0*Q)*abs(real(Kpend))) #viscous dissipation is added

# Suspension thermal noise is calculated using fluctuation-dissipation theorem. The equations of motion are written
# in a martix formalism Ah*Xh=B. The admittance is calculated. VIR-0.15A-09
# sec. 1.1-1.2-1.3

    Ah = np.zeros([len(f),2,2],dtype =  "complex")
    Ah[:,0,0]=kp1+Kpend-m1*w**2
    Ah[:,0,1]=-Kpend
    Ah[:,1,0]=-Kpend
    Ah[:,1,1]=Kpend-m*w**2
    B=np.array([0,1])

    Xh=np.zeros([2,len(w)],dtype =  "complex")
    def solver_AXB(Ah_matrix,B_matrix):
        for i in range(len(w)):
            Xh[:,i]=np.linalg.solve(Ah_matrix[i,:,:],np.transpose(B_matrix))

    Xh[:,i] = solver_AXB(Ah,B)

    #admitt=1j*w*Xh[1,:]
    def admittance(impedance):
        admittance = 1j*w*impedance
        return admittance
    def FDT(Temperature, admittance, w):
        s_density=4*kb*Temperature*np.real(admittance)/(w**2)*4
        return s_density

    admittance_h = admittance(Xh[1,:])
    # Horizontal Suspension Thermal Noise PDS
    horizontal_PSD = FDT(Temp, admittance_h,w)




 #***********Vertical motion*********

# Surface contribution for circular cross section fibres (vertical motion)
    mu_v=2/r
    phiv=phi0*(1 + mu_v*ds)

    wv1=2*pi*0.4 # 0.4Hz is the resonance frequency of the marionette stage due to the anti-spring system
#Vertical resonance frequency for mirror wires (from Hook's law)
    wv=sqrt(4*E0*S0/L/m)

#Vertical viscous Q's
    Qv=ifo.Suspension.Stage[0].Qvv
    Qv1=ifo.Suspension.Stage[2].Qvv


# Total loss angle internal+viscous
    phiv=phiv+w/(Qv*wv)
    phiv1=phi1+w/(Qv1*wv1)


# Elastic constants for vertical motion
    kv1=wv1**2*M1*(1+1j*phiv1)
    kv=wv**2*m*(1+1j*phiv)

# Suspension thermal noise is calculated using fluctuation-dissipation theorem. The equations of motion are written
# in a martix formalism Ah*Xh=B. The admittance is calculated. VIR-0.15A-09
# sec. 1.1-1.2-1.3
    Xv = np.zeros([2,len(f)],dtype = "complex")
    Av = np.zeros([len(f),2,2],dtype = "complex")
    Av[:,0,0]=kv1+kv-m1*w**2
    Av[:,0,1]=-kv
    Av[:,1,0]=-kv
    Av[:,1,1]=kv-m*w**2
    B=[0,1]

    for i in range(len(w)):
        Xv[:,i]=np.linalg.solve(Av[i,:,:],np.transpose(B))


    admitt_v = admittance(Xv[1,:])
    Vertical_PSD = FDT(Temp, admitt_v,w)                         # Vertical suspension thermal noise PSD

    VerticaltoH_PSD = thetaVH**2*Vertical_PSD;                       # Vertical coupled to horizontal noise
    noise=horizontal_PSD + VerticaltoH_PSD;                          # Total th noise PSD
    return noise, horizontal_PSD, VerticaltoH_PSD




def Optimisedfibres(f,phi,ifo):
    WireMat0=ifo.Suspension.Stage[0].WireMaterial
    WireMat1=ifo.Suspension.Stage[1].WireMaterial
    WireMat2=ifo.Suspension.Stage[2].WireMaterial


    #ifo=ifopar;

    N=ifo.Suspension.Stage[0].NWires         # Number of suspension's fibres

    E0=ifo.Suspension[WireMat0].Y
    E=E0*(1+1j*phi);                         # Fused Silica Young Modulus {Pa]
    rho=ifo.Suspension[WireMat0].Rho            # FS density [kg/m^3]

    g=const.g                                # Gravitational field [N/kg]
    m=ifo.Suspension.Stage[0].Mass           # Mirror mass [kg]


    L1=ifo.Suspension.Stage[0].Length1
    L2=ifo.Suspension.Stage[0].Length2
    L3=ifo.Suspension.Stage[0].Length3

    d1=ifo.Suspension.Stage[0].WireDiameter1
    d2=ifo.Suspension.Stage[0].WireDiameter2
    d3=ifo.Suspension.Stage[0].WireDiameter3

    r1=d1/2
    r2=d2/2
    r3=d3/2

    S1=cross_section_area(r1)
    I1=cross_section_momentum(r1)
    S2=cross_section_area(r2)
    I2=cross_section_momentum(r2)
    S3=cross_section_area(r3)
    I3=cross_section_momentum(r3)

    T=tension_per_fibre(m,N,g)        # Fibre tension
    delta=1e-3                        # shift at the end of the fibre

    w=2*pi*f

# The solutions for the 3 segments contain 4 integration constants each (see VIR-015A-09
# eq(50)). *coeff* is a vector that contains the twelve integration constants.
# Twelve equation are needed. They are the 4 boundary conditions and the 8 joining conditions (x, x',force, moment)
# in the diameter discontinuity points.

    lambda1=sqrt((T+sqrt(T**2+4*E*I1*w**2*rho*S1))/(2*E*I1))
    p1=sqrt((-T+sqrt(T**2+4*E*I1*w**2*rho*S1))/(2*E*I1))

    lambda2=sqrt((T+sqrt(T**2+4*E*I2*w**2*rho*S2))/(2*E*I2));
    p2=sqrt((-T+sqrt(T**2+4*E*I2*w**2*rho*S2))/(2*E*I2))

    lambda3=sqrt((T+sqrt(T**2+4*E*I3*w**2*rho*S3))/(2*E*I3));
    p3=sqrt((-T+sqrt(T**2+4*E*I3*w**2*rho*S3))/(2*E*I3))

    AA=np.zeros([12,12],dtype = "complex")

    ex1=exp(-lambda1*L1)
    ex3=exp(-lambda3*L3)
    sn1=sin(p1*L1)
    cn1=cos(p1*L1)
    sn2=sin(p2*L2)
    cn2=cos(p2*L2)
    sn3=sin(p3*L3)
    cn3=cos(p3*L3)

# eq1: x1(0)=0
    AA[0,:]=[1,ex1,1,0,0,0,0,0,0,0,0,0]

# eq2: x1'(0)=0
    AA[1,:]=[-lambda1,lambda1*ex1,0,p1,0,0,0,0,0,0,0,0]

# eq3: x1(L1)=x2(0)
    AA[2,:]=[ex1,1,cn1,sn1,-1,0,-1,0,0,0,0,0]

# eq4: x1'(L1)=x2'(0)
    AA[3,:]=[-lambda1*ex1,lambda1,-p1*sn1,p1*cn1,lambda2,0,0,-p2,0,0,0,0]

# eq5: F1(L1)=F2(0)
    AA[4,:]=[I1*(-lambda1**3*ex1),I1*lambda1**3,I1*p1**3*sn1,I1*(-p1**3*cn1),I2*lambda2**3,0,0,I2*p2**3,0,0,0,0]

# eq6: M1(L1)=M2(0)
    AA[5,:]=[I1*lambda1**2*ex1,I1*lambda1**2,I1*(-p1**2*cn1), I1*(-p1**2*sn1),-I2*lambda2**2,0,I2*p2**2,0,0,0,0,0]

# eq7: M2(L2)=M3(0)
    AA[6,:]=[0,0,0,0,0,I2*lambda2**2,-I2*p2**2*cn2,-I2*p2**2*sn2,-I3*lambda3**2,-I3*lambda3**2*ex3,I3*p3**2,0 ]

# eq8: F2(L2)=F3(0)
    AA[7,:]=[0,0,0,0,0,I2*lambda2**3,I2*p2**3*sn2,I2*(-p2**3*cn2),I3*lambda3**3,-I3*lambda3**3*ex3,0,I3*p3**3]

# eq9: x2'(L2)=x3'(0)
    AA[8,:]=[0,0,0,0,0,lambda2,-p2*sn2,p2*cn2,lambda3,-lambda3*ex3,0,-p3]

# eq10: x2(L2)=x3(0)
    AA[9,:]=[0,0,0,0,0,1,cn2,sn2,-1,-ex3,-1,0]

# eq11: x3'(L3)=0
    AA[10,:]=[0,0,0,0,0,0,0,0,-lambda3*ex3,lambda3,-p3*sn3,p3*cn3]

# eq12: x3(L3)'=delta
    AA[11,:]=[0,0,0,0,0,0,0,0,ex3,1,cn3,sn3]

#Force
    BB=[0,0,0,0,0,0,0,0,0,0,0,delta]

    coef=np.linalg.solve(AA,BB)

# first fibre's segment
    z1=np.linspace(0,L1,1000)

    x1 = coef[0]*exp(-lambda1*z1)+coef[1]*exp(-lambda1*(L1-z1))+coef[2]*cos(p1*z1)+coef[3]*sin(p1*z1)
    x12 = coef[0]*lambda1**2*exp(-lambda1*z1)+coef[1]*lambda1**2*exp(-lambda1*(L1-z1))+coef[2]*(-p1**2)*cos(p1*z1)+coef[3]*(-p1**2)*sin(p1*z1)

    Energy1=0.5*E*I1*x12**2*z1[1];   #bending energy E=1/2 E I x''^2

# second fibre's segment
    z2=np.linspace(0,L2,1000)
    x2 = coef[4]*exp(-lambda2*z2)+coef[5]*exp(-lambda2*(L2-z2))+coef[6]*cos(p2*z2)+coef[7]*sin(p2*z2)
    x22 = coef[4]*lambda2**2*exp(-lambda2*z2)+coef[5]*lambda2**2*exp(-lambda2*(L2-z2))+coef[6]*(-p2**2)*cos(p2*z2)+coef[7]*(-p2**2)*sin(p2*z2)
    Energy2=0.5*E*I2*x22**2*z2[1]    #bending energy E=1/2 E I x''^2


# third fibre's segment
    z3=np.linspace(0,L3,1000)
    x3 = coef[8]*exp(-lambda3*z3)+coef[9]*exp(-lambda3*(L3-z3))+coef[10]*cos(p3*z3)+coef[11]*sin(p3*z3)
    x32 = coef[8]*lambda3**2*exp(-lambda3*z3)+coef[9]*lambda3**2*exp(-lambda3*(L3-z3))+coef[10]*(-p3**2)*cos(p3*z3)+coef[11]*(-p3**2)*sin(p3*z3)

    x33 = coef[8]*(-lambda3**3)*ex3+coef[9]*lambda3**3+coef[10]*p3**3*sn3+coef[11]*(-p3**3)*cn3
    Energy3=0.5*E*I3*x32**2*z3[1]    #bending energy E=1/2 E I x''^2

    z=[z1,z2+z1[-1],z3+z2[-1]+z1[-1]]
    x=[x1,x2,x3]

    Fend=-E*I3*x33                   #Force at the end of the fibre
    Energy=[Energy1,Energy2,Energy3] #Energy distribution along the fibre
    return z,x,Fend,Energy
