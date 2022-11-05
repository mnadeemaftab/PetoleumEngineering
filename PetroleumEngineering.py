# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 23:57:48 2022

@author: Muhammad Nadeem Afta
"""

# functions to calculate gas deviation (Z) and gas formation volume factors

# first fuction to correct gas gravity for impurities
import math

def SG_HC(SG_gas=0.75, yN2=0, yCO2=0.5, yH2S=1):
    """Calculates the specific gravity of hydrocorbon with measured specific gravity of gas
    and composition of impurities with Standing-Katz Gas (1977). SG_gas measured sepcific gravity of gas,
    yN2, yCO2, yH2S compositions of N2, CO2 and H2S in percentage"""
    
    #convert composition from percentage to fraction
    yN2 = yN2/100
    yCO2 = yCO2/100
    yH2S = yH2S/100
    
    
    return (SG_gas-0.9672*yN2-1.5195*yCO2-
            1.1765*yH2S)/(1-yN2-yCO2-yH2S)

#function to calcuate compressiblity factor Z

def Z_DPR(SG_gas=0.75, yN2=0, yCO2=0.5, yH2S=1, pres=500, temp = 150):
     
    """ calculates the Z deviation factor with Dranchuk Purvis Robinson (1974) correlations
    and CO2 and H2S correction with Wichert-Aziz (1972) Correlations. SG_gas measured sepcific gravity of gas,
    yN2, yCO2, yH2S compositions of N2, CO2 and H2S in percentage. Pressure and temperature are in psia and deg F"""
    
    # calculate specific gravity of hydrocorbon by calling function SG_HC
    
    SG_HC1 = SG_HC(SG_gas, yN2, yCO2, yH2S)
    #print(SG_HC1)
    
    #convert composition from percentage to fraction
    yN2 = yN2/100
    yCO2 = yCO2/100
    yH2S = yH2S/100
    # convert temperature from deg F to deg R
    temp = temp + 459.67
    
    #Calculate Pseudo-critical temperature and pressure for the hydrocarbon fraction
    
    TpcHC = 168+325*SG_HC1-12.5*math.pow(SG_HC1,2)
    PpcHC = 677+15.0*SG_HC1-37.5*math.pow(SG_HC1,2)
    #print(PpcHC)
    
    #pseudo-critical temperature and pressure for the whole gas
    
    TpcM = (1-yN2-yCO2-yH2S)*TpcHC+227.3*yN2+547.6*yCO2+672.4*yH2S
    PpcM = (1-yN2-yCO2-yH2S)*PpcHC+493.0*yN2+1071*yCO2+1306*yH2S
    #print(PpcM)
    #Wichert-Aziz correction for CO2 and H2S 
    
    eita = 120*((yCO2+yH2S)**0.9 - (yCO2+yH2S)**1.6)+15*(yH2S**0.5-yH2S**4)
    TpcM_st = TpcM-eita
    PpcM_st = PpcM*(TpcM-eita)/(TpcM+yH2S*(1-yH2S)*eita)
    #print(eita)
    #Calculate the Tpr and Ppr 
    Tpr = temp/TpcM_st
    Ppr = pres/PpcM_st
    #print(Tpr)
    #print(PpcM_st)
    #print(Ppr)
    A1, A2, A3, A4, A5, A6, A7, A8 = 0.31506237, -1.0467099, -0.57832729, 0.53530771, -0.61232032, -0.10488813, 0.68157001, 0.68446549
    
    A = A5*A6
    B = A4*Tpr+A5
    C = A1*Tpr+A2+A3/Tpr**2
    D = Tpr
    E = A7/Tpr**2
    F = A8
    G = 0.27*Ppr
    #intial reduced density
    den_red = 0.27 * Ppr/Tpr 
    # function with above reduced density
    Fp = A*den_red**6+B*den_red**3+C*den_red**2+D*den_red+E*den_red**3*(1+F*den_red**2\
        )*math.exp(-F*den_red**2)-G
    # Derivative of above function    
    Fp_der = 6*A*den_red**5+3*B*den_red**2+2*C*den_red+D+E*den_red**2*(3+F*den_red**2*\
             (3-2*F*den_red**2))*math.exp(-F*den_red**2)
    #iteration until funtion Fp becomes zero
    while Fp != 0:
        den_red = den_red - Fp/Fp_der
        
        Fp = round(A*den_red**6+B*den_red**3+C*den_red**2+D*den_red+E*den_red**3*(1+F*den_red**2\
        )*math.exp(-F*den_red**2)-G,6)
        Fp_der = 6*A*den_red**5+3*B*den_red**2+2*C*den_red+D+E*den_red**2*(3+F*den_red**2*\
             (3-2*F*den_red**2))*math.exp(-F*den_red**2)
        
        #print(den_red)
        #print(Fp)
        #print(Fp_der)
        
    return 0.27*Ppr/Tpr/den_red

    


# funtion to calcuate the Gas formation volume factor

def Bg(SG_gas=0.75, yN2=0, yCO2=0.5, yH2S=1, pres=500, temp = 150):
    """ Calcuate the Gas formation volume factor in the units of cu.ft/SCF. SG_gas measured sepcific gravity of gas,
    yN2, yCO2, yH2S compositions of N2, CO2 and H2S in percentage. Pressure and temperature are in psia and deg F"""
           
    # call gas deviation factor to calculate Z factor.
    Z_DPR1 = Z_DPR(SG_gas, yN2, yCO2, yH2S, pres, temp)
    # convert temperature from deg F to deg R
    temp = temp + 459.67
       
    return 0.0283*Z_DPR1* temp/pres
       

      
# funtion to calculate oil formation volume factor

def Bo(oil_API=30, SG_gas=0.75, pres=500, temp=150):
    
    """ calculate the oil formation volume factor with Katz (SG_gst, 1942) and Standing (Rst, 1977) and Standing (Bo, 1977) correlations,
    oil_API, SG_gas, pres and temp are oil API at 60F, Gas specific gravity, pressure in psia and temperature in deg F"""
    #convert oil API to oil SG
    oil_SG = 141.5/(oil_API + 131.5)
    #estimate tank GOR
    GOR_tank = pres*(1.797/oil_SG - 1.838)
    #calcualte SG of stock tank gas
    SG_gst = oil_API*(0.02-3.57*10**-6*GOR_tank)+0.25
    # calculate Rs
    Rs = SG_gst*((pres/18.2+1.4)*10**(0.0125*oil_API-0.00091*temp))**1.20482
    
    #iteration process to find SG_gst and Rs.
    while (Rs - GOR_tank) != 0:
        GOR_tank = Rs
        SG_gst = oil_API*(0.02-3.57*10**-6*GOR_tank)+0.25
        Rs = SG_gst*((pres/18.2+1.4)*10**(0.0125*oil_API-0.00091*temp))**1.20482
        #print(Rs)
     
    Bo = 0.9759+12*10**-5*(Rs*math.sqrt(SG_gst/oil_SG)+1.25*temp)**1.2
    #print(oil_SG)
    
    return Bo



def Rs(oil_API=30, SG_gas=0.75, pres=500, temp=150):
    
    """ calculate the dissolved gas in oil at given pressure and temperature with Katz (SG_gst, 1942) and Standing (Rst, 1977) and Standing (Bo, 1977) correlations,
    oil_API, SG_gas, pres and temp are oil API at 60F, Gas specific gravity, pressure in psia and temperature in deg F"""
    #convert oil API to oil SG
    oil_SG = 141.5/(oil_API + 131.5)
    #estimate tank GOR
    GOR_tank = pres*(1.797/oil_SG - 1.838)
    #calcualte SG of stock tank gas
    SG_gst = oil_API*(0.02-3.57*10**-6*GOR_tank)+0.25
    # calculate Rs
    Rs = SG_gst*((pres/18.2+1.4)*10**(0.0125*oil_API-0.00091*temp))**1.20482
    
    #iteration process to find SG_gst and Rs.
    while (Rs - GOR_tank) != 0:
        GOR_tank = Rs
        SG_gst = oil_API*(0.02-3.57*10**-6*GOR_tank)+0.25
        Rs = SG_gst*((pres/18.2+1.4)*10**(0.0125*oil_API-0.00091*temp))**1.20482
        #print(Rs)
     
    #Bo = 0.9759+12*10**-5*(Rs*math.sqrt(SG_gst/oil_SG)+1.25*temp)**1.2 (not required for Rs)
    #print(oil_SG)
    
    return Rs



def Pb(oil_API=30, SG_gas=0.75, Rs=500, temp=150):
    """calculate bubble point pressure with Standing (1977) correlations"""
    #convert oil API to oil SG
    #oil_SG = 141.5/(oil_API + 131.5) not required here
    #correlation constants
    a1, a2, a3,  a4, a5 = 18.2, 0.83, 0.00091, 0.0125, 1.4
    #correlation
    X = a3*temp-a4*oil_API
    Pb = a1*((Rs/SG_gas)**a2*10**X - a5)
      
    return Pb

def ql_gilbert(choke=64, pwh=500, Rs=500):
    """calculate liquid rate by Gilbert Equation"""
      
    Gilbert = [3.86E-03, 0.546, 1.89]
    
    n_choke = choke/64.0
    n_pwh = pwh
    n_Rs = Rs
    
    ql_gilbert = (n_pwh * n_choke**Gilbert[2]) / (Gilbert[0] * n_Rs**Gilbert[1])
    
    return ql_gilbert



def ql_rose(choke=64, pwh=500, Rs=500):
    """calculate liquid rate by Rose Equation"""
    n_choke = choke/64.0
    n_pwh = pwh
    n_Rs = Rs

    Rose = [4.26E-03, 0.500, 2.00]
    ql_rose = (n_pwh * n_choke**Rose[2]) / (Rose[0] * n_Rs**Rose[1])
    
    return ql_rose



def ql_baxendell(choke=64, pwh=500, Rs=500):
    """calculate liquid rate by Baxendlell Equation"""
    n_choke = choke/64.0
    n_pwh = pwh
    n_Rs = Rs

    Baxendell = [3.12E-03, 0.546, 1.93]
    ql_baxendell = (n_pwh * n_choke**Baxendell[2]) / (Baxendell[0] * n_Rs**Baxendell[1])
    
    return ql_baxendell



def ql_achong(choke=64, pwh=500, Rs=500):
    """calculate liquid rate by Achong Equation"""
    n_choke = choke/64.0
    n_pwh = pwh
    n_Rs = Rs
    
    Achong = [1.54E-03, 0.650, 1.88]
    ql_achong = (n_pwh * n_choke**Achong[2]) / (Achong[0] * n_Rs**Achong[1])
    
    return ql_achong



def ql_aftab(choke=64, pwh=500, Rs=500):
    """calculate liquid rate by Aftab Equation"""
    n_choke = choke/64.0
    n_pwh = pwh
    n_Rs = Rs
    
    Aftab = [2.97E-03, 0.5605, 1.925]
    ql_aftab = (n_pwh * n_choke**Aftab[2]) / (Aftab[0] * n_Rs**Aftab[1])
    
    return ql_aftab


def ql_kargarpour(choke=64, pwh=500, GLR=500, pd=250, SG=0.85):
    """Calcuate liquid rate by Mohammad Ali Kargarpour equation based on critical and sub critical flow
    where SG is SG of liquid and pd is downstream pressure"""
    
    d = choke/64.0
    P1 = pwh
    P2 = pd
    PR = P2/P1
    #condition if P2 is equal or greater than P1 gives exception but greater than to zero
    if PR < 1.0 and P2 > 0:
        Eq1 = math.sqrt(P1)
        Eq2 = 552*math.sqrt((1-PR)/SG)
        Eq3 = GLR/(65554*math.sqrt((PR)**1.5625*(1-(PR)**0.21875)))
        Eq4 = GLR/14387
        ql_subcritical = P1*d**2*(Eq1/Eq2+Eq3)**(-1)
        ql_critical = P1*d**2*(Eq1/Eq2 + Eq4)**(-1)
        #print(ql_subcritical, ql_critical)
    
        #divide flow based on critical and non critical conditions (i.e. P2/P1)            
        if PR <= 0.55:
            ql_kargarpour = ql_critical
        else:
            ql_kargarpour = ql_subcritical
    # condition if P2 is qual to zero or negative
    elif P2 == 0 or P2 <0:
        Eq1 = math.sqrt(P1)
        Eq2 = 552*math.sqrt((1-PR)/SG)
        #Eq3 = GLR/(65554*math.sqrt((PR)**1.5625*(1-(PR)**0.21875)))
        Eq4 = GLR/14387
        #ql_subcritical = P1*d**2*(Eq1/Eq2+Eq3)**(-1)
        ql_critical = P1*d**2*(Eq1/Eq2 + Eq4)**(-1)
        ql_kargarpour = ql_critical
    
    # exception error incase pd is equal or higher than pwh
    else:
        raise Exception("pd cannot be equal or higher than pwh")
             
    return ql_kargarpour


# just change choke, pwh, GLR, pd and SG in in the below list to calculate rates.
#x = [choke, pwh, GLR, pd, SG]

x = [32,350,365,100,0.85]

print("Flow rate by Aftab:")
print(ql_aftab(x[0],x[1],x[2]))

print()
print("Average flow rate by Gilbert, Rose, Baxendell and Achong:")
print((ql_gilbert(x[0],x[1],x[2]) + ql_rose(x[0],x[1],x[2]) + ql_baxendell(x[0],x[1],x[2]) + ql_achong(x[0],x[1],x[2]))/4.0)

print()

print("Flow rate by Kargarpour considering downstream pressure:")
print(ql_kargarpour(x[0],x[1],x[2],x[3],x[4]))

print()

#calculate shrinkage factor, GOR2, Pb, Z facor and expansion factor
#x = [SG_gas, yN2, yCO2, yH2S, pres, temp]
x = [0.85, 0, 1.2, 0.5, 600, 130] # for Bg and Z factor
# for Rs and Bo
#y = [oil_API, SG_gas, pres, temp]
y = [35, 0.80, 600, 130]
# for Pb
#z = [oil_API, SG_gas, Rs, temp]
z = [35, 0.80, 300, 180]
print(1/Bo(y[0],y[1],y[2],y[3]), Rs(y[0],y[1],y[2],y[3]), Pb(z[0],z[1],z[2],z[3]),\
      Z_DPR(x[0],x[1],x[2],x[3],x[4],x[5]), 1/Bg(x[0],x[1],x[2],x[3],x[4],x[5]))

