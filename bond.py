import numpy as np
import math
import scipy.optimize
"""
    spot=vector of spot rates
    b=vector of sd of spot rates
    face=face value of the bond
    c=annual cupon rate

    
    """
a=np.array([3.00,3.04,3.07,3.09,3.10,3.10,2.84,2.77,2.70,2.63])
b=np.array([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1])
spot=np.array([3,3.1,3.2,3.3,3.4,3.5,3.55,3.6,3.65,3.7])
def bdt_aux(a,b,spot):
    
    bdt=np.zeros((spot.size,spot.size))
    bdt[0][0]=a[0]
    for i in range(1,spot.size):
        for j in range(0,i+1):
            bdt[i][j]=a[i]*math.exp(b[i]*(i-j))
    q=0.5
    #fazer com os elementary prices
    element=np.zeros((spot.size+1,spot.size+1))
    element[0][0]=1
    for i in range(1,spot.size+1):
        for j in range(0,i+1):
            if(j==i):
                element[i][j]=q*element[i-1][j-1]/(1+bdt[i-1][j-1]/100)
            elif(j==0):
                element[i][j]=(1-q)*element[i-1][j]/(1+bdt[i-1][j]/100)
            else:
                element[i][j]=q*element[i-1][j-1]/(1+bdt[i-1][j-1]/100)+(1-q)*element[i-1][j]/(1+bdt[i-1][j]/100)
    bdtprice=np.zeros(spot.size)
    for i in range(1,spot.size+1):
        aux=0
        for j in range(0,i+1):
            aux=aux+element[i][j]
        bdtprice[i-1]=aux

    bdtspot=np.zeros(spot.size)
    func=0
    for i in range(0,spot.size):
        bdtspot[i]=100*((1/bdtprice[i])**(1/(i+1))-1)
        func=func+(bdtspot[i]-spot[i])**2
    return func
def interest_tree(spot,b):
    aux=np.zeros((spot.size))
    auxx=scipy.optimize.minimize(bdt_aux,aux,args=(b,spot))
    a=auxx.x
    bdt=np.zeros((spot.size,spot.size))
    bdt[0][0]=a[0]
    for i in range(1,spot.size):
        for j in range(0,i+1):
            bdt[i][j]=a[i]*math.exp(b[i]*(i-j))
    return bdt
def price_tree(face,spot,b,c):
    
    bdt=interest_tree(spot,b)
    ep=np.zeros((spot.size+1,spot.size+1))
    for i in range(0,spot.size+1):
        ep[spot.size][i]=face*(1+c)
    for i in range(spot.size-1,-1,-1):
        for j in range(0,i+1):
            ep[i][j]=100*c+(0.5*ep[i+1][j]+(0.5)*ep[i+1][j+1])/(1+(bdt[i,j]/100)) 
    return ep
def cap(spot,b,k,exp):
    t=exp-1
    tree=interest_tree(spot,b)
    cap=np.zeros((t+1,t+1))
    for i in range(0,t+1):
        cap[t][i]=max(0,(tree[t][i]/100-k/100)/(1+tree[t][i]/100))

    for i in range(t-1,-1,-1):
        for j in range(0,i+1):
            cap[i][j]=(0.5*cap[i+1,j]+0.5*cap[i+1,j+1])/(1+tree[i][j]/100)
    print(cap)
    return(cap[0][0])

def floor(spot,b,k,exp):
    t=exp-1
    tree=interest_tree(spot,b)
    floor=np.zeros((t+1,t+1))
    for i in range(0,t+1):
        floor[t][i]=max(0,(k/100-tree[t][i]/100)/(1+tree[t][i]/100))

    for i in range(t-1,-1,-1):
        for j in range(0,i+1):
            floor[i][j]=(0.5*floor[i+1,j]+0.5*floor[i+1,j+1])/(1+tree[i][j]/100)
    return(floor[0][0])

def swap(spot,b,fixed,exp):
    t=exp-1
    tree=interest_tree(spot,b)
    swap=np.zeros((t+1,t+1))
    for i in range(0,t+1):
        swap[t][i]=max(0,(tree[t][i]/100-fixed/100)/(1+tree[t][i]/100))

    for i in range(t-1,-1,-1):
        for j in range(0,i+1):
            swap[i][j]=((tree[i][j]/100-fixed/100)+0.5*swap[i+1,j]+0.5*swap[i+1,j+1])/(1+tree[i][j]/100)
    return(swap[0][0])



def bond_opt(spot,b,face,k,t,c,european=True,call=True):
    tree=price_tree(face,spot,b,c)
    interest=interest_tree(spot,b)
    option=np.zeros((t+1,t+1))
    if european:
        if call:
            for i in range(0,t+1):
                option[t][i]=max(0,(tree[t][i]-k))

            for i in range(t-1,-1,-1):
                for j in range(0,i+1):
                    option[i][j]=(0.5*option[i+1][j]+0.5*option[i+1][j+1])/(1+interest[i][j]/100)
        else:
            for i in range(0,t+1):
                option[t][i]=max(0,(k-tree[t][i]))

            for i in range(t-1,-1,-1):
                for j in range(0,i+1):
                    option[i][j]=(0.5*option[i+1][j]+0.5*option[i+1][j+1])/(1+interest[i][j]/100)
        
    else:
        if call:
            for i in range(0,t+1):
                option[t][i]=max(0,(tree[t][i]-k))

            for i in range(t-1,-1,-1):
                for j in range(0,i+1):
                    option[i][j]=max(tree[i,j]-k,(0.5*option[i+1][j]+0.5*option[i+1][j+1])/(1+interest[i][j]/100))
        else:
            for i in range(0,t+1):
                option[t][i]=max(0,(tree[t][i]-k))

            for i in range(t-1,-1,-1):
                for j in range(0,i+1):
                    option[i][j]=max(k-tree[i,j],(0.5*option[i+1][j]+0.5*option[i+1][j+1])/(1+interest[i][j]/100))     
            
   
        
    return(option[0][0])
