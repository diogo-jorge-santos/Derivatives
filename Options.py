import numpy as np
from scipy.stats import norm

def binomial_pricing(s0,T,sigma,r,c,n,k,european,call):
    """
    s0=current price
    T=number of years
    sigma=anualized volatility
    r=risk free rate, compounded continuously
    c=dividend rate
    n=number of steps
    k=strike price
    european=True if european option , False if american option
    call=True if call option, False if put option 
    u = factor change of upstate
    d = factor change of downstate
    """
    u=np.exp(sigma*((T/n)**0.5))
    d=1/u
    
    #risk neutral
    q=(np.exp((r-c)*T/n)-d)/(u-d)
    
    #prices lattice
    lattice=np.zeros((n+1,n+1))
    
    for i in range(n + 1):
        for j in range(i + 1):
            lattice[j, i] = s0*(u**(i - j))*(d**j)
    option=np.zeros((n+1,n+1))
    #option pricing
    if call==True:
        
        for i in range(n+1):
            option[i,n]=max(lattice[i,n]-k,0)
    
        if european==True:
            for i in range(n-1,-1,-1):
                for j in range(0,i+1):
                    option[j,i]=np.exp(-T*r/n)*(q*option[j,i+1]+(1-q)*option[j+1,i+1])
                    
        #allowing exercising if it is american option
        if european==False:
            for i in range(n-1,-1,-1):
                for j in range(0,i+1):
                    option[j,i]=max((np.exp(-T*r/n)*(q*option[j,i+1]+(1-q)*option[j+1,i+1])),lattice[j,i]-k)
    if call==False:
        
        for i in range(n+1):
            option[i,n]=max(k-lattice[i,n],0)
            
        if european==True:
            for i in range(n-1,-1,-1):
                for j in range(0,i+1):
                    option[j,i]=np.exp(-T*r/n)*(q*option[j,i+1]+(1-q)*option[j+1,i+1])
                    
        #allowing exercising if it is american option            
        if european==False:
            for i in range(t-1,-1,-1):
                for j in range(0,i+1):
                    option[j,i]=max((np.exp(-T*r/n)*(q*option[j,i+1]+(1-q)*option[j+1,i+1])),max((k-lattice[j,i]),0))
    return (option[0,0])

def bs_pricing(s0,T,sigma,r,c,k,call=True):
    d1=(np.log(s0/k)+(r-c+sigma*sigma/2)*T)/(sigma*T**0.5)
    d2=d1-sigma*T**0.5
    c0=s0*np.exp(-c*T)*norm.cdf(d1)-k*np.exp(-r*T)*norm.cdf(d2)
    if call==True:
        return c0
    else:
        return  c0+k*np.exp(-r*T)-s0*np.exp(-c*T)

def greeks(s0,T,sigma,r,c,n,k,european,call,delta_finite):
    """
    for european optios it will use analytic greeks
    for american option it will use finite difference
    will assume that t=0
    """
    if european==True:
        d1=(np.log(s0/k)+(r-c+sigma*sigma/2)*T)/(sigma*T**0.5)
        d2=d1-sigma*T**0.5
        if call==True:
            delta=norm.cdf(d1)
            theta=-(s0*norm.pdf(d1)*sigma)/(2*T**0.5)-r*k*np.exp(-r*T)*norm.cdf(d2)
            rho=k*T*np.exp(-r*T)*norm.cdf(d2)
        if call==False:
            delta=norm.cdf(d1)-1
            theta=-(s0*norm.pdf(d1)*sigma)/(2*T**0.5)+r*k*np.exp(-r*T)*norm.cdf(-d2)
            rho=-k*T*np.exp(-r*T)*norm.cdf(-d2)
        gamma=norm.pdf(d1)/(s0*sigma*T**0.5)
        vega=s0*norm.pdf(d1)*T**0.5
    if european==False:
        c0=binomial_pricing(s0,T,sigma,r,c,n,k,european,call)
        delta=(binomial_pricing(s0+delta_finite,T,sigma,r,c,n,k,european,call)-c0)/delta_finite
        vega=(binomial_pricing(s0,T,sigma+delta_finite,r,c,n,k,european,call)-c0)/delta_finite
        theta=-((binomial_pricing(s0,T+delta_finite,sigma,r,c,n,k,european,call)-c0)/delta_finite)
        rho=(binomial_pricing(s0,T,sigma,r+delta_finite,c,n,k,european,call)-c0)/delta_finite
        gamma=(binomial_pricing(s0+delta_finite,T,sigma,r,c,n,k,european,call)-2*c0+binomial_pricing(s0-delta_finite,T,sigma,r,c,n,k,european,call))/(delta_finite**2)
    print(f"delta: {delta}")
    print(f"gamma: {gamma}") 
    print(f"vega: {vega}") 
    print(f"rho: {rho}") 
    print(f"theta: {theta}") 
        
