def bedrockchannel(tend,uplift,kappa1,kappa2,deltaz):

    import numpy as np

    # Initial Topography
    nx=200
    dx=10
    xgrid=np.arange(0,nx*dx,dx)            # Grid
    area=np.zeros(nx)
    area=500+0.5*xgrid**2                  # Hack's law relating drainage area and stream length

    # Channel Width
    # width=8*(area*1e-6).^(0.5)

    # Define Slope of Steady State Channel Profile Based on Stream Power Erosion Law
    S=np.zeros(nx);
    S=(uplift/kappa1)*area**(-1/2)

    # Define Steady State Topography Based on Slope
    topo=np.zeros(nx)
    for i in range(nx-2,-1,-1):
        topo[i]=topo[i+1]+0.5*(S[i]+S[i+1])*dx;
        
    topo=topo-min(topo);

    slope=np.zeros(nx)
    kappa=np.ones(nx)
    halfway=round(0.5*nx)
    kappa[0:halfway]=kappa1*kappa[0:halfway]
    kappa[halfway:nx]=kappa2*kappa[halfway:nx]
    topoold=np.zeros(nx)

    t=0
    dt=0.05*dx/(max(kappa)*(area[nx-1])**(1/2))
    while t<tend:
    
        topoold=topo[:]
    
        slope[0:nx-1]=1/dx*abs((topo[1:nx]-topo[0:nx-1]))
        slope[nx-1]=slope[nx-2]
        erode=kappa*slope*area**(1/2)
    
        topo=topo+dt*uplift-dt*erode
        topo[nx-1]=topoold[nx-1]
    
        if (t>0.001) and (deltaz>0):
            topo[nx-1]=topo[nx-1]-deltaz
            deltaz=0;
    
        t=t+dt

    return(xgrid,area,topo,slope)


