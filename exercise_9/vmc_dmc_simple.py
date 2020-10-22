"""
Simple VMC and DMC code for FYS-4096 course at TUT

- Fill in more comments based on lecture slides
- Follow assignment instructions  


By Ilkka Kylanpaa
"""

from numpy import *
from matplotlib.pyplot import *

class Walker:
    def __init__(self,*args,**kwargs):
        self.Ne = kwargs['Ne']
        self.Re = kwargs['Re']
        self.spins = kwargs['spins']
        self.Nn = kwargs['Nn']
        self.Rn = kwargs['Rn']
        self.Zn = kwargs['Zn']
        self.sys_dim = kwargs['dim']

    def w_copy(self):
        return Walker(Ne=self.Ne,
                      Re=self.Re,
                      spins=self.spins,
                      Nn=self.Nn,
                      Rn=self.Rn,
                      Zn=self.Zn,
                      dim=self.sys_dim)
    

def vmc_run(Nblocks,Niters,Delta,Walkers_in,Ntarget):
    Eb = zeros((Nblocks,))
    Accept=zeros((Nblocks,))

    vmc_walker=Walkers_in[0] # just one walker needed
    Walkers_out=[]
    for i in range(Nblocks):
        for j in range(Niters):
            # moving only electrons
            for k in range(vmc_walker.Ne):
                R = vmc_walker.Re[k]
                Psi_R = wfs(vmc_walker)

                # move the particle
                vmc_walker.Re[k] = R + Delta*(random.rand(vmc_walker.sys_dim)-0.5)
                
                # calculate wave function at the new position
                Psi_Rp = wfs(vmc_walker)

                # calculate the sampling probability
                A_RtoRp = min((Psi_Rp/Psi_R)**2,1.0)
                
                # Metropolis
                if (A_RtoRp > random.rand()):
                    Accept[i] += 1.0/vmc_walker.Ne
                else:
                    vmc_walker.Re[k]=R
                #end if
            #end for
            Eb[i] += E_local(vmc_walker)
        #end for
        if (len(Walkers_out)<Ntarget):
            Walkers_out.append(vmc_walker.w_copy())
        Eb[i] /= Niters
        Accept[i] /= Niters
        print('Block {0}/{1}'.format(i+1,Nblocks))
        print('    E   = {0:.5f}'.format(Eb[i]))
        print('    Acc = {0:.5f}'.format(Accept[i]))
    
    return Walkers_out, Eb, Accept


def dmc_run(Nblocks,Niters,Walkers,tau,E_T,Ntarget):
    
    max_walkers = 2*Ntarget
    lW=len(Walkers)
    while len(Walkers)<Ntarget:
        Walkers.append(Walkers[max(1,int(lW*random.rand()))].w_copy())

    Eb=zeros((Nblocks,))
    Accept=zeros((Nblocks,))
    AccCount=zeros((Nblocks,))
    
    obs_interval=5
    mass=1

    for i in range(Nblocks):
        EbCount=0
        for j in range(Niters):
            Wsize = len(Walkers)
            Idead = []
            for k in range(Wsize):
                for np in range(Walkers[k].Ne):
                    # initiate accuracy
                    Acc = 0.0
                    # move only the electrons
                    R = Walkers[k].Re[np]
                    Psi_R = wfs(Walkers[k])
                    DriftPsi_R = 2*Gradient(Walkers[k],np)/Psi_R*tau/2/mass
                    E_L_R=E_local(Walkers[k])

                    DeltaR=random.randn(Walkers[k].sys_dim)
                    logGf=-0.5*dot(DeltaR,DeltaR)
                    # move the particles
                    Walkers[k].Re[np]=R+DriftPsi_R+DeltaR*sqrt(tau/mass)
                    # calculate wave function at the new position
                    Psi_Rp = wfs(Walkers[k])
                    DriftPsi_Rp = 2*Gradient(Walkers[k],np)/Psi_Rp*tau/2/mass
                    E_L_Rp = E_local(Walkers[k])
                    
                    DeltaR = R-Walkers[k].Re[np]-DriftPsi_Rp
                    logGb = -dot(DeltaR,DeltaR)/2/tau*mass
                    #calculate the sampling probability
                    A_RtoRp = min(1, (Psi_Rp/Psi_R)**2*exp(logGb-logGf))
                    #Metropolis
                    if (A_RtoRp > random.rand()):
                        Acc += 1.0/Walkers[k].Ne
                        Accept[i] += 1
                    else:
                        Walkers[k].Re[np]=R
                    
                    AccCount[i] += 1
                
                tau_eff = Acc*tau
                GB = exp(-(0.5*(E_L_R+E_L_Rp) - E_T)*tau_eff)
                MB = int(floor(GB + random.rand()))
                
                if MB>1:
                    for n in range(MB-1):
                        if (len(Walkers)<max_walkers):
                            Walkers.append(Walkers[k].w_copy())
                elif MB==0:
                    Idead.append(k)
 
            Walkers = DeleteWalkers(Walkers,Idead)

            # Calculate observables every now and then
            if j % obs_interval == 0:
                EL = Observable_E(Walkers)
                Eb[i] += EL
                EbCount += 1
                E_T += 0.01/tau*log(Ntarget/len(Walkers))
                

        Nw = len(Walkers)
        dNw = Ntarget-Nw
        for kk in range(abs(dNw)):
            ind=int(floor(len(Walkers)*random.rand()))
            if (dNw>0):
                Walkers.append(Walkers[ind].w_copy())
            elif dNw<0:
                Walkers = DeleteWalkers(Walkers,[ind])

        
        Eb[i] /= EbCount
        Accept[i] /= AccCount[i]
        print('Block {0}/{1}'.format(i+1,Nblocks))
        print('    E   = {0:.5f}'.format(Eb[i]))
        print('    Acc = {0:.5f}'.format(Accept[i]))


    return Walkers, Eb, Accept

def DeleteWalkers(Walkers,Idead):
    if (len(Idead)>0):
        if (len(Walkers)==len(Idead)):
            Walkers = Walkers[0]
        else:
            Idead.sort(reverse=True)   
            for i in range(len(Idead)):
                del Walkers[Idead[i]]

    return Walkers

def H_1s(r1,r2):
    return exp(-sqrt(sum((r1-r2)**2)))
     
def wfs(Walker):
    # H2 approx
    f = H_1s(Walker.Re[0],Walker.Rn[0])+H_1s(Walker.Re[0],Walker.Rn[1])
    f *= (H_1s(Walker.Re[1],Walker.Rn[0])+H_1s(Walker.Re[1],Walker.Rn[1]))

    J = 0.0
    # Jastrow e-e
    for i in range(Walker.Ne-1):
        for j in range(i+1,Walker.Ne):
           r = sqrt(sum((Walker.Re[i]-Walker.Re[j])**2))          
           if (Walker.spins[i]==Walker.spins[j]):
            #Jastrrow factor used here: a1 = 0.25, b1 =1.0
                J += 0.25*r/(1.0+1.0*r)
           else:
            #Jastrrow factor used here: a1 = 0.5, b1 =1.0
               J += 0.5*r/(1.0+1.0*r)
    
    # Jastrow e-Ion
    for i in range(Walker.Ne):
        for j in range(Walker.Nn):
           r = sqrt(sum((Walker.Re[i]-Walker.Rn[j])**2))
           #Jastrrow factor used here: a1 = Zn, b1 =100
           J -= Walker.Zn[j]*r/(1.0+100.0*r)
       

    return f*exp(J)

def potential(Walker):
    V = 0.0
    r_cut = 1.0e-4
    # e-e
    for i in range(Walker.Ne-1):
        for j in range(i+1,Walker.Ne):
            r = sqrt(sum((Walker.Re[i]-Walker.Re[j])**2))
            V += 1.0/max(r_cut,r)

    # e-Ion
    for i in range(Walker.Ne):
        for j in range(Walker.Nn):
            r = sqrt(sum((Walker.Re[i]-Walker.Rn[j])**2))
            V -= Walker.Zn[j]/max(r_cut,r)

    # Ion-Ion
    for i in range(Walker.Nn-1):
        for j in range(i+1,Walker.Nn):
            r = sqrt(sum((Walker.Rn[i]-Walker.Rn[j])**2))
            V += 1.0/max(r_cut,r)

    return V

def Local_Kinetic(Walker):
    # laplacian -0.5 \nabla^2 \Psi / \Psi
    h = 0.001
    h2 = h*h
    K = 0.0
    Psi_R = wfs(Walker)
    for i in range(Walker.Ne):
        for j in range(Walker.sys_dim):
            Y=Walker.Re[i][j]
            Walker.Re[i][j]-=h
            wfs1 = wfs(Walker)
            Walker.Re[i][j]+=2.0*h
            wfs2 = wfs(Walker)
            K -= 0.5*(wfs1+wfs2-2.0*Psi_R)/h2
            Walker.Re[i][j]=Y
    return K/Psi_R

def Gradient(Walker,particle):
    h=0.001
    dPsi = zeros(shape=shape(Walker.Re[particle]))
    for i in range(Walker.sys_dim):
        Y=Walker.Re[particle][i]
        Walker.Re[particle][i]-=h
        wfs1=wfs(Walker)
        Walker.Re[particle][i]+=2.0*h
        wfs2=wfs(Walker)
        dPsi[i] = (wfs2-wfs1)/2/h
        Walker.Re[particle][i]=Y

    return dPsi

def E_local(Walker):
    return Local_Kinetic(Walker)+potential(Walker)

def Observable_E(Walkers):
    E=0.0
    Nw = len(Walkers)
    for i in range(Nw):
        E += E_local(Walkers[i])
    E /= Nw
    return E

def main():
    Walkers=[]
    Walkers.append(Walker(Ne=2,
                          Re=[array([0.5,0,0]),array([-0.5,0,0])],
                          spins=[0,1],
                          Nn=2,
                          Rn=[array([-0.7,0,0]),array([0.7,0,0])],
                          Zn=[1.0,1.0],
                          dim=3))

    vmc_only = True
    Ntarget=100
    vmc_time_step = 2.0

    Walkers, Eb, Acc = vmc_run(100,50,vmc_time_step,Walkers,Ntarget)
    n = np.array(linspace(1, 100, 100))
    Eb_mean = zeros(100,)
    Eb_std = zeros(100,)
    err_l = zeros(100,)
    err_u = zeros(100,)
    for i in range(100):
        a = int(n[i])
        Eb_mean[i] = mean(Eb[0:a])
        Eb_std[i] =std(Eb[0:a])
    figure(),
    plot(linspace(1,100,100), Eb,'b-', linspace(1,100,100), (mean(Eb)-std(Eb))*ones(100,), 'k--',linspace(1,100,100),  mean(Eb)*ones(100,), 'r-',linspace(1,100,100), (mean(Eb)+std(Eb))*ones(100,), 'k--')
    legend(['Energy','Standard Deviation','Total Energy'])
    title('Energy')
    savefig('Energy.pdf',dpi = 200)
    figure(),
    errorbar(linspace(1,100,100),Eb_mean, yerr=Eb_std)
    title('Mean Energy ')
    savefig('Mean_Energy.pdf',dpi = 200)
    figure(),
    plot(Acc)
    title('Acceptance Ratio')
    savefig('Acceptance.pdf',dpi = 200)
    if not vmc_only:
        Walkers, Eb_dmc, Accept_dmc = dmc_run(10,10,Walkers,0.05,mean(Eb),Ntarget)
    
if __name__=="__main__":
    main()
        
