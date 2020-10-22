"""
Simple Molecular Dynamics code for course FYS-4096 Computational Physics

Problem 1:
- Make the code to work and solve H2 using the Morse potential.
- Modify and comment especially at parts where it reads "# ADD"
- Follow instructions on ex12.pdf

Problems 2:
- Add observables: temperature, distance, and heat capacity.
- Follow instructions on ex12.pdf

Problem 3:
- Add Lennard-Jones capability etc.
- Follow instructions on ex12.pdf

Problem 4:
- Temperature dependent data and analysis
- Follow instructions on ex12.pdf

"""




from numpy import *
from matplotlib.pyplot import *

class Atom:
    def __init__(self,index,mass,dt,dims):
        self.index = index
        self.mass = mass
        self.dt = dt
        self.dims = dims
        self.LJ_epsilon = None
        self.LJ_sigma = None
        self.R = None
        self.v = None
        self.force = None

    def set_LJ_parameters(self,LJ_epsilon,LJ_sigma):
        self.LJ_epsilon = LJ_epsilon
        self.LJ_sigma = LJ_sigma

    def set_position(self,coordinates):
        self.R = coordinates
        
    def set_velocity(self,velocity):
        self.v = velocity

    def set_force(self,force):
        self.force = force

class Observables:
    def __init__(self):
        self.E_kin = []
        self.E_pot = []
        self.distance = []
        self.Temperature = []
        self.Capacity = []
        self.E2  = []

def calculate_energetics(atoms):
    N = len(atoms)
    V = 0.0
    E_kin = 0.0
    # ADD calculation of kineti and potential energy
    for i in range(N):
        E_kin += 0.5*(atoms[i].mass*sum(atoms[i].v**2))
        for j in range(i+1,N):
            V += pair_potential(atoms[i],atoms[j]) 
    
    
    
    return E_kin, V

def calculate_force(atoms):
    # ADD comments, e.g., what are the elements of the return function F
    # 
    N = len(atoms)
    ij_map = zeros((N,N),dtype=int)
    Fs = []
    ind = 0
    # Get the map for index with upper and lower triangle order. 
    for i in range(0,N-1):
        for j in range(i+1,N):
            Fs.append(pair_force(atoms[i],atoms[j]))
            ij_map[i,j] = ind
            ij_map[j,i] = ind
            ind += 1
    F = []
    for i in range(N):
        f = zeros(shape=shape(atoms[i].R))
        for j in range(N):
            ind = ij_map[i,j]
            if i<j:
                f += Fs[ind] #for rji, when j>i, sum Fji
            elif i>j:
                f -= Fs[ind] #for rji, when j<i, the paired -Fji is sumed as well
        F.append(f) 
    F = array(F) 
    return F # get the array of total force F

def pair_force(atom1,atom2):
    return Morse_force(atom1,atom2)

def pair_potential(atom1,atom2):
    return Morse_potential(atom1,atom2)

def Morse_potential(atom1,atom2):
    # H2 parameters given here
    De = 0.1745
    re = 1.40
    a = 1.0282
    # ADD comments and calculation of Morse potential
    r = sqrt(sum((atom1.R-atom2.R)**2))
    V = De*((1-exp(-a*(r-re)))**2)-De
    #
    return V

def Morse_force(atom1,atom2):
    # H2 parameters
    De = 0.1745
    re = 1.40
    a = 1.0282
    # ADD comments and calculation of Morse force
    r = sqrt(sum((atom1.R-atom2.R)**2))
    r_hat = (atom1.R-atom2.R)/r
    F = -2*De*a*exp(-a*(r-re))*(1-exp(-a*(r-re)))*r_hat
    return F

def lennard_jones_potential(atom1,atom2):
    epsilon = sqrt(atom1.LJ_epsilon*atom2.LJ_epsilon)
    sigma = (atom1.LJ_sigma+atom2.LJ_sigma)/2
    # ADD --> in problem 3: comments and calculation of LJ potential
    r = sqrt(sum((atom1.R-atom2.R)**2))
    V = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    #
    return V

def lennard_jones_force(atom1,atom2):
    epsilon = sqrt(atom1.LJ_epsilon*atom2.LJ_epsilon)
    sigma = (atom1.LJ_sigma+atom2.LJ_sigma)/2
    # ADD --> in problem 3: comments and calculation of LJ force
    r = sqrt(sum((atom1.R-atom2.R)**2))
    r_hat = (atom1.R-atom2.R)/r
    F = -4*epsilon*((-12/r)*(sigma/r)**12+(6/r)*(sigma/r)**6)*r_hat
    return F

def velocity_verlet_update(atoms):
    # ADD comments and integration algoritm according to function name
    dt = atoms[0].dt
    dt2 = dt**2
    for i in range(len(atoms)):
        atoms[i].R += dt*atoms[i].v+dt2/(2*atoms[i].mass)*atoms[i].force
    Fnew = calculate_force(atoms)
    for i in range(len(atoms)):
        atoms[i].v += dt/(2*atoms[i].mass)*(Fnew[i]+atoms[i].force)
        atoms[i].force = Fnew[i] # update force
    return atoms
    
def initialize_positions(atoms):
    # diatomic case
    atoms[0].set_position(array([-0.8,0.0,0.0]))
    atoms[1].set_position(array([0.7,0.0,0.0]))

def initialize_velocities(atoms):
    # ADD --> in problem 4: comment on what is v_max
    # diatomic case
    dims = atoms[0].dims
    kB=3.16e-6 # in hartree/Kelvin
    for i in range(len(atoms)):
        v_max = sqrt(3.0/atoms[i].mass*kB*10.0)
        atoms[i].set_velocity(array([1.0,0.0,0.0])*v_max)
    atoms[1].v = -1.0*atoms[0].v

def initialize_force(atoms):
    F=calculate_force(atoms)
    for i in range(len(atoms)):
        atoms[i].set_force(F[i])

def Temperature():
    # Boltzmann constant in Hartree/Kelvin
    kB = 3.16e-6
    # ADD calculation of temperature in problem 2
#    T = 2*observables.E_kin
 #   observables.tamperature.append(T)
    #
 #   return T

def calculate_observables(atoms,observables):
            # ADD calculation of temperature in problem 2

    
    E_k, E_p = calculate_energetics(atoms)
    
#    C = (E**2-)
    E_2 = dot(E_k+E_p,E_k+E_p)
    E_1 = E_k+E_p
#    print(E_2)
    
#    print(dot(mean(E_k+E_p),mean(E_k+E_p)))
#    print('Tot E')
#    print(E_2)
    dims = 3
    kB = 3.16e-6
    T = 2*mean(observables.E_kin)/(kB*dims*len(atoms))
#    C = (E_2-dot((E_k+E_p),mean(E_k+E_p)))/(kB*T)
#    print(C)
#    print('Temperature')
#    print(T)
#    print('E2')
#    print(E_2)

    for i in range(0,len(atoms)-1):
        for j in range(i+1, len(atoms)):
#            print(atoms[i].R-atoms[j].R)
#            print(sum(dot(atoms[i].R-atoms[j].R,atoms[i].R-atoms[j].R)))
            d = sqrt(sum(dot(atoms[i].R-atoms[j].R,atoms[i].R-atoms[j].R)))
#    print('distance')
#    print(d)
    observables.E_kin.append(E_k)
    observables.E_pot.append(E_p)
#    observables.Capacity.append(C)
    observables.distance.append(d)
    observables.Temperature.append(T)
    # ADD calculation of observables in problem 2
    
    #
    #
    return observables,d, T,E_1,E_2

def main():
    N_atoms = 2
    dims = 3
    dt = 0.1
    mass = 1860.0

    # Initialize atoms
    atoms = []    
    for i in range(N_atoms):
        atoms.append(Atom(i,mass,dt,dims))
        atoms[i].set_LJ_parameters(0.1745,1.25)
    # Initialize observables
    observables = Observables()

    # Initialize positions, velocities, and forces
    #T = observables.Temperature
    initialize_positions(atoms)
    initialize_velocities(atoms)
    initialize_force(atoms)

    E_k, E_p = calculate_energetics(atoms)
    
    T = [];
    E_1 = [];
    E_2 = [];
    d = [];
    for i in range(100000):
        atoms = velocity_verlet_update(atoms)
        if ( i % 10 == 0):
            observables,D,TT,E1,E2 = calculate_observables(atoms,observables)            
            T.append(TT)
            E_1.append(E1)
            E_2.append(E2)
            d.append(D) 

#            observables= calculate_observables(atoms,observables)            
    kB = 3.16e-6
    # Print energies
    E_kin = array(observables.E_kin)
    E_pot = array(observables.E_pot)
#    capacity = array(observables.Capacity)
#    distance = array(observables.distance)
    E_tot = E_kin+E_pot
    print(mean(E_tot))
    print(mean(E_2))
    C = (mean(E_2)-E_tot**2)/dot(kB,T)
#    T = array(observables.Temperature)  
    print('C',C)
    print('E_kin',mean(E_kin))
    print('E_pot',mean(E_pot))
    print('E_tot',mean(E_tot))
    print('T', mean(T))
    print('T_std', std(T))
    print("Capacity",mean(C))
    print("Distance",mean(d))
    
    figure()
    plot(E_tot)
    title('E_tot')
    savefig('E_tot.pdf',dpi = 200)
    figure()
    plot(E_kin)
    title('E_kint')
    savefig('E_kin.pdf',dpi = 200)
    figure()
    plot(E_pot)
    title('E_pot')
    savefig('E_pot.pdf',dpi = 200)
    figure()
    plot(T)
    title('Temperature with standard:',std(T))
    savefig('Temperature.pdf',dpi = 200)
    figure()
    plot(d)
    title('Distance')
    savefig('Distance.pdf',dpi = 200) 
    figure()
    plot(C)
    title('Capacity')
    savefig('Capacity.pdf',dpi = 200) 
    show()

if __name__=="__main__":
    main()
        
