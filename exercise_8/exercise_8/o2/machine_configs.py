#! /usr/bin/env python
from nexus import job

def general_configs(machine):
    if machine=='taito':
        jobs = get_taito_configs()
    else:
        print 'Using taito as defaul machine'
        jobs = get_taito_configs()
    return jobs

def get_taito_configs():
    scf_presub = '''
    module purge
    module load gcc
    module load openmpi
    module load openblas
    module load hdf5-serial
    '''

    pmc_presub = '''
    module purge
    module load gcc
    module load mkl
    module load intelmpi
    module load hdf5-par
    module load fftw
    module load boost
    module load cmake   
    '''

    qe='pw.x'
    qe2='qmcpack_taito_cpu_cpu_comp_SoA'
    # 4 processes
    scf  = job(cores=4,minutes=30,user_env=False,presub=scf_presub,app=qe)
    # 24 processes (1 node = 24 processors at taito)
    #scf  = job(nodes=1,hours=1,user_env=False,presub=scf_presub,app=qe)
    vmc = job(cores = 16, minutes = 10, user_env +False,presub = qmc_presub,qpp =qe2)
    
    jobs = {'scf' : scf, 'vmc': vmc}

    return jobs
