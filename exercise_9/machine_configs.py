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
    qmc_presub ='''
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
    qe2='qmcpack.x'
    qe3='pw2qmcpack.x'
    
    # 4 processes
    scf  = job(cores=16,minutes=30,user_env=False,presub=scf_presub,app=qe)
    # 24 processes (1 node = 24 processors at taito)
    #scf  = job(nodes=1,hours=1,user_env=False,presub=scf_presub,app=qe)
    qmc = job(cores=16,minutes=30,threads=4,user_env=False,presurb=scf_presub,app=qe2)
    conv = job(cores=1, minutes=30,user_env=False,presub=scf_presub,app=qe3)
    
    jobs = {'scf' : scf, 'qmc': qmc, 'conv': conv}

    return jobs
