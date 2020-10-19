import numpy as np
from gwf import Workflow
import os, sys
import itertools

gwf = Workflow()

def nquireRef(ploidy=2,depth=50,L=5000,nInd=1,random_seed=1):
    outFolder='./ploidy_ref/ploidy_'+str(ploidy)
    inputs=[]
    outputs=[outFolder+'/sim.DP'+str(depth)+'.genolikes']
    options = {
        'cores': 1,
        'memory': '24g',
        'walltime': '8:00:00'
    }

    spec = '''
    mkdir -p {outFolder} 
    #source activate HMMploidyNew
    echo {ploidy}x{nInd} > {outFolder}/ploidy_{ploidy}.file
    ./HMMploidy/simulationScript.sh -p {outFolder}/ploidy_{ploidy}.file -d {depth} -l {L} -g {random_seed} -o {outFolder}/sim
    echo ./{outFolder}/sim.DP{depth}.mpileup.gz > {outFolder}/file.list
    if [ {ploidy}==1 ]; then
        python3 ./HMMploidy/Genotype_Likelihoods.py {outFolder}/file.list -m 0
    else
        python3 ./HMMploidy/Genotype_Likelihoods.py {outFolder}/file.list
    fi
    gunzip {outFolder}/sim.DP{depth}.genolikes.gz
    '''.format(ploidy=ploidy, nInd=nInd, outFolder=outFolder, depth=depth, L=L, random_seed=random_seed)
    
    return inputs, outputs, options, spec

def simulate(ploidy=2, nInd=1, depth=10, L=5000, sim=1, random_seed=1, walltime='1:00:00', memory='8g'):
    inputs = []
    outFolder='ploidy_'+str(ploidy)+'/depth_'+str(depth)+'/nind_'+str(nInd)+'/run_'+str(sim)
    outputs = [outFolder+'/sim.DP'+str(depth)+'.genolikes']
    options = {
        'cores': 1,
        'memory': memory,
        'walltime': walltime
    }

    outFolder='ploidy_'+str(ploidy)+'/depth_'+str(depth)+'/nind_'+str(nInd)+'/run_'+str(sim)
    
    spec = '''
    mkdir -p {outFolder} 
    #source activate HMMploidyNew
    echo {ploidy}x{nInd} > {outFolder}/ploidy.file
    ./HMMploidy/simulationScript.sh -p {outFolder}/ploidy.file -d {depth} -l {L} -g {random_seed} -o {outFolder}/sim
    echo ./{outFolder}/sim.DP{depth}.mpileup.gz > {outFolder}/file.list
    if [ {ploidy}==1 ]; then
        python3 ./HMMploidy/Genotype_Likelihoods.py {outFolder}/file.list -m 0
    else
        python3 ./HMMploidy/Genotype_Likelihoods.py {outFolder}/file.list
    fi
    gunzip {outFolder}/sim.DP{depth}.genolikes.gz
    '''.format(ploidy=ploidy, nInd=nInd, outFolder=outFolder, depth=depth, L=L, random_seed=random_seed)

    return inputs, outputs, options, spec

def runHMM(ploidy=2, nInd=1, depth=10, sim=1, walltime='1:00:00', memory='8g'):
    inFolder='ploidy_'+str(ploidy)+'/depth_'+str(depth)+'/nind_'+str(nInd)+'/run_'+str(sim)
    inputs = [inFolder+'/sim.DP'+str(depth)+'.genolikes']
    outputs = [inFolder+'/sim.DP'+str(depth)+'.HMMploidy',
               inFolder+'/sim.DP'+str(depth)+'.pdf',
               inFolder+'/sim.DP'+str(depth)+'.testRate']
    options = {
        'cores': 1,
        'memory': memory,
        'walltime': walltime
    }

    spec = '''
    mkdir -p {outFolder} 
    #source activate HMMploidyNew   
    Rscript ./HMMploidy/HMMploidy.R file={outFolder}/sim.DP{depth} maxPloidy=5 truePl={ploidy} wind=250
    '''.format(outFolder=inFolder, depth=depth, ploidy=ploidy)

    return inputs, outputs, options, spec

def runPloidyNGS(ploidy=2, nInd=1, depth=10, sim=1, walltime='1:00:00', memory='8g'):
    inFolder='ploidy_'+str(ploidy)+'/depth_'+str(depth)+'/nind_'+str(nInd)+'/run_'+str(sim)
    inputs = [inFolder+'/sim.DP'+str(depth)+'.genolikes']
    outputs = [inFolder+'/sim.DP'+str(depth)+'.ploidyNGS']
    options = {
        'cores': 1,
        'memory': memory,
        'walltime': walltime
    }

    spec = '''
    mkdir -p {outFolder} 
    #source activate HMMploidyNew   
    Rscript ./HMMploidy/tools/ploidyNGS.R {outFolder}/sim.DP{depth} {NIND} {ploidy}
    '''.format(outFolder=inFolder, depth=depth, ploidy=ploidy, NIND=nInd)

    return inputs, outputs, options, spec

def runNquire(ploidy=2, nInd=1, depth=10, sim=1, walltime='2:00:00', memory='16g'):
    inFolder='ploidy_'+str(ploidy)+'/depth_'+str(depth)+'/nind_'+str(nInd)+'/run_'+str(sim)
    inputs = [inFolder+'/sim.DP'+str(depth)+'.genolikes']
    outputs = [inFolder+'/sim.DP'+str(depth)+'.noisy.gmmu']
    options = {
        'cores': 1,
        'memory': memory,
        'walltime': walltime
    }

    spec = '''
    mkdir -p {outFolder} 
    #source activate HMMploidyNew   
    Rscript ./HMMploidy/tools/nQuire.R {outFolder}/sim.DP{depth} {NIND} {ploidy}
    '''.format(outFolder=inFolder, depth=depth, ploidy=ploidy, NIND=nInd)

    return inputs, outputs, options, spec




PLOIDY = [1, 2, 3, 4, 5]
for p in PLOIDY:
    gwf.target_from_template('NQref_'+str(p), nquireRef(ploidy=p,
                                                nInd=1,
                                                depth=50,
                                                L=5000,
                                                random_seed=p))



DEPTH = [0.5, 1, 2, 5]
IND = [1, 2, 5, 10, 20]
PLOIDY = [1, 2, 3, 4, 5]
RUN = range(0,100) 
parameter_space = itertools.product(RUN,DEPTH,IND,PLOIDY)
for run, depth, ind, ploidy in parameter_space:
    name = 'depth_'+str(depth)+'_ind_'+str(ind)+'_pl_'+str(ploidy)+'_sim_'+str(run)
    gwf.target_from_template(name, simulate(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                L=5000,
                                                sim=run,
                                                random_seed=run,
                                                walltime='8:00:00',
                                                memory='64g'))
    gwf.target_from_template('HMM_'+name, runHMM(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='6:00:00'))    
    gwf.target_from_template('ploidyNGS_'+name, runPloidyNGS(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='6:00:00'))
    gwf.target_from_template('nQuire_'+name, runNquire(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='3:00:00'))




DEPTH = [10, 20]
IND = [1, 2, 5]
PLOIDY = [1, 2, 3]
RUN = range(0,100) 
parameter_space = itertools.product(RUN,DEPTH,IND,PLOIDY)
for run, depth, ind, ploidy in parameter_space:
    name = 'depth_'+str(depth)+'_ind_'+str(ind)+'_pl_'+str(ploidy)+'_sim_'+str(run)
    gwf.target_from_template(name, simulate(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                L=5000,
                                                sim=run,
                                                random_seed=run,
                                                walltime='8:00:00',
                                                memory='64g'))
    gwf.target_from_template('HMM_'+name, runHMM(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='14:00:00'))
    gwf.target_from_template('ploidyNGS_'+name, runPloidyNGS(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='6:00:00'))
    gwf.target_from_template('nQuire_'+name, runNquire(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='3:00:00'))



DEPTH = [10, 20]
IND = [1, 2, 5]
PLOIDY = [4, 5]
RUN = range(0,100) 
parameter_space = itertools.product(RUN,DEPTH,IND,PLOIDY)
for run, depth, ind, ploidy in parameter_space:
    name = 'depth_'+str(depth)+'_ind_'+str(ind)+'_pl_'+str(ploidy)+'_sim_'+str(run)
    gwf.target_from_template(name, simulate(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                L=5000,
                                                sim=run,
                                                random_seed=run,
                                                walltime='20:00:00',
                                                memory='64g'))
    gwf.target_from_template('HMM_'+name, runHMM(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='14:00:00'))
    gwf.target_from_template('ploidyNGS_'+name, runPloidyNGS(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='6:00:00'))
    gwf.target_from_template('nQuire_'+name, runNquire(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='3:00:00'))

    


DEPTH = [10, 20]
IND = [10, 20]
PLOIDY = [1, 2, 3]
RUN = range(0,100) 
parameter_space = itertools.product(RUN,DEPTH,IND,PLOIDY)
for run, depth, ind, ploidy in parameter_space:
    name = 'depth_'+str(depth)+'_ind_'+str(ind)+'_pl_'+str(ploidy)+'_sim_'+str(run)
    gwf.target_from_template(name, simulate(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                L=5000,
                                                sim=run,
                                                random_seed=run,
                                                memory='64g',
                                                walltime='20:00:00'))
    gwf.target_from_template('HMM_'+name, runHMM(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='14:00:00'))
    gwf.target_from_template('ploidyNGS_'+name, runPloidyNGS(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='6:00:00'))
    gwf.target_from_template('nQuire_'+name, runNquire(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='6:00:00'))




DEPTH = [10, 20]
IND = [10, 20]
PLOIDY = [4, 5]
RUN = range(0,100) 
parameter_space = itertools.product(RUN,DEPTH,IND,PLOIDY)
for run, depth, ind, ploidy in parameter_space:
    name = 'depth_'+str(depth)+'_ind_'+str(ind)+'_pl_'+str(ploidy)+'_sim_'+str(run)
    gwf.target_from_template(name, simulate(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                L=5000,
                                                sim=run,
                                                random_seed=run,
                                                memory='64g',
                                                walltime='24:00:00'))
    gwf.target_from_template('HMM_'+name, runHMM(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='14:00:00'))
    gwf.target_from_template('ploidyNGS_'+name, runPloidyNGS(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='6:00:00'))
    gwf.target_from_template('nQuire_'+name, runNquire(ploidy=ploidy,
                                                nInd=ind,
                                                depth=depth,
                                                sim=run,
                                                memory='24g',
                                                walltime='2:00:00'))
