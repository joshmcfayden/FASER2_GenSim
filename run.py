from generate_forsee_events import ForeseeGenerator
import numpy as np
import os
import gc # import garbage collector interface
import shutil
import subprocess
from array import array

import utils,setups

runmode="run combine"
#runmode="run"
#runmode="combine"
#runmode="G4"
#runmode="eff"
#runmode="plotsep"

do_hepmc=False
#do_hepmc=True

# Set model inputs
model="DarkPhoton"
#model="DarkHiggs"
masses,couplings,decays=utils.get_model_setup(model)

# Set run parameters
energy="14"
t0=0
notime=False
lumi=3000.*1000. #3000 fb-1 (sample is in pb)
nevents=10000

randomSeed=1
suffix=f"s{randomSeed}"

# Set paths
currdir=os.getcwd()
outdir="DUMMY"
#path="/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/FORESEE/"
path="/its/home/jam58/FASER2/FASER2_GenSim/FORESEE/"
G4path="../G4_test"

# Get details of setups
setup_dict=setups.setup_dict



# Loop over setups

ndettot=len(setup_dict)

for ndet,setup_name in enumerate(setup_dict):

    outdir=setup_name+"_"+model+"/"
    setup=setup_dict[setup_name]

    # make empty csv file (will append to this later
    utils.clear_csvs(runmode,setup,currdir,outdir,energy)
   
    print(f"\n\n=== Detector setup {setup_name} === ({ndet+1}/{ndettot})") 


    if not os.path.exists(outdir):
        os.makedirs(outdir)

    
    m_c_nevents=[]

    nmasstot=len(masses)

    for nmass,mass in enumerate(masses):

        c_nevents=[]
        ncouptot=len(couplings)
        
        for ncoup,coup in enumerate(couplings):
    
            nsignal=0.



            ndectot=len(decays)
            
            for ndec,decay in enumerate(decays):
                
                print(f"\nGenerating {model} events at Ecom = {energy}, mass = {mass} GeV ({nmass+1}/{nmasstot}), coupling = {coup} ({ncoup+1}/{ncouptot}), decay = {decay} ({ndec+1}/{ndectot}) [{(nmass*ncouptot*ndectot)+(ncoup*ndectot)+(ndec)+1}/{nmasstot*ncouptot*ndectot}]")
                               
                npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}.npy"
                hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}_{suffix}.hepmc"
                

                if "run" in runmode:
                    f = None

                    # Only run if npy file doesn't already exist
                    if os.path.exists(npname):
                        print(f"File {npname} found, moving to next step...")

                    else:
                        print(f"   Working on {npname}")
                        f = ForeseeGenerator(model, energy, mass, [coup], decay, outdir = outdir, path = path, randomSeed = randomSeed, t0 = t0, notime = notime, suffix = suffix, selection=setup["selection"], length=setup["length"], distance=setup["distance"])
                        f.write()
                       
                    # Run HepMC creation
                    if do_hepmc:
                        if os.path.exists(hepname):
                            print(f"File {hepname} found, skipping to next")
                        else:
                            print(f"   Working on {hepname}")

                            skip=False
                            if os.path.exists(npname):
                                np_arr=np.load(npname)
                                energies=np_arr[0]
                                thetas=np_arr[1]
                                weights=np_arr[2]
                                if not len(weights) or max(weights) == 0.:
                                    print("WARNING: Max weight = 0, not trying HepMC generation")
                                    skip=True

                            if not skip:
                                if not f:
                                    f = ForeseeGenerator(model, energy, mass, [coup], decay, outdir = outdir, path = path, randomSeed = randomSeed, t0 = t0, notime = notime, suffix = suffix, selection=setup["selection"], length=setup["length"], distance=setup["distance"])

                                try:
                                    hepname=f.write_hepmc(nevents)
                                    hepname=outdir+"/"+hepname
                                    print("saving hepmc file:",hepname)
                                except:
                                    print("Error in HepMC output creation")
                            

                # Sum up all the signal yields from all the decay modes and put in one .npy file (this is what plot_reach() expects)
                if "combine" in runmode:

                    if not os.path.exists(npname):
                        print(f"WARNING: File {npname} not found, skipping to next")
                        continue
                    
                    np_arr=np.load(npname)
                    energies=np_arr[0]
                    thetas=np_arr[1]
                    weights=np_arr[2]
                
                    nsignal+=sum(weights)*lumi
                    print(f"   nsignal: {nsignal}")

                    # trying to reduce memory usage
                    del np_arr
                    #np_arr.close()
                    gc.collect()
                



                # Run G4 simulation for how every many setups are defined
                if runmode=="G4":                       

                    if "G4" not in setup:
                        print(f"Error: No G4 configuration found for {setup_name}")
                        continue
                    
                    if os.path.exists(hepname):
                        print(f"Running G4 on {hepname}...")

                        # loop over setups (could be e.g. different tracker positions)
                        for G4setup in setup["G4"]:

                            outroot=f"{currdir}/{hepname}".replace('.hepmc',f'_{G4setup}.root')
                            if os.path.exists(outroot):
                                print(f"File {outroot} found, skipping to next")
                                continue

                            
                            G4dir=G4path+"/"+G4setup+"-build"
                            print("Running G4, moving to",G4dir)
                            os.chdir(G4dir)

                            
                            infilename="foresee_hepmc_ascii_1M_tmp.in"
                            infiletext=f'''/control/execute vis.mac
/vis/scene/add/magneticField 100
/generator/select hepmcAscii
/generator/hepmcAscii/open {currdir}/{hepname}
/generator/hepmcAscii/verbose 0
/run/beamOn 1000000     
'''                     
                            infile = open(infilename,'w')
                            infile.write(infiletext+'\n')
                            infile.close()
                        
                            with open(f'{currdir}/{outdir}/G4_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}_{suffix}_{G4setup}.log', 'w') as f:
                                process = subprocess.Popen(["./"+G4setup, infilename], stdout=f)
                                process.wait()
                                print("   ...G4 done!")
                        
                                print(f"Copying output.root to {outroot}")
                                shutil.copyfile("output.root",outroot)

                            os.chdir(currdir)
                    
                    else:
                        print(f"Not running G4, {hepname} not found")


                # Plot particle separations for each of the three stations
                if runmode=="plotsep":
               
                    for G4setup in setup["G4"]:
                        
                        utils.plot_seps(currdir,outdir,energy,mass,coup,decay,G4setup,setup_name)

            c_nevents.append(nsignal)
                        

            # Read G4 output files to get efficiencies for different cuts in each configuration
            if runmode=="eff":

                if "G4" not in setup:
                    print(f"Error: No G4 configuration found for {setup_name}")
                    continue
                
                # Do this for each G4 setup
                for G4setup in setup["G4"]:

                    # Go station by station
                    for station in setup["stations"]:
                    
                        # Calculate efficiencies from cuts
                        for eff in setup["effs"]:
    
                            sumall=0. # denominator for eff calc
                            sumeffs=[0. for x in setup["effs"][eff]] # numerator(s) for eff calc
                            
                            for decay in decays:
                                hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}_{suffix}.hepmc"
                                outroot=f"{currdir}/{hepname}".replace('.hepmc',f'_{G4setup}.root')
    
                                sumall,sumeffs=utils.get_effs(hepname,sumall,sumeffs,outroot,decay,station,setup["effs"][eff])
                                
                            # open output file to append results to
                            effname=f"{currdir}/{outdir}/eff_{energy}TeV_{G4setup}_station{station}_{eff}.csv"
                            print(f"Opening eff file {effname}")
                            efffile = open(effname,'a')
    
                            for n,(effval,efftitle,effstring) in enumerate(setup["effs"][eff]):
                                effcalc=sumeffs[n]/sumall if sumall > 0 else 1.0
                                print(f'{mass},{coup},{effval},{effcalc}\n')
                                efffile.write(f'{mass},{coup},{effval},{effcalc}\n')
                                




                                  

        m_c_nevents.append(c_nevents)
    
    # Only want one file per setup
    if "combine" in runmode:
        outdir="FORESEE/Models/"+model+"/model/results/"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        outfile=outdir+"/"+energy+"TeV_"+setup_name+".npy"
        
        print("\nWriting output file:",outfile)
        print("m:",len(masses),"c:",len(couplings),"n:",len(m_c_nevents[0]))
        np.save(outfile,np.array([masses,couplings,m_c_nevents],dtype=object))    

    
    # Go back to starting directory 
    os.chdir(currdir)
