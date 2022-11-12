from generate_forsee_events import ForeseeGenerator
import numpy as np
import os
import gc # import garbage collector interface
import shutil
import subprocess



#args: model, Ecom, mass, couplings, pid1, pid2, outdir, path, randomSeed, t0, notime, suffix
#args: DarkPhoton 14 0.1 [0.0001] 11 -11 mytest /Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/FORESEE/ 0 0 False test1

mode="run"
#mode="G4"
#mode="runcombine"
mode="combine"
#mode="eff"

do_hepmc=False
do_hepmc=True


model="DarkPhoton"
energy="14"
outdir="mytest2"
path="/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/FORESEE/"
t0=0
notime=False

lumi=3000.*1000. #3000 fb-1 (sample is in pb)

nevents=10000

randomSeed=1
suffix=f"s{randomSeed}"


G4path="../G4_test"

#masses = [ 0.01  ,  0.1,  1.0, 10. ]
masses = [ 
    0.01  ,
    0.0501,
    0.1585,
    0.3548,
    0.6457,
    0.7586,
    0.8913,
    1.2589,
    2.8184
]

#masses = [ 
#    0.01  ,  0.0126,  0.0158,  0.02  ,  0.0251,  0.0316,  0.0398,
#    0.0501,  0.0631,  0.0794,  0.1   ,  0.1122,  0.1259,  0.1413,
#    0.1585,
#


#    0.1778,  0.1995,  0.2239,  0.2512,  0.2818,  0.3162,
#    0.3548,  0.3981,  0.4467,  0.5012,  0.5623,  0.6026,  0.631 ,
#    0.6457,  0.6607,  0.6761,  0.6918,  0.7079,  0.7244,  0.7413,
#    0.7586,  0.7762,  0.7943,  0.8128,  0.8318,


#    0.8511,  0.871 ,
#    0.8913,  0.912 ,  0.9333,  0.955 ,  0.9772,  1.    ,  1.122 ,
#    1.2589,  1.4125,  1.5849,  1.7783,  1.9953,  2.2387,  2.5119,
#    2.8184,  3.1623,  3.9811,  5.0119,  6.3096,  7.9433, 10.    
#]


#couplings=np.logspace(-8,-3,12)
couplings=np.logspace(-8,-3,20)


#pids=[[-11, 11]]
pids=[[-11, 11],[-13,13],[999,999]]
#pids=[[-13,13],[999,999]]



setup_dict={
    "F2-default":{
        "name":"FASER2 Orig",# (Default)",
        #"color":"maroon",
        "color":"firebrick",
        "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
        "length":5,
        "distance":480,
        "channels": None,
        "G4":["FASER2_HepMC_v4_FASER2_Default_3rdTrkStation"],
        "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_y-em_y)>0.1'),
                       (1,'sep>1mm','abs(ep_y-em_y)>1'),
                       (5,'sep>5mm','abs(ep_y-em_y)>5'),
                       (10,'sep>10mm','abs(ep_y-em_y)>10'),
                       (100,'sep>100mm','abs(ep_y-em_y)>10')]}
    },

    "S2-L10-D2":{
        #"name":"S2 L=10m D=2m",
        "name":"Old Baseline",
        "color":"royalblue",
        #"color":"darkgreen",
        #"color":"orange",
        "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
        "length":10,
        "distance":615,
        "channels": None
    },

        "R1-L10-R0p5x2":{
            #"name":"R1 L=10m X=2m Y=0.5m",
            "name":"New Baseline (X=2m Y=0.5m)",
            #"color":"lightsteelblue",
            "color":"pink",
            "style":"dotted",
            "selection":"(np.sqrt(x.x**2)<0.25) * (np.sqrt(x.y**2)<1.0)",
            "length":10,
            "distance":615,
            "channels": None
        },

    
        "R1-L10-R1x3":{
            #"name":"R1 L=10m X=3m Y=1m",
            "name":"New Baseline (X=3m Y=1m)",
            #"color":"cornflowerblue",
            "color":"forestgreen",
            #"color":"limegreen",
            "style":"dashed",
            "selection":"(np.sqrt(x.x**2)<0.5) * (np.sqrt(x.y**2)<1.5)",
            "length":10,
            "distance":615,
            "channels": None
        },
        "R1-L10-R0p5x3":{
            #"name":"R1 L=10m X=3m Y=0.5m",
            "name":"New Baseline (X=3m Y=0.5m)",
            #"color":"lightsteelblue",
            "color":"orchid",
            #"color":"limegreen",
            "style":"dotted",
            "selection":"(np.sqrt(x.x**2)<0.25) * (np.sqrt(x.y**2)<1.5)",
            "length":10,
            "distance":615,
            "channels": None
        },

    
}


currdir=os.getcwd()

#clear existing eff files




for setup_name in setup_dict:

    outdir=setup_name+"_"+model+"/"
    setup=setup_dict[setup_name]

    if mode=="eff":
        if "G4" not in setup:
            continue
        for G4setup in setup["G4"]:
            for eff in setup["effs"]:
                effname=f"{currdir}/{outdir}/eff_{energy}TeV_{G4setup}_{eff}.csv"
                print(f"Clearing eff file {effname}")
                efffile = open(effname,'w')
                efffile.close()

    
    print(f"\n\n=== Detector setup {setup_name} ===") 

        
    m_c_nevents=[]
    for mass in masses:
        c_nevents=[]
        for coup in couplings:
    
            nsignal=0.




            # Read G4 output files to get efficiencies for each configuration
            if mode=="eff":

                import ROOT
                import pyhepmc

                if "G4" not in setup:
                    print(f"Error: No G4 configuration found for {setup_name}")
                    continue
                

                for G4setup in setup["G4"]:
                    
                    for eff in setup["effs"]:


                        sumall=0.
                        sumeffs=[0. for x in setup["effs"][eff]]
                        
                        for pid1,pid2 in pids:
                            hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}_{suffix}.hepmc"
                            outroot=f"{currdir}/{hepname}".replace('.hepmc',f'_{G4setup}.root')

                            # get xs from hepmc
                            xs=0.
                            with pyhepmc.open("/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/mytest_G4_mini/events_14TeV_m0.1GeV_c1e-05to_11_-11_test1.hepmc") as f:
                                event = f.read()
                                xs=event.cross_section.xsec()
                                print("Found cross section:",xs)
        
                            if not os.path.exists(outroot):
                                print(f"ROOT file {outroot} not found - skipping")
                                continue
                            
                            f_eff=ROOT.TFile.Open(outroot,"READ")
                            t_eff = f_eff.Get("Hits")

                            var='ep_y'
                            if abs(pid1) == 13:
                                var='mp_y'
                            elif abs(pid1) == 999:
                                var='hp_y'

                            print(f"Using var {var}")


                            hall=ROOT.TH1F("hall","hall",1000,0,10000)                                
                            t_eff.Draw(f"{var}>>hall");
                            sumall+=hall.Integral()*xs
                            print(f"pid1 = {pid1}, all =",hall.Integral(),xs,hall.Integral()*xs)
                                
                            for n,(effval,efftitle,effstring) in enumerate(setup["effs"][eff]):

                                heff=ROOT.TH1F("heff","heff",1000,0,10000)
                                t_eff.Draw(f"{var}>>heff",effstring)

                                print(f"pid1 = {pid1}, effs[{n}] =",heff.Integral(),xs,heff.Integral()*xs)
                                sumeffs[n]+=heff.Integral()*xs

                                
                        # open output file to append results to
                        effname=f"{currdir}/{outdir}/eff_{energy}TeV_{G4setup}_{eff}.csv"
                        print(f"Opening eff file {effname}")
                        efffile = open(effname,'a')

                        for n,(effval,efftitle,effstring) in enumerate(setup["effs"][eff]):
                            effcalc=sumeffs[n]/sumall if sumall > 0 else 1.0
                            print(f'{mass},{coup},{effval},{effcalc}\n')
                            efffile.write(f'{mass},{coup},{effval},{effcalc}\n')
                             
                        efffile.close()

            
            for pid1,pid2 in pids:
                
                #print(f"Generating {model} events at Ecom = {energy}") 
                #print(f"   mother mass = {mass} GeV")
                #print(f"   decay = {pid1} {pid2}")
                #print(f"   couplings = {coup}")    

                print(f"\nGenerating {model} events at Ecom = {energy}, mass = {mass} GeV, decay = {pid1} {pid2}, coupling = {coup}")
                
                npname=""
                
                #break
                if "run" in mode:
                    npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}.npy"
                    if os.path.exists(npname):
                        print(f"File {npname} found, skipping to next")
                        continue
                    
                    f = ForeseeGenerator(model, energy, mass, [coup], pid1, pid2, outdir = outdir, path = path, randomSeed = randomSeed, t0 = t0, notime = notime, suffix = suffix, selection=setup["selection"], length=setup["length"], distance=setup["distance"])
                
                    npname=f.write()

                if "combine" in mode:
                    if "run" not in mode:
                        npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}.npy"

                        print("   Using file:",npname)


                    np_arr=np.load(npname)
                    energies=np_arr[0]
                    thetas=np_arr[1]
                    weights=np_arr[2]
                
                    nsignal+=sum(weights)*lumi
                    print(f"   nsignal: {nsignal}")
                
                    del np_arr
                    #np_arr.close()
                    gc.collect()
                
                #except:
                #    print("Error in numpy output creation")
                    
                hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}_{suffix}.hepmc"

                if do_hepmc:
                    try:
                        hepname=f.write_hepmc(nevents)
                        hepname=outdir+"/"+hepname
                        print("saving hepmc file:",hepname)
                    except:
                        print("Error in HepMC output creation")


                # Run G4 simulation for how every many setups are defined
                if mode=="G4":                       

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
                        
                            with open(f'{currdir}/{outdir}/G4_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}_{suffix}_{G4setup}.log', 'w') as f:
                                process = subprocess.Popen(["./"+G4setup, infilename], stdout=f)
                                process.wait()
                                print("   ...G4 done!")
                        
                                print(f"Copying output.root to {outroot}")
                                shutil.copyfile("output.root",outroot)

                            os.chdir(currdir)
                    
                    else:
                        print(f"Not running G4, {hepname} not found")


                        



            c_nevents.append(nsignal)
                                  


        m_c_nevents.append(c_nevents)
    

    if "combine" in mode:
        outdir="FORESEE/Models/"+model+"/model/results/"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        outfile=outdir+"/"+energy+"TeV_"+setup_name+".npy"
        
        print("\nWriting output file:",outfile)
        #np.save(outfile,[masses,couplings,m_c_nevents])
        np.save(outfile,np.array([masses,couplings,m_c_nevents],dtype=object))    

    
    # Go back to starting directory 
    os.chdir(currdir)
