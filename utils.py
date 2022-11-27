import numpy as np
import ROOT
import pyhepmc
import os

def get_model_setup(model,setup=None):

    masses = couplings = decays = None

    if model=="DarkPhoton":
    
        #pids=[[-11, 11],[-13,13],[999,999]]
        decays=["e_e","mu_mu","had"]
        
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
        couplings=np.logspace(-8,-3,20)
        
        
        
        
    elif model=="DarkHiggs":
        #decays=["e_e","mu_mu","pi_pi","4pi","K_K","other"]
        decays=["e_e","mu_mu","pi_pi","K_K"]
        
        
        masses = [   
            0.1   ,  
            0.2239,  
            0.5012,  
            0.6918,  
            0.8128,  
            0.955 ,
            1.5849,  
            3.5   ,
            11.22 ,  
            25.119
        ]

        couplings=np.logspace(-6,-3,20)
        
        #masses = [0.6918]
        #couplings = [0.001]
        
        # Use full mass-coupling fidelity (same as original FORESEE reach plots)
        if setup and setup == "full":
            masses = [
                0.1   ,  0.1122,  0.1259,  0.1413,  0.1585,  0.1778,  0.1995,  
                0.2239,  0.2512,  0.2818,  0.3162,  0.3548,  0.3981,  0.4467,  
                0.5012,  0.5623,  0.6026,  0.631 ,  0.6457,  0.6607,  0.6761,  
                0.6918,  0.7079,  0.7244,  0.7413,  0.7586,  0.7762,  0.7943,  
                0.8128,  0.8318,  0.8511,  0.871 ,  0.8913,  0.912 ,  0.9333,  
                0.955 ,  0.9772,  1.    ,  1.122 ,  1.2589,  1.4125,  1.5   ,
                1.5849,  1.7783,  1.9953,  2.2387,  2.5119,  2.8184,  3.1623,  
                3.5   ,  3.7   ,  3.9811,  5.0119,  6.3096,  7.9433,  10.   ,
                11.22 ,  12.589,  14.125,  15.849,  17.783,  19.953,  22.387,  
                25.119,  28.184,  31.623,  39.811,  50.119,  55.000,  60.000,
                63.096,  79.430,  99.9
            ]
            couplings=np.logspace(-8,-3,20)

    else:
        print(f"ERROR: Couldn't find setup for model {model}")

    return masses,couplings,decays



def clear_csvs(runmode,setup,currdir,outdir,energy):
    if runmode=="eff":
        if "G4" not in setup:
            return
        for G4setup in setup["G4"]:
            for station in setup["stations"]:
                for eff in setup["effs"]:
                    effname=f"{currdir}/{outdir}/eff_{energy}TeV_{G4setup}_station{station}_{eff}.csv"
                    print(f"Clearing eff file {effname}")
                    efffile = open(effname,'w')
                    efffile.close()


def get_effs(hepname,sumall,sumeffs,outroot,decay,station,effs):
          
    # get xs from hepmc
    if not os.path.exists(hepname):
        print(f"HepMC file {hepname} not found - skipping")
        return sumall,sumeffs
    
    xs=0.
    print("JOSH2",hepname)
    with pyhepmc.open(hepname) as f:
        print("JOSH3")
        event = f.read()
        print("JOSH4")
        xs=event.cross_section.xsec()
        print("Found cross section:",xs)
        print("JOSH5")
        if not os.path.exists(outroot):
            print(f"ROOT file {outroot} not found - skipping")
            return sumall,sumeffs

        print(f"Working on ROOT file {outroot}")
        
        try:
            f_eff=ROOT.TFile.Open(outroot,"READ")
        except:
            print(f"ERROR: ROOT file {outroot} not opened properly")
            return sumall,sumeffs
            
        t_eff = f_eff.Get(f"Hits{station}")
        if not t_eff:
            print(f"ERROR: TTree Hits{station} not opened properly")
            return sumall,sumeffs

        
        varp=''
        varm=''
        if decay=="e_e":
            varp='ep'
            varm='em'
        elif decay=="mu_mu":
            varp='mp'
            varm='mm'
        else:
            varp='hp'
            varm='hm'
            
        # Just use this to make a histogram, doesn't really matter what the variable is
        vary=varp+'_y'                                
        print(f"Using var {vary}")
        
        # Get total number of events
        hall=ROOT.TH1F("hall","hall",1000,0,10000)                                
        t_eff.Draw(f"{vary}>>hall");
        # Normalise to cross section
        sumall+=hall.Integral()*xs
        print(f"decay = {decay}, all =",hall.Integral(),xs,hall.Integral()*xs)
        
        # Get number of events passing each efficiency cut
        for n,(effval,efftitle,effstring) in enumerate(effs):
            # Correct eff string to be for relevant object in root file
            effstring=effstring.replace('ep',varp).replace('em',varm)
            heff=ROOT.TH1F("heff","heff",1000,0,10000)
            t_eff.Draw(f"{vary}>>heff",effstring)
            
            print(f"decay = {decay},",effstring,",",effs[n],"=",heff.Integral(),xs,heff.Integral()*xs)
            sumeffs[n]+=heff.Integral()*xs

        f_eff.Close()
        
    return sumall,sumeffs
