import numpy as np
from foresee import Foresee, Utility, Model
#import foresee
import os


scan_name="pick"
scan_search={"pick":["F2-default","S1-L1p5-D2","S2-L2-D2","S3-L10-D1","S3-L10-D2"]}
scan_search={"pick":["F2-default","S1-L1p5-D1","S1-L1p5-D2","S2-L10-D1","S2-L10-D2"]}
scan_search=None

setup_dict={
    "R1-L10-R1x3":{
        #"name":"R1 L=10m X=3m Y=1m",
        "name":"New Baseline (X=3m Y=1m)",
        #"color":"cornflowerblue",
        "color":"forestgreen",
        #"color":"limegreen",
        "style":"dashed",
        "selection":"(np.sqrt(x.x**2)<1.5) * (np.sqrt(x.y**2)<0.5)",
        "length":10,
        "distance":615,
        "channels": None,
#        "G4":["FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_AllTrkStations"],
#        "effs":{
#            "sep":[(0.1,'sep>0.1mm','abs(ep_x-em_x)>0.1'),
#                   (1,'sep>1mm','abs(ep_x-em_x)>1'),
#                   (5,'sep>5mm','abs(ep_x-em_x)>5'),
#                   (10,'sep>10mm','abs(ep_x-em_x)>10'),
#                   (100,'sep>100mm','abs(ep_x-em_x)>10')],
#        }
    },
}

print("INFO   : Initialise FORESEE")
foresee = Foresee()


print("INFO   : Setting up Dark Higgs model")
energy = "14"
modelname = "DarkHiggs"
model = Model(modelname,path="./FORESEE/Models/DarkHiggs/")


print("INFO   :   - Adding production modes")

## 2-body decays
model.add_production_2bodydecay(
    pid0 = "5",
    pid1 = "321",
    br = "5.7 * coupling**2 * pow(1.-pow(mass/5.279,2),2)",
    generator = "Pythia8",
    energy = energy,
    nsample = 10,
)

model.add_production_2bodydecay(
    pid0 = "-5",
    pid1 = "321",
    br = "5.7 * coupling**2 * pow(1.-pow(mass/5.279,2),2)",
    generator = "Pythia8",
    energy = energy,
    nsample = 10,
)

#    model.add_production_2bodydecay(
#        pid0 = "25",
#        pid1 = "0",
#        br = "2*0.05",
#        generator = "Pythia8",
#        energy = energy,
#        nsample = 100,
#        scaling = 0,
#    )
#
#    # 3-body
#
#    model.add_production_3bodydecay(
#        label= "5_di",
#        pid0 = "5",
#        pid1 = "321",
#        pid2 = "0",
#        br = "7.37e-10*np.sqrt(1-4*mass**2/q**2)*(1-q**2/4.5**2)**2",
#        generator = "Pythia8",
#        energy = energy,
#        nsample = 10,
#        scaling = 0, 
#    )
#    
#    model.add_production_3bodydecay(
#        label= "-5_di",
#        pid0 = "-5",
#        pid1 = "321",
#        pid2 = "0",
#        br = "7.37e-10*np.sqrt(1-4*mass**2/q**2)*(1-q**2/4.5**2)**2",
#        generator = "Pythia8",
#        energy = energy,
#        nsample = 10,
#        scaling = 0, 
#    )


print("INFO   :   - Setting lifetimes")
## Lifetime
model.set_ctau_1d(
    filename="model/ctau.txt", 
    coupling_ref=1
)

print("INFO   :   - Setting branching fractions")
## Branching ratio
#model.set_br_1d(
#    modes=["e_e", "mu_mu"],
#    filenames=["files/models/"+modelname+"/br/e_e.txt","files/models/"+modelname+"/br/mu_mu.txt"]
#)

allmodes=[
    "4pi",
    #        "K+_K-",
    "K_K",
    #        "Kl_Kl",
    #        "Ks_Ks",
    "b_b",
    "c_c",
    "e_e",
    "g_g",
    "mu_mu",
    #        "pi+_pi-",
    #        "pi+_pi-_pi+_pi-",
    #        "pi+_pi-_pi0_pi0",
    #        "pi0_pi0",
    #        "pi0_pi0_pi0_pi0",
    "pi_pi",
    "s_s",
    "tau_tau"
]

allmodes=[
    "e_e",
    "mu_mu",
    "pi_pi",
    "4pi",
    #        "K+_K-",
    "K_K",
    #        "Kl_Kl",
    #        "Ks_Ks",
    "b_b",
    "c_c",
    "g_g",
    #        "pi+_pi-",
    #        "pi+_pi-_pi+_pi-",
    #        "pi+_pi-_pi0_pi0",
    #        "pi0_pi0",
    #        "pi0_pi0_pi0_pi0",
    "s_s",
    "tau_tau"
]

## Branching ratio
model.set_br_1d(
    #modes=["e_e", "mu_mu"],
    #filenames=["files/models/"+modelname+"/br/e_e.txt","files/models/"+modelname+"/br/mu_mu.txt"]
    modes=allmodes,
    filenames=["model/br/"+mode+".txt" for mode in allmodes]
    
)


## Set model just created
foresee.set_model(model=model)


print("INFO   :   - Plotting reach")

#Now let's plot the results. We first specify all detector setups for which we want to show result (filename in model/results directory, label, color, linestyle, opacity alpha for filled contours, required number of events).

setups=[]
for setup in setup_dict:
    if not scan_search or setup in scan_search[scan_name]:
        setups.append(["14TeV_%s.npy"%setup, setup_dict[setup]["name"],setup_dict[setup]["color"], "solid", 0., 3, None])
        
        
print(f"INFO   :     - Found {len(setups)} setups")
#Then we specify all the existing bounds (filename in model/bounds directory, label, label position x, label position y, label rotation)

bounds = [ 
    ["bounds_1508.04094.txt", "LHCb $B^0$"  , 0.430, 2.2*10**-3, 90 ],
    ["bounds_1612.08718.txt", "LHCb $B^+$"  , 0.330, 2.2*10**-3, 90 ],
    ["bounds_1612.08718.txt", "LHCb $B^+$"  , 2.500, 2.2*10**-3, 90 ],
    ["bounds_LSND.txt"      , "LSND"        , 0.250, 9.0*10**-5, 90 ],
    ["bounds_CHARM.txt"     , "CHARM"       , 0.250, 4.0*10**-4, 90 ],
    ["bounds_MicroBoone.txt", "$\mu$BooNE"  , 0.138, 2.6*10**-4, 90 ],
    ["bounds_E949.txt"      , "E949"        , 0.102, 1.5*10**-4, 90 ],
    ["bounds_2011.11329.txt", "NA62 $K^+$"  , 0.170, 6.2*10**-4, 90 ],
    ["bounds_2010.07644.txt", "NA62 $\pi^+$", 0.125, 2.4*10**-3, 90 ],
]



#We then specify other projected sensitivitities (filename in model/bounds directory, color, label, label position x, label position y, label rotation)
    
projections = []
#projections = [
#    ["limits_SHiP.txt",       "teal",         "SHiP"    , 2.700, 3.2*10**-5, 0  ],
#    ["limits_MATHUSLA.txt",   "dodgerblue",   "MATHUSLA", 0.120, 5.0*10**-6, 0  ],
#    ["limits_CodexB.txt",     "deepskyblue",  "CodexB"  , 1.700, 2.0*10**-5, 0  ],
#    ["limits_LHCb.txt",       "cyan",         "LHCb"    , 3.800, 1.0*10**-4, 0  ],
#]

# Finally, we can plot everything using foresee.plot_reach(). It returns a matplotlib instance, to which we can add further lines and which we can show or save. Below, we add the dark matter relict target line for a specific benchmark.
plot = foresee.plot_reach(
    setups=setups,
    bounds=bounds,
    projections=projections,
    title="Dark Higgs", 
    xlims=[0.1,10], 
    ylims=[10**-5.5,10**-2.0],
    xlabel=r"Dark Higgs Mass $m_{\phi}$ [GeV]", 
    ylabel=r"Mixing $\theta$",
    legendloc=(0.43,0.25),
    figsize=(8,6),
)

plot.subplots_adjust(left=0.12, right=0.97, bottom=0.10, top=0.95)
plot.savefig("DarkHiggs_Pythia-Reach%s_TEST.pdf"%("_"+scan_name if scan_search else ""))
