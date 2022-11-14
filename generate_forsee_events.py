import os,sys
from foresee import Foresee, Model, Utility
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

class ForeseeGenerator(object):
    """
    Generate LLP particles within FASER acceptance from FORESEE
    """
    
    def __init__(self, modelname, energy, mass, couplings, daughter1_pid, daughter2_pid, outdir = None, path = '.', randomSeed = 12345, t0 = 0, notime = False, suffix = "", selection="np.sqrt(x.x**2 + (x.y - 0.06)**2)< 0.1", length=1.5, distance=480.):

        self.modelname = modelname
        self.energy = energy
        self.mass = mass
        self.couplings = [couplings] if isinstance(couplings, (str, int, float)) else couplings
        self.daughter1_pid = daughter1_pid
        self.daughter2_pid = daughter2_pid
        self.outdir = outdir
        self.path = path
        self.version = 3  # Forsee "version"
        self.seed = randomSeed
        self.t0 = t0
        self.notime = notime
        self.suffix = f"_{suffix}" if suffix else ""
        self.nbinsample = 1
        self.selection = selection
        self.length = length
        self.distance = distance

        # Set decay mode ...
        
        self.pid_map = { 
            (-11, 11) : ["e_e"],
            (11, -11) : ["e_e"],                       
            (-13, 13) : ["mu_mu"],
            (13, -13) : ["mu_mu"],            
            (22, 22) : ["gamma_gamma"],
            (999,999) : ["K_K",
                         "K_K_pi",
                         "c_c",
                         "d_d",
                         "other_hadrons",
                         "pi+_pi-",
                         "pi+_pi-_pi+_pi-",
                         "pi+_pi-_pi0",
                         "pi+_pi-_pi0_pi0",
                         "pi0_gamma",
                         "s_s",
                         "u_u"]

            }

        self.mode = self.pid_map.get((self.daughter1_pid, self.daughter2_pid), None)
        if self.mode is None:
            sys.exit(f"Undefined decay to {self.daughter1_pid} + {self.daughter2_pid} for {self.modelname}")

        # Set detector ...
        if self.version == 1:
            self.foresee = Foresee()
        else:
            self.foresee = Foresee(path = self.path)

        # Generate 6 cm high to account for translation from ATLAS to FASER coord. system
        # TODO: relax this a bit as daughters may enter even if mother doesn't
        self.foresee.set_detector(selection=self.selection,
                                  channels=self.mode, distance=float(self.distance), length=float(self.length) ,
                                  luminosity=1/1000.) # 1 pb-1        

        # Set model ...
        if self.version == 1:
            self.model = Model(self.modelname)
        else:
            self.model = Model(self.modelname, path = f"{self.path}/Models/{self.modelname}/")
        
        if self.modelname == "DarkPhoton":
            self.data = self.darkphoton()
        elif self.modelname == "ALP-W":
            self.data = self.alp_W()
        else:
            sys.exit(f"Unknown model {self.modelname}")

        return

    def darkphoton(self):

        self.nbinsample = 100 # resample bins to help with asymmetric detector

        # Production modes
        self.model.add_production_2bodydecay(
            pid0 =  "111",
            pid1 = "22",
            br = "2.*0.99 * coupling**2 * pow(1.-pow(mass/self.masses('111'),2),3)",
            generator = "EPOSLHC",
            energy = self.energy,
            nsample = 10)
    
        self.model.add_production_2bodydecay(
            pid0 = "221",
            pid1 = "22",
            br = "2.*0.39 * coupling**2 * pow(1.-pow(mass/self.masses('221'),2),3)",
            generator = "EPOSLHC",
            energy = self.energy,
            nsample = 10)

        self.model.add_production_mixing(
            pid = "113",
            mixing = "coupling * 0.3/5. * 0.77545**2/abs(mass**2-0.77545**2+0.77545*0.147*1j)",
            generator = "EPOSLHC",
            energy = self.energy,
            )

        if self.version == 1:
            self.model.add_production_direct(
                label = "Brem",
                energy = self.energy,
                condition = "p.pt<1",
                coupling_ref=1,
                )
            
            self.model.add_production_direct(
                label = "DY",
                energy = self.energy,
                coupling_ref=1,
                massrange=[1.5, 10.]
                )

            self.model.set_br_1d(
                modes = self.mode,
                filenames=[f"files/models/{self.modelname}/br/{self.mode}.txt"] 
                )

        else:
            masses_brem = [ 
                0.01  ,  0.0126,  0.0158,  0.02  ,  0.0251,  0.0316,  0.0398,
                0.0501,  0.0631,  0.0794,  0.1   ,  0.1122,  0.1259,  0.1413,
                0.1585,  0.1778,  0.1995,  0.2239,  0.2512,  0.2818,  0.3162,
                0.3548,  0.3981,  0.4467,  0.5012,  0.5623,  0.6026,  0.631 ,
                0.6457,  0.6607,  0.6761,  0.6918,  0.7079,  0.7244,  0.7413,
                0.7586,  0.7762,  0.7943,  0.8128,  0.8318,  0.8511,  0.871 ,
                0.8913,  0.912 ,  0.9333,  0.955 ,  0.9772,  1.    ,  1.122 ,
                1.2589,  1.4125,  1.5849,  1.7783,  1.9953,  2.2387,  2.5119,
                2.8184,  3.1623,  3.9811,  5.0119,  6.3096,  7.9433, 10.    
                ]
            
            self.model.add_production_direct(
                label = "Brem",
                energy = self.energy,
                condition = "p.pt<1",
                coupling_ref=1,
                masses = masses_brem,
                )
            
            masses_dy = [ 
                1.5849,  1.7783,  1.9953,  2.2387,  2.5119, 2.8184,  3.1623,  3.9811,  5.0119,  6.3096,  7.9433, 10.    
                ]
            
            self.model.add_production_direct(
                label = "DY",
                energy = self.energy,
                coupling_ref=1,
                masses = masses_dy,
                )

            #self.model.set_br_1d(
            #    modes = ["e_e", "mu_mu"],
            #    finalstates = [[11, -11], [13, -13]],
            #    filenames=["model/br/e_e.txt", "model/br/mu_mu.txt"],
            #    #filenames=[f"model/br/all.txt"] 
            #    )

            allmodes=[
                "K_K",
                "K_K_pi",
                "c_c",
                "d_d",
                "e_e",
                "mu_mu",
                #"other_hadrons",
                "pi+_pi-",
                "pi+_pi-_pi+_pi-",
                "pi+_pi-_pi0",
                "pi+_pi-_pi0_pi0",
                "pi0_gamma",
                "s_s",
                "tau_tau",
                "u_u"
            ]

            allfinalstates=[
                [321, -321],
                [321, -321, 111],
                [4, -4],
                [1,-1],
                [11, -11],
                [13, -13],
#                "other_hadrons",
                [211, -211],
                [211, -211, 211, -211],
                [211, -211, 111],
                [211, -211, 111, 111],
                [111, 22],
                [3, -3],
                [15, -15],
                [2,-2]
            ]
            #finalstates = [[11, -11], [13, -13]],
            
            ## Branching ratio
            self.model.set_br_1d(
                #modes=["e_e", "mu_mu"],
                #filenames=["files/models/"+modelname+"/br/e_e.txt","files/models/"+modelname+"/br/mu_mu.txt"]
                modes=allmodes,
                finalstates = allfinalstates,
                filenames=["model/br/"+mode+".txt" for mode in allmodes]
            )

            
        return self.decays()


    def alp_W(self):

        self.nbinsample = 100 # resample bins to help with smoothness

        self.model.add_production_2bodydecay(
            pid0 = "5",
            pid1 = "321",
            br = "2.2e4 * coupling**2 * np.sqrt((1-(mass+0.495)**2/5.279**2)*(1-(mass-0.495)**2/5.279**2))",
            generator = "Pythia8",
            energy = self.energy,
            nsample = 500, # Vary over phi and theta -> increase number to improve smoothness for B
            )
            
        self.model.add_production_2bodydecay(
            pid0 = "-5",
            pid1 = "321",
            br = "2.2e4 * coupling**2 * np.sqrt((1-(mass+0.495)**2/5.279**2)*(1-(mass-0.495)**2/5.279**2))",
            generator = "Pythia8",
            energy = self.energy,
            nsample = 500,
            )
        
        self.model.add_production_2bodydecay(
            pid0 = "130",
            pid1 = "111",
            br = "4.5 * coupling**2 * np.sqrt((1-(mass+0.135)**2/0.495**2)*(1-(mass-0.135)**2/0.495**2))",
            generator = "EPOSLHC",
            energy = self.energy,
            nsample = 50,
            )

        self.model.add_production_2bodydecay(
            pid0 = "321",
            pid1 = "211",
            br = "10.5 * coupling**2 * np.sqrt((1-(mass+0.135)**2/0.495**2)*(1-(mass-0.135)**2/0.495**2))",
            generator = "EPOSLHC",
            energy = self.energy,
            nsample = 50,
            )

        self.model.add_production_2bodydecay(
            pid0 = "-321",
            pid1 = "211",
            br = "10.5 * coupling**2 * np.sqrt((1-(mass+0.135)**2/0.495**2)*(1-(mass-0.135)**2/0.495**2))",
            generator = "EPOSLHC",
            energy = self.energy,
            nsample = 50,
            )


        if self.version == 1:
            self.model.set_br_1d(
                modes = self.mode,
                filenames=[f"files/models/{self.modelname}/br/{self.mode}.txt"] 
                )
        else:
            self.model.set_br_1d(
                modes = ["gamma_gamma"],
                finalstates = [[22, 22]],
                filenames=["model/br/gamma_gamma.txt"]
                #filenames=[f"model/br/all.txt"] 
                )   

        return self.decays()


    def decays(self):
        # Set up liftime and BRs

        if self.version == 1:
            self.model.set_ctau_1d(
                filename=f"files/models/{self.modelname}/ctau.txt", 
                coupling_ref=1
                )        
        else:
            self.model.set_ctau_1d(
                filename=f"model/ctau.txt", 
                coupling_ref=1
                )
            
        # Get LLP spectrum
        self.foresee.set_model(model=self.model)
        # This is just a reference coupling 
        plt = self.foresee.get_llp_spectrum(self.mass, coupling=1, do_plot=True)  
        plt.savefig(f"{self.modelname}_m{self.mass}.png")
        plt.close()

        def flatten(l):
            return [i for sublist in l for i in sublist]

        # Get list of events within detector
        output = self.foresee.get_events(mass=self.mass, energy=self.energy, couplings=self.couplings)        

        if self.version >= 3:
            coups, ctaus, nsigs, momenta, weights = output
            fmomenta  = flatten(momenta)
            fweights  = flatten(weights)            
            fenergies = [p.e for p in fmomenta] # Now in MeV ?
            fthetas   = [p.pt/p.pz for p in fmomenta]            

            self.plot(fthetas, fenergies, fweights)

            # Return energy, theta and weights
            return [fenergies, fthetas, fweights]             
        else:
            coups, ctaus, nsigs, energies, weights, thetas = output
            self.plot(flatten(thetas), flatten(energies), flatten(weights))

            # Return energy (converting to MeV), theta and weights
            return [[e*1000 for e in flatten(energies)], flatten(thetas), flatten(weights)] 

    def plot(self, thetas, energies, weights):
        # Plot the results in Forsee format

        #print("Nwgts:",len(weights))
        #print("min wgts:",min(weights))
        #print("max wgts:",max(weights))
        if max(weights) == 0:
            print("Max wgts = 0 - skipping plots")
            return
        
        t = np.array(thetas)
        p = np.sqrt(np.array(energies)**2 - self.mass**2)

        prange=[[-6, 0, 120],[ 0, 5, 50]]
        tmin, tmax, tnum = prange[0]
        pmin, pmax, pnum = prange[1]
        t_edges = np.logspace(tmin, tmax, num=tnum+1)
        p_edges = np.logspace(pmin, pmax, num=pnum+1)  

        ticks = np.array([[np.linspace(10**(j),10**(j+1),9)] for j in range(-7,6)]).flatten()
        ticks = [np.log10(x) for x in ticks]
        ticklabels = np.array([[r"$10^{"+str(j)+"}$","","","","","","","",""] for j in range(-7,6)]).flatten()
        matplotlib.rcParams.update({'font.size': 15})

        fig = plt.figure(figsize=(8,5.5))
        ax = plt.subplot(1,1,1)
        h=ax.hist2d(x=np.log10(t),y=np.log10(p),weights=weights,
                    bins=[tnum,pnum],range=[[tmin,tmax],[pmin,pmax]],
                    norm=matplotlib.colors.LogNorm(), cmap="hsv",
                    )
        
        #fig.colorbar(h[3], ax=ax)
        ax.set_xlabel(r"angle wrt. beam axis $\theta$ [rad]")
        ax.set_ylabel(r"momentum $p$ [GeV]")
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabels)
        ax.set_xlim(tmin, tmax)
        ax.set_ylim(pmin, pmax)
        plt.savefig(f"{self.modelname}_m{self.mass}_acc.png")
        plt.close()

    def write(self):
        # Write LLP results to a file
        
        energies, thetas, weights = self.data

        if self.outdir is None:
            if self.version == 1:
                self.outdir = f"files/models/{self.modelname}/events"
            else:
                self.outdir = f"{self.foresee.dirpath}/Models/{self.modelname}/model/events"

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        if len(self.couplings) == 1:
            filename = f"{self.outdir}/events_{self.energy}TeV_m{self.mass}GeV_c{self.couplings[0]}to_{self.daughter1_pid}_{self.daughter2_pid}.npy"
        else:
            filename = f"{self.outdir}/events_{self.energy}TeV_m{self.mass}GeV_to_{self.daughter1_pid}_{self.daughter2_pid}.npy"

        print(f"Generated {len(thetas)} events")
        print(f"save data to file: {filename}")
        np.save(filename,[energies,thetas, weights])

        cfgname = filename.replace(".npy", ".cfg")
        print(f"save config to file: {cfgname}")
        with open(cfgname, "w") as f:
            f.write(" ".join(sys.argv))

        return filename

    def write_hepmc(self, nevents):

        if self.outdir is None:
            self.outdir = "model/events/"
        elif not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

            
        #filename =  f"{self.outdir}/events_{self.energy}TeV_m{self.mass}GeV_c{self.couplings[0]}to_{self.daughter1_pid}_{self.daughter2_pid}{self.suffix}.hepmc"
        filename =  f"events_{self.energy}TeV_m{self.mass}GeV_c{self.couplings[0]}to_{self.daughter1_pid}_{self.daughter2_pid}{self.suffix}.hepmc"
        print("JOSH0")
        self.foresee.write_events(self.mass, self.couplings[0], self.energy, filename=filename, numberevent=nevents, zfront = -1.5, seed = self.seed, decaychannels = self.mode, notime = self.notime, t0 = self.t0, nsample = self.nbinsample, outdir=self.outdir)

        #cfgname = f"{self.foresee.dirpath}/Models/{self.modelname}/" + filename.replace(".hepmc", ".cfg")
        cfgname = self.outdir+"/"+filename.replace(".hepmc", ".cfg")
        print(f"save config to file: {cfgname}")
        with open(cfgname, "w") as f:
            f.write(" ".join(sys.argv))

        return filename


def setup_foresee(path):

    if path is None:
        return

    # Add foresee to python path
    path = os.path.expandvars(os.path.expanduser(path))
    os.sys.path.append(f"{path}/FORESEE/src")

    # Symlink foresee files/Models dirs to current dir
    #if not os.path.exists("files"):
    #    os.symlink(os.path.expandvars(f"{path}/FORESEE/files"), "files")
    #if not os.path.exists("Models"):
    #    os.symlink(os.path.expandvars(f"{path}/FORESEE/Models"), "files")        

    # Install scikit-hep if needed.

    try:
        from skhep.math.vectors import LorentzVector, Vector3D
    except ModuleNotFoundError:
        os.system("pip install scikit-hep --user")
        try:
            from skhep.math.vectors import LorentzVector, Vector3D
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Unable to find skhep.  Please install the scikit-hep package")
        
    return

def add_to_python_path(path):
    if path in sys.path: return
    path = os.path.expandvars(os.path.expanduser(path))
    os.sys.path.append(path)
    return

def parse_couplings(data, write_hepMC = False):    

    if write_hepMC:
        if len(data) == 1:
            couplings = float(data[0])
        else:
            sys.exit("Only a single coupling allowed when writing HEPMC events")

    try:
        couplings = [float(d) for d in data]
    except ValueError:
        try:
            couplings = np.logspace(*eval(data[0]))
        except:
            sys.exit("Unable to parse couplings")

    return couplings

def main():

    import argparse
    
    parser = argparse.ArgumentParser(description="Run FORSEE generation")
    parser.add_argument("model", help = "Name of foresee model")
    parser.add_argument("--mass", "-m", required = True, type = float, help = "Mass of mother [GeV]")
    parser.add_argument("--couplings", "-c", required = True, nargs = "+", help = "Couplings of mother (either single/mulitple values or tuple to pass to np.logspace)")
    parser.add_argument("--pid1", required = True, type = int, help = "PID of daughter 1")
    parser.add_argument("--pid2", default = None, type = int, help = "PID of daughter 2 (if not set then will be -PID1)")
    parser.add_argument("--Ecom", default = "14", help = "Center of mass energy [TeV]")
    parser.add_argument("--outdir", "-o", default = None, help = "Output path")    
    parser.add_argument("--path", default = ".", help = "Path to foresee installation")
    parser.add_argument("--hepmc", action = "store_true", help = "Write HepMC events")
    parser.add_argument("--nevents", "-n", default = 10, type = int, help = "Number of HepMC events ")
    parser.add_argument("--randomSeed", "-s", default = 1234, type = int, help = "Random seed for HepMC generation")
    parser.add_argument("--t0", "-t", default = 0, type = int, help = "Time offset for start of decay volume")
    parser.add_argument("--notime", action = "store_true", help = "Set all vertex times to 0 rather than calculating from start of decay volume")
    parser.add_argument("--suffix", default = "", help = "Filename suffix")
    parser.add_argument("--selection", default = "np.sqrt(x.x**2 + (x.y - 0.06)**2)< 0.1", help = "Detector transverse selection")
    parser.add_argument("--length", default = 1.5, help = "Detector length")
    parser.add_argument("--distance", default = 480., help = "Detector distance from IP")

    args = parser.parse_args()

    add_to_python_path(f"{args.path}/src")


    # Create PIDs
    if args.pid2 is None:
        args.pid2 = -args.pid1
    
    couplings = parse_couplings(args.couplings, args.hepmc)

    print(f"Generating {args.model} events at Ecom = {args.Ecom}") 
    print(f"   mother mass = {args.mass} GeV")
    print(f"   decay = {args.pid1} {args.pid2}")
    print(f"   couplings = {couplings}")    

    f = ForeseeGenerator(args.model, args.Ecom, args.mass, couplings, args.pid1, args.pid2, outdir = args.outdir, path = args.path, randomSeed = args.randomSeed, t0 = args.t0, notime = args.notime, suffix = args.suffix, selection=args.selection, length=args.length, distance=args.distance)

    #print("args: model, Ecom, mass, couplings, pid1, pid2, outdir, path, randomSeed, t0, notime, suffix")
    #print("args:",args.model, args.Ecom, args.mass, couplings, args.pid1, args.pid2, args.outdir, args.path, args.randomSeed, args.t0, args.notime, args.suffix)

    if args.hepmc:
        f.write_hepmc(args.nevents)
    else:
        f.write()
    
        

if __name__ == "__main__":
    main()
