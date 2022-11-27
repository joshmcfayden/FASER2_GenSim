from generate_forsee_events import ForeseeGenerator
import numpy as np
import os
import gc # import garbage collector interface
import shutil
import subprocess
from array import array
import utils,setups



#args: model, Ecom, mass, couplings, pid1, pid2, outdir, path, randomSeed, t0, notime, suffix
#args: DarkPhoton 14 0.1 [0.0001] 11 -11 mytest /Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/FORESEE/ 0 0 False test1

runmode="run"
runmode="combine"
#runmode="G4"
#runmode="eff"
#runmode="plotsep"

do_hepmc=False
do_hepmc=True

# Set model inputs
model="DarkPhoton"
model="DarkHiggs"
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
path="/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/FORESEE/"
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
        
    m_c_nevents=[]

    nmasstot=len(masses)

    for nmass,mass in enumerate(masses):

        c_nevents=[]
        ncouptot=len(couplings)
        
        for ncoup,coup in enumerate(couplings):
    
            nsignal=0.

            # Read G4 output files to get efficiencies for each configuration
            if runmode=="eff":

                if "G4" not in setup:
                    print(f"Error: No G4 configuration found for {setup_name}")
                    continue
                

                for G4setup in setup["G4"]:

                    for station in setup["stations"]:
                    
                        # Calculate efficiencies 
                        for eff in setup["effs"]:
    
                            sumall=0. # denominator for eff calc
                            sumeffs=[0. for x in setup["effs"][eff]] # numerator(s) for eff calc
                            
                            for decay in decays:
                                #hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}_{suffix}.hepmc"
                                hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}_{suffix}.hepmc"
                                outroot=f"{currdir}/{hepname}".replace('.hepmc',f'_{G4setup}.root')
    
                                sumall,sumeffs=utils.get_effs(hepname,sumall,sumeffs,outroot,decay,station,setup["effs"][eff])
                                print("JOSH1:",sumall,sumeffs)
                                
                            # open output file to append results to
                            effname=f"{currdir}/{outdir}/eff_{energy}TeV_{G4setup}_station{station}_{eff}.csv"
                            print(f"Opening eff file {effname}")
                            efffile = open(effname,'a')
    
                            for n,(effval,efftitle,effstring) in enumerate(setup["effs"][eff]):
                                effcalc=sumeffs[n]/sumall if sumall > 0 else 1.0
                                print(f'{mass},{coup},{effval},{effcalc}\n')
                                efffile.write(f'{mass},{coup},{effval},{effcalc}\n')
                                
                            #efffile.close()


            ndectot=len(decays)
            
            for ndec,decay in enumerate(decays):
                
                #print(f"Generating {model} events at Ecom = {energy}") 
                #print(f"   mother mass = {mass} GeV")
                #print(f"   decay = {pid1} {pid2}")
                #print(f"   couplings = {coup}")    

                print(f"\nGenerating {model} events at Ecom = {energy}, mass = {mass} GeV ({nmass+1}/{nmasstot}), coupling = {coup} ({ncoup+1}/{ncouptot}), decay = {decay} ({ndec+1}/{ndectot}) [{(nmass*ncouptot*ndectot)+(ncoup*ndectot)+(ndec)+1}/{nmasstot*ncouptot*ndectot}]")
                
                npname=""
                
                npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}.npy"
                hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}_{suffix}.hepmc"
                
                #break

            

                if "run" in runmode:
                    f = None
                    #npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}.npy"
                    #npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}.npy"
                    if os.path.exists(npname):
                        #if do_hepmc and os.path.exists(hepname):
                        #    print(f"Files {npname} and {hepname} found, skipping to next mass/coup/decay point")
                        #    continue
                        #else:
                        print(f"File {npname} found, moving to next step...")
                    #f = ForeseeGenerator(model, energy, mass, [coup], pid1, pid2, outdir = outdir, path = path, randomSeed = randomSeed, t0 = t0, notime = notime, suffix = suffix, selection=setup["selection"], length=setup["length"], distance=setup["distance"])
                    else:

                        print(f"   Working on {npname}")
                        f = ForeseeGenerator(model, energy, mass, [coup], decay, outdir = outdir, path = path, randomSeed = randomSeed, t0 = t0, notime = notime, suffix = suffix, selection=setup["selection"], length=setup["length"], distance=setup["distance"])
                        f.write()
                       
                    # Run HepMC creation
                    if do_hepmc:
                        if os.path.exists(hepname):
                            print(f"File {hepname} found, skipping to next")
                            #continue
                        else:

                            print(f"   Working on {hepname}")

                            skip=False
                            if os.path.exists(npname):
                                np_arr=np.load(npname)
                                energies=np_arr[0]
                                thetas=np_arr[1]
                                weights=np_arr[2]
                                if max(weights) == 0.:
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
                            


                if "combine" in runmode:
                    if "run" not in runmode:
                        #npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}.npy"
                        #npname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}.npy"

                        print("   Using file:",npname)


                    if not os.path.exists(npname):
                        print(f"WARNING: File {npname} not found, skipping to next")
                        continue
                    
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
                    
                #hepname=f"{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}_{suffix}.hepmc"




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

                if runmode=="plotsep":
                    import ROOT
                
                    for G4setup in setup["G4"]:

                        ROOT.gStyle.SetOptStat(0)
                        ROOT.gStyle.SetOptTitle(0)
                        
                        
                        #cout << "../"+outdir+"/events_"+energy+"TeV_m"+mass+"GeV_c"+coupling+"to_"+pid1+"_"+pid2+"_s1_"+G4setup+".root" << endl

                        infile=f"{currdir}/{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}_s1_{G4setup}.root"
                        if not os.path.exists(infile):
                            print(f"ROOT file {infile} not found - skipping")
                            continue
                        
                        f_L1 = ROOT.TFile.Open(infile)
                        t_L1 = f_L1.Get("Hits1")
                        
                        f_L2 = ROOT.TFile.Open(infile)
                        t_L2 = f_L2.Get("Hits2")
                        
                        f_L3 = ROOT.TFile.Open(infile)
                        t_L3 = f_L3.Get("Hits3")
       
                        nbins = 50
                        xmin = 1e-3
                        xmax = 1e5
                        logxmin = ROOT.TMath.Log10(xmin)
                        logxmax = ROOT.TMath.Log10(xmax)
                        binwidth = (logxmax-logxmin)/nbins
                        xbins = [xmin]
                        for i in range(1,nbins+1) :
                            xbins.append(float(xmin + ROOT.TMath.Power(10,logxmin+i*binwidth)))
                            
  
                        c1 = ROOT.TCanvas("c1","c1")
                        
                        #print(xmin,array('d',xbins))
                        
                        h_dy_L1 = ROOT.TH1D("h_dy_L1","h_dy_L1",nbins,array('d',xbins))
                        print(t_L1,h_dy_L1)
                        t_L1.Draw("abs(ep_y-em_y)>>h_dy_L1")
                        print(h_dy_L1)
                        
                        h_dy_L2 = ROOT.TH1D("h_dy_L2","h_dy_L2",nbins,array('d',xbins))
                        t_L2.Draw("abs(ep_x-em_x)>>h_dy_L2")
                        
                        h_dy_L3 = ROOT.TH1D("h_dy_L3","h_dy_L3",nbins,array('d',xbins))
                        t_L3.Draw("abs(ep_x-em_x)>>h_dy_L3")
                        
                        if h_dy_L1.Integral(): h_dy_L1.Scale(1./h_dy_L1.Integral())
                        if h_dy_L2.Integral(): h_dy_L2.Scale(1./h_dy_L2.Integral())
                        if h_dy_L3.Integral(): h_dy_L3.Scale(1./h_dy_L3.Integral())

  
                        h_dy_L1.SetMaximum(1.1*ROOT.TMath.Max(h_dy_L1.GetMaximum(),h_dy_L3.GetMaximum()))
                        h_dy_L1.GetXaxis().SetTitle("Separation [mm]")
                        h_dy_L1.GetXaxis().SetTitleSize(0.045)
                        
                        h_dy_L1.SetLineColor(ROOT.kBlue+1)
                        h_dy_L1.SetLineWidth(3)
                        h_dy_L1.Draw("hist")
                        
                        h_dy_L2.SetLineStyle(ROOT.kDashed)
                        h_dy_L2.SetLineColor(ROOT.kBlue+1)
                        h_dy_L2.SetLineWidth(3)
                        h_dy_L2.Draw("histsame")

                        h_dy_L3.SetLineStyle(ROOT.kDotted)
                        h_dy_L3.SetLineColor(ROOT.kBlue+1)
                        h_dy_L3.SetLineWidth(3)
                        h_dy_L3.Draw("histsame")

                        
                        c1.SetLogx()
                        c1.SetTickx()
                        c1.SetTicky()
                        
                        leg = ROOT.TLegend(0.7,0.85,0.85,0.6)
                        leg.SetBorderSize(0)
                        leg.SetFillColor(0)
                        leg.SetTextSize(0.04)
                        leg.AddEntry(h_dy_L1,"Station 1","l")
                        leg.AddEntry(h_dy_L2,"Station 2","l")
                        leg.AddEntry(h_dy_L3,"Station 3","l")
      
                        leg.Draw()
                        
                        latex = ROOT.TLatex()
                        latex.SetNDC()
                        latex.SetTextFont(42)
                        latex.SetTextSize(0.06)
                        latex.DrawLatex(0.15,0.8,"#bf{#it{FASER2}}")
                        latex.SetTextSize(0.05)
                        latex.DrawLatex(0.32,0.8,setup_name)
                        latex.SetTextSize(0.04)

                        latex2 = ROOT.TLatex()
                        latex2.SetNDC()
                        latex2.SetTextFont(42)
                        latex2.SetTextSize(0.04)
                        latex2.DrawLatex(0.15,0.75,f"m={mass} GeV")
                        latex2.DrawLatex(0.15,0.7,f"#varepsilon={coup:.5g}")
                        latex2.SetTextSize(0.04)

                        c1.SaveAs(f"{currdir}/{outdir}/plot_sep_stations_{energy}TeV_m{mass}GeV_c{coup}to_{pid1}_{pid2}_s1_{G4setup}.pdf")

  




            c_nevents.append(nsignal)
                                  


        m_c_nevents.append(c_nevents)
    

    if "combine" in runmode:
        outdir="FORESEE/Models/"+model+"/model/results/"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        outfile=outdir+"/"+energy+"TeV_"+setup_name+".npy"
        
        print("\nWriting output file:",outfile)
        #np.save(outfile,[masses,couplings,m_c_nevents])
        np.save(outfile,np.array([masses,couplings,m_c_nevents],dtype=object))    

    
    # Go back to starting directory 
    os.chdir(currdir)
