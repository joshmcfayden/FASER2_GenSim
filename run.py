from generate_forsee_events import ForeseeGenerator
import numpy as np
import os
import gc # import garbage collector interface
import shutil
import subprocess
from array import array




#args: model, Ecom, mass, couplings, pid1, pid2, outdir, path, randomSeed, t0, notime, suffix
#args: DarkPhoton 14 0.1 [0.0001] 11 -11 mytest /Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/FORESEE/ 0 0 False test1

mode="run"
mode="G4"
#mode="runcombine"
#mode="combine"
#mode="eff"
mode="plotsep"

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

#    "F2-default":{
#        "name":"FASER2 Orig",# (Default)",
#        #"color":"maroon",
#        "color":"firebrick",
#        "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#        "length":5,
#        "distance":480,
#        "channels": None,
#        "G4":["FASER2_HepMC_v4_FASER2_Default_3rdTrkStation"],
#        "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_y-em_y)>0.1'),
#                       (1,'sep>1mm','abs(ep_y-em_y)>1'),
#                       (5,'sep>5mm','abs(ep_y-em_y)>5'),
#                       (10,'sep>10mm','abs(ep_y-em_y)>10'),
#                      (100,'sep>100mm','abs(ep_y-em_y)>100')]}
#    },
#
#    "S2-L10-D2":{
#        #"name":"S2 L=10m D=2m",
#        "name":"Old Baseline",
#        "color":"royalblue",
#        #"color":"darkgreen",
#        #"color":"orange",
#        "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#        "length":10,
#        "distance":615,
#        "channels": None,
#        "G4":["FASER2_HepMC_v4_FASER2_Cavern_3rdTrkStation"],
#        "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_y-em_y)>0.1'),
#                       (1,'sep>1mm','abs(ep_y-em_y)>1'),
#                       (5,'sep>5mm','abs(ep_y-em_y)>5'),
#                       (10,'sep>10mm','abs(ep_y-em_y)>10'),
#                       (100,'sep>100mm','abs(ep_y-em_y)>100')]}
#
#    },
#
    
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
             "G4":["FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_AllTrkStations"],
            "effs":{
                "sep":[(0.1,'sep>0.1mm','abs(ep_x-em_x)>0.1'),
                       (1,'sep>1mm','abs(ep_x-em_x)>1'),
                       (5,'sep>5mm','abs(ep_x-em_x)>5'),
                       (10,'sep>10mm','abs(ep_x-em_x)>10'),
                       (100,'sep>100mm','abs(ep_x-em_x)>10')],
#                "env":[(-250,'env<-250mm','abs(ep_x)<(1500-250)&&abs(em_x)<(1500-250)&&abs(ep_y)<(500-250)&&abs(em_y)<(500-250)'),
#                       (-100,'env<-100mm','abs(ep_x)<(1500-100)&&abs(em_x)<(1500-100)&&abs(ep_y)<(500-100)&&abs(em_y)<(500-100)'),
#                       (0,'env<0mm','abs(ep_x)<(1500)&&abs(em_x)<(1500)&&abs(ep_y)<(500)&&abs(em_y)<(500)'),
#                       (100,'env<100mm','abs(ep_x)<(1500+100)&&abs(em_x)<(1500+100)&&abs(ep_y)<(500+100)&&abs(em_y)<(500+100)'),
#                       (250,'env<250mm','abs(ep_x)<(1500+250)&&abs(em_x)<(1500+250)&&abs(ep_y)<(500+250)&&abs(em_y)<(500+250)'),
#                       (500,'env<500mm','abs(ep_x)<(1500+500)&&abs(em_x)<(1500+500)&&abs(ep_y)<(500+500)&&abs(em_y)<(500+500)'),
#                       (1000,'env<1000mm','abs(ep_x)<(1500+1000)&&abs(em_x)<(1500+1000)&&abs(ep_y)<(500+1000)&&abs(em_y)<(500+1000)')],
#                "envx":[(-1000,'env<-1000mm','abs(ep_x)<(1500-1000)&&abs(em_x)<(1500-1000)'),
#                        (-500,'env<-500mm','abs(ep_x)<(1500-500)&&abs(em_x)<(1500-500)'),
#                        (-250,'env<-250mm','abs(ep_x)<(1500-250)&&abs(em_x)<(1500-250)'),
#                        (-100,'env<-100mm','abs(ep_x)<(1500-100)&&abs(em_x)<(1500-100)'),
#                        (0,'env<0mm','abs(ep_x)<(1500)&&abs(em_x)<(1500)'),
#                        (100,'env<100mm','abs(ep_x)<(1500+100)&&abs(em_x)<(1500+100)'),
#                        (250,'env<250mm','abs(ep_x)<(1500+250)&&abs(em_x)<(1500+250)'),
#                        (500,'env<500mm','abs(ep_x)<(1500+500)&&abs(em_x)<(1500+500)'),
#                        (1000,'env<1000mm','abs(ep_x)<(1500+1000)&&abs(em_x)<(1500+1000)')],

            }
        },

#      "R1-L10-R0p5x3":{
#          #"name":"R1 L=10m X=3m Y=0.5m",
#          "name":"New Baseline (X=3m Y=0.5m)",
#          #"color":"lightsteelblue",
#          "color":"orchid",
#          #"color":"limegreen",
#          "style":"dotted",
#          "selection":"(np.sqrt(x.x**2)<1.5) * (np.sqrt(x.y**2)<0.25)",
#          "length":10,
#          "distance":615,
#          "channels": None,
#          "G4":["FASER2_HepMC_v4_FASER2_Cavern_Rect_KEKRect_3rdTrkStation"],
#          "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_x-em_x)>0.1'),
#                         (1,'sep>1mm','abs(ep_x-em_x)>1'),
#                         (5,'sep>5mm','abs(ep_x-em_x)>5'),
#                         (10,'sep>10mm','abs(ep_x-em_x)>10'),
#                         (100,'sep>100mm','abs(ep_x-em_x)>100')]}
#      },

#    "R1-L10-R0p5x2":{
#        #"name":"R1 L=10m X=2m Y=0.5m",
#        "name":"New Baseline (X=2m Y=0.5m)",
#        #"color":"lightsteelblue",
#        "color":"pink",
#        "style":"dotted",
#        "selection":"(np.sqrt(x.x**2)<1.0) * (np.sqrt(x.y**2)<0.25)",
#        "length":10,
#        "distance":615,
#        "channels": None,
#        "G4":["FASER2_HepMC_v4_FASER2_Cavern_Rect_KEKCircle_3rdTrkStation"],
#        "effs":{"sep":[(0.1,'sep>0.1mm','abs(ep_x-em_x)>0.1'),
#                       (1,'sep>1mm','abs(ep_x-em_x)>1'),
#                       (5,'sep>5mm','abs(ep_x-em_x)>5'),
#                       (10,'sep>10mm','abs(ep_x-em_x)>10'),
#                       (100,'sep>100mm','abs(ep_x-em_x)>100')]}
#    },
    

    
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


                            varp='ep'
                            varm='em'
                            if abs(pid1) == 13:
                                varp='mp'
                                varm='mm'
                            elif abs(pid1) == 999:
                                varp='hp'
                                varm='hm'

                            vary=varp+'_y'                                
                            print(f"Using var {vary}")


                            hall=ROOT.TH1F("hall","hall",1000,0,10000)                                
                            t_eff.Draw(f"{vary}>>hall");
                            sumall+=hall.Integral()*xs
                            print(f"pid1 = {pid1}, all =",hall.Integral(),xs,hall.Integral()*xs)
                                
                            for n,(effval,efftitle,effstring) in enumerate(setup["effs"][eff]):

                                effstring=effstring.replace('ep',varp).replace('em',varm)
                                heff=ROOT.TH1F("heff","heff",1000,0,10000)
                                t_eff.Draw(f"{vary}>>heff",effstring)

                                print(f"pid1 = {pid1},",effstring,", effs[{n}] =",heff.Integral(),xs,heff.Integral()*xs)
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

                # Run HepMC creation
                if mode=="run" and do_hepmc:
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

                if mode=="plotsep":
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
