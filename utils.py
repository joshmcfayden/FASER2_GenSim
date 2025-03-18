import numpy as np
import ROOT
import pyhepmc
import os
from array import array

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
        
        
#        masses=[0.1585]
#        couplings=[1.2742749857031322e-06]

        
        
    elif model=="DarkHiggs":
        #decays=["e_e","mu_mu","pi_pi","4pi","K_K","other"]
        #decays=["e_e","mu_mu","pi_pi","K_K"]
        #decays=["e_e"]
        decays=["mu_mu"]
        
        
        masses = [   
            0.1   ,  
            0.2239,  
            0.5012,  
#            0.6918,  
#            0.8128,  
#            0.955 ,
#            1.5849,  
#            3.5   ,
#            11.22 ,  
#            25.119
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


def plot_seps(currdir,outdir,energy,mass,coup,decay,G4setup,setup_name):
    infile=f"{currdir}/{outdir}/events_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}_s1_{G4setup}.root"
    if not os.path.exists(infile):
        print(f"ROOT file {infile} not found - skipping")
        return
                        
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
 
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

    c1.SaveAs(f"{currdir}/{outdir}/plot_sep_stations_{energy}TeV_m{mass}GeV_c{coup}_to_{decay}_s1_{G4setup}.pdf")
