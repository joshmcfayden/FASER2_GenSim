import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 15})

#m,c,effcut,eff = np.loadtxt('/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkPhoton/eff_14TeV_FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_3rdTrkStation_envx.csv',delimiter=',',unpack=True)

xlims = None
ylims= None
xlabel= None
ylabel= None
legendloc= None
figsize= None

title="Dark Photons"
title="Dark Higgs"

if title=="Dark Photon":
    xlims = [0.01,3]
    ylims=[10**-7,0.002]
    xlabel=r"Dark Photon Mass $m_{A'}$ [GeV]"
    ylabel=r"Kinetic Mixing $\epsilon$"
    legendloc=(1.02,0.72)
    figsize=(10,6)
#csvfile='/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkPhoton/eff_14TeV_FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_3rdTrkStation_envx.csv'



if title=="Dark Higgs":
    xlims=[0.1,10]
    ylims=[10**-5.5,10**-2.0]
    xlabel=r"Dark Higgs Mass $m_{\phi}$ [GeV]"
    ylabel=r"Mixing $\theta$"
    legendloc=(0.43,0.25)
    figsize=(8,6)
      
#csvfile='/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkHiggs//eff_14TeV_FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_AllTrkStations_station3_sep.csv'
csvfile='/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkHiggs//eff_14TeV_FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_AllTrkStations_station3_env.csv'
csvfile='/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkHiggs//eff_14TeV_FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_AllTrkStations_station1_sep.csv'

df= pd.read_csv(csvfile,names=['masses', 'couplings', 'effcuts','effs'],header=None)

unique_effcuts=df.effcuts.unique()
min_eff=df.effs.min()



for pick_effcut in unique_effcuts:
    #pick_effcut=-500
    pick_df = df[df['effcuts']==pick_effcut]
    
    # = [x,y,z for x,y,tmpeffcut,z in zip(m,c,effcut,eff) if tmpeffcut == pick_effcut]
    #
    
    unique_m=pick_df.masses.unique()
    unique_c=pick_df.couplings.unique()
    nx=len(unique_m)
    ny=len(unique_c)
    print("nx:",nx,"ny:",ny)
    
    #print(pick_df.masses)#.values.reshape(nx,ny).T,
    #print(pick_df['masses'])#.values.reshape(nx,ny).T,
    #print(pick_df['masses'].values)#.reshape(nx,ny).T,
    print(pick_df['masses'].values.reshape(nx,ny).T)
    print(pick_df['couplings'].values.reshape(nx,ny).T)
    print(pick_df['effs'].values.reshape(nx,ny).T)
    
    
    
    
    fig, ax = plt.subplots(figsize=figsize)
    #ax.pcolormesh(np.log(pick_df['masses'].values).reshape(nx,ny).T,
    #              np.log(pick_df['couplings'].values).reshape(nx,ny).T,
    #              pick_df['effs'].values.reshape(nx,ny).T)
    
    im=ax.pcolormesh(pick_df['masses'].values.reshape(nx,ny).T,
                     pick_df['couplings'].values.reshape(nx,ny).T,
                     pick_df['effs'].values.reshape(nx,ny).T,
                     cmap=plt.cm.viridis,
                     vmin=min_eff, vmax=1.0)

    fig.colorbar(im, ax=ax)
    
    
    ax.set_title(title)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xlims[0],xlims[1])
    ax.set_ylim(ylims[0],ylims[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

        
    #plt.show()

    print(pick_effcut)
    plt.savefig(csvfile.replace('.csv','')+'_'+str(pick_effcut)+'.pdf')
    #plt.savefig('test.pdf')
