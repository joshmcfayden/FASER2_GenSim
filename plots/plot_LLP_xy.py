import pyhepmc
import matplotlib.pyplot as plt
import os


run="R1-L10-R1x3_DarkHiggs"
mass="1.5849"
#mass="0.1"
#mass="0.5012"
#coupling="2.06913808111479e-06"
#coupling="0.0001623776739188721"
coupling="1.274274985703132e-05"
decayname="mu_mu"
decays=["mu_mu"]
#decayname="all"
#decays=["mu_mu","K_K","pi_pi"]#,"4pi"]

#run="R1-L10-R1x3_DarkPhoton"
#mass="0.1585"
#coupling="1.2742749857031322e-06"
#decaysname="e_e"
##decays=["-11_11"]
#decays=["e_e"]



x=[]
y=[]
xs=[]

# pyhepmc.open can read most HepMC formats using auto-detection
#with pyhepmc.open("/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/mytest_G4_mini/events_14TeV_m0.1GeV_c1e-05to_11_-11_test1.hepmc") as f:
#with pyhepmc.open("/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkPhoton//events_14TeV_m0.01GeV_c1e-08to_-11_11_s1.hepmc") as f:
#with pyhepmc.open("R1-L10-R1x3_DarkHiggs//events_14TeV_m3.5GeV_c2.976351441631319e-06_to_4pi_s1.hepmc") as f:
#with pyhepmc.open("R1-L10-R1x3_DarkHiggs//events_14TeV_m0.955GeV_c7.847599703514606e-05_to_mu_mu_s1.hepmc") as f:
#with pyhepmc.open("R1-L10-R1x3_DarkHiggs//events_14TeV_m1.5849GeV_c2.06913808111479e-06_to_pi_pi_s1.hepmc") as f:

files=[f"../{run}/events_14TeV_m{mass}GeV_c{coupling}_to_{decay}_s1.hepmc" for decay in decays]
print(files)

tmpxs=None
minxs=None

for infile in files:
    if not os.path.exists(infile):
        print(f"WARNING: File {infile} not found")
        continue
    
    with pyhepmc.open(infile) as f:
        #print(f)
        #print(f.read())
        #event = f.read()
        for n,event in enumerate(f):
            if n==0:
                tmpxs=event.cross_section.xsec()
                print(tmpxs)
                if not minxs or tmpxs < minxs:
                    minxs=tmpxs
                    
            #print(event)
            #print(event.cross_section.xsec())
            #print()
            for v in event.vertices:
                #print(v)
                #if tmpxs: > 1e-10:
                x.append(v.position.x)
                y.append(v.position.y)
                xs.append(tmpxs)
            #for p in event.particles:
            #    print(p)

        f.close()

print("minxs:",minxs)
print("xs:",xs[0],xs[-1])
wgt=[wxs/minxs for wxs in xs]
print("wgt:",wgt[0],wgt[-1])
#wgt=xs
#plt.hist2d(y, x, bins=(300, 100), range=[[-1000, 1000], [-1000, 1000]], cmap=plt.cm.jet)
#plt.hist2d(x, y, bins=(49, 51), range=[[-3000, 3000], [-3000, 3000]], cmap=plt.cm.viridis,cmin=1e-9)
#plt.hist2d(x, y, bins=(50, 50), range=[[-3000, 3000], [-3000, 3000]], cmap=plt.cm.viridis,cmin=1e-9)
plt.hist2d(x, y, bins=(50, 50), range=[[-3000, 3000], [-3000, 3000]], cmap=plt.cm.viridis,cmin=1e-9,weights=wgt)
#plt.show()
        
plt.subplots_adjust(left=0.12, right=0.97, bottom=0.10, top=0.95)
plt.savefig(f"../{run}/LLP_xy_m{mass}GeV_c{coupling}_to_{decayname}.pdf")
