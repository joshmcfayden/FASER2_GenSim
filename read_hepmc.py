import pyhepmc
import matplotlib.pyplot as plt

#run="R1-L10-R1x3_DarkHiggs"
#mass="1.5849"
#coupling="2.06913808111479e-06"
#decay="pi_pi"


run="F2-default_DarkPhoton"
mass="0.01"
coupling="1e-08"
decay="e_e"


run="F2-default_DarkHiggs"
mass="0.2239"
coupling="1.8329807108324375e-05"
decay="mu_mu"



x=[]
y=[]

# pyhepmc.open can read most HepMC formats using auto-detection
#with pyhepmc.open("/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/mytest_G4_mini/events_14TeV_m0.1GeV_c1e-05to_11_-11_test1.hepmc") as f:
#with pyhepmc.open("/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkPhoton//events_14TeV_m0.01GeV_c1e-08to_-11_11_s1.hepmc") as f:
#with pyhepmc.open("R1-L10-R1x3_DarkHiggs//events_14TeV_m3.5GeV_c2.976351441631319e-06_to_4pi_s1.hepmc") as f:
#with pyhepmc.open("R1-L10-R1x3_DarkHiggs//events_14TeV_m0.955GeV_c7.847599703514606e-05_to_mu_mu_s1.hepmc") as f:
#with pyhepmc.open("R1-L10-R1x3_DarkHiggs//events_14TeV_m1.5849GeV_c2.06913808111479e-06_to_pi_pi_s1.hepmc") as f:
with pyhepmc.open(f"{run}/events_14TeV_m{mass}GeV_c{coupling}_to_{decay}_s1.hepmc") as f:
    #print(f)
    #print(f.read())
    #event = f.read()
    for event in f:
        #print(event)
        #print(event.cross_section.xsec())
        for v in event.vertices:
            #print(v.position.x,v.position.y)
            x.append(v.position.x)
            y.append(v.position.y)



#plt.hist2d(y, x, bins=(300, 100), range=[[-1000, 1000], [-1000, 1000]], cmap=plt.cm.jet)
plt.hist2d(x, y, bins=(50, 50), range=[[-3000, 3000], [-3000, 3000]], cmap=plt.cm.viridis,cmin=1e-9)
plt.show()
        
plt.subplots_adjust(left=0.12, right=0.97, bottom=0.10, top=0.95)
plt.savefig(f"{run}/events_14TeV_m{mass}GeV_c{coupling}_to_{decay}_s1.png")
