import pyhepmc
import matplotlib.pyplot as plt

x=[]
y=[]

# pyhepmc.open can read most HepMC formats using auto-detection
with pyhepmc.open("/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/mytest_G4_mini/events_14TeV_m0.1GeV_c1e-05to_11_-11_test1.hepmc") as f:
#with pyhepmc.open("/Users/mcfayden/Work/FASER/FASER2/FASER_FORESEE/R1-L10-R1x3_DarkPhoton//events_14TeV_m0.01GeV_c1e-08to_-11_11_s1.hepmc") as f:
    #event = f.read()
    for event in f:
        #print(event.cross_section.xsec())
        for v in event.vertices:
            #print(v.position.x,v.position.y)
            x.append(v.position.x)
            y.append(v.position.y)



#plt.hist2d(y, x, bins=(300, 100), range=[[-1000, 1000], [-1000, 1000]], cmap=plt.cm.jet)
plt.hist2d(x, y, bins=(300, 100), range=[[-3000, 3000], [-3000, 3000]], cmap=plt.cm.viridis,cmin=1e-9)
plt.show()
        
