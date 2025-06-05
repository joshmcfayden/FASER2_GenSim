import pyhepmc
import matplotlib.pyplot as plt
import sys
import parse
  

#print(sys.argv)

if len(sys.argv) != 2:
    print("ERROR: Usage = python read_hepmc.py <filename>")

infile=sys.argv[1]


template = "{run}/events_14TeV_m{mass}GeV_c{coupling}_to_{decay}_s1.hepmc"

params = parse.parse(template, infile).named
print("INFO: found parameters:", params)

run=params["run"]
mass=params["mass"]
coupling=params["coupling"]
decay=params["decay"]



x=[]
y=[]

with pyhepmc.open(infile) as f:

    for event in f:
        #print(event)
        #print(event.cross_section.xsec())

        for v in event.vertices:
            #print(v.position.x,v.position.y)
            x.append(v.position.x)
            y.append(v.position.y)




plt.hist2d(x, y, bins=(50, 50), range=[[-3000, 3000], [-3000, 3000]], cmap=plt.cm.viridis,cmin=1e-9)
#plt.show()
plt.colorbar()

plt.subplots_adjust(left=0.12, right=0.97, bottom=0.10, top=0.95)
outfile=f"{run}/events_14TeV_m{mass}GeV_c{coupling}_to_{decay}_s1.png"
print(f"INFO: Writing outfile: {outfile}")
plt.savefig(outfile)
