# FASER2_GenSim

## Setup

Checkout this package:

```bash
git clone https://github.com/joshmcfayden/FASER2_GenSim.git
cd FASER2_GenSim/
```

Checkout FORESEE and apply patch:
```bash
git clone https://github.com/KlingFelix/FORESEE.git
cd FORESEE
git checkout d359a4d5b56cff9b797a2df8c2cd2de3deb6da47
patch src/foresee.py ../foresee.patch
cd ..
```

Setup environment for ROOT and G4 on MacOS:
Edit `$pathtoG4` and `$G4version` in `setup.sh` then source:
```bash
source setup.sh
```

Setup environment for ROOT and G4 on lxplus:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc10-opt/setup.sh
```

Setup G4 simulations:
Follow the instructions here: https://github.com/joshmcfayden/FASER2_G4/tree/faser2_newbaseline_KEK/ to compile G4 simulations

## Running all the steps

Setup the correct paths for your installations in `run.py`, in particular the path to FORESEE, `path` and the path to the G4 setups, `G4path`.

Then set the mode you want:
- `"run"`: Run FORESEE event generation, outputs `.npy` files.  If `do_hepmc=True` this will also output `.hepmc` files to pass to G4.
- `"G4"`: Run Geant4 simulation based on `"G4"` entries in `setup_dict`
- `"combine"`: Combines together the signal yields for e.g. different decay modes and puts into new `.npy` in the format expected by FORESEE's `plot_reach` function.
- `"eff"`: Calculate efficiencies based on criterial in `"effs"` entries in `setup_dict` using ROOT file output from running G4 step.
- `"plotsep"`: Plot the particle separations using ROOT file output from running G4 step.


This will then loop over all the `masses` and `couplings` values and execute other plots/calculations as indicated from the `mode` and `setup_dict` entries.

```bash
python3 run.py
```
