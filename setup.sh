export PYTHONPATH=$PYTHONPATH:$PWD/FORESEE/src/
#source /usr/local/Cellar/root/6.24.06/bin/thisroot.sh
source /usr/local/Cellar/root/6.26.06_2//bin/thisroot.sh

pathtoG4="../G4_test"
G4version="geant4.10.07.p02"

cd ${pathtoG4}/${G4version}-install/bin
source geant4.sh
cd -
