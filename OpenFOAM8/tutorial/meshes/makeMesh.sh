# source $HOME/foam/foam-extend-4.0/etc/bashrc

#copire il file geo con il nome untitled.geo prima di lanciare lo script

rm -rf *.msh
gmsh mesh.geo -3 -o untitled.msh
rm -rf case
cp -r config case
gmshToFoam untitled.msh -case case
rm -rf *.msh

sed --in-place '22s/patch/empty/' $(pwd)/case/constant/polyMesh/boundary
# sed --in-place '43s/patch/empty/' $(pwd)/case/constant/polyMesh/boundary
# sed --in-place '40 s/patch/empty/' $(pwd)/case/constant/polyMesh/boundary
# sed --in-place '50s/patch/wall/' $(pwd)/case/constant/polyMesh/boundary
# sed --in-place '57s/patch/wall/' $(pwd)/case/constant/polyMesh/boundary
checkMesh -case case

# polyDualMesh 180 -case case
# rm *.obj
# checkMesh -case case
