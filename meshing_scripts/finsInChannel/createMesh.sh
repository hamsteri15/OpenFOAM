python createCylinders.py > geometry/cyls

cd geometry


#gmsh -2 pin.geo

gmsh -2 pin.geo

wait

gmsh -2 pin.geo -o pin.vtk



#wait
#NekMesh -v pin.msh mesh.xml


#mv mesh.xml ../base/
