def read_vtk(filename) :
	import vtk
	import numpy as np
	reader=vtk.vtkUnstructuredGridReader()
	reader.SetFileName(filename)
	reader.Update()
	nbvv=reader.GetNumberOfVectorsInFile()
	nbss=reader.GetNumberOfScalarsInFile()
	m1=reader.GetOutput()
	numpoints=m1.GetNumberOfPoints()
	points=m1.GetPoints()
	xyz=np.zeros((numpoints,3))
	for k in range(numpoints):
		xyz[k,0],xyz[k,1],xyz[k,2]=points.GetPoint(k)
	pointData=m1.GetPointData()
	nba=pointData.GetNumberOfArrays()
	S=np.zeros((numpoints,3*nba*nbvv+nbss*nba))
	l=0
	for i in range(nbss*nba):
		a=pointData.GetArrayName(i)
		ve=pointData.GetScalars(a)
		for k in range(numpoints):
			S[k,0+l]=ve.GetTuple1(k)
		l=l+1
	ll=0
	for i in range(nbvv*nba):
		a=pointData.GetArrayName(i)
		ve=pointData.GetVectors(a)
		for k in range(numpoints):
			S[k,l+ll*3],S[k,l+1+ll*3],S[k,l+2+ll*3]=ve.GetTuple3(k)
		ll=ll+1

	return xyz[:numpoints/2,:],S[:numpoints/2,:],nba*nbvv,nba*nbss

  #######################################################################"
def ecrivtkbis(filename,PHI,sortie,change_domaine,nom,nom_champ):
# Pour ecrire sur le maillage de code_saturne (ou sur un maillage Ensight en general)
#Ecriture en VTK d un fichier EnsightGoldBinaries
# en entree
# filename : l adresse du fichier ( ex :   "/Calcul/eliberge/DATA/240_SNAP/CHR.case")
# PHI : variable a ecrire PHI(nombre de cellules, nombre de mode, dimension) ou PHI(nombre de cellules, dimension) (ou PHI(nombre de cellules, nombre de modes) ) ou PHI(nombre de cellules)
# sortie : repertoire de sortie
#variable : nom a donner a la variable

	import vtk
	#~ del(grid)
	#~ del(vectors)
	#~ print variable
	#~ print PHI[0,:,0]
	#~ print PHI[0,:,1]
	if change_domaine==1:
		vtkfile=vtk.vtkUnstructuredGridReader()
		vtkfile.SetFileName('%s/mask.vtk' %(sortie))
		vtkfile.Update()
		grid=vtkfile.GetOutput()

	else :
		reader=vtk.vtkEnSightGoldBinaryReader()
		reader.SetCaseFileName(filename)
		reader.Update()
		mstructure=reader.GetOutput().GetBlock(0)
		grid = vtk.vtkUnstructuredGrid()
		grid.SetPoints(mstructure.GetPoints())
		grid.SetCells(mstructure.GetCellTypesArray(),mstructure.GetCellLocationsArray(), mstructure.GetCells())
		volume=mstructure.GetCellData().GetArray('Volume_cellule')

	if len(PHI.shape) == 3 :
		nbve=PHI[:,0,0].size
		nmode=PHI[0,:,0].size
		ndim=PHI[0,0,:].size
		for i in range(nmode):
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(3)
			vectors.SetNumberOfTuples(nbve)
			vectors.SetName("%s%d" %(nom_champ,i))
			for k in range(nbve):
				vectors.InsertTuple3(k, PHI[k,i,0],PHI[k,i,1],PHI[k,i,2])
			grid.GetCellData().AddArray(vectors)
			grid.GetCellData().SetActiveVectors("%s%d" %(nom_champ,i))
		#~ grid.GetCellData().AddArray(volume)
		#~ grid.GetCellData().SetActiveScalars("Volume_cellule")
		vtkfile=vtk.vtkUnstructuredGridWriter()
		vtkfile.SetInput(grid)
		vtkfile.SetFileName('%s/%s.vtk' %(sortie,nom))
		vtkfile.Write()
	else :

		if len(PHI.shape) == 2 :
			if PHI[0,:].size == 3 :
				nbve=PHI[:,0].size
				ndim=PHI[0,:].size
				vectors = vtk.vtkFloatArray()
				vectors.SetNumberOfComponents(3)
				vectors.SetNumberOfTuples(nbve)
				vectors.SetName("%s" %(nom_champ))
				for k in range(nbve):
					vectors.InsertTuple3(k, PHI[k,0],PHI[k,1],PHI[k,2])
				grid.GetCellData().AddArray(vectors)
				grid.GetCellData().SetActiveVectors("%s" %(nom_champ))
				#~ grid.GetCellData().AddArray(volume)
				#~ grid.GetCellData().SetActiveScalars("volume")
				vtkfile=vtk.vtkUnstructuredGridWriter()
				vtkfile.SetInput(grid)
				vtkfile.SetFileName('%s/%s.vtk' %(sortie,nom))
				vtkfile.Write()
			else :
				nbve=PHI[:,0].size
				nmode=PHI[0,:].size
				for i in range(nmode):
					vectors = vtk.vtkFloatArray()
					vectors.SetNumberOfComponents(1)
					vectors.SetNumberOfTuples(nbve)
					vectors.SetName("%s%d" %(nom_champ,i))
					for k in range(nbve):
						vectors.InsertTuple1(k, PHI[k,i])
					grid.GetCellData().AddArray(vectors)
					grid.GetCellData().SetActiveScalars("%s%d" %(nom_champ,i))
				#~ grid.GetCellData().AddArray(volume)
				#~ grid.GetCellData().SetActiveScalars("volume")
				vtkfile=vtk.vtkUnstructuredGridWriter()
				vtkfile.SetInput(grid)
				vtkfile.SetFileName('%s/%s.vtk' %(sortie,nom))
				vtkfile.Write()
		else :
			nbve=PHI[:].size
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(1)
			vectors.SetNumberOfTuples(nbve)
			vectors.SetName("%s" %(nom_champ))
			for k in range(nbve):
				vectors.InsertTuple1(k, PHI[k])
			grid.GetCellData().AddArray(vectors)
			grid.GetCellData().SetActiveScalars("%s" %(nom_champ))
			#~ grid.GetCellData().AddArray(volume)
			#~ grid.GetCellData().SetActiveScalars("volume")
			vtkfile=vtk.vtkUnstructuredGridWriter()
			vtkfile.SetInput(grid)
			vtkfile.SetFileName('%s/%s.vtk' %(sortie,nom))
			vtkfile.Write()


#######################################################################"

def ecrivtk2(nx,ny,nz,F,PHI,sortie):
#Ecriture sur un maillage cartesien faussement 2D
# en entree :
# nx :nombre de noeuds selon x
# ny :nombre de noeuds selon y
# nz :nombre de noeuds selon z
# F : coordonnee des noeuds du maillage : x1 y1
#                                                                x2 y2
#
#                                                                xn yn
#F doit etre construit ainsi :
#les nx premieres lignes voient seulemnt varier les xi
# PHI la variable a ecrire :
#PHI(nombre de noeuds (en 2D), nombre de mode, dimension) ou PHI(nombre de noeuds (en 2D), dimension)
# sortie :repertoire de sortie

	import vtk
	import scipy as sci
	numpoint=F[:,0].size

	numpoints=numpoint*2
	points=vtk.vtkPoints()
	points.SetNumberOfPoints(numpoints)
	#~ CO=np.zeros((numpoints,4))
	for i in range(numpoint):
		points.InsertPoint(i,F[i,0],F[i,1],0.0)
		points.InsertPoint(i+numpoint,F[i,0],F[i,1],1.0)
		#~ CO[i,:]=[i,F[i,0],F[i,1],0.0]
		#~ CO[i+numpoint,:]=[i+numpoint,F[i,0],F[i,1],1.0]


	grid=vtk.vtkUnstructuredGrid()
	grid.SetPoints(points)

	numberofcells=(nx-1)*(ny-1)*(nz-1)
	for lz in range(nz-1):
		for ly in range(ny-1):
			for lx in range(nx-1):
				i=lx+nx*ly+(nx*ny)*lz
				ids=vtk.vtkIdList()
				ids.InsertNextId(i)
				ids.InsertNextId(i+1)
				ids.InsertNextId(i+1+nx)
				ids.InsertNextId(i+nx)
				ids.InsertNextId(i+nx*ny)
				ids.InsertNextId(i+nx*ny+1)
				ids.InsertNextId(i+nx*ny+1+nx)
				ids.InsertNextId(i+nx*ny+nx)
				grid.InsertNextCell(vtk.VTK_HEXAHEDRON, ids)

	if len(PHI.shape) == 3 :
		nbve=PHI[:,0,0].size
		nmode=PHI[0,:,0].size
		ndim=PHI[0,0,:].size

		listemode=sci.arange(nmode)

		for i1 in range(nmode):
			i=listemode[(i1+1) %(nmode)]
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(3)
			vectors.SetNumberOfTuples(nbve)
			vectors.SetName("phi%d" %(i))
			for k in range(nbve):
				vectors.InsertTuple3(k, PHI[k,i,0],PHI[k,i,1],PHI[k,i,2])
				vectors.InsertTuple3(k+numpoint, PHI[k,i,0],PHI[k,i,1],PHI[k,i,2])
			grid.GetPointData().AddArray(vectors)
			grid.GetPointData().SetActiveVectors("phi%d" %(i))

		vtkfile=vtk.vtkUnstructuredGridWriter()
		vtkfile.SetInput(grid)
		vtkfile.SetFileName('%s/basepod.vtk' %(sortie))
		vtkfile.Write()
	else :
		nbve=PHI[:,0].size
		ndim=PHI[0,:].size
		vectors = vtk.vtkFloatArray()
		vectors.SetNumberOfComponents(3)
		vectors.SetNumberOfTuples(nbve)
		vectors.SetName("vm")
		for k in range(nbve):
			vectors.InsertTuple3(k, PHI[k,0],PHI[k,1],PHI[k,2])
			vectors.InsertTuple3(k+numpoint, PHI[k,0],PHI[k,1],PHI[k,2])
		grid.GetPointData().AddArray(vectors)
		grid.GetPointData().SetActiveVectors("vm")

		vtkfile=vtk.vtkUnstructuredGridWriter()
		vtkfile.SetInput(grid)
		vtkfile.SetFileName('%s/vm.vtk' %(sortie))
		vtkfile.Write()
######################################################################################
def ecrivtk3(nx,ny,nz,F,PHI,sortie):
#Ecriture sur un maillage cartesien 3D
# en entree :
# nx :nombre de noeuds selon x
# ny :nombre de noeuds selon y
# nz :nombre de noeuds selon z
# F : coordonnee des noeuds du maillage : x1 y1 z1
#                                                                x2 y2 z2
#
#                                                                xn yn zn
#F doit etre construit ainsi :
#les nx premieres lignes voient seulemnt varier les xi
#sur les nx*ny premieres lignes les zi ne varient pas
# ensuite on passe au deuxieme "etage" et ainsi de suite
# PHI la variable a ecrire :
#PHI(nombre de noeuds (en 2D), nombre de mode, dimension) ou PHI(nombre de noeuds (en 2D), dimension)
# sortie :repertoire de sortie
	import vtk
	numpoints=F[:,0].size


	points=vtk.vtkPoints()
	points.SetNumberOfPoints(numpoints)
	#~ CO=np.zeros((numpoints,4))
	for i in range(numpoints):
		points.InsertPoint(i,F[i,0],F[i,1],F[i,2])

		#~ CO[i,:]=[i,F[i,0],F[i,1],0.0]
		#~ CO[i+numpoint,:]=[i+numpoint,F[i,0],F[i,1],1.0]


	grid=vtk.vtkUnstructuredGrid()
	grid.SetPoints(points)

	numberofcells=(nx-1)*(ny-1)*(nz-1)
	for lz in range(nz-1):
		for ly in range(ny-1):
			for lx in range(nx-1):
				i=lx+nx*ly+(nx*ny)*lz
				ids=vtk.vtkIdList()
				ids.InsertNextId(i)
				ids.InsertNextId(i+1)
				ids.InsertNextId(i+1+nx)
				ids.InsertNextId(i+nx)
				ids.InsertNextId(i+nx*ny)
				ids.InsertNextId(i+nx*ny+1)
				ids.InsertNextId(i+nx*ny+1+nx)
				ids.InsertNextId(i+nx*ny+nx)
				grid.InsertNextCell(vtk.VTK_HEXAHEDRON, ids)

	if len(PHI.shape) == 3 :
		nbve=PHI[:,0,0].size
		nmode=PHI[0,:,0].size
		ndim=PHI[0,0,:].size
		listemode=sci.arange(nmode)

		for i1 in range(nmode):
			i=listemode[(i1+1) %(nmode)]
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(3)
			vectors.SetNumberOfTuples(nbve)
			vectors.SetName("phi%d" %(i))
			for k in range(nbve):
				vectors.InsertTuple3(k, PHI[k,i,0],PHI[k,i,1],PHI[k,i,2])
			grid.GetPointData().AddArray(vectors)
			grid.GetPointData().SetActiveVectors("phi%d" %(i))

		vtkfile=vtk.vtkUnstructuredGridWriter()
		vtkfile.SetInput(grid)
		vtkfile.SetFileName('%s/basepod.vtk' %(sortie))
		vtkfile.Write()
	else :
		nbve=PHI[:,0].size
		ndim=PHI[0,:].size
		vectors = vtk.vtkFloatArray()
		vectors.SetNumberOfComponents(3)
		vectors.SetNumberOfTuples(nbve)
		vectors.SetName("vm")
		for k in range(nbve):
			vectors.InsertTuple3(k, PHI[k,0],PHI[k,1],PHI[k,2])

		grid.GetPointData().AddArray(vectors)
		grid.GetPointData().SetActiveVectors("vm")

		vtkfile=vtk.vtkUnstructuredGridWriter()
		vtkfile.SetInput(grid)
		vtkfile.SetFileName('%s/vm.vtk' %(sortie))
		vtkfile.Write()
#########################################################################################
def xyztodata2(nx,ny,nz,F,v,sortie,variable) :
#Ecriture sur un maillage cartesien faussement 2D
# en entree :
# nx :nombre de noeuds selon x
# ny :nombre de noeuds selon y
# nz :nombre de noeuds selon z
# F : coordonnee des noeuds du maillage : x1 y1
#                                                                x2 y2
#
#                                                                xn yn
#F doit etre construit ainsi :
#les nx premieres lignes voient seulemnt varier les xi
# v la variable a ecrire :
#v(nombre de cliches, nombre de noeuds (nx*ny), dimension) ou v(nx*ny, valeurs nodales)
# UN FICHIER EST CREER PAR SNAPSHOT
# sortie :repertoire de sortie
# nom de la variable
	import vtk
	if len(v.shape) == 3 :
		M=v[:,0,0].size

		numpoint=F[:,0].size

		numpoints=numpoint*2
		points=vtk.vtkPoints()
		points.SetNumberOfPoints(numpoints)
		#~ CO=np.zeros((numpoints,4))
		for i in range(numpoint):
			points.InsertPoint(i,F[i,0],F[i,1],0.0)
			points.InsertPoint(i+numpoint,F[i,0],F[i,1],1.0)
		#~ CO[i,:]=[i,F[i,0],F[i,1],0.0]
		#~ CO[i+numpoint,:]=[i+numpoint,F[i,0],F[i,1],1.0]


		grid=vtk.vtkUnstructuredGrid()
		grid.SetPoints(points)

		numberofcells=(nx-1)*(ny-1)*(nz-1)
		for lz in range(nz-1):
			for ly in range(ny-1):
				for lx in range(nx-1):
					i=lx+nx*ly+(nx*ny)*lz
					ids=vtk.vtkIdList()
					ids.InsertNextId(i)
					ids.InsertNextId(i+1)
					ids.InsertNextId(i+1+nx)
					ids.InsertNextId(i+nx)
					ids.InsertNextId(i+nx*ny)
					ids.InsertNextId(i+nx*ny+1)
					ids.InsertNextId(i+nx*ny+1+nx)
					ids.InsertNextId(i+nx*ny+nx)
					grid.InsertNextCell(vtk.VTK_HEXAHEDRON, ids)
			#~ print "cellule"
#Creation maillage surfacique
		for t in range(M):
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(3)
			vectors.SetNumberOfTuples(numpoints)
			vectors.SetName("vint" )
			for i in range(numpoint):
				vectors.InsertTuple3(i, v[t,i,0],v[t,i,1],0)
				vectors.InsertTuple3(i+numpoint, v[t,i,0],v[t,i,1],0)
			grid.GetPointData().AddArray(vectors)
			grid.GetPointData().SetActiveVectors("vint")
			vtkfile=vtk.vtkUnstructuredGridWriter()
			vtkfile.SetInput(grid)
			vtkfile.SetFileName('%s/%sint%d.vtk' %(sortie,variable,t))
			vtkfile.Write()
	if len(v.shape) == 2 :
		M=v[:,0].size

		numpoint=F[:,0].size

		numpoints=numpoint*2
		points=vtk.vtkPoints()
		points.SetNumberOfPoints(numpoints)
		#~ CO=np.zeros((numpoints,4))
		for i in range(numpoint):
			points.InsertPoint(i,F[i,0],F[i,1],0.0)
			points.InsertPoint(i+numpoint,F[i,0],F[i,1],1.0)
		#~ CO[i,:]=[i,F[i,0],F[i,1],0.0]
		#~ CO[i+numpoint,:]=[i+numpoint,F[i,0],F[i,1],1.0]


		grid=vtk.vtkUnstructuredGrid()
		grid.SetPoints(points)

		numberofcells=(nx-1)*(ny-1)*(nz-1)
		for lz in range(nz-1):
			for ly in range(ny-1):
				for lx in range(nx-1):
					i=lx+nx*ly+(nx*ny)*lz
					ids=vtk.vtkIdList()
					ids.InsertNextId(i)
					ids.InsertNextId(i+1)
					ids.InsertNextId(i+1+nx)
					ids.InsertNextId(i+nx)
					ids.InsertNextId(i+nx*ny)
					ids.InsertNextId(i+nx*ny+1)
					ids.InsertNextId(i+nx*ny+1+nx)
					ids.InsertNextId(i+nx*ny+nx)
					grid.InsertNextCell(vtk.VTK_HEXAHEDRON, ids)
			#~ print "cellule"
#Creation maillage surfacique
		for t in range(M):
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(1)
			vectors.SetNumberOfTuples(numpoints)
			vectors.SetName("fc" )
			for i in range(numpoint):
				vectors.InsertTuple1(i, v[t,i])
				vectors.InsertTuple1(i+numpoint, v[t,i])
			grid.GetPointData().AddArray(vectors)
			grid.GetPointData().SetActiveScalars("fc")
			vtkfile=vtk.vtkUnstructuredGridWriter()
			vtkfile.SetInput(grid)
			vtkfile.SetFileName('%s/%sint%d.vtk' %(sortie,variable,t))
			vtkfile.Write()

###############################################################
def xyztodata3(nx,ny,nz,F,v,sortie,variable) :
#Ecriture sur un maillage cartesien 3D
# en entree :
# nx :nombre de noeuds selon x
# ny :nombre de noeuds selon y
# nz :nombre de noeuds selon z
# F : coordonnee des noeuds du maillage : x1 y1 z1
#                                                                x2 y2 z2
#
#                                                                xn yn zn
#F doit etre construit ainsi :
#les nx premieres lignes voient seulemnt varier les xi
#sur les nx*ny premieres lignes les zi ne varient pas
# ensuite on passe au deuxieme "etage" et ainsi de suite
# v la variable a ecrire : [V(t,x).x1 V(t,x).x2 V(t,x).x3]  ou [f(t,x)]
#v(nombre de cliches, nombre de noeuds (nx*ny*nz), dimension) ou v((nx*ny*nz), valeurs nodales)
# UN FICHIER EST CREER PAR SNAPSHOT
# sortie :repertoire de sortie
# variable : nom de la variable
	import vtk
	if len(v.shape) == 3 :
		M=v[:,0,0].size

		numpoints=F[:,0].size


		points=vtk.vtkPoints()
		points.SetNumberOfPoints(numpoints)
		#~ CO=np.zeros((numpoints,4))
		for i in range(numpoints):
			points.InsertPoint(i,F[i,0],F[i,1],F[i,2])
		#~ CO[i,:]=[i,F[i,0],F[i,1],0.0]
		#~ CO[i+numpoint,:]=[i+numpoint,F[i,0],F[i,1],1.0]


		grid=vtk.vtkUnstructuredGrid()
		grid.SetPoints(points)

		numberofcells=(nx-1)*(ny-1)*(nz-1)
		for lz in range(nz-1):
			for ly in range(ny-1):
				for lx in range(nx-1):
					i=lx+nx*ly+(nx*ny)*lz
					ids=vtk.vtkIdList()
					ids.InsertNextId(i)
					ids.InsertNextId(i+1)
					ids.InsertNextId(i+1+nx)
					ids.InsertNextId(i+nx)
					ids.InsertNextId(i+nx*ny)
					ids.InsertNextId(i+nx*ny+1)
					ids.InsertNextId(i+nx*ny+1+nx)
					ids.InsertNextId(i+nx*ny+nx)
					grid.InsertNextCell(vtk.VTK_HEXAHEDRON, ids)
			#~ print "cellule"
#Creation maillage surfacique
		for t in range(M):
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(3)
			vectors.SetNumberOfTuples(numpoints)
			vectors.SetName("vint" )
			for i in range(numpoints):
				vectors.InsertTuple3(i, v[t,i,0],v[t,i,1],v[t,i,2])
			grid.GetPointData().AddArray(vectors)
			grid.GetPointData().SetActiveVectors("vint")
			vtkfile=vtk.vtkUnstructuredGridWriter()
			vtkfile.SetInput(grid)
			vtkfile.SetFileName('%s/%sint%d.vtk' %(sortie,variable,t))
			vtkfile.Write()
	if len(v.shape) == 2 :
		M=v[:,0].size

		numpoints=F[:,0].size


		points=vtk.vtkPoints()
		points.SetNumberOfPoints(numpoints)
		#~ CO=np.zeros((numpoints,4))
		for i in range(numpoints):
			points.InsertPoint(i,F[i,0],F[i,1],F[i,2])

		#~ CO[i,:]=[i,F[i,0],F[i,1],0.0]
		#~ CO[i+numpoint,:]=[i+numpoint,F[i,0],F[i,1],1.0]


		grid=vtk.vtkUnstructuredGrid()
		grid.SetPoints(points)

		numberofcells=(nx-1)*(ny-1)*(nz-1)
		for lz in range(nz-1):
			for ly in range(ny-1):
				for lx in range(nx-1):
					i=lx+nx*ly+(nx*ny)*lz
					ids=vtk.vtkIdList()
					ids.InsertNextId(i)
					ids.InsertNextId(i+1)
					ids.InsertNextId(i+1+nx)
					ids.InsertNextId(i+nx)
					ids.InsertNextId(i+nx*ny)
					ids.InsertNextId(i+nx*ny+1)
					ids.InsertNextId(i+nx*ny+1+nx)
					ids.InsertNextId(i+nx*ny+nx)
					grid.InsertNextCell(vtk.VTK_HEXAHEDRON, ids)
			#~ print "cellule"
#Creation maillage surfacique
		for t in range(M):
			vectors = vtk.vtkFloatArray()
			vectors.SetNumberOfComponents(1)
			vectors.SetNumberOfTuples(numpoints)
			vectors.SetName("fc" )
			for i in range(numpoints):
				vectors.InsertTuple1(i, v[t,i])
			grid.GetPointData().AddArray(vectors)
			grid.GetPointData().SetActiveScalars("fc")
			vtkfile=vtk.vtkUnstructuredGridWriter()
			vtkfile.SetInput(grid)
			vtkfile.SetFileName('%s/%sint%d.vtk' %(sortie,variable,t))
			vtkfile.Write()
