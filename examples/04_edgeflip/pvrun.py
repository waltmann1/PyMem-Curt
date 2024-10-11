import os, glob
from paraview.simple import *


def paraviewvis(viewindex, facefiles):
	print(viewindex, facefiles)
	if viewindex==0: 
		renderView=GetActiveViewOrCreate('RenderView')
	else:
		renderView=CreateRenderView()
	layout1 = GetLayoutByName("Layout #1")
	print("viewindex:",viewindex)
	AssignViewToLayout(view=renderView, layout=layout1, hint=viewindex)
	facereader = XMLPolyDataReader(FileName=facefiles)
	animationScene = GetAnimationScene()
	animationScene.UpdateAnimationUsingDataTimeSteps()
	programmableFilter = ProgrammableFilter(Input=facereader)
	programmableFilter.Script = """import numpy as np
input = inputs[0]
Numfaces = input.GetNumberOfCells()
neighbor = np.zeros(input.GetNumberOfPoints())
for i in range(Numfaces):
       cell = input.GetCell(i)
       neighbor[cell.GetPointId(0)] += 1
       neighbor[cell.GetPointId(1)] += 1
       neighbor[cell.GetPointId(2)] += 1
output.PointData.append(neighbor, "neighbor")"""
	programmabledisplay = Show(programmableFilter,renderView)
	ColorBy(programmabledisplay, ('POINTS', 'neighbor'))
	vertexTypeLUT = GetColorTransferFunction('neighbor')
	vertexTypeLUT.RGBPoints = [5.0, 0.231373, 0.298039, 0.752941, 6.0, 0.705882, 0.0156863, 0.14902, 7.0, 0.17647058823529413, 0.596078431372549, 0.13333333333333333]
	renderView.OrientationAxesVisibility = 0
	programmabledisplay.SetScalarBarVisibility(renderView, False)

	pythonAnnotation = PythonAnnotation(Input=facereader)
	pythonAnnotation.ArrayAssociation = 'Point Data'
	pythonAnnotation.Expression = " 'Nv = %d' % len(inputs[0].PointData['Id'])"
	pythonAnnotationDisplay = Show(pythonAnnotation, renderView)
	pythonAnnotationDisplay.WindowLocation = 'LowerCenter'
	pythonAnnotationDisplay.FontSize = 16

	renderView.ResetCamera()
	Render()
	return renderView


def getint(filename):
	try:
		return int(filename.split(".vtp")[0].split("_")[-1])
	except:
		return -1

path = os.getcwd()
for viewindex, Nv in enumerate([32, 42, 72, 132]):
	facefiles = sorted(glob.glob('*edge_flip_Nv_'+ str(Nv) + '*.vtp'), key = getint)
	paraviewvis(viewindex, [os.path.join(path, facefile) for facefile in facefiles])

servermanager.SaveState("movie_paraview.pvsm")
