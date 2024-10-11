import os, glob
from paraview.simple import *


def paraviewvis(viewindex, edgefile):
	print(viewindex, edgefile)
	if viewindex==0: 
		renderView=GetActiveViewOrCreate('RenderView')
	else:
		renderView=CreateRenderView()
	layout1 = GetLayoutByName("Layout #1")
	print("viewindex:",viewindex)
	AssignViewToLayout(view=renderView, layout=layout1, hint=viewindex)
	edgereader = XMLPolyDataReader(FileName=edgefile)
	edgedisplay=Show(edgereader,renderView)
	edgedisplay.LineWidth = 5.
	edgedisplay.DiffuseColor=[1.,1.,1.]
	ColorBy(edgedisplay, ('CELLS', 'edge_type'))
	edgeTypeLUT = GetColorTransferFunction('edge_type')
	edgeTypeLUT.ApplyPreset('X Ray', True)
	edgeTypeLUT.RGBPoints = [0.0, 0.0862745,0.0862745,0.65098, 1.0, 0.705882, 0.0156863, 0.14902]
	edgedisplay.SetScalarBarVisibility(renderView, False)
	renderView.OrientationAxesVisibility = 0
	text=Text()
	text.Text=edgefile.split("/")[-1]
	textdisplay=Show(text,renderView)
	textdisplay.FontSize = 10
	renderView.ResetCamera()
	Render()
	return renderView



path = os.getcwd()
edgefiles = glob.glob('*fraction*_edges.vtp')

for viewindex, edgefile in enumerate(edgefiles):
	paraviewvis(viewindex, os.path.join(path, edgefile))

servermanager.SaveState("snapshots_paraview.pvsm")
