import os, glob
from paraview.simple import *


def paraviewvis(viewindex, facefile):
	print(viewindex, facefile)
	if viewindex==0: 
		renderView=GetActiveViewOrCreate('RenderView')
	else:
		renderView=CreateRenderView()
	layout1 = GetLayoutByName("Layout #1")
	print("viewindex:",viewindex)
	AssignViewToLayout(view=renderView, layout=layout1, hint=viewindex)
	facereader = XMLPolyDataReader(FileName=facefile)
	facedisplay=Show(facereader,renderView)
	ColorBy(facedisplay, ('POINTS', 'VertexType'))
	vertexTypeLUT = GetColorTransferFunction('VertexType')
	vertexTypeLUT.ApplyPreset('X Ray', True)
	vertexTypeLUT.RGBPoints = [0.0, 1.,0.988235,0.968627, 1.0, 0.47451, 0.580392, 0.85098, 2.0, 0.25098,0.172549,0.690196]
	renderView.OrientationAxesVisibility = 0
	facedisplay.SetScalarBarVisibility(renderView, False)
	text = Text()
	text.Text = facefile.split("/")[-1]
	textdisplay = Show(text,renderView)
	textdisplay.FontSize = 10
	renderView.ResetCamera()
	Render()
	return renderView



path = os.getcwd()
facefiles = glob.glob('*fraction*.vtp')

for viewindex, facefile in enumerate(facefiles):
	paraviewvis(viewindex, os.path.join(path, facefile))

servermanager.SaveState("snapshots_paraview.pvsm")
