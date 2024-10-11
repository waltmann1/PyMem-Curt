#create a system 
system = mb.System()
#read the mesh
vertex_file = 'vertices_R1.0.inp'
face_file = 'faces_R1.0.inp'
system.read_mesh_from_files(files = {'vertices': vertex_file, 'faces': face_file})