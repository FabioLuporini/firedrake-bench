2D_MESH = box_cylinder_2d
2D_SCALES = 1.0 0.707 0.5 0.354 0.25
2D_MESHES := $(foreach scale, $(2D_SCALES), $(2D_MESH)_$(scale).msh)

3D_MESH = box_cylinder_3d
3D_SCALES = 1.0 0.8 0.62 0.48 0.38 0.3
3D_MESHES := $(foreach scale, $(3D_SCALES), $(3D_MESH)_$(scale).msh)

# SCALES = 1.0 0.707 0.5 0.354 0.25 0.177 0.125
meshes: 2d 3d

2d: $(2D_MESHES)

3d: $(3D_MESHES)

$(2D_MESH)_%.msh: $(2D_MESH).geo
	gmsh -2 -clscale $* -o $@ $<

$(3D_MESH)_%.msh: $(3D_MESH).geo
	gmsh -3 -clscale $* -o $@ $<
