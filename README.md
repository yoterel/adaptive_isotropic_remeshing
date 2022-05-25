This is an unofficial implementation of the paper:

[Adaptive Remeshing for Real-Time Mesh Deformation](https://diglib.eg.org/handle/10.2312/conf.EG2013.short.029-032)

The compilation of this code (C++) produces a python library as well (using numpyeigen).

<table>
  <tr>
    <td>Nefratiti compressed from ~100k faces to ~36k using this implementation. Notice the mesh is much more isotropic.</td>
  </tr>
 <tr>
<td align="center">
<img src="https://github.com/yoterel/adaptive_isotropic_remeshing/blob/master/resource/nefratiti.png" alt="0" width = 400px>
</td>
</tr>
</table>

Based mostly on a previous implementation from [here](https://github.com/sgsellan/botsch-kobbelt-remesher-libigl), which implements [A Remeshing Approach to Multiresolution Modeling](https://dl.acm.org/doi/10.1145/1057432.1057457).

The key contribution of this paper is introducing a "sizing field" which is a scalar-per-vertex dependant on the local curvature of the mesh, and dictates the target edge length (for which the vertex is one of its endpoints) - high curvature areas recieve a smaller sizing field).

This implementation degenerates into the original one with the flag -na, forcing a constant sizing field upon all the mesh.

## Installation

1. Clone the repository with all submodules (libigl, numpyeigen)

`git clone --recurse-submodules https://github.com/yoterel/adaptive_isotropic_remeshing.git`

2. The usual cmake & make:

```
cd adaptive_isotropic_remeshing
mkdir build
cd build
cmake ..
make
```

3. launch:

`./adaptive_remesh`

or for python (requires libigl, but you can use any other library to load some mesh):

```
import igl
import numpy as np
from py_ada_remesh import adaptive_remesh_botsch
v, f = igl.read_triangle_mesh("./../resource/nefratiti_orig.ply")
new_v, new_f = adaptive_remesh_botsch(v, np.int32(f), 0.0008, 5, False, True)
igl.write_obj("nefratiti_compressed.obj", new_v, new_f)
```

where v is a Nx3 numpy array, f is a Mx3 numpy array, and for the other parameters see below.

## Options

-i  | Number of iterations to run remeshing

-e  | Epsilon (for adaptive) / Target Edge Length (For non adaptive), controls the resolution of output

-p  | Project result vertices of each iteration on original mesh

-na | Not Adaptive, will use original botsch & kobbelt remeshing algorithm
