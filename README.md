This is an unofficial implementation of the paper:

[Adaptive Remeshing for Real-Time Mesh Deformation](https://diglib.eg.org/handle/10.2312/conf.EG2013.short.029-032)

<table>
  <tr>
    <td>Nefratiti compressed from ~1 million faces to ~360k using this implementation</td>
  </tr>
 <tr>
<td align="center">
<img src="https://github.com/yoterel/adaptive_isotropic_remeshing/blob/master/resource/nefratiti.png" alt="0" width = 400px>
</td>
</tr>
</table>

Based mostly on a previous implementation from [here](https://github.com/sgsellan/botsch-kobbelt-remesher-libigl), which implements [A Remeshing Approach to Multiresolution Modeling](https://dl.acm.org/doi/10.1145/1057432.1057457).

The paper introduces a sizing field which is a scalar dependant on the local curvature of the mesh, and dictates how the target edge length changes (smaller for higher curvature areas).

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

## Options

-i  | Number of iterations to run remeshing

-e  | Epsilon (for adaptive) / Target Edge Length (For non adaptive), controls the resolution of output

-p  | Project result vertices of each iteration on original mesh

-na | Not Adaptive, will use original botsch & kobbelt remeshing algorithm
