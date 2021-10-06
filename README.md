# SurfMeshSmooth
Surface triangle mesh smoothing, an extension of basic laplacian smoothing that takes vertex normals into account in order to preserve surface curvature during smoothing.

Might not properly handle sharp concave corners on 1D mesh boundary, if any.

Runs parallel on OpenMP or CUDA if available.

Performance (on a quad core Intel Xeon E5-1620 v4 + NVIDIA Quadro M4000): 1000 smoothing iterations on a 1 million vertex mesh with ~2 million triangles in about ~2 seconds with CUDA, and about ~27 sec with OpenMP.

![1a](https://github.com/codearxiv/SurfMeshSmooth/blob/main/images/smooth0.PNG)
