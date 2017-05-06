
# Mesh Segmentation

A tool for cutting a closed surface (represented as triangular mesh) at a user defined position. It is an implementation of the paper *Mesh Decomposition with Cross-Boundary Brushes* by Youyi Zheng and Chiew-Lan Tai (Eurographics 2010).

## Libraries

Qt, OpenGL, SuiteSparse (for solving linear system, replacable by Eigen), meshlib (courtesy of David Gu)

## Algorithm

A quasi-harmonic scalar field is computed on the mesh, with boundary conditions defined by the user stroke. Then isolines that intersect the user stroke are found. The optimal one of them, with respect to the metric of local shape information, is selected as the cut.

## Source File descriptions

### MeshViewer.h/.cpp

The framework. Contains methods for rendering the scene (draw the mesh in points mode, wireframe mode, etc.), and methods responding to toolbar button pushes (show the harmonic field, show the isolines, show the final cut).

### isoline.h/.cpp

The core of the algorithm. Defines how an isoline is represented and found. Implements the meandering triangles algorithm (variation of marching cube).

### MainWindow.h/.cpp, checkableAction.h/.cpp, RenderOptionsWindow.h/.cpp

Defines toolbar buttons and display options.

### functions.h/.cpp

A collection of some commonly used methods.
