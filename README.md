# Interpolation Workflow

## Building
### Supported Operating Systems

* Windows (Visual Studio)

* Linux (g++ / clang)

* macOS (XCode)

### System Requirements

* C++17 Compiler with OpenMP support

### Dependencies

* [ViennaPS](https://github.com/ViennaTools/ViennaPS)
* LAPACK
* ([fmt](https://github.com/fmtlib/fmt))

## Components

* __InterpolationWorkflow__: The main executable of this example project.
* __CreateData__: The executable used to create the _data.csv_ file.
* __GeometryReconstruction__: A simple example which first runs a physical deposition on a trench geometry. Afterwards the dimensions of the final geometry are extracted and subsequently used to geometrically reconstruct the trench geometry.
* __TestSpline__: A small example program showcasing the Spline interpolation functionality in 1D and ND.
