# Interpolation Workflow

## Building
### Supported Operating Systems

* Windows (Visual Studio)
* Linux (g++ / clang)
* macOS (XCode)

### System Requirements

* C++17 Compiler with OpenMP support

### Dependencies
* Required:
    * [ViennaLS](https://github.com/ViennaTools/ViennaLS)
    * LAPACK
* Optional:
  * [ViennaPS](https://github.com/ViennaTools/ViennaPS)

## Components

* __InterpolationWorkflow__: The main executable of this example project.

* If _ViennaPS_ installation is available:
  * __CreateData__: The executable used to create the _data.csv_ file.
  * __GeometryReconstruction__: A simple example which first runs a physical deposition on a trench geometry. Afterwards the dimensions of the final geometry are extracted and subsequently used to geometrically reconstruct the trench geometry.

## Building

```bash
git clone https://github.com/yozoon/ViennaPS-InterpolationWorkflow.git
cd ViennaPS-InterpolationWorkflow/
mkdir build && cd build
cmake ..
make
```

## Usage

The _data.csv_ file contains data sampled from a regular grid of the parameter space with the following values:

| Parameter | Values |
| --- | --- |
| _taperAngle_ | -10°,-5°,0°, 5°,10°,15° |
| _stickingProbability:_ | 1.0, 0.7, 0.4, 0.1 |
| _processTime:_ | 0, 1, 2, 3, 4, 5 |

When the _InterpolationWorkflow_ executable is invoked, two _.vtp_ files are generated: 1. _IW_initial.vtp_ which is the trench geometry before the deposition process and 2. _IW_interpolated.vtp_ which is the geometry after the deposition process, created through cubic interpolation of geometric features based on the provided data samples (_data.csv_ file).

```bash
./InterpolationWorkflow data.csv config.txt
```

You can update the values in the _config.txt_-file in your build directory to match your desired configuration.
