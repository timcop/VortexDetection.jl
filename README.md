# VortexDetection3D.jl

A Julia package for detecting and analyzing vortices in 3D systems.

## Requirements

This package has been tested and verified to work with:
- Julia 1.10.4

Earlier versions may work but are not officially supported.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/timcop/VortexDetection3D.jl.git
cd VortexDetection3D.jl
```

2. Start Julia and activate the project:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will install all required dependencies specified in the `Project.toml` file, including:
- BenchmarkTools
- ColorSchemes
- DSP
- FourierGPE
- GLMakie
- GraphMakie
- JLD2
- and more

## Running Simulations

The package includes simulation capabilities in the `sims` directory. If the simulation data is not present in the `data` directory, you can generate it by running the simulation notebook:

1. Open VS Code
2. Navigate to `sims/simulation_64.jl` and open it
3. Click "Open With..." in the top right corner and select "Julia Notebook"
4. Run the notebook cells sequentially

This will create the simulation data files (`sim_64.jld2` and `sol_64.jld2`) in the `data` directory.

## Example Usage

An example demonstrating the vortex detection algorithm is provided in `examples/quench_64_3d_algo.jl`. To run the example:

1. Open VS Code
2. Navigate to `examples/quench_64_3d_algo.jl` and open it
3. Click "Open With..." in the top right corner and select "Julia Notebook"
4. Run the notebook cells sequentially

This example demonstrates:
- Loading simulation data
- Processing 3D vortex structures
- Visualizing the results using GLMakie and GraphMakie

## Project Structure

- `src/core/`: Core functionality
  - `ccma.jl`: Curvature corrected moving average smoothing functions 
  - `utils.jl`: Core 3d detection algorithm functions
  - `utils_plots.jl`: Plotting functions
- `src/data/`: Simulation data storage
- `src/examples/`: Example scripts
- `src/sims/`: Simulation scripts

## Dependencies

The project uses several key Julia packages:
- GLMakie and GraphMakie for visualization
- FourierGPE for GPE calculations
- JLD2 for data storage
- DSP for signal processing
- Various other utilities for graphs, interpolation, and distributions

## Contributing

Feel free to open issues or submit pull requests with improvements or bug fixes.
