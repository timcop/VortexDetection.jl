## Imports 
using Revise
using FourierGPE
using VortexDistributions
using JLD2
using Graphs, GraphPlot, GraphMakie
using BenchmarkTools
using SparseArrays
using StaticArrays
using FLoops
using PyCall
using Interpolations

## Include utils functions
includet("../core/utils_plots.jl") # includet is short hand for include and then track the changes with revise
includet("../core/utils.jl")
includet("../core/ccma.jl")


# Load data, if it doesn't exist run the src/sims/simulation_64.jl file to generate it
@load "src/data/sol_64.jld2" sol;
@load "src/data/sim_64.jld2" sim;


## Params
t = 20 # Choose a t to analyse, smaller t will be more turbulent (t = 0-10) and larger t will be more relaxed (t > 10)
psi = sol(t);
X = sim.X;
x = X[1]; y = X[2]; z = X[3];
x = round.(x, digits=3); y = round.(y, digits=3); z = round.(z, digits=3); # Rounding to avoid floating point issues
dx = x[2] - x[1]; dy = y[2] - y[1]; dz = z[2] - z[1];

## Plot the isosurface before detection to look at the state
fig, lscene = plot_iso(
    psi, 
    X, 
    show_axis=true, 
    isovalue=0.5, isorange=0.15, # Adjust isorange to change the thickness of the isosurface
    is_128=false, # Hack for bug with 128 grid
    visible=true # Set to false if you don't want to see the isosurface but want to keep the scene for later
)
fig

##
n_interpolate = 4 # Number of times to interpolate the wavefunction before vortex detection, higher = better resolution but slower. 1 = no interpolation
@time g, vort_lines, vort_loops, vort_rings, vorts_coords = full_algorithm(psi, x, y, z, n_itp = n_interpolate);
graphplot!(g, layout = (adj) -> vorts_coords, edge_width=1, node_size=2, show_axis=false)

## Smoothing
vort_lines_coords = [reduce(vcat, transpose.(vorts_coords[v_line])) for v_line in vort_lines];
vort_loops_coords = [reduce(vcat, transpose.(vorts_coords[v_loop])) for v_loop in vort_loops];
vort_rings_coords = [reduce(vcat, transpose.(vorts_coords[v_ring])) for v_ring in vort_rings];

n_itr = 10

@time vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j = vortex_ccma_j(vort_lines_coords, vort_loops_coords, vort_rings_coords, w_ma = 3, w_cc = 2, distrib = "hanning");

# Optionally perform smoothing multiple times to get smoother lines
for i in 1:n_itr
    vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j = vortex_ccma_j(vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j, w_ma = 3, w_cc = 2, distrib = "hanning");
end


# Colour code based on rings/loops/lines
plot_unconnected_vorts_mat(vort_lines_ccma_j, vort_loops_ccma_j, vort_rings_ccma_j, linewidth=3, mono=true, repeat_ring_ends=true, color=:blue)