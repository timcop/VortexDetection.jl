using FLoops

function findvortices_planes_threaded(ψ; n_itr = 1)

    ϕ = angle.(ψ)

    Δϕx = zeros(size(ψ))
    Δϕy = zeros(size(ψ)) 
    Δϕz = zeros(size(ψ))

    @floop for i in eachindex(ψ[:, 1, 1])
        Δϕx[i, :, :] = findvortices_jumps_plane(ϕ[i, :, :])
    end

    @floop for j in eachindex(ψ[1, :, 1])
        Δϕy[:, j, :] = findvortices_jumps_plane(ϕ[:, j, :])
    end

    @floop for k in eachindex(ψ[1, 1, :])
        Δϕz[:, :, k] = findvortices_jumps_plane(ϕ[:, :, k])
    end

    return Δϕx, Δϕy, Δϕz
end

function link_graph_vorts(g; repeat_end_of_ring=true)
    g_temp = deepcopy(g)
    visited = Set()
    deg_1 = Set(findall(x -> x == 1, degree(g_temp)))

    # Link vortex lines by starting at degree 1 vortices (end points) and terminating at degree != 2 vertices
    vort_lines = []
    while length(deg_1) > 0
        current_vort = []
        vc = pop!(deg_1)
        
        push!(current_vort, vc)
        vc_neighbors = neighbors(g_temp, vc)

        # setdiff!(vc_neighbors, visited)
        while length(vc_neighbors) <= 2 # vc isn't a reconnection
            push!(visited, vc)
            setdiff!(vc_neighbors, visited)
            # print(length(vc_neighbors))
            if length(vc_neighbors) == 0
                break
            end
            vc = pop!(vc_neighbors)
            push!(current_vort, vc)
            vc_neighbors = neighbors(g_temp, vc)  
        end
        push!(vort_lines, current_vort)
        setdiff!(deg_1, visited)
    end

    # Only rings are left, now link rings that have reconnection points
    g_temp = deepcopy(g)
    deg_2g = Set(findall(x -> x > 2, degree(g_temp)))
    vort_loops = []

    while length(deg_2g) > 0
        current_vort = []
        vi = pop!(deg_2g)
        push!(current_vort, vi)
        push!(visited, vi)
        
        vi_neighbors = neighbors(g_temp, vi)

        # Now traverse a neighbor that is not a reconnection point and not visited
        vi = nothing
        for v in vi_neighbors
            if v ∉ visited
                if length(neighbors(g_temp, v)) <= 2 # This should only be 2 but just in case 
                    vi = v
                    break
                end
            end
        end

        while !(isnothing(vi))
            push!(current_vort, vi)
            push!(visited, vi)
            vi_neighbors = neighbors(g_temp, vi)
            setdiff!(vi_neighbors, visited)
            if length(vi_neighbors) == 1
                vi = pop!(vi_neighbors)
                # Check not reconnection point
                if length(neighbors(g_temp, vi)) > 2
                    push!(current_vort, vi)
                    push!(visited, vi)
                    vi = nothing
                end
            else
                vi = nothing
            end
        end

        # Now we've traversed the ring, remove visited vertices
        push!(vort_loops, current_vort)
        setdiff!(deg_2g, visited)
    end
        
    # Now what's left are rings that don't have reconnection points
    g_temp = deepcopy(g)

    deg_2 = Set(findall(x -> x == 2, degree(g_temp)))
    setdiff!(deg_2, visited)
    vort_rings = []
    while length(deg_2) > 0
        current_vort = []
        vi = pop!(deg_2)
        push!(current_vort, vi)
        push!(visited, vi)
        vi_neighbors = neighbors(g_temp, vi)
        setdiff!(vi_neighbors, visited)
        while length(vi_neighbors) != 0
            vc = pop!(vi_neighbors)
            push!(current_vort, vc)
            push!(visited, vc)
            vi_neighbors = neighbors(g_temp, vc)
            setdiff!(vi_neighbors, visited)
        end
        push!(vort_rings, current_vort)
        # put start point at end for ring plotting
        if repeat_end_of_ring
            push!(current_vort, current_vort[1])
        end
        deg_2 = setdiff!(deg_2, visited)
    end

    return vort_lines, vort_loops, vort_rings
    # return vort_lines, vort_loops
end

# Returns an array A where A[i] is a connected vortex filament built up of individual lines
function connect_vortex_ends(vorts, vort_lines, vort_loops, vort_rings, X)
    x = X[1]; y = X[2]; z = X[3];
    dx = x[2]-x[1]; dy = y[2]-y[1]; dz = z[2]-z[1];

    vorts_lines_loops = vcat(vort_lines, vort_loops) # v[i] are lines/loops of indicies referencing vorts array
    line_ends = []
    for i in eachindex(vorts_lines_loops)
        v = vorts_lines_loops[i]
        push!(line_ends, [vorts[v[1]], vorts[v[end]]])
    end

    connection_array = zeros(Bool, length(vorts_lines_loops), length(vorts_lines_loops))
    for i in eachindex(vorts_lines_loops)
        for j in 1:2
            vij = line_ends[i][j]
    
            is_x_vort = vij[1] % dx ≈ 0
            is_y_vort = vij[2] % dy ≈ 0
            is_z_vort = vij[3] % dz ≈ 0
    
            x_dist_1 = is_x_vort ? dx : dx/2
            y_dist_1 = is_y_vort ? dy : dy/2
            z_dist_1 = is_z_vort ? dz : dz/2
    
            for ii in eachindex(vorts_lines_loops)
                for jj in 1:2
                    vjj = line_ends[ii][jj]
    
                    x_dist = abs(vij[1] - vjj[1])
                    y_dist = abs(vij[2] - vjj[2])
                    z_dist = abs(vij[3] - vjj[3])
    
                    if (x_dist <= x_dist_1 || x_dist >= 16 - x_dist_1) && (y_dist <= y_dist_1 || y_dist >= 16 - y_dist_1) && (z_dist <= z_dist_1 || z_dist >= 16 - z_dist_1)
                        connection_array[i, ii] = true
                    end
                end
            end
        end
    end
    
    
    for i in axes(connection_array, 1)
        connection_array[i, i] = false
    end

    connection_graph = SimpleGraph(connection_array)

    connected_vorts_idxs = connected_components(connection_graph)

    connected_vorts = []
    for i in eachindex(connected_vorts_idxs)
        push!(connected_vorts, vorts_lines_loops[connected_vorts_idxs[i]])
    end

    for i in eachindex(vort_rings)
        push!(connected_vorts, [vort_rings[i]])
    end

    return connected_vorts
end

function findvortices_jumps_plane(phase)
    # phase = angle.(ψ);

    Δϕx, Δϕy = phase_jumps(phase,1),phase_jumps(phase,2)

    circshift!(phase,Δϕx,(0,1))
    Δϕx .-= phase; Δϕx .-= Δϕy
    circshift!(phase,Δϕy,(1,0))
    Δϕx .+= phase

    return abs.(Δϕx)
end



function vortex_coords(vorts_x, vorts_y, vorts_z, x, y, z)
    dx = x[2] - x[1]; dy = y[2] - y[1]; dz = z[2] - z[1];
    vorts_x_coords = [[x[v[1]], y[v[2]] - dy/2 , z[v[3]] - dz/2] for v in vorts_x]
    vorts_y_coords = [[x[v[1]] - dx/2, y[v[2]], z[v[3]] - dz/2] for v in vorts_y]
    vorts_z_coords = [[x[v[1]] - dx/2, y[v[2]] - dy/2, z[v[3]]] for v in vorts_z]
    vorts_coords = vcat(vorts_x_coords, vorts_y_coords, vorts_z_coords);
    return vorts_coords
end

function neighbour_vort(vorts, i, Δϕd, vorts_map_d, Δneighbour)
    v_idx = vorts[i]
    v_idx_i = v_idx .+ Δneighbour
    v_idx_i_mod = mod1.(v_idx_i, size(Δϕd))
    if Δϕd[v_idx_i_mod[1], v_idx_i_mod[2], v_idx_i_mod[3]] > 0.0
        if v_idx_i == v_idx_i_mod
            return true, false, vorts_map_d[v_idx_i_mod]
        else
            return true, true, vorts_map_d[v_idx_i_mod]
        end
    else
        return false, false, nothing
    end
end

function vortex_edge_list_threaded(Δϕx, Δϕy, Δϕz)
    vorts_x = Tuple.(findall(Δϕx .> 0.0)) .|> collect;
    vorts_y = Tuple.(findall(Δϕy .> 0.0)) .|> collect;
    vorts_z = Tuple.(findall(Δϕz .> 0.0)) .|> collect;

    vorts_x_len = length(vorts_x)
    vorts_y_len = length(vorts_y)
    vorts_z_len = length(vorts_z)

    vorts_x_map = Dict((vorts_x[i], i) for i in eachindex(vorts_x))
    vorts_y_map = Dict((vorts_y[i], i + vorts_x_len) for i in eachindex(vorts_y))
    vorts_z_map = Dict((vorts_z[i], i + vorts_x_len + vorts_y_len) for i in eachindex(vorts_z))

    edge_list_x = [[] for _ in 1:Threads.nthreads()]
    edge_list_periodic_x = [[] for _ in 1:Threads.nthreads()]
 
    @floop for i in eachindex(vorts_x)
        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕx, vorts_x_map, [-1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [0, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [0, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕx, vorts_x_map, [1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕy, vorts_y_map, [1, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_x, i, Δϕz, vorts_z_map, [1, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_x[Threads.threadid()], (i, v_idx)) : push!(edge_list_x[Threads.threadid()], (i, v_idx)))
    end

    edge_list_y = [[] for _ in 1:Threads.nthreads()]
    edge_list_periodic_y = [[] for _ in 1:Threads.nthreads()]

    @floop for i in eachindex(vorts_y)
        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕy, vorts_y_map, [0, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [-1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕy, vorts_y_map, [0, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕz, vorts_z_map, [0, 1, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [0, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_y, i, Δϕx, vorts_x_map, [-1, 1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_y[Threads.threadid()], (i + vorts_x_len, v_idx)) : push!(edge_list_y[Threads.threadid()], (i + vorts_x_len, v_idx)))

    end

    edge_list_z = [[] for _ in 1:Threads.nthreads()]
    edge_list_periodic_z = [[] for _ in 1:Threads.nthreads()]

    @floop for i in eachindex(vorts_z)
        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕz, vorts_z_map, [0, 0, -1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [-1, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, 0, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, -1, 0])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕz, vorts_z_map, [0, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [0, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕx, vorts_x_map, [-1, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, 0, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))

        is_neighbour, is_periodic, v_idx = neighbour_vort(vorts_z, i, Δϕy, vorts_y_map, [0, -1, 1])
        is_neighbour && (is_periodic ? push!(edge_list_periodic_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)) : push!(edge_list_z[Threads.threadid()], (i + vorts_x_len + vorts_y_len, v_idx)))
    end
    
    edge_list_x = reduce(vcat, edge_list_x)
    edge_list_y = reduce(vcat, edge_list_y)
    edge_list_z = reduce(vcat, edge_list_z)
    edge_list_periodic_x = reduce(vcat, edge_list_periodic_x)
    edge_list_periodic_y = reduce(vcat, edge_list_periodic_y)
    edge_list_periodic_z = reduce(vcat, edge_list_periodic_z)

    return vcat(edge_list_x, edge_list_y, edge_list_z), vcat(edge_list_periodic_x, edge_list_periodic_y, edge_list_periodic_z), vorts_x, vorts_y, vorts_z
end

function smooth_vortex_rings(vort_rings, max_itr=1)
    vort_rings_prev = deepcopy(vort_rings)
    vort_rings_next = deepcopy(vort_rings_prev)

    for itr in 1:max_itr
        for i in eachindex(vort_rings_prev)
            v_ring_length = length(vort_rings_prev[i])
            for j in eachindex(vort_rings_prev[i])
                vort_rings_next[i][j] = ((0.8/5) * sum(vort_rings_prev[i][mod1.(j-2:j+2, v_ring_length)])) + (0.2 * vort_rings_prev[i][j])
                # vort_rings_coords_smoothed_next[i][j] = ((0.8/3) * sum(vort_rings_coords_smoothed_prev[i][mod1.(j-1:j+1, v_ring_length)])) + (0.2 * vort_rings_coords_smoothed_prev[i][j]) # 3 point window
            end
        end
        vort_rings_prev = deepcopy(vort_rings_next)
    end

    return vort_rings_next
end




function full_algorithm(psi, x, y, z; n_itp = 1)
    if n_itp > 1
        range_x = collect(range(start=0, stop=length(x), length=length(x)*n_itp + 1))[2:end]
        range_y = collect(range(start=0, stop=length(y), length=length(y)*n_itp + 1))[2:end]
        range_z = collect(range(start=0, stop=length(z), length=length(z)*n_itp + 1))[2:end]

        x_itp = interpolate(x, BSpline(Linear()))
        y_itp = interpolate(y, BSpline(Linear()))
        z_itp = interpolate(z, BSpline(Linear()))

        psi_itp = interpolate(psi, BSpline(Linear()))

        x_etp = extrapolate(x_itp, Line())
        y_etp = extrapolate(y_itp, Line())
        z_etp = extrapolate(z_itp, Line())

        psi_etp = extrapolate(psi_itp, Line())

        x = x_etp(range_x); y = y_etp(range_y); z = z_etp(range_z)
        psi = psi_etp(range_x, range_y, range_z)
    end


    # @time Δϕx, Δϕy, Δϕz = find_vortices_circshift(psi, x, y, z);
    println("==============================")
    print("Finding vortex points on planes:")
    @time Δϕx, Δϕy, Δϕz = findvortices_planes_threaded(psi);

    print("Creating edge list:")
    @time edge_list, edge_list_periodic, vorts_x, vorts_y, vorts_z = vortex_edge_list_threaded(Δϕx, Δϕy, Δϕz);
    
    print("Creating graph:")
    @time g = SimpleGraph(Edge.(edge_list))

    print("Linking vortices:")
    @time vort_lines, vort_loops, vort_rings = link_graph_vorts(g, repeat_end_of_ring=true);

    print("Vort coords:")
    @time vorts_coords = vortex_coords(vorts_x, vorts_y, vorts_z, x, y, z)
    println("")
    println("Number of vortex points: $(length(vorts_coords))")
    println("Number of vortex lines: $(length(vort_lines))")
    println("Number of vortex loops: $(length(vort_loops))")
    println("Number of vortex rings: $(length(vort_rings))")
    println("")
    println("Number of vortices: $(length(vort_lines) + length(vort_loops) + length(vort_rings))")
    # @time connected_vorts = connect_vortex_ends(vorts_coords, vort_lines, vort_loops, vort_rings, X);
    println("==============================")
    return g, vort_lines, vort_loops, vort_rings, vorts_coords
end

function vortex_moving_avg(vorts_coords, vort_lines, vort_loops, vort_rings; n_itr = 10)
    vort_lines = [vorts_coords[line] for line in vort_lines]
    v_lines_prev = copy(vort_lines)
    vort_lines_avg = []
    for itr in 1:n_itr

        vort_lines_avg = []
        for v_line in v_lines_prev
            if length(v_line) > 5
                v_line_avg = []
                # Handing first 2 points
                
                mi = v_line[1]
                push!(v_line_avg, mi)
        
        
                mi = (0.8/3)*(v_line[1] + v_line[2] + v_line[3]) + 0.2*v_line[2]
                push!(v_line_avg, mi)
        
        
                for i in 3:length(v_line)-2
                    mi = (0.8/5)*(v_line[i-2] + v_line[i-1] + v_line[i] + v_line[i+1] + v_line[i+2]) + 0.2*v_line[i]
                    push!(v_line_avg, mi)
                end
        
                mi = (0.8/3)*(v_line[end-2] + v_line[end-1] + v_line[end]) + 0.2*v_line[end-1]
                push!(v_line_avg, mi)
        
                mi = v_line[end]
                push!(v_line_avg, mi)
                push!(vort_lines_avg, v_line_avg)
            else
                v_line_avg = v_line
                push!(vort_lines_avg, v_line_avg)
            end
        end
        v_lines_prev = copy(vort_lines_avg)
    end
    
    vort_loops = [vorts_coords[loop] for loop in vort_loops]
    v_loops_prev = copy(vort_loops)
    vort_loops_avg = []
    for itr in 1:n_itr

        vort_loops_avg = []
        for v_loop in v_loops_prev
            if length(v_loop) > 5
                v_loop_avg = []
                # Handing first 2 points
                
                mi = v_loop[1]
                push!(v_loop_avg, mi)
        
        
                mi = (0.8/3)*(v_loop[1] + v_loop[2] + v_loop[3]) + 0.2*v_loop[2]
                push!(v_loop_avg, mi)
        
        
                for i in 3:length(v_loop)-2
                    mi = (0.8/5)*(v_loop[i-2] + v_loop[i-1] + v_loop[i] + v_loop[i+1] + v_loop[i+2]) + 0.2*v_loop[i]
                    push!(v_loop_avg, mi)
                end
        
                mi = (0.8/3)*(v_loop[end-2] + v_loop[end-1] + v_loop[end]) + 0.2*v_loop[end-1]
                push!(v_loop_avg, mi)
        
                mi = v_loop[end]
                push!(v_loop_avg, mi)
                push!(vort_loops_avg, v_loop_avg)
            else
                v_loop_avg = v_loop
                push!(vort_loops_avg, v_loop_avg)
            end
        end
        v_loops_prev = copy(vort_loops_avg)
    end

    vort_rings = [vorts_coords[ring] for ring in vort_rings]
    v_rings_prev = copy(vort_rings)
    vort_rings_avg = []
    # vort_rings_next[i][j] = ((0.8/5) * sum(vort_rings_prev[i][mod1.(j-2:j+2, v_ring_length)])) + (0.2 * vort_rings_prev[i][j])
    for itr in 1:n_itr

        vort_rings_avg = []
        for v_ring in v_rings_prev
            if length(v_ring) > 5
                v_ring_avg = []
                # Handing first 2 points
                
                
                
                v_ring_length = length(v_ring)
                for i in eachindex(v_ring)
                    mi = (0.8/5) * sum(v_ring[mod1.(i-2:i+2, v_ring_length)]) + (0.2 * v_ring[i])
                    push!(v_ring_avg, mi)
                end
                
                push!(vort_rings_avg, v_ring_avg)
            else
                v_ring_avg = v_ring
                push!(vort_rings_avg, v_ring_avg)
            end
        end
        v_rings_prev = copy(vort_rings_avg)
    end

   
    return vort_lines_avg, vort_loops_avg, vort_rings_avg
end

includet("ccma.jl")
function vortex_ccma(ccma, vort_lines_mat, vort_loops_mat, vort_rings_mat; distrib="hanning", w_ma=2, w_cc=3, itr=1)
    # Fix single vortex error
    ccma_inst = ccma.CCMA(w_ma, w_cc, distrib=distrib)
    # vort_lines_ccma = []
    # for v_l in vort_lines_mat
    #     try
    #         v_l_ccma = ccma_inst.filter(v_l, mode="fill_boundary")
    #         push!(vort_lines_ccma, v_l_ccma)
    #     catch err
    #         println(err)
    #         push!(vort_lines_ccma, v_l)
    #     end
    # end

    # vort_loops_ccma = []
    # for v_l in vort_loops_mat
    #     try
    #         v_l_ccma = ccma_inst.filter(v_l)
    #         push!(vort_loops_ccma, v_l_ccma)
    #     catch
    #         push!(vort_loops_ccma, v_l)
    #     end
    # end

    # vort_rings_ccma = []
    # for v_l in vort_rings_mat
    #     try
    #         v_l_ccma = ccma_inst.filter(v_l, mode="wrapping")
    #         push!(vort_rings_ccma, v_l_ccma)
    #     catch
    #         push!(vort_rings_ccma, v_l)
    #     end
    # end
    # return vort_lines_ccma, vort_loops_ccma, vort_rings_ccma

    vort_lines_ccma = Vector{Any}(undef, length(vort_lines_mat))
    for (i, v_l) in enumerate(vort_lines_mat)
        try
            v_l_ccma = ccma_inst.filter(v_l, mode="fill_boundary")

            vort_lines_ccma[i] = v_l_ccma
        catch err
            println("Error at vort_lines: $i")
            println("")
            # println(err)
            # showerror(stdout, err, catch_backtrace())
            vort_lines_ccma[i] = v_l
        end
    end

    # vort_loops_ccma = []
    vort_loops_ccma = Vector{Any}(undef, length(vort_loops_mat))
    for (i, v_l) in enumerate(vort_loops_mat)
        try
            v_l_ccma = ccma_inst.filter(v_l)
            vort_loops_ccma[i] = v_l_ccma
        catch
            println("Error at vort_loops: $i")
            println("")
            println(err)
            # push!(vort_loops_ccma, v_l)
            vort_loops_ccma[i] = v_l
        end
    end

    # vort_rings_ccma = []
    vort_rings_ccma = Vector{Any}(undef, length(vort_rings_mat))
    for (i, v_l) in enumerate(vort_rings_mat)
        try
            v_l_ccma = ccma_inst.filter(v_l, mode="wrapping")
            # push!(vort_rings_ccma, v_l_ccma)
            vort_rings_ccma[i] = v_l_ccma
        catch
            println("Error at vort_rings: $i")
            println("")
            println(err)
            # push!(vort_rings_ccma, v_l)
            vort_rings_ccma[i] = v_l
        end
    end
    return vort_lines_ccma, vort_loops_ccma, vort_rings_ccma
end

function vortex_ccma_j(vort_lines_mat, vort_loops_mat, vort_rings_mat; distrib="hanning", w_ma=2, w_cc=3, itr=1, cc_mode=true)
    ccma_inst = CCMA(w_ma = w_ma, w_cc = w_cc, distrib=distrib)
    vort_lines_ccma = Vector{Any}(undef, length(vort_lines_mat))

    @floop for (i, v_l) in enumerate(vort_lines_mat)
        try
            v_l_ccma = filter(ccma_inst, v_l, mode="fill_boundary", cc_mode=cc_mode)
            # push!(vort_lines_ccma, v_l_ccma)
            vort_lines_ccma[i] = v_l_ccma
        catch err
            println("Error at vort_lines: $i")
        #     println("")
        #     # println(err)
        #     showerror(stdout, err, catch_backtrace())
        #     # push!(vort_lines_ccma, v_l)
            vort_lines_ccma[i] = v_l
        end
    end

    vort_loops_ccma = Vector{Any}(undef, length(vort_loops_mat))
    @floop for (i, v_l) in enumerate(vort_loops_mat)
        try
            v_l_ccma = filter(ccma_inst, v_l, mode="fill_boundary", cc_mode=cc_mode)
            # push!(vort_loops_ccma, v_l_ccma)
            vort_loops_ccma[i] = v_l_ccma
        catch
            println("Error at vort_loops: $i")
            println("")
            # println(err)
            # push!(vort_loops_ccma, v_l)
            vort_loops_ccma[i] = v_l
        end
    end

    # vort_rings_ccma = []
    vort_rings_ccma = Vector{Any}(undef, length(vort_rings_mat))
    @floop for (i, v_l) in enumerate(vort_rings_mat)
        try
            v_l_ccma = filter(ccma_inst, v_l, mode="wrapping", cc_mode=cc_mode)
            # push!(vort_rings_ccma, v_l_ccma)
            vort_rings_ccma[i] = v_l_ccma
        catch
            println("Error at vort_rings: $i")
            println("")
            # println(err)
            # push!(vort_rings_ccma, v_l)
            vort_rings_ccma[i] = v_l
        end
    end
    return vort_lines_ccma, vort_loops_ccma, vort_rings_ccma
end