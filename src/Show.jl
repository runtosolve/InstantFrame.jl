module Show

using GLMakie, InstantFrame, LinearAlgebra


function deformed_shape(node, element, u, scale)

    #Plot solution.
    X = [node.coordinates[i][1] for i in eachindex(node.list)]
    Y = [node.coordinates[i][2] for i in eachindex(node.list)]
    Z = [node.coordinates[i][3] for i in eachindex(node.list)]

    ΔX = u[1:6:end]
    ΔY = u[2:6:end]
    ΔZ = u[3:6:end]
    θX = u[4:6:end]
    θY = u[5:6:end]
    θZ = u[6:6:end]

    X_range = maximum(X) - minimum(X)
    Y_range = maximum(Y) - minimum(Y)
    Z_range = maximum(Z) - minimum(Z)

    f = Figure()
    ax = Axis3(f[1,1])
    ax.aspect = (1.0, Y_range/X_range, Z_range/X_range)

    for i in eachindex(element.list)

        index_i = findfirst(x->x==element.start_node[i], node.list)
        index_j = findfirst(x->x==element.end_node[i], node.list)

        scatterlines!([X[index_i], X[index_j]], [Y[index_i], Y[index_j]], [Z[index_i], Z[index_j]], markersize = 5)
        scatterlines!([X[index_i] + scale[1] * ΔX[index_i], X[index_j] + scale[1] * ΔX[index_j]], [Y[index_i] + scale[2] * ΔY[index_i], Y[index_j] + scale[2] * ΔY[index_j]], [Z[index_i] + scale[3] * ΔZ[index_i], Z[index_j] + scale[3] * ΔZ[index_j]], markersize = 5,  linestyle=:dash)

    end

    f

end



function beam_shape_function(q1, q2, q3, q4, L, x)
    
    a0 = q1
    a1 = q2
    a2 = 1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
    a3 = 1/L^3*(2*q1+q2*L-2*q3+q4*L)
    w = a0 + a1*x + a2*x^2 + a3*x^3

    return w

end


function column_shape_function(a1, a2, L, x)

    a_x = a1 + (a2-a1)/L * x

    return a_x

end

function get_element_deformed_shape(u_local_e, L, Γ, n)

    x = range(0.0, L, n)

    #local x-y plane
    dof = [2, 6, 8, 12]
    q1, q2, q3, q4 = u_local_e[dof]
    w_xy = beam_shape_function.(q1, q2, q3, q4, L, x)

    #local x-z plane
    dof = [3, 5, 9, 11]
    q1, q2, q3, q4 = u_local_e[dof]
    w_xz = beam_shape_function.(q1, q2, q3, q4, L, x)

    #axial deformation along local x-axis
    dof = [1, 7]
    a1, a2 = u_local_e[dof]
    a_x = column_shape_function.(a1, a2, L, x)

    #local element deformation
    δ = [zeros(Float64, 6) for i in eachindex(x)]

    for i in eachindex(x)

        δ[i][1] = a_x[i]
        δ[i][2] = w_xy[i]
        δ[i][3] = w_xz[i]

    end

    #global element deformation
    Δ = [Γ'[1:6,1:6] * δ[i] for i in eachindex(x)]

    return δ, Δ, x

end



function discretized_element_global_coords(node_i_coords, Γ, x)

    local_element_discretized_coords = [zeros(Float64, 3) for i in eachindex(x)]

    [local_element_discretized_coords[i][1] = x[i] for i in eachindex(x)]

    global_element_discretized_coords = [Γ'[1:3, 1:3] * (local_element_discretized_coords[i]) .+  node_i_coords for i in eachindex(x)]

    return global_element_discretized_coords

end


function get_display_coords_element(u_local_e, node_i_coords, L, Γ, n)

    δe, Δe, x = get_element_deformed_shape(u_local_e, L, Γ, n)
    discretized_element_global_coords = Show.discretized_element_global_coords(node_i_coords, Γ, x)

    return discretized_element_global_coords, Δe

end


function get_display_coords(element, node, properties, u_local_e, n)

    display_coords = Array{Array{Array{Float64, 1}, 1}}(undef, length(properties.L))
    display_Δ = Array{Array{Array{Float64, 1}}}(undef, length(properties.L))

    for i=1:length(properties.L)

        index = findfirst(num->num==element.nodes[i][1], node.numbers)
        node_i_coords = node.coordinates[index]

        display_coords[i], display_Δ[i] = get_display_coords_element(u_local_e[i], node_i_coords, properties.L[i], properties.Γ[i], n[i])

        # get_display_coords_element(u_local_e[i], node_i_coords, properties.L[i], properties.Γ[i], n[i])


    end

    return display_coords, display_Δ

end



function show_element_deformed_shape!(ax, element_XYZ, Δ, scale, color)

    X = [element_XYZ[i][1] for i in eachindex(element_XYZ)]
    Y = [element_XYZ[i][2] for i in eachindex(element_XYZ)]
    Z = [element_XYZ[i][3] for i in eachindex(element_XYZ)]

    ΔX = [Δ[i][1] for i in eachindex(element_XYZ)]
    ΔY = [Δ[i][2] for i in eachindex(element_XYZ)]
    ΔZ = [Δ[i][3] for i in eachindex(element_XYZ)]


    for i=1:(length(X)-1)

        # scatterlines!(ax, [X[i], X[i+1]], [Y[i], Y[i+1]], [Z[i], Z[i+1]], markersize = 5)
        lines!(ax, [X[i] + scale[1] * ΔX[i], X[i+1] + scale[1] * ΔX[i+1]], [Y[i] + scale[2] * ΔY[i], Y[i+1] + scale[2] * ΔY[i+1]], [Z[i] + scale[3] * ΔZ[i], Z[i+1] + scale[3] * ΔZ[i+1]], linestyle=:dash, color=color, linewidth=2)

    end

    # return ax

end


function define_global_element_displacements(u, global_dof, element, element_connections)

    u_global_e = [zeros(Float64, 12) for i in eachindex(global_dof)]

    for i in eachindex(global_dof)

        u_global_e[i] = u[global_dof[i]]

        #update element deformations to consider partially restrained connections
        index = findfirst(num->num==element.numbers[i], element_connections.elements)

        if !isnothing(index)

            u_global_e[i][[5, 6, 11, 12]] .= element_connections.displacements[index][[5, 6, 11, 12]]

        end
        
    end

    return u_global_e

end

function show_model_deformed_shape(nodal_displacements, element_connections, element, node, properties, scale)

    n = vec(ones(Int64, size(properties.L, 1)) * 11)

    u = Array{Float64}(undef, 0)

    for i in eachindex(nodal_displacements)

        u = [u; nodal_displacements[i]]

    end

    u_global_e = define_global_element_displacements(u, properties.global_dof, element, element_connections)

    u_local_e = [properties.Γ[i]*u_global_e[i] for i in eachindex(u_global_e)]

    element_display_coords, element_display_Δ = InstantFrame.UI.get_display_coords(element, node, properties, u_local_e, n)

    figure = Figure()
    ax = Axis3(figure[1,1])

    for i in eachindex(element_display_coords)
        ax = show_element_deformed_shape(element_display_coords[i], element_display_Δ[i], scale, ax)
    end

    # ax.azimuth[] = 3π/2
    # ax.elevation[] = π/2
    # ax.aspect = (1.0, Y_range/X_range, Z_range/X_range)

    return figure

end


######new stuff



function define_element_nodal_start_end_coordinates(element, node)

    element_nodal_coords = Vector{NamedTuple{(:start_node, :end_node), Tuple{Tuple{Float64, Float64, Float64}, Tuple{Float64, Float64, Float64}}}}(undef, length(element.numbers))

    for i in eachindex(element.numbers)

        start_node_index = findfirst(num->num==element.nodes[i][1], node.numbers)
        end_node_index = findfirst(num->num==element.nodes[i][2], node.numbers)

        element_nodal_coords[i] = (start_node=node.coordinates[start_node_index], end_node=node.coordinates[end_node_index])

    end

    return element_nodal_coords

end



function get_node_XYZ(node)

    X = [node.coordinates[i][1] for i in eachindex(node.coordinates)]
    Y = [node.coordinates[i][2] for i in eachindex(node.coordinates)]
    Z = [node.coordinates[i][3] for i in eachindex(node.coordinates)]

    return X, Y, Z

end

function elements!(ax, element_nodal_coords, color)

    for i in eachindex(element_nodal_coords)
        lines!(ax, [element_nodal_coords[i].start_node[1], element_nodal_coords[i].end_node[1]], [element_nodal_coords[i].start_node[2], element_nodal_coords[i].end_node[2]], [element_nodal_coords[i].start_node[3], element_nodal_coords[i].end_node[3]], color = color)

    end

end

function nodes!(ax, X, Y, Z, markersize, color)

    scatter!(ax, X, Y, Z, markersize=markersize, color = color)

end


function point_loads!(ax, point_load, node, arrow_scale, arrow_head_scale, unit_arrow_head_size, arrowcolor, linecolor, linewidth)

    if !isnothing(point_load.nodes)

        for i in eachindex(point_load.nodes)

            index = findfirst(num->num==point_load.nodes[i], node.numbers)
            tail_location = node.coordinates[index]

            FX = point_load.magnitudes.FX[i] 
            FY = point_load.magnitudes.FY[i] 
            FZ = point_load.magnitudes.FZ[i] 

            unit_arrow_vector = [FX, FY, FZ] / norm([FX, FY, FZ])
        
            arrow_vector = unit_arrow_vector .* arrow_scale

            arrow_head_size = unit_arrow_head_size * arrow_head_scale

            arrow_size_vector = arrow_head_scale * unit_arrow_vector

            arrows!(ax, [tail_location[1]-arrow_vector[1] - arrow_size_vector[1]], [tail_location[2]-arrow_vector[2] - arrow_size_vector[2]], [tail_location[3]-arrow_vector[3]-arrow_size_vector[3]], [arrow_vector[1]], [arrow_vector[2]], [arrow_vector[3]], arrowsize = arrow_head_size, linewidth=linewidth,
                arrowcolor = arrowcolor, linecolor = linecolor)

        end

    end

end


function uniform_loads!(ax, uniform_load, element, node, unit_arrow_head_size, arrow_head_scale, arrow_scale, linewidth, arrowcolor, linecolor)

    for i in eachindex(uniform_load.elements)

        index = findfirst(num->num==uniform_load.elements[i], element.numbers)
        element_nodes = element.nodes[index]

        qX = uniform_load.magnitudes.qX[i] 
        qY = uniform_load.magnitudes.qY[i] 
        qZ = uniform_load.magnitudes.qZ[i] 

        #start node

        index = findfirst(num->num==element_nodes[1], node.numbers)
        tail_location = node.coordinates[index]

        if qX != 0.0
            unit_arrow_vector = [qX, 0.0, 0.0] / norm([qX, 0.0, 0.0])
            arrow_vector = unit_arrow_vector .* arrow_scale
            define_load_arrow!(ax, unit_arrow_head_size, arrow_head_scale, unit_arrow_vector, arrow_vector, tail_location, linewidth, arrowcolor, linecolor)
        end

        if qY != 0.0
            unit_arrow_vector = [0.0, qY, 0.0] / norm([0.0, qY, 0.0])
            arrow_vector = unit_arrow_vector .* arrow_scale
            define_load_arrow!(ax, unit_arrow_head_size, arrow_head_scale, unit_arrow_vector, arrow_vector, tail_location, linewidth, arrowcolor, linecolor)
        end

        if qZ != 0.0
            unit_arrow_vector = [0.0, 0.0, qZ] / norm([0.0, 0.0, qZ])
            arrow_vector = unit_arrow_vector .* arrow_scale
            define_load_arrow!(ax, unit_arrow_head_size, arrow_head_scale, unit_arrow_vector, arrow_vector, tail_location, linewidth, arrowcolor, linecolor)
        end


        #end node 

        index = findfirst(num->num==element_nodes[2], node.numbers)
        tail_location = node.coordinates[index]

        if qX != 0.0
            unit_arrow_vector = [qX, 0.0, 0.0] / norm([qX, 0.0, 0.0])
            arrow_vector = unit_arrow_vector .* arrow_scale
            define_load_arrow!(ax, unit_arrow_head_size, arrow_head_scale, unit_arrow_vector, arrow_vector, tail_location, linewidth, arrowcolor, linecolor)
        end

        if qY != 0.0
            unit_arrow_vector = [0.0, qY, 0.0] / norm([0.0, qY, 0.0])
            arrow_vector = unit_arrow_vector .* arrow_scale
            define_load_arrow!(ax, unit_arrow_head_size, arrow_head_scale, unit_arrow_vector, arrow_vector, tail_location, linewidth, arrowcolor, linecolor)
        end

        if qZ != 0.0
            unit_arrow_vector = [0.0, 0.0, qZ] / norm([0.0, 0.0, qZ])
            arrow_vector = unit_arrow_vector .* arrow_scale
            define_load_arrow!(ax, unit_arrow_head_size, arrow_head_scale, unit_arrow_vector, arrow_vector, tail_location, linewidth, arrowcolor, linecolor)
        end

    end

end


function define_load_arrow!(ax, unit_arrow_head_size, arrow_head_scale, unit_arrow_vector, arrow_vector, tail_location, linewidth, arrowcolor, linecolor)
    
 
    arrow_head_size = unit_arrow_head_size * arrow_head_scale

    arrow_size_vector = arrow_head_scale * unit_arrow_vector


    arrows!(ax, [tail_location[1]-arrow_vector[1] - arrow_size_vector[1]], [tail_location[2]-arrow_vector[2] - arrow_size_vector[2]], [tail_location[3]-arrow_vector[3]-arrow_size_vector[3]], [arrow_vector[1]], [arrow_vector[2]], [arrow_vector[3]], arrowsize = arrow_head_size, linewidth=linewidth,
        arrowcolor = arrowcolor, linecolor = linecolor)

end


function element_local_axes!(ax, element, node, model, unit_arrow_head_size, arrow_head_scale, arrow_scale, arrowcolor, linecolor, linewidth)

    for i in eachindex(element.numbers)

        element_nodes = element.nodes[i]

        #start node

        start_index = findfirst(num->num==element_nodes[1], node.numbers)
        end_index = findfirst(num->num==element_nodes[2], node.numbers)

        Δ = node.coordinates[end_index] .- node.coordinates[start_index]

        tail_location = node.coordinates[start_index] .+ Δ./2


        unit_vector_Y = [0.0, 1.0, 0.0]
        local_Y = model.properties.Γ[i][1:3,1:3]' * unit_vector_Y
        # unit_arrow_vector = global_Y
        arrow_vector = local_Y .* arrow_scale
        arrow_head_size = unit_arrow_head_size * arrow_head_scale
        arrows!(ax, [tail_location[1]], [tail_location[2]], [tail_location[3]], [arrow_vector[1]], [arrow_vector[2]], [arrow_vector[3]], arrowsize = arrow_head_size, linewidth=linewidth,
            arrowcolor = arrowcolor, linecolor = linecolor)

        unit_vector_X = [1.0, 0.0, 0.0]
        local_X = model.properties.Γ[i][1:3,1:3]' * unit_vector_X
        arrow_vector = local_X .* arrow_scale
        arrow_head_size = unit_arrow_head_size * arrow_head_scale
        arrows!(ax, [tail_location[1]], [tail_location[2]], [tail_location[3]], [arrow_vector[1]], [arrow_vector[2]], [arrow_vector[3]], arrowsize = arrow_head_size, linewidth=linewidth,
            arrowcolor = arrowcolor, linecolor = linecolor)

    end

end

function deformed_shape!(ax, nodal_displacements, global_dof, element, node, properties, element_connections, n, scale, linecolor)

    u = Array{Float64}(undef, 0)

    for i in eachindex(nodal_displacements)

        u = [u; nodal_displacements[i]]

    end

    u_global_e = InstantFrame.Show.define_global_element_displacements(u, global_dof, element, element_connections)
    u_local_e = [properties.Γ[i]*u_global_e[i] for i in eachindex(u_global_e)]

    element_display_coords, element_display_Δ = InstantFrame.Show.get_display_coords(element, node, properties, u_local_e, n)

    for i in eachindex(element_display_coords)
        InstantFrame.Show.show_element_deformed_shape!(ax, element_display_coords[i], element_display_Δ[i], scale, linecolor)
    end

end

function node_numbers!(ax, node, fontsize, color)

    text!(ax,
        [Point3f(node.coordinates[i][1], node.coordinates[i][2], node.coordinates[i][3]) for i in eachindex(node.coordinates)],
        text = [string(node.numbers[i]) for i in eachindex(node.numbers)],
        # rotation = [i / 7 * 1.5pi for i in 1:7],
        color = color,
        # align = (:left, :baseline),
        fontsize = fontsize,
        # markerspace = :data
    )

end

function element_numbers!(ax, element, node, fontsize, color)


    text_location = Vector{Tuple{Float64, Float64, Float64}}(undef, size(element.numbers, 1))
    for i in eachindex(element.numbers)

        start_index = findfirst(num->num==element.nodes[i][1], node.numbers)
        end_index = findfirst(num->num==element.nodes[i][2], node.numbers)

        Δ = node.coordinates[end_index] .- node.coordinates[start_index]

        text_location[i] = node.coordinates[start_index] .+ Δ./2

    end

    text!(ax,
        [Point3f(text_location[i][1], text_location[i][2], text_location[i][3]) for i in eachindex(text_location)],
        text = [string(element.numbers[i]) for i in eachindex(element.numbers)],
        # rotation = [i / 7 * 1.5pi for i in 1:7],
        color = color,
        # align = (:left, :baseline),
        fontsize = fontsize,
        # markerspace = :data
    )
  

end

end #module

