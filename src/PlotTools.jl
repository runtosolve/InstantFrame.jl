module PlotTools

using CairoMakie


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

# function show_model_deformed_shape(nodal_displacements, element_connections, element, node, properties, scale)

#     n = vec(ones(Int64, size(properties.L, 1)) * 11)

#     u = Array{Float64}(undef, 0)

#     for i in eachindex(nodal_displacements)

#         u = [u; nodal_displacements[i]]

#     end

#     u_global_e = define_global_element_displacements(u, properties.global_dof, element, element_connections)

#     u_local_e = [properties.Γ[i]*u_global_e[i] for i in eachindex(u_global_e)]

#     element_display_coords, element_display_Δ = InstantFrame.UI.get_display_coords(element, node, properties, u_local_e, n)

#     figure = Figure()
#     ax = Axis3(figure[1,1])

#     for i in eachindex(element_display_coords)
#         ax = show_element_deformed_shape(element_display_coords[i], element_display_Δ[i], scale, ax)
#     end

#     # ax.azimuth[] = 3π/2
#     # ax.elevation[] = π/2
#     # ax.aspect = (1.0, Y_range/X_range, Z_range/X_range)

#     return figure

# end


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

function get_XYZ_element_ij(element_nodal_coords)

    Xij=Float64[]
    Yij=Float64[]
    Zij=Float64[]
    for i in eachindex(element_nodal_coords)

        Xij = push!(Xij, element_nodal_coords[i].start_node[1])
        Xij = push!(Xij, element_nodal_coords[i].end_node[1])
        Yij = push!(Yij, element_nodal_coords[i].start_node[2])
        Yij = push!(Yij, element_nodal_coords[i].end_node[2])
        Zij = push!(Zij, element_nodal_coords[i].start_node[3])
        Zij = push!(Zij, element_nodal_coords[i].end_node[3])

    end

    return Xij, Yij, Zij

end


function get_text_location_on_elements(element, node)

    text_location = Vector{Tuple{Float64, Float64, Float64}}(undef, size(element.numbers, 1))
    for i in eachindex(element.numbers)

        start_index = findfirst(num->num==element.nodes[i][1], node.numbers)
        end_index = findfirst(num->num==element.nodes[i][2], node.numbers)

        Δ = node.coordinates[end_index] .- node.coordinates[start_index]

        text_location[i] = node.coordinates[start_index] .+ Δ./2

    end

    return text_location

end

end #module