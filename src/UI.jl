module UI

using GLMakie, InstantFrame


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

    #local element deformation
    δ = [zeros(Float64, 6) for i in eachindex(x)]

    for i in eachindex(x)

        δ[i][2] = w_xy[i]
        δ[i][3] = w_xz[i]

    end

    #global element deformation
    Δ = [Γ[1:6,1:6] * δ[i] for i in eachindex(x)]

    return δ, Δ, x

end



function discretized_element_global_coords(node_i_coords, Γ, x)

    local_element_discretized_coords = [zeros(Float64, 3) for i in eachindex(x)]

    [local_element_discretized_coords[i][1] = x[i] for i in eachindex(x)]

    global_element_discretized_coords = [Γ'[1:3, 1:3] * (local_element_discretized_coords[i] .+  node_i_coords) for i in eachindex(x)]

    return global_element_discretized_coords

end


function get_display_coords_element(u_local_e, node_i_coords, L, Γ, n)

    δe, Δe, x = get_element_deformed_shape(u_local_e, L, Γ, n)
    discretized_element_global_coords = UI.discretized_element_global_coords(node_i_coords, Γ, x)

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



function display_element_deformed_shape(element_XYZ, Δ, scale, ax)

    X = [element_XYZ[i][1] for i in eachindex(element_XYZ)]
    Y = [element_XYZ[i][2] for i in eachindex(element_XYZ)]
    Z = [element_XYZ[i][3] for i in eachindex(element_XYZ)]

    ΔX = [Δ[i][1] for i in eachindex(element_XYZ)]
    ΔY = [Δ[i][2] for i in eachindex(element_XYZ)]
    ΔZ = [Δ[i][3] for i in eachindex(element_XYZ)]


    for i=1:(length(X)-1)

        scatterlines!(ax, [X[i], X[i+1]], [Y[i], Y[i+1]], [Z[i], Z[i+1]], markersize = 5)
        scatterlines!(ax, [X[i] + scale[1] * ΔX[i], X[i+1] + scale[1] * ΔX[i+1]], [Y[i] + scale[2] * ΔY[i], Y[i+1] + scale[2] * ΔY[i+1]], [Z[i] + scale[3] * ΔZ[i], Z[i+1] + scale[3] * ΔZ[i+1]], markersize = 5,  linestyle=:dash)

    end

    return ax

end

function display_model_deformed_shape(u, element, node, properties, scale)

    n = vec(ones(Int64, size(properties.L, 1)) * 11)

    u_global_e = InstantFrame.define_global_element_displacements(u, properties.global_dof)

    u_local_e = [properties.Γ[i]*u_global_e[i] for i in eachindex(u_global_e)]

    element_display_coords, element_display_Δ = InstantFrame.UI.get_display_coords(element, node, properties, u_local_e, n)

    figure = Figure()
    ax = Axis3(figure[1,1])

    for i in eachindex(element_display_coords)
        ax = display_element_deformed_shape(element_display_coords[i], element_display_Δ[i], scale, ax)
    end

    # ax.azimuth[] = 3π/2
    # ax.elevation[] = π/2
    # ax.aspect = (1.0, Y_range/X_range, Z_range/X_range)

    return figure

end

end #module

