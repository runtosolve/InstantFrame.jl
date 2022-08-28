module UI

using GLMakie


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

end #module

