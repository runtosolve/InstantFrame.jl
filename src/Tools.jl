module Tools

using DataFrames, CSV


function report_nodal_inputs(node, support, point_load, save_filename)

    #add node numbers
    input = DataFrame(node_numbers=node.numbers)

    #nodal coordinates
    nodal_coords = Matrix{Float64}(undef, size(node.coordinates, 1), 3)
    for i in eachindex(node.coordinates)

        nodal_coords[i,:] = [node.coordinates[i][1], node.coordinates[i][2], node.coordinates[i][3]]

    end

    input[!, :X] = nodal_coords[:,1]
    input[!, :Y] = nodal_coords[:,2]
    input[!, :Z] = nodal_coords[:,3]

    #add support info
    input[!, :uX] = zeros(size(node.numbers, 1))
    input[!, :uY] = zeros(size(node.numbers, 1))
    input[!, :uZ] = zeros(size(node.numbers, 1))
    input[!, :rX] = zeros(size(node.numbers, 1))
    input[!, :rY] = zeros(size(node.numbers, 1))
    input[!, :rZ] = zeros(size(node.numbers, 1))

    for i in eachindex(support.nodes)

        index = findfirst(num->num==support.nodes[i], node.numbers)

        input[index, :uX] = support.stiffness.uX[i]
        input[index, :uY] = support.stiffness.uY[i]
        input[index, :uZ] = support.stiffness.uZ[i]
        input[index, :rX] = support.stiffness.rX[i]
        input[index, :rY] = support.stiffness.rY[i]
        input[index, :rZ] = support.stiffness.rZ[i]


    end

    #add point loads
    input[!, :FX] = zeros(size(node.numbers, 1))
    input[!, :FY] = zeros(size(node.numbers, 1))
    input[!, :FZ] = zeros(size(node.numbers, 1))
    input[!, :MX] = zeros(size(node.numbers, 1))
    input[!, :MY] = zeros(size(node.numbers, 1))
    input[!, :MZ] = zeros(size(node.numbers, 1))

    if point_load != InstantFrame.PointLoad(nothing)

        for i in eachindex(support.nodes)

            index = findfirst(num->num==support.nodes[i], node.numbers)

            input[index, :FX] = point_load.magnitudes.FX[i]
            input[index, :FY] = point_load.magnitudes.FY[i]
            input[index, :FZ] = point_load.magnitudes.FZ[i]
            input[index, :MX] = point_load.magnitudes.MX[i]
            input[index, :MY] = point_load.magnitudes.MY[i]
            input[index, :MZ] = point_load.magnitudes.MZ[i]

        end

    end

    CSV.write(save_filename, input)

    return input

end

end #module

