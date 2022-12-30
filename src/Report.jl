module Report

using DataFrames, CSV

using ..InstantFrame

function nodal_inputs(node, support; filename)

    #add node numbers
    report = DataFrame(node_numbers=node.numbers)

    #nodal coordinates
    nodal_coords = Matrix{Float64}(undef, size(node.coordinates, 1), 3)
    for i in eachindex(node.coordinates)

        nodal_coords[i,:] = [node.coordinates[i][1], node.coordinates[i][2], node.coordinates[i][3]]

    end

    report[!, :X] = nodal_coords[:,1]
    report[!, :Y] = nodal_coords[:,2]
    report[!, :Z] = nodal_coords[:,3]

    #add support info
    report[!, :uX] = zeros(size(node.numbers, 1))
    report[!, :uY] = zeros(size(node.numbers, 1))
    report[!, :uZ] = zeros(size(node.numbers, 1))
    report[!, :rX] = zeros(size(node.numbers, 1))
    report[!, :rY] = zeros(size(node.numbers, 1))
    report[!, :rZ] = zeros(size(node.numbers, 1))

    for i in eachindex(support.nodes)

        index = findfirst(num->num==support.nodes[i], node.numbers)

        report[index, :uX] = support.stiffness.uX[i]
        report[index, :uY] = support.stiffness.uY[i]
        report[index, :uZ] = support.stiffness.uZ[i]
        report[index, :rX] = support.stiffness.rX[i]
        report[index, :rY] = support.stiffness.rY[i]
        report[index, :rZ] = support.stiffness.rZ[i]


    end

    # #add point loads
    # report[!, :FX] = zeros(size(node.numbers, 1))
    # report[!, :FY] = zeros(size(node.numbers, 1))
    # report[!, :FZ] = zeros(size(node.numbers, 1))
    # report[!, :MX] = zeros(size(node.numbers, 1))
    # report[!, :MY] = zeros(size(node.numbers, 1))
    # report[!, :MZ] = zeros(size(node.numbers, 1))

    # if point_load != InstantFrame.PointLoad(nothing)

    #     for i in eachindex(point_load.nodes)

    #         index = findfirst(num->num==point_load.nodes[i], node.numbers)

    #         report[index, :FX] = point_load.magnitudes.FX[i]
    #         report[index, :FY] = point_load.magnitudes.FY[i]
    #         report[index, :FZ] = point_load.magnitudes.FZ[i]
    #         report[index, :MX] = point_load.magnitudes.MX[i]
    #         report[index, :MY] = point_load.magnitudes.MY[i]
    #         report[index, :MZ] = point_load.magnitudes.MZ[i]

    #     end

    # end

    CSV.write(filename, report)

    return report

end


function element_inputs(element, properties; filename)

    #add element numbers
    report = DataFrame(element_numbers=element.numbers)

    #add element info
    node_i = [element.nodes[i][1] for i in eachindex(element.nodes)]
    node_j = [element.nodes[i][2] for i in eachindex(element.nodes)]
    report[!, :node_i] = node_i
    report[!, :node_j] = node_j
    report[!, :L] = properties.L
    report[!, :local_axis_rotation] = element.orientation

    #add local element axes unit vectors

    local_xX = [properties.local_axis_directions[i].x[1] for i in eachindex(properties.local_axis_directions)]    
    local_xY = [properties.local_axis_directions[i].x[2] for i in eachindex(properties.local_axis_directions)]  
    local_xZ = [properties.local_axis_directions[i].x[3] for i in eachindex(properties.local_axis_directions)]  

    local_yX = [properties.local_axis_directions[i].y[1] for i in eachindex(properties.local_axis_directions)]    
    local_yY = [properties.local_axis_directions[i].y[2] for i in eachindex(properties.local_axis_directions)]  
    local_yZ = [properties.local_axis_directions[i].y[3] for i in eachindex(properties.local_axis_directions)]  

    local_zX = [properties.local_axis_directions[i].z[1] for i in eachindex(properties.local_axis_directions)]    
    local_zY = [properties.local_axis_directions[i].z[2] for i in eachindex(properties.local_axis_directions)]  
    local_zZ = [properties.local_axis_directions[i].z[3] for i in eachindex(properties.local_axis_directions)]  

    report[!, :xX] = local_xX 
    report[!, :xY] = local_xY 
    report[!, :xZ] = local_xZ 

    report[!, :yX] = local_yX 
    report[!, :yY] = local_yY 
    report[!, :yZ] = local_yZ 

    report[!, :zX] = local_zX 
    report[!, :zY] = local_zY 
    report[!, :zZ] = local_zZ 


    #add material properties
    report[!, :E] = properties.E
    report[!, :Poissons_ratio] = properties.ν
    report[!, :G] = properties.G
    report[!, :mass_density] = properties.ρ

    #add section properties
    report[!, :A] = properties.A
    report[!, :Iyy] = properties.Iy
    report[!, :Izz] = properties.Iz
    report[!, :Io] = properties.Io
    report[!, :J] = properties.J

    #add element end conditions
    report[!, :ux_i] = fill(Inf, size(element.numbers, 1))
    report[!, :uy_i] = fill(Inf, size(element.numbers, 1))
    report[!, :uz_i] = fill(Inf, size(element.numbers, 1))
    report[!, :rx_i] = fill(Inf, size(element.numbers, 1))
    report[!, :ry_i] = fill(Inf, size(element.numbers, 1))
    report[!, :rz_i] = fill(Inf, size(element.numbers, 1))

    report[!, :ux_j] = fill(Inf, size(element.numbers, 1))
    report[!, :uy_j] = fill(Inf, size(element.numbers, 1))
    report[!, :uz_j] = fill(Inf, size(element.numbers, 1))
    report[!, :rx_j] = fill(Inf, size(element.numbers, 1))
    report[!, :ry_j] = fill(Inf, size(element.numbers, 1))
    report[!, :rz_j] = fill(Inf, size(element.numbers, 1))

    for i in eachindex(properties.start_connection)
        report.ux_i[i] = properties.start_connection[i].ux
        report.uy_i[i] = properties.start_connection[i].uy
        report.uz_i[i] = properties.start_connection[i].uz
        report.rx_i[i] = properties.start_connection[i].rx
        report.ry_i[i] = properties.start_connection[i].ry
        report.rz_i[i] = properties.start_connection[i].rz
    end

    for i in eachindex(properties.end_connection)
        report.ux_j[i] = properties.end_connection[i].ux
        report.uy_j[i] = properties.end_connection[i].uy
        report.uz_j[i] = properties.end_connection[i].uz
        report.rx_j[i] = properties.end_connection[i].rx
        report.ry_j[i] = properties.end_connection[i].ry
        report.rz_j[i] = properties.end_connection[i].rz
    end

    # #add element uniform loads
    # report[!, :qX] = zeros(size(element.numbers, 1))
    # report[!, :qY] = zeros(size(element.numbers, 1))
    # report[!, :qZ] = zeros(size(element.numbers, 1))
    # report[!, :mX] = zeros(size(element.numbers, 1))
    # report[!, :mY] = zeros(size(element.numbers, 1))
    # report[!, :mZ] = zeros(size(element.numbers, 1))

    # if uniform_load != InstantFrame.UniformLoad(nothing)

    #     for i in eachindex(uniform_load.elements)

    #         index = findfirst(num->num==uniform_load.elements[i], element.numbers)

    #         report[index, :qX] = uniform_load.magnitudes.qX[i]
    #         report[index, :qY] = uniform_load.magnitudes.qY[i]
    #         report[index, :qZ] = uniform_load.magnitudes.qZ[i]
    #         report[index, :mX] = uniform_load.magnitudes.mX[i]
    #         report[index, :mY] = uniform_load.magnitudes.mY[i]
    #         report[index, :mZ] = uniform_load.magnitudes.mZ[i]

    #     end

    # end

    CSV.write(filename, report)

    return report

end

function nodal_output(node, support, displacements, reactions; filename)

    #add node numbers
    report = DataFrame(node_numbers=node.numbers)

    #add nodal displacements

    report[!, :uX] = zeros(Float64, size(node.numbers, 1))
    report[!, :uY] = zeros(Float64, size(node.numbers, 1))
    report[!, :uZ] = zeros(Float64, size(node.numbers, 1))
    report[!, :rX] = zeros(Float64, size(node.numbers, 1))
    report[!, :rY] = zeros(Float64, size(node.numbers, 1))
    report[!, :rZ] = zeros(Float64, size(node.numbers, 1))

    for i in eachindex(displacements)

        report[i, :uX] = displacements[i][1]
        report[i, :uY] = displacements[i][2]
        report[i, :uZ] = displacements[i][3]
        report[i, :rX] = displacements[i][4]
        report[i, :rY] = displacements[i][5]
        report[i, :rZ] = displacements[i][6]

    end

    #add nodal reactions
    report[!, :RX] = zeros(Float64, size(node.numbers, 1))
    report[!, :RY] = zeros(Float64, size(node.numbers, 1))
    report[!, :RZ] = zeros(Float64, size(node.numbers, 1))
    report[!, :MX] = zeros(Float64, size(node.numbers, 1))
    report[!, :MY] = zeros(Float64, size(node.numbers, 1))
    report[!, :MZ] = zeros(Float64, size(node.numbers, 1))

    for i in eachindex(support.nodes)

        index = findfirst(num-> num==support.nodes[i], node.numbers)
        report[index, :RX] = reactions[i][1]
        report[index, :RY] = reactions[i][2]
        report[index, :RZ] = reactions[i][3]
        report[index, :MX] = reactions[i][4]
        report[index, :MY] = reactions[i][5]
        report[index, :MZ] = reactions[i][6]

    end

    CSV.write(filename, report)

    return report

end

function element_output(element, forces, connections; filename)

    #add element numbers
    report = DataFrame(element_numbers=element.numbers)

    #add element forces
    report[!, :Px_i] = zeros(Float64, size(element.numbers, 1))
    report[!, :Py_i] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :Pz_i] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :Mx_i] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :My_i] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :Mz_i] = zeros(Float64, size(element.numbers, 1))  

    report[!, :Px_j] = zeros(Float64, size(element.numbers, 1))
    report[!, :Py_j] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :Pz_j] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :Mx_j] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :My_j] = zeros(Float64, size(element.numbers, 1)) 
    report[!, :Mz_j] = zeros(Float64, size(element.numbers, 1))  

    for i in eachindex(element.numbers)

        report[i, :Px_i] = forces[i][1]
        report[i, :Py_i] = forces[i][2]
        report[i, :Pz_i] = forces[i][3]
        report[i, :Mx_i] = forces[i][4]
        report[i, :My_i] = forces[i][5]
        report[i, :Mz_i] = forces[i][6]

        report[i, :Px_j] = forces[i][7]
        report[i, :Py_j] = forces[i][8]
        report[i, :Pz_j] = forces[i][9]
        report[i, :Mx_j] = forces[i][10]
        report[i, :My_j] = forces[i][11]
        report[i, :Mz_j] = forces[i][12]

    end

    #add element end displacements for partially restrained connections

    if isempty(connections.elements)

        report[!, :ux_i] = fill("", size(element.numbers, 1))
        report[!, :uy_i] = fill("", size(element.numbers, 1))
        report[!, :uz_i] = fill("", size(element.numbers, 1))
        report[!, :rx_i] = fill("", size(element.numbers, 1))
        report[!, :ry_i] = fill("", size(element.numbers, 1))
        report[!, :rz_i] = fill("", size(element.numbers, 1))

        report[!, :ux_j] = fill("", size(element.numbers, 1))
        report[!, :uy_j] = fill("", size(element.numbers, 1))
        report[!, :uz_j] = fill("", size(element.numbers, 1))
        report[!, :rx_j] = fill("", size(element.numbers, 1))
        report[!, :ry_j] = fill("", size(element.numbers, 1))
        report[!, :rz_j] = fill("", size(element.numbers, 1))

    else

        report[!, :ux_i] = fill(0.0, size(element.numbers, 1))
        report[!, :uy_i] = fill(0.0, size(element.numbers, 1))
        report[!, :uz_i] = fill(0.0, size(element.numbers, 1))
        report[!, :rx_i] = fill(0.0, size(element.numbers, 1))
        report[!, :ry_i] = fill(0.0, size(element.numbers, 1))
        report[!, :rz_i] = fill(0.0, size(element.numbers, 1))

        report[!, :ux_j] = fill(0.0, size(element.numbers, 1))
        report[!, :uy_j] = fill(0.0, size(element.numbers, 1))
        report[!, :uz_j] = fill(0.0, size(element.numbers, 1))
        report[!, :rx_j] = fill(0.0, size(element.numbers, 1))
        report[!, :ry_j] = fill(0.0, size(element.numbers, 1))
        report[!, :rz_j] = fill(0.0, size(element.numbers, 1))

        for i in eachindex(connections.elements)

            index = findfirst(num->num==connections.elements[i], element.numbers)

            report[index, :ux_i] = connections.displacements[i][1]
            report[index, :uy_i] = connections.displacements[i][2] 
            report[index, :uz_i] = connections.displacements[i][3] 
            report[index, :rx_i] = connections.displacements[i][4] 
            report[index, :ry_i] = connections.displacements[i][5] 
            report[index, :rz_i] = connections.displacements[i][6]  

            report[index, :ux_j] = connections.displacements[i][7]
            report[index, :uy_j] = connections.displacements[i][8] 
            report[index, :uz_j] = connections.displacements[i][9] 
            report[index, :rx_j] = connections.displacements[i][10] 
            report[index, :ry_j] = connections.displacements[i][11] 
            report[index, :rz_j] = connections.displacements[i][12]
            
        end

    end

    CSV.write(filename, report)

    return report

end

function uniform_loads(uniform_load; filename)

    #add element numbers
    report = DataFrame(element_numbers=uniform_load.elements)

    #add load case labels
    report[!, :load_case] = uniform_load.labels

    #add element uniform loads
    report[!, :qX] = uniform_load.magnitudes.qX
    report[!, :qY] = uniform_load.magnitudes.qY
    report[!, :qZ] = uniform_load.magnitudes.qZ
    report[!, :mX] = uniform_load.magnitudes.mX
    report[!, :mY] = uniform_load.magnitudes.mY
    report[!, :mZ] = uniform_load.magnitudes.mZ

    CSV.write(filename, report)

    return report

end


function point_loads(point_load; filename)

    #add node numbers
    report = DataFrame(element_numbers=point_load.nodes)

    #add load case labels
    report[!, :load_case] = point_load.labels

    #add nodal point loads
    report[!, :FX] = point_load.magnitudes.FX
    report[!, :FY] = point_load.magnitudes.FY
    report[!, :FZ] = point_load.magnitudes.FZ
    report[!, :MX] = point_load.magnitudes.MX
    report[!, :MY] = point_load.magnitudes.MY
    report[!, :MZ] = point_load.magnitudes.MZ

    CSV.write(filename, report)

    return report

end


end #module

