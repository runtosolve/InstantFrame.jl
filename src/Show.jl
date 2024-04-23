module Show

using CairoMakie, LinearAlgebra
using ..PlotTools


#3D
function elements!(ax, element_nodal_coords, color)

    for i in eachindex(element_nodal_coords)
        lines!(ax, [element_nodal_coords[i].start_node[1], element_nodal_coords[i].end_node[1]], [element_nodal_coords[i].start_node[2], element_nodal_coords[i].end_node[2]], [element_nodal_coords[i].start_node[3], element_nodal_coords[i].end_node[3]], color = color)

    end

end

#2D
function elements!(ax, Xij, Yij, attributes)

    linesegments!(ax, Xij, Yij, color = attributes.color, linewidth = attributes.linewidth)

end



function nodes!(ax, X, Y, Z, markersize, color)

    scatter!(ax, X, Y, Z, markersize=markersize, color = color)

end

function nodes!(ax, X, Y, attributes)

    scatter!(ax, X, Y, markersize=attributes.size, color = attributes.color)

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

    u_global_e = InstantFrame.PlotTools.define_global_element_displacements(u, global_dof, element, element_connections)
    u_local_e = [properties.Γ[i]*u_global_e[i] for i in eachindex(u_global_e)]

    element_display_coords, element_display_Δ = InstantFrame.PlotTools.get_display_coords(element, node, properties, u_local_e, n)

    for i in eachindex(element_display_coords)
        InstantFrame.PlotTools.show_element_deformed_shape!(ax, element_display_coords[i], element_display_Δ[i], scale, linecolor)
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

#2D
function axial_force!(ax, Xij, Yij, axial_forces, attributes)

    num_elements = Int(size(Xij)[1]/2)
    force_color = Vector{String}(undef, num_elements)
    linewidths = Vector{Float64}(undef, num_elements)

    max_force = maximum(abs.(axial_forces))

    for i=1:num_elements

        if axial_forces[i] >= 0.0
            force_color[i] = attributes.tension_color
        else
            force_color[i] = attributes.compression_color
        end

        linewidths[i] = abs(axial_forces[i])/max_force * attributes.scale

    end

    linesegments!(ax, Xij, Yij, color = force_color, linewidth = linewidths)

end

function axial_force_magnitude!(ax, element, node, axial_forces, active_element_index, attributes)


    text_location = PlotTools.get_text_location_on_elements(element, node)

    three_dimensional = isempty(findall(coord->coord == 0.0, [text_location[i][3] for i in eachindex(text_location)]))

    if three_dimensional

        text!(ax,
        [Point3f(text_location[i][1], text_location[i][2], text_location[i][3]) for i in eachindex(text_location)],
        text = [string(axial_forces[i]) for i in eachindex(axial_forces)],
        # rotation = [i / 7 * 1.5pi for i in 1:7],
        color = attributes.color,
        # align = (:left, :baseline),
        fontsize = attributes.fontsize,
        # markerspace = :data)
        )
  
    else

        text!(ax,
        [Point2f(text_location[active_element_index[i]][1], text_location[active_element_index[i]][2]) for i in eachindex(active_element_index)],
        text = [string(axial_forces[active_element_index[i]]) for i in eachindex(active_element_index)],
        # rotation = [i / 7 * 1.5pi for i in 1:7],
        color = attributes.color,
        # align = (:left, :baseline),
        fontsize = attributes.fontsize,
        # markerspace = :data)
        )
  
    end

end
       

end #module

