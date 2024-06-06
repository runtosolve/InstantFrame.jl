module InstantFrame

using LinearAlgebra, Rotations, Parameters, NonlinearSolve, StaticArrays


export PlotTools
include("PlotTools.jl")
using .PlotTools

export Show
include("Show.jl")
using .Show


# export Report
# include("Report.jl")
# using .Report

# export Export
# include("Export.jl")
# using .Export


@with_kw mutable struct Material

    names::Array{String}
    E::Array{Float64}
    ν::Array{Float64}
    ρ::Array{Float64}

end

@with_kw mutable struct CrossSection

    names::Array{String}
    A::Array{Float64}
    Iy::Array{Float64}
    Iz::Array{Float64}
    J::Array{Float64}

end

@with_kw mutable struct Connection

    names::Array{String}
    stiffness::NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Vector{Float64}}}

end

@with_kw mutable struct Node

    numbers::Array{Int64}
    coordinates::Array{Tuple{Float64, Float64, Float64}, 1}

end

@with_kw mutable struct Element

    types::Array{String}
    numbers::Array{Int64}
    nodes::Array{Tuple{Int64, Int64}}
    orientation::Array{Float64}
    connections::Array{Tuple{String, String}, 1}
    cross_section::Array{String}
    material::Array{String}

end

@with_kw mutable struct Support

    nodes::Array{Int64}
    stiffness::NamedTuple{(:uX, :uY, :uZ, :rX, :rY, :rZ), NTuple{6, Vector{Float64}}}

end

@with_kw mutable struct UniformLoad

    labels::Union{Array{String}, Nothing}
    elements::Union{Array{Int64}, Nothing}
    magnitudes::Union{NamedTuple{(:qX, :qY, :qZ, :mX, :mY, :mZ), NTuple{6, Vector{Float64}}}, Nothing}

end

UniformLoad(nothing) = UniformLoad(labels=nothing, elements=nothing, magnitudes=nothing)

@with_kw mutable struct PointLoad

    labels::Union{Array{String}, Nothing}
    nodes::Union{Array{Int64}, Nothing}
    magnitudes::Union{NamedTuple{(:FX, :FY, :FZ, :MX, :MY, :MZ), NTuple{6, Vector{Float64}}}, Nothing}

end

PointLoad(nothing) = PointLoad(labels=nothing, nodes=nothing, magnitudes=nothing)

@with_kw mutable struct ElementProperties

    L::Array{Float64}
    A::Array{Float64}
    Iz::Array{Float64}
    Iy::Array{Float64}
    Io::Array{Float64}
    J::Array{Float64}
    E::Array{Float64}
    ν::Array{Float64}
    G::Array{Float64}
    ρ::Array{Float64}

    rotation_angles::Array{NamedTuple{(:Y, :Z, :X), NTuple{3, Float64}}}
    Γ::Array{Array{Float64, 2}}
    local_axis_directions::Array{NamedTuple{(:x, :y, :z), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}
    global_dof::Array{Array{Int64, 1}}

    start_connection::Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}
    end_connection::Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}

end

@with_kw struct ElementConnections
    elements::Array{Int64, 1}
    displacements::Array{Array{Float64, 1}}
end

@with_kw struct ElasticSupport

    global_dof::Array{Int64}
    global_stiffness::Array{Float64}

end



@with_kw struct FirstOrderEquations
  
    free_dof::Array{Int64}
    fixed_dof::Array{Int64}
    elastic_supports::ElasticSupport
    ke_local::Array{Array{Float64, 2}}
    ke_global::Array{Array{Float64, 2}}
    Ke::Array{Float64, 2}
  
end

@with_kw struct SecondOrderEquations
  
    free_dof::Array{Int64}
    fixed_dof::Array{Int64}
    elastic_supports::ElasticSupport
    ke_local::Array{Array{Float64, 2}}
    ke_global::Array{Array{Float64, 2}}
    Ke::Array{Float64, 2}
    kg_local::Array{Array{Float64, 2}}
    kg_global::Array{Array{Float64, 2}}
    Kg::Array{Float64, 2}
  
end


@with_kw struct ModalEquations
  
    free_dof::Array{Int64}
    ke_local::Array{Array{Float64, 2}}
    ke_global::Array{Array{Float64, 2}}
    Ke::Array{Float64, 2}

    m_local::Array{Array{Float64, 2}}
    m_global::Array{Array{Float64, 2}}
    M::Array{Float64, 2}
  
end

@with_kw struct Reactions

    nodes::Vector{Int64}
    magnitudes::Vector{Float64}

end

@with_kw struct Solution

    displacements::Array{Array{Float64, 1}}
    forces::Array{Array{Float64, 1}}
    connections::ElementConnections
    reactions::Array{Array{Float64, 1}}
    u::Vector{Float64}
    uf::Vector{Float64}

end


# @with_kw struct SecondOrderSolution

#     nodal_displacements::Array{Array{Float64, 1}}
#     element_forces::Array{Array{Float64, 1}}
#     element_connections::ElementConnections

# end

@with_kw struct ModalSolution

    ωn::Array{Float64}
    ϕ::Array{Array{Float64, 1}}

end


@with_kw struct Forces

    global_uniform_nodal::Array{Array{Float64}}
    local_fixed_end_forces::Array{Array{Float64}}
    global_dof_nodal_forces_uniform_load::Array{Float64}

    global_dof_point_loads::Array{Float64}

    F::Array{Float64}

end

@with_kw struct Inputs

    node::Node
    cross_section::CrossSection
    material::Material
    connection::Connection
    element::Element
    support::Support
    uniform_load::UniformLoad
    point_load::PointLoad
    analysis_type::String
    solution_tolerance::Union{Float64, Nothing}

end

Inputs(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type) = Inputs(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type, nothing)

@with_kw struct Model

    inputs::Inputs
    properties::ElementProperties
    forces::Union{Forces, Nothing}
    equations::Union{FirstOrderEquations, SecondOrderEquations, ModalEquations}
    solution::Union{Solution, ModalSolution}

end




function define_local_elastic_stiffness_matrix(Iy, Iz, A, J, E, ν, L)


    G = E/(2*(1+ν))

    ke = zeros(Float64, (12, 12))

    ke[1, 1] = E*A/L
    ke[1, 7] = -E*A/L
    ke[2, 2] = 12*E*Iz/L^3
    ke[2, 6] = 6*E*Iz/L^2
    ke[2, 8] = -12*E*Iz/L^3
    ke[2, 12] = 6*E*Iz/L^2
    ke[3, 3] = 12*E*Iy/L^3
    ke[3, 5] = -6*E*Iy/L^2
    ke[3, 9] = -12*E*Iy/L^3
    ke[3, 11] = -6*E*Iy/L^2
    ke[4, 4] = G*J/L
    ke[4, 10] = -G*J/L
    ke[5, 5] = 4*E*Iy/L
    ke[5, 9] = 6*E*Iy/L^2
    ke[5, 11] = 2*E*Iy/L
    ke[6, 6] = 4*E*Iz/L
    ke[6, 8] = -6*E*Iz/L^2
    ke[6, 12] = 2*E*Iz/L
    ke[7, 7] = E*A/L
    ke[8, 8] = 12*E*Iz/L^3
    ke[8, 12] = -6*E*Iz/L^2
    ke[9, 9] = 12*E*Iy/L^3
    ke[9, 11] = 6*E*Iy/L^2
    ke[10, 10] = G*J/L
    ke[11, 11] = 4*E*Iy/L
    ke[12, 12] = 4*E*Iz/L

    for i = 1:12

        for j = 1:12

            ke[j, i] = ke[i, j]

        end
        
    end

    return ke

end


function define_local_geometric_stiffness_matrix(P, L)


    kg = zeros(Float64, (12, 12))

    kg[1, 1] = P/L
    kg[1, 7] = -P/L
    kg[2, 2] = 6/5*P/L
    kg[2, 6] = P/10
    kg[2, 8] = -6/5*P/L
    kg[2, 12] = P/10
    kg[3, 3] = 6/5*P/L
    kg[3, 5] = -P/10
    kg[3, 9] = -6/5*P/L 
    kg[3, 11] = -P/10
    kg[5, 5] = 2/15*P*L 
    kg[5, 9] =  P/10
    kg[5, 11] = -P*L/30
    kg[6, 6] = 2/15*P*L
    kg[6, 8] = -P/10
    kg[6, 12] = -P*L/30
    kg[7, 7] = P/L
    kg[8, 8] = 6/5*P/L
    kg[8, 12] = - P/10
    kg[9, 9] = 6/5 * P/L
    kg[9, 11] = P/10
    kg[11, 11] = 2/15*P*L
    kg[12, 12] = 2/15*P*L
    
    for i = 1:12

        for j = 1:12

            kg[j, i] = kg[i, j]

        end
        
    end
    

    return kg

end


# function define_rotation_matrix(A, B, β)

#     #ρ is x-z plane,  χ is x-y plane , ω is y-z plane
#     AB = B - A

#     length_AB = norm(AB)

#     χ = π/2 - acos(AB[2]/length_AB)

#     if (AB[3]==0.0) & (AB[1]==0.0)  #avoid 0/0 case
#         ρ = 0.0    
#     else
#         ρ = -atan(AB[3]/AB[1])
#     end

#     ω = β

#     γ = RotXZY(-ω, -χ, -ρ)

#     Γ = zeros(Float64, (12, 12))

#     Γ[1:3, 1:3] .= γ
#     Γ[4:6, 4:6] .= γ
#     Γ[7:9, 7:9] .= γ
#     Γ[10:12, 10:12] .= γ

#     return Γ

# end


function define_rotation_matrix(A, B, β)

    AB = B - A

    ΔX = AB[1]
    ΔZ = AB[3]
    ΔY = AB[2]

    #Rotation global coordinate system about Y axis
    ρ = atan(-ΔZ, ΔX)

    #Rotation revised global coordinate system about Z axis
    proj_AB_xz = sqrt(ΔX^2 + ΔZ^2)
    χ = atan(ΔY, proj_AB_xz)

    #Rotation revised global coordinate system about X axis
    # current_local_y_axis = RotZ(-χ) * RotY(-ρ) * [0.0, 1.0, 0.0]  #where Y is pointing after Y and Z rotations 
    # ω = acos(dot(current_local_y_axis, [0.0, 1.0, 0.0])/ norm(current_local_y_axis))

    # # matrix of direction cosines
    # γ = RotX(-(ω+β)) * RotZ(-χ) * RotY(-ρ) # add β here to rotate local y-axis to orient cross-section in global coordinate system 

    ω = β

    γ = RotYZX(ρ, χ, ω)' #transpose!

    Γ = zeros(Float64, (12, 12))

    Γ[1:3, 1:3] .= γ
    Γ[4:6, 4:6] .= γ
    Γ[7:9, 7:9] .= γ
    Γ[10:12, 10:12] .= γ

    angles = (Y=ρ, Z=χ, X=ω)

    return Γ, angles 

end


function define_element_properties(node, cross_section, material, element, connection)

    num_elem = length(element.numbers)

    L = Array{Float64}(undef, num_elem)
    Γ = Array{Array{Float64, 2}}(undef, num_elem)
    global_dof = Array{Array{Int64, 1}}(undef, num_elem)
    Iz = Array{Float64}(undef, num_elem)
    Iy = Array{Float64}(undef, num_elem)
    A = Array{Float64}(undef, num_elem)
    Io = Array{Float64}(undef, num_elem)
    J = Array{Float64}(undef, num_elem)
    E = Array{Float64}(undef, num_elem)
    ν = Array{Float64}(undef, num_elem)
    G = Array{Float64}(undef, num_elem)
    ρ = Array{Float64}(undef, num_elem)
    start_connection = Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}(undef, num_elem)
    end_connection = Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}(undef, num_elem)
    rotation_angles = Array{NamedTuple{(:Y, :Z, :X), NTuple{3, Float64}}}(undef, num_elem)
    local_axis_directions = Array{NamedTuple{(:x, :y, :z), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}(undef, num_elem)

    for i in eachindex(element.numbers)

        #nodal coordinates
        node_i_index = findfirst(node_num->node_num == element.nodes[i][1], node.numbers)
        node_j_index = findfirst(node_num->node_num == element.nodes[i][2], node.numbers)
        node_i = collect(node.coordinates[node_i_index])
        node_j = collect(node.coordinates[node_j_index])

        #cross section
        index = findfirst(name->name == element.cross_section[i], cross_section.names)
        Iz[i] = cross_section.Iz[index]
        Iy[i] = cross_section.Iy[index]
        A[i] = cross_section.A[index]
        J[i] = cross_section.J[index]
        Io[i] = J[i]   #I think this is correct?
        
        #material
        index = findfirst(name->name == element.material[i], material.names)
        E[i] = material.E[index]
        ν[i] = material.ν[index]
        ρ[i] = material.ρ[index]
        G[i] = E[i]/(2*(1+ν[i]))

        #element end conditions
        index = findfirst(name->name == element.connections[i][1], connection.names)
        start_connection[i] = (ux=connection.stiffness.ux[index], uy=connection.stiffness.uy[index], uz=connection.stiffness.uz[index], rx=connection.stiffness.rx[index], ry=connection.stiffness.ry[index], rz=connection.stiffness.rz[index])
        index = findfirst(name->name == element.connections[i][2], connection.names)
        end_connection[i] = (ux=connection.stiffness.ux[index], uy=connection.stiffness.uy[index], uz=connection.stiffness.uz[index], rx=connection.stiffness.rx[index], ry=connection.stiffness.ry[index], rz=connection.stiffness.rz[index])
    
        #element length
        L[i] = norm(node_j - node_i)

        #rotation matrix
        Γ[i], rotation_angles[i] = define_rotation_matrix(node_i, node_j, element.orientation[i])

        #element local axis directions 
        local_axis_directions[i] = calculate_element_local_axis_directions(Γ[i])

        #global dof for each element
        num_dof_per_node = 6 #hard code this for now
        node_i_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_i_index-1)
        node_j_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_j_index-1) 
        global_dof[i] = [node_i_dof; node_j_dof]

    end

    element_properties = ElementProperties(L, A, Iz, Iy, Io, J, E, ν, G, ρ, rotation_angles, Γ, local_axis_directions, global_dof, start_connection, end_connection)

    return element_properties

end


function assemble_global_matrix(k, dof)

    total_dof = maximum([maximum(dof[i]) for i in eachindex(dof)])

    K = zeros(Float64, (total_dof , total_dof))
    
    for i in eachindex(k)
       K[dof[i], dof[i]] += k[i]
    end

    return K

end

function calculate_element_internal_forces(properties, ke_local, element, uniform_load, local_fixed_end_forces, u)

    P_element_local = Array{Array{Float64, 1}}(undef, length(properties.L))

    for i in eachindex(properties.L)
        u_element_global = u[properties.global_dof[i]]
        u_element_local = properties.Γ[i] * u_element_global
        P_element_local[i] = ke_local[i] * u_element_local

        elem_num = element.numbers[i]

        if !isnothing(uniform_load.elements)

            index = findall(num->num==elem_num, uniform_load.elements)

            if !isnothing(index) #add in fixed end forces

                for j in eachindex(index)  #if there are multiple load assignments on one member, need a loop 

                    P_element_local[i][1] += local_fixed_end_forces[index[j]][1]
                    P_element_local[i][2] += local_fixed_end_forces[index[j]][2]
                    P_element_local[i][3] += local_fixed_end_forces[index[j]][3]
                    P_element_local[i][4] += local_fixed_end_forces[index[j]][4]

                    #deal with internal moment sign switch for partially rigid connections
                    # if properties.start_connection[i].ry != Inf
                    #     P_element_local[i][5] -= local_fixed_end_forces[index[j]][5]
                    # else
                        P_element_local[i][5] += local_fixed_end_forces[index[j]][5]
                    # end

                    # if properties.start_connection[i].rz != Inf
                    #     P_element_local[i][6] -= local_fixed_end_forces[index[j]][6]
                    # else
                        P_element_local[i][6] += local_fixed_end_forces[index[j]][6]
                    # end

                    P_element_local[i][7] += local_fixed_end_forces[index[j]][7]
                    P_element_local[i][8] += local_fixed_end_forces[index[j]][8]
                    P_element_local[i][9] += local_fixed_end_forces[index[j]][9]
                    P_element_local[i][10] += local_fixed_end_forces[index[j]][10]


                    # if properties.end_connection[i].ry != Inf
                    #     P_element_local[i][11] -= local_fixed_end_forces[index[j]][11]
                    # else
                        P_element_local[i][11] += local_fixed_end_forces[index[j]][11]
                    # end

                    # if properties.end_connection[i].rz != Inf
                    #     P_element_local[i][12] -= local_fixed_end_forces[index[j]][12]
                    # else
                        P_element_local[i][12] += local_fixed_end_forces[index[j]][12]
                    # end

                end

            end

        end

    end

    return P_element_local

end


function residual(u, p)
    Kff, Ff = p
    Kff * u - Ff
end

function nonlinear_solution(Kff, Ff, u1f)

    p = [Kff, Ff]

    u1f = MVector{length(u1f)}(u1f)
    # u1f_S = zeros(Float64, length(u1f))
    # u1fSS = [u1f_S[i] for i in eachindex(u1f)]
    probN = NonlinearSolve.NonlinearProblem{false}(residual, u1f, p)
    u2f = NonlinearSolve.solve(probN, NewtonRaphson(), tol = 1e-9)

    return u2f

end

function define_nodal_displacements(node, u)

    nodal_displacements = Array{Array{Float64, 1}}(undef, length(node.numbers))
    # nodal_displacements = Vector{Vector{Float64, 1}}(undef, length(node.numbers))

    num_dof_per_node = 6 #hard code this for now
    for i in eachindex(node.numbers)

        nodal_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (i-1)
        nodal_displacements[i] = u[nodal_dof]

    end

    return nodal_displacements

end



function first_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load)

    element_properties = InstantFrame.define_element_properties(node, cross_section, material, element, connection)

    ke_local = [InstantFrame.define_local_elastic_stiffness_matrix(element_properties.Iy[i], element_properties.Iz[i], element_properties.A[i], element_properties.J[i], element_properties.E[i], element_properties.ν[i], element_properties.L[i]) for i in eachindex(element_properties.L)]

    ke_local = InstantFrame.modify_element_local_connection_stiffness(element_properties, ke_local, element)

    ke_global = [element_properties.Γ[i]'*ke_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]

    Ke = InstantFrame.assemble_global_matrix(ke_global, element_properties.global_dof)

    free_global_dof, fixed_global_dof, elastic_supports = InstantFrame.define_free_global_dof(node, support)

    for i in eachindex(elastic_supports.global_dof)  #add springs 

        Ke[elastic_supports.global_dof[i], elastic_supports.global_dof[i]] += elastic_supports.global_stiffness[i]

    end

    equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load = InstantFrame.calculate_nodal_forces_from_uniform_loads(uniform_load, element, node, element_properties)

    global_dof_point_loads = InstantFrame.define_global_dof_point_loads(node, point_load)

    F = global_dof_point_loads .+ global_dof_nodal_forces_uniform_load

    forces = InstantFrame.Forces(equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load, global_dof_point_loads, F)

    
    #calculate nodal displacements
    Ff = F[free_global_dof]
    Ke_ff = Ke[free_global_dof, free_global_dof]
    uf = Ke_ff \ Ff

    u = zeros(Float64, size(Ke,1))
    u[free_global_dof] = uf

    nodal_displacements = define_nodal_displacements(node, u)

    #### calculate reaction ####
    reactions = zeros(Float64, length(node.numbers)*6)

    #get reactions from rigid supports
    Ke_sf = Ke[fixed_global_dof, free_global_dof]
    R = Ke_sf * uf
    reactions[fixed_global_dof] = R

    #get reactions from elastic supports
    elastic_support_reactions = -u[elastic_supports.global_dof] .* elastic_supports.global_stiffness
    reactions[elastic_supports.global_dof] = elastic_support_reactions

    #add point loads to rigid support reactions 
    reactions[fixed_global_dof] += -forces.global_dof_point_loads[fixed_global_dof]

    #add uniform load equivalent nodal forces to rigid support reactions  
    reactions[fixed_global_dof] += -forces.global_dof_nodal_forces_uniform_load[fixed_global_dof] 
    # print(fixed_global_dof)

    # reactions[fixed_global_dof] -= -forces.global_dof_nodal_forces_uniform_load[fixed_global_dof]

    # reactions[fixed_global_dof] = reactions[fixed_global_dof] - 

    #package up reactions at each node, these are reactions from elastic supports and rigid supports
    nodal_reactions = Array{Array{Float64, 1}}(undef, length(support.nodes))

    num_dof_per_node = 6 #hard code this for now
    for i in eachindex(support.nodes)

        nodal_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (support.nodes[i]-1)
        nodal_reactions[i] = reactions[nodal_dof]

    end

    ######

    #calculate element internal forces
    element_forces = InstantFrame.calculate_element_internal_forces(element_properties, ke_local, element, uniform_load, local_fixed_end_forces, u)

    #calculate connection deformations
    element_connections = calculate_element_connection_deformations(element_properties, element, node, nodal_displacements, element_forces)


    equations = FirstOrderEquations(free_global_dof, fixed_global_dof, elastic_supports, ke_local, ke_global, Ke)

    solution = Solution(nodal_displacements, element_forces, element_connections, nodal_reactions, u, uf)

    inputs = Inputs(node, cross_section, material, connection, element, support, uniform_load, point_load, "1st order")

    model = Model(inputs, element_properties, forces, equations, solution)

    return model

end


function second_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load, solution_tolerance)

    element_properties = InstantFrame.define_element_properties(node, cross_section, material, element, connection)

    ke_local = [InstantFrame.define_local_elastic_stiffness_matrix(element_properties.Iy[i], element_properties.Iz[i], element_properties.A[i], element_properties.J[i], element_properties.E[i], element_properties.ν[i], element_properties.L[i]) for i in eachindex(element_properties.L)]

    ke_local = InstantFrame.modify_element_local_connection_stiffness(element_properties, ke_local, element)

    ke_global = [element_properties.Γ[i]'*ke_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]

    Ke = InstantFrame.assemble_global_matrix(ke_global, element_properties.global_dof)

    free_global_dof, fixed_global_dof, elastic_supports = InstantFrame.define_free_global_dof(node, support)

    for i in eachindex(elastic_supports.global_dof)  #add springs 

        Ke[elastic_supports.global_dof[i], elastic_supports.global_dof[i]] += elastic_supports.global_stiffness[i]

    end

    equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load = InstantFrame.calculate_nodal_forces_from_uniform_loads(uniform_load, element, node, element_properties)

    global_dof_point_loads = InstantFrame.define_global_dof_point_loads(node, point_load)

    F = global_dof_point_loads .+ global_dof_nodal_forces_uniform_load

    forces = InstantFrame.Forces(equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load, global_dof_point_loads, F)


    Ff = F[free_global_dof]
    Ke_ff = Ke[free_global_dof, free_global_dof]
    u1f = Ke_ff \ Ff

    u1 = zeros(Float64, size(Ke,1))
    u1[free_global_dof] = u1f

    # P1 = InstantFrame.calculate_element_internal_forces(element_properties, ke_local, u1)
    
    #calculate element internal forces
    P1 = InstantFrame.calculate_element_internal_forces(element_properties, ke_local, element, uniform_load, local_fixed_end_forces, u1)
    P1_axial = [P1[i][7] for i in eachindex(P1)]
    kg_local = [InstantFrame.define_local_geometric_stiffness_matrix(P1_axial[i], element_properties.L[i]) for i in eachindex(element_properties.L)]
    kg_global = [element_properties.Γ[i]'*kg_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]
    Kg = InstantFrame.assemble_global_matrix(kg_global, element_properties.global_dof)
    Kg_ff = Kg[free_global_dof, free_global_dof]

    Kff = Ke_ff + Kg_ff

    p = [Kff, Ff]
    u1f = MVector{length(u1f)}(u1f)
    probN = NonlinearSolve.NonlinearProblem{false}(residual, u1f, p)
    u2f = NonlinearSolve.solve(probN, NewtonRaphson(), reltol = solution_tolerance)

    u2 = zeros(Float64, size(Ke,1))
    u2[free_global_dof] = u2f

    nodal_displacements = define_nodal_displacements(node, u2)

    #### calculate reaction ####
    reactions = zeros(Float64, length(node.numbers)*6)

    #get reactions from rigid supports
    Ke_sf = Ke[fixed_global_dof, free_global_dof]
    R = Ke_sf * u2f
    reactions[fixed_global_dof] = R

    #get reactions from elastic supports
    elastic_support_reactions = -u2[elastic_supports.global_dof] .* elastic_supports.global_stiffness
    reactions[elastic_supports.global_dof] = elastic_support_reactions

    #add point loads to rigid support reactions 
    reactions[fixed_global_dof] += -forces.global_dof_point_loads[fixed_global_dof]

    #add uniform load equivalent nodal forces to rigid support reactions  
    reactions[fixed_global_dof] += -forces.global_dof_nodal_forces_uniform_load[fixed_global_dof]

    #package up reactions at each node, these are reactions from elastic supports and rigid supports
    nodal_reactions = Array{Array{Float64, 1}}(undef, length(support.nodes))

    num_dof_per_node = 6 #hard code this for now
    for i in eachindex(support.nodes)

        nodal_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (support.nodes[i]-1)
        nodal_reactions[i] = reactions[nodal_dof]

    end

    ######

    #calculate element internal forces
    element_forces = InstantFrame.calculate_element_internal_forces(element_properties, ke_local, element, uniform_load, local_fixed_end_forces, u2)

    #calculate connection deformations
    element_connections = calculate_element_connection_deformations(element_properties, element, node, nodal_displacements, element_forces)


 
    equations = SecondOrderEquations(free_global_dof, fixed_global_dof, elastic_supports, ke_local, ke_global, Ke, kg_local, kg_global, Kg)

    solution = Solution(nodal_displacements, element_forces, element_connections, nodal_reactions, u2, u2f)

    inputs = Inputs(node, cross_section, material, connection, element, support, uniform_load, point_load, "2nd order", solution_tolerance)

    model = Model(inputs, element_properties, forces, equations, solution)

    return model

end



function modal_vibration_analysis(node, cross_section, material, connection, element, support)

    element_properties = InstantFrame.define_element_properties(node, cross_section, material, element, connection)

    # free_global_dof = InstantFrame.define_free_global_dof(node, support)

    free_global_dof, fixed_global_dof, elastic_supports = InstantFrame.define_free_global_dof(node, support)

    ke_local = [InstantFrame.define_local_elastic_stiffness_matrix(element_properties.Iy[i], element_properties.Iz[i], element_properties.A[i], element_properties.J[i], element_properties.E[i], element_properties.ν[i], element_properties.L[i]) for i in eachindex(element_properties.L)]

    ke_local = InstantFrame.modify_element_local_connection_stiffness(element_properties, ke_local, element)

    ke_global = [element_properties.Γ[i]'*ke_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]

    Ke = InstantFrame.assemble_global_matrix(ke_global, element_properties.global_dof)

    m_local = [InstantFrame.define_local_3D_mass_matrix(element_properties.A[i], element_properties.L[i], element_properties.Io[i], element_properties.ρ[i]) for i in eachindex(element_properties.L)]
    m_global = [element_properties.Γ[i]'*m_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]
    M = InstantFrame.assemble_global_matrix(m_global, element_properties.global_dof)
    
    Ke_ff = Ke[free_global_dof, free_global_dof]
    Mff = M[free_global_dof, free_global_dof]

    ωn_squared = real.(eigvals(Ke_ff, Mff))
    ωn = sqrt.(ωn_squared)
    ϕff = eigvecs(Ke_ff, Mff)

    ϕ = [zeros(Float64, size(M, 1)) for i in eachindex(ωn)]

    for i in eachindex(ωn)
        ϕ[i][free_global_dof] .= real.(ϕff[:, i])
    end


    equations = ModalEquations(free_global_dof, ke_local, ke_global, Ke, m_local, m_global, M)

    solution = ModalSolution(ωn, ϕ)

    model = Model(element_properties, nothing, equations, solution)

    return model

end


function solve(node, cross_section, material, connection, element, support, uniform_load, point_load; analysis_type, solution_tolerance)

    if analysis_type == "first order"

        model = first_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load)

    elseif analysis_type == "modal vibration"

        model = modal_vibration_analysis(node, cross_section, material, connection, element, support)

    elseif analysis_type == "second order"

        model = second_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load, solution_tolerance)

    end

    return model

end

# function solve(node, cross_section, material, connection, element, support, uniform_load, point_load; analysis_type, solution_tolerance)

#     if analysis_type == "second order"

#         model = second_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load, solution_tolerance)

#     end

#     return model 

# end



function calculate_local_element_fixed_end_forces_partial_restraint(wx_local, wy_local, wz_local, L, k1z, k2z, k1y, k2y, E, Iy, Iz)

    local_fixed_end_forces = calculate_local_element_fixed_end_forces(wx_local, wy_local, wz_local, L)

    k1z, k2z = InstantFrame.connection_zeros_and_inf(k1z, k2z, E, Iz)

    M_fixed = [local_fixed_end_forces[6], local_fixed_end_forces[12]]
    F_fixed = [local_fixed_end_forces[2], local_fixed_end_forces[8]]
    M_spring, F_spring = fixed_end_forces_partial_restraint(k1z, k2z, E, Iz, L, M_fixed, F_fixed)

    local_fixed_end_forces[6] = M_spring[1]
    local_fixed_end_forces[12] = -M_spring[2]
    local_fixed_end_forces[2] = F_spring[2] #tricky here
    local_fixed_end_forces[8] = F_spring[1] #tricky here
  
    k1y, k2y = InstantFrame.connection_zeros_and_inf(k1y, k2y, E, Iy)

    M_fixed = [local_fixed_end_forces[5], local_fixed_end_forces[11]]
    F_fixed = [local_fixed_end_forces[3], local_fixed_end_forces[9]]
    M_spring, F_spring = fixed_end_forces_partial_restraint(k1y, k2y, E, Iy, L, M_fixed, F_fixed)

    local_fixed_end_forces[3] = F_spring[1]
    local_fixed_end_forces[5] = -M_spring[1]
    local_fixed_end_forces[9] = F_spring[2]
    local_fixed_end_forces[11] = M_spring[2]
   
    return local_fixed_end_forces

end



function calculate_local_element_fixed_end_forces_partial_restraint_xy(local_fixed_end_forces, L, k1z, k2z, E, Iz)

    k1z, k2z = InstantFrame.connection_zeros_and_inf(k1z, k2z, E, Iz)

    M_fixed = [local_fixed_end_forces[6], local_fixed_end_forces[12]]
    F_fixed = [local_fixed_end_forces[2], local_fixed_end_forces[8]]
    M_spring, F_spring = fixed_end_forces_partial_restraint(k1z, k2z, E, Iz, L, M_fixed, F_fixed)

    local_fixed_end_forces[6] = M_spring[1]
    local_fixed_end_forces[12] = -M_spring[2]
    local_fixed_end_forces[2] = F_spring[2] #tricky here
    local_fixed_end_forces[8] = F_spring[1] #tricky here
  
    return local_fixed_end_forces

end




function calculate_local_element_fixed_end_forces_partial_restraint_xz(local_fixed_end_forces, L, k1y, k2y, E, Iy)
  
    k1y, k2y = InstantFrame.connection_zeros_and_inf(k1y, k2y, E, Iy)

    M_fixed = [local_fixed_end_forces[5], local_fixed_end_forces[11]]
    F_fixed = [local_fixed_end_forces[3], local_fixed_end_forces[9]]
    M_spring, F_spring = fixed_end_forces_partial_restraint(k1y, k2y, E, Iy, L, M_fixed, F_fixed)

    local_fixed_end_forces[3] = F_spring[1]
    local_fixed_end_forces[5] = -M_spring[1]
    local_fixed_end_forces[9] = F_spring[2]
    local_fixed_end_forces[11] = M_spring[2]
   
    return local_fixed_end_forces

end


function calculate_local_element_fixed_end_forces(wx_local, wy_local, wz_local, L)


    local_fixed_end_forces = zeros(Float64, 12)

    local_fixed_end_forces[1] = -wx_local*L/2
    local_fixed_end_forces[7] = -wx_local*L/2

    local_fixed_end_forces[2] = -wy_local*L/2
    local_fixed_end_forces[6] = -wy_local*L^2/12
    local_fixed_end_forces[8] = -wy_local*L/2
    local_fixed_end_forces[12] = +wy_local*L^2/12


    local_fixed_end_forces[3] = -wz_local*L/2
    local_fixed_end_forces[5] = +wz_local*L^2/12
    local_fixed_end_forces[9] = -wz_local*L/2
    local_fixed_end_forces[11] = -wz_local*L^2/12

    return local_fixed_end_forces

end




function fixed_end_forces_partial_restraint(k1, k2, E, I, L, M_fixed, F_fixed)

    α1 = k1/(E*I/L)
    α2 = k2/(E*I/L)

    #MGZ Example 13.7

    k_moment = E*I/L*([4+α1  2.0
            2.0  4+α2])

    θi = k_moment^-1*M_fixed

    # M_spring = θi .* [k1, k2]

    M_spring = abs.(θi) .* [k1, k2]

    k_shear = (E*I/L)*[6/L  6/L
            -6/L -6/L]

    # F_spring = F_fixed - k_shear * θi
    F_spring = F_fixed + k_shear * θi


    return M_spring, F_spring

end


function calculate_nodal_forces_from_uniform_loads(uniform_load, element, node, element_properties)

    nodal_forces_uniform_load = zeros(Float64, length(node.numbers) * 6)  #hard coded for now

    if !isnothing(uniform_load.elements)

        element_global_nodal_forces_uniform_load = Array{Array{Float64}}(undef, length(uniform_load.elements))
        local_fixed_end_forces = Array{Array{Float64}}(undef, length(uniform_load.elements))

        for i in eachindex(uniform_load.elements)

            elem_num = uniform_load.elements[i]
            elem_index = findfirst(num->num==elem_num, element.numbers)

            global_element_uniform_loads = [uniform_load.magnitudes.qX[i], uniform_load.magnitudes.qY[i], uniform_load.magnitudes.qZ[i], uniform_load.magnitudes.mX[i], uniform_load.magnitudes.mY[i], uniform_load.magnitudes.mZ[i]]

            local_element_uniform_loads = element_properties.Γ[elem_index][1:6, 1:6] * global_element_uniform_loads

            wx_local = local_element_uniform_loads[1]
            wy_local = local_element_uniform_loads[2]
            wz_local = local_element_uniform_loads[3]

            k1z = element_properties.start_connection[elem_index].rz
            k2z = element_properties.end_connection[elem_index].rz
            k1y = element_properties.start_connection[elem_index].ry
            k2y = element_properties.end_connection[elem_index].ry

            local_fixed_end_forces[i] = calculate_local_element_fixed_end_forces(wx_local, wy_local, wz_local, element_properties.L[elem_index])

            #now modify fixed end forces if there are partial end restraints 
            # partial_restraint_index = findall(x->x!=Inf, [k1z, k2z])
            if (k1z < Inf) | (k2z < Inf)
                local_fixed_end_forces[i] = calculate_local_element_fixed_end_forces_partial_restraint_xy(local_fixed_end_forces[i], element_properties.L[elem_index], k1z, k2z, element_properties.E[elem_index], element_properties.Iz[elem_index])
            end

            # partial_restraint_index = findall(x->x!=Inf, [k1y, k2y])
            if (k1y < Inf) | (k2y < Inf)
                local_fixed_end_forces[i] = calculate_local_element_fixed_end_forces_partial_restraint_xz(local_fixed_end_forces[i], element_properties.L[elem_index], k1y, k2y, element_properties.E[elem_index], element_properties.Iy[elem_index])
            end

            global_fixed_end_forces = element_properties.Γ[elem_index]' * local_fixed_end_forces[i]

            element_global_nodal_forces_uniform_load[i] = -global_fixed_end_forces

            nodal_forces_uniform_load[element_properties.global_dof[elem_index]] += element_global_nodal_forces_uniform_load[i]

        end

    else

        element_global_nodal_forces_uniform_load = Array{Array{Float64}}(undef, 0)
        local_fixed_end_forces = Array{Array{Float64}}(undef, 0)

    end

    return element_global_nodal_forces_uniform_load, local_fixed_end_forces, nodal_forces_uniform_load


end


function define_local_elastic_element_torsional_stiffness_matrix_partially_restrained(J, G, L, k1, k2)

    k1, k2 = connection_zeros_and_inf(k1, k2, G, J)

    α1 = k1/(G*J/L)
    α2 = k2/(G*J/L)

    α = 1 + 1/α1 + 1/α2

    kt = G*J/L * [1/α -1/α
                  -1/α 1/α]

    return kt

end

function define_local_elastic_element_axial_stiffness_matrix_partially_restrained(A, E, L, k1, k2)

    k1, k2 = connection_zeros_and_inf(k1, k2, E, A)

    α1 = k1/(A*E/L)
    α2 = k2/(A*E/L)

    α = 1 + 1/α1 + 1/α2

    ka = A*E/L * [1/α -1/α
                  -1/α 1/α]

    return ka

end



function modify_element_local_connection_stiffness(element_properties, ke_local, element)

    for i in eachindex(ke_local)

        start_index = findall(u->u!=Inf, element_properties.start_connection[i])
        end_index = findall(u->u!=Inf, element_properties.end_connection[i])
    
        if (!isempty(start_index)) | (!isempty(end_index))
    
                index = findfirst(elem_num->elem_num==element.numbers[i], element.numbers)
    
                ke_local_xy = InstantFrame.define_local_elastic_element_stiffness_matrix_partially_restrained(element_properties.Iz[index], element_properties.E[index], element_properties.L[index], element_properties.start_connection[index][6], element_properties.end_connection[index][6])
                ke_local_xz = InstantFrame.define_local_elastic_element_stiffness_matrix_partially_restrained(element_properties.Iy[index], element_properties.E[index], element_properties.L[index], element_properties.start_connection[index][5], element_properties.end_connection[index][5])

                #signs are opposite in xz plane for these terms
                ke_local_xz[1, 2] = -ke_local_xz[1, 2]
                ke_local_xz[1, 4] = -ke_local_xz[1, 4]
                ke_local_xz[2, 1] = -ke_local_xz[2, 1]
                ke_local_xz[2, 3] = -ke_local_xz[2, 3]
                ke_local_xz[3, 2] = -ke_local_xz[3, 2]
                ke_local_xz[3, 4] = -ke_local_xz[3, 4]
                ke_local_xz[4, 1] = -ke_local_xz[4, 1]
                ke_local_xz[4, 3] = -ke_local_xz[4, 3]

                ke_local_torsion = define_local_elastic_element_torsional_stiffness_matrix_partially_restrained(element_properties.J[index], element_properties.G[index], element_properties.L[index], element_properties.start_connection[index][4], element_properties.end_connection[index][4])

                ke_local_axial = define_local_elastic_element_axial_stiffness_matrix_partially_restrained(element_properties.A[index], element_properties.E[index], element_properties.L[index], element_properties.start_connection[index][1], element_properties.end_connection[index][1])

                #flexural dof
                dof = [2, 6, 8, 12]

                for i in eachindex(dof)
                    for j in eachindex(dof)
                        ke_local[index][dof[i], dof[j]] = ke_local_xy[i,j]
                    end
                end

                dof = [3, 5, 9, 11]

                for i in eachindex(dof)
                    for j in eachindex(dof)

                        ke_local[index][dof[i], dof[j]] = ke_local_xz[i,j]
                    end
                end

               
                # ke_local[index][dof, dof] .= ke_local_xz

                #torsion dof
                dof = [4, 10]
                for i in eachindex(dof)
                    for j in eachindex(dof)
                        ke_local[index][dof[i], dof[j]] = ke_local_torsion[i,j]
                    end
                end

                dof = [1, 7]
                for i in eachindex(dof)
                    for j in eachindex(dof)
                        ke_local[index][dof[i], dof[j]] = ke_local_axial[i,j]
                    end
                end


                # ke_local[index][dof, dof] .= ke_local_torsion
                
        end
    
    end

    return ke_local

end


function connection_zeros_and_inf(k1, k2, E, I)

    if k1 == Inf  #need big number instead of Inf because 0.0*Inf = NaN
        k1 = E*I*10.0^20
    elseif k1 == 0.0
        k1 = E*I*10.0^-30
    end

    if k2 == Inf #need big number instead of Inf
        k2 = E*I*10.0^20
    elseif k2 == 0.0
        k2 = E*I*10.0^-30
    end

    return k1, k2

end


function define_local_elastic_element_stiffness_matrix_partially_restrained(I, E, L, k1, k2)

    #from McGuire et al. (2000) Example 13.6

    k1, k2 = connection_zeros_and_inf(k1, k2, E, I)
    
    α1 = k1/(E*I/L)  #start node

    α2 = k2/(E*I/L)  #end node

    α = (α1*α2)/(α1*α2 + 4*α1 + 4*α2 + 12)

    ke = zeros(Float64, 4, 4)

    ke[1,1] = 12/L^2*(1 + (α1+α2)/(α1*α2))
    ke[1,2] = 6/L * (1+2/α2)
    ke[1,3] = -12/L^2*(1+(α1+α2)/(α1*α2))
    ke[1,4] = 6/L*(1+2/α1)
    ke[2,2] = 4*(1+3/α2)
    ke[2,3] = -6/L*(1+2/α2)
    ke[2,4] = 2
    ke[3,3] = 12/L^2*(1+(α1+α2)/(α1*α2))
    ke[3,4] = -6/L*(1+2/α1)
    ke[4,4] = 4*(1+3/α1)

    for i = 1:4

        for j = 1:4

            ke[j, i] = ke[i, j]

        end
        
    end

    ke = α * (E * I / L) .* ke

    return ke

end



function define_local_3D_mass_matrix(A, L, Io, ρ)

    #https://pdfs.semanticscholar.org/2361/90b717f2b950078e8c1ddc7da16080d601eb.pdf

    m = zeros(Float64, (12, 12))

    m[1, 1] = 140
    m[1, 7] = 70
    m[2, 2] = 156
    m[2, 6] = 22L
    m[2, 8] = 54
    m[2, 12] = -13L
    m[3, 3] = 156
    m[3, 5] = -22L
    m[3, 9] = 54 
    m[3, 11] = 13L
    m[4, 4] = 140Io/A
    m[4, 10] = 70Io/A
    m[5, 5] = 4L^2 
    m[5, 9] =  -13L
    m[5, 11] = -3L^2
    m[6, 6] = 4L^2
    m[6, 8] = 13L
    m[6, 12] = -3L^2
    m[7, 7] = 140
    m[8, 8] = 156
    m[8, 12] = -22L
    m[9, 9] = 156
    m[9, 11] = 22L
    m[10, 10] = 140Io/A
    m[11, 11] = 4L^2
    m[12, 12] = 4L^2
    
    for i = 1:12

        for j = 1:12

            m[j, i] = m[i, j]

        end
        
    end
    
    m = m .* (ρ*A*L)/420

    return m

end

function define_free_global_dof(node, support)

    num_dof_per_node = 6
    dof_support_stiffness = zeros(Float64, length(node.numbers)*num_dof_per_node)

    nodal_support_stiffness = reduce(hcat, collect(support.stiffness))

    for i in eachindex(support.nodes)

        node_index = findfirst(node_num->node_num == support.nodes[i], node.numbers)
        node_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_index-1)

        dof_support_stiffness[node_dof] = nodal_support_stiffness[i, :]

    end

    free_global_dof = vec(findall(dof->dof!=Inf, dof_support_stiffness))

    fixed_global_dof = vec(findall(dof->dof==Inf, dof_support_stiffness))

    elastic_support_dof = vec(findall(dof->(dof!=Inf)&(dof!=0.0), dof_support_stiffness))

    elastic_support_dof = vec(findall(dof->(dof!=Inf)&(dof!=0.0), dof_support_stiffness))

    elastic_support_stiffness = dof_support_stiffness[elastic_support_dof]

    elastic_support = ElasticSupport(global_dof=elastic_support_dof, global_stiffness = elastic_support_stiffness)

    return free_global_dof, fixed_global_dof, elastic_support

end


function define_global_dof_point_loads(node, point_load)


    num_dof_per_node = 6
    global_dof_point_loads = zeros(Float64, length(node.numbers)*num_dof_per_node)

    if !isnothing(point_load.nodes)

        nodal_point_loads = reduce(hcat, collect(point_load.magnitudes))

        for i in eachindex(point_load.nodes)

            node_index = findfirst(node_num->node_num == point_load.nodes[i], node.numbers)
            node_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_index-1)

            global_dof_point_loads[node_dof] = nodal_point_loads[i, :]

        end

    end

    return global_dof_point_loads

end


function calculate_element_local_axis_directions(Γ)

    unit_vector_X = [1.0, 0.0, 0.0]
    local_x = Γ[1:3,1:3]' * unit_vector_X

    unit_vector_Y = [0.0, 1.0, 0.0]
    local_y = Γ[1:3,1:3]' * unit_vector_Y

    unit_vector_Z = [0.0, 0.0, 1.0]
    local_z = Γ[1:3,1:3]' * unit_vector_Z

    element_local_axes = (x=local_x, y=local_y, z=local_z)

    return  element_local_axes

end

function calculate_element_connection_deformations(properties, element, node, nodal_displacements, element_forces)

    u_element_connection = Array{Array{Float64, 1}}(undef, 0)
    elements_with_connections = Array{Int64}(undef, 0)

    #go through each element with a connection
    for i in eachindex(properties.L)

        start_index = findall(u->u!=Inf, properties.start_connection[i])
        end_index = findall(u->u!=Inf, properties.end_connection[i])

        if (!isempty(start_index)) | (!isempty(end_index))

                u_connection_local = zeros(Float64, 12)

                index = findfirst(elem_num->elem_num==element.numbers[i], element.numbers)

                push!(elements_with_connections, element.numbers[i])

                E = properties.E[index]
                A = properties.A[index]
                I = properties.Iz[index]
                L = properties.L[index]
                G = properties.G[index]
                J = properties.J[index]

                node_i_num = element.nodes[index][1]
                node_j_num = element.nodes[index][2]

                node_i_index = findfirst(num->num==node_i_num, node.numbers)
                node_j_index = findfirst(num->num==node_j_num, node.numbers)

                u_global_element = [nodal_displacements[node_i_index]; nodal_displacements[node_j_index]]
                u_local_element = properties.Γ[index] * u_global_element

                #x-y plane
                k1 = properties.start_connection[index][6]
                k2 = properties.end_connection[index][6]
                k1, k2 = connection_zeros_and_inf(k1, k2, E, I)
                v1 = u_local_element[2]
                θ1 = u_local_element[6]
                v2 = u_local_element[8]
                θ2 = u_local_element[12]
                Mi = element_forces[index][6]
                Mj = element_forces[index][12]
                θij = calculate_connection_rotation(k1, k2, E, I, L, v1, θ1, v2, θ2, Mi, Mj)  
                u_connection_local[6] = θij[1]
                u_connection_local[12] = θij[2]

                #x-z plane
                k1 = properties.start_connection[index][5]
                k2 = properties.end_connection[index][5]
                k1, k2 = connection_zeros_and_inf(k1, k2, E, I)
                v1 = u_local_element[3]
                θ1 = u_local_element[5]
                v2 = u_local_element[9]
                θ2 = u_local_element[11]
                Mi = element_forces[index][5]
                Mj = element_forces[index][11]
                θij = calculate_connection_rotation(k1, k2, E, I, L, v1, θ1, v2, θ2, Mi, Mj)  
                u_connection_local[5] = θij[1]
                u_connection_local[11] = θij[2]

                #torsion 
                k1 = properties.start_connection[index][4]
                k2 = properties.end_connection[index][4]
                k1, k2 = connection_zeros_and_inf(k1, k2, G, J)
                θx1 = u_local_element[4]
                θx2 = u_local_element[10]
                Ti = element_forces[index][4]
                Tj = element_forces[index][10]
                θij = calculate_connection_rotation_torsion(k1, k2, G, J, L, θx1, θx2, Ti, Tj)  
                u_connection_local[4] = θij[1]
                u_connection_local[10] = θij[2]

                #axial
                k1 = properties.start_connection[index][1]
                k2 = properties.end_connection[index][1]
                k1, k2 = connection_zeros_and_inf(k1, k2, E, A)
                Δx1 = u_local_element[1]
                Δx2 = u_local_element[7]
                Pi = element_forces[index][1]
                Pj = element_forces[index][7]
                Δij = calculate_connection_deformation_axial(k1, k2, E, A, L, Δx1, Δx2, Pi, Pj)  
                u_connection_local[1] = Δij[1]
                u_connection_local[7] = Δij[2]

                u_connection_global = properties.Γ[index]' * u_connection_local  #convert from local to global

                push!(u_element_connection, u_connection_global)

        end

    end



    element_connections = ElementConnections(elements_with_connections, u_element_connection)

    return element_connections

end

function calculate_connection_rotation(k1, k2, E, I, L, v1, θ1, v2, θ2, Mi, Mj)

    #Use MZG Example 13.7 equations to solve for rotation at a partially rigid or hinged connection.

    α1 = k1/(E*I/L)
    α2 = k2/(E*I/L)

    kB = (E*I/L)*[4+α1  2
                2   4+α2]

    kC = (E*I/L)*[6/L  -α1 -6/L  0
                  6/L  0    -6/L -α2]

    u = [v1, θ1, v2, θ2]

    M = [Mi, Mj]

    θ = kB^-1*(M - kC*u)

    return θ

end



function calculate_connection_rotation_torsion(k1, k2, G, J, L, θ1, θ2, Ti, Tj)

    α1 = k1/(G*J/L)
    α2 = k2/(G*J/L)

    kB = (G*J/L)*[(1-α1)/2  -1/2
                -1/2   (1-α2)/2]

    kC = (G*J/L)*[α1/2  0
                  0    α2/2]

    u = [θ1, θ2]

    T = [Ti, Tj]

    θ = kB^-1*(T - kC*u)

    return θ

end


function calculate_connection_deformation_axial(k1, k2, E, A, L, Δ1, Δ2, Pi, Pj)

    α1 = k1/(E*A/L)
    α2 = k2/(E*A/L)

    kB = (E*A/L)*[(1-α1)/2  -1/2
                -1/2   (1-α2)/2]

    kC = (E*A/L)*[α1/2  0
                  0    α2/2]

    u = [Δ1, Δ2]

    P = [Pi, Pj]

    Δ = kB^-1*(P - kC*u)

    return Δ

end



end # module



    

#     ke = zeros(Float64, 6, 6)

#     c = [2;3;5;6]

#     ke[c, c] = ke_condensed

#     ke[1, 1] = E*A/L
#     ke[1, 4] = -E*A/L
#     ke[4, 1] = -E*A/L
#     ke[4, 4] = E*A/L


#     return ke

# end


#     # α = (α1 * α2) / (α1*α2 + 4*α1 + 4*α2 + 12)

#     # ke = zeros(Float64, 6, 6)

#     # ke[1, 1] = E*A/L
#     # ke[1, 4] = E*A/L
#     # ke[2, 2] = α*E*I/L * 12/L^2 * (1 + (α1 + α2)/(α1*α2))
#     # ke[2, 3] = α*E*I/L * 6/L * (1 + 2/α2)
#     # ke[2, 5] = α*E*I/L * -12/L^2 * (1 + (α1 + α2)/(α1*α2))
#     # ke[2, 6] = α*E*I/L * 6/L * (1 + 2/α1)
#     # ke[3, 3] = α*E*I/L * 4 * (1 + 3/α2)
#     # ke[3, 5] = α*E*I/L * -6/L * (1 + 2/α2)
#     # ke[3, 6] = α*E*I/L * 2 
#     # ke[5, 5] = α*E*I/L * 12/L^2 * (1 + (α1 + α2)/(α1*α2))
#     # ke[5, 6] = α*E*I/L * -6/L * (1 + 2/α1)
#     # ke[6, 6] = α*E*I/L * 4 * (1 + 3/α1)

#     # for i = 1:6

#     #     for j = 1:6

#     #         ke[j, i] = ke[i, j]

#     #     end
        
#     # end

# #     return ke

# # end


# function define_local_3D_mass_matrix(A, L, ρ)


#     m = zeros(Float64, (12, 12))

#     m[1, 1] = 140
#     m[1, 7] = 70
#     m[2, 2] = 156
#     m[2, 6] = 22L
#     m[2, 8] = 54
#     m[2, 12] = -13L
#     m[3, 3] = 156
#     m[3, 5] = -22L
#     m[3, 9] = 54 
#     m[3, 11] = 13L
#     m[4, 4] = 140Io/A
#     m[4, 10] = 70Io/A
#     m[5, 5] = 4L^2 
#     m[5, 9] =  -13L
#     m[5, 11] = -3L^2
#     m[6, 6] = 4L^2
#     m[6, 8] = 13L
#     m[6, 12] = -3L^2
#     m[7, 7] = 140
#     m[8, 8] = 156
#     m[8, 12] = -22L
#     m[9, 9] = 156
#     m[9, 11] = 22L
#     m[10, 10] = 140Io/A
#     m[11, 11] = 4L^2
#     m[12, 12] = 4L^2
    
#     for i = 1:12

#         for j = 1:12

#             m[j, i] = m[i, j]

#         end
        
#     end
    
#     m = m .* (ρ*A*L)/420

#     return m

# end

