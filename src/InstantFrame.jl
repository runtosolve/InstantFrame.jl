module InstantFrame

using SparseArrays, StaticArrays, LinearAlgebra, Rotations, Parameters, NonlinearSolve

export UI
include("UI.jl")
using .UI


@with_kw struct Material

    names::Array{String}
    E::Array{Float64}
    ν::Array{Float64}
    ρ::Array{Float64}

end

@with_kw struct CrossSection

    names::Array{String}
    A::Array{Float64}
    Iy::Array{Float64}
    Iz::Array{Float64}
    J::Array{Float64}

end

@with_kw struct Connection

    names::Array{String}
    stiffness::NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Vector{Float64}}}

end

@with_kw struct Node

    numbers::Array{Int64}
    coordinates::Array{Tuple{Float64, Float64, Float64}, 1}

end

@with_kw struct Element

    numbers::Array{Int64}
    nodes::Array{Tuple{Int64, Int64}}
    orientation::Array{Float64}
    connections::Array{Tuple{String, String}, 1}
    cross_section::Array{String}
    material::Array{String}

end

@with_kw struct Support

    nodes::Array{Int64}
    stiffness::NamedTuple{(:uX, :uY, :uZ, :rX, :rY, :rZ), NTuple{6, Vector{Float64}}}

end

@with_kw struct UniformLoad

    labels::Union{Array{String}, Nothing}
    elements::Union{Array{Int64}, Nothing}
    loads::Union{NamedTuple{(:qX, :qY, :qZ, :mX, :mY, :mZ), NTuple{6, Vector{Float64}}}, Nothing}

end

UniformLoad(nothing) = UniformLoad(labels=nothing, elements=nothing, loads=nothing)

@with_kw struct PointLoad

    labels::Union{Array{String}, Nothing}
    nodes::Union{Array{Int64}, Nothing}
    loads::Union{NamedTuple{(:FX, :FY, :FZ, :MX, :MY, :MZ), NTuple{6, Vector{Float64}}}, Nothing}

end

PointLoad(nothing) = PointLoad(labels=nothing, nodes=nothing, loads=nothing)

@with_kw struct ElementProperties

    L::Array{Float64}
    A::Array{Float64}
    Iz::Array{Float64}
    Iy::Array{Float64}
    Io::Array{Float64}
    J::Array{Float64}
    E::Array{Float64}
    ν::Array{Float64}
    ρ::Array{Float64}

    Γ::Array{Array{Float64, 2}}
    global_dof::Array{Array{Int64, 1}}

    start_connection::Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}
    end_connection::Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}

end

@with_kw struct FirstOrderEquations
  
    free_dof::Array{Int64}
    ke_local::Array{Array{Float64, 2}}
    ke_global::Array{Array{Float64, 2}}
    Ke::Array{Float64, 2}
  
end

@with_kw struct SecondOrderEquations
  
    free_dof::Array{Int64}
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


@with_kw struct FirstOrderSolution

    u1::Array{Float64}
    P1::Array{Array{Float64, 1}}

end

@with_kw struct SecondOrderSolution

    u1::Array{Float64}
    P1::Array{Array{Float64, 1}}
    u2::Array{Float64}
    P2::Array{Array{Float64, 1}}

end

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


@with_kw struct Model

    properties::ElementProperties
    forces::Union{Forces, Nothing}
    equations::Union{FirstOrderEquations, SecondOrderEquations, ModalEquations}
    solution::Union{FirstOrderSolution, SecondOrderSolution, ModalSolution}

end


function define_local_elastic_stiffness_matrix(Iy, Iz, A, J, E, ν, L)


    G = E/(2*(1+ν))

    ke = zeros(Float64, (12, 12))

    ke[1, 1] = E*A/L
    ke[1, 6] = -E*A/L
    ke[2, 2] = 12*E*Iz/L^3
    ke[2, 6] = 6*E*Iz/L^3
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


function define_rotation_matrix(A, B, β)
    #ρ is x-z plane,  χ is x-y plane , ω is y-z plane
    AB = B - A

    length_AB = norm(AB)

    χ = π/2 - acos(AB[2]/length_AB)

    if (AB[3]==0.0) & (AB[1]==0.0)  #avoid 0/0 case
        ρ = 0.0    
    else
        ρ = -atan(AB[3]/AB[1])
    end

    ω = β

    γ = RotXZY(-ω, -χ, -ρ)

    Γ = zeros(Float64, (12, 12))

    Γ[1:3, 1:3] .= γ
    Γ[4:6, 4:6] .= γ
    Γ[7:9, 7:9] .= γ
    Γ[10:12, 10:12] .= γ

    return Γ

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
    ρ = Array{Float64}(undef, num_elem)
    start_connection = Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}(undef, num_elem)
    end_connection = Array{NamedTuple{(:ux, :uy, :uz, :rx, :ry, :rz), NTuple{6, Float64}}}(undef, num_elem)

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

        #element end conditions
        index = findfirst(name->name == element.connections[i][1], connection.names)
        start_connection[i] = (ux=connection.stiffness.ux[index], uy=connection.stiffness.uy[index], uz=connection.stiffness.uz[index], rx=connection.stiffness.rx[index], ry=connection.stiffness.ry[index], rz=connection.stiffness.rz[index])
        index = findfirst(name->name == element.connections[i][2], connection.names)
        end_connection[i] = (ux=connection.stiffness.ux[index], uy=connection.stiffness.uy[index], uz=connection.stiffness.uz[index], rx=connection.stiffness.rx[index], ry=connection.stiffness.ry[index], rz=connection.stiffness.rz[index])
    
        #element length
        L[i] = norm(node_j - node_i)

        #rotation matrix
        Γ[i] = define_rotation_matrix(node_i, node_j, element.orientation[i])

        #global dof for each element
        num_dof_per_node = 6 #hard code this for now
        node_i_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_i_index-1)
        node_j_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_j_index-1) 
        global_dof[i] = [node_i_dof; node_j_dof]

    end

    element_properties = ElementProperties(L, A, Iz, Iy, Io, J, E, ν, ρ, Γ, global_dof, start_connection, end_connection)

    return element_properties

end


function assemble_global_matrix(k, dof)

    # total_dof = sum([Int(size(k[i], 1)/2) for i in eachindex(k)])
    total_dof = maximum([maximum(dof[i]) for i in eachindex(dof)])

    K = zeros(Float64, (total_dof , total_dof))
    
    for i in eachindex(k)
       K[dof[i], dof[i]] += k[i]
    end

    return K

end

function calculate_element_internal_forces(properties, ke_local, u)

    P_element_local = Array{Array{Float64, 1}}(undef, length(properties.L))

    for i in eachindex(properties.L)
        u_element_global = u[properties.global_dof[i]]
        u_element_local = properties.Γ[i] * u_element_global
        P_element_local[i] = ke_local[i] * u_element_local
    end

    return P_element_local

end


function residual(u, p)
    Kff, Ff = p
    Kff * u - Ff
end

function nonlinear_solution(Kff, Ff, u1f)

    p = [Kff, Ff]

    u1f = SVector{length(u1f)}(u1f)
    # u1f_S = zeros(Float64, length(u1f))
    # u1fSS = [u1f_S[i] for i in eachindex(u1f)]
    probN = NonlinearSolve.NonlinearProblem{false}(residual, u1f, p)
    u2f = NonlinearSolve.solve(probN, NewtonRaphson(), tol = 1e-9)

    return u2f

end


function first_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load)

    element_properties = InstantFrame.define_element_properties(node, cross_section, material, element, connection)

    free_global_dof = InstantFrame.define_free_global_dof(node, support)

    ke_local = [InstantFrame.define_local_elastic_stiffness_matrix(element_properties.Iy[i], element_properties.Iz[i], element_properties.A[i], element_properties.J[i], element_properties.E[i], element_properties.ν[i], element_properties.L[i]) for i in eachindex(element_properties.L)]

    ke_local = InstantFrame.modify_element_local_connection_stiffness(element_properties, ke_local, element)

    ke_global = [element_properties.Γ[i]'*ke_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]

    Ke = InstantFrame.assemble_global_matrix(ke_global, element_properties.global_dof)

    equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load = InstantFrame.calculate_nodal_forces_from_uniform_loads(uniform_load, element, node, element_properties)

    global_dof_point_loads = InstantFrame.define_global_dof_point_loads(node, point_load)

    F = global_dof_point_loads .+ global_dof_nodal_forces_uniform_load

    forces = InstantFrame.Forces(equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load, global_dof_point_loads, F)

    Ff = F[free_global_dof]
    Ke_ff = Ke[free_global_dof, free_global_dof]
    u1f = Ke_ff \ Ff

    u1 = zeros(Float64, size(Ke,1))
    u1[free_global_dof] = u1f

    P1 = InstantFrame.calculate_element_internal_forces(element_properties, ke_local, u1)

    equations = FirstOrderEquations(free_global_dof, ke_local, ke_global, Ke)

    solution = FirstOrderSolution(u1, P1)

    model = Model(element_properties, forces, equations, solution)

    return model

end


function second_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load)

    element_properties = InstantFrame.define_element_properties(node, cross_section, material, element, connection)

    free_global_dof = InstantFrame.define_free_global_dof(node, support)

    ke_local = [InstantFrame.define_local_elastic_stiffness_matrix(element_properties.Iy[i], element_properties.Iz[i], element_properties.A[i], element_properties.J[i], element_properties.E[i], element_properties.ν[i], element_properties.L[i]) for i in eachindex(element_properties.L)]

    ke_local = InstantFrame.modify_element_local_connection_stiffness(element_properties, ke_local, element)

    ke_global = [element_properties.Γ[i]'*ke_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]

    Ke = InstantFrame.assemble_global_matrix(ke_global, element_properties.global_dof)

    equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load = InstantFrame.calculate_nodal_forces_from_uniform_loads(uniform_load, element, node, element_properties)

    global_dof_point_loads = InstantFrame.define_global_dof_point_loads(node, point_load)

    F = global_dof_point_loads .+ global_dof_nodal_forces_uniform_load

    forces = InstantFrame.Forces(equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load, global_dof_point_loads, F)

    Ff = F[free_global_dof]
    Ke_ff = Ke[free_global_dof, free_global_dof]
    u1f = Ke_ff \ Ff

    u1 = zeros(Float64, size(Ke,1))
    u1[free_global_dof] = u1f

    P1 = InstantFrame.calculate_element_internal_forces(element_properties, ke_local, u1)
    P1_axial = [P1[i][7] for i in eachindex(P1)]
    kg_local = [InstantFrame.define_local_geometric_stiffness_matrix(P1_axial[i], element_properties.L[i]) for i in eachindex(element_properties.L)]
    kg_global = [element_properties.Γ[i]'*kg_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]
    Kg = InstantFrame.assemble_global_matrix(kg_global, element_properties.global_dof)
    Kg_ff = Kg[free_global_dof, free_global_dof]

    Kff = Ke_ff + Kg_ff

    p = [Kff, Ff]
    u1f = SVector{length(u1f)}(u1f)
    probN = NonlinearSolve.NonlinearProblem{false}(residual, u1f, p)
    u2f = NonlinearSolve.solve(probN, NewtonRaphson(), tol = 1e-9)

    u2 = zeros(Float64, size(Ke,1))
    u2[free_global_dof] = u2f
    P2 = InstantFrame.calculate_element_internal_forces(element_properties, ke_local, u2)

    equations = SecondOrderEquations(free_global_dof, ke_local, ke_global, Ke, kg_local, kg_global, Kg)

    solution = SecondOrderSolution(u1, P1, u2, P2)

    model = Model(element_properties, forces, equations, solution)

    return model

end



function modal_vibration_analysis(node, cross_section, material, connection, element, support)

    element_properties = InstantFrame.define_element_properties(node, cross_section, material, element, connection)

    free_global_dof = InstantFrame.define_free_global_dof(node, support)

    ke_local = [InstantFrame.define_local_elastic_stiffness_matrix(element_properties.Iy[i], element_properties.Iz[i], element_properties.A[i], element_properties.J[i], element_properties.E[i], element_properties.ν[i], element_properties.L[i]) for i in eachindex(element_properties.L)]

    ke_local = InstantFrame.modify_element_local_connection_stiffness(element_properties, ke_local, element)

    ke_global = [element_properties.Γ[i]'*ke_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]

    Ke = InstantFrame.assemble_global_matrix(ke_global, element_properties.global_dof)

    m_local = [InstantFrame.define_local_3D_mass_matrix(element_properties.A[i], element_properties.L[i], element_properties.Io[i], element_properties.ρ[i]) for i in eachindex(element_properties.L)]
    m_global = [element_properties.Γ[i]'*m_local[i]*element_properties.Γ[i] for i in eachindex(element_properties.L)]
    M = InstantFrame.assemble_global_matrix(m_global, element_properties.global_dof)
    
    Ke_ff = Ke[free_global_dof, free_global_dof]
    Mff = M[free_global_dof, free_global_dof]

    ωn_squared = eigvals(Ke_ff, Mff)
    ωn = sqrt.(ωn_squared)
    ϕff = eigvecs(Ke_ff, Mff)

    ϕ = [zeros(Float64, size(M, 1)) for i in eachindex(ωn)]

    for i in eachindex(ωn)
        ϕ[i][free_global_dof] .= ϕff[:, i]
    end


    equations = ModalEquations(free_global_dof, ke_local, ke_global, Ke, m_local, m_global, M)

    solution = ModalSolution(ωn, ϕ)

    model = Model(element_properties, nothing, equations, solution)

    return model

end


function solve(node, cross_section, material, connection, element, support, uniform_load, point_load; analysis_type)

    if analysis_type == "first order"

        first_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load)

    elseif analysis_type == "second order"

        second_order_analysis(node, cross_section, material, connection, element, support, uniform_load, point_load)

    elseif analysis_type == "modal vibration"

        modal_vibration_analysis(node, cross_section, material, connection, element, support)

    end

end


function calculate_local_element_fixed_end_forces(wx_local, wy_local, wz_local, L)

    local_fixed_end_forces = zeros(Float64, 12)

    local_fixed_end_forces[1] = -wx_local*L/2
    local_fixed_end_forces[7] = -wx_local*L/2

    local_fixed_end_forces[3] = -wy_local*L/2
    local_fixed_end_forces[6] = -wy_local*L^2/12
    local_fixed_end_forces[9] = -wy_local*L/2
    local_fixed_end_forces[12] = wy_local*L^2/12

    local_fixed_end_forces[2] = -wz_local*L/2
    local_fixed_end_forces[5] = wz_local*L^2/12
    local_fixed_end_forces[8] = -wz_local*L/2
    local_fixed_end_forces[11] = -wz_local*L^2/12

    return local_fixed_end_forces

end


function calculate_nodal_forces_from_uniform_loads(uniform_load, element, node, element_properties)

    element_global_nodal_forces_uniform_load = Array{Array{Float64}}(undef, length(uniform_load.elements))
    local_fixed_end_forces = Array{Array{Float64}}(undef, length(uniform_load.elements))
    nodal_forces_uniform_load = zeros(Float64, length(node.numbers) * 6)  #hard coded for now

    if !isnothing(uniform_load.elements)

        for i in eachindex(uniform_load.elements)

            elem_num = uniform_load.elements[i]
            elem_index = findfirst(num->num==elem_num, element.numbers)

            global_element_uniform_loads = [uniform_load.loads.qX[elem_index], uniform_load.loads.qY[elem_index], uniform_load.loads.qZ[elem_index], uniform_load.loads.mX[elem_index], uniform_load.loads.mY[elem_index], uniform_load.loads.mZ[elem_index]]

            local_element_uniform_loads = element_properties.Γ[elem_index][1:6, 1:6] * global_element_uniform_loads

            wx_local = local_element_uniform_loads[1]
            wy_local = local_element_uniform_loads[2]
            wz_local = local_element_uniform_loads[3]

            local_fixed_end_forces[i] = calculate_local_element_fixed_end_forces(wx_local, wy_local, wz_local, element_properties.L[elem_index])

            global_fixed_end_forces = element_properties.Γ[elem_index]' * local_fixed_end_forces[i]

            element_global_nodal_forces_uniform_load[i] = -global_fixed_end_forces

            nodal_forces_uniform_load[element_properties.global_dof[elem_index]] += element_global_nodal_forces_uniform_load[i]

        end

    end

    return element_global_nodal_forces_uniform_load, local_fixed_end_forces, nodal_forces_uniform_load

end

function modify_element_local_connection_stiffness(element_properties, ke_local, element)

    for i in eachindex(ke_local)

        start_index = findall(u->u!=Inf, element_properties.start_connection[i])
        end_index = findall(u->u!=Inf, element_properties.end_connection[i])
    
        if (!isempty(start_index)) | (!isempty(end_index))
    
                index = findfirst(elem_num->elem_num==element.numbers[i], element.numbers)
    
                ke_local_xy = InstantFrame.define_local_elastic_element_stiffness_matrix_partially_restrained(element_properties.Iz[index], element_properties.E[index], element_properties.L[index], element_properties.start_connection[index][6], element_properties.end_connection[index][6])
                ke_local_xz = InstantFrame.define_local_elastic_element_stiffness_matrix_partially_restrained(element_properties.Iy[index], element_properties.E[index], element_properties.L[index], element_properties.start_connection[index][5], element_properties.end_connection[index][5])
    
                dof = [2, 6, 8, 12]
                ke_local[index][dof, dof] .= ke_local_xy
    
                dof = [3, 5, 9, 11]
                ke_local[index][dof, dof] .= ke_local_xz
    
        end
    
    end

    return ke_local

end


function define_local_elastic_element_stiffness_matrix_partially_restrained(I, E, L, k1, k2)

    #from McGuire et al. (2000) Example 13.6
    
    α1 = k1/(E*I/L)  #start node

    α2 = k2*(E*I/L)  #end node

    Kbb = E*I/L * [4+α1     2
                   2        4+α2]

    Kbc = E*I/L *[6/L   -α1     -6/L    0.0
                  6/L   0.0     -6/L    -α2]

    Kcb = Kbc'

    Kcc = E*I/L * [12/L^2   0.0     -12/L^2     0.0
                   0.0      α1      0.0         0.0
                   -12/L^2  0.0     12/L^2      0.0
                   0.0      0.0     0.0         α2]

    
    ke = Kcc - Kcb * Kbb^-1 * Kbc

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

    return free_global_dof

end


function define_global_dof_point_loads(node, point_load)


    num_dof_per_node = 6
    global_dof_point_loads = zeros(Float64, length(node.numbers)*num_dof_per_node)

    if !isnothing(point_load.nodes)

        nodal_point_loads = reduce(hcat, collect(point_load.loads))

        for i in eachindex(point_load.nodes)

            node_index = findfirst(node_num->node_num == point_load.nodes[i], node.numbers)
            node_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_index-1)

            global_dof_point_loads[node_dof] = nodal_point_loads[i, :]

        end

    end

    return global_dof_point_loads

end


function define_global_element_displacements(u, global_dof)

    u_global_e = [zeros(Float64, 12) for i in eachindex(global_dof)]

    for i in eachindex(global_dof)

        u_global_e[i] = u[global_dof[i]]
        
    end

    return u_global_e

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


function define_local_3D_mass_matrix(A, L, ρ)


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