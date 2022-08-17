using SparseArrays, StaticArrays, LinearAlgebra, Rotations, Parameters, NonlinearSolve

# export UI
# include("UI.jl")
# using .UI


@with_kw struct Node

    list::Array{Int64}
    coordinates::Array{Array{Float64}}

end

@with_kw struct CrossSection

    names::Array{String}
    A::Array{Float64}
    Iy::Array{Float64}
    Iz::Array{Float64}
    J::Array{Float64}

end

@with_kw struct Material

    names::Array{String}
    E::Array{Float64}
    ν::Array{Float64}

end

@with_kw struct Connection

    names::Array{String}
    rx::Array{Float64}
    ry::Array{Float64}
    rz::Array{Float64} 

end

@with_kw struct Element

    list::Array{Int64}
    start_node::Array{Int64}
    end_node::Array{Int64}
    β::Array{Float64}
    start_connection::Array{String}
    end_connection::Array{String}
    cross_section::Array{String}
    material::Array{String}

end

@with_kw struct Properties

    L::Array{Float64}
    A::Array{Float64}
    Iz::Array{Float64}
    Iy::Array{Float64}
    J::Array{Float64}
    E::Array{Float64}
    ν::Array{Float64}

    Γ::Array{Array{Float64, 2}}
    global_dof::Array{Array{Int64, 1}}

end

@with_kw struct Equations
  
    free_dof::Array{Int64}
    ke_local::Array{Array{Float64, 2}}
    ke_global::Array{Array{Float64, 2}}
    Ke::Array{Float64, 2}
    kg_local::Array{Array{Float64, 2}}
    kg_global::Array{Array{Float64, 2}}
    Kg::Array{Float64, 2}  
  
end

@with_kw struct Solution

    u1::Array{Float64}
    P1::Array{Float64}
    u2::Array{Float64}
    P2::Array{Float64}

end


@with_kw struct Model

    properties::Properties
    equations::Equations
    solution::Solution

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


function define_element_properties(node, cross_section, material, element)

    num_elem = length(element.list)

    L = Array{Float64}(undef, num_elem)
    Γ = Array{Array{Float64, 2}}(undef, num_elem)
    global_dof = Array{Array{Int64, 1}}(undef, num_elem)
    Iz = Array{Float64}(undef, num_elem)
    Iy = Array{Float64}(undef, num_elem)
    A = Array{Float64}(undef, num_elem)
    J = Array{Float64}(undef, num_elem)
    E = Array{Float64}(undef, num_elem)
    ν = Array{Float64}(undef, num_elem)


    for i in eachindex(element.list)

        #nodal coordinates
        node_i_index = findfirst(node_num->node_num == element.start_node[i], node.list)
        node_j_index = findfirst(node_num->node_num == element.end_node[i], node.list)
        node_i = node.coordinates[node_i_index]
        node_j = node.coordinates[node_j_index]

        #cross section
        index = findfirst(name->name == element.cross_section[i], cross_section.names)
        Iz[i] = cross_section.Iz[index]
        Iy[i] = cross_section.Iy[index]
        A[i] = cross_section.A[index]
        J[i] = cross_section.J[index]
        
        #material
        index = findfirst(name->name == element.material[i], material.names)
        E[i] = material.E[index]
        ν[i] = material.ν[index]
    
        #element length
        L[i] = norm(node_j - node_i)
    
        #global dof for each element
        num_dof_per_node = Int(size(ke_local[i], 1)/2)
        node_i_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_i_index-1)
        node_j_dof = range(1, num_dof_per_node) .+ num_dof_per_node * (node_j_index-1) 
        global_dof[i] = [node_i_dof; node_j_dof]

    end

    properties = FrameProperties(L, A, Iz, Iy, J, E, ν, Γ, global_dof)

    return properties

end


function assemble_stiffness_matrix(k, dof)

    total_dof = sum([Int(size(k[i], 1)/2) for i in eachindex(k)])

    K = zeros(Float64, (total_dof , total_dof))
    
    for i in eachindex(k)
       K[dof[i], dof[i]] += k[i]
    end

    return K

end

function calculate_element_internal_forces(properties, u)

    P_element_local = Array{Array{Float64, 1}}(undef, length(properties.L))

    for i in eachindex(properties.L)
        u_element_global = u[properties.global_dof[i]]
        u_element_local = properties.Γ[i] * u_element_global
        P_element_local[i] = properties.ke_local[i] * u_element_local
    end

    return P_element_local

end


function residual(u, p)
    Kff, Ff = p
    Kff * u - Ff
end

function nonlinear_solution(Kff, Ff, u1f)

    p = [Kff, Ff]
    u1f = @SVector [u1f[i] for i in eachindex(u1f)]
    probN = NonlinearSolve.NonlinearProblem{false}(residual, u1fs, p)
    u2f = NonlinearSolve.solve(probN, NewtonRaphson(), tol = 1e-9)

    return u2f

end


function solve(node, cross_section, material, element, supports, loads)

    properties = InstantFrame.define_element_properties(node, cross_section, material, element)

    free_dof = vec(findall(dof->dof==0, supports))

    ke_local = [InstantFrame.InstantFrame.define_local_elastic_stiffness_matrix(properties.Iy[i], properties.Iz[i], properties.A[i], properties.J[i], properties.E[i], properties.ν[i], properties.L[i]) for i in eachindex(properties.L)]
    ke_global = [properties.Γ[i]'*ke_local[i]*properties.Γ[i] for i in eachindex(properties.L)]

    Ke = InstantFrame.assemble_stiffness_matrix(ke_global, properties.global_dof)
        
    Ff = loads[free_dof]
    Ke_ff = Ke[free_dof, free_dof]
    u1f = Ke_ff \ Ff

    u1 = zeros(Float64, size(Ke,1))
    u1[free_dof] = u1f

    P1 = InstantFrame.calculate_element_internal_forces(properties, u1)
    P1_axial = [P_local_1[i][7] for i in eachindex(P1)]
    kg_local = [InstantFrame.define_local_geometric_stiffness_matrix(P1_axial[i], properties.L[i]) for i in eachindex(properties.L)]
    kg_global = [properties.Γ[i]'*kg_local[i]*properties.Γ[i] for i in eachindex(properties.L)]
    Kg = InstantFrame.assemble_stiffness_matrix(kg_global, properties.global_dof)
    Kg_ff = Kg[free_dof, free_dof]

    Kff = Ke_ff + Kg_ff

    # u2f = nonlinear_solution(Kff, Ff, u1f)

    # u2 = zeros(Float64, size(Ke,1))
    # u2[free_dof] = u2f
    # P2 = InstantFrame.calculate_element_internal_forces(properties, u2)

    # equations = Equations(free_dof, ke_local, ke_global, Ke, kg_local, kg_global, Kg)

    # solution = Solution(u1, P1, u2, P2)

    # model = Model(properties, equations, solution)

    return Kff

end