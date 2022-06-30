module InstantFrame

using SparseArrays, StaticArrays, LinearAlgebra, Rotations

# struct Nodes

#     x::Array{Float64}
#     y::Array{Float64}
#     z::Array{Float64}

# end

struct CrossSections

    name::Array{String}
    A::Array{Float64}
    Ixx::Array{Float64}
    Iyy::Array{Float64}
    J::Array{Float64}

end

struct Materials

    name::Array{String}
    E::Array{Float64}
    ν::Array{Float64}

end

struct Connections

    name::Array{String}
    rx::Array{Float64}
    ry::Array{Float64}
    rz::Array{Float64} 

end

struct Members

    type::Array{String}
    start_node::Array{Int64}
    end_node::Array{Int64}
    β::Array{Float64}
    start_connection::Array{String}
    end_connection::Array{String}
    cross_section::Array{String}
    material::Array{String}
  
end


function define_local_elastic_element_stiffness_matrix(I, A, E, L)


    ke = [E*A/L   0.          0.          -E*A/L    0.          0.
         0.     12E*I/L^3    6E*I/L^2     0.      -12E*I/L^3   6E*I/L^2
         0.     6E*I/L^2     4E*I/L       0.      -6E*I/L^2    2E*I/L
         -E*A/L  0.          0.          E*A/L    0.          0.
         0.     -12E*I/L^3   -6E*I/L^2    0.      12E*I/L^3    -6E*I/L^2
         0.     6E*I/L^2     2E*I/L       0.      -6E*I/L^2    4E*I/L]

    return ke

end


function define_local_elastic_stiffness_matrix(Ix, Iy, A, J, E, ν, L)


    G = E/(2*(1+ν))

    ke = zeros(Float64, (12, 12))

    ke[1, 1] = E*A/L
    ke[1, 6] = -E*A/L
    ke[2, 2] = 12*E*Iy/L^3
    ke[2, 6] = 6*E*Iy/L^3
    ke[2, 8] = -12*E*Iy/L^3
    ke[2, 12] = 6*E*Iy/L^2
    ke[3, 3] = 12*E*Ix/L^3
    ke[3, 5] = -6*E*Ix/L^2
    ke[3, 9] = -12*E*Ix/L^3
    ke[3, 11] = -6*E*Ix/L^2
    ke[4, 4] = G*J/L
    ke[4, 10] = -G*J/L
    ke[5, 5] = 4*E*Ix/L
    ke[5, 9] = 6*E*Ix/L^2
    ke[5, 11] = 2*E*Ix/L
    ke[6, 6] = 4*E*Iy/L
    ke[6, 8] = -6*E*Iy/L^2
    ke[6, 12] = 2*E*Iy/L
    ke[7, 7] = E*A/L
    ke[8, 8] = 12*E*Iy/L^3
    ke[8, 12] = -6*E*Iy/L^2
    ke[9, 9] = 12*E*Ix/L^3
    ke[9, 11] = 6*E*Ix/L^2
    ke[10, 10] = G*J/L
    ke[11, 11] = 4*E*Ix/L
    ke[12, 12] = 4*E*Iy/L

    for i = 1:12

        for j = 1:12

            ke[j, i] = ke[i, j]

        end
        
    end

    return ke

end



function define_local_geometric_stiffness_matrix(P, L)


    kg = P/L * [0. 0.     0.          0.  0.      0.
                0. 6/5    L/10        0.  -6/5    L/10
                0. L/10   2/15*L^2    0.  -L/10   -L^2/30
                0. 0.     0.          0.  0.      0.
                0. -6/5   -L/10       0.  6/5     -L/10
                0. L/10   -L^2/30     0. -L/10    2/15*L^2]

    return kg

end

function define_local_3D_geometric_stiffness_matrix(P, L)


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




# function define_global_transformation(βo)

#     T = [ cos(βo)     sin(βo)     0.  0.          0.      0.
#                     -sin(βo)    cos(βo)     0.  0.          0.      0.
#                     0.          0.          1.  0.          0.      0.
#                     0.          0.          0.  cos(βo)     sin(βo) 0.
#                     0.          0.          0.  -sin(βo)    cos(βo) 1.]
                    
#     return T

# end


function define_rotation_matrix(A, B, β)

    AB = B - A

    length_AB = norm(AB)

    χ = π/2 - acos(AB[2]/length_AB)
    ρ = -atan(AB[3]/AB[1])
    ω = β

    γ = RotXZY(-ω, -χ, -ρ)

    Γ = zeros(Float64, (12, 12))

    Γ[1:3, 1:3] .= γ
    Γ[4:6, 4:6] .= γ
    Γ[7:9, 7:9] .= γ
    Γ[10:12, 10:12] .= γ

    return Γ

end



function define_rotation_matrix(A, B)

    AB = B - A

    χ = atan(AB[2], AB[1])

    γ = Angle2d(-χ)  #need negative sign here since z-axis is pointing in opposite direction in Rotations.jl

    Γ = zeros(Float64, (6, 6))

    Γ[1:2, 1:2] .= γ
    Γ[3, 3] = 1.0
    Γ[4:5, 4:5] .= γ
    Γ[6, 6] = 1.0

    return Γ

end



# function second_order_analysis_residual!(R, U, K, F)

#     for i=1:length(F)
 
#        R[i] = transpose(K[i,:]) * U - F[i]
 
#     end
 
#     return R
 
#  end

 function beam_shape_function(q1,q2,q3,q4,L,x, offset)

    a0=q1
    a1=q2
    a2=1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
    a3=1/L^3*(2*q1+q2*L-2*q3+q4*L)

    w = a0 .+a1.*x .+a2.*x.^2 .+a3.*x.^3 .+ offset

end


# function define_local_element_mass_matrix(A, L, ρ)

#     m=zeros(Float64, (6,6))

#     m[1,1] = 140
#     m[1,4] = 70
#     m[2,2] = 156
#     m[2,3] = 22L
#     m[2,5] = 54
#     m[2,6] = -13L
#     m[3,3] = 4L^2
#     m[3,5] = 13L
#     m[3,6] = -3L^2
#     m[4,4] = 140
#     m[5,5] = 156
#     m[5,6] = -22L
#     m[6,6] = 4L^2

#     m[2:6,1] = m[1,2:6]
#     m[3:6,2] = m[2,3:6]
#     m[4:6,3] = m[3,4:6]
#     m[5:6,4] = m[4,5:6]
#     m[6,5] = m[5,6]

#     m =(ρ*A*L/420) * m

#     return m

# end



function get_load_deformation_response(I, A, E, L, P, F, free_dof)

    #Define the frame elastic stiffness.
    Ke = InstantFrame.define_local_elastic_stiffness_matrix(I, A, E, L)
    Kg = InstantFrame.define_local_geometric_stiffness_matrix(P, L)

    #Partition external force vector and elastic stiffness matrix.
    Ff = F[free_dof]
    Ke_ff = Ke[free_dof, free_dof]
    Kg_ff = Kg[free_dof, free_dof]

    #Get the total stiffness matrix, elastic + geometric.
    Kff = Ke_ff + Kg_ff

    #Define the deformation initial guess for the nonlinear solver.
    deformation_guess = Ke_ff \ Ff

    #Solve for the beam deformations.
    solution = nlsolve((R,U) ->InstantFrame.second_order_analysis_residual!(R, U, Kff, Ff), deformation_guess)

    #Get the frame deformed shape.
    u = zeros(Float64, 6)
    u[free_dof] .= solution.zero

    return u

end

function define_global_elastic_stiffness_matrix(nodes, cross_sections, materials, members)

    num_dof_per_node = 3
    num_nodes = length(nodes)
    Ke = zeros(Float64, (num_dof_per_node * num_nodes, num_dof_per_node * num_nodes))
    
    num_elements = length(members.type)
    
    for i=1:num_elements
    
        node_i = members.start_node[i]
        node_j = members.end_node[i]
    
        cross_section_label = members.cross_section[i]
        index = findfirst(name->name == cross_section_label, cross_sections.name)
    
        I = cross_sections.Ixx[index]
        A = cross_sections.A[index]
    
        material_label = members.material[i]
        index = findfirst(name->name == material_label, materials.name)
    
        E = materials.E[index]
        ν = materials.ν[index]
    
        L = norm(nodes[node_j] - nodes[node_i])
    
        ke_local = InstantFrame.define_local_elastic_element_stiffness_matrix(I, A, E, L)
    
        A = nodes[node_i]
        B = nodes[node_j]
        Γ = InstantFrame.define_rotation_matrix(A, B)
    
        ke_global = Γ' * ke_local * Γ
    
        node_i_dof = range(1, num_dof_per_node) .+ (node_i - 1) * num_dof_per_node
        node_j_dof = range(1, num_dof_per_node) .+ (node_j - 1) * num_dof_per_node
    
        global_dof = [node_i_dof; node_j_dof]
    
        Ke[global_dof, global_dof] = Ke[global_dof, global_dof] .+ ke_global
    
    end

    return Ke

end


function calculate_internal_forces(nodes, cross_sections, materials, members, u)

    num_dof_per_node = 3
    num_elements = length(members.type)
    P = Array{Array{Float64}}(undef, num_elements)
    
    for i=1:num_elements
    
        node_i = members.start_node[i]
        node_j = members.end_node[i]
    
        cross_section_label = members.cross_section[i]
        index = findfirst(name->name == cross_section_label, cross_sections.name)
    
        I = cross_sections.Ixx[index]
        A = cross_sections.A[index]
    
        material_label = members.material[i]
        index = findfirst(name->name == material_label, materials.name)
    
        E = materials.E[index]
        ν = materials.ν[index]
    
        L = norm(nodes[node_j] - nodes[node_i])
    
        ke_local = InstantFrame.define_local_elastic_element_stiffness_matrix(I, A, E, L)

        A = nodes[node_i]
        B = nodes[node_j]
        Γ = InstantFrame.define_rotation_matrix(A, B)

        node_i_dof = range(1, num_dof_per_node) .+ (node_i - 1) * num_dof_per_node
        node_j_dof = range(1, num_dof_per_node) .+ (node_j - 1) * num_dof_per_node

        global_dof = [node_i_dof; node_j_dof]

        ue_local = Γ * u[global_dof]
    
        P[i] = ke_local * ue_local
    
    end

    return P
    
end
    

    function define_global_geometric_stiffness_matrix(nodes, members, P)

        num_dof_per_node = 3
        num_nodes = length(nodes)
        Kg = zeros(Float64, (num_dof_per_node * num_nodes, num_dof_per_node * num_nodes))
        
        num_elements = length(members.type)
        
        for i=1:num_elements
        
            node_i = members.start_node[i]
            node_j = members.end_node[i]
            
            L = norm(nodes[node_j] - nodes[node_i])
        
            element_axial_force = P[i][4]  #this needs to change for 3D!!! 
            kg_local = define_local_geometric_stiffness_matrix(element_axial_force, L)
        
            A = nodes[node_i]
            B = nodes[node_j]
            Γ = InstantFrame.define_rotation_matrix(A, B)
        
            kg_global = Γ' * kg_local * Γ
        
            node_i_dof = range(1, num_dof_per_node) .+ (node_i - 1) * num_dof_per_node
            node_j_dof = range(1, num_dof_per_node) .+ (node_j - 1) * num_dof_per_node
        
            global_dof = [node_i_dof; node_j_dof]
        
            Kg[global_dof, global_dof] = Kg[global_dof, global_dof] .+ kg_global
        
        end
    
        return Kg
    
    end


end # module
