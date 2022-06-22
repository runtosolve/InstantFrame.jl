module InstantFrame

struct Nodes

    x::Array{Float64}
    y::Array{Float64}
    z::Array{Float64}

end

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
    start_connection::Array{String}
    end_connection::Array{String}
    cross_section::Array{String}
    material::Array{String}

end


function define_local_elastic_stiffness_matrix(I, A, E, L)


    ke = [E*A/L   0.          0.          E*A/L    0.          0.
         0.     12E*I/L^3    6E*I/L^2     0.      -12E*I/L^3   6E*I/L^2
         0.     6E*I/L^2     4E*I/L       0.      -6E*I/L^2    2E*I/L
         -E*A/L  0.          0.          E*A/L    0.          0.
         0.     -12E*I/L^3   -6E*I/L^2    0.      12E*I/L^3    -6E*I/L^2
         0.     6E*I/L^2     2E*I/L       0.      -6E*I/L^2    4E*I/L]

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



function define_global_transformation(βo)

    T = [ cos(βo)     sin(βo)     0.  0.          0.      0.
                    -sin(βo)    cos(βo)     0.  0.          0.      0.
                    0.          0.          1.  0.          0.      0.
                    0.          0.          0.  cos(βo)     sin(βo) 0.
                    0.          0.          0.  -sin(βo)    cos(βo) 1.]
                    
    return T

end


function second_order_analysis_residual!(R, U, K, F)

    for i=1:length(F)
 
       R[i] = transpose(K[i,:]) * U - F[i]
 
    end
 
    return R
 
 end

 function beam_shape_function(q1,q2,q3,q4,L,x, offset)

    a0=q1
    a1=q2
    a2=1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
    a3=1/L^3*(2*q1+q2*L-2*q3+q4*L)

    w = a0 .+a1.*x .+a2.*x.^2 .+a3.*x.^3 .+ offset

end


function define_local_element_mass_matrix(A, L, ρ)

    m=zeros(Float64, (6,6))

    m[1,1] = 140
    m[1,4] = 70
    m[2,2] = 156
    m[2,3] = 22L
    m[2,5] = 54
    m[2,6] = -13L
    m[3,3] = 4L^2
    m[3,5] = 13L
    m[3,6] = -3L^2
    m[4,4] = 140
    m[5,5] = 156
    m[5,6] = -22L
    m[6,6] = 4L^2

    m[2:6,1] = m[1,2:6]
    m[3:6,2] = m[2,3:6]
    m[4:6,3] = m[3,4:6]
    m[5:6,4] = m[4,5:6]
    m[6,5] = m[5,6]

    m =(ρ*A*L/420) * m

    return m

end



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

# function define_local_transformation(δΔy, Lo)

#     T_to_local = [ -1.  δΔy/Lo  0.  1.  -δΔy/Lo.    0.
#                     0.  1/Lo    1.  0.  -1/Lo       0.
#                     0.  1/Lo    1.  0.  1/Lo        1.]


#     return T_to_local

# end

# function define_geometric_stiffness_matrix(P, Lo)

#     Kg = zeros(Float64, (6, 6))

#     Kg[2, 2] = P/Lo
#     Kg[2, 4] = -P/Lo
#     Kg[5, 2] = -P/Lo
#     Kg[5, 4] =  P/Lo

#     return Kg

# end


#calculate element orientations

#define element elastic stiffness matrices

#define element geometric stiffness matrices

#define element mass matrices

#apply member end conditions

#calculate element global stiffness matrices

#assemble Ke and Kg

#assemble M

#apply boundary conditions

#solve for displacement field

#solve for internal forces

#solve for vibration modes

end # module
