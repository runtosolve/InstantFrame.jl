using InstantFrame 


#Apply a uniform load to a simply-supported frame element.

material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[0.627], Iy=[0.627], Iz=[0.627], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "end release"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 230000.0]))

node = InstantFrame.Node(numbers=[1, 2, 3, 4, 5, 6], coordinates=[(0.0, 0.0, 0.0), (96.0, 0.0, 0.0), (0.0, 48.0, 0.0), (0.0, -48.0, 0.0), (96.0, 48.0, 0.0), (96.0, -48.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2, 3, 4, 5], nodes=[(1,2), (1, 3), (1, 4), (2, 5), (2, 6)], orientation=[0.0, 0.0, 0.0, 0.0, 0.0], connections=[("end release", "end release"), ("rigid", "rigid"), ("rigid", "rigid"), ("rigid", "rigid"), ("rigid", "rigid")], cross_section=["beam", "beam", "beam", "beam", "beam"], material=["steel", "steel", "steel", "steel", "steel"])

support = InstantFrame.Support(nodes=[3, 4, 5, 6], stiffness=(uX=[Inf,Inf,Inf,Inf], uY=[Inf,Inf,Inf,Inf], uZ=[Inf,Inf,Inf,Inf], rX=[Inf,Inf,Inf,Inf], rY=[Inf,Inf,Inf,Inf], rZ=[Inf,Inf,Inf,Inf]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1], loads=(qX=[0.0], qY=[-12.5], qZ=[0.0], mX=[0.0], mY=[0.0], mZ=[0.0]))

# uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(nothing)

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

# scale = (1.0, 1.0, 1.0)
# figure = InstantFrame.UI.display_model_deformed_shape(model.solution.u1, element, node, model.properties, scale)




# E*I/L


# ke_pr = InstantFrame.define_local_elastic_element_stiffness_matrix_partially_restrained(I, E, L, k1, k2)

# dof = [2; 6; 8; 12]

# ke = model.equations.ke_local[1][dof, dof]

# free_dof = [1; 3]
# u = [1.0, -1.0]


# ke_ff = ke[free_dof, free_dof]

# M_fixed = ke_ff * u

# ke_pr_ff = ke_pr[free_dof, free_dof]
# M_pr = ke_pr_ff * u


# 41.5118   1992.56         -41.5118   1992.56
# 1992.56        1.67198e5  -1992.56    24088.3
#  -41.5118  -1992.56          41.5118  -1992.56
# 1992.56    24088.3        -1992.56        1.67198e5



# 246.623   11837.9          -246.623   11837.9
# 11837.9    757625.0        -11837.9         3.78812e5
#  -246.623  -11837.9           246.623  -11837.9
# 11837.9         3.78812e5  -11837.9    757625.0


I = 0.627
E = 29000000.0
L = 96.0
k1 = 230000.0
k2 = 2000.0



M_fixed = [9600.0, -9600.0]
F_fixed = [600.0, 600.0]      

function fixed_end_forces_partial_restraint(k1, k2, E, I, L, M_fixed, F_fixed)

    α1 = k1/(E*I/L)
    α2 = k2/(E*I/L)

    #MGZ Example 13.7

    k_moment = E*I/L*([4+α1  2.0
            2.0  4+α2])

    θi = k_moment^-1*M_fixed

    M_spring = θi .* [k1, k2]

    k_shear = (E*I/L)*[6/L  6/L
            -6/L -6/L]

    F_spring = F_fixed - k_shear * θi

    return M_spring, F_spring

end


wx_local = 0.0
wy_local = -12.5
wz_local = -12.5

k1z = 230000.0
k2z = 2000.0
k1y = 230000.0
k2y = 230000.0

local_fixed_end_forces = zeros(Float64, 12)

    local_fixed_end_forces[1] = -wx_local*L/2
    local_fixed_end_forces[7] = -wx_local*L/2

    # local_fixed_end_forces[3] = -wy_local*L/2
    # local_fixed_end_forces[6] = -wy_local*L^2/12
    # local_fixed_end_forces[9] = -wy_local*L/2
    # local_fixed_end_forces[12] = wy_local*L^2/12

    local_fixed_end_forces[2] = -wy_local*L/2
    local_fixed_end_forces[6] = -wy_local*L^2/12
    local_fixed_end_forces[8] = -wy_local*L/2
    local_fixed_end_forces[12] = +wy_local*L^2/12

    k1z, k2z = InstantFrame.connection_zeros_and_inf(k1z, k2z, E, I)

    M_fixed = [local_fixed_end_forces[6], local_fixed_end_forces[12]]
    F_fixed = [local_fixed_end_forces[2], local_fixed_end_forces[8]]
    M_spring, F_spring = fixed_end_forces_partial_restraint(k1z, k2z, E, I, L, M_fixed, F_fixed)

    local_fixed_end_forces[6] = M_spring[1]
    local_fixed_end_forces[12] = M_spring[2]
    local_fixed_end_forces[2] = F_spring[1]
    local_fixed_end_forces[8] = F_spring[2]
    
    local_fixed_end_forces[3] = -wz_local*L/2
    local_fixed_end_forces[5] = +wz_local*L^2/12
    local_fixed_end_forces[9] = -wz_local*L/2
    local_fixed_end_forces[11] = -wz_local*L^2/12

    k1y, k2y = InstantFrame.connection_zeros_and_inf(k1y, k2y, E, I)

    M_fixed = [local_fixed_end_forces[5], local_fixed_end_forces[11]]
    F_fixed = [local_fixed_end_forces[3], local_fixed_end_forces[9]]
    M_spring, F_spring = fixed_end_forces_partial_restraint(k1y, k2y, E, I, L, M_fixed, F_fixed)

    local_fixed_end_forces[5] = M_spring[1]
    local_fixed_end_forces[11] = M_spring[2]
    local_fixed_end_forces[3] = F_spring[1]
    local_fixed_end_forces[9] = F_spring[2]