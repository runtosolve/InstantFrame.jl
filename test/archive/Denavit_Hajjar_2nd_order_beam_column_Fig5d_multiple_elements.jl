using InstantFrame, StaticArrays, NLsolve, LineSearches

#Denavit and Hajjar (2013) Figure 5
#http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf


material = InstantFrame.Material(names=["steel"], E=[29000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[9.12], Iy=[37.1], Iz=[110.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

tip_disp = Vector{Float64}(undef, 10)
base_moment = Vector{Float64}(undef, 10)

for i =2:11

    num_nodes = i
    x_range = range(0.0, 180.0, num_nodes)
    node = InstantFrame.Node(numbers=range(1, num_nodes, num_nodes), coordinates=[(x_range[i], 0.0, 0.0) for i in eachindex(x_range)])

    num_elem = num_nodes - 1
    element = InstantFrame.Element(numbers=range(1, num_elem, num_elem), nodes=[(i,i+1) for i=1:num_elem], orientation=fill(0.0, num_elem), connections=fill(("rigid", "rigid"), num_elem), cross_section=fill("beam", num_elem), material=fill("steel", num_elem))

    support = InstantFrame.Support(nodes=[1, num_nodes], stiffness=(uX=[Inf,0.0], uY=[Inf,Inf], uZ=[Inf,0.0], rX=[Inf,Inf], rY=[Inf,0.0], rZ=[Inf,Inf]))

    # uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1], magnitudes=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))
    uniform_load = InstantFrame.UniformLoad(nothing)

    point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[num_nodes], magnitudes=(FX=[-50.0], FY=[0.0], FZ=[1.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

    analysis_type = "second order"
    solution_tolerance = 1E-9
    model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type, solution_tolerance)

    tip_disp[i-1] = model.solution.displacements[end][3]
    base_moment[i-1] = model.solution.reactions[1][5]

end

using Plots
plot(1:10, tip_disp, markershape = :o)

plot(1:10, base_moment, markershape = :o)


#global to local

node.coordinates = [node.coordinates[i] .+ tuple(model.solution.displacements[i][1:3]...) for i in eachindex(node.coordinates)]

analysis_type = "first order"
model_update = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type)



u_elem_g = model.solution.u

InstantFrame.define_nodal_displacements(node, model.solution.u)



rotation_angles, Γ = update_rotation_angles(node, element, model.properties, model.solution.u)

u_elem_l = model.properties.Γ[1] * model.solution.u

u_elem_l = Γ[1]' * u_elem_l

P_rotated = model.equations.ke_local[1] * u_elem_l

#convert local to new local 






model.properties.Γ = Γ

P = InstantFrame.calculate_element_internal_forces(model.properties, model.equations.ke_local, element, uniform_load, model.forces.local_fixed_end_forces, model.solution.u)



#########

properties = InstantFrame.define_element_properties(node, cross_section, material, element, connection)

ke_local = [InstantFrame.define_local_elastic_stiffness_matrix(properties.Iy[i], properties.Iz[i], properties.A[i], properties.J[i], properties.E[i], properties.ν[i], properties.L[i]) for i in eachindex(properties.L)]

ke_local = InstantFrame.modify_element_local_connection_stiffness(properties, ke_local, element)

ke_global = [properties.Γ[i]'*ke_local[i]*properties.Γ[i] for i in eachindex(properties.L)]

Ke = InstantFrame.assemble_global_matrix(ke_global, properties.global_dof)

free_global_dof, fixed_global_dof, elastic_supports = InstantFrame.define_free_global_dof(node, support)

for i in eachindex(elastic_supports.global_dof)  #add springs 

    Ke[elastic_supports.global_dof[i], elastic_supports.global_dof[i]] += elastic_supports.global_stiffness[i]

end

equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load = InstantFrame.calculate_nodal_forces_from_uniform_loads(uniform_load, element, node, properties)

global_dof_point_loads = InstantFrame.define_global_dof_point_loads(node, point_load)

F = global_dof_point_loads .+ global_dof_nodal_forces_uniform_load

forces = InstantFrame.Forces(equiv_global_nodal_forces_uniform_load, local_fixed_end_forces, global_dof_nodal_forces_uniform_load, global_dof_point_loads, F)


Ff = F[free_global_dof]
Ke_ff = Ke[free_global_dof, free_global_dof]
u1f = Ke_ff \ Ff

u1 = zeros(Float64, size(Ke,1))
u1[free_global_dof] = u1f



# P1 = InstantFrame.calculate_element_internal_forces(properties, ke_local, u1)

#calculate element internal forces
P1 = InstantFrame.calculate_element_internal_forces(properties, ke_local, element, uniform_load, local_fixed_end_forces, u1)
P1_axial = [P1[i][7] for i in eachindex(P1)]
kg_local = [InstantFrame.define_local_geometric_stiffness_matrix(P1_axial[i], properties.L[i]) for i in eachindex(properties.L)]
kg_global = [properties.Γ[i]'*kg_local[i]*properties.Γ[i] for i in eachindex(properties.L)]
Kg = InstantFrame.assemble_global_matrix(kg_global, properties.global_dof)
Kg_ff = Kg[free_global_dof, free_global_dof]

Kff = Ke_ff + Kg_ff

# u = fill(0.0, 100)
# u[free_global_dof] = convert(Vector{Float64}, u1f)

# length(node.numbers)*6
# u_all = SVector{num_dof}(zeros(Float64, num_dof))

# function residual(uf, p)
#     Ke_ff, Kg_ff, Ff, node, element, properties, free_global_dof, fixed_global_dof = p
#     num_dof = length(node.numbers)*6

#     u = Vector{Any}(undef, num_dof)
#     u[free_global_dof] = uf 
#     u[fixed_global_dof] .= 0.0
#     # u = zeros(Any, num_dof)
#     # u_all = MVector{num_dof}(zeros(Float64, num_dof))
#     # u_all = deepcopy(uf[1])
#     #convert(Vector{Float64}, uf)
#     print("type is " * string(typeof(uf)))
#     print("          ")
#     Ke_ff = update_Ke(node, element, properties, u)
#     Kff = Ke_ff + Kg_ff

#     Kff * uf - Ff
# end

# p = [Ke_ff, Kg_ff, Ff, node, element, properties, free_global_dof, fixed_global_dof]
# u1f = SVector{length(u1f)}(u1f)
# probN = NonlinearSolve.NonlinearProblem{false}(residual, u1f, p)
# solution_tolerance = 0.01
# u2f = NonlinearSolve.solve(probN, NewtonRaphson(), reltol = solution_tolerance)

# u1f = 0.01*u1f

solution = nlsolve((R,uf) ->residual!(R, uf, Ke_ff, Kg_ff, Ff, free_global_dof, ke_local, uniform_load, local_fixed_end_forces, properties, node), u1f, method = :newton)



function residual!(R, uf, Ke_ff, Kg_ff, Ff, free_global_dof, ke_local, uniform_load, local_fixed_end_forces, properties, node)

    # num_dof = length(node.numbers) * 6
    # u = zeros(Float64, num_dof)
    # u[free_global_dof] = uf 
    # Kg_ff = update_Kg(element, properties, ke_local, uniform_load, local_fixed_end_forces, u)

    # print(Ke_ff)
    # print("             ")

    Kff = Ke_ff + Kg_ff 

    for i in eachindex(Ff)
 
       R[i] = transpose(Kff[i,:]) * uf - Ff[i]
 
    end
 
    return R
 
 end



function update_Kg(element, properties, ke_local, uniform_load, local_fixed_end_forces, u)

    Γ = update_rotation_angles(node, element, properties, u)
    properties.Γ =  deepcopy(Γ)
    # properties.Γ =  [0.0]

    P1 = InstantFrame.calculate_element_internal_forces(properties, ke_local, element, uniform_load, local_fixed_end_forces, u)
    P1_axial = [P1[i][7] for i in eachindex(P1)]
    # P1_axial = [-50.0]
    kg_local = [InstantFrame.define_local_geometric_stiffness_matrix(P1_axial[i], properties.L[i]) for i in eachindex(properties.L)]
    kg_global = [Γ[i]'*kg_local[i]*Γ[i] for i in eachindex(properties.L)]
    Kg = InstantFrame.assemble_global_matrix(kg_global, properties.global_dof)
    Kg_ff = Kg[free_global_dof, free_global_dof]

    print(P1_axial)
    print("        ")

    return Kg_ff

end



function update_rotation_angles(node, element, properties, u)

    num_elem = length(properties.L)

    Γ = Array{Array{Float64, 2}}(undef, num_elem)
    rotation_angles = Array{NamedTuple{(:Y, :Z, :X), NTuple{3, Float64}}}(undef, num_elem)

    nodal_displacements = InstantFrame.define_nodal_displacements(node, u)

    node_i = collect(node.coordinates[1])

    for i in eachindex(element.numbers)

        #nodal coordinates
        node_i_index = findfirst(node_num->node_num == element.nodes[i][1], node.numbers)
        node_j_index = findfirst(node_num->node_num == element.nodes[i][2], node.numbers)
        node_i = collect(node.coordinates[node_i_index]) + nodal_displacements[node_i_index][1:3]  #update geometry
        node_j = collect(node.coordinates[node_j_index]) + nodal_displacements[node_j_index][1:3]  #update geometry

        #rotation matrix
        Γ[i], rotation_angles[i] = InstantFrame.define_rotation_matrix(node_i, node_j, element.orientation[i])

        # print(Γ[i])
        # print("                ")

    end

    return rotation_angles, Γ

end



##############################

##tests
#deflection
isapprox(model.solution.nodal_displacements[2][3], 4.597, rtol=0.05)

#cantilever moment 
isapprox(model.solution.element_forces[1][5], 409.8, rtol=0.08)  #tolerance is a little loose here



α = 1/50
L = 180.0
P = 50.0
E = 29000.0
I = 37.1

Δ = α*L * (tan(sqrt(P/(E*I))*L)/(sqrt(P/(E*I))*L) - 1)
