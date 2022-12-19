using InstantFrame 


#Apply a point load to a fixed-fixed beam with partial restraint (rotational springs) at each end.
#Discretize it with 2 elements.

material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[1.0], Iy=[0.627], Iz=[0.55967], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "PR"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 230000.0]))

beam_span = 96.0
x = range(0.0, beam_span, 3) 
coordinates = [(x[i], 0.0, 0.0) for i in eachindex(x)]

node = InstantFrame.Node(numbers=1:length(x), coordinates=coordinates)

num_elem = length(x)-1
element_connectivity = [(i, i+1) for i=1:num_elem]
element = InstantFrame.Element(numbers=1:length(element_connectivity), nodes=element_connectivity, orientation=zeros(Float64, length(element_connectivity)), connections=[("PR", "rigid"), ("rigid", "PR")], cross_section=["beam" for i in eachindex(element_connectivity)], material=["steel" for i in eachindex(element_connectivity)])

support = InstantFrame.Support(nodes=[1, length(x)], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[Inf,Inf]))

# uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=1:num_elem, loads=(qX=zeros(Float64, num_elem), qY=ones(Float64, num_elem)*10.0, qZ=zeros(Float64, num_elem), mX=zeros(Float64, num_elem), mY=zeros(Float64, num_elem), mZ=zeros(Float64, num_elem)))
uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[-1000.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")




element_nodal_coords = InstantFrame.UI.define_element_nodal_start_end_coordinates(element, node)

X, Y, Z = InstantFrame.UI.get_node_XYZ(node)

X_range = abs(maximum(X) - minimum(X))
Y_range = abs(maximum(Y) - minimum(Y))
Z_range = abs(maximum(Z) - minimum(Z))


using GLMakie, LinearAlgebra
figure = Figure()
ax = Axis3(figure[1,1])
ax.aspect = (1.0, 1.0, 1.0)
ax.yticks = WilkinsonTicks(2)
# ylims!(ax, 0.0, 10.0)

color = :gray
InstantFrame.UI.show_elements!(ax, element_nodal_coords, color)
figure

markersize = 1000
color = :blue
InstantFrame.UI.show_nodes!(ax, X, Y, Z, markersize, color)
figure

unit_arrow_head_size = [1.0, 1.0, 1.0]
arrow_head_scale = 0.01
arrow_scale = 5.0
arrowcolor = :orange 
linecolor = :orange 
linewidth = 0.01

InstantFrame.UI.show_element_local_y_axis!(ax, element, node, model, unit_arrow_head_size, arrow_head_scale, arrow_scale, arrowcolor, linecolor, linewidth)
figure






# scale = (1000.0, 1000.0, 1000.0)
# figure = InstantFrame.UI.display_model_deformed_shape(model.solution.nodal_displacements, model.solution.element_connections, element, node, model.properties, scale)



# ##test

# #Calculate deflection at the center of the beam.
# wY = uniform_load.loads.qY[1]
# L = sum(model.properties.L)
# E = material.E[1]
# Iz = cross_section.Iz[1]
# Δ_midspan = wY*L^4/(384*E*Iz)


# #midspan beam deflection
# Δ_midspan ≈ model.solution.nodal_displacements[Int(num_elem/2)+1][2]



# elem_num = 1
# start_node = element.nodes[elem_num][1]
# end_node = element.nodes[elem_num][2]

# element_global_disp = [model.solution.nodal_displacements[start_node]; model.solution.nodal_displacements[end_node]]

# element_local_disp = model.properties.Γ[elem_num] * element_global_disp

# #local x-z plane
# q1 =  element_local_disp[2]
# q2 =  element_local_disp[6]
# q3 =  element_local_disp[8]
# q4 =  element_local_disp[12]
# L = model.properties.L[elem_num]
# E = material.E[1]
# I = cross_section.Iz[1]
# x = range(0.0, L, 5)

# function calculate_local_element_moment(q1, q2, q3, q4, L, x, E, I)

#     a2 = 1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
#     a3 = 1/L^3*(2*q1+q2*L-2*q3+q4*L)

#     w_xx = 2*a2 + 6*a3*x
#     M = E * I * w_xx

#     return M

# end

# M = calculate_local_element_moment.(q1, q2, q3, q4, L, x, E, I)
