using InstantFrame 


#Apply a uniform load to a fixed-fixed beam.
#Discretize it with 2 elements.

material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[1.0], Iy=[1.0], Iz=[1.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "end release"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 0.0]))

beam_span = 120.0
x = range(0.0, beam_span, 9) 
coordinates = [(x[i], 0.0, 0.0) for i in eachindex(x)]

node = InstantFrame.Node(numbers=1:length(x), coordinates=coordinates)

num_elem = length(x)-1
element_connectivity = [(i, i+1) for i=1:num_elem]
element = InstantFrame.Element(numbers=1:length(element_connectivity), nodes=element_connectivity, orientation=zeros(Float64, length(element_connectivity)), connections=[("rigid", "rigid") for i in eachindex(element_connectivity)], cross_section=["beam" for i in eachindex(element_connectivity)], material=["steel" for i in eachindex(element_connectivity)])

support = InstantFrame.Support(nodes=[1, length(x)], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[Inf,Inf]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=1:num_elem, loads=(qX=zeros(Float64, num_elem), qY=ones(Float64, num_elem)*10.0, qZ=zeros(Float64, num_elem), mX=zeros(Float64, num_elem), mY=zeros(Float64, num_elem), mZ=zeros(Float64, num_elem)))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

scale = (1000.0, 1000.0, 1000.0)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.nodal_displacements, model.solution.element_connections, element, node, model.properties, scale)



##test

#Calculate deflection at the center of the beam.
wY = uniform_load.loads.qY[1]
L = sum(model.properties.L)
E = material.E[1]
Iz = cross_section.Iz[1]
Δ_midspan = wY*L^4/(384*E*Iz)


#midspan beam deflection
Δ_midspan ≈ model.solution.nodal_displacements[Int(num_elem/2)+1][2]



elem_num = 1
start_node = element.nodes[elem_num][1]
end_node = element.nodes[elem_num][2]

element_global_disp = [model.solution.nodal_displacements[start_node]; model.solution.nodal_displacements[end_node]]

element_local_disp = model.properties.Γ[elem_num] * element_global_disp

#local x-z plane
q1 =  element_local_disp[2]
q2 =  element_local_disp[6]
q3 =  element_local_disp[8]
q4 =  element_local_disp[12]
L = model.properties.L[elem_num]
E = material.E[1]
I = cross_section.Iz[1]
x = range(0.0, L, 5)

function calculate_local_element_moment(q1, q2, q3, q4, L, x, E, I)

    a2 = 1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
    a3 = 1/L^3*(2*q1+q2*L-2*q3+q4*L)

    w_xx = 2*a2 + 6*a3*x
    M = E * I * w_xx

    return M

end

M = calculate_local_element_moment.(q1, q2, q3, q4, L, x, E, I)
