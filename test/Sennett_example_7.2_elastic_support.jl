using InstantFrame 


#Sennett Example 7.2
#Elastic support at the center of a fixed-fixed beam


material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[8.0], Iy=[100.0], Iz=[100.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

beam_span = 240.0
x = range(0.0, beam_span, 3) 
coordinates = [(x[i], 0.0, 0.0) for i in eachindex(x)]

node = InstantFrame.Node(numbers=1:length(x), coordinates=coordinates)

num_elem = length(x)-1
element_connectivity = [(i, i+1) for i=1:num_elem]
element = InstantFrame.Element(numbers=1:length(element_connectivity), nodes=element_connectivity, orientation=zeros(Float64, length(element_connectivity)), connections=[("rigid", "rigid") for i in eachindex(element_connectivity)], cross_section=["beam" for i in eachindex(element_connectivity)], material=["steel" for i in eachindex(element_connectivity)])

support = InstantFrame.Support(nodes=[1, length(x)], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[Inf,Inf]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[2], loads=(qX=zeros(Float64, 1), qY=zeros(Float64, 1), qZ=ones(Float64, 1)*-1000.0, mX=zeros(Float64, 1), mY=zeros(Float64, 1), mZ=zeros(Float64, 1)))

# point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))
point_load = InstantFrame.PointLoad(nothing)

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")



######

element_nodal_coords = InstantFrame.UI.define_element_nodal_start_end_coordinates(element, node)

X, Y, Z = InstantFrame.UI.get_node_XYZ(node)

X_range = abs(maximum(X) - minimum(X))
Y_range = abs(maximum(Y) - minimum(Y))
Z_range = abs(maximum(Z) - minimum(Z))

Y_range = 100.0
Z_range = 100.0

using GLMakie, LinearAlgebra
figure = Figure()
ax = Axis3(figure[1,1])
ax.aspect = (1.0, Y_range/X_range, Z_range/X_range)
ax.yticks = WilkinsonTicks(2)
ylims!(ax, -50.0, 50.0)
zlims!(ax, -50.0, 50.0)

color = :gray
InstantFrame.UI.show_elements!(ax, element_nodal_coords, color)
figure

markersize = 10
color = :blue
InstantFrame.UI.show_nodes!(ax, X, Y, Z, markersize, color)
figure

unit_arrow_head_size = [1.0, 1.0, 1.0]
arrow_head_scale = 3.0
arrow_scale = 10.0
linewidth = 1.0
arrowcolor = :green 
linecolor = :green
InstantFrame.UI.show_uniform_loads!(ax, uniform_load, element, node, unit_arrow_head_size, arrow_head_scale, arrow_scale, linewidth, arrowcolor, linecolor)



figure 

