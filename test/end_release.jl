using InstantFrame, LinearAlgebra

#McGuire et. al(2000) Example 4.9 with end release at point 'a'.


material = InstantFrame.Material(names=["steel"], E=[200.0], ν=[0.3], ρ=[8000/1000^3])  #ρ = kg/mm^3

cross_section = InstantFrame.CrossSection(names=["beam ab", "beam bc"], A=[6E3, 4E3], Iy=[700E6, 540E6], Iz=[200E6, 50E6], J=[300E3, 100E3])

connection = InstantFrame.Connection(names=["rigid", "hinge"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 0.0]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (8000.0, 0.0, 0.0), (13000.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2], nodes=[(1,2), (2,3)], orientation=[0.0, 0.0], connections=[("hinge", "rigid"), ("rigid", "rigid")], cross_section=["beam ab", "beam bc"], material=["steel", "steel"])

support = InstantFrame.Support(nodes=[1, 2, 3], stiffness=(uX=[Inf,0.0, 0.0], uY=[Inf,Inf,0.0], uZ=[Inf,Inf,Inf], rX=[Inf,Inf,Inf], rY=[Inf,Inf,Inf], rZ=[Inf,0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1, 2], loads=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[3], loads=(FX=[5cos(-π/4)], FY=[5sin(-π/4)], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

scale = (30.0, 30.0, 30.0)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.nodal_displacements, model.solution.element_connections, element, node, model.properties, scale)


##tests

#element local stiffness matrix with hinge

#Use Sennett Eq. 7.39
E = model.properties.E[1]
I = model.properties.Iz[1]
L = model.properties.L[1]

ke = zeros(Float64, 4, 4)

#define upper triangular portion of stiffness matrix with Eq. 7.39
ke[1, 1] = 3*E*I/L^3
ke[1, 2] = 0.0
ke[1, 3] = -3*E*I/L^3
ke[1, 4] = 3*E*I/L^2
ke[2, 2] = 0.0
ke[2, 3] = 0.0
ke[2, 4] = 0.0
ke[3, 3] = 3*E*I/L^3
ke[3, 4] = -3*E*I/L^2
ke[4, 4] = 3*E*I/L


#compare to local stiffness matrix from model
dof = [2, 6, 8, 12]
isapprox(triu(ke), triu(model.equations.ke_local[1])[dof, dof])


