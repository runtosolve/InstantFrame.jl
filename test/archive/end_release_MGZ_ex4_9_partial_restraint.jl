using InstantFrame, LinearAlgebra

#McGuire et. al(2000) Example 4.9 with partial restraint in member ab, b side.


material = InstantFrame.Material(names=["steel"], E=[200.0], ν=[0.3], ρ=[8000/1000^3])  #ρ = kg/mm^3

cross_section = InstantFrame.CrossSection(names=["beam ab", "beam bc"], A=[6E3, 4E3], Iy=[700E6, 540E6], Iz=[200E6, 50E6], J=[300E3, 100E3])

connection = InstantFrame.Connection(names=["rigid", "partial restraint"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 1.0]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (8000.0, 0.0, 0.0), (13000.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2], nodes=[(1,2), (2,3)], orientation=[0.0, 0.0], connections=[("rigid", "partial restraint"), ("rigid", "rigid")], cross_section=["beam ab", "beam bc"], material=["steel", "steel"])

support = InstantFrame.Support(nodes=[1, 2, 3], stiffness=(uX=[Inf,0.0, 0.0], uY=[Inf,Inf,0.0], uZ=[Inf,Inf,Inf], rX=[Inf,Inf,Inf], rY=[Inf,Inf,Inf], rZ=[Inf,0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1, 2], loads=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[3], loads=(FX=[5cos(-π/4)], FY=[5sin(-π/4)], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

scale = (30.0, 30.0, 30.0)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.nodal_displacements, model.solution.element_connections, element, node, model.properties, scale)


##tests

#element local stiffness matrix with partial restrained connection

#Use McGuire et al. (2000) Eq. 13.31
E = model.properties.E[1]
I = model.properties.Iz[1]
L = model.properties.L[1]

k1 = E*I*10E20
k2 = connection.stiffness.rz[2]

α1 = k1/(E*I/L)  #start node

α2 = k2*(E*I/L)  #end node

α = (α1*α2)/(α1*α2 + 4*α1 + 4*α2 + 12)

ke = zeros(Float64, 4, 4)

ke[1,1] = 12/L^2*(1 + (α1+α2)/(α1*α2))
ke[1,2] = 6/L * (1+2/α2)
ke[1,3] = -12/L^2*(1+(α1+α2)/(α1*α2))
ke[1,4] = 6/L*(1+2/α1)
ke[2,2] = 4*(1+3/α2)
ke[2,3] = -6/L*(1+2/α2)
ke[2,4] = 2
ke[3,3] = 12/L^2*(1+(α1+α2)/(α1*α2))
ke[3,4] = -6/L*(1+2/α1)
ke[4,4] = 4*(1+3/α1)

ke = α * (E * I / L) .* ke

#compare to local stiffness matrix from model
dof = [2, 6, 8, 12]
isapprox(triu(ke), triu(model.equations.ke_local[1])[dof, dof])



