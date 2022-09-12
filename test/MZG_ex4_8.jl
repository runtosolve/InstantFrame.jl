using InstantFrame 

#McGuire et. al(2000) Example 4.8


material = InstantFrame.Material(names=["steel"], E=[200.0], ν=[0.3], ρ=[8000/1000^3])  #ρ = kg/mm^3

cross_section = InstantFrame.CrossSection(names=["beam ab", "beam bc"], A=[6E3, 4E3], Iy=[700E6, 540E6], Iz=[200E6, 50E6], J=[300E3, 100E3])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (8000.0, 0.0, 0.0), (13000.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2], nodes=[(1,2), (2,3)], orientation=[0.0, 0.0], connections=[("rigid", "rigid"), ("rigid", "rigid")], cross_section=["beam ab", "beam bc"], material=["steel", "steel"])

support = InstantFrame.Support(nodes=[1, 2, 3], stiffness=(uX=[Inf,0.0, 0.0], uY=[Inf,Inf,0.0], uZ=[Inf,Inf,0.0], rX=[Inf,Inf,Inf], rY=[Inf,Inf,0.0], rZ=[Inf,0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1, 2], loads=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

# uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[3], loads=(FX=[5cos(-π/4)], FY=[5sin(-π/4)], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")


# dof=[1, 2, 4, 6, 7, 8, 10, 12]
# model.equations.ke_local[2][dof, dof] ./ 200

##tests

#deflections and rotations
isapprox(model.solution.u1[7], 0.024, rtol=0.05)
isapprox(model.solution.u1[13], 0.046, rtol=0.05)
isapprox(model.solution.u1[12], -0.00088, rtol=0.05)
isapprox(model.solution.u1[18], -0.00530, rtol=0.05)
isapprox(model.solution.u1[14], -19.15, rtol=0.05)

scale = (1.0, 1.0, 1.0)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.u1, element, node, model.properties, scale)

