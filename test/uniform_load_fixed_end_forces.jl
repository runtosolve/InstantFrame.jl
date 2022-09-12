using InstantFrame 


material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[1.0], Iy=[1.0], Iz=[1.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "end release"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 0.0]))

node = InstantFrame.Node(numbers=[1, 2], coordinates=[(0.0, 0.0, 0.0), (120.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1], nodes=[(1,2)], orientation=[0.0, 0.0], connections=[("rigid", "rigid")], cross_section=["beam"], material=["steel"])

support = InstantFrame.Support(nodes=[1, 2], stiffness=(uX=[Inf,0.0], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,0.0], rY=[0.0,0.0], rZ=[0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1], loads=(qX=[0.0], qY=[10.0], qZ=[0.0], mX=[0.0], mY=[0.0], mZ=[0.0]))

# uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

scale = (1.0, 1.0, 1.0)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.u1, element, node, model.properties, scale)

wY = uniform_load.loads.qY[1]
L = sum(model.properties.L)
E = material.E[1]
Iz = cross_section.Iz[1]
θ_start = wY*L^3/(24*E*Iz)

##test

#rotation at left end of beam
θ_start ≈ model.solution.u1[6]

#add nodal displacements to output