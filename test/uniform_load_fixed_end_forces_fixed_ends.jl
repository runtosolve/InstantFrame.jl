using InstantFrame 


#Apply a uniform load to a simply-supported frame element.

material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[1.0], Iy=[1.0], Iz=[1.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "end release"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 0.0]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (120.0, 0.0, 0.0), (240.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2], nodes=[(1,2), (2,3)], orientation=[0.0, 0.0], connections=[("rigid", "rigid"), ("rigid", "rigid")], cross_section=["beam", "beam"], material=["steel", "steel"])

support = InstantFrame.Support(nodes=[1, 3], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[Inf,Inf]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1, 2], loads=(qX=[0.0, 0.0], qY=[-10.0,-10.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

# uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

# scale = (1.0, 1.0, 1.0)
# figure = InstantFrame.UI.display_model_deformed_shape(model.solution.u1, element, node, model.properties, scale)

w=-10.0
L = 240.0
E = 29000000.0
I = 1.0

Δ = w*L^4/(384*E*I)

##test

#Calculate rotation at the left end of the beam.

wY = uniform_load.loads.qY[1]
L = sum(model.properties.L)
E = material.E[1]
Iz = cross_section.Iz[1]
θ_start = wY*L^3/(24*E*Iz)


#rotation at left end of beam
θ_start ≈ model.solution.nodal_displacements[1][6]

