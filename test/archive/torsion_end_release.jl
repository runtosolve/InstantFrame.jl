using InstantFrame 


#Apply a point torsion at midspan of a beam that is twist fixed.  The beam is discretized with two elements.  The first element has a torsion release at its right end (midspan).

material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[1.0], Iy=[1.0], Iz=[1.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "end release"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, 0.0], ry=[Inf, Inf], rz=[Inf, Inf]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (12.0, 0.0, 0.0), (24.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2], nodes=[(1,2), (2,3)], orientation=[0.0, 0.0, 0.0], connections=[("rigid", "end release"), ("rigid", "rigid")], cross_section=["beam", "beam"], material=["steel", "steel"])

support = InstantFrame.Support(nodes=[1, 3], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[Inf,Inf]))

# uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1], loads=(qX=[0.0], qY=[10.0], qZ=[0.0], mX=[0.0], mY=[0.0], mZ=[0.0]))

uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[1000.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

##test

#Twist at the right end of the first element should be zero.
isapprox(model.solution.element_connections.end_displacements[1][4], 0.0, atol=0.0001)

#Twist at the left end of the second element should be Θ = TL/GJ where L=12.0 in.
isapprox(model.solution.nodal_displacements[2][4], (1000.0*12.0)/(model.properties.G[1]*model.properties.J[1]), rtol=0.01)

