using InstantFrame 


#Sennett Example 7.4 - portal frame with corner hinge, 2D, 1st order

material = InstantFrame.Material(names=["steel"], E=[29E6], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[10.0], Iy=[100.0], Iz=[100.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "hinge"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 0.0]))

node = InstantFrame.Node(numbers=[1, 2, 3, 4], coordinates=[(0.0, 0.0, 0.0), (0.0, 120.0, 0.0), (120.0, 120.0, 0.0), (120.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2, 3], nodes=[(1,2), (2,3), (3, 4)], orientation=[0.0, 0.0, 0.0], connections=[("rigid", "rigid"), ("hinge", "rigid"), ("rigid", "rigid")], cross_section=["beam", "beam", "beam"], material=["steel", "steel", "steel"])

support = InstantFrame.Support(nodes=[1, 2, 3, 4], stiffness=(uX=[Inf,0.0,0.0,Inf], uY=[Inf,0.0,0.0,Inf], uZ=[Inf,Inf, Inf,Inf], rX=[Inf,Inf, Inf, Inf], rY=[Inf,Inf, Inf,Inf], rZ=[0.0,0.0, 0.0, 0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1, 2], loads=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[1000.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

scale = (30.0, 30.0, 30.0)
# figure = InstantFrame.UI.display_model_deformed_shape(model.solution.nodal_displacements, element, node, model.properties, scale)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.nodal_displacements, model.solution.element_connections, element, node, model.properties, scale)

##tests
#rotation at hinge
isapprox(model.solution.element_connections.end_displacements[1][6], 0.0008, rtol=0.05)  #hinge rotation using Sennett Eq. 7.33



