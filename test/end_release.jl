using InstantFrame 


material = InstantFrame.Material(names=["steel"], E=[29000000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[1.0], Iy=[1.0], Iz=[1.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

node = InstantFrame.Node(numbers=[1, 2], coordinates=[(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1], nodes=[(1,2)], orientation=[0.0], connections=[("rigid", "rigid")], cross_section=["beam"], material=["steel"])

support = InstantFrame.Support(nodes=[1, 2], stiffness=(uX=[Inf,0.0], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,0.0], rY=[0.0,0.0], rZ=[0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(elements=[1], loads=(qX=[0.0], qY=[1.0], qZ=[0.0], mX=[0.0], mY=[0.0], mZ=[0.0]))

point_load = InstantFrame.PointLoad(nodes=[1], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))


model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "modal vibration")



