using InstantFrame 

material = InstantFrame.Material(names=["steel"], E=[29500.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000])  ##ρ = kips * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["top chord", "bottom chord", "post"], A=[1.0, 2.0, 3.0], Iy=[3.1, 4.2, 2.3], Iz= [2.4, 6.5, 4.6], J=[0.001, 0.002, 0.005])

connection = InstantFrame.Connection(names=["rigid", "pinned"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, 0.0], rz=[Inf, 0.0]))

node = InstantFrame.Node(numbers=[1, 2, 3, 4], coordinates=[(0.0, 0.0, 0.0), (300.0, 0.0, 0.0), (600.0, 0.0, 0.0), (300.0, -30.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2, 3, 4, 5], nodes = [(1,2), (2,3), (1,4), (3,4), (2,4)], orientation = [0.0, 0.0, 0.0, 0.0, 0.0], connections=[("rigid", "rigid"), ("rigid", "rigid"), ("pinned", "pinned"), ("pinned", "pinned"), ("rigid", "pinned")], cross_section= ["top chord", "top chord", "bottom chord", "bottom chord", "post"], material = ["steel", "steel", "steel", "steel", "steel"], types=["frame", "frame", "frame", "frame", "frame"])

support = InstantFrame.Support(nodes=[1, 3], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["snow", "snow"], elements=[1, 2], magnitudes=(qX=[0.0, 0.0], qY=[-0.100, -0.100], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

point_load = InstantFrame.PointLoad(labels = ["lights"], nodes=[4], magnitudes=(FX=[0.0], FY=[-0.500], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type="first order")