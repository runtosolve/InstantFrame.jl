using InstantFrame 

#Denavit and Hajjar (2013) Figure 5
#http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf


material = InstantFrame.Material(names=["steel"], E=[29000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[9.12], Iy=[37.1], Iz=[110.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (90.0, 0.0, 0.0), (180.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2], nodes=[(1,2), (2,3)], orientation=[0.0, 0.0], connections=[("rigid", "rigid"), ("rigid", "rigid")], cross_section=["beam", "beam"], material=["steel", "steel"])

support = InstantFrame.Support(nodes=[1, 2, 3], stiffness=(uX=[Inf,0.0,0.0], uY=[Inf,0.0,0.0], uZ=[Inf,0.0,0.0], rX=[Inf,Inf,Inf], rY=[Inf,0.0,0.0], rZ=[Inf,0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1, 2], loads=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

# uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[3], loads=(FX=[-50.0], FY=[1.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "second order")


##tests

#deflection
isapprox(model.solution.u2[8], 0.765, rtol=0.05)

#cantilever moment 
isapprox(-model.solution.P2[1][6], 218.3, rtol=0.05)

scale = (1.0, 1.0, 1.0)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.u1, element, node, model.properties, scale)

