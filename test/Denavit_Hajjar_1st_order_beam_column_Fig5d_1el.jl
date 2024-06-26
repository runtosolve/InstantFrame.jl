using InstantFrame 

#Denavit and Hajjar (2013) Figure 5
#http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf


material = InstantFrame.Material(names=["steel"], E=[29000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[9.12], Iy=[37.1], Iz=[110.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

node = InstantFrame.Node(numbers=[1, 2], coordinates=[(0.0, 0.0, 0.0), (180.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1], nodes=[(1,2)], orientation=[0.0], connections=[("rigid", "rigid")], cross_section=["beam"], material=["steel"])

support = InstantFrame.Support(nodes=[1, 2], stiffness=(uX=[Inf,0.0], uY=[Inf,Inf], uZ=[Inf,0.0], rX=[Inf,Inf], rY=[Inf,0.0], rZ=[Inf,Inf]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1], magnitudes=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], magnitudes=(FX=[-50.0], FY=[0.0], FZ=[1.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

analysis_type = "first order"
model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type)


##tests
#deflection
isapprox(model.solution.nodal_displacements[2][3], 1.807, rtol=0.05)

#cantilever moment 
isapprox(model.solution.element_forces[1][5], 180.0, rtol=0.05) 



