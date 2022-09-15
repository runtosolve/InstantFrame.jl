using InstantFrame 

#Denavit and Hajjar (2013) Figure 5
#http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf

#Vibration modes 

material = InstantFrame.Material(names=["steel"], E=[29000.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[9.12], Iy=[37.1], Iz=[110.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

node = InstantFrame.Node(numbers=[1, 2], coordinates=[(0.0, 0.0, 0.0), (180.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1], nodes=[(1,2)], orientation=[0.0], connections=[("rigid", "rigid")], cross_section=["beam"], material=["steel"])

support = InstantFrame.Support(nodes=[1, 2], stiffness=(uX=[Inf,0.0], uY=[Inf,0.0], uZ=[Inf,0.0], rX=[Inf,0.0], rY=[Inf,0.0], rZ=[Inf,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1], loads=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "modal vibration")


##tests

#first mode natural frequency for a cantilever
#https://roymech.org/Useful_Tables/Vibrations/Natural_Vibrations.html

A = cross_section.A[1]
E = material.E[1]
I_weak = cross_section.Iy[1]
I_strong = cross_section.Iz[1]
L = 180.0
ρ = material.ρ[1]
m = ρ * A  #mass per unit length 

ωn_weak = 1.875^2 * sqrt((E*I_weak)/(m*L^4))
ωn_strong = 1.875^2 * sqrt((E*I_strong)/(m*L^4))


isapprox(model.solution.ωn[1], ωn_weak, rtol=0.05)
isapprox(model.solution.ωn[2], ωn_strong, rtol=0.05)




