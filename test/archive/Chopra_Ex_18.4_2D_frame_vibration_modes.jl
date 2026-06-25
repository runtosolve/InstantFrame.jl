using InstantFrame, GLMakie 


#Chopra Example 18.4 - 2D portal frame vibration modes 

material = InstantFrame.Material(names=["steel"], E=[29E6], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam", "column"], A=[10.0, 5.0], Iy=[100.0, 50.0], Iz=[100.0, 50.0], J=[0.001, 0.0005])

connection = InstantFrame.Connection(names=["rigid", "hinge"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 0.0]))

node = InstantFrame.Node(numbers=[1, 2, 3, 4], coordinates=[(0.0, 0.0, 0.0), (0.0, 120.0, 0.0), (240.0, 120.0, 0.0), (240.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2, 3], nodes=[(1,2), (2,3), (3, 4)], orientation=[0.0, 0.0, 0.0], connections=[("rigid", "rigid"), ("rigid", "rigid"), ("rigid", "rigid")], cross_section=["column", "beam", "column"], material=["steel", "steel", "steel"])

support = InstantFrame.Support(nodes=[1, 2, 3, 4], stiffness=(uX=[Inf,0.0,0.0,Inf], uY=[Inf,0.0,0.0,Inf], uZ=[Inf,Inf, Inf,Inf], rX=[Inf,Inf, Inf, Inf], rY=[Inf,Inf, Inf,Inf], rZ=[Inf,0.0, 0.0, Inf]))

uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(nothing)

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "modal vibration")


#Chopra modal frequencies

E = 29000000.0
I = 50.0
h = 120.0
A = 5.0
m = material.ρ[1] * A
ω1 = 1.9004*sqrt(E*I/(m*h^4))
ω2 = 4.6609*sqrt(E*I/(m*h^4))
ω3 = 14.5293*sqrt(E*I/(m*h^4))

#Chopra mode shape for mode 1
ϕ1 = [0.4636, -0.2741/h, -0.2741/h]

#test
#small differences from including axial stiffness?
isapprox(model.solution.ωn[1], ω1, rtol=0.01)
isapprox(model.solution.ωn[2], ω2, rtol=0.02)
isapprox(model.solution.ωn[3], ω3, rtol=0.04)

isapprox(model.solution.ϕ[1, ][[7, 12, 18]] .* 0.4636, ϕ1, atol = 0.01)