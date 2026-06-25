using InstantFrame 


#Sennett Example 3.1 - fixed-fixed beam with two elements, 2D

material = InstantFrame.Material(names=["steel"], E=[10^7], ν=[0.3], ρ=[492.0 / 32.17 / 12^4])  #ρ = lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[1.0], Iy=[1.0], Iz=[1.0], J=[0.001])

connection = InstantFrame.Connection(names=["rigid", "end release"], stiffness=(ux=[Inf, Inf], uy=[Inf, Inf], uz=[Inf, Inf], rx=[Inf, Inf], ry=[Inf, Inf], rz=[Inf, 0.0]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2], nodes=[(1,2), (2,3)], orientation=[0.0, 0.0], connections=[("rigid", "rigid"), ("rigid", "rigid")], cross_section=["beam", "beam"], material=["steel", "steel"])

support = InstantFrame.Support(nodes=[1, 3], stiffness=(uX=[Inf,Inf], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[Inf,Inf]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1, 2], loads=(qX=[0.0, 0.0], qY=[0.0, 0.0], qZ=[0.0, 0.0], mX=[0.0,0.0], mY=[0.0,0.0], mZ=[0.0,0.0]))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[-1.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")

scale = (1.0, 1.0, 1.0)
figure = InstantFrame.UI.display_model_deformed_shape(model.solution.u1, element, node, model.properties, scale)


##tests

P = point_load.loads.FY[1]

#deflections
model.solution.u1[8]≈P/24E7  #midspan deflection
model.solution.u1[12]≈0.0  #midspan rotation

#internal forces
dof = [2, 6, 8, 12]
model.solution.P1[1][dof] == [-P/2,-P/4, P/2, -P/4]
model.solution.P1[2][dof] == [P/2,P/4, -P/2, P/4]




