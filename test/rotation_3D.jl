
using InstantFrame, LinearAlgebra

#McGuire et al. (2000) Example 5.5

material = InstantFrame.Material(names=["steel"], E=[200.0], ν=[0.3], ρ=[8000/1000^3])  #ρ = kg/mm^3

cross_section = InstantFrame.CrossSection(names=["beam ab"], A=[6E3], Iy=[700E6], Iz=[200E6], J=[300E3])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))

node = InstantFrame.Node(numbers=[1, 2], coordinates=[(0.0, 0.0, 0.0), (10.0, 5.0, 3.0)])

element = InstantFrame.Element(numbers=[1], nodes=[(1,2)], orientation=[π/6], connections=[("rigid", "rigid")], cross_section=["beam ab"], material=["steel"])

support = InstantFrame.Support(nodes=[1, 2], stiffness=(uX=[Inf,0.0], uY=[Inf,Inf], uZ=[Inf,Inf], rX=[Inf,Inf], rY=[Inf,Inf], rZ=[Inf,0.0]))

uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[1], loads=(qX=[0.0], qY=[0.0], qZ=[0.0], mX=[0.0], mY=[0.0], mZ=[0.0]))

point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[5cos(-π/4)], FY=[5sin(-π/4)], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")



##test

#compare MGZ rotation matrix from Eq. 5.7 to model.

γ_MGZ = [0.8638     0.4319  0.2591
         -0.5019    0.7811  0.3714
         -0.0420    -0.4510 0.8915]

Γ_MGZ = zeros(Float64, 12, 12)

Γ_MGZ[1:3, 1:3] = γ_MGZ 
Γ_MGZ[4:6, 4:6] = γ_MGZ 
Γ_MGZ[7:9, 7:9] = γ_MGZ 
Γ_MGZ[10:12, 10:12] = γ_MGZ 

isapprox(Γ_MGZ, model.properties.Γ[1], rtol=0.01)




A = [0.0, 0.0, 0.0]
B = [10.0, 5.0, 3.0]
β = 0.0


# AB = B - A

# ΔX = AB[1]
# ΔZ = AB[3]
# ΔY = AB[2]

# ρ = atan(-ΔZ, ΔX)

# rad2deg(ρ)

# proj_AB_xz = sqrt(ΔX^2 + ΔZ^2)

# χ = atan(ΔY, proj_AB_xz)
# rad2deg(χ)

# current_local_y_axis = RotZ(-χ) * RotY(-ρ) * [0.0, 1.0, 0.0]  #where y-local is pointing after Y and Z rotations 

# ω = acos(dot(current_local_y_axis, [0.0, 1.0, 0.0])/ norm(current_local_y_axis))

# rad2deg(ω)



AB = B - A

ΔX = AB[1]
ΔZ = AB[3]
ΔY = AB[2]

ρ = atan(-ΔZ, ΔX)

proj_AB_xz = sqrt(ΔX^2 + ΔZ^2)

χ = atan(ΔY, proj_AB_xz)

RotY(ρ)'

RotZ(χ)'

RotX(π/2 - π/3)'

RotYZX(ρ, χ, π/2 - π/3)'

ω = β

x_axis = RotYZX(ρ, χ, ω) * [1.0, 0.0, 0.0]
y_axis = RotYZX(ρ, χ, ω) * [0.0, 1.0, 0.0]
z_axis = RotYZX(ρ, χ, ω) * [0.0, 0.0, 1.0]

RotYZX(ρ, χ, ω)' * [1.0, 0.0, 0.0]