using Rotations
#Consider Zeimian and McGuire Example 5.5

#global 

A = [0, 0, 0]
B = [10, 5, 3]

AB = B - A

ρ = deg2rad(-16.699)
χ = deg2rad(25.589)
ω = deg2rad(30.0)

RotY(-ρ)

RotZ(-χ)

RotX(-ω)



RotXZY(-ω, -χ, -ρ)

RotX(0.1) * RotY(0.2) * RotZ(0.3) === RotXYZ(0.1, 0.2, 0.3)





#local 

# A_prime = [0, 0, 0]
# B_prime = [length_AB, 0.0, 0.0]


using LinearAlgebra, Rotations
#[x'] = [β][x]

β = deg2rad(30.0)

A = [0, 0, 0]
B = [10, 5, 3]


AB = B - A

length_AB = norm(AB)

χ = π/2 - acos(AB[2]/length_AB)


ρ = -atan(AB[3]/AB[1])




ω = β

rad2deg(χ)
rad2deg(ρ)
rad2deg(ω)

RotXZY(-ω, -χ, -ρ)
