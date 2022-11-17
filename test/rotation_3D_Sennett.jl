using LinearAlgebra, Rotations 


#Sennett Example 6.4


function rotation_matrix(A, B, β)

    AB = B - A

    # A_prime = [A[1], A[2]+1.0, A[3]]

    # AA_prime = A_prime - A

    # AK = RotX(β) * AA_prime

    # AK = K - A

    x_cosines = AB/norm(AB)

    # AK = A + x_cosines * norm(AB)/2 + [0.0, 1.0, 0.0]

    AK = x_cosines * norm(AB)/2 + RotX(-β) * [0.0, 1.0, 0.0]

    z_cosines = cross(AB, AK) / norm(cross(AB, AK))

    y_cosines = cross(z_cosines, x_cosines)/norm(cross(z_cosines, x_cosines))

    ℓ = [x_cosines'
        y_cosines'
        z_cosines']

    return ℓ

end







#Member 1

A = [180.0, 0.0, 180.0]
B = [0.0, 0.0, 180.0]
β = π/2


AB = B - A

ΔX = AB[1]
ΔZ = AB[3]
ΔY = AB[2]

ρ = atan(-ΔZ, ΔX)

proj_AB_xz = sqrt(ΔX^2 + ΔZ^2)

χ = atan(ΔY, proj_AB_xz)

RotYZ(ρ, χ)

x_axis = RotYZ(ρ, χ) * [1.0, 0.0, 0.0]
y_axis = RotYZ(ρ, χ) * [0.0, 1.0, 0.0]
z_axis = RotYZ(ρ, χ) * [0.0, 0.0, 1.0]

ω = π

RotYZ(-ρ, -χ)

# RotZY(ρ, χ) * [1.0, 0.0, 0.0]
# RotZY(ρ, χ) * [0.0, 1.0, 0.0]

# current_local_y_axis = RotZ(-χ) * RotY(-ρ) * [0.0, 1.0, 0.0]  #where y-local is pointing after Y and Z rotations 

# ω = acos(dot(current_local_y_axis, [0.0, 1.0, 0.0])/ norm(current_local_y_axis))

RotX(-ω) * RotZ(-χ) * RotY(-ρ) * [0.0, 0.0, 1.0]



#Member 2


A = [180.0, 0.0, 180.0]
B = [180.0, 180.0, 180.0]
β = -π/2


AB = B - A

ΔX = AB[1]
ΔZ = AB[3]
ΔY = AB[2]

ρ = atan(-ΔZ, ΔX)

proj_AB_xz = sqrt(ΔX^2 + ΔZ^2)

χ = atan(ΔY, proj_AB_xz)


RotYZX(ρ, χ, β)

ω = β

x_axis = RotYZX(ρ, χ, ω) * [1.0, 0.0, 0.0]
y_axis = RotYZX(ρ, χ, ω) * [0.0, 1.0, 0.0]
z_axis = RotYZX(ρ, χ, ω) * [0.0, 0.0, 1.0]

RotYZX(ρ, χ, ω)'

# x_axis = RotYZ(ρ, χ) * [1.0, 0.0, 0.0]
# y_axis = RotYZ(ρ, χ) * [0.0, 1.0, 0.0]
# z_axis = RotYZ(ρ, χ) * [0.0, 0.0, 1.0]


# current_local_y_axis = RotZ(-χ) * RotY(-ρ) * [0.0, 1.0, 0.0]  #where y-local is pointing after Y and Z rotations 

# ω = acos(dot(current_local_y_axis, [0.0, 1.0, 0.0])/ norm(current_local_y_axis))


# proj_AB_xy = sqrt(ΔX^2 + ΔY^2)

# ω = atan(ΔZ, proj_AB_xy)

# ω = π/2

# ρ = 0.0
# χ = π/2
# ω = π/2

# RotX(-ω) * RotZ(-χ) * RotY(-ρ) * [0.0, 0.0, 1.0]


#Member 3

A = [180.0, 0.0, 180.0]
B = [180.0, 0.0, 0.0]
β = -π/2


AB = B - A

ΔX = AB[1]
ΔZ = AB[3]
ΔY = AB[2]

ρ = atan(-ΔZ, ΔX)

proj_AB_xz = sqrt(ΔX^2 + ΔZ^2)

χ = atan(ΔY, proj_AB_xz)


RotYZX(ρ, χ, β)

ω = β

x_axis = RotYZX(ρ, χ, ω) * [1.0, 0.0, 0.0]
y_axis = RotYZX(ρ, χ, ω) * [0.0, 1.0, 0.0]
z_axis = RotYZX(ρ, χ, ω) * [0.0, 0.0, 1.0]

RotYZX(ρ, χ, ω)' * [1.0, 0.0, 0.0]
RotYZX(ρ, χ, ω)' * [0.0, 1.0, 0.0]
RotYZX(ρ, χ, ω)' * [0.0, 0.0, 1.0]

RotYZX(ρ, χ, ω)' * [1.0, 0.0, 0.0]





# current_local_y_axis = RotZ(-χ) * RotY(-ρ) * [0.0, 1.0, 0.0]  #where y-local is pointing after Y and Z rotations 

# ω = acos(dot(current_local_y_axis, [0.0, 1.0, 0.0])/ norm(current_local_y_axis))


# ω = 0.0

# # ρ = π/2
# # χ = 0.0
# # ω = 0.0

# RotX(-ω) * RotZ(-χ) * RotY(-ρ) * [0.0, 0.0, 1.0]
