using InstantFrame

#use modified version of Sennett Eq. 7-54 to obtain element stiffness matrix including partial restraint at the left end.

I = 100.0
A = 10.0
E = 29E6
L = 120.0

ke = InstantFrame.define_local_elastic_element_stiffness_matrix(I, A, E, L)




#Insert a rotational spring at the left end of the frame element.

# p = [1; 2; 4; 5; 6]
# s = [3]  #spring at left end

kϕ1 = 300000.0 #lb-in./rad  this is consistent with a typical rack connection.

k_spring_left = zeros(Float64, 6, 6)
k_spring_left[3, 3] = kϕ1
k_spring_left[3, 6] = -kϕ1
k_spring_left[6, 3] = -kϕ1
k_spring_left[6, 6] = kϕ1

k_sys = zeros(Float64, 9 ,9)

k_sys[4:9, 4:9] = ke

k_sys[1:6, 1:6] = k_sys[1:6, 1:6] + k_spring_left

s = [4;5;6]
p = [1;2;3;7;8;9]


# k_spring = zeros(Float64, 6, 6)
# k_spring[s, s] .= kϕ_left

k_pp = k_sys[p, p]
k_ps = k_sys[p, s]
k_ss = k_sys[s, s]
k_sp = k_sys[s, p]

# k_spring_s = k_spring[s, s]

ke_condensed = k_pp - k_ps*k_ss^-1*k_sp

ke_updated = zeros(Float64, 6, 6)

ke_updated[p, p] = ke_condensed   


#From Zeimian 13.31...
k1 = kϕ_left
α1 = k1/(E*I/L)

α2 = k1*1000000000000.0

α = (α1 * α2) / (α1*α2 + 4*α1 + 4*α2 + 12)

α*E*I/L * 12/L^2*(1 + (α1 + α2)/(α1*α2))

#This matches the Sennett static condensation, 5081=5081.

#Now analyze a beam with a hinge and vertical roller at the left end, and fixed at the right end.   Apply a vertical force at the left end.

Ke = ke_updated

#Define the frame loading as an external force vector. 
F = [0., 1000.0, 0., 0., 0., 0.]

#Define the frame boundary conditions.
free_dof = [2]

#Partition external force vector and elastic stiffness matrix.
Ff = F[free_dof]
Ke_ff = Ke[free_dof, free_dof]

#Solution
uf = Ke_ff \ Ff

u = zeros(Float64, 12)

u[free_dof] = uf

#This matches Mastan2, 0.1968 in.


