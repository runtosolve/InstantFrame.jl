using InstantFrame

#use Sennett Eq. 7-54 to obtain element stiffness matrix including a hinge at the left end.

I = 100.0
A = 10.0
E = 29E6
L = 120.0

ke = InstantFrame.define_local_elastic_element_stiffness_matrix(I, A, E, L)

p = [1; 2; 4; 5; 6]
s = [3]  #hinge at left end

k_pp = ke[p, p]
k_ps = ke[p, s]
k_ss = ke[s, s]
k_sp = ke[s, p]

ke_condensed = k_pp - k_ps*k_ss^-1*k_sp

ke_updated = zeros(Float64, 6, 6)

ke_updated[p, p] = ke_condensed   #This matches Sennett book Eq. 7-38.




