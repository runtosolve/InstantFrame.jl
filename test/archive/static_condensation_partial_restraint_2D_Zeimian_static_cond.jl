using InstantFrame

#use modified version of Sennett Eq. 7-54 to obtain element stiffness matrix including partial restraint at the left end.

I = 100.0
A = 10.0
E = 29E6
L = 120.0

ke = InstantFrame.define_local_elastic_element_stiffness_matrix(I, A, E, L)


k1 = 300000.0 #lb-in./rad  this is consistent with a typical rack connection.
k2 = 300000.0 #lb-in./rad  this is consistent with a typical rack connection.


ke = InstantFrame.define_local_elastic_element_stiffness_matrix_condensed(I, A, E, L, k1, k2)


#From Zeimian 13.31...

α1 = k1/(E*I/L)
α2 = k2/(E*I/L)

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


