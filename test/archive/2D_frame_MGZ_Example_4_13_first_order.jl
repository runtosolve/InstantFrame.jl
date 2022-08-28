using InstantFrame, NonlinearSolve, Plots, LinearAlgebra, BenchmarkTools, StaticArrays, LinearAlgebra, BenchmarkTools

#McGuire, W., Gallagher, R. H., & Ziemian, R. D. (2000). Matrix Structural Analysis, John Wiley and Sons. Inc., New York.

#Consider the study summarized in Example 4.13 of McGuire et al. (2000), a 2D rigid frame.


#Define structural system.
nodes = [[0.0, 0.0], [8000.0, 0.0], [8000.0, -5000.0]]

cross_sections = InstantFrame.CrossSections(["member_ab", "member_bc"], [6E3, 4E3], [200E6, 50E6], [], [])

materials = InstantFrame.Materials(["steel"], [200000.], [0.3])

connections = InstantFrame.Connections(["semi-rigid", "rigid"], [30000.0, Inf], [Inf, Inf], [Inf, Inf])

members = InstantFrame.Members(["frame", "frame"], [1, 2], [2, 3], [0.0, 0.0], ["rigid", "rigid"], ["rigid", "rigid"], ["member_ab", "member_bc"], ["steel", "steel"])

#Calculate global elastic stiffness matrix.
Ke = InstantFrame.define_global_elastic_stiffness_matrix(nodes, cross_sections, materials, members)
 
#Define the frame loading as an external force vector. 
F = [0., 0., 0., 100000.0/sqrt(2), -100000.0/sqrt(2), 50000000., 0., 0., 0.]

#Define the frame boundary conditions.
free_dof = [4; 5; 6]

#Partition external force vector and elastic stiffness matrix.
Ff = F[free_dof]
Ke_ff = Ke[free_dof, free_dof]

#Solution
u = Ke_ff \ Ff
