using InstantFrame, NonlinearSolve, Plots, LinearAlgebra, BenchmarkTools, StaticArrays, LinearAlgebra, BenchmarkTools

#McGuire, W., Gallagher, R. H., & Ziemian, R. D. (2000). Matrix Structural Analysis, John Wiley and Sons. Inc., New York.

#Consider the study summarized in Example 5.4 of McGuire et al. (2000), a 3D space truss.


#Define structural system.
nodes = [[2000.0, 4000.0, 8000.0], [0.0, 0.0, 0.0], [8000.0, 0.0, 0.0], [8000.0, 6000.0, 0.0], [0.0, 6000.0, 0.0]]

cross_sections = InstantFrame.CrossSections(["member_ab", "member_ac", "member_ad", "member_ae"], [20E3, 30E3, 40E3, 30E3], [], [], [])

materials = InstantFrame.Materials(["steel"], [200000.], [0.3])

connections = InstantFrame.Connections(["semi-rigid", "rigid", "pinned, twist fixed"], [30000.0, Inf, Inf], [Inf, Inf, 0.0], [Inf, Inf, 0.0])

members = InstantFrame.Members(["frame", "frame"], [1, 1, 1, 1], [2, 3, 4, 5], [0.0, 0.0, 0.0, 0.0, 0.0], ["pinned, twist fixed", "pinned, twist fixed", "pinned, twist fixed", "pinned, twist fixed"], ["pinned, twist fixed", "pinned, twist fixed", "pinned, twist fixed", "pinned, twist fixed"], ["member_ab", "member_ac", "member_ad", "member_ae"], ["steel", "steel", "steel", "steel"])

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
