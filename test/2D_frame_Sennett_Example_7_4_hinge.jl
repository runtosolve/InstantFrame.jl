using InstantFrame, NonlinearSolve, Plots, LinearAlgebra, BenchmarkTools, StaticArrays, LinearAlgebra, BenchmarkTools

#Sennett, R. E. (2000). Matrix analysis of structures. Waveland Press.

#Consider the study summarized in Example 7-4 of Sennett (2000), a 2D portal frame with a hinge.


#Define structural system.
nodes = [[0.0, 0.0], [0.0, 120.0], [120.0, 120.0], [120.0, 0.0]]

cross_sections = InstantFrame.CrossSections(["beam"], [10.0], [100.0], [], [])

materials = InstantFrame.Materials(["steel"], [29e6], [0.3])

connections = InstantFrame.Connections(["rigid", "hinge"], [Inf, 0.0], [], [])

members = InstantFrame.Members(["frame", "frame", "frame"], [1, 2, 3], [2, 3, 4], [0.0, 0.0, 0.0], ["rigid", "hinge", "rigid"], ["hinge", "rigid", "rigid"], ["beam", "beam", "beam"], ["steel", "steel", "steel"])

#Calculate global elastic stiffness matrix.
Ke = InstantFrame.define_global_elastic_stiffness_matrix(nodes, cross_sections, materials, members)
 
#Define the frame loading as an external force vector. 
F = [0., 0., 0., 1000.0, 0., 0., 0., 0., 0., 0., 0., 0.]

#Define the frame boundary conditions.
free_dof = [3; 4; 5; 6; 7; 8; 9; 12]

#Partition external force vector and elastic stiffness matrix.
Ff = F[free_dof]
Ke_ff = Ke[free_dof, free_dof]

#Solution
uf = Ke_ff \ Ff

u = zeros(Float64, 12)

u[free_dof] = uf




