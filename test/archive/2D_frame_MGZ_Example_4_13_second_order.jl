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

#First order solution
u_1f = Ke_ff \ Ff

#Define global system displacement vector for the first order analysis.
num_nodes = size(nodes, 1)
num_dof_per_node = 3
u_1 = zeros(Float64, num_nodes * num_dof_per_node)
u_1[free_dof] .= u_1f

#Solve for first order internal forces.  Need axial forces for second order analysis.
@btime P = InstantFrame.calculate_internal_forces(nodes, cross_sections, materials, members, u_1)

#Calculate global geometric stiffness matrix.
@btime Kg = InstantFrame.define_global_geometric_stiffness_matrix(nodes, members, P)

#Partition the geometric stiffness matrix.
Kg_ff = Kg[free_dof, free_dof]

#Get the total stiffness matrix, elastic + geometric.
Kff = Ke_ff + Kg_ff

p = [Kff, Ff]

function residual(u, p)

    Kff, Ff = p

    Kff * u - Ff

end

u0 = u_1f
u0 = @SVector [u0[i] for i in eachindex(u0)]
probN = NonlinearProblem{false}(residual, u0, p)
@btime u_2f = solve(probN, NewtonRaphson(), tol = 1e-9)

