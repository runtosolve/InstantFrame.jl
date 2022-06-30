using InstantFrame, NonlinearSolve, Plots, LinearAlgebra, BenchmarkTools, StaticArrays

#Denavit, M. D., & Hajjar, J. F. (2013). Description of geometric nonlinearity for beam-column analysis in OpenSees.

#Consider the study summarized in Figure 5 of Denavit and Hajjar (2013), a cantilever under axial and transverse loading.

#Define frame inputs.
E = 29000.0  #ksi
I = 37.1 #in^4
A = 9.12 #in^2
L = 180.0 #in
P = -50.0 #kips
H = -1.0 #kips
ρ = 492.0 / 32.17 / 1000.0 / 12^4 #kips * s^2 / in^4
 
#Define the frame loading as an external force vector. There is an axial load of 50 kips and a lateral load of 1 kip.
F = [0., 0., 0., P, H, 0.0]

#Define the frame boundary conditions.
free_dof = [4; 5; 6]

#Define the frame elastic stiffness.
Ke = InstantFrame.define_local_elastic_stiffness_matrix(I, A, E, L)
Kg = InstantFrame.define_local_geometric_stiffness_matrix(P, L)

K=Ke+Kg


Ku=F 

Ku-F = 0 



#Partition external force vector and elastic stiffness matrix.
Ff = F[free_dof]
Ke_ff = Ke[free_dof, free_dof]
Kg_ff = Kg[free_dof, free_dof]

#Get the total stiffness matrix, elastic + geometric.
Kff = Ke_ff + Kg_ff

#Define the deformation initial guess for the nonlinear solver.
deformation_guess = Ke_ff \ Ff


p = [Kff, Ff]

function residual(u, p)

    Kff, Ff = p

    Kff * u - Ff

end



u0 = deformation_guess
u0 = @SVector [u0[i] for i in eachindex(u0)]
probN = NonlinearProblem{false}(residual, u0, p)
@btime solver = solve(probN, NewtonRaphson(), tol = 1e-9)








#Solve for the beam deformations.
@btime solution = nlsolve((R,U) ->InstantFrame.second_order_analysis_residual!(R, U, Kff, Ff), deformation_guess)

#Newton-Raphson is faster here, than trust region...


#Get the frame deformed shape.
u = zeros(Float64, 6)
u[free_dof] .= solution.zero

#Show the frame deformed shape.
x = range(0., L, 20)
offset = 0.
w = beam_shape_function(u[2],u[3],u[5],u[6], L, x, offset)
plot(x, w, linewidth = 10, linecolor = :gray, legend=false)


#Plot load-deformation response.
num_steps = 20
P_range = range(0., P, num_steps)
H_range = range(0., H, num_steps)
u_range = Array{Array}(undef, num_steps)

for i = 1:num_steps
    F = [0., 0., 0., P_range[i], H_range[i], 0.0]
    u_range[i] = get_load_deformation_response(I, A, E, L, P_range[i], F, free_dof)
end

u_load = [u_range[i][5] for i = 1:num_steps]
plot(-u_load, -P_range ./ -P, markershape = :o, xlabel="cantilever disp. [in.]", ylabel="load factor", label="varying P")

u_range_constant_P = Array{Array}(undef, num_steps)
for i = 1:num_steps
    F = [0., 0., 0., P, H_range[i], 0.0]
    u_range_constant_P[i] = get_load_deformation_response(I, A, E, L, P, F, free_dof)
end

u_load_constant_P = [u_range_constant_P[i][5] for i = 1:num_steps]
plot!(-u_load_constant_P, -P_range ./ -P, markershape = :o, xlabel="cantilever disp. [in.]", ylabel="load factor", label="constant P", legend=:bottomright)



#Calculate the vibration natural frequencies and mode shapes.
M = InstantFrame.define_local_element_mass_matrix(A,L,ρ)
Mff = m[free_dof, free_dof]

ωn_squared=eigvals(Ke_ff, Mff)
ωn=sqrt.(ωn_squared)

f = ωn/(2π)

#And the natural periods of vibration...
T = 1 ./ f

mode_shapes=eigvecs(Ke_ff,Mff)

#And show the vibration modes.
mode_shape = zeros(Float64, 6)
mode = 3
mode_shape[free_dof] .= modeshapes[:, mode] ./ maximum(abs.(modeshapes[:, mode]))
w = beam_shape_function(mode_shape[2],mode_shape[3],mode_shape[5],mode_shape[6], L, x, offset)
plot(x, w, linewidth = 10, linecolor = :gray, legend=false)







