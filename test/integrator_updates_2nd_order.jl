using NonlinearSolve

p = 2.0
function residual(u, p) 
    
    p = 1.0*u    
    u .* u .- 2.0*p

end


u0 = 1.5
probB = NonlinearProblem(f, u0, p)
cache = init(probB, NewtonRaphson()); # Can iterate the solver object
solver = solve!(cache)

step!(cache)

for i in take(cache,3) end