using ApproxFun, Plots


E = 29000  #ksi
I = 37.1 #in^4
A = 9.12 #in^2
L = 180.0 #in
P = -50.0 #kips
H = 1.0  #kips


x = Fun(identity, 0..L)
# x = Fun(Taylor(), [1,2,3])
y₀ = 0.0x # initial guess
N = y -> [y(0), y'(0), y''(L), y'''(L) - H/(E*I), E*I*y'''' + P*y'']
y = newton(N, y₀) 
plot(y, legend = false)

rad2deg(y'(L))




s = Fun(identity, 0..L)
θ₀ = 0.0s # initial guess
N = θ -> [θ(0), θ'(L), E*I*θ'' + P*sin(θ) + H*L - H*x]
θ = newton(N, θ₀) 
plot(θ, legend = false)

y(L)


x = Fun(identity, 0..L)
y₀ = 0.0x # initial guess
N = y -> [y(0), y'(0), E*I*y'' + P*y + H*L - H*x]
y = newton(N, y₀) 
plot(y, legend = false)



α = sqrt(-P*L^2/(E*I))

y_theor = H*L^3/(3*E*I) * (3 * (tan(α) - α))/α^3

plot(y''')

plot(y''', legend=false)


H/(E*I)

y(L)





L, E, Iy, P, qx = 12.0, 1.0, 1.0, 2.0, 0.0001
d = 0..L
z = Fun(identity, d)
B = [[Evaluation(d,leftendpoint,k) for k=0:1]... ; [Evaluation(d,rightendpoint,k) for k=2:3]... ;]
v = [B; E*Iy*Derivative()^4 + P*Derivative()^2] \ [ zeros(4)..., qx*one(z)]
func_name = zip([v, v', v'', v'''], ["Deflection", "Angle", "Moment", "Shear"])
plot([plot(z, f, title=n, label="") for (f,n) in func_name]..., lw=3)


x = Fun()
B = Evaluation(Chebyshev(),-1)

A = [B      0;
            B*𝒟    0;
            0      B;
            𝒟^2-I  2I;
            I      𝒟+I]

 u,v = A\[0;0;0;exp(x);cos(x)]

 u(-1),u'(-1),v(-1)





Δ = qx * L^4/(8 * E * Iy)

maximum(u)

u.(0:0.1:12)


x = Fun(identity, 0..1)
N = (u1,u2) -> [u1'(0) - 0.5*u1(0)*u2(0);
                u2'(0) + 1;
                u1(1) - 1;
                u2(1) - 1;
                u1'' + u1*u2;
                u2'' - u1*u2]

u10 = one(x)
u20 = one(x)
u1,u2 = newton(N, [u10,u20])

plot(u1, label="u1")
plot!(u2, label="u2")