using InstantFrame








using ApproxFun, Plots
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



x = Fun(identity, 0..L)
u₀ = 0.1x # initial guess
L, E, Iy, P, qx = 12.0, 1.0, 2.0, -0.1, 0.0001
N = u -> [u(0), u'(0), u''(L), u'''(L), E*Iy*u'''' + P*u'' - qx]
u = newton(N, u₀) # perform Newton iteration in function space
plot(u, legend = false)

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