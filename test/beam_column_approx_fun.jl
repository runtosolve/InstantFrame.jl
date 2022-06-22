using InstantFrame


nodes = InstantFrame.Nodes([0.0, 10.0, 20.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

cross_sections = InstantFrame.CrossSections(["W16x26"], [7.68], [301.0], [9.59], [0.262])

materials = InstantFrame.Materials(["steel"], [29000.0], [0.3])

connections = InstantFrame.Connections(["semi-rigid", "rigid"], [30000.0, Inf], [Inf, Inf], [Inf, Inf])

members = InstantFrame.Members(["frame", "frame"], [1, 2], [2, 3], ["rigid", "rigid"], ["rigid", "rigid"], ["W16x26", "W16x26"], ["steel", "steel"])













using ApproxFun, Plots
L, E, I = 12.0, 1.0, 1.0
d = 0..L
z = Fun(identity, d)
B = [[Evaluation(d,leftendpoint,k) for k=0:1]... ; [Evaluation(d,rightendpoint,k) for k=2:3]... ;]
v = [B; E*I*Derivative()^4] \ [ zeros(4)..., one(z)]
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
