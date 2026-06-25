using InstantFrame

I = 0.55967
L = 48.0
E = 29000000.0
k1 = 230000.0
k2 = Inf 

ke_left_spring = InstantFrame.define_local_elastic_element_stiffness_matrix_partially_restrained(I, E, L, k1, k2)

using InstantFrame

I = 0.55967
L = 48.0
E = 29000000.0
k1 = Inf
k2 = 230000.0

ke_right_spring = InstantFrame.define_local_elastic_element_stiffness_matrix_partially_restrained(I, E, L, k1, k2)


I = 0.55967
L = 48.0
E = 29000000.0
k1 = 230000.0
k2 = 230000.0 * 10^30


α1 = k1/(E*I/L)  #start node

α2 = k2/(E*I/L)  #end node


Kbb = E*I/L* [4 + α1      2
              2        4+α2]

Kbc = E*I/L * [6/L  -α1   -6/L   0.0
               6/L   0    -6/L   -α2]

Kcb = Kbc'


Kcc = E*I/L * [12/L^2  0   -12/L^2  0
       0     α1    0     0
       -12/L^2  0   12/L^2  0
       0   0  0  α2 ]


Kcc_hat = Kcc - Kcb*Kbb^-1*Kbc 


α = (α1*α2)/(α1*α2 + 4*α1 + 4*α2 + 12)

k34 = α*E*I/L * (-6/L * (1 + 2/α1))

k44 = α*E*I/L * 4 * (1 + 3/α1)