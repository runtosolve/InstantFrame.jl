using InstantFrame

#use modified version of Sennett Eq. 7-54 to obtain element stiffness matrix including partial restraint at the left end.

I = 100.0
A = 10.0
E = 29E6
L = 120.0

ke = InstantFrame.define_local_elastic_element_stiffness_matrix(I, A, E, L)


k1 = 0.0 #lb-in./rad  this is consistent with a typical rack connection.
k2 = 3000000000000000000000.0 #lb-in./rad  


ke = InstantFrame.define_local_elastic_element_stiffness_matrix_condensed(I, A, E, L, k1, k2)


#From Zeimian Example 13.7

    α1 = k1/(E*I/L)
    α2 = k2/(E*I/L)

    Kbb = E*I/L*[4+α1   2
                2      4+α2]

    Kbc = E*I/L * [6/L  -α1     -6/L    0
                6/L  0       -6/L    -α2]


    Kcc = E*I/L*[12/L^2     0   -12/L^2     0
                0          α1  0       0
                -12/L^2    0   12/L^2  0
                0          0   0       α2]


    Kcc_bar = Kcc - Kbc'*Kbb^-1*Kbc


    #Check against Sennett Eq. 7.37

    3*E*I/L^3  #Kcc_bar[1,1] matches 




