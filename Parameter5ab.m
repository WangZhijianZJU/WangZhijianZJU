 syms a b x1 x2 x3 x4 x5 real
 
% b=1;
rM5ab=[ a-b 1+b 2 3 4
        4 a-b 1+b 2 3
        3 4 a-b 1+b 2
        2 3 4 a-b 1+b
        1+b 2 3 4 a-b];
payoff_matrix = rM5ab
       Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
         mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
        V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U); 
 Jac = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
            diff(V_Eq_0,'x4') diff(V_Eq_0,'x5') ] ;  
                D_Eq_at_NE =  subs(subs(subs(subs(subs(Jac,'x1',0.2),'x2',0.2),'x3',0.2),'x4',0.2),'x5',0.2);
               [eigen_vector, eigen_value] = eig(D_Eq_at_NE); 
               d_eigen_value =diag(eigen_value);
               
               
               