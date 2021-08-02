R =[]; V =[]; Ne=[]; 
        s1=perms(1:5);
for t=1:10000 
    
            t=randi([1 120],[5 1]);
            payoff_matrix = s1(t,:);
         s = bimat(payoff_matrix,payoff_matrix');
     if abs(sum(s)-1) < 0.001 & prod(s) > 0
           syms x1 x2 x3 x4 x5   real
           Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
             mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
            V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
            Jac = jacobian(V_Eq_0); 

            S=vpasolve(V_Eq_0);
            A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);

                s3=prod(A,2);
                id = find(s3>0);
            if size(id,1)==1
                    x1=A(id,1);x2=A(id,2);x3=A(id,3);x4=A(id,4);x5=A(id,5);
                    D_Eq_at_NE = eval(Jac);
                   [eigen_vector, eigen_value] = eig(D_Eq_at_NE); 

                 if prod(s,2) > 0
                              R=[R; payoff_matrix];
                              V=[V; diag(eigen_value)'];
                              Ne=[Ne; s [x1 x2 x3 x4 x5]];
                      save('manifold0718try2.mat', 'R','V','Ne', '-append')
                 end 
            end
     end


end

