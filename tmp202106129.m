%payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0 ]
%payoff_matrix_1 = [0 -1 1 1 -1; 1 0 -1 -1 1;-1 1 0 1 -1;-1 1 -1 0 1;1 -1 1 -1 0 ]

R =[]; V =[]; Ne=[]; 
for t=1:100
    syms x1 x2 x3 x4 x5 real
        payoff_matrix =randi([1 5],5,5) ;
    %    payoff_matrix = [1,3,1,5,5;2,2,4,2,3;5,5,2,4,1;5,1,2,1,3;1,4,2,2,5]
    %    payoff_matrix = [1,2,3,4,5;1,3,4,5,1;4,1,2,5,1;5,1,4,1,4;2,5,5,1,1]
%                          [1,4,1,5,3;1,2,4,4,4;3,4,4,2,1;4,2,3,1,5;3,2,4,3,2];
%                          [5,1,2,1,5;3,1,3,4,1;1,3,2,3,4;4,5,2,3,2;1,1,3,3,3];
%                          [3,5,2,2,4;5,1,3,1,4;2,1,5,4,4;5,4,2,1,1;2,5,2,5,2] 
         s = bimat(payoff_matrix,payoff_matrix');
         if abs(sum(s)-1) < 0.001 & prod(s) > 0
                Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
                 mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
                V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
                Jac = jacobian(V_Eq_0); 

                S=vpasolve(V_Eq_0);
                A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);

                    s3=prod(A,2);
                    id = find(s3>0);
                    if size(id,1)==1
                            x1=A(id,1); x2=A(id,2); x3=A(id,3); x4=A(id,4); x5=A(id,5);
                            D_Eq_at_NE = eval(Jac);
                           [eigen_vector, eigen_value] = eig(D_Eq_at_NE);

                %         VI1 = diag(eigen_value);
                %         VI2 = sum(abs(sign(imag(VI1)))');
                %             if VI2 > 3
                         if prod(s,2) > 0
                                      R=[R; payoff_matrix];
                                      V=[V; diag(eigen_value)'];
                                      Ne=[Ne; s];
                         end 
                    end
         end
end
save('manifold0718try2.mat', 'R','Ne' )
