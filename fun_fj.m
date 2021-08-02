function r = fun_fj(M)
syms a x1 x2 x3 x4 x5 real
  A=[0 a 1 -1 -a;
    -a 0 a 1 -1;
    -1 -a 0 a 1;
    1 -1 -a 0 a;
    a 1 -1 -a 0]
 
% A=  [[ 0 0 a 0 1 ]
%      [ 1 0 0 a 0 ]
%      [ 0 1 0 0 a ]
%      [ a 0 1 0 0 ]
%      [ 0 a 0 1 0 ]]
    eigen_cycle=[];
    payoff_matrix=A;
       Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
         mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
        V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
        Jac =  [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
                 diff(V_Eq_0,'x4') diff(V_Eq_0,'x5') ]
%         S=solve(V_Eq_0);
%         A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);

%             s3=prod(A,2);
%             id = find(s3>0);
%  
%                 x1=A(id,1);x2=A(id,2);x3=A(id,3);x4=A(id,4);x5=A(id,5);
                x1=0.2;x2=0.2;x3=0.2;x4=0.2;x5=0.2;
                D_Eq_at_NE = vpa(eval(Jac));
               [eigen_vector, eigen_value] = eig(D_Eq_at_NE); 
               d_eigen_value =diag(eigen_value);
%                save('F:\tmpdata\result_ana202107\2centerTheory.mat')
                   for kv=1:5
                       Lmn= from_eigenvector_out_am(eigen_vector(:,kv));
                       eigen_cycle = [eigen_cycle Lmn];
                   end 
end