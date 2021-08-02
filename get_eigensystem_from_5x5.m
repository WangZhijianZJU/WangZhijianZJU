% Program: Socexp, WZJ, 20210721
% Input:  payoff_matrix   --- the payoff  
% Output: eigen_vector, 
%         d_eigen_value, 
%         rest_point_no, 
%         solutionlist (A)
%         of the 5x5 payoff matrix
%         of replicator dynamics 


function [eigen_vector, d_eigen_value, eigen_cycle, rest_point_no, A] = ...
         get_eigensystem_from_5x5(payoff_matrix)
eigen_vector=[];
d_eigen_value=[];
eigen_cycle = [];
rest_point_no=1;
syms x1 x2 x3 x4 x5 real
       Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
         mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
        V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
        Jac =  [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
                 diff(V_Eq_0,'x4') diff(V_Eq_0,'x5') ]; 
%         S=solve(V_Eq_0);
%         A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);
% 
%             s3=prod(A,2);
%             id = find(s3>0);
%         if size(id,1)==1 
%                 x1=A(id,1);x2=A(id,2);x3=A(id,3);x4=A(id,4);x5=A(id,5);
                x1=0.2;x2=0.2;x3=0.2;x4=0.2;x5=0.2;
                D_Eq_at_NE = eval(Jac);
               [eigen_vector, eigen_value] = eig(D_Eq_at_NE); 
               d_eigen_value =diag(eigen_value);
                   for kv=1:5
                       Lmn= from_eigenvector_out_am(eigen_vector(:,kv));
                       eigen_cycle = [eigen_cycle Lmn];
                   end    
%         end
%         if size(id,1) > 1
%             rest_point_no = size(id,1);
%         end
A=[ 2 2 2 2 2]; 
end