% R =[]; V =[]; Ne=[]; 
%         s=perms(1:5);
% for t=1:10000 
syms x1 x2 x3 x4 x5   real
  payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0 ]
% a= 1 
% payoff_matrix = [ 0 0 0 0 a; 1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ]
%    payoff_matrix  = [0 -1  1  1 -1; 
%                      1  0 -1 -1  1;
%                     -1  1  0  1 -1;
%                     -1  1 -1  0  1;
%                      1 -1  1 -1  0 ]
%         t=randi([1 120],[5 1])
%         payoff_matrix = s(t,:)
 Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
 mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
 V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
 Jac = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
            diff(V_Eq_0,'x4') diff(V_Eq_0,'x5') ] ; 

    S=solve(V_Eq_0);
    A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);
        s=prod(A,2);
        id = find(s~= 0);
        if size(id,1)==1
                x1=A(id,1); x2=A(id,2); x3=A(id,3); x4=A(id,4); x5=A(id,5); 
%                 x1=0.2; x2=0.2; x3=0.2; x4=0.2; x5=0.2; 
                Jac_at_NE = eval(Jac)
                [eigen_vector eigen_value] = eig(Jac_at_NE)
%                           R=[R; payoff_matrix];
%                           V=[V; diag(eigen_value)'];
% %                           Ne=[Ne; s];
%                   figure
%                 v=diag(eigen_vector(:,1));hold on
%                 quiver([1:5]'*0,[1:5]'*0,real(v),imag(v),1)
%                 text(real(v)*1.05,imag(v)*1.05, num2str([1:5]'),'fontsize',15);axis square;
         end
% end
        v = eigen_vector(:,3)
        [Lmn Tmn]= from_eigenvector_out_am(v)
% save('manifold5try.mat', 'R','V','Ne' )