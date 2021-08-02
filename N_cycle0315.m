%%
%
%
%%
%
syms x1 x2 x3 x4 a real

%% ֧������
% payoff_matrix = [0 -1 0 1; 1 0 -1 0; 0 1 0 -1; -1 0 1 0];
payoff_matrix = [0 0 0 a;
                 1 0 0 0;
                 0 1 0 0;
                 0 0 1 0];

%% �ռ�������֧������
Payoff_vector_field_F = payoff_matrix * [x1 x2 x3 x4]';

%% �����֧����ֵ
mean_U = [x1 x2 x3 x4] * Payoff_vector_field_F;

%% ����ѧ����
V_F = [x1 x2 x3 x4]'.*(Payoff_vector_field_F - mean_U);

%% ȫ΢�־���
D_V_F = [diff(V_F,'x1') diff(V_F,'x2') diff(V_F,'x3') diff(V_F,'x4')];

%% ��0��
S=solve(V_F);

%% ������֪0��
W = subs(subs(subs(subs(D_V_F,'x3','x2'),'x4','x2'),'x1','x2/a'),'x2','a/(3*a + 1)')
U=simplify(W)
[v d]=eig(U)
%x1=1/4;x2=1/4;x3=1/4;x4=1/4;
% % % % eval(subs([v; d],'a',2))
% % % % % �����Ȿ��ֵ������������ʲ����֮��Ĺ�ϵ��
aK =[0.1:0.1:1 1:1:9]
a_NE_eigV_3col=[]
for K=1:length(aK)
%     for a=1:9
    a=aK(K) %1:9
        A=eval([S.x1 S.x2 S.x3 S.x4]);
        %A=A(2,:)
        for k=1:length(A(:,1))
            if prod(A(k,:)) > 0
                x1=A(k,1); x2=A(k,2); x3=A(k,3); x4=A(k,4);
        % x1=1.5/4;x2=.5/4;x3=1.5/4;x4=.5/4;

        %% ������0���ϵı���ֵ�ͱ�������
                [eigen_vector eigen_value] = eig(eval(D_V_F))
                [(payoff_matrix) A(k,:)']
                
                  b=eigen_vector(:,2);
                  normal_b = b*exp(-atan(imag(b(4))/real(b(4)))*i);
                  normal_b = normal_b*sign(real(b(4)));
                
          ret = [ eval(payoff_matrix(:,4)) A(k,:)' ... 
                    [eigen_value(1,1) eigen_value(2,2) ...
                     eigen_value(3,3) eigen_value(4,4)]'  b  normal_b]

            else
                'No inner mixed strategy'
            end
        end
        a_NE_eigV_3col = [a_NE_eigV_3col; ret]
end
%% �۲챾��������ֲ�
%   for j=1:4; 
%       [(b(:,j).*b(:,j))   b(:,1).*conj(b(:,1))]
%   end
%   c=ones(4)*8;
%   for j=1:4; 
%   for k=1:4; 
%        c(j,k)=(b(:,j)'*b(:,k))  
%   end  
%   end  

%% ��ʾ���ű��������ϵ�ƫ�룬ϵͳ�˶�������
  c2=eigen_vector(:,2)'; 
  % ���Ǿ��������ŵڶ������������ϣ�ƫ���� 0.01�� ϵͳ��������
  
  t2=real([[A(1,:) + 0.01*c2] 0])
  step=0.01
for k=0:1:10000
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*step;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*step;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*step;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*step;
    
  t2=[t2; real([x1 x2 x3 x4]) k];
    
end
figure(1)
    plot(real(t2(2:end,[1 2 3 4])),'linewidth',4)
figure(2)
    plot3(t2(:,[1]),t2(:,[2]),t2(:,[3]));hold on
    plot3(1/4,1/4,1/4,'r*')
    plot3(t2(1:100:1000,[1]),t2(1:100:1000,[2]),t2(1:100:1000,[3]),'r.-')
    text(t2(100,[1]),t2(100,[2]),t2(100,[3]),'100')
    xlabel('x');ylabel('y');zlabel('z');
% eigen_vector =
% 
%	-0.7071     -0.5	 -0.5	0.1692
%         0      0.5i	 -0.5i	0.6866
%	-0.7071      0.5      0.5 	0.1692
%         0     -0.5i	  0.5i	0.6866
%				
%				
% eigen_value =			
%				
%	-0.25           0        0       0
%       0        0.25i       0       0
%       0           0	 -0.25i      0
%       0           0        0	 -0.25

% û����������û����λ�Ĳ�ֵ�� ��1,4̬��13��24ͬ�࣬����23̬��13��24����

%% ���������ϵĽ��
% figure(3) 
%                 e1=[1 0 0];
%                 e2=[1/2 sqrt(3/4) 0];
%                        phi3=(pi/2-acos(...
%                            ((sqrt(3)/2)^2 *2 - 1)/ ((sqrt(3)/2)^2 *2) ...
%                                 ));
%                 e3=1/3*(e2+e1) + [0 0 cos(phi3)*sqrt(3/4)];
%  t3=[];
%          for m=1:5001
%                 tmp = real(t2(m,1:3))* [e1; e2; e3]';
%                 t3=[t3; tmp 1-sum(tmp) t2(m,5)];
%          end
%      plot3(t3(:,[1]),t3(:,[2]),t3(:,[3]));hold on
%      Nsah = [1/4,1/4,1/4]*[e1; e2; e3]';
%     plot3(Nsah(1),Nsah(2),Nsah(3),'r*')
%     plot3(t3(1:100:1000,[1]),t3(1:100:1000,[2]),t3(1:100:1000,[3]),'r.')
%     text(t3(100,[1]),t3(100,[2]),t3(100,[3]),'100')
%     xlabel('x');ylabel('y');zlabel('z'); 


 