
% D:\MATLAB\R2016a\bin\tmp202106122.m

% % 6.5.2 A simple HIV model 
% syms alpha c1 p1 beta1 beta2 c2 p2
%   c1dot  = -alpha*c1 + p1*(c1+beta1*c2) ;
%   p1dot  = p1*(1-c1-beta1*c2) ;
%   c2dot  = -alpha*c2 + p2*(c2+beta2*c1) ;
%   p2dot  = p2*(1-c2-beta2*c1) ;
% A=solve([c1dot ,p1dot ,c2dot ,p2dot ])
% % reach There are four equilibria. The first one corresponds to extinction of the entire
% % population, the second and third to extinction of one subpopulation or the other,
% % and the last to coexistence of both subpopulations.
% % We will now determine the stability of these four equilibria. Note that to calculate
% % the Jacobian, we have to order the variables, and maintain the same order
% % throughout the calculation. 
% jacobian([c1dot,p1dot,c2dot,p2dot],[c1,p1,c2,p2])  
% %   
% % syms x_1 x_2 x_3 x_4 a a_1 a_2
% % x_1  =  a/(3*a+1) + a_1 * (a - 1)/4 + a_2
% % x_2  =  a/(3*a+1) + a_1 * (- a - 1)/2
% % x_3  =  a/(3*a+1) + a_1 * (a - 1)/4 - a_2
% % x_4  =  1/(3*a+1) + a_1
% % 
% % x_2_2 = subs(x_2,'a_1','x_4 - 1/(3*a+1)')  % a/(3*a + 1) - ((x_4 - 1/(3*a + 1))*(a + 1))/2
% % x_1_2 = subs(x_1,'a_1','x_4 - 1/(3*a+1)')  % x_1/2 - x_3/2 + a/(3*a + 1) + ((x_4 - 1/(3*a + 1))*(a - 1))/4
% % x_1_3 = subs(x_1_2,'a_2','(x_1-x_3)/2') 
%   
% 
% x2 = a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2
% x1 = 2*(- x3/2 + a/(3*a + 1) + ((_4 - 1/(3*a + 1))*(a - 1))/4)

function [T,Y] = tmp202106122()
global a
Ym=[];
Yr=[];
noise_Amp =0.01
        for i=1:1000
            [T,Y] = tmp202106122old(noise_Amp);
            [ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(Y,mean(Y));
            Yr=[Yr  Yret3];
            Ym=[Ym; mean(Y)];
        end
end


 


function [T,Y] = tmp202106122old(noise_Amp) 

rr=[];
tspan = [0 100];  
initial_pos = [1 1 1 1 2 1 1 1 ]/4 + noise_Amp * rand(1,8);
initial_pos = [initial_pos(1:4)/sum(initial_pos(1:4))] % initial_pos(5:8)/sum(initial_pos(5:8)) ];

for c=1:1;
    b=[0.25 1 4];
    a=b(c);
% [T,Y] = ode15s(@osc,tspan,initial_pos);
% [T,Y] = ode23(@osc,tspan,initial_pos);
[T,Y] = ode45(@osc,tspan,initial_pos);
%         plot(T,Y(:,1),'.-')
%         figure;plotmatrix(Y); title(num2str(a)) 
end
        function dydt = osc(t,y) 
            dydt = [ -(y(1)*(y(1)*y(2) - a*y(4) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4)))/(y(1)*y(2) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4))
                     -(y(2)*(y(1)*y(2) - y(1) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4)))/(y(1)*y(2) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4))
                     -(y(3)*(y(1)*y(2) - y(2) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4)))/(y(1)*y(2) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4))
                     -(y(4)*(y(1)*y(2) - y(3) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4)))/(y(1)*y(2) + y(2)*y(3) + y(3)*y(4) + a*y(1)*y(4))]; 
        end
end 


%% III. Solving systems of first-order ODEs
% function [T,Y] = call_osc()
% tspan = [0 30000];
% y1_0 = 2;
% y2_0 = 0;
% [T,Y] = ode15s(@osc,tspan,[y1_0 y2_0]);
% plot(T,Y(:,1),'o')
% 
%     function dydt = osc(t,y)
%     dydt = [y(2)
%            1000*(1 - y(1)^2)*y(2) - y(1)];
%     %Still y(1) is y1 and y(2) is y2, and dydt(1)
%     %is dy1/dt and dydt(2) is dy2/dt.
%     end
% 
% end
% 

%% II. Solving first-order ODEs
% function [t,y] = call_dstate()
% % https://web.mit.edu/voigtlab/BP205/Notes/BP205_Matlab_slides.pdf
%     tspan = [0 9]; % set time interval
%     y0 = 10; % set initial condition
%     % dstate evaluates r.h.s. of the ode
%     [t,y] = ode45( @dstate ,tspan ,y0);
%     plot(t,y)
%     disp([t,y]) % displays t and y(t)
%             function dydt = dstate (t,y)
%             alpha = 2; gamm a= 0.0001;
%             dydt = alpha* y - gamma *y^2;
%             end
% end