% code for Strategy Space Collapse: Experiment and Theory
% Yijia Wang* and Zhijian Wang†
% April 2019
% http://gtcenter.org/Archive/2019/Conf/Wang3311.pdf
% 
% for K=1:3;Ctime(K);end
function r=Ctime() %(Exp1Dyn2Sim3)
    P2=[];
    for Exp1Dyn2Sim3=1:3;
        if Exp1Dyn2Sim3==1; load('PE26.mat'); end %Experiment  
        if Exp1Dyn2Sim3==2; load('PD26.mat');  Tm=[]; for k=1:12; Tm=[Tm; PS];end;PS=Tm; end %Dynamics
        if Exp1Dyn2Sim3==3; load('PS26.mat'); end %ABED Simulation
        [CrossPoint4Dim T6accu T8accu T10accu] = getCrossPoint4(PS);
        matrix2latex(CrossPoint4Dim,strcat('tmp',num2str(Exp1Dyn2Sim3),'.tex'));
           P1=[];  %output distrribution
           for Pa=1:3
                 if Pa==1;   r=T6accu;  Accu_P1=r(:,1:8);Accu_P2=r(:,9:16);  end
                 if Pa==2;   r=T8accu;  Accu_P1=r(:,1:8);Accu_P2=r(:,9:16);  end
                 if Pa==3;   r=T10accu; Accu_P1=r(:,1:8);Accu_P2=r(:,9:16);  end
                        subFigId = 1+(Pa-1)*2;
                        subplot(3,6, subFigId + (Exp1Dyn2Sim3-1)*6);plot(Accu_P1,'linewidth',2);title(strcat('X: ',num2str(Pa))) 
                            legend('1','2','3','4','5','6','7','8','Location','northwest');xlim([0,600]);ylim([0,300]); axis square;  
                                if subFigId ~= 1; set(gca,'ytick',[]); end
                                if subFigId ~= 6;  legend off; end ; set(gca,'xtick',[]); 
                            xlabel('Period','FontSize',15), ylabel('Accumulated Frequency','FontSize',15)
                        subplot(3,6, subFigId+1+ (Exp1Dyn2Sim3-1)*6);plot(Accu_P2,'linewidth',2);title(strcat('Y:  ',num2str(Pa))) 
                            legend('1','2','3','4','5','6','7','8','Location','eastoutside');xlim([0,600]);ylim([0,300]); axis square;  
                            xlabel('Period','FontSize',15) %, ylabel('Accumulated Frequency')
                                 if subFigId+1 ~= 1; set(gca,'ytick',[]); end; set(gca,'xtick',[]); 
                                 legend off;
                  P1=[P1; [(r(1000,1:8)-r(500,1:8)); (r(1000,9:16)-r(500,9:16))]'/500];
           end
                P2 = [P2 P1];
    end
                NE = FindNash8x8_20200901;
                P2 = [[ones(1,8) ones(1,8)*2 ones(1,8)*3]' [1:8 1:8 1:8]' roundn(P2,-3) roundn(NE,-3)];
                matrix2latex(P2,'Exp_Dyn_Sim_NE_distribution.tex')
end

function [CrossPoint4Dim T6accu T8accu T10accu] = getCrossPoint4(PS)
        T6=PS(find(PS(:,1)==6),:);
        T8=PS(find(PS(:,1)==8),:);
        T10=PS(find(PS(:,1)==10),:); 
%% 6
Le=length(T6(:,1))/12;
T6poll=zeros(Le,16); T6accu=zeros(Le,16); 
for i=1:Le;
    a=T6(i:Le:Le*12,[4:11 13:20]);
    T6poll(i,:)=mean(a);
    T6accu(i,:)=sum(T6poll(1:i,:));
end 
    T6accu(1,:)=(T6poll(1,:));

    Crossover_si_sj_time_strength=[]; 
    for i=1:16; 
        for j=1:16; 
            tmp=T6accu(:,i)-T6accu(:,j);
               for k=length(T6accu(:,1)):-1:2
                   if tmp(k)*tmp(k-1) < 0 && (i-8.1)*(j-8.1)>0  
                       Crossover_si_sj_time_strength= ...
                            [Crossover_si_sj_time_strength; ...
                            i j k T6accu(k,i) sign(tmp(k))];
                       break
                   end  
               end 
        end;
    end

A6 = sortrows(Crossover_si_sj_time_strength,-4) ;



%% 8
T8poll=zeros(Le,16); T8accu=zeros(Le,16); 
for i=1:Le;
    a=T8(i:Le:Le*12,[4:11 13:20]);
    T8poll(i,:)=mean(a);
    T8accu(i,:)=sum(T8poll(1:i,:));
end 
    T8accu(1,:)=(T8poll(1,:));

    Crossover_si_sj_time_strength=[]; 
    for i=1:16; 
        for j=1:16; 
            tmp=T8accu(:,i)-T8accu(:,j);
               for k=length(T8accu(:,1)):-1:2
                   if tmp(k)*tmp(k-1) < 0 && (i-8.1)*(j-8.1)>0  
                       Crossover_si_sj_time_strength= ...
                            [Crossover_si_sj_time_strength; ...
                            i j k T8accu(k,i) sign(tmp(k))];
                       break
                   end  
               end 
        end;
    end

A8 = sortrows(Crossover_si_sj_time_strength,-4) ;



%% 10
T10poll=zeros(Le,16); T10accu=zeros(Le,16); 
for i=1:Le;
    a=T10(i:Le:Le*12,[4:11 13:20]);
    T10poll(i,:)=mean(a);
    T10accu(i,:)=sum(T10poll(1:i,:));
end 
    T10accu(1,:)=(T10poll(1,:));

    Crossover_si_sj_time_strength=[]; 
    for i=1:16; 
        for j=1:16; 
            tmp=T10accu(:,i)-T10accu(:,j);
               for k=length(T10accu(:,1)):-1:2
                   if tmp(k)*tmp(k-1) < 0 && (i-8.1)*(j-8.1)>0  
                       Crossover_si_sj_time_strength= ...
                            [Crossover_si_sj_time_strength; ...
                            i j k T10accu(k,i) sign(tmp(k))];
                       break
                   end  
               end 
        end;
    end
A10 = sortrows(Crossover_si_sj_time_strength,-4) ;

Sig=A6;
Si6_X = Sig(find((Sig(:,2)==2 | Sig(:,2)==6)&(Sig(:,1)~=2 & Sig(:,1)~=6 )) ,:)
Si6_Y = Sig(find((Sig(:,2)==2+8 | Sig(:,2)==4+8)&(Sig(:,1)~=2+8 & Sig(:,1)~=4+8 )) ,:)
Sig=A8;
Si8_X = Sig(find((Sig(:,2)==1)&(Sig(:,1)~=1)) ,:)
Si8_Y = Sig(find((Sig(:,2)==4+8)&(Sig(:,1)~=4+8)) ,:)
Sig=A10;
Si10_X = Sig(find((Sig(:,2)==1 | Sig(:,2)==2 | Sig(:,2)==5)&(Sig(:,1)~=1 & Sig(:,1)~=2 & Sig(:,1)~=5)) ,:)
Si10_Y = Sig(find((Sig(:,2)==2+8 | Sig(:,2)==4+8)&(Sig(:,1)~=2+8 & Sig(:,1)~=4+8 )) ,:)
CrossPoint4Dim = [Si6_X;Si6_Y;Si8_X;Si8_Y;Si10_X;Si10_Y;];
end
