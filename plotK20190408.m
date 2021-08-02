function r = plotK20190408()
clf
Para = [2 1; 3 2; 4 2]; 
     r32=[];  PS=[];
     dt=0.01;Lambda=50;totalStep =1000; subjectsNumber = 12;
    for Pa=1:3

    % 定下参数   Lambda=50 dt=0.02
    %         r = logit_dyn_0909(payoff88(Para(Pa,1),Para(Pa,2)),Lambda,totalStep);
    %         r = logit_dyn_090A(payoff88(Para(Pa,1),Para(Pa,2)),Lambda,totalStep);
             r = logit_dyn_090B20190408(payoff88(Para(Pa,1),Para(Pa,2)),Lambda,totalStep,dt);
                Accu_P1=r(:,17:24)*subjectsNumber;
                Accu_P2=r(:,25:32)*subjectsNumber;
        %     Accu_P1=r(:,17-16:24-16)*subjectsNumber;
        %     Accu_P2=r(:,25-16:32-16)*subjectsNumber;
                subFigId = 1+(Pa-1)*2;
                subplot(1,6, subFigId);plot(Accu_P1,'linewidth',2);title(strcat('X: ',num2str(Pa))) 
                    legend('1','2','3','4','5','6','7','8','Location','northwest');xlim([0,600]);ylim([0,3200]); axis square; grid on
                        if subFigId ~= 1; set(gca,'ytick',[]); end
                        if subFigId ~= 6;  legend off; end ; set(gca,'xtick',[]);
                    %xlabel('Period','FontSize',15), ylabel('Accumulated Frequency','FontSize',15)
                subplot(1,6, subFigId+1);plot(Accu_P2,'linewidth',2);title(strcat('Y:  ',num2str(Pa))) 
                    legend('1','2','3','4','5','6','7','8','Location','eastoutside');xlim([0,600]);ylim([0,3200]); axis square;  grid on
%                     xlabel('Period','FontSize',15) %, ylabel('Accumulated Frequency')
                         if subFigId+1 ~= 1; set(gca,'ytick',[]); end; set(gca,'xtick',[]);
                         legend off;
        r32 =[r32; r(1000,:) Para(Pa,:) Lambda];
        TreatID=[6 8 10]
        PS = [PS; TreatID(Pa)*ones(1000,1) ones(1000,1) [1:1000]' ,r(:,1:8),[1:1000]',r(:,9:end)];
    end
    
                            set(gcf,'papersize',[20 3],'paperposition',[3,3,20,3]);hold on; 
                            saveas(gcf,'D:\MATLAB\R2016a\bin\html\wzj20180131\Theo600.png') 
                            
                            save('F:\cDownload\PokerAnan20200716\8x8abedata\PD26.mat','PS')
end


function anan=logit_dyn_090B20190408(A,Lambda,repeat,wlogit)
%%logit_dyn_0909(payoff88(2,1),6,100)  Pa=6;Lambda=0.48;totalStep=100;logit_dyn_090B(payoff88(Para(Pa,1),Para(Pa,2)),Lambda,totalStep,0.01);
%%logit_dyn_090B(payoff88(2,1),6,100)
%%anan=[p1-8 q9-16 p_sum17-24 q_sum25-32];

pay_a=A/6;
pay_b=-A'/6;
%%
%
%初始
p=[0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125];
q=[0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125];

Ua=zeros(1,8);
Ub=zeros(1,8);

p_sum=zeros(1,8);
q_sum=zeros(1,8); 

anan=[];

    for ri=1:repeat
        for i=1:8
            Ua(i)=sum(q.*pay_a(i,:));
            Ub(i)=sum(p.*pay_b(i,:));	
        end
        p_sum=p_sum+p;
        q_sum=q_sum+q;		

        c_a=exp(Lambda*Ua);
        c_b=exp(Lambda*Ub);

    %     p=c_a./sum(c_a)*0.1 + p.*0.9;
    %     q=c_b./sum(c_b)*0.1 + q.*0.9;
        p=c_a./sum(c_a)*wlogit + p.*(1-wlogit);
        q=c_b./sum(c_b)*wlogit + q.*(1-wlogit);

        anan=[anan;p q p_sum q_sum];
    end
end

 
function P=payoff8820190408(n,m)

P=zeros(8,8);

Stra=zeros(8,3);

for i=0:7
    tem=dec2bin(i,3);
    Stra(i+1,1)=str2num(tem(1));
    Stra(i+1,2)=str2num(tem(2));
    Stra(i+1,3)=str2num(tem(3));
end

situa=[];
situa_n=0;
for i=1:3
    for j=1:3
        if i==j
            continue
        end
        situa=[situa;i j];
        situa_n=situa_n+1;
    end
end

for i=1:8%for Ann %算Anan的收益
    for j=1:8 %for Bob
        for k=1:situa_n
            if Stra(i,situa(k,1))==0
                P(i,j)=P(i,j)+sign(situa(k,1)-situa(k,2))*m;
            else
                if Stra(j,situa(k,2))==0
                    P(i,j)=P(i,j)+1;
                else
                    P(i,j)=P(i,j)+sign(situa(k,1)-situa(k,2))*n;
                end
            end
        end
    end
end


end

   
