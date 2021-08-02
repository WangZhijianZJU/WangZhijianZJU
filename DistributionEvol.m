%% 实验结果 

% C:\Users\Think\Downloads\【批量下载】8x8(82).rar等3个文件\DistributionEvol.m
%  
YiJiaDesign = [6 8 10];
ParaID=1;
DisEvolution =[];
for ParaID=1:3
    Para = YiJiaDesign(ParaID);
    [n8x8list,stateList,R] = Get8x8Data(Para);
%% 静态分布
Disoll = zeros(8,2); Dis1st = zeros(8,2); Dis2nd = zeros(8,2);Dis200Last = zeros(8,2);
    for I=1:8
        Disoll(I,1) = length(find(R(:,6)==I)); %T001 distribution
        Disoll(I,2) = length(find(R(:,7)==I)); %T002 distribution
        Dis1st(I,1) = length(find(R(:,6)==I & R(:,3) < 501 )); %T001 distribution
        Dis1st(I,2) = length(find(R(:,7)==I & R(:,3) < 501 )); %T002 distribution 
        Dis2nd(I,1) = length(find(R(:,6)==I & R(:,3) > 500 )); %T001 distribution
        Dis2nd(I,2) = length(find(R(:,7)==I & R(:,3) > 500 )); %T002 distribution 
        Dis200Last(I,1) = length(find(R(:,6)==I & R(:,3) > 800 ));
        Dis200Last(I,2) = length(find(R(:,7)==I & R(:,3) > 800 ));
    end
    DisEvolution = [DisEvolution; ones(8,1)*ParaID Disoll/12000 Dis1st/6000 Dis2nd/6000 Dis200Last/2400];
%% 累积分布-时间依赖
Accu_P1 = zeros(1000,8); %Accu(1,R(1,6))=1;
Accu_P2 = zeros(1000,8); %Accu(1,R(1,6))=1;
    for I = 1:1000
        for J=1:8
            Accu_P1(I,J) = length(find(R(:,6)==J & R(:,3) <= I));
            Accu_P2(I,J) = length(find(R(:,7)==J & R(:,3) <= I));
        end
    end
% 安安想用中位数，屏蔽个例
% 中位数的算法有二种，
%  Algorithm 1 -- given t, for subject j and strategy j, calc 
%  p_{i,j}(t-25:t+25), then calc
%  x = median(p_{i,1:12}(t-25:t+25))

%  Algorithm 2 -- given t, for subject j and strategy j, calc 
%  p_{i,t} = median(p_{i,1:12}) -- result must be 0 or 1, then calc
%  p_{i,t} = mean(p_{i,t-25:t+25})
%
%     Accu_P1
%     Accu_P2
   subFigId = 1+(ParaID-1)*2;
    subplot(1,6,subFigId);plot(Accu_P1,'linewidth',2);%title(strcat('X:',  num2str(ParaID))); 
%     legend('1','2','3','4','5','6','7','8','Location','northwest');
    xlim([0,600]);ylim([0,3200]);axis square
                         if subFigId ~= 1; set(gca,'ytick',[]); end
%     xlabel('Period','FontSize',15), ylabel('Accumulated Frequency','FontSize',15)
    subplot(1,6,subFigId+1);plot(Accu_P2,'linewidth',2);%title(strcat('Y:',num2str(ParaID)))
%     legend('1','2','3','4','5','6','7','8','Location','northwest');
    xlim([0,600]);ylim([0,3200]);;axis square
                         if subFigId+1 ~= 1; set(gca,'ytick',[]); end
%     xlabel('Period','FontSize',15) %, ylabel('Accumulated Frequency')
    
end
                            set(gcf,'papersize',[17 3],'paperposition',[3,3,20,3]);hold on; 
                            saveas(gcf,'D:\MATLAB\R2016a\bin\html\wzj20180131\Expe600.png')

