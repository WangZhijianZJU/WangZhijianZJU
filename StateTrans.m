
%% 计算转移矩阵 TransMatrixOrig， TransMatrixAsym
%   + 通过TransMatrixAsym中净转移次数 以constTran下限，求出三元圈 Loop3 
%   过滤Loop3 中的等价圈，求得 Loop_cleaned 圈数（含对应的态） 
  YiJiaDesign = [6 8 10];
  constTran = 5; 
ret = []; 
for ParaID=1:3
    Para = YiJiaDesign(ParaID);          
    [n8x8list,stateList,R] = Get8x8Data(Para);
                     % count 1st and 2nd half vortex  start                       
%                          R=R(find(R(:,3) > 500),:);   
%                          R=R(find(R(:,3) < 501),:); 
                        R=sortrows(R,[1,2,3]);
                     % count 1st and 2nd half N-flux  end 
    TransMatrixOrig = zeros(64);
    for I = 1:64
            for J = 1:64
                TransMatrixOrig(I,J) = length(R(find(R(:,10) == stateList(I) & R(:,11) == stateList(J)),1)); 
            end
    end
                TransMatrixAsym = TransMatrixOrig - TransMatrixOrig';
                [sum(sum(abs(TransMatrixAsym))) sum(sum(TransMatrixAsym))]

                Loop3=[];
    for I=1:64
        for J=1:64
            for K = 1:64
                if TransMatrixAsym(I,J) > constTran & TransMatrixAsym(J,K) > constTran & TransMatrixAsym(K,I) > constTran 
                    Loop3=[Loop3; stateList(I) stateList(J) stateList(K) ...
                        TransMatrixAsym(I,J) TransMatrixAsym(J,K) TransMatrixAsym(K,I) ...
                        TransMatrixAsym(I,J)+TransMatrixAsym(J,K)+TransMatrixAsym(K,I) ...
                         999 Para constTran];
                end
            end
        end
    end


    Loop_cleaned = [Loop3(1,1:7)  min(Loop3(1,4:7)) 999 Para constTran ];
    L6 = [Loop3(:,1:3) Loop3(:,1:7)];
    for I = 2:length(L6(:,1)) 
        temp=0;
          for K=1:length(Loop_cleaned(:,1))
              for J=2:4
                  if L6(I,J:J+2) == Loop_cleaned(K,1:3)
                      L6(I,1:6) = 999*ones(1,6); 
                      temp=1;
                  end
              end
          end
         if temp==0; Loop_cleaned = [Loop_cleaned; L6(I,4:10) min(L6(I,7:9)) 999 Para constTran ];end
    end
    Loop4 = Loop_cleaned;
    ret = [ret; Loop_cleaned]
end

fLatex = [ret(:,10) ret(:,1:8)]
% 
% fuction [n8x8list R06 R08 R10] = Get8x8Data()
% n8x8list = load('C:\Users\Think\Downloads\n8x8.csv');
% n8x8list(:,1) = n8x8list(:,1) * 100 - 8000000;
% 
% %构造 State_t0,State_tp1(+1,后一轮),State_tm1(-1,前一轮)
% n8x8list(:,10) = n8x8list(:,6)*10 + n8x8list(:,7);
% stateList=round(unique(n8x8list(:,10)));
% 
% n8x8list(:,11:12) = 999;
% for L = 1:72000-1; n8x8list(L,11) = n8x8list(L+1,10);end
% for L = 1+1:72000; n8x8list(L,12) = n8x8list(L-1,10);end
% 
% s=round(unique(n8x8list(:,1)));
% %分参数，每个参数有24000记录
% R06 = n8x8list(find(n8x8list(:,9) == 6),:);
% R08 = n8x8list(find(n8x8list(:,9) == 8),:);
% R10 = n8x8list(find(n8x8list(:,9) == 10),:);
% end


