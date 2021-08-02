
function [f4set disMarkov]=TheoVortecInMarkov20190419()
Pa=[2 1;3 2;4 2]; f4set=[]; f3set=[]; 
    for K=1:3
                     A = payoff88(Pa(K,1),Pa(K,2));
                     [fLatex3 fLatex4 disMarkov]=TheoVortecInMarkov20190419emb(A);
                     f3set=[f3set;fLatex3];
                     f4set=[f4set;fLatex4];
    end
    matrix2latex(f3set,'Dyn_vortex3.tex')
    matrix2latex(f4set,'Dyn_vortex4.tex')
end
function [fLatex3 fLatex4 disMarkov]=TheoVortecInMarkov20190419emb(A)
                 Lambda=98;  
                 wlogit=0.1; 
                 repeat= 100; 
%     [M64x64 v disMarkov]= GenTheoMarkoW(A,Lambda,repeat,wlogit);
    [M64x64 v disMarkov]= GenTheoMarkoW20190408(A,Lambda,repeat,wlogit);
                f=v*12000;
                Para=A(8,7);
                constTran=5;
    T64 = zeros(64); for I=1:64; T64(I,:) = f(I)*M64x64(I,:);end

    TransMatrixOrig = T64;
    TransMatrixAsym = TransMatrixOrig - TransMatrixOrig'; 

% find 3 points vortex  
statelist=stateList();ret3=[];
                Loop3 = f3vortex(TransMatrixAsym,stateList,Para,constTran);  
                Loop3_cleaned =  CleanedLoop3(Loop3,Para,constTran);
                ret3 = [ret3; Loop3_cleaned];
                   fLatex3 = [ret3(:,10) ret3(:,1:8)];
                   fLatex3 = sortrows(fLatex3,-9);
% find 4 points vortex       
statelist=stateList();ret4=[];
                Loop4 = f4vortex(TransMatrixAsym,statelist,Para,constTran);  
                Loop4_cleaned =  CleanedLoop4(Loop4,Para,constTran);
                ret4 = [ret4; Loop4_cleaned];

                fLatex4 = [ret4(:,12) ret4(:,1:10)];
                fLatex4 = sortrows(fLatex4,-11);
                fLatex4 = fLatex4(1:end,:);
end



function ret4vortex = f4vortex(TransMatrixAsym,stateList,Para,constTran)
                Loop4=[];
    for I=1:64
        for J=1:64
            for K = 1:64
            for L = 1:64
                if TransMatrixAsym(I,J) > constTran & TransMatrixAsym(J,K) > constTran ...
                        & TransMatrixAsym(K,L) > constTran & TransMatrixAsym(L,I) > constTran 
                    Loop4=[Loop4; stateList(I) stateList(J) stateList(K) stateList(L) ...
                        TransMatrixAsym(I,J) TransMatrixAsym(J,K) TransMatrixAsym(K,L) TransMatrixAsym(L,I) ...
                        TransMatrixAsym(I,J)+TransMatrixAsym(J,K)+TransMatrixAsym(K,L)+TransMatrixAsym(L,I) ...
                        min([TransMatrixAsym(I,J) TransMatrixAsym(J,K) TransMatrixAsym(K,L) TransMatrixAsym(L,I)]) ...
                         999 Para constTran];
                end
            end
            end
        end
    end
                ret4vortex = Loop4;
end


function retCleanedLoop4 = CleanedLoop4(Loop4,Para,constTran)
Loop_cleaned = [];
    if  isempty(Loop4) == 0; 
            Loop_cleaned = [Loop4(1,1:9)  min(Loop4(1,5:8)) 999 Para constTran ];
            L6 = [Loop4(:,1:4) Loop4(:,1:9)];
            for I = 2:length(L6(:,1)) 
                temp=0;
                  for K=1:length(Loop_cleaned(:,1))
                      for J=2:5
                          if L6(I,J:J+3) == Loop_cleaned(K,1:4)
                              L6(I,1:8) = 999*ones(1,8); 
                              temp=1;
                          end
                      end
                  end
                 if temp==0; 
                     Loop_cleaned = [Loop_cleaned; L6(I,5:13) min(L6(I,9:12)) 999 Para constTran ];
                 end
            end
    end
 retCleanedLoop4 = Loop_cleaned;
end

function retCleanedLoop3 = CleanedLoop3(Loop3,Para,constTran)
Loop_cleaned = [];
    if  isempty(Loop3) == 0;
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
    end
 retCleanedLoop3 = Loop_cleaned;
end

function ret3vortex = f3vortex(TransMatrixAsym,stateList,Para,constTran)
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
ret3vortex = Loop3;

end


function r=stateList()
r=[11;12;13;14;15;16;17;18;21;22;23;24;25;26;27;28; ...
    31;32;33;34;35;36;37;38;41;42;43;44;45;46;47;48; ...
    51;52;53;54;55;56;57;58;61;62;63;64;65;66;67;68; ...
    71;72;73;74;75;76;77;78;81;82;83;84;85;86;87;88];
end