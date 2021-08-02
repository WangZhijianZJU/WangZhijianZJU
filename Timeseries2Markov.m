

load('PS26.mat')
stateList=[11;12;13;14;15;16;17;18;21;22;23;24;25;26;27;28; ...
    31;32;33;34;35;36;37;38;41;42;43;44;45;46;47;48; ...
    51;52;53;54;55;56;57;58;61;62;63;64;65;66;67;68; ...
    71;72;73;74;75;76;77;78;81;82;83;84;85;86;87;88];

 TransMatrixOrig = zeros(64);r3vor=[];r4vor=[];
  YiJiaDesign = [6 8 10];  
    for ParaID=1:3
        Para = YiJiaDesign(ParaID);    
T8=PS((find(abs(PS(:,3)-40)>0.001 & PS(:,1)==Para)),[4:11 13:20]);
T8(:,17:32)=99;
for i=1:12000; for j=1:8; T8(i,16+j)=sum(T8(i,1:j));T8(i,24+j)=sum(T8(i,9:8+j));T8(i,33:34)=rand(1,2);end;end
for i=1:12000; T8(i,35:36) = [min(find(T8(i,17:24)>T8(i,33))) min(find(T8(i,25:32)>T8(i,34)))] ;end;
for i=1:12000; T8(i,37) =  T8(i,35)*10+ T8(i,36);end;
for i=1:12000-1; T8(i,38) =  T8(i+1,37);end;


    for I = 1:64
            for J = 1:64
                TransMatrixOrig(I,J) = length(T8(find(T8(:,37) == stateList(I) & T8(:,38) == stateList(J)),1)); 
            end
    end
             M=TransMatrixOrig;
             subplot(2,3,ParaID);Mf = flipud(M);pcolor(Mf);figure(gcf);colorbar  %('Ticks',[0:500:2500])
             subplot(2,3,ParaID+3); A=(abs(M-M')+(M-M'))/2;Af=flipud(A);pcolor(Af);figure(gcf);colorbar

        [r3vortex,r4vortex] = Markov2Vortex(TransMatrixOrig)
        r3vortex(:,1)=ones(length(r3vortex(:,1)),1)*Para;
        r4vortex(:,1)=ones(length(r4vortex(:,1)),1)*Para;
        r3vor=[r3vor;r3vortex];
        r4vor=[r4vor;r4vortex]; 
    end
    r3vor