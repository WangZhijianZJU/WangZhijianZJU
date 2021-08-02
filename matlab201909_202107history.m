%-- 2019/9/18 18:47 --%
V2_OSA
[S, C] = graphconncomp(Grig_M, 'Directed', true)
[MaxFlow, FlowMatrix, Cut] = graphmaxflow(G, 12, 14)
[MaxFlow, FlowMatrix, Cut] = graphmaxflow(Grig_M, 12, 14)
P = [ 0   0  1/2 1/4 1/4  0   0 ;
0   0  1/3  0  2/3  0   0 ;
0   0   0   0   0  1/3 2/3;
0   0   0   0   0  1/2 1/2;
0   0   0   0   0  3/4 1/4;
1/2 1/2  0   0   0   0   0 ;
1/4 3/4  0   0   0   0   0 ];
mc = dtmc(P);
openExample('econ/CreateMarkovChainFromRandomTransitionMatrixExample')
%-- 2019/9/20 14:35 --%
[1/3 1/9 0; 0 1/3 1/9; 1/9 0 0]
a=[1/3 1/9 0; 0 1/3 1/9; 1/9 0 0]
eig(a)
[v d]eig(a)
[v d]=eig(a)
[v d]=eig(a')
v'*v
inv(v)*v
a=[1/3 1/9 0; 0 1/3 1/9; 1/9 0 0].*([4/9 4/9 1/9]'.*ones(1,3))
([4/9 4/9 1/9]'.*ones(1,3))
([4/9 4/9 1/9]'*ones(1,3))
a=[1/3 1/9 0; 0 1/3 1/9; 1/9 0 0].*([4/9 4/9 1/9]'*ones(1,3))
a=[1/3 1/9 0; 0 1/3 1/9; 1/9 0 0].\([4/9 4/9 1/9]'*ones(1,3))
a=[1/3 1/9 0; 0 1/3 1/9; 1/9 0 0]./([4/9 4/9 1/9]'*ones(1,3))
a2=[    0.7500    0.2500         0
0    0.7500    0.2500
0.98         0.01         0.01]
[v d]=eig(a2')
[v d]=eig(a')
%-- 2019/9/24 16:21 --%
textCurrent_V2_Selten12
plot(v(:,2))
textCurrent_V2_Selten12
DD2=round(DD,2)
saveas(strcat('.\fig-try\gid=',num2str(gid),'.png'),'png')
strcat('.\fig-try\gid=',num2str(gid),'.png')
saveas(gcf,strcat('.\fig-try\gid=',num2str(gid),'.png'),'png')
saveas(gcf,strcat('./fig-try/gid=',num2str(gid),'.png'),'png')
saveas(gcf,strcat('fig-try/gid=',num2str(gid),'.png'),'png')
strcat('fig-try/gid=',num2str(gid),'.png')
saveas(gcf,strcat('C:/Users/Think/Desktop/V2b/fig-try/gid=',num2str(gid),'.png'),'png')
close all
textCurrent_V2_Selten12
close all
sum(sum(M'))
(sum(M))
clear all
V2_5005_Exp
textCurrent_V2_Selten12
testv2_25x25
S=[0,0;0,1;0,2;0,3;0,4;1,0;1,1;1,2;1,3;1,4;2,0;2,1;2,2;2,3;2,4;3,0;3,1;3,2;3,3;3,4;4,0;4,1;4,2;4,3;4,4]
textCurrent_V2_Selten12
cloae all
close all
clear all
textCurrent_V2_Selten12
close all
clear all
textCurrent_V2_Selten12
saveas(gcf,strcat('C:/Users/Think/Desktop/V2b/fig-try/WZJ-flux-Complex.png'),'png')
%-- 2019/9/25 21:30 --%
textCurrent_V2_Selten12
sum(sum(T))
textCurrent_V2_Selten12
saveas(gcf,strcat('C:/Users/Think/Desktop/V2b/fig-try/Tid=',num2str(gid),'.png'),'png')
saveas(gcf,strcat('C:/Users/Think/Desktop/V2b/fig-try/Tid-V1.png'),'png')
textCurrent_V2_Selten12
%-- 2019/10/23 15:09 --%
a=[-1 -1 ones(1,40)]
[p h]=ttest(a)
%-- 2019/11/27 20:37 --%
main
%-- 2019/12/22 16:46 --%
createIcosahedron
varargout = drawPolygon3d(nodes)
demoDrawPoint3d
createIcosahedron
figure('color','w')
props.Marker='o';
props.MarkerEdgeColor='r';
props.MarkerFaceColor='g';
props.LineStyle='none';
drawPoint3d(gca,nodes,props,'MarkerSize',5)
line(nodes(1,:),nodes(2,:))
line3(nodes(1,:),nodes(2,:))
line(nodes(1,:),nodes(2,:))
plot3(nodes(1,:),nodes(2,:))
nodes(1,:)
plot3(nodes(1,1),nodes(1,2),nodes(1,3),nodes(2,1),nodes(2,2),nodes(2,3))
drawPoint3d(gca,nodes,props,'MarkerSize',5)
hold on
plot3(nodes(1,1),nodes(1,2),nodes(1,3),nodes(2,1),nodes(2,2),nodes(2,3))
plot3(nodes(1,1),nodes(1,2),nodes(1,3),nodes(2,1),nodes(2,2),nodes(2,3),'--')
line3(nodes(1,1),nodes(1,2),nodes(1,3),nodes(2,1),nodes(2,2),nodes(2,3),'--')
line(nodes(1,1),nodes(1,2),nodes(1,3),nodes(2,1),nodes(2,2),nodes(2,3),'--')
drawLine3d( [nodes(1,1),nodes(1,2),nodes(1,3),nodes(2,1),nodes(2,2),nodes(2,3)],'k')
drawLine3d( [nodes(1,1) nodes(1,2) nodes(1,3) nodes(2,1) nodes(2,2) nodes(2,3)],'k')
quiver3
hold on;plot3([nodes(1,:); nodes(2,:)])
hold on;plot3([nodes(1,:); nodes(2,:)],'k')
hold on;plot3([0 0 0; 0 0 1],'k')
hold on;plot3(nodes(1:3),'k')
A=nodes(1,:);B=nodes(2,:);plot3(A(:,1),A(:,2),A(:,3),B(:,1),B(:,2),B(:,3));
%-- 2019/12/22 18:38 --%
drawPoint3d(gca,nodes,props,'MarkerSize',5)
createIcosahedron
drawPolygon3d(nodes)
figure('color','w')
props.Marker='o';
props.MarkerEdgeColor='r';
props.MarkerFaceColor='g';
props.LineStyle='none';
drawPoint3d(gca,nodes,props,'MarkerSize',5)
A=nodes(1,:);B=nodes(2,:);plot3(A(:,1),A(:,2),A(:,3),B(:,1),B(:,2),B(:,3));
drawPoint3d(gca,nodes,props,'MarkerSize',5)
hold on
A=nodes(1,:);B=nodes(2,:);plot3(A(:,1),A(:,2),A(:,3),B(:,1),B(:,2),B(:,3));
A=nodes(1,:);B=nodes(2,:);plot3(A(:,1),A(:,2),A(:,3),B(:,1),B(:,2),B(:,3),'k');
plot3([A(:,1),B(:,1)],[A(:,2),B(:,2)],[A(:,3),B(:,3)],'k');
for i=1:30; A=nodes(edges(i,1),:);B=nodes(edges(i,2),:);plot3(A(:,1),A(:,2),A(:,3),B(:,1),B(:,2),B(:,3),'k');hold on;end
plot3([A(:,1),B(:,1)],[A(:,2),B(:,2)],[A(:,3),B(:,3)],'r');
drawPoint3d(gca,nodes,props,'MarkerSize',5)
hold on
for i=1:30; A=nodes(edges(i,1),:);B=nodes(edges(i,2),:);plot3(A(:,1),A(:,2),A(:,3),B(:,1),B(:,2),B(:,3),'r');hold on;end
PlatomicSolid
%-- 2019/12/23 16:57 --%
PlatomicSolid
plot3(V)
plot3(V(:,1),V(:,2),V(:,3))
plot3(V(:,1),V(:,2),V(:,3), 'o')
for k=1:60;plot3(V(k,1),V(k,2),V(k,3), 'o');hold on;end
for k=1:60;plot3(V(k,1),V(k,2),V(k,3), 'o');text3(V(k,1),V(k,2),V(k,3), num2str(k));hold on;end
for k=1:60;plot3(V(k,1),V(k,2),V(k,3), 'o');text(V(k,1),V(k,2),V(k,3), num2str(k));hold on;end
PlatomicSolid
%-- 2019/12/24 9:22 --%
PlatomicSolid
[B,V] = bucky;
G = graph(B);
p = plot(G);
axis equal
PlatomicSolid
[B,V] = bucky;
G = graph(B);
p = plot(G);
axis equal
PlatomicSolid
sqrt(5)-1
(sqrt(5)+1)/2
PlatomicSolid
%-- 2019/12/24 19:31 --%
PlatomicSolid
(sqrt(5)+1)/2
PlatomicSolid
plot(V(:,end))
PlatomicSolid
plot(V(:,end))
V(:,4)=sqrt( V(:,1).*V(:,1)  )
V(:,4)=sqrt( V(:,1).*V(:,1) +   V(:,2).*V(:,2) )
V(:,5)= V(:,1).*V(:,4)
V(:,5)= V(:,1)./V(:,4)
V(:,6)= V(:,2)./V(:,4)
V(:,7)= V(:,3)
plot(V(:,5:6),'DisplayName','V(:,5:6)')
scatter(getcolumn(V(:,5:6),1),getcolumn(V(:,5:6),2))
V(:,8) = atan(V(:,5)./V(:,6))
plot(V(:,end))
V(:,8) = atan(V(:,6)./V(:,5))
plot(V(:,end))
%-- 2020/1/31 19:18 --%
anan20200131
plot(a1)
anan20200131
plot(acc)
anan20200131
plot(acc)
anan20200131
%-- 2020/6/8 1:44 --%
createTetrahedron4RPSD
axis off
createTetrahedron4RPSD
[xState e123]=createTetrahedron4RPSD(4)
[1,1,0],[-4,4,0], [2,-2,0]
[Omega3XYZ180 ...
AngleMonentum3XYZ ...
Raidus ...
Raidus_Sign ...
] ...
= ...
AngleFunction3DThetaR( [1,1,0],[-4,4,0], [2,-2,0] )
AngleFunction3DThetaR( [1,1,0],[-4,4,0], [2,-2,0] )
[Omega3XYZ180 ...
AngleMonentum3XYZ ...
Raidus ...
Raidus_Sign ...
] ...
= ...
AngleFunction3DThetaR( [1,1,0],[-4,4,0], [2,-2,0] )
%-- 2020/6/8 22:53 --%
see5strategyPatternYQM
%-- 2020/6/8 23:31 --%
see5strategyPatternYQM
close all
see5strategyPatternYQM
%-- 2020/6/9 0:58 --%
see5strategyPatternYQM
plot(a,'DisplayName','a')
see5strategyPatternYQM_a
AngleFunction3DThetaR
see5strategyPatternYQM_a
Angu=Angu+AngleMonentum3XYZ(3);Raid=Raid+Raidus;Raid_sign=Raid_sign+Raidus_Sign(3);
see5strategyPatternYQM_a
xlabel=num2str(Angu)
text(0,0,num2str(Angu))
see5strategyPatternYQM_a
%-- 2020/6/9 15:35 --%
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i]
B=[]; for k = 0:0.2:pi/3; plot(A(1)*exp(i*k),A(2)*exp(i*k),'b');end
B=[]; for k = 0:0.2:pi/3; plot(A(1)*exp(i*k),A(2)*exp(i*k),'b');hold on;end
B=[]; for k = 0:0.2:pi/3; plot(real(A(1)*exp(i*k)),real(A(2)*exp(i*k)),'b');hold on;end
-
B=[]; for k = 0:0.2:pi/3; scatter(real(A(1)*exp(i*k)),real(A(2)*exp(i*k)),'b-');hold on;end
B=[]; for k = 0:0.2:pi/3; scatter(real(A(1)*exp(i*k)),real(A(2)*exp(i*k)),'*');hold on;end
B=[]; for k = 0:0.2:2*pi; scatter(real(A(1)*exp(i*k)),real(A(2)*exp(i*k)),'*');hold on;end
B=[]; for k = 0:0.2:2*pi; scatter(real(A(1)*exp(i*k)),real(A(2)*exp(i*k)),'*');hold on;B=[B; real(A(1)*exp(i*k)),real(A(2)*exp(i*k))]; end
B=[]; for k = 0:0.2:2*pi; scatter(real(A(2)*exp(i*k)),real(A(4)*exp(i*k)),'*');hold on;B=[B; real(A(1)*exp(i*k)),real(A(2)*exp(i*k))]; end
B=[]; for k = 0:0.2:2*pi; scatter(real(A(2)*exp(i*k)),real(A(4)*exp(i*k)),'*');hold on;B=[B; real(A(2)*exp(i*k)),real(A(4)*exp(i*k))]; end
see5strategyPatternYQM_a
AngleFunction3DThetaR
see5strategyPatternYQM_a
aaa=mean(a)
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
for I=2:5
for j=1:I
if I~=j;
B=[]; for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B=[B; cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k))], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd))]) ];
end
end
end
end
B
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
for I=2:5
for j=1:I
if I~=j;
B=[]; for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B=[B; cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0]) ];
end
end
end
end
B
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
B=[]
for I=2:5
for j=1:I
if I~=j;
B1=0;
for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0])  ;
end
B=[B;B1]
end
end
end
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
B=[]
for I=2:5
for j=1:I
if I~=j;
B1=0;
for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0])  ;
end
B=[B;j I B1]
end
end
end
B
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
B=[]
for I=2:5
for j=1:I
if I~=j;
B1=0;
for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0]-[real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0])  ;
end
B=[B;j I B1]
end
end
end
B
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
B=[]
i=sqrt(-1);
for I=2:5
for j=1:I
if I~=j;
B1=0;
for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0]-[real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0])  ;
end
B=[B;j I B1]
end
end
end
B
[B AngRaidwSign]
zz=[B AngRaidwSign]
scatter(getcolumn(zz(:,[5,8]),1),getcolumn(zz(:,[5,8]),2))
scatter(getcolumn(zz(:,[8,end]),1),getcolumn(zz(:,[8,end]),2))
scatter(getcolumn(zz(:,[5,end]),1),getcolumn(zz(:,[5,end]),2))
scatter(getcolumn(zz(:,[5,8]),1),getcolumn(zz(:,[5,8]),2))
zz(:,5).*zz(:,8)
zz=[zz zz(:,5).*zz(:,8)]
B=[]
i=sqrt(-1);
for I=2:5
for j=1:I
if I~=j;
B1=0;
for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k) A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd) A(I)*exp(i*kd)),0]-[real(A(j)*exp(i*k) A(I)*exp(i*k)),0])  ;
end
B=[B;j I B1]
end
end
end
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
B=[]
i=sqrt(-1);
for I=2:5
for j=1:I
if I~=j;
B1=0;
for k = 0:0.2:2*pi; kd=k+0.2;
%                         scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0]-[real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0])  ;
end
B=[B;j I B1]
end
end
end
zz=[B AngRaidwSign]
zz=[zz zz(:,5).*zz(:,8)]
see5strategyPatternYQM_a
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
B=[];kf=0;
i=sqrt(-1);
for I=2:5
for j=1:I
if I~=j;
B1=0;
for k = 0:0.2:2*pi; kd=k+0.2;
kf=kf+1;subplot(12,kf) ;                        scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0]-[real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0])  ;
end
B=[B;j I B1]
end
end
end
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
B=[];kf=0;
i=sqrt(-1);
for I=2:5
for j=1:I
if I~=j;
B1=0;
kf=kf+1;subplot(3,4,kf) ;
for k = 0:0.2:2*pi; kd=k+0.2;
scatter(real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0]-[real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0])  ;
end
B=[B;j I B1]
end
end
end
see5strategyPatternYQM_a
A
A'
conj(A)
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
A=conj(A);
B=[];kf=0;
i=sqrt(-1);
for I=2:5
for j=1:I
if I~=j;
B1=0;
kf=kf+1;subplot(3,4,kf) ;
for k = 0:0.2:2*pi; kd=k+0.2;
scatter(real(A(j)*exp(i*k))*0.1+NE(j),real(A(I)*exp(i*k))*0.1+NE(I),'*');hold on;
B1=B1+ cross([real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0], ...
[real(A(j)*exp(i*kd)),real(A(I)*exp(i*kd)),0]-[real(A(j)*exp(i*k)),real(A(I)*exp(i*k)),0])  ;
end
B=[B;j I B1]
end
end
end
see5strategyPatternYQM_a
zz=[B AngRaidwSign]
zz=[zz zz(:,5).*zz(:,8)]
see5strategyPatternYQM_a
zz=[B AngRaidwSign]
zz=[zz zz(:,5).*zz(:,8)]
B=[];kf=0;
i=sqrt(-1);
for I=2:5
for j=1:I
if I~=j;
B1=0;
kf=kf+1;subplot(3,4,kf) ;
for k = 0:0.2:2*pi; kd=k+0.2;
scatter(real(A(j)*exp(ilambda3*k))*0.5+NE(j),real(A(I)*exp(ilambda3*k))*0.5+NE(I),'*');hold on;
B1=B1+ cross([real(A(j)*exp(ilambda3*k)),real(A(I)*exp(ilambda3*k)),0], ...
[real(A(j)*exp(ilambda3*kd)),real(A(I)*exp(ilambda3*kd)),0]-[real(A(j)*exp(ilambda3*k)),real(A(I)*exp(ilambda3*k)),0])  ;
end
B=[B;j I B1]
end
end
end
see5strategyPatternYQM_a
%-- 2020/6/10 18:47 --%
A=[0.3147+0.2844i
0.3939-0.3464i
-0.6364
-0.1077+0.2811i
0.0355-0.2191i];
plot(A(1))
compass(A)
feather(A)
compass(A)
see5strategyPatternYQM_a
AngleFunction3DThetaR
see5strategyPatternYQM_a
%-- 2020/6/10 20:42 --%
see5strategyPatternYQM_a
dataset=load('Binmore8.csv');
stdr=std(dataset);                      %求个变量的标准差
[n,m]=size(dataset);                    %定义矩阵行列数
sddata=dataset./stdr(ones(n,1),:);      %将原始数据采集标准化
sddata                                  %输出标准化数据
[p,princ,eigenvalue,t2]=princomp(sddata);%调用前三个主成分系数
p3=p(:,1:3);                            %提取前三个主成分得分系数,通过看行可以看
dataset=load('Binmore8.csv');
dataset=[:,[1:3 5:7]];
dataset=load('Binmore8.csv');
dataset=[:,[1 2 3 5 6 7]]
dataset=load('Binmore8.csv');
dataset=dataset(:,[1 2 3 5 6 7])
dataset=load('Binmore8.csv');
dataset=dataset(:,[1 2 3 5 6 7]);
stdr=std(dataset);                      %求个变量的标准差
[n,m]=size(dataset);                    %定义矩阵行列数
sddata=dataset./stdr(ones(n,1),:);      %将原始数据采集标准化
sddata                                  %输出标准化数据
[p,princ,eigenvalue,t2]=princomp(sddata);%调用前三个主成分系数
p3=p(:,1:3);                            %提取前三个主成分得分系数,通过看行可以看
sc=princ(:,1:3)
e=eigenvalue(1:3)
sc                                      %输出前三个主成分得分值
e=eigenvalue(1:3)';                     %提取前三个特征根并转置
M=e(ones(m,1),:).^0.5;                  %输出前三个特征根并转置
compmat=p3.*M;                          %利用特征根构造变换矩阵
per=100*eigenvalue/sum(eigenvalue);     %求出成分载荷矩阵的前三列
per
%求出各主成分的贡献率
cumsum(per);                            %列出各主成分的累积贡献率
figure(1)
pareto(per);                            %将贡献率绘成直方图
t2
figure(2)
%输出各省与平局距离
plot(eigenvalue,'r+');                  %绘制方差贡献散点图
hold on
%保持图形
plot(eigenvalue,'g-');                  %绘制方差贡献山麓图
figure(3)
%关闭图形
plot(princ(:,1),princ(:,2),'+');        %绘制2维成份散点图
%gname
%,(rowname)                         %标示个别散点代表的省data市
[st2,index]=sort(t2);
%st2=flipud(st2);
%index=flipud(index);
%extreme=index(1);
%-- 2020/6/14 1:20 --%
see5strategyPatternYQM_a
contour(a(:,1:3))
contour(a(:,1),a(:,2))
hist2(a(:,1),a(:,2))
hist(a(:,1),a(:,2))
histc(a(:,1),a(:,2))
histc(a(:,1))
histc(a(:,1),100)
histc(a(:,1),0:1)
histc(a(:,1),0:.1:1)
histc
histogram2(a(:,1),a(:,2))
histogram2(a(:,1),a(:,2),'DisplayStyle','tile','ShowEmptyBins','on')
histcounts2(a(:,1),a(:,2),'DisplayStyle','tile','ShowEmptyBins','on')
histcounts2(a(:,1),a(:,2))
histogram2(a(:,1),a(:,2),'DisplayStyle','tile','ShowEmptyBins','on')
histogram2(a(:,4),a(:,5),'DisplayStyle','tile','ShowEmptyBins','on')
for iI=2:5
for j=1:iI
if iI~=j;
k=k+1;
subplot(3,4,k) ;histogram2(a(:,j),a(:,I),'DisplayStyle','tile','ShowEmptyBins','on')
hold on;             title(strcat('Mean = [',num2str(mean(0)),'] Obs=',num2str(L),'  tick/s=',num2str([j iI])));
axis square;xlim([-0.1 0.5]);ylim([-0.1 0.5]);box on; end;end;end
k=0
for iI=2:5
for j=1:iI
if iI~=j;
k=k+1;
subplot(3,4,k) ;histogram2(a(:,j),a(:,I),'DisplayStyle','tile','ShowEmptyBins','on')
hold on;             title(strcat('Mean = [',num2str(mean(0)),'] Obs=',num2str(L),'  tick/s=',num2str([j iI])));
axis square;xlim([-0.1 0.5]);ylim([-0.1 0.5]);box on; end;end;end
for iI=2:5
for j=1:iI
if iI~=j;
k=k+1;
subplot(3,4,k) ;histogram2(a(:,j),a(:,iI),'DisplayStyle','tile','ShowEmptyBins','on')
hold on;             title(strcat('Mean = [',num2str(mean(0)),'] Obs=',num2str(L),'  tick/s=',num2str([j iI])));
axis square;xlim([-0.1 0.5]);ylim([-0.1 0.5]);box on; end;end;end
L=14236
for iI=2:5
for j=1:iI
if iI~=j;
k=k+1;
subplot(3,4,k) ;histogram2(a(:,j),a(:,iI),'DisplayStyle','tile','ShowEmptyBins','on')
hold on;             title(strcat('Mean = [',num2str(mean(0)),'] Obs=',num2str(L),'  tick/s=',num2str([j iI])));
axis square;xlim([-0.1 0.5]);ylim([-0.1 0.5]);box on; end;end;end
k=0
for iI=2:5
for j=1:iI
if iI~=j;
k=k+1;
subplot(3,4,k) ;histogram2(a(:,j),a(:,iI),'DisplayStyle','tile','ShowEmptyBins','on')
hold on;             title(strcat('Mean = [',num2str(mean(0)),'] Obs=',num2str(L),'  tick/s=',num2str([j iI])));
axis square;xlim([-0.1 0.5]);ylim([-0.1 0.5]);box on; end;end;end
k=0
for iI=2:5
for j=1:iI
if iI~=j;
k=k+1;
subplot(3,4,k) ;histogram2(a(:,j),a(:,iI),'DisplayStyle','tile','ShowEmptyBins','on')
hold on;             title(strcat('Mean = [',num2str(mean(0)),'] Obs=',num2str(L),'  tick/s=',num2str([j iI])));
axis square;box on; end;end;end
see5strategyPatternYQM_a
k=0
for iI=2:5
for j=1:iI
if iI~=j;
k=k+1;
subplot(3,4,k) ;histogram2(a(:,j),a(:,iI),'DisplayStyle','tile','ShowEmptyBins','on')
hold on;             title(strcat('Mean = [',num2str(mean(0)),'] Obs=',num2str(L),'  tick/s=',num2str([j iI])));
axis square;box on; end;end;end
%-- 2020/6/14 23:59 --%
vol3d
plot(map,'DisplayName','map')
squeeze(D)
load spine
figure
image(X)
colormap(map)
%Make the random volume
mat = rand(50,50,100);
%Place high values in particular parts of the volume
sigCoors.rows = [23:33];
sigCoors.columns = [40:45];
sigCoors.time = [55:85];
mat(sigCoors.rows, sigCoors.columns, sigCoors.time) = 10.*rand(length(sigCoors.rows),   length(sigCoors.columns), length(sigCoors.time));
%Visualize the volume:
plot3
VaryMarkerSizeExample
scatter3(x,y,z,s)
%-- 2020/6/15 0:59 --%
a = rand(1000,3);       % Create random matrix, use your data here
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.2;              % Tolerance for (squared) distance to count as "nearby"
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance
n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;
a_0=load('yqm-20200607.csv');
%     a_0=load('mp8（从界面历史策略整理而来的数据）2.csv');
a_0=[a_0(:,1)-a_0(:,2) a_0(:,2)-a_0(:,3) a_0(:,3)-a_0(:,4) a_0(:,4)-a_0(:,5) a_0(:,5)];
L_0=length(a_0(:,1));
see5strategyPatternYQM_a
a_0=load('yqm-20200607.csv');
%     a_0=load('mp8（从界面历史策略整理而来的数据）2.csv');
a_0=[a_0(:,1)-a_0(:,2) a_0(:,2)-a_0(:,3) a_0(:,3)-a_0(:,4) a_0(:,4)-a_0(:,5) a_0(:,5)];
L_0=length(a_0(:,1));
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.2;              % Tolerance for (squared) distance to count as "nearby"
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance
n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;
a=a_0(:,3:5)
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.2;              % Tolerance for (squared) distance to count as "nearby"
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance
n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.01;              % Tolerance for (squared) distance to count as "nearby"
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance
n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.005;              % Tolerance for (squared) distance to count as "nearby"
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance
n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.001;              % Tolerance for (squared) distance to count as "nearby"
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance
n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;
%-- 2020/6/15 21:05 --%
a_0=load('yqm-20200607.csv');
%     a_0=load('mp8（从界面历史策略整理而来的数据）2.csv');
a_0=[a_0(:,1)-a_0(:,2) a_0(:,2)-a_0(:,3) a_0(:,3)-a_0(:,4) a_0(:,4)-a_0(:,5) a_0(:,5)];
L_0=length(a_0(:,1));
tick_second = 1;  %+round(7*tick);
a = a_0(1:tick_second:L_0,:);
L=length(a(:,1));
NE=[ 444/2987,22/103, 641/2987, 288/2987, 976/2987];
see5strategyPatternYQM_a
a=a_0(:,3:5)
n = zeros(size(a,1),1); % Set up array for number of nearby points
tol = 0.001;              % Tolerance for (squared) distance to count as "nearby"
sz = size(a,1);         % Shorthand for size of data
% Loop over every point
for ii = 1:sz;
dists = sum((repmat(a(ii,:), sz, 1) - a).^2, 2); % Get standard Euclidean distance
n(ii) = nnz(dists < tol); % Count number of points within tolerance
end
% Plot, colouring by an nx3 RGB array, in this case just
% scaling the red and having no green or blue.
scatter3(a(:,1), a(:,2), a(:,3), [], [n./max(n), zeros(numel(n),2)], 'filled');
grid on;
view(30,93)
view(130,93)
view(130,33)
%-- 2020/6/17 1:34 --%
anan20200131
matrix_lower_tri_to_surf
a_0=load('yqm-20200607.csv');
see5strategyPatternYQM_a
AngleFunction3DThetaR
see5strategyPatternYQM_a
A=[sum(a(:,1:3)')' a(:,4) a(:,5)];
matrix_lower_tri_to_surf(A)
plot(A,'DisplayName','A')
A=[sum(a(:,1:3)')' a(:,4) a(:,5)];
contour(A)
surf(A)
histogram(A)
hist3(A)
histogram3(A)
histogram(A)
histc(A)
contour(A(:,1:2))
pie(A(:,1:2))
scatter(getcolumn(A(:,1:2),1),getcolumn(A(:,1:2),2))
scatter(getcolumn(A(:,1:2),1),getcolumn(A(:,1:2),2))
histogram(A(:,1:2))
bar(A(:,2:end),'DisplayName','A(:,2:end)')
contour(A(:,2:end))
mesh(A(:,2:end))
histogram(A(:,1:2))
histogram(A(:,2:end))
plotmatrix(A(:,2:end))
plotmatrix(A(:,2:end))
plotmatrix(a)
hist3(A)
hist3(A(:,1:2))
hist3(A(:,1:2),100)
hist3(A(:,1:2))
hist3(A(:,1:2),'CDataMode','auto','FaceColor','interp')
histogram2(A(:,1:2),'CDataMode','auto','FaceColor','interp')
histogram2(A(:,1:2),20,'CDataMode','auto','FaceColor','interp')
histogram2(A(:,1:2),20)
histogram2(A(:,1),A(:,2),20)
histogram2(A(:,1),A(:,2),120)
Xedges = linspace(0,1);
Yedges = linspace(0,1);
binScatterPlot(A(:,1),A(:,2),Xedges,Yedges)
Y = tall(randn(1e5,1));
Y = all(randn(1e5,1));
%-- 2020/7/12 20:16 --%
x=[2,2,1,3];p=[1/8,1/8,2/8,4/8];
ExpectationProb
%-- 2020/7/18 16:09 --%
ZSJ20200718_HalfSpaceProduction
a_0=load('yqm-20200607.csv');
see5strategyPatternYQM_a
%-- 2020/8/1 21:50 --%
Strategy5
payoff_matrix =  [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0];
Strategy5
[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
BJ1983
xdot1 + ydot1 + zdot1
simple(xdot1 + ydot1 + zdot1)
simplex(xdot1 + ydot1 + zdot1)
simplify(xdot1 + ydot1 + zdot1)
BJ1983
(eval([RPpartialXY(1) RPpartialXY(2) ; RPpartialXY(3) RPpartialXY(4)]))
factor(eval([RPpartialXY(1) RPpartialXY(2) ; RPpartialXY(3) RPpartialXY(4)]))
JacobianRP2x2=  (eval([RPpartialXY(1) RPpartialXY(2) ; RPpartialXY(3) RPpartialXY(4)])) ;
LatexJacobianRP2x2=strcat( '$J_{|_{x_0,y_0}}',latex(JacobianRP2x2) ,'$' );
LatexJacobianRP2x2
BJ1983
eval(payoffM)
syms p1 p2 p3 p4 y1
y1o2 = y1/2;
Mp=[p1 0 0 0;0 p2 0 0; 0 0 p3 0; 0 0 0 p4];
Mq=[1/2 1-y1o2 1-y1o2 1-y1o2; y1o2 1/2  1-y1o2  1-y1o2; y1o2  y1o2 1/2  1-y1o2; y1o2  y1o2 1-y1o2 1/2 ];
payoffM = Mp*Mq;
p1=1;p2=2;p3=3;p4=4;y1=0.3; M=eval(payoffM);N=M';  %          [A,B,a,b,iterations,err,ms]=bimat(M,N);
[NashVectorA,NashVectorB,NashEarn_a,NashEarn_b,iterations,err,ms]=bimat(M,N);
[NashVectorA 999 NashEarn_a;NashVectorB 999 NashEarn_b]
[NashVectorA,NashVectorB,NashEarn_a,NashEarn_b,iterations,err,ms]=bimat(M,N)
%-- 2020/8/2 11:03 --%
a=1:8
b=1:8
c=a.*10 + b
for a=1:8
for b=1:8
c=a*10 + b;end;end
c=[];for a=1:8;for b=1:8
c=[c;a*10 + b];end;end
c=[];r=1;k=0
for a=1:8;for b=1:8
c=[c;a*10 + b];
plot(r*cos(2*pi/64*k),r*sin(2*pi/64*k),'ro')
text(r*cos(2*pi/64*k),r*sin(2*pi/64*k),num2str(c))
k=k+1;
end;end
anan8x8cyclicstate
%-- 2020/8/3 11:54 --%
m = 8;
n = 6;
A = randi(5,[m n]);
imagesc(A);
m = 12;
n = 8;
A = randi(5,[m n]);
Arot = flipud(A);
Arot = [ Arot; Arot(end,:) ];
Arot = [ Arot, Arot(:,end) ];
pcolor(Arot);figure(gcf);
m = 12;
n = 8;
A = randi(5,[m n]);
Arot = flipud(A);
pcolor(Arot);figure(gcf);
m = 64;
n = 64;
A = randi(5,[m n]);
Arot = flipud(A);
pcolor(Arot);figure(gcf);
axis square
Get8x8Data202008
p=find(stateList==RPara(i,10));
M=zeros(64);
for i=1:12000-1;
p=find(stateList==RPara(i,10));
q=find(stateList==RPara(i,11));
M(p,q) = M(p,q)+1;
end;M
Arot = flipud(M);
pcolor(Arot);figure(gcf);
Get8x8Data202008
colormap(parula(5))
colormap(parula(2))
colormap(parula(12))
colormap(parula(120))
colorbar('Ticks',[0,1000,2000,3000],)
colorbar('Ticks',[0,1000,2000,3000])
Get8x8Data202008
%-- 2020/8/3 17:03 --%
bimat
readc
b=csvread(filename,7,19)
b=fread(filename,7,19)
b=fopen(filename,7,19)
b=fopen(filename)
fid = fopen(filename);
tline = fgetl(fid);
while ischar(tline)
disp(tline)
tline = fgetl(fid);
end
fclose(fid);
readc
strcat('E:/AbedData/01/',oldname)
readc
strcat('E:/AbedData/01/',oldname)
strcat('E:/AbedData/01/',newname)
readc
strcat('E:/AbedData/01/',oldname)
strcat('E:/AbedData/01/',newname)
readc
tline
command = ['xcopy' 64 strcat('E:/AbedData/01/',oldname) 64 strcat('E:/AbedData/01/',newname)];
status = dos(command);
if status == 0
disp([oldname, ' 已被重命名为 ', newname])
else
disp([oldname, ' 重命名失败!'])
end
command = ['xcopy'  strcat('E:/AbedData/01/',oldname) strcat('E:/AbedData/01/',newname)];
status = dos(command);
if status == 0
disp([oldname, ' 已被重命名为 ', newname])
else
disp([oldname, ' 重命名失败!'])
end
readc
cell2str(status)
cell2str(command)
cellstr(command)
char(command)
status = dos(char(command));
char(command)
status = dos(char(command'));
eval(command)
command = ['!rename' 32 strcat('E:/AbedData/01/',oldname) 32 strcat('E:/AbedData/01/', newname) ];
eval(command)
strcat('E:/AbedData/01/',oldname)
char(strcat('E:/AbedData/01/',oldname))
cell(strcat('E:/AbedData/01/',oldname))
strcell(strcat('E:/AbedData/01/',oldname))
strel(strcat('E:/AbedData/01/',oldname))
fileName=strcat('E:/AbedData/01/',oldname);fileName2 = strcat('E:/AbedData/01/',newname);  eval(['!rename' 32 fileName 32 fileName2]);
fileName=strcat('E:/AbedData/01/',oldname);fileName2 = strcat('E:/AbedData/01/',newname);  exec(['!rename' 32 fileName 32 fileName2]);
fileName=strcat('E:/AbedData/01/',oldname);fileName2 = strcat('E:/AbedData/01/',newname);  dos(['!rename' 32 fileName 32 fileName2]);
fileName=strcat('E:/AbedData/01/',oldname);fileName2 = strcat('E:/AbedData/01/',newname);  dos(['rename' 32 fileName 32 fileName2]);
readc
pwd
path('newpath')
clear
%-- 2020/8/4 4:01 --%
readc
userpath( 'E:/AbedData/01/')
files = dir('*.csv')
readc
pathtool
readc
pwd
readc
eval(['!rename' 32 '8x8_Y_s1xxx0.3645220705315265.csv' 32 '8x8_Y_s1xxx0.3645220705315265.ccc']);
fid = fopen(filename);
tline = fgetl(fid);k=1
while ischar(tline)
disp(tline)
if strfind(tline,'[8 -8]')
%             copy(filename,strcat(filename,'8'))
newname = strcat('8',num2str(k),'.aaa');
k=k+1;
end
tline = fgetl(fid);
end
files = dir('*.csv')
len=length(files);
readc
files(i).name
fid = fopen(files(i).name)
fid = fopen(files(i).name);
tline = fgetl(fid);k=1
while ischar(tline)
disp(tline)
if strfind(tline,'[8 -8]')
%             copy(filename,strcat(filename,'8'))
newname = strcat('8',num2str(k),'.aaa');
k=k+1;
end
tline = fgetl(fid);
end
fclose(fid);
eval(['!rename' 32 oldname  32 newname])
readc
%-- 2020/8/4 13:10 --%
readc
strfind(tline,'[6 -6] [0 0]]]')
strfind(tline,'[8 -8] [0 0]]]')
strfind(tline,'[10 -10] [0 0]]]')
strfind(tline,'[6 -6] [0 0]]]')
strfind(tline,'[8 -8] [0 0]]]')
strfind(tline,'[8 -8] [0 0]] ]')
readc
files(i).name
char(files(i).name,1)
char(files(i).name,2)
char(files(i).name,[1:2])
readc
[num,txt,raw] = xlsread('phase_exp Strategy distributions_3.csv');
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps = [v.u1 v.u2-v.u6 v.u6-v.u10 v.u10-v.u14 v.u14];
x = ps(:,1);
[num,txt,raw] = xlsread('X_10_35.aaa');
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps = [v.u1 v.u2-v.u6 v.u6-v.u10 v.u10-v.u14 v.u14];
x = ps(:,1);
readc
[num,txt,raw] = xlsread('X_10_35.xls');
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps = [v.u1 v.u2-v.u6 v.u6-v.u10 v.u10-v.u14 v.u14];
x = ps(:,1);
readc
clear
%-- 2020/8/4 14:57 --%
readc
oldnamf = replace(oldname,'X_','Y_')
oldnamf = strrep(oldname,'X_','Y_')
readc
[num,txt,raw] = xlsread('X_10_35.xls');
[num,txt,raw] = load('X_10_35.xls');
readc
[num,txt,raw] = load('X_10_35.xls');
[num,txt,raw] = xlsread('X_10_135.xls');
readc
[num,txt,raw] = xlsread('X_10_35.csv')
[num,txt,raw] = xlsread('X_10_135.csv')
readc
close all
%-- 2020/8/4 15:18 --%
readc
[num,txt,raw] = xlsread('X_10_135.csv')
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps = [v.u1 v.u2-v.u6 v.u6-v.u10 v.u10-v.u14 v.u14];
ps = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22 v.u26 v.u30]
plot(ps(:,2:end),'DisplayName','ps(:,2:end)')
readc
save(PS)
save('PS.csv',PS)
save('PS.csv','PS')
save('PS')
load('PS2.mat')
clear
readc
load(PS22)
load(PS22.mat)
PS22.mat
(PS22.mat)
load('PS22.mat')
max(PS(3,4:11))
[M,I] =max(PS(3,4:11))
[M,I] =max(PS(16,4:11))
max
PS(16,4:11)
[M,I] =max(PS(16,4:11))
find(PS(16,4:11)==M)
Mf=find(PS(16,4:11)==M); randi(len(Mf))
Mf=find(PS(16,4:11)==M); randi(length(Mf))
[M,I] =max(PS(16,4:11)); Mf=find(PS(16,4:11)==M); tmp=randi(length(Mf)); Mf(tmp)
[M,I] =max(PS(16,4:11)); Mf=find(PS(16,4:11)==M); tmp=randi(length(Mf)); chooseStrategy=Mf(tmp)
Markov=[];for m=1:length(PS(:,1)); [M,I] =max(PS(m,4:11)); Mf=find(PS(m,4:11)==M); tmp=randi(length(Mf)); chooseStrategyX=Mf(tmp);
[M,I] =max(PS(m,13:20)); Mf=find(PS(m,13:20)==M); tmp=randi(length(Mf)); chooseStrategyY=Mf(tmp);
Markov=[Markov; chooseStrategyX chooseStrategyY];end
%-- 2020/8/4 19:22 --%
load('PS22.mat')
Markov=[Markov; chooseStrategyX chooseStrategyY];end
Markov=[];for m=1:length(PS(:,1)); [M,I] =max(PS(m,4:11)); Mf=find(PS(m,4:11)==M); tmp=randi(length(Mf)); chooseStrategyX=Mf(tmp);
[M,I] =max(PS(m,13:20)); Mf=find(PS(m,13:20)==M); tmp=randi(length(Mf)); chooseStrategyY=Mf(tmp);
Markov=[Markov; chooseStrategyX chooseStrategyY];end
%-- 2020/8/4 19:22 --%
T=rowsorts(PS, [1 2 3])
a
T=rawsorts(PS, [1 2 3])
T=sortraws(PS, [1 2 3])
T=sortrows(PS, [1 2 3])
T6=PS(find(PS(:,1)==6),:)
T8=PS(find(PS(:,1)==8),:)
T10=PS(find(PS(:,1)==10),:)
plot(T6(:,4:11),'DisplayName','T6(:,4:11)')
T6poll=zeros(501,16); for i=1:501;  a=mean(i:501:6012,[4:11 13 20]);end
T6poll=zeros(501,16); for i=1:501;  a=mean(i:501:6012,[4:11 13:20]);T6poll(i,:)=a;end
T6poll=zeros(501,16); for i=1:501;  a=T6(i:501:6012,[4:11 13:20]);T6poll(i,:)=mean(a);end
plot(T6poll(:,1:8),'DisplayName','T6poll(:,1:8)')
plot(T6poll(:,9:end),'DisplayName','T6poll(:,9:end)')
a=T8(i:501:6012,[4:11 13:20])
eeee
clear
eeee
figure(6)
plot(T6poll(:,1:8),'DisplayName','T6poll(:,1:8)')
plot(T6poll(:,1+8:8+8),'DisplayName','T6poll(:,1:8)')
figure(6)
plot(T6poll(:,1:8) )
plot(T6poll(:,1+8:8+8) )
figure(8)
plot(T8poll(:,1:8) )
plot(T8poll(:,1+8:8+8) )
figure(10)
plot(T10poll(:,1:8) )
plot(T10poll(:,1+8:8+8) )
eeee
plot(T6accu(:,1:8),'DisplayName','T6accu(:,1:8)')
eeee
plot(T6accu(:,9:end),'DisplayName','T6accu(:,9:end)')
eeee
readc
eeee
readc
eeee
readc
eeee
readc
eeee
figure(5)
subplot(1,2,1);plot(T6accu(1:600,1:8) ); axis square;ylim([0 250])
subplot(1,2,2);plot(T6accu(1:600,1+8:8+8) ); axis square;ylim([0 250])
figure(7)
subplot(1,2,1);plot(T6accu(1:600,1:8) ); axis square;ylim([0 250])
subplot(1,2,2);plot(T6accu(1:600,1+8:8+8) ); axis square;ylim([0 250])
eeee
%-- 2020/8/5 18:54 --%
eigenVV8x8
S=solve(V_Eq_0)
format short
S=vpasolve(V_Eq_0, x1, [0, 1]);
S
S=vpasolve(V_Eq_0, x2, [0, 1]);
S
S=fsolve(V_Eq_0);
%-- 2020/8/27 18:47 --%
eeee
%-- 2020/8/27 19:34 --%
eeee
CTime=zeros(16);for i=1:16; for j=1:16; tmp=T6accu(:,i)-T6accu(:,j);[va id]=min(abs(tmp));CTime(i,j)=id;end;end
eeee
save('CrossoveTime26.mat')
Ctime
A=CTime(:,:,1)
Ctime
A=CTime(:,:,1)
A2=CTime(:,:,2)
A1=CTime(:,:,1)
Ctime
%-- 2020/8/28 8:38 --%
Impulsive
[h.p]=ttest(a(:,4),a(:.6))
[h.p]=ttest(a(:,4),a(:,6))
[h,p]=ttest(a(:,4),a(:,6))
Impulsive
Sih = Sig(find(Sig(:,2)==4 | Sig(:,2)==6 ),:)
Sih = Sig(find((Sig(:,2)==2 | Sig(:,2)==6 )&(Sig(:,1)!=2 | Sig(:,1)!=6 )) ,:)
Sih = Sig(find((Sig(:,2)==2 | Sig(:,2)==6 )&(Sig(:,1)~=2 | Sig(:,1)~=6 )) ,:)
Sih = Sig(find((Sig(:,2)==2 | Sig(:,2)==6 )&(Sig(:,1)~=2 & Sig(:,1)~=6 )) ,:)
A = sortrows(Sig,1)
A = sortrows(Sih,1)
Sih_Y = Sig(find((Sig(:,2)==2+8 | Sig(:,2)==4+8)&(Sig(:,1)~=2+8 & Sig(:,1)~=4+8 )) ,:)
Sih_X = Sig(find((Sig(:,2)==2 | Sig(:,2)==6)&(Sig(:,1)~=2 & Sig(:,1)~=6 )) ,:)
Sih_Y = Sig(find((Sig(:,2)==2+8 | Sig(:,2)==4+8)&(Sig(:,1)~=2+8 & Sig(:,1)~=4+8 )) ,:)
Impulsive
Si8_X = Sig(find((Sig(:,2)==1)&(Sig(:,1)~=1)) ,:)
% r = [207,:]
Si8_Y = Sig(find((Sig(:,2)==4+8)&(Sig(:,1)~=4+8)) ,:)
Impulsive
%-- 2020/8/28 14:13 --%
Impulsive
save('ImpulseSignal_Stat')
unique(Si10_X(:,1))
unique(Si8_X(:,1))
unique(Si6_X(:,1))
ImpulseSignal=[661;unique(Si6_X(:,1));662;;881;unique(Si8_X(:,1));882;;101;unique(Si10_X(:,1));102;]
unique(Si6_X(:,2))
unique(Si6_Y(:,1))
ImpulseSignal=[661;unique(Si6_X(:,1));662;unique(Si6_Y(:,1));881;unique(Si8_X(:,1));882;;101;unique(Si10_X(:,1));102;]
ImpulseSignal=[661;unique(Si6_X(:,1));662;unique(Si6_Y(:,1)); ...
881;unique(Si8_X(:,1));882;unique(Si8_Y(:,1));...
101;unique(Si10_X(:,1));102;;unique(Si10_Y(:,1))];
clear
Ctime
%-- 2020/8/28 21:31 --%
Ctime
clear
Ctime
clear
load('PS26.mat')
Get8x8Data202008
DistributionEvol
R=sorteows(R,[1 2 3])
R=sortrows(R,[1 2 3])
eeee
T3=[];
[n8x8list,stateList,R] = Get8x8Data(Para);
% experimentid	subjectid	periodid	groupid	myNo	myStrategy	opStrategy	myPay	Para	t0	tp1	tm1
for m=1:length(n8x8list(:,1))
if n8x8list(:,1)==1
x8=zeros(1,8); x8(n8x8list(m,6))=1;
x8=zeros(1,8); x8(n8x8list(m,7))=1;
T3 = [T3; n8x8list(m,9) n8x8list(m,2:3) x8 n8x8list(m,3) y8 n8x8list(m,1)];
end
end
T3=[];
[n8x8list,stateList,R] = Get8x8Data(Para);
% experimentid	subjectid	periodid	groupid	myNo	myStrategy	opStrategy	myPay	Para	t0	tp1	tm1
for m=1:length(n8x8list(:,1))
if n8x8list(:,5)==1
x8=zeros(1,8); x8(n8x8list(m,6))=1;
y8=zeros(1,8); y8(n8x8list(m,7))=1;
T3 = [T3; n8x8list(m,9) n8x8list(m,2:3) x8 n8x8list(m,3) y8 n8x8list(m,1)];
end
end
T3
T3=[];
YiJiaDesign = [6 8 10];
ParaID=1;
DisEvolution =[];
for ParaID=1:3
Para = YiJiaDesign(ParaID);
[n8x8list,stateList,R] = Get8x8Data(Para);
% experimentid	subjectid	periodid	groupid	myNo	myStrategy	opStrategy	myPay	Para	t0	tp1	tm1
for m=1:length(n8x8list(:,1))
if n8x8list(:,5)==1
x8=zeros(1,8); x8(n8x8list(m,6))=1;
y8=zeros(1,8); y8(n8x8list(m,7))=1;
T3 = [T3; n8x8list(m,9) n8x8list(m,2:3) x8 n8x8list(m,3) y8 n8x8list(m,1)];
end
end
end
clear
expdata20
PS = T3;
save('PE26.mat','PS')
Impulsive
%-- 2020/8/29 0:20 --%
Impulsive
%-- 2020/8/29 13:15 --%
ImpulExp
clear
ImpulExp
ImpulseSignal
%-- 2020/8/29 22:42 --%
ImpulExp
ImpulseSignal_unique=[661;unique(Si6_X(:,1));662;unique(Si6_Y(:,1)); ...
881;unique(Si8_X(:,1));882;unique(Si8_Y(:,1));...
101;unique(Si10_X(:,1));102;;unique(Si10_Y(:,1))];
ImpulseSignal_all=[661; (Si6_X(:,1));662; (Si6_Y(:,1)); ...
881; (Si8_X(:,1));882; (Si8_Y(:,1));...
101; (Si10_X(:,1));102; (Si10_Y(:,1))];
ImpulseSignal_all
ImpulseSignal_all=[ (Si6_X(:,1)); ; (Si6_Y(:,1)); ...
; (Si8_X(:,1)); ; (Si8_Y(:,1));...
; (Si10_X(:,1)); ; (Si10_Y(:,1))]
ImpulseSignal_all=[ (Si6_X ); ; (Si6_Y ); ...
; (Si8_X(:,1)); ; (Si8_Y );...
; (Si10_X ); ; (Si10_Y )]
Impulsive
DistributionEvol
Ctime
A10 = sortrows(Crossover_si_sj_time_strength,-4) ;
Ctime
plotK20190408
%-- 2020/8/30 12:08 --%
Ctime
A10 = sortrows(Crossover_si_sj_time_strength,-4)
A6(:,6) = A(:,4)./ A(:,3)
A6(:,6) = A6(:,4)./ A6(:,3)
sortrows(A6(:,6))
A7=sortrows(A6(:,6))
A7=sortrows(A6,6)
%-- 2020/8/30 16:54 --%
YQM20200830
A=eval([S.x1 S.x2 S.x3  S.y1 S.y2 S.y3])
t=8/15, 2/15, 1/3, 8/15, 2/15, 1/3
t=[8/15, 2/15, 1/3, 8/15, 2/15, 1/3]
x1=t(1);x2=t(2);x(3)=t(3);y1=t(4);y2=t(5);y3=t(6);
eval(V_Eq_0)
x1=t(1);x2=t(2);x3=t(3);y1=t(4);y2=t(5);y3=t(6);
eval(V_Eq_0)
[S.x1 S.x2 S.x3 S.y1 S.y2 S.y3]
t=[1/5, 2/5, 2/5, 8/15, 2/15, 1/3]
x1=t(1);x2=t(2);x3=t(3);y1=t(4);y2=t(5);y3=t(6);
eval(V_Eq_0)
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2 x3]'
YQM20200830
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
YQM20200830
eigen_value
YQM20200830
A=A(1:9,:)
clear
YQM20200830
A=A(1:9,:)
A
YQM20200830
round(eigen_value,2)
YQM20200830
%-- 2020/8/30 20:21 --%
bimat
M=[2 -1; -1 1]; N=[1 -1; -1 2];[A,B,a,b,iterations,err,ms]=bimat(M,N);
bimat
M=[2 -1; -1 1]; N=[1 -1; -1 2];[A,B,a,b,iterations,err,ms]=bimat(M,N);
M=[2 -1; -1 1]; N=[1 -1; -1 2];[A,B,a,b,iterations,err,ms]=bimat(M,N)
M=[-2 3 -3; -1 -3 0; 3 -1 1]; N=-M;[A,B,a,b,iterations,err,ms]=bimat(M,N)
YQM20200830
[-2 3 -3; -1 -3 0; 3 -1 1]
YQM20200830
A=A(1:15,:)
x1=A(length(A(:,1)),1); x2=A(length(A(:,1)),2);
x3=A(length(A(:,1)),3);
y1=A(length(A(:,1)),4); y2=A(length(A(:,1)),5);
y3=A(length(A(:,1)),6);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
payoff_matrix_1 = [-2 3 -3; -1 -3 0; 0 -1 1]
%-- 2020/8/31 10:00 --%
YQM20200830
Ctime
Si6_X = Sig(find((Sig(:,2)==2 | Sig(:,2)==6)&(Sig(:,1)~=2 & Sig(:,1)~=6 )) ,:)
% r = [207,:]
Si6_Y = Sig(find((Sig(:,2)==2+8 | Sig(:,2)==4+8)&(Sig(:,1)~=2+8 & Sig(:,1)~=4+8 )) ,:)
% r = []
Sig=A6;
Si6_X = Sig(find((Sig(:,2)==2 | Sig(:,2)==6)&(Sig(:,1)~=2 & Sig(:,1)~=6 )) ,:)
Si6_Y = Sig(find((Sig(:,2)==2+8 | Sig(:,2)==4+8)&(Sig(:,1)~=2+8 & Sig(:,1)~=4+8 )) ,:)
Sig=A8;
Si8_X = Sig(find((Sig(:,2)==1)&(Sig(:,1)~=1)) ,:)
Si8_Y = Sig(find((Sig(:,2)==4+8)&(Sig(:,1)~=4+8)) ,:)
Sig=A10;
Si10_X = Sig(find((Sig(:,2)==1 | Sig(:,2)==2 | Sig(:,2)==5)&(Sig(:,1)~=1 & Sig(:,1)~=2 & Sig(:,1)~=5)) ,:)
Si10_Y = Sig(find((Sig(:,2)==2+8 | Sig(:,2)==4+8)&(Sig(:,1)~=2+8 & Sig(:,1)~=4+8 )) ,:)
CrossExperimen = [Si6_X;Si6_Y;Si8_X;Si8_Y;Si10_X;Si10_Y;]
matrix2latex(CrossExperimen,'tmp.tex')
Ctime
Ctime(1)
for K=1:3;CtimeExp1Dyn2Sim3(K);end
for K=1:3;Ctime(K);end
%-- 2020/9/1 0:18 --%
Ctime(1)
for K=1:3;Ctime(K);end
matrix2latex
for K=1:3;Ctime(K);end
Ctime(1)
for K=1:3;Ctime(K);end
Ctime(1)
clear
for K=1:3;Ctime(K);end
Ctime(1)
Ctime()
%-- 2020/9/1 14:42 --%
Get8x8Data202008
Get8x8Data
DistributionEvol
matrix2latex(DisEvolution,'Exp_distribution.tex')
roundn(DistributionEvol,3)
roundn(DisEvolution,3)
roundn(DisEvolution,-3)
r3 = roundn(DisEvolution,-3)
matrix2latex(r3,'Exp_distribution.tex')
%-- 2020/9/1 20:55 --%
Ctime()
matrix2latex
Ctime()
P1=[mean(r(1000,1:8) - r(500,1:8))]'/500
P1=[ (r(1000,1:8) - r(500,1:8))]'/500
Ctime()
P1=[(r(1000,1:8)-r(500,1:8)) (r(1000,9:16)-r(500,9:16))]'/500
P1=[(r(1000,1:8)-r(500,1:8)); (r(1000,9:16)-r(500,9:16))]'/500;
Ctime()
P2 = [[ones(1:8) ones(1:8) ones(1:8)]' [1:8 1:8 1:8]' P2]
P2 = [[ones(1,8) ones(1,8) ones(1,8)]' [1:8 1:8 1:8]' P2]
matrix2latex(P2,'Exp_Dyn_Sim_distribution.tex')
P2 = [[ones(1,8) ones(1,8)*2 ones(1,8)*3]' [1:8 1:8 1:8]' P2]
matrix2latex(P2,'Exp_Dyn_Sim_distribution.tex')
Ctime()
matrix2latex
Ctime()
r = fpayoff88(2,1)
r = payoff88(2,1)
FindNash8x8_20200901
Ctime()
StateTrans
fLatex = [ret(:,10) ret(:,1:8)]
matrix2latex(fLatex,'Exp_vortex.tex')
A = payoff88(2,1); Lambda=128; repeat=100;  wlogit=0.02; %6
[rM64x64 v disMarkov]= GenTheoMarkoW20190408(A,Lambda,repeat,wlogit)
TheoVortecInMarkov20190419
matrix2latex(f4set,'Dyn_vortex4.tex')
TheoVortecInMarkov20190419
%-- 2020/9/2 16:41 --%
TheoVortecInMarkov20190419
GenTheoMarkoW20190408
TheoVortecInMarkov20190419
sum(sum(M64x64))
sum(v)
sum(sum(T64))
TheoVortecInMarkov20190419
%-- 2020/9/2 20:19 --%
plotK20190408
plot(r(:,1:8),'DisplayName','r(:,1:8)')
plotK20190408
plot(r(:,1:8),'DisplayName','r(:,1:8)')
plotK20190408
plot(r(:,1:8),'DisplayName','r(:,1:8)')
plotK20190408
plot(r(:,1:8),'DisplayName','r(:,1:8)')
plotK20190408
plot(r(:,1:8),'DisplayName','r(:,1:8)')
plotK20190408
plot(r(:,1:8),'DisplayName','r(:,1:8)')
plotK20190408
%-- 2020/9/3 17:48 --%
A=[0.707	0.051	0.85	0.026	0.928	0.019	1	0.033
0.17	0.186	0.083	0.238	0.032	0.244	0	0.186
0.012	0.02	0	0.054	0	0.015	0	0.041
0.011	0.636	0	0.498	0.001	0.659	0	0.573
0.039	0.026	0.061	0.006	0.036	0.001	0	0.02
0.016	0.025	0.006	0.054	0.002	0.005	0	0.049
0.025	0.011	0	0.012	0	0.002	0	0.023
0.019	0.045	0	0.112	0	0.055	0	0.074
0.912	0.039	0.887	0.079	0.9	0.046	1	0.058
0.049	0.411	0.073	0.411	0.062	0.498	0	0.422
0	0.022	0	0.061	0	0.017	0	0.043
0.002	0.4	0	0.319	0	0.404	0	0.33
0.031	0.042	0.037	0.012	0.036	0.001	0	0.025
0.001	0.03	0.003	0.061	0.001	0.009	0	0.055
0.001	0.01	0	0.009	0	0.001	0	0.021
0.004	0.047	0	0.048	0	0.025	0	0.046
];
[p h]=signrank(A(:,2),A(:,8))
scatter(A(:,2),A(:,8))
[p h]=ttest(A(:,2),A(:,8))
%-- 2020/9/4 0:10 --%
plotK20190408
plot(r(:,1:8),'DisplayName','r(:,1:8)')
function NE = FindNash8x8()
NE=[];
M=payoff8820190408(2,1); N = -M; f = latex2MxWithMxPrecision(M, 0)
[A,B,a,b,iterations,err,ms]=bimat(M,N)
NE=[NE; A;B];
M=payoff8820190408(3,2); N = -M;
[A,B,a,b,iterations,err,ms]=bimat(M,N)
NE=[NE; A;B];
M=payoff8820190408(4,2); N = -M;
[A,B,a,b,iterations,err,ms]=bimat(M,N)
NE=[NE; A;B];
g = latex2MxWithMxPrecision(NE, 3)
end
M=payoff8820190408(2,1); N = -M; f = latex2MxWithMxPrecision(M, 0)
[A,B,a,b,iterations,err,ms]=bimat(M,N)
FindNash8x8_20200901
bimat
FindNash8x8_20200901
%-- 2020/9/4 21:46 --%
[]=eig([0.459939117	0.047377516	0.06189479	0.050662668
0.047377516	0.360108065	0.156140561	0.125943673
0.06189479	0.156140561	0.273543883	0.165534118
0.050662668	0.125943673	0.165534118	0.344487148
])
[v d]=eig([0.459939117	0.047377516	0.06189479	0.050662668
0.047377516	0.360108065	0.156140561	0.125943673
0.06189479	0.156140561	0.273543883	0.165534118
0.050662668	0.125943673	0.165534118	0.344487148
])
[v d]=eig([0.4229	0.1789	0.0761	0.2175
0.1789	0.4656	0.0675	0.1929
0.0761	0.0675	0.5959	0.0818
0.2175	0.1929	0.0818	0.3905
])
[v d]=eig([0.4229	0.1789	0.0761
0.1789	0.4656	0.0675
0.0761	0.0675	0.5959
])
[v d]=eig([0.4656	0.0675	0.1929
0.0675	0.5959	0.0818
0.1929	0.0818	0.3905
])
v(:,1)'*v(:,2)
[v d]=eig([0.459939117	0.047377516	0.06189479	0.050662668
0.047377516	0.360108065	0.156140561	0.125943673
0.06189479	0.156140561	0.273543883	0.165534118
0.050662668	0.125943673	0.165534118	0.344487148
])
v(:,1)'*v(:,2)
v(:,1)'*v(:,3)
v(:,1)'*v(:,4)
%-- 2020/9/5 14:04 --%
Ctime()
matrix2latex
Ctime()
A=[0.061
0.188
0.03
0.03
0.03
0.02
0.063
];[t p] = ttest(A,0,579)
];[t p] = ttest(A-0,579)
[t p] = ttest(A-0,579)
b.*[1 1]
b = [b b b b]
A=[0.2
0.3
0.4
0.3
0.3
]
[t p] = ttest(A-0.579)
A=[];
A=[0.2
0.3
0.7
0.3
0.3
];
[t p] = ttest(A-0.9)
[t p] = ttest(A-0.8)
A
[t p] = ttest(A-0.579)
A=[0.2
0.3
0.5
0.3
0.3
];
[t p] = ttest(A-0.579)
Ctime()
clear
Ctime()
%-- 2020/9/5 19:01 --%
Ctime()
matrix2latex
Ctime()
load('PS26.mat')
A=PS(find(:,1)==6,:);
A=PS(find(PS(:,1)==6),:);
plot(A(:,4:11),'DisplayName','A(:,4:11)')
plot(A(:,7))
plot(A(:,[7,11]),'DisplayName','A(:,[7,11])')
plot(A(:,4:11),'DisplayName','A(:,4:11)')
plot(A(:,[7,11]),'DisplayName','A(:,[7,11])')
plot(A(:,13:end),'DisplayName','A(:,13:end)')
plot(A(:,4:11),'DisplayName','A(:,4:11)')
A=PS(find(PS(:,1)==10),:);
plot(A(:,4:11),'DisplayName','A(:,4:11)')
plot(A(:,13:19),'DisplayName','A(:,13:19)')
%-- 2020/9/6 0:26 --%
load('PS26.mat')
T=PS((find(abs(PS(:,3)-40)>0.001 & PS(:,1)==6)),[4:11 13:20]);
T(:,17:32)=99;
for i=1:12000; for j=1:8; T(i,16+j)=sum(T(i,1:j));T(i,24+j)=sum(T(i,9:8+j));T(i,33:34)=rand(1,2);end;end
histogram(T(:,end))
min(find(T(1,17:24)>T(1,33)))
[min(find(T(1,17:24)>T(1,33))) min(find(T(1,25:32)>T(1,34)))]
for i=1:12000; T(i,35:36) = [min(find(T(i,17:24)>T(i,33))) min(find(T(i,25:32)>T(i,34)))] ;end;
for i=1:12000; T(i,37) =  T(i,35)*10+ T(i,36);end;
for i=1:12000-1; T(i,38) =  T(i+1,37);end;
T8=PS((find(abs(PS(:,3)-40)>0.001 & PS(:,1)==8)),[4:11 13:20]);
T8(:,17:32)=99;
for i=1:12000; for j=1:8; T8(i,16+j)=sum(T8(i,1:j));T8(i,24+j)=sum(T8(i,9:8+j));T8(i,33:34)=rand(1,2);end;end
for i=1:12000; T8(i,35:36) = [min(find(T8(i,17:24)>T8(i,33))) min(find(T8(i,25:32)>T8(i,34)))] ;end;
for i=1:12000; T8(i,37) =  T8(i,35)*10+ T8(i,36);end;
for i=1:12000-1; T8(i,38) =  T8(i+1,37);end;
T8=PS((find(abs(PS(:,3)-40)>0.001 & PS(:,1)==10)),[4:11 13:20]);
T8(:,17:32)=99;
for i=1:12000; for j=1:8; T8(i,16+j)=sum(T8(i,1:j));T8(i,24+j)=sum(T8(i,9:8+j));T8(i,33:34)=rand(1,2);end;end
for i=1:12000; T8(i,35:36) = [min(find(T8(i,17:24)>T8(i,33))) min(find(T8(i,25:32)>T8(i,34)))] ;end;
for i=1:12000; T8(i,37) =  T8(i,35)*10+ T8(i,36);end;
for i=1:12000-1; T8(i,38) =  T8(i+1,37);end;
%-- 2020/9/6 10:58 --%
TheoVortecInMarkov20190419
StateTrans
save('Exp_count_Markov.mat',TransMatrixOrig)
save('Exp_count_Markov.mat','TransMatrixOrig')
clear
load('Exp_count_Markov.mat')
[r3vortex r4votex] = Markov2Vortex(TransMatrixOrig)
%-- 2020/9/6 14:32 --%
load('Exp_count_Markov.mat')
[r3vortex r4votex] = Markov2Vortex(TransMatrixOrig)
Markov2Vortex
[r3vortex r4votex] = Markov2Vortex(TransMatrixOrig)
Timeseries2Markov
[ret4(:,12) ret4(:,1:10)]
Timeseries2Markov
[ret4(:,12) ret4(:,1:10)]
Timeseries2Markov
Markov2Vortex
Timeseries2Markov
sum(sum(TransMatrixOrig))
clear
Timeseries2Markov
A=[0.3101; 0.4455;  0.4476; 0.2011; 0.6816]
A/sum(A)
[-0.2235; -0.3988;-0.3526;0.1781;0.7968]
C= [0.0612;-0.4021;0.4647;-0.6146;0.4908]
B=[-0.2235; -0.3988;-0.3526;0.1781;0.7968]
B'*C
B2=B/sqrt(B'*B)
C2=C/sqrt(C'*C)
C2'*B2
B2-B
B2'*A
a=[0.3147+0.2844i; 0.3939-0.3464i; -0.6364; -0.1077+0.2811i; 0.0355-0.2191i]
y=[];for t=0:0.01:1; x=a'*exp(0.67+0.43i);y=[y;x];end
y=[];for t=0:0.01:1; x=a'*exp((0.67+0.43i)*t);y=[y;x];end
yr=real(y)
plot(yr,'DisplayName','yr')
grid on
y=[];for t=0:0.01:10; x=a'*exp((0.67+0.43i)*t);y=[y;x];end
yr=real(y)
plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:1; x=a'*exp((-0.67+0.43i)*t);y=[y;x];end
yr=real(y)
plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:10; x=a'*exp((-0.67+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
grid on
y=[];for t=0:0.01:40; x=a'*exp((-0.67+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
grid on
ylim([-0.05 0.05])
ylim([-0.005 0.005])
%-- 2020/9/6 20:44 --%
y=[];for t=0:0.01:40; x=a'*exp((-0.67+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
a=[0.3147+0.2844i; 0.3939-0.3464i; -0.6364; -0.1077+0.2811i; 0.0355-0.2191i]
y=[];for t=0:0.01:40; x=a'*exp((-0.67+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:40; x=a'*exp((-0.0067+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
grid on
y=[];for t=0:0.01:40; x=a'*exp((-0.67+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:40; x=a'*exp((-0.0067+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:40; x=a'*exp((-0.0000+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:40; x=a'*exp((-0.0067+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:40; x=a'*exp((-0.067+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:40; x=a'*exp((-0.67+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
y=[];for t=0:0.01:40; x=a'*exp((-0.0067+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
grid on
y=[];for t=0:0.01:40; x=a'*exp((-0.0000+0.43i)*t);y=[y;x];end;yr=real(y);plot(yr,'DisplayName','yr')
grid on
%-- 2020/9/7 0:32 --%
Timeseries2Markov
Get8x8Data202008
figure(8)
Timeseries2Markov
%-- 2020/9/7 1:21 --%
Timeseries2Markov
StateTrans
%-- 2020/9/7 14:05 --%
Timeseries2Markov
%-- 2020/9/7 15:09 --%
NE=[];a=xlsread('F:\cDownload\abed-1pop-master\parameter-files\pop1-5strategies.xlsx')
a=a(1:14200,:);
NE=mean(a)
c=[];for i=1:length(a(:,1)); b=a(i,:)-NE;c=[c;b];end
d=sign(c);
d=(sign(c)+1)/2;
d(:,6)=(sign(c)+1)/2+900000;
d(:,6)=0
d(:,6)=(sign(c)+1)/2+900000;
c=[];for i=1:length(a(:,1)); d(i,6)=1000*d(i,2)+100*d(i,3)+10*d(i,4)+1*d(i,5);end
c=[];for i=1:length(a(:,1)); d(i,6)=10000*d(i,1)+1000*d(i,2)+100*d(i,3)+10*d(i,4)+1*d(i,5)+900000;end
NE=[0.148644125
0.213592233
0.214596585
0.096417811
0.326749247]
c=[];for i=1:length(a(:,1)); b=a(i,:)-NE';c=[c;b];end
d=(sign(c)+1)/2;
f=unique(d)
f=unique(d,[1 5])
f=unique(d,[1:5])
f=unique(d(:,[1:5]))
f=unique(d,'rows')
e=d(2:end,:);
g=d(1:end-1,:)
M=zeros(30);for i=1:30;forj=1:30;
M=length(T8(find(e(:,1:5) == f(i,:) & g(:,38) == f(j,:)),1));
end
end
for i=1:lengthe(e(:,1)); m=find(f==e(i,:)); n=find(f==g(i,:));M(m,n)=M(m,n)+1;end
for i=1:length(e(:,1)); m=find(f==e(i,:)); n=find(f==g(i,:));M(m,n)=M(m,n)+1;end
Untitledddd
aM=M-M'
aM=(abs(M-M')+(M-M'))/2
Untitledddd
aM=(oM-oM')/2
%-- 2020/9/8 0:13 --%
Timeseries2Markov
plot(T8(:,1:8),'DisplayName','T8(:,1:8)')
plot(T8(:,[2,10]),'DisplayName','T8(:,[2,10])')
Timeseries2Markov
plot(T8(:,[2,10]),'DisplayName','T8(:,[2,10])')
Timeseries2Markov
scatter(getcolumn(T8(1:1000,[2,10]),1),getcolumn(T8(1:1000,[2,10]),2))
surf(T8(1:1000,[2,10]))
plot(T8(:,[2,10]),'DisplayName','T8(:,[2,10])')
scatter(getcolumn(T8(1:1000,[2,10]),1),getcolumn(T8(1:1000,[2,10]),2))
plot(getcolumn(T8(1:1000,[2,10]),1),getcolumn(T8(1:1000,[2,10]),2))
for i=1:900;plot(T8(i,[2,10]));hold on;wait(0.1);end
for i=1:900;plot(T8(i,[2,10]));hold on;palse(0.1);end
for i=1:900;plot(T8(i,[2,10]));hold on;pause(0.1);end
for i=1:900;scatter(T8(i,[2,10]));hold on;pause(0.1);end
for i=1:900;scatter(T8(i,2),T8(i,10));hold on;pause(0.1);end
for i=1:900;scatter(T8(i,2),T8(i,10));hold on;pause(0.1);end;xlim([0 1]);ylim([0 1])
for i=1:900;scatter(T8(i,2),T8(i,10));hold on;pause(0.1);xlim([0 1]);ylim([0 1]);end;
for i=1:900;scatter(T8(i,2),T8(i,10),'.');hold on;pause(0.1);xlim([0 1]);ylim([0 1]);end;
%-- 2020/9/9 15:16 --%
load('PE26.mat');
f100 = PS(find(PS(:,3)<100),:)
f100 = PS(find(PS(:,3)<=100),:)
f100 = PS(find(PS(:,3)<=100 & PS(:,1)<=6),:)
mean(f100)'
mean(f100(:,[4 20]))'
mean(f100(:,[4:20]))
mean(f100(:,[4:11]))
f100 = PS(find(PS(:,3)<=200 & PS(:,1)<=6),:)
mean(f100(:,[4:11]))
f200 = PS(find(PS(:,3)<=200 & PS(:,1)<=6),:);
f100 = PS(find(PS(:,3)<=100 & PS(:,1)==6),:)
f100 = PS(find(PS(:,3)<=100 & PS(:,1)==6),:);
f200 = PS(find(PS(:,3)<=200 & PS(:,1)==6),:);
mean(f100(:,[4:11]))
s = [mean(f100(:,[4:11])) mean(f100(:,[4:11]))]'
s = [mean(f100(:,[4:11])); mean(f200(:,[4:11]))]'
f200 = PS(find(PS(:,3)<=200 & PS(:,3)>100 & PS(:,1)==6),:);
s = [mean(f100(:,[4:11])); mean(f200(:,[4:11]))]'
f300 = PS(find(PS(:,3)<=300 & PS(:,3)>200 & PS(:,1)==6),:);
s = [mean(f100(:,[4:11])); mean(f200(:,[4:11])); mean(f300(:,[4:11]))]'
f100 = PS(find(PS(:,3)<=100 & PS(:,1)==6),:)f200 = PS(find(PS(:,3)<=200 & PS(:,3)>100 & PS(:,1)==6),:);f300 = PS(find(PS(:,3)<=300 & PS(:,3)>200 & PS(:,1)==6),:);s = [mean(f100(:,[4:11])); mean(f200(:,[4:11])); mean(f300(:,[4:11]))]'
f100 = PS(find(PS(:,3)<=100 & PS(:,1)==6),:);f200 = PS(find(PS(:,3)<=200 & PS(:,3)>100 & PS(:,1)==6),:);f300 = PS(find(PS(:,3)<=300 & PS(:,3)>200 & PS(:,1)==6),:);s = [mean(f100(:,[4:11])); mean(f200(:,[4:11])); mean(f300(:,[4:11]))]'
f100 = PS(find(PS(:,3)<=100 & PS(:,1)==8),:);f200 = PS(find(PS(:,3)<=200 & PS(:,3)>100 & PS(:,1)==8),:);f300 = PS(find(PS(:,3)<=300 & PS(:,3)>200 & PS(:,1)==8),:);s = [mean(f100(:,[4:11])); mean(f200(:,[4:11])); mean(f300(:,[4:11]))]'
t=10;f100 = PS(find(PS(:,3)<=100 & PS(:,1)==t),:);f200 = PS(find(PS(:,3)<=200 & PS(:,3)>100 & PS(:,1)==t),:);f300 = PS(find(PS(:,3)<=300 & PS(:,3)>200 & PS(:,1)==t),:);s = [mean(f100(:,[4:11])); mean(f200(:,[4:11])); mean(f300(:,[4:11]))]'
%-- 2020/9/9 18:41 --%
Ctime()
matrix2latex
Ctime()
plot(Accu_P1,'DisplayName','Accu_P1')
xlim([0 200])
xlim([0 200]);ylim([0 300])
xlim([0 200]);ylim([0 30])
plot(Accu_P2,'DisplayName','Accu_P1')
xlim([0 200]);ylim([0 30])
plot(Accu_P1,'DisplayName','Accu_P1')
xlim([0 200]);ylim([0 30])
%-- 2020/9/9 22:15 --%
StateTrans
aa=[6	7	16	9	10	11	12	14		55	80
24		8			11
14	8	9	10	22	12		15
]
aa=[6	7	16	9	10	11	12	14	0	55	80
24	0	8	0	0	11	0	0	0	0	0
0	14	8	9	10	22	12	0	15	0	0
]
aa'>
a=aa'
[h p]=ttest(a(:,1)-a(:,3),0)
[h p]=ttest(a(:,1)-a(:,2),0)
[h p]=signrank(a(:,1),a(:,2))
[h p]=ttest(a(:,1)-a(:,3),0)
[h p]=ttest(a(:,2)-a(:,3),0)
a(:,1)-a(:,3)
b=a(:,1)-a(:,3)
[h p]=ttest(b)
stdev(a)
std(a)
mean(a)
StateTrans
Timeseries2Markov
Markov2Vortex
%-- 2020/9/11 19:01 --%
run('F:\cDownload\PokerAnan20200716\8x8abedata\eigenVV8x8.m')
eigenVV8x8
D_Eq_0
clear
clean
D_Eq_0
eigenVV8x8
D_Eq_0
payoff_matrix_1 = [0,0,0,0,0,0,0,0
-2,-2,1,1,1,1,4,4
2,-3,2,-3,5,0,5,0
0,-5,3,-2,6,1,9,4
6,1,1,-4,6,1,1,-4
4,-1,2,-3,7,2,5,0
8,-2,3,-7,11,1,6,-4
6,-4,4,-6,12,2,10,0];
payoff_matrix_2 = -payoff_matrix_1
syms x(1:8) y(1:8)  real
syms x(8) y(8)  real
help syms
eigenVV8x8
save('eigenVV8x8.mat' )
clear
load('eigenVV8x8.mat' )
v=V_Eq_0;
p=0.125*ones(16,1)
x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);	y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16);
eval(v)
sum(v)
sum(eval(v))
eval(v(1))
T=[];p=0.125*ones(16,1);t=0:0.001:3; p=p+eval(v)*0.001;T=[T;eval(v)];end
T=[];p=0.125*ones(16,1);for t=0:0.001:3; p=p+eval(v)*0.001;T=[T;eval(v)];end
T=[];p=0.125*ones(16,1);for t=0:0.001:1; p=p+eval(v)*0.001;T=[T;p];end
T=[];p=0.125*ones(16,1);for t=0:0.001:1; p=p+eval(v)*0.001;T=[T;p'];end
plot(T(:,1:8),'DisplayName','T(:,1:8)')
v
sum(v(1:8))
simply(sum(v(1:8)))
simplify(sum(v(1:8)))
T=[];p=0.125*ones(16,1);for t=1:100; p = p + eval(v)*0.001;T=[T;p'];end
plot(T(:,9:end),'DisplayName','T(:,9:end)')
T=[];p=0.125*ones(16,1);for t=1:100; x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);	y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16); p = p + eval(v)*0.001;T=[T;p'];end
plot(T(:,9:end),'DisplayName','T(:,9:end)')
T=[];p=0.125*ones(16,1);for t=1:100; x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);	y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16); p = p + eval(v)*0.005;T=[T;p'];end
plot(T(:,9:end),'DisplayName','T(:,9:end)')
plot(T(:,1:8),'DisplayName','T(:,1:8)')
T=[];p=0.125*ones(16,1);for t=1:100; x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);	y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16); p = p + eval(v)*0.01;T=[T;p'];end
plot(T(:,1:8),'DisplayName','T(:,1:8)')
T=[];p=0.125*ones(16,1);for t=1:1000; x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);	y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16); p = p + eval(v)*0.01;T=[T;p'];end
plot(T(1:108,1:8),'DisplayName','T(1:108,1:8)')
plot(T(:,9:end),'DisplayName','T(:,9:end)')
plot(T(:,1:8),'DisplayName','T(:,1:8)')
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
eigenVV8x8
payoff88
eigenVV8x8
T=[];p=0.125*ones(16,1);for t=1:1000; x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);	y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16); p = p + eval(v)*0.01;T=[T;p'];end
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
eigenVV8x8
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
plot(T(:,9:end),'DisplayName','T(:,9:end)')
plot(T(:,1:8),'DisplayName','T(:,1:8)')
plot(T(:,9:end),'DisplayName','T(:,9:end)')
clear
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
eigenVV8x8
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
v=V_Eq_0;
T=[];p=0.125*ones(16,1);for t=1:1000;
x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);
y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16);
p = p + eval(v)*0.01;T = [T;p'];
end
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
%-- 2020/9/11 23:22 --%
eigenVV8x8
Timeseries2Markov
eigenVV8x8
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
v=V_Eq_0;
T=[];p=0.125*ones(16,1);for t=1:1000;
x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);
y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16);
p = p + eval(v)*0.001;T = [T;p'];
end
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
eigenVV8x8
plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
new figure;plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
figure new;plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
figure  ;plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
eigenVV8x8
figure  ;  plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
figure  ;  plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
eigenVV8x8
figure  ;  plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
figure  ;  plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
eigenVV8x8
figure  ;  plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
figure  ;  plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
%-- 2020/9/12 23:44 --%
a=[(-0.0974
0.75064
-1.6532
1)
;(-0.00232+0.65265 I
0.16909+0.18737 I
-1.1668-0.8400 I
1);
(-0.00232-0.65265 I
0.16909-0.18737 I
-1.1668+0.8400 I
1)]
a=[ -0.0974 0.75064 -1.6532 1 ;
-0.00232+0.65265i 0.16909+0.18737i -1.1668-0.8400i 1;
-0.00232-0.65265i 0.16909-0.18737i -1.1668+0.8400i 1 ]
a=a'
a(:,1)'*a(:,2)
dot(a(:,1),a(:,2))
dot(a(:,1),a(:,3))
dot(a(:,2),a(:,3))
dot(a(:,2),a(:,2))
dot(a(:,2)',a(:,3))
a(:,2)'
a
%-- 2020/9/13 15:30 --%
%双人群3策略稳定博弈
%本征值、本征向量计算
syms x1 x2 x3  y1 y2 y3  real
payoff_matrix_1 = [4 10 12 ; 15 0 15 ; 18 0 14  ]
payoff_matrix_2 = [0 12 16 ; 15 6 10 ; 10 12 8  ]
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2 y3]'
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2 x3]'
mean_U_1 = [x1 x2 x3 ] * Payoff_vector_field_F_1
mean_U_2 = [y1 y2 y3] * Payoff_vector_field_F_2
V_Eq_1 = [x1 x2 x3]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_2 = [y1 y2 y3]'.*(Payoff_vector_field_F_2 - mean_U_2);
V_Eq_0 = [V_Eq_1; V_Eq_2]
D_Eq_0 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3')  ...
]
for i = 1:6
'column-i of the charactor matrix)'
D_Eq_0(:,i)
end
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3  S.y1 S.y2 S.y3 ]
%下面将纳什均衡填入
A=eval([S.x1 S.x2 S.x3  S.y1 S.y2 S.y3 ]);
A=A(1:12,:)
x1=A(length(A(:,1)),1); x2=A(length(A(:,1)),2);
x3=A(length(A(:,1)),3);
y1=A(length(A(:,1)),4); y2=A(length(A(:,1)),5);
y3=A(length(A(:,1)),6);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
matrix2latex(eigen_value,'tmp.txt')
matrix2latex(eigen_vector,'tmp.txt')
latex2MxWithMxPrecision(eigen_vector,3)
round(eigen_vector,3)
round(eigen_vector,2)
round(eigen_vector,3)
a = round(eigen_vector,3)
matrix2latex(a,'tmp.txt')
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny')
columnLabels=['$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_6$' ]
matrix2latex(a,'tmp.txt',  'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny')
columnLabels={'$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_6$' }
matrix2latex(a,'tmp.txt',  'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny')
rowLabels={'$\mu_1$' '$\mu_1$' '$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_6$' }
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny')
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
rowLabels={'$\mu_1$' '$\mu_2$' '$\mu_1$' '$\mu_1$' '$\mu_1$' '$\mu_6$' }
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
rowLabels={'$\mu_1$' '$\mu_2$' '$\mu_3$' '$\mu_4$' '$\mu_5$' '$\mu_6$' }
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
rowLabels={'$\mu_1$' '$\mu_1$' '$\xi_1$' '$\xi_1$' '$\xi_1$' '$\xi_6$' }
columnLabels={'$\xi_1$' '$\xi_2$' '$\xi_3$' '$\xi_4$' '$\xi_5$' '$\xi_6$' }
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
diag(eigen_value)
diag(eigen_value)'
matrix2latex(diag(eigen_value)','tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
columnLabels={'$\xi_1$' '$\xi_2$' '$\xi_3$' '$\xi_4$' '$\xi_5$' '$\xi_6$' }
rowLabels={'$\mu_1$' '$\mu_2$' '$\mu_3$' '$\mu_4$' '$\mu_5$' '$\mu_6$' }
matrix2latex(a,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
matrix2latex( eigen_vector,'tmp.txt', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small')
%-- 2020/9/13 21:19 --%
f='C:\Users\Think\Desktop\abed-1pop5replicator Strategy distributions (complete history).csv'
data=readAbed(csvfile)
data=readAbed(f)
f='C:\Users\Think\Desktop\abed-1pop5replicator Strategy distributions (complete history).csv'
data=readAbed(f)
plot(data,'DisplayName','data')
NE=[0.3101;   0.4455;   0.4476;  0.2011;  0.6816]/sum[0.3101;   0.4455;   0.4476;  0.2011;  0.6816]
NE=[0.3101;   0.4455;   0.4476;  0.2011;  0.6816]/sum([0.3101;   0.4455;   0.4476;  0.2011;  0.6816])
b=[];for k=1:5378;b=[b;NE'];end
c=data-b;
plot(c(61,:))
plot(c,'DisplayName','c')
f='C:\Users\Think\Desktop\abed-1pop5replicator Strategy distributions (complete history).csv'
data=readAbed(f)
b=[];for k=1:length(data(:,1));b=[b;NE'];end
b=[];for k=1:length(data(:,1));b=[b;NE'];end;;c=data-b;
plot(c,'DisplayName','c')
plot(c,'DisplayName','c','LineWidth',2)
histogram(c)
scatter(getcolumn(c(:,[2,end]),1),getcolumn(c(:,[2,end]),2))
scatter(getcolumn(c(:,4:end),1),getcolumn(c(:,4:end),2))
scatter(getcolumn(c(:,2:3),1),getcolumn(c(:,2:3),2))
scatter(getcolumn(c(:,4:end),1),getcolumn(c(:,4:end),2))
scatter(getcolumn(c(:,1:2),1),getcolumn(c(:,1:2),2))
scatter(getcolumn(c(:,[1,3]),1),getcolumn(c(:,[1,3]),2))
scatter(getcolumn(c(:,[1,4]),1),getcolumn(c(:,[1,4]),2))
scatter(getcolumn(c(:,[1,end]),1),getcolumn(c(:,[1,end]),2))
data=readAbed(f)
b=[];for k=1:length(data(:,1));b=[b;NE'];end;;c=data-b;
plot(c,'DisplayName','c','LineWidth',2)
bar(c,'DisplayName','c')
histogram(c)
histogram(c(:,1))
bar(c(:,1:2),'DisplayName','c(:,1:2)')
histogram(c(:,1:2),'DisplayName','c(:,1:2)')
histogram(c(:,[1 2]),'DisplayName','c(:,1:2)')
histogram(c(:,[1]),'DisplayName','c(:,1:2)')
hold on;histogram(c(:,[2]),'DisplayName','c(:,1:2)')
for k=1:5;hold on;histogram(c(:,[k]));end
for k=1:5;subplot(2,3,k);hold on;histogram(c(:,[k]));end
f='C:\Users\Think\Desktop\abed-1pop5replicator Strategy distributions (complete history).csv'
data=readAbed(f)
% 构造 NE
NE=[0.3101;   0.4455;   0.4476;  0.2011;  0.6816]';
NE = mean(data);
NE = NE/sum(NE);
b=[];for k=1:length(data(:,1));b=[b;NE];end;;
% 构造差值
c=data-b;
% 绘制围绕NE振动图（1 时间序列）
plot(c,'DisplayName','c','LineWidth',2)
% 绘制围绕NE振动图（2 偏差的直方图(Histogram)）
for k=1:5;subplot(2,3,k);hold on;histogram(c(:,[k]));end
clear
clf
f='C:\Users\Think\Desktop\abed-1pop5replicator Strategy distributions (complete history).csv'
data=readAbed(f)
% 构造 NE
NE=[0.3101;   0.4455;   0.4476;  0.2011;  0.6816]';
NE = mean(data);
NE = NE/sum(NE);
b=[];for k=1:length(data(:,1));b=[b;NE];end;;
% 构造差值
c=data-b;
% 绘制围绕NE振动图（1 时间序列）
plot(c,'DisplayName','c','LineWidth',2)
% 绘制围绕NE振动图（2 偏差的直方图(Histogram)）
for k=1:5;subplot(2,3,k);hold on;histogram(c(:,[k]));end
plot(c,'DisplayName','c','LineWidth',2)
grid on
for k=1:5;subplot(2,3,k);hold on;histogram(c(:,[k]));end
f='C:\Users\Think\Desktop\abed-1pop5replicator Strategy distributions (complete history).csv'
data=readAbed(f)
% 构造 NE
NE=[0.3101;   0.4455;   0.4476;  0.2011;  0.6816]';
NE = mean(data);
NE = NE/sum(NE);
b=[];for k=1:length(data(:,1));b=[b;NE];end;;
% 构造差值
c=data-b;
% 绘制围绕NE振动图（1 时间序列）
plot(c,'DisplayName','c','LineWidth',2)
% 绘制围绕NE振动图（2 偏差的直方图(Histogram)）
for k=1:5;subplot(2,3,k);hold on;histogram(c(:,[k]));end
plot(c,'DisplayName','c','LineWidth',2)
figure;plot(c,'DisplayName','c','LineWidth',2)
figure;for k=1:5;subplot(2,3,k);hold on;histogram(c(:,[k]));end
[a ]=std(c)
std(c)
R = corrcoef(c)
R = corrcoef(c(:,[4:5]))
mean(c)
f='C:\Users\Think\Desktop\abed-1pop5replicator Strategy distributions (complete history).csv'
data=readAbed(f)
% 构造 NE
NE=[0.3101;   0.4455;   0.4476;  0.2011;  0.6816]';
%NE = mean(data);
NE = NE/sum(NE);
b=[];for k=1:length(data(:,1));b=[b;NE];end;;
% 构造差值
c=data-b;
% 绘制围绕NE振动图（1 时间序列）
plot(c,'DisplayName','c','LineWidth',2)
% 绘制围绕NE振动图（2 偏差的直方图(Histogram)）
figure;for k=1:5;subplot(2,3,k);hold on;histogram(c(:,[k]));end
R = corrcoef(c)
c(:,6)=c(:,2)+c(:,4);
R = corrcoef(c)
c(:,6)=(c(:,1)+c(:,4));
R = corrcoef(c)
std(c)
plot(c(:,4:end),'DisplayName','c(:,4:end)')
plot(c(:,1:3),'DisplayName','c(:,1:3)')
v=[0.3147-0.2844*i
0.3939+0.3464*i
-0.6364
-0.1077-0.2811*i
0.0355+0.2191*i]
plot(v)
plot(v,'.')
plot(v,'o')
u=[c(:,3) (c(:,2)+c(:,5))  (c(:,4)+c(:,1)) ]
plot(u,'DisplayName','u')
%-- 2020/9/18 0:51 --%
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 4 ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ]
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]'
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1]
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ]
for i = 1:4
'column-i of the charactor matrix)'
D_Eq_0(:,i)
end
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3 S.x4]
%下面将纳什均衡填入
x1=4/13;x2=4/13;x3=4/13;x4=1/13;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
ZSJ4x44
[S.x1 S.x2 S.x3 S.x4]
ZSJ4x44
[S.x1 S.x2 S.x3 S.x4]
S.x1
ZSJ4x44
%-- 2020/9/20 19:24 --%
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 4 ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ]
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]'
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1]
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ]
for i = 1:4
'column-i of the charactor matrix)'
D_Eq_0(:,i)
end
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3 S.x4]
%下面将纳什均衡填入
x1=4/13;x2=4/13;x3=4/13;x4=1/13;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 1/4 ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ]
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]'
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1]
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ]
for i = 1:4
'column-i of the charactor matrix)'
D_Eq_0(:,i)
end
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3 S.x4]
%下面将纳什均衡填入
x1=4/13;x2=4/13;x3=4/13;x4=1/13;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
a = [1/4 1 4]
for K=1:3
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 a(K) ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ];
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]';
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1;
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1];
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ];
%下面将纳什均衡填入
x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
end
a = [1/4 1 4]
for K=1:3
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0.0001 0 0 a(K) ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ];
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]';
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1;
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1];
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ];
%下面将纳什均衡填入
x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
end
a = [1/4 1 4]
for K=1:3
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 a(K) ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ];
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]';
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1;
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1];
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ];
%下面将纳什均衡填入
x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
end
help eig
a = [1/4 1 4]
for K=1:3
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 a(K) ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ];
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]';
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1;
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1];
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ];
%下面将纳什均衡填入
x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
end
a = [1/4 1 4]
for K=1:3
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 a(K) ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ];
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]';
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1;
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1];
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ];
%下面将纳什均衡填入
x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE,'qz')
end
a = [1/4 1 4]
for K=1:3
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 a(K) ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ];
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]';
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1;
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1];
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ];
%下面将纳什均衡填入
x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE,'balance')
end
a = [1/4 1 4]
for K=1:3
syms x1 x2 x3 x4 real
payoff_matrix_1 = [ 0 0 0 a(K) ; 1 0 0 0  ; 0 1 0 0
0 0 1 0 ];
Payoff_vector_field_F_1 = payoff_matrix_1 *[x1 x2 x3 x4 ]';
mean_U_1 = [x1 x2 x3 x4 ] * Payoff_vector_field_F_1;
V_Eq_1 = [x1 x2 x3 x4 ]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_0 = [V_Eq_1];
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') ...
diff(V_Eq_0,'x4') ];
%下面将纳什均衡填入
x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE,'vector')
end
maple('convert(sin(x),`exp`)')
convert(sin(x),`exp`)
convert(sin(x),'exp')
syms x
convert(sin(x),'exp')
convert(sin(3),'exp')
angle(eigen_vector)
1.7976 + 1.3440
2.6880/pi
1.3440/pi
1.7976/pi
v=eigen_vector;
v(2,4)/v(4,4)
v
-0.1446/0.2374
abs(v)
z = abs(v)
z2=z(:,2)
sum(z2(1:3))
sum(z2(1:3))/z2(4)
%-- 2020/9/21 1:05 --%
replicator4x4
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3 S.x4]
S.x1(12)
x1=S.x1(12);x2=S.x2(12);x3=S.x3(12);x4=S.x4(12);
c=S(12)
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE,'vector')
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
replicator4x4
c=[S.x1 S.x2 S.x3 S.x4]
c(:,2)
c(1,:)
[ 0, 0, 0, 0]
S(9)
S(9,:)
A=eval(S(9,:))
A=eval(c(9,:))
Sx0=eval([S.x1 S.x2 S.x3 S.x4])
%下面将纳什均衡填入
%   x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
LI = Sx0(9,:);
x1=Sx(1);x2=Sx(2);x3=Sx(3);x4=Sx(4);
Sx0=eval([S.x1 S.x2 S.x3 S.x4])
%下面将纳什均衡填入
%   x1=a(K)/(3*a(K)+1);x2=x1;x3=x1;x4=1/(3*a(K)+1);
Sx = Sx0(9,:);
x1=Sx(1);x2=Sx(2);x3=Sx(3);x4=Sx(4);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
abs(eigen_value)
replicator4x4
D_Eq_at_NE^(-1)
D_Eq_at_NE^(-1) * D_Eq_at_NE
1/D_Eq_at_NE * D_Eq_at_NE
1./D_Eq_at_NE * D_Eq_at_NE
replicator4x4
D_Eq_at_NE^(-1) * D_Eq_at_NE
D_Eq_at_NE^(-1)
%-- 2020/9/21 11:06 --%
replicator4x4
payoff_matrix_1 = [ 0 0 0 4 ; 1 0 0 0  ; 0 1 0 0; 0 0 1 0 ]
replicator4x4
Sx0=eval([S.x1 S.x2 S.x3 S.x4])
Sx = Sx0(9,:)
Sx = Sx0(9,:);
x1=Sx(1);x2=Sx(2);x3=Sx(3);x4=Sx(4);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
replicator4x4
Sx = Sx0(2,:);
x1=Sx(1);x2=Sx(2);x3=Sx(3);x4=Sx(4);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
D_Eq_0
D_Eq_at_NE
sum([?0.56487, 5.6666, ?2.8375, 1])
sum([-0.56487, 5.6666, -2.8375, 1])
sum([[ 0.75041 + 3.2512i, -2.501 + 8.1447e-5i, 0.75041 - 3.2515i, 1]])
abs(eigen_vector)
abs( eigen_value)
%-- 2020/9/21 22:56 --%
replicator4x4
D_Eq_0
d=subs(D_Eq_0,'x4','1-x1-x2-x3')
sum(d)
simplify(sum(d))
replicator4x4
D_Eq_0=subs(D_Eq_0,'x4','1-x1-x2-x3')
D_Eq_0 =  [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2')  ...
diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ]
D_Eq_0=subs(D_Eq_0,'x4','1-x1-x2-x3')
diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ]
D_Eq_0 =  [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2')  ...
diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ]
D_Eq_0
D_Eq_0=subs(D_Eq_0,'x4','1-x1-x2-x3')
simplify((d))
replicator4x4
D_Eq_0 =  [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2')  ...
diff(V_Eq_0,'x3')];
D_Eq_0
V_Eq_0 = V_Eq_0(1:3,:)
D_Eq_0 =  [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2')  ...
diff(V_Eq_0,'x3')]
S=solve(V_Eq_0)
Sx0=eval([S.x1 S.x2 S.x3])
Sx = Sx0(2,:)
x1=Sx(1);x2=Sx(2);x3=Sx(3);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
D_Eq_at_NE = eval(D_Eq_0)'
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
D_Eq_at_NE = eval(D_Eq_0 + D_Eq_0')
[eigen_vector eigen_value W] = eig(D_Eq_at_NE)
D_Eq_at_NE = eval(D_Eq_0 )
[eigen_vector eigen_value ] = eig(D_Eq_at_NE)
replicator4x4
abs(eigen_vector)
abs(eigen_vector(:,4))
abs(eigen_vector(:,4)^2)
abs(eigen_vector(:,4).^2)
replicator4x4
A= [1 1 1 -2; 0 1 0 -1; 0 0 1 1; 0 0 0 1]
[v d]=eig(A)
replicator4x4
0.1446/0.2374
abs(eigen_vector(:,4).^2)
0.9496 - 0.1446
0.9496/0.707
0.9496 * 0.707
eigen_vector
0.7071 / 4
0.7071 / 0.1446
%-- 2020/9/26 11:13 --%
a=[1 0 1 0]; b=[0 4 0 1];
syms x real
c = a*x + b; d= c/sum(c);
a=[1 0 1 0]; b=[0 4 0 1];
syms x real
c = a*x + b; d= c/sum(c);
E = d.* log(d)
f = diff(E)
S=solve(f)
S=solve(f,'x')
S.x
f'
f*(2*x + 5)^2
simpliff*(2*x + 5)^2)y(
simplify(f*(2*x + 5)^2)
E = d * log(d)
E = d * log(d)'
f = diff(E)
S=solve(f)
x=S
eval(d)
(2*2^(3/5))/(4*2^(3/5) + 5)
[ (2*2^(3/5))/(4*2^(3/5) + 5), 4/(4*2^(3/5) + 5), (2*2^(3/5))/(4*2^(3/5) + 5), 1/(4*2^(3/5) + 5)]
1/13
a=[1 0 1 0]; b=[0 4 0 1];
syms x real
c = a*x + b; d= c/sum(c);
E = d * log(d)'
V = 0; for k =1:4; V= d(k)*log(k)+V; end; V
solve(V)
-18729944304496077/4947709893870346
log(10)
V = 0; for k =1:4; V= d(k)*log(d(k))+V; end; V
V*(2*x + 5)
simplify(V*(2*x + 5))
zz= simplify(V*(2*x + 5))
s=solve(zz)
s.x
s=nsolve(zz)
s=vpasolve(zz)
for t=-10:0.2:10; x=t;plot(x,eval(zz));hold on; end
for t=0:0.2:10; x=t;plot(x,eval(zz));hold on; end
for t=0:0.2:10; x=t;plot(x,eval(zz).'r*');hold on; end
for t=0:0.2:10; x=t;plot(x,eval(zz),'r*');hold on; end
log(256)
for t=0:0.2:10; x=t;plot(x,eval(zz),'r*');hold on; end
for t=0:0.2:10; x=t;plot(x,eval(V),'r*');hold on; end
for t=0:0.02:10; x=t;plot(x,eval(V),'r.');hold on; end
for t=1:0.02:5; x=t;plot(x,eval(V),'r.');hold on; end
for t=2.5:0.02:3.5; x=t;plot(x,eval(V),'r.');hold on; end
eval(d)
sum(d)
x=4;
rval(d)
eval(d)
eval(V)
for t=1:0.02:5; x=t;plot(x,eval(V),'r.');hold on; end
V
latex(V)
a=[1 0 1 0]; b=[-2.5 4 -2.5 1];
syms x real
c = a*x + b; d= c/sum(c);
V = 0; for k =1:4; V= d(k)*log(d(k))+V; end; V
for t=2.5:0.02:3.5; x=t;plot(x,eval(V),'r.');hold on; end
for t=2:0.02:5; x=t;plot(x,eval(V),'r.');hold on; end
for t=2:0.02:10; x=t;plot(x,eval(V),'r.');hold on; end
for t=2:0.02:10; x=t;plot(x,eval(V),'r.');cc=[cc; x eval(V)];hold on; end
cc=[];for t=2:0.02:10; x=t;plot(x,eval(V),'r.');cc=[cc; x eval(V)];hold on; end
x=5.54; eval(d)
b=[0.06920.24140.06920.9655]
b = [0.06920 .24140 .06920 .9655]
syms x real
c = a*x + b; d= c/sum(c);
V = 0; for k =1:4; V= d(k)*log(d(k))+V; end; V
cc=[];for t=2:0.02:10; x=t;plot(x,eval(V),'r.');cc=[cc; x eval(V)];hold on; end
cc=[];for t=2:0.02:20; x=t;plot(x,eval(V),'r.');cc=[cc; x eval(V)];hold on; end
a
b
cc=[];for t=2:0.2:40; x=t;plot(x,eval(V),'r.');cc=[cc; x eval(V)];hold on; end
cc=[];for t=2:1:140; x=t;plot(x,eval(V),'r.');cc=[cc; x eval(V)];hold on; end
b=[-0.6631 0.3369 -0.6631 0.0842]
syms x real
c = a*x + b; d= c/sum(c);
V = 0; for k =1:4; V= d(k)*log(d(k))+V; end; V
cc=[];for t=2:1:140; x=t;plot(x,eval(V),'r.');cc=[cc; x eval(V)];hold on; end
d
sum(s)
sum(d)
simplify(sum(d))
%-- 2020/10/6 12:25 --%
payoff_matrix_1=payoff88(4,2)
payoff_matrix_2 = -payoff_matrix_1;
syms x1 x2 x3 x4 x5 x6 x7 x8 y1 y2 y3 y4 y5 y6 y7 y8 real
% 求二个群体的各个策略在给定点的收益
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2 y3 y4 y5 y6 y7 y8]'
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2 x3 x4 x5 x6 x7 x8]'
% 求各个群体分别的收益平均值
mean_U_1 = [x1 x2 x3 x4 x5 x6 x7 x8] * Payoff_vector_field_F_1
mean_U_2 = [y1 y2 y3 y4 y5 y6 y7 y8] * Payoff_vector_field_F_2
% 求二个群体的各个策略在给定点的演化速度
V_Eq_1 = [x1 x2 x3 x4 x5 x6 x7 x8]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_2 = [y1 y2 y3 y4 y5 y6 y7 y8]'.*(Payoff_vector_field_F_2 - mean_U_2);
V_Eq_0 = [V_Eq_1; V_Eq_2]
D_Eq_0 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ...
diff(V_Eq_0,'x5') diff(V_Eq_0,'x6') diff(V_Eq_0,'x7') diff(V_Eq_0,'x8') ...
diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3') diff(V_Eq_0,'y4') ...
diff(V_Eq_0,'y5') diff(V_Eq_0,'y6') diff(V_Eq_0,'y7') diff(V_Eq_0,'y8') ...
]
% 显示特征矩阵各列
% 绘制演化时间序列
v=V_Eq_0;
T=[];p=0.125*ones(16,1);for t=1:1000;
x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);
y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16);
p = p + eval(v)*0.01;T = [T;p'];
end
figure  ;  plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
figure  ;  plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
payoff_matrix_1=payoff88(4,2)
payoff_matrix_2 = -payoff_matrix_1;
syms x1 x2 x3 x4 x5 x6 x7 x8 y1 y2 y3 y4 y5 y6 y7 y8 real
% 求二个群体的各个策略在给定点的收益
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2 y3 y4 y5 y6 y7 y8]'
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2 x3 x4 x5 x6 x7 x8]'
% 求各个群体分别的收益平均值
mean_U_1 = [x1 x2 x3 x4 x5 x6 x7 x8] * Payoff_vector_field_F_1
mean_U_2 = [y1 y2 y3 y4 y5 y6 y7 y8] * Payoff_vector_field_F_2
% 求二个群体的各个策略在给定点的演化速度
V_Eq_1 = [x1 x2 x3 x4 x5 x6 x7 x8]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_2 = [y1 y2 y3 y4 y5 y6 y7 y8]'.*(Payoff_vector_field_F_2 - mean_U_2);
V_Eq_0 = [V_Eq_1; V_Eq_2]
% D_Eq_0 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ...
% diff(V_Eq_0,'x5') diff(V_Eq_0,'x6') diff(V_Eq_0,'x7') diff(V_Eq_0,'x8') ...
% diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3') diff(V_Eq_0,'y4') ...
% diff(V_Eq_0,'y5') diff(V_Eq_0,'y6') diff(V_Eq_0,'y7') diff(V_Eq_0,'y8') ...
% ]
% 显示特征矩阵各列
% 绘制演化时间序列
v=V_Eq_0;
T=[];p=0.125*ones(16,1);
for t=1:1000;
x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);
y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16);
p = p + eval(v)*0.01;T = [T;p'];
end
figure  ;  plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3)
figure  ;  plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3)
grod on
grid on
ylim([-0.1 1])
ylim([-0.1 1.1])
saveas(gcf,'F:\8x8code\forZhengjie2020\X3.png','png')
saveas(gcf,'F:\8x8code\forZhengjie2020\Y3.png','png')
payoff_matrix_1=payoff88(3,2)
payoff_matrix_2 = -payoff_matrix_1;
syms x1 x2 x3 x4 x5 x6 x7 x8 y1 y2 y3 y4 y5 y6 y7 y8 real
% 求二个群体的各个策略在给定点的收益
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2 y3 y4 y5 y6 y7 y8]'
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2 x3 x4 x5 x6 x7 x8]'
% 求各个群体分别的收益平均值
mean_U_1 = [x1 x2 x3 x4 x5 x6 x7 x8] * Payoff_vector_field_F_1
mean_U_2 = [y1 y2 y3 y4 y5 y6 y7 y8] * Payoff_vector_field_F_2
% 求二个群体的各个策略在给定点的演化速度
V_Eq_1 = [x1 x2 x3 x4 x5 x6 x7 x8]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_2 = [y1 y2 y3 y4 y5 y6 y7 y8]'.*(Payoff_vector_field_F_2 - mean_U_2);
V_Eq_0 = [V_Eq_1; V_Eq_2]
% 绘制演化时间序列
v=V_Eq_0;
T=[];p=0.125*ones(16,1);
for t=1:1000;
x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);
y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16);
p = p + eval(v)*0.01;T = [T;p'];
end
figure  ;  plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3);ylim([-0.1 1.1]);grid on;saveas(gcf,'F:\8x8code\forZhengjie2020\X2.png','png')
figure  ;  plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3);ylim([-0.1 1.1]);grid on;saveas(gcf,'F:\8x8code\forZhengjie2020\Y2.png','png')
payoff_matrix_1=payoff88(2,1)
payoff_matrix_2 = -payoff_matrix_1;
syms x1 x2 x3 x4 x5 x6 x7 x8 y1 y2 y3 y4 y5 y6 y7 y8 real
% 求二个群体的各个策略在给定点的收益
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2 y3 y4 y5 y6 y7 y8]'
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2 x3 x4 x5 x6 x7 x8]'
% 求各个群体分别的收益平均值
mean_U_1 = [x1 x2 x3 x4 x5 x6 x7 x8] * Payoff_vector_field_F_1
mean_U_2 = [y1 y2 y3 y4 y5 y6 y7 y8] * Payoff_vector_field_F_2
% 求二个群体的各个策略在给定点的演化速度
V_Eq_1 = [x1 x2 x3 x4 x5 x6 x7 x8]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_2 = [y1 y2 y3 y4 y5 y6 y7 y8]'.*(Payoff_vector_field_F_2 - mean_U_2);
V_Eq_0 = [V_Eq_1; V_Eq_2]
% 绘制演化时间序列
v=V_Eq_0;
T=[];p=0.125*ones(16,1);
for t=1:1000;
x1=p(1);	x2=p(2);	x3=p(3);	x4=p(4);	x5=p(5);	x6=p(6);	x7=p(7);	x8=p(8);
y1=p(9);	y2=p(10);	y3=p(11);	y4=p(12);	y5=p(13);	y6=p(14);	y7=p(15);	y8=p(16);
p = p + eval(v)*0.01;T = [T;p'];
end
figure  ;  plot(T(:,1:8),'DisplayName','T(:,1:8)','LineWidth',3);ylim([-0.1 1.1]);grid on;saveas(gcf,'F:\8x8code\forZhengjie2020\X1.png','png')
figure  ;  plot(T(:,9:end),'DisplayName','T(:,9:end)','LineWidth',3);ylim([-0.1 1.1]);grid on;saveas(gcf,'F:\8x8code\forZhengjie2020\Y1.png','png')
saveas(gcf,'F:\8x8code\forZhengjie2020\X1.png','png')
%-- 2020/10/7 14:57 --%
?a =?[0.1381 + 0.3222i
0.4603 + 0.0000i
0.1381-0.3222i]
b=a*i;
cross(a,b)
a=?[0.1381 + 0.3222i
0.4603 + 0.0000i
0.1381-0.3222i]
b=a*i;
0.1381 + 0.3222i
0.4603 + 0.0000i
0.1381-0.3222i
a=[0.1381 + 0.3222i
0.4603 + 0.0000i
0.1381-0.3222i]
b=a*i;
cross(a,b)
b
a=[0.4603-0.0000i;0.1381 + 0.3222i ;-0.7365 + 0.0000i]
cross(a,a*i)
a=[1 -i]
cross(a,a*i)
a=[1 -i 0]
cross(a,a*i)
a*i
syms x
c=int(e^x)
c=int(exp(x))
c=int(exp(x)*exp(x))
syms x lambda
c=int(exp(lambda*x)*exp(lambda*x))
syms x lambda t
c=int(exp(lambda*t)*exp(lambda*t))
a=[1 exp(-pi/2) 0]; b=a*i; cross(a,b)
a=[1 exp(-pi/2) 0]; b=a*i; cross(real(a),real(b))
real(b)
b
exp(-pi/2)
a=[1 exp(-pi/2*i) 0]; b=a*i; cross(real(a),real(b))
a=[1 exp(-pi/2*i) 0]; b=a*i; cross(a,b)
cos(x)*sin(x)
simplify(cos(x)*sin(x))
simplify(sin(x)*sin(x))
int(sin(x)^2)
int(sin(x)*cos(x+t))
a=[0.1381 + 0.3222i
0.4603 + 0.0000i
0.1381-0.3222i]
syms x lambda
int(cross(real(a*exp(x*0.2i)),real(a*exp(x*0.2i))))
int(cross(real(a*exp(x*0.2i)),real(i*a*exp(x*0.2i))))
integral(cross(real(a*exp(x*0.2i)),real(i*a*exp(x*0.2i))))
integral(cross(real(a*exp(x*0.2i)),real(i*a*exp(x*0.2i))),0,10000)
cross(real(a*exp(x*0.2i)),real(i*a*exp(x*0.2i)))
fun=cross(real(a*exp(x.*0.2i)),real(i*a*exp(x.*0.2i)))
q=vap(int(fun,x,0,10000))
q=vpa(int(fun,x,0,10000))
a=[1;-1/2+sqrt(3)/2;-1/2-sqrt(3)/2;]
a=[1;-1/2+sqrt(3)/2*i;-1/2-sqrt(3)/2*i;]
fun=cross(real(a*exp(x.*0.2i)),real(i*a*exp(x.*0.2i)))
q=vpa(int(fun,x,0,10000))
lambda = 0.2i; fun=cross(real(a*exp(x.*lambda)),real(lambda*a*exp(x.*lambda)))
q=vpa(int(fun,x,0,10000))
a=[1;-i;-1/2-sqrt(3)/2*i;]
q=vpa(int(fun,x,0,10000))
a=[1;-i;-1/2-sqrt(3)/2*i;]
fun=cross(real(a*exp(x.*0.2i)),real(i*a*exp(x.*0.2i)))
q=vpa(int(fun,x,0,10000))
a=[1;-i;-1/2-sqrt(3)/2*i;]
fun=cross(real(a*exp(x.*0.2i)),real(i*a*exp(x.*0.2i)))
q=vpa(int(fun,x,0,10000))
lambda = 0.2i; fun=cross(real(a*exp(x.*lambda)),real(lambda*a*exp(x.*lambda)))
q=vpa(int(fun,x,0,10000))
clear
syms x lambda
a=[1;-i;-1/2-sqrt(3)/2*i;]
lambda = 0.2i; fun=cross(real(a*exp(x.*lambda)),real(lambda*a*exp(x.*lambda)))
q=vpa(int(fun,x,0,10000))
%-- 2020/10/12 20:47 --%
payoff_matrix_1 = [1 -1; -1 1]
payoff_matrix_2 = -payoff_matrix_1;
syms x1 x2 x3 x4 x5 x6 x7 x8 y1 y2 y3 y4 y5 y6 y7 y8 real
% 求二个群体的各个策略在给定点的收益
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2]'
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2]'
% 求各个群体分别的收益平均值
mean_U_1 = [x1 x2] * Payoff_vector_field_F_1
mean_U_2 = [y1 y2] * Payoff_vector_field_F_2
% 求二个群体的各个策略在给定点的演化速度
V_Eq_1 = [x1 x2]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_2 = [y1 y2]'.*(Payoff_vector_field_F_2 - mean_U_2);
V_Eq_0 = [V_Eq_1; V_Eq_2]
t1 = subs(V_Eq_0,'y2','1-y1')
t1 = subs(t1,'x2','1-x1')
simplify(t1)
L=[]; for x1=0:0.125:1; for y1=0:0.125:1; L=[L; eval([x1 y1 vpa(t1(1)) vpa(t1(3))])];end;end
quiver(L(:,1),L(:,2),L(:,3),L(:,4),1)
t2 = [(t1(1)); (t1(3))]
t3 = [diff(t2,'x1') diff(t2,'y1')  ] * t2
L=[]; for x1=0:0.125/2:1; for y1=0:0.125/2:1; L=[L; eval([x1 y1 vpa(t3(1)) vpa(t3(2))])];end;end
quiver(L(:,1),L(:,2),L(:,3),L(:,4),1)
clear
acceleration5005
figure
acceleration5005
figure
acceleration5005
%-- 2020/10/19 15:35 --%
accelerationRPS
L=[]; for x1=0:0.125/4:1; for x2=0:0.125/4:1; L=[L; eval([x1 x2 vpa(t1(1)) vpa(t1(2))])];end;end
quiver(L(:,1),L(:,2),L(:,3),L(:,4),1);hold on
L=[]; for x1=0:0.125:1; for x2=0:0.125:1; L=[L; eval([x1 x2 vpa(t1(1)) vpa(t1(2))])];end;end
quiver(L(:,1),L(:,2),L(:,3),L(:,4),1);hold on
L=[]; for x1=0:0.125:1; for x2=0:0.125:1-x1; L=[L; eval([x1 x2 vpa(t1(1)) vpa(t1(2))])];end;end
quiver(L(:,1),L(:,2),L(:,3),L(:,4),1);hold on
figure
t2 = [(t1(1)); (t1(2))]
t3 = [diff(t2,'x1') diff(t2,'x2')] * t2
L=[]; for x1=0:0.125:1; for x2=0:0.125:1-x1; L=[L; eval([x1 x2 vpa(t3(1)) vpa(t3(2))])];end;end
quiver(L(:,1),L(:,2),L(:,3),L(:,4),1,'red')
figure
t2 = [(t1(1)); (t1(2))]
t3 = [diff(t2,'x1') diff(t2,'x2')] * t2
L=[]; for x1=0:0.125/4:1; for x2=0:0.125/4:1-x1; L=[L; eval([x1 x2 vpa(t3(1)) vpa(t3(2))])];end;end
quiver(L(:,1),L(:,2),L(:,3),L(:,4),1,'red')
accelerationRPS
%-- 2020/11/12 17:21 --%
A=[12	-4.81E-17	0	10.8	-2.0278	6.4
13	-4.81E-17	0	-17.6	0.9556	-34.4
14	-4.81E-17	0	6.8	1.0722	28
16	0.3927	0	29.8	6.75	206.8
17	0.3927	0	16.4	15.1167	166.6
18	0.3927	0	50.8	6.8944	191
23	0	-0.036	6.8	-1.5222	16.6
24	0	0.036	4	-0.5056	-10.2
25	0.3927	0	30.2	20.2722	183.4
26	-0.1309	-0.0774	-8	-0.75	-46.8
27	-0.1309	0.2733	-18.2	-11.7333	-97.4
28	-0.1309	-0.1958	-4	-7.7889	-39.2
34	0	-0.0361	-10.8	-0.5667	-17.8
35	0.3927	0	47.6	3.5667	191.2
36	-0.1309	0.2733	-1.8	-2.3111	-41.4
37	-0.1309	-0.9813	-3	1.1222	-16
38	-0.1309	0.7079	-42.8	-2.3778	-133.8
45	0.3927	0	19.2	4.9222	189.8
46	-0.1309	-0.1958	-20	-3.6889	-118.6
47	-0.1309	0.7079	4.8	-4.5056	-53.2
48	-0.1309	-0.512	-4	3.2722	-18
56	-4.81E-17	0	-7.2	-3.4389	-4
57	-4.81E-17	0	30.4	3.2889	21.8
58	-4.81E-17	0	-23.2	0.15	-17.8
67	0	-0.036	-31.2	-4.9556	-12.6
68	0	0.036	24	1.5167	8.6
78	0	-0.0361	-0.8	-1.6667	9.2
];
plotmatrix(A(:,4:end))
grid on
plotmatrix(A(:,[2 4:end]))
plotmatrix(A(:,4:end))
scatter(getcolumn(unnamed(:,[10,end]),1),getcolumn(unnamed(:,[10,end]),2))
grid on
scatter((b(:,4), b(:,12))
scatter((b(:,4),b(:,12))
scatter(b(:,4),b(:,12))
scatter(b(:,4),b(:,12));grid on
text(b(:,4),b(:,12),num2str(b(:,1)));grid on
scatter(b(:,4),b(:,12));grid on
text(b(:,4),b(:,12),num2str(b(:,1)));grid on
figure
text(b(:,4),b(:,11),num2str(b(:,1)));grid on
scatter(b(:,4),b(:,11));grid on
text(b(:,4),b(:,11),num2str(b(:,1)));grid on
figure
scatter(b(:,4),b(:,10));grid on
text(b(:,4),b(:,10),num2str(b(:,1)));grid on
%-- 2020/11/13 11:14 --%
scatter(b(:,4),b(:,12));grid on
text(b(:,4),b(:,12),num2str(b(:,1)));grid on
figure
scatter(b(:,4),b(:,10));grid on
text(b(:,4),b(:,10),num2str(b(:,1)));grid on
c=b
b=c(find(c(:,15)==150),:)
scatter(b(:,4),b(:,12));grid on
text(b(:,4),b(:,12),num2str(b(:,1)));grid on
figure
scatter(b(:,4),b(:,10));grid on
text(b(:,4),b(:,10),num2str(b(:,1)));grid on
b=c;
c=b
b=c(find(c(:,15)==150),:)
scatter(b(:,4),b(:,12));grid on
text(b(:,4),b(:,12),num2str(b(:,1)));grid on
figure
scatter(b(:,4),b(:,10));grid on
text(b(:,4),b(:,10),num2str(b(:,1)));grid on
close all
b=c
scatter(b(:,4),b(:,12));grid on
text(b(:,4),b(:,12),num2str(b(:,1)));grid on
figure
scatter(b(:,4),b(:,10));grid on
text(b(:,4),b(:,10),num2str(b(:,1)));grid on
%-- 2020/11/16 0:16 --%
y=[7613.51  7850.91  8381.86  9142.81 10813.6 8631.43 8124.94 9429.79 10230.81 10163.61 9737.56 8561.06 7781.82 7110.97];
x1=[7666 7704 8148 8571 8679 7704 6471 5870 5289 3815 3335 2927 2758 2591];
x2=[16.22 16.85 17.93 17.28 17.23 17 19 18.22 16.3 13.37 11.62 10.36 9.83 9.25];
X = [ones(size(y)) x1.^2 x2.^2 x1 x2 x1.*x2];
[b,bint,r,rint,stats] = regress(y,X);
x1.^2
ones(size(y))
y=[7613.51  7850.91  8381.86  9142.81 10813.6 8631.43 8124.94 9429.79 10230.81 10163.61 9737.56 8561.06 7781.82 7110.97]
y=[7613.51  7850.91  8381.86  9142.81 10813.6 8631.43 8124.94 9429.79 10230.81 10163.61 9737.56 8561.06 7781.82 7110.97]';
x1=[7666 7704 8148 8571 8679 7704 6471 5870 5289 3815 3335 2927 2758 2591]';
x2=[16.22 16.85 17.93 17.28 17.23 17 19 18.22 16.3 13.37 11.62 10.36 9.83 9.25]';
X = [ones(size(y)) x1.^2 x2.^2 x1 x2 x1.*x2];
[b,bint,r,rint,stats] = regress(y,X);
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);
bb= [b bint]
y=data(:,8); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);
bb= [b bint]
bb= [b bint stats]
bb= [b bint stats']
r = []; for k=8:12; y=data(:,k); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']; r=[r; bb];end
rr = []; for k=8:12; y=data(:,k); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']; rr=[rr; bb];end
data(:,13)=data(:,10)+data(:,11)+data(:,12)+data(:,8)
data(:,13)=data(:,10)+data(:,11)+data(:,12)-data(:,8)
y=data(:,13); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']
-data(:,8)
data(:,14)=data(:,10)+data(:,11)
y=data(:,14); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']
y=data(:,13); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']
y=data(:,14); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2  ];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint ]
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x3 ];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint ]
y=data(:,8); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x3 ];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint ]
y=data(:,8); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 ];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint ]
rr = []; for k=8:12; y=data(:,k); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']; rr=[rr; bb];end
corr(data(:,8:12))
%-- 2020/11/16 19:23 --%
MultiLRYQM
y=data(:,15); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']
bb(3,4)
rint
fi = 'F:\cDownload\Binmore_2.csv'
angY = in_8colExp_out_am(fi)
[aa1 aa2 ] =  in_8colExp_out_am(fi)
fi = 'E:\AbedData\01\OGame1.csv';[aa1 aa2 ] =  in_8colExp_out_am(fi)
y=data(:,15); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);
scatter(data(:,15),data(:,6),'ro')
scatter(data(:,6),data(:,15),'ro')
scatter(data(:,14),data(:,15),'ro')
scatter(data(:,12),data(:,15),'ro')
sum(data)
y=data(:,15); x1=data(:,6);   X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X)
y=data(:,15); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);
fi = 'E:\AbedData\01\OGame2.csv';[aa1 aa2 ] =  in_8colExp_out_am(fi)
y=data(:,16); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);
y=data(:,15); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X)
scatter(getcolumn(data(:,[6,end]),1),getcolumn(data(:,[6,end]),2))
a(:,17) = a(:,2) + a(:.4)
a(:,17) = a(:,2) + a(:,4)
data(:,17) = data(:,2) + data(:,4)
y=data(:,16); x1=data(:,6); x2=data(:,17) X = [ones(size(y)) x1 x2 ];[b,bint,r,rint,stats] = regress(y,X);
y=data(:,16); x1=data(:,6); x2=data(:,17); X = [ones(size(y)) x1 x2 ];[b,bint,r,rint,stats] = regress(y,X);
y=data(:,16); x1=data(:,6); x2=data(:,17); X = [ones(size(y)) x1 x2 ];[b,bint,r,rint,stats] = regress(y,X)
%-- 2020/11/16 23:35 --%
fi = 'E:\AbedData\01\OGame1.csv';[aa1 aa2 ] =  in_8colExp_out_am(fi)
data =[26	0.900211091	-0.900211091	0.077437117	-0.0774	-0.1309	0.130868728	-2.2	-0.749999997	10.6	-33.2	-16.2
27	-0.143086216	0.143086216	-0.27330414	0.2733	-0.1309	0.130868728	-3.8	-2.31111111	-30.2	13	-62
28	0.741983761	-0.741983761	0.195831955	-0.1958	-0.1309	0.130868728	-19	-3.688888886	-42.2	6	1
36	-0.143086216	0.143086216	-0.27330414	0.2733	-0.1309	0.130868728	-23.2	-11.73333332	-8.2	-6.2	-25.2
37	0.028412595	-0.028412595	0.981336815	-0.9813	-0.1309	0.130868728	-8.8	1.12222222	-13	-11	11
38	-0.131098343	0.131098343	-0.707913246	0.7079	-0.1309	0.130868728	4	-4.505555552	-24	-17	-50
46	-0.756919893	0.756919893	0.195831955	-0.1958	-0.1309	0.130868728	0.8	-7.788888883	-56.2	-17.2	-25.2
47	0.114632307	-0.114632307	-0.707913246	0.7079	-0.1309	0.130868728	-42.8	-2.377777774	-18	-21	-19
48	-0.610696192	0.610696192	0.511996947	-0.512	-0.1309	0.130868728	-4	3.27222222	22	-12	-24
]
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X);
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X)
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y))  x2 x3];[b,bint,r,rint,stats] = regress(y,X)
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1  x3];[b,bint,r,rint,stats] = regress(y,X);
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1  x3];[b,bint,r,rint,stats] = regress(y,X)
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x2  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x2  x3];[b,bint,r,rint,stats] = regress(y,X)
fi = 'E:\AbedData\01\OGame1.csv';[aa1 aa2 ] =  in_8colExp_out_am(fi)
y=data(:,18); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X)
y=data(:,14); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];[b,bint,r,rint,stats] = regress(y,X)
data2=data(find(abs(data(:,6)-0.1309)<0.01),:)
data2=data(find(abs(data(:,6)+0.1309)<0.01),:)
y=data2(:,12); x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x2  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y)) x2  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y)) z1 x2  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y)) x1 x2  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x2  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x2  x3];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x2 ];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x3 ];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,14); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x3 ];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,14); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x2 ];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,14); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x3];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,12); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x3];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,11); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x3];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,8); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x3];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,8); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x2];[b,bint,r,rint,stats] = regress(y,X)
y=data2(:,9); x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x2];[b,bint,r,rint,stats] = regress(y,X)
%-- 2020/11/17 16:37 --%
MultiLRYQM
y=data(:,10) + data(:,11) + data(:,12);
x1=data2(:,6); x2=data2(:,2);x3=data2(:,4); X = [ones(size(y))  x2];[b,bint,r,rint,stats] = regress(y,X)
x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x2  x2];[b,bint,r,rint,stats] = regress(y,X)
x1=data(:,6); x2=data(:,2);x3=data(:,4); X = [ones(size(y)) x1 x2  x3];[b,bint,r,rint,stats] = regress(y,X)
scatter(y,x1)
scatter(x1,y)
lsline(x1,y)
lsline
data(:,15) = y
for k=10:15;subplot(3,2,k-9);scatter(x1,data(:,k-9));lsline;end
x1=data(:,6);
for k=10:15;subplot(3,2,k-9);scatter(x1,data(:,k-9));lsline;end
for k=10:15;subplot(3,2,k-9);scatter(x1,data(:,k-9));end
for k=10:15;subplot(3,2,k-9);scatter(x1,data(:,k));end
for k=10:15;subplot(3,2,k-9);scatter(x1,data(:,k));lsline;end
polyfit(x1,y)
polyfit(x1,y,1)
[h p]=ttest([])
[h p]=ttest([10.6
-30.2
-42.2
-8.2
-13
-24
-56.2
-18
22
])
[h p]=ttest([])10.6
-30.2
-42.2
-8.2
-13
-24
-56.2
-18
22
[h p]=ttest([  10.6
-30.2
-42.2
-8.2
-13
-24
-56.2
-18
22
])
x=[8	0.75	10.6	-33.2	-16.2
18.2	11.7333	-30.2	13	-62
4	7.7889	-42.2	6	1
1.8	2.3111	-8.2	-6.2	-25.2
3	-1.1222	-13	-11	11
42.8	2.3778	-24	-17	-50
20	3.6889	-56.2	-17.2	-25.2
-4.8	4.5056	-18	-21	-19
4	-3.2722	22	-12	-24
]
[h p]=ttest(x(:,1))
[h p]=ttest(x(:,2))
[h p]=ttest(x(:,3))
[h p]=ttest(x(:,4))
[h p]=ttest(x(:,5))
[h p]=ttest(x(:,1:5))
[h p]=mean(x(:,1:5))
mean(x(:,1:5))
[h p]=mean(x(:,1);x(:,2))
[h p]=mean([x(:,1);x(:,2)])
[h p]=ttest([x(:,1);x(:,2)])
pp = [0.0544    0.0670    0.0615    0.0462    0.0146]
samp=[2650 1950 3000 3000 4700]
scatter(pp,samp)
data = [23	0	-0.071439963	0.036010289	29	4.955555547	-27.2	2	19.2	16	901	2	3
24	0	0.071549996	-0.036019622	-21	-1.516666671	3.8	-6	-15.2	3	901	2	4
26	-0.1309	0.900211091	0.077437117	-2.2	-0.749999997	10.6	-33.2	4	-20.2	911	2	6
27	-0.1309	-0.143086216	-0.27330414	-3.8	-2.31111111	-30.2	13	-28	-33.8	911	2	7
28	-0.1309	0.741983761	0.195831955	-19	-3.688888886	-42.2	6	-11.2	12	911	2	8
34	0	-0.071441185	0.036059535	-1.33227E-15	1.666666663	-17	1	-2.4	12	901	3	4
36	-0.1309	-0.143086216	-0.27330414	-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2	911	3	6
37	-0.1309	0.028412595	0.981336815	-8.8	1.12222222	-13	-11	5.8	5.2	911	3	7
38	-0.1309	-0.131098343	-0.707913246	4	-4.505555552	-24	-17	-25.4	-24	911	3	8
46	-0.1309	-0.756919893	0.195831955	0.8	-7.788888883	-56.2	-17.2	-7.8	-17.2	911	4	6
47	-0.1309	0.114632307	-0.707913246	-42.8	-2.377777774	-18	-21	-13.8	-5.8	911	4	7
48	-0.1309	-0.610696192	0.511996947	-4	3.27222222	22	-12	-2	-22	911	4	8
67	0	-0.071439963	0.036010289	-8.6	1.522222224	-9.8	1.2	11	16.4	901	6	7
68	0	0.16585302	-0.036019622	-7.8	0.505555554	6.2	3.2	-14.2	-10.8	901	6	8
78	0	0.032521274	0.036059535	14.8	0.56666667	1	15	-0.2	-5.2	901	7	8
]
rr = [];
for k=5:10; y=data(:,k);
x1=data(:,2); x2=data(:,3);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];
[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']; rr=[rr; bb];
end
rr = []; bb = []; for k=5:10; y=data(:,k);
x1=data(:,2); x2=data(:,3);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];
[b,bint,r,rint,stats] = regress(y,X);bb= [bb;b bint ]; rr=[rr; stats];
end
rr = []; bb = []; for k=5:10; y=data(:,k);
x1=data(:,2); x2=data(:,3);x3=data(:,4); X = [ones(size(y)) x1  x3];
[b,bint,r,rint,stats] = regress(y,X);bb= [bb;b bint ]; rr=[rr; stats];
end
rr = []; bb = []; for k=5:10; y=data(:,k);
x1=data(:,2); x2=data(:,3);x3=data(:,4); X = [ones(size(y)) x1];
[b,bint,r,rint,stats] = regress(y,X);bb= [bb;b bint ]; rr=[rr; stats];
end
data(:,14) = sum(data(:,5:10)')'
rr = []; bb = []; y=data(:,14);
x1=data(:,2); x2=data(:,3);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];
[b,bint,r,rint,stats] = regress(y,X);bb= [bb;b bint ]; rr=[rr; stats];
data(:,15) =  data(:,14) -  data(:,6)
rr = []; bb = []; y=data(:,15);
x1=data(:,2); x2=data(:,3);x3=data(:,4); X = [ones(size(y)) x1 x2 x3];
[b,bint,r,rint,stats] = regress(y,X);bb= [bb;b bint ]; rr=[rr; stats];
r=[];
for i1 = 0:6
for i2 = 0:6-i1
for i3=0:6-i1-i2
i4 = 6-i1-i2-i3
if i4 >= 0
r=[r;i1 i2 i3 i4]
end
end
end
end
84*84
r2=[];
for i1 = 0:6
for i2 = 0:6
for i3=0:6
i4 = 0:6
if i1+i2+i3+i4 == 6
r2=[r2;i1 i2 i3 i4]
end
end
end
end
r2=[];
for i1 = 0:6
for i2 = 0:6
for i3=0:6
for i4 = 0:6
if i1+i2+i3+i4 == 6
r2=[r2;i1 i2 i3 i4]
end
end
end
end
end
%-- 2020/11/18 16:02 --%
m=3;n=2;lissa(m,n)
m=3;n=2;lissa(m,n,1/4*pi)
m=1;n=1;lissa(m,n,1/4*pi)
lissa
x=1+i
rewrite(x,'cos')
lissa2Complex2areaPlot
vim=1;vjm=1+i;lambda_i=8i;lambda_j=4i;lissa2Complex2areaPlot(vim,vjn,lambda_i,lambda_j)
vim=1;vjn=1+i;lambda_i=8i;lambda_j=4i;lissa2Complex2areaPlot(vim,vjn,lambda_i,lambda_j)
real(vim)
image(vim)
Im(vim)
imag(vim)
vim=1;vjn=1+i;lambda_i=8i;lambda_j=4i;lissa2Complex2areaPlot(vim,vjn,lambda_i,lambda_j)
abs(vjn)
angle(vjn)*180/pi
vim=1;vjn=1+i;lambda_i=8i;lambda_j=4i;lissa2Complex2areaPlot(vim,vjn,lambda_i,lambda_j)
pi/4
vim=1;vjn=1-i;lambda_i=8i;lambda_j=4i;lissa2Complex2areaPlot(vim,vjn,lambda_i,lambda_j)
vi=1;vj=1-i;lambda_i=8i;lambda_j=4i;lissa2Complex2areaPlot(vi,vj,lambda_i,lambda_j)
YQM20200830
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
OneillYMQ20201118
v=eigen_vector(:,1)
ae=[];for m=1:7;for n=m:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j]end;end
ae=[];for m=1:7;for n=m:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end
ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end
diag(eigen_value)
diag(eigen_value)'
ag=round(af,4)
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end;evt=round(af,4);evl=round(diag(eigen_value)',4);
diag(eigen_value
diag(eigen_value)
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end;evt=round(af,4);evl=round(diag(eigen_value),4);
MultiLRYQM
e4=round(eigen_vector,3)
plot(e4(:,5))
spy(e4(:,5))
plot(e4(:,5))
plot(e4(:,7))
figure
plot(e4(:,5))
plot(e4(:,7))
plot(real(e4(:,5)),imag(e4(:,5)))
scatter(real(e4(:,5)),imag(e4(:,5)))
scatter(real(e4(:,5)),imag(e4(:,5)));axis([-0.5 0.5])
scatter(real(e4(:,5)),imag(e4(:,5)));axis([-0.5 0.5],[-0.5 0.5])
scatter(real(e4(:,5)),imag(e4(:,5)));axis([-0.5 0.5 -0.5 0.5])
figure
scatter(real(e4(:,7)),imag(e4(:,7)));axis([-0.5 0.5 -0.5 0.5])
grid on
yv=[  -0.0000 - 0.6124i  -0.0000 + 0.6124i  -0.7559 + 0.0000i  -0.0000 + 0.0000i  -0.0000 - 0.0000i  -0.0000 + 0.0000i  -0.0000 + 0.0000i  -0.0000 - 0.0000i
0.0000 + 0.2041i   0.0000 - 0.2041i  -0.3780 + 0.0000i   0.0000 + 0.0000i  -0.0000 + 0.5353i  -0.0000 - 0.5353i   0.1557 - 0.0205i   0.1557 + 0.0205i
0.0000 + 0.2041i   0.0000 - 0.2041i  -0.3780 + 0.0000i  -0.0000 + 0.0000i  -0.0425 - 0.0851i  -0.0425 + 0.0851i  -0.5589 + 0.0000i  -0.5589 + 0.0000i
0.0000 + 0.2041i   0.0000 - 0.2041i  -0.3780 + 0.0000i   0.0000 + 0.0000i   0.0425 - 0.4501i   0.0425 + 0.4501i   0.4032 + 0.0205i   0.4032 - 0.0205i
0.6124 + 0.0000i   0.6124 + 0.0000i  -0.0000 + 0.0000i   0.7559 + 0.0000i   0.0000 - 0.0000i   0.0000 + 0.0000i  -0.0000 + 0.0000i  -0.0000 - 0.0000i
-0.2041 - 0.0000i  -0.2041 + 0.0000i  -0.0000 + 0.0000i   0.3780 + 0.0000i   0.5353 + 0.0000i   0.5353 + 0.0000i  -0.0205 - 0.1557i  -0.0205 + 0.1557i
-0.2041 + 0.0000i  -0.2041 - 0.0000i  -0.0000 + 0.0000i   0.3780 + 0.0000i  -0.0851 + 0.0425i  -0.0851 - 0.0425i  -0.0000 + 0.5589i  -0.0000 - 0.5589i
-0.2041 + 0.0000i  -0.2041 - 0.0000i   0.0000 + 0.0000i   0.3780 + 0.0000i  -0.4501 - 0.0425i  -0.4501 + 0.0425i   0.0205 - 0.4032i   0.0205 + 0.4032i]
v4(:,5).*v4(:,7)
e4(:,5).*e4(:,7)
e4(:,5)*e4(:,7)
dot(e4(:,5),e4(:,7))
abs(e4(:,5))
angle(e4(:,5))
angle(e4(:,5))/pi
angle(e4(:,7))/pi
angle(e4(:,6))/pi
dot(e4(:,5),e4(:,7))
dot(e4(:,5),e4(:,7))*8
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end;evt=round(af,4);evl=round(diag(eigen_value),4);
dot(eigen_vector(:,5),eigen_vector(:,7))
cross(eigen_vector(:,5),eigen_vector(:,7))
e5 = [abs(eigen_vector(:,5)),angle(eigen_vector(:,7))]
e5 = [abs(eigen_vector(:,5)),angle(eigen_vector(:,5))]
e5 = [abs(eigen_vector(:,5)),angle(eigen_vector(:,5)) abs(eigen_vector(:,7)),angle(eigen_vector(:,7))]
%-- 2020/11/19 1:42 --%
YQM20200830
eqq = subs(D_Eq_0,'x4','1-x1-x2-x3')
OneillYMQ20201118
D_Eq_at_NE = eval(D_Eq_0)
OneillYMQ20201118
OneillYMQ20201118
OneillYMQ20201118
V_Eq_0
OneillYMQ20201118
D_Eq_0
sum(D_Eq_0(:,1:4))
simplify(sum(D_Eq_0(:,1:4)))
simplify(sum(V_Eq_0(:,1:4)))
simplify(sum(V_Eq_0))
OneillYMQ20201118
A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]);
x1=A(length(A(:,1))-2,1); x2=A(length(A(:,1))-2,2);
x3=A(length(A(:,1))-2,3); %x4=A(length(A(:,1))-2,4);
y1=A(length(A(:,1))-2,5); y2=A(length(A(:,1))-2,6);
y3=A(length(A(:,1))-2,7); %y4=A(length(A(:,1))-2,8);
% A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]);
A=eval([S.x1 S.x2 S.x3  S.y1 S.y2 S.y3 ]);
x1=A(length(A(:,1))-2,1); x2=A(length(A(:,1))-2,2);
x3=A(length(A(:,1))-2,3); %x4=A(length(A(:,1))-2,4);
y1=A(length(A(:,1))-2,5); y2=A(length(A(:,1))-2,6);
y3=A(length(A(:,1))-2,7); %y4=A(length(A(:,1))-2,8);
% D_Eq_at_NE = eval(D_Eq_0)
A=eval([S.x1 S.x2 S.x3  S.y1 S.y2 S.y3 ]);
x1=A(length(A(:,1))-2,1); x2=A(length(A(:,1))-2,2);
x3=A(length(A(:,1))-2,3); %x4=A(length(A(:,1))-2,4);
y1=A(length(A(:,1))-2,5-1); y2=A(length(A(:,1))-2,6-1);
y3=A(length(A(:,1))-2,7-1); %y4=A(length(A(:,1))-2,8);
OneillYMQ20201118
epp = subs(V_Eq_0, x4 , 1-x1-x2-x3 )
OneillYMQ20201118
A=eval([S.x1 S.x2 S.x3  S.y1 S.y2 S.y3 ]);
INs = -4
x1=A(INs,1); x2=A(INs,2);
x3=A(INs,3); %x4=A(length(A(:,1))-2,4);
y1=A(INs,5-1); y2=A(INs,6-1);
y3=A(INs,7-1); %y4=A(length(A(:,1))-2,8);
% D_Eq_at_NE = eval(D_Eq_0)
INs = length(A(:,1))-4
x1=A(INs,1); x2=A(INs,2);
x3=A(INs,3); %x4=A(length(A(:,1))-2,4);
y1=A(INs,5-1); y2=A(INs,6-1);
y3=A(INs,7-1); %y4=A(length(A(:,1))-2,8);
A6= D_Eq_0([1:3 5:7],[1:3 5:7])
D_Eq_at_NE = eval(A6)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
af=[]; for k=1:6;v=eigen_vector(:,k);ae=[];for m=1:5;for n=m+1:6;
vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];
end;end;af=[af ae];end;
evt=round(af,4); evl=round(diag(eigen_value),4);
MultiLRYQM
a=[-2.2	-0.749999997	10.6	-33.2	4	-20.2
-3.8	-2.31111111	-30.2	13	-28	-33.8
-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2
-8.8	1.12222222	-13	-11	5.8	5.2
]
>
[h p]=ttest(a)
[h p]=ttest([a(:1);a(:2);a(:3);a(:4);a(:5);a(:6)])
[h p]=ttest([a(:,1);a(:,2);a(:,3);a(:,4);a(:,5);a(:,6)])
mean(z)
mean(a)
b=[-4.8088E-17	0	0	0.003047619	0.001763533	-0.0078	-0.001333333	0.001687764	0.008016878	0.000897077
-4.8088E-17	0	0	-0.011047619	-0.00168661	0.0034	-0.000333333	-0.009113924	-0.001687764	-0.003411542
0.3927	0	0	0.009371429	0.010396011	0.017933333	0.018866667	0.007594937	0.020506329	0.014111451
0.3927	0	0	0.021104762	0.00182906	0.0204	0.006333333	0.015189873	0.014514768	0.013228633
0	-0.071439963	0.036010289	0.011047619	0.002541311	-0.009066667	0.000666667	0.008101266	0.006751055	0.003340208
0.3927	0	0	0.00952381	0.003461538	0.0206	0.004733333	0.014852321	0.017721519	0.01181542
-0.1309	0.900211091	0.077437117	-0.000838095	-0.000384615	0.003533333	-0.011066667	0.001687764	-0.008523207	-0.002598581
-0.1309	-0.143086216	-0.27330414	-0.001447619	-0.001185185	-0.010066667	0.004333333	-0.011814346	-0.014261603	-0.005740348
0.3927	0	0	0.010666667	0.007752137	0.015066667	0.0114	0.014261603	0.012658228	0.01196755
-0.1309	-0.143086216	-0.27330414	-0.008838095	-0.006017094	-0.002733333	-0.002066667	-0.005991561	-0.004725738	-0.005062081
-0.1309	0.028412595	0.981336815	-0.003352381	0.000575499	-0.004333333	-0.003666667	0.002447257	0.002194093	-0.001022589
-4.8088E-17	0	0	-0.006247619	0.001039886	-0.0012	0.001466667	-0.001350211	0.002362869	-0.000654735
-4.8088E-17	0	0	0.008914286	-0.000490028	0.0036	0.0046	-0.004725738	-0.009113924	0.000464099
0	-0.071439963	0.036010289	-0.00327619	0.000780627	-0.003266667	0.0004	0.00464135	0.006919831	0.001033159
-1.1782	0	0	-0.037714286	-0.014749288	-0.053066667	-0.032866667	-0.03907173	-0.049367089	-0.037805954
]
x1=b(:,1); x2=b(:, 2);x3=b(:,3);y=b(:,10); %wi8i4
X = [ones(size(y)) x1 x2 x3];
[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']
cc=[1	2	0	0	0
1	3	0	0	0
1	6	0.4284	0	0
1	7	0.4284	0	0
2	3	0	-0.172	-0.3518
2	5	0.4284	0	0
2	6	-0.1428	0.6841	0.5758
2	7	-0.1428	0.7596	-0.6702
3	5	0.4284	0	0
3	6	-0.1428	0.7596	-0.6702
3	7	-0.1428	0.8867	0.995
5	6	0	0	0
5	7	0	0	0
6	7	0	-0.172	-0.3518
1	5	-1.1781	0	0
]
b=[-4.8088E-17	0	0	0.003047619	0.001763533	-0.0078	-0.001333333	0.001687764	0.008016878	0.000897077
-4.8088E-17	0	0	-0.011047619	-0.00168661	0.0034	-0.000333333	-0.009113924	-0.001687764	-0.003411542
0.3927	0	0	0.009371429	0.010396011	0.017933333	0.018866667	0.007594937	0.020506329	0.014111451
0.3927	0	0	0.021104762	0.00182906	0.0204	0.006333333	0.015189873	0.014514768	0.013228633
0	-0.071439963	0.036010289	0.011047619	0.002541311	-0.009066667	0.000666667	0.008101266	0.006751055	0.003340208
0.3927	0	0	0.00952381	0.003461538	0.0206	0.004733333	0.014852321	0.017721519	0.01181542
-0.1309	0.900211091	0.077437117	-0.000838095	-0.000384615	0.003533333	-0.011066667	0.001687764	-0.008523207	-0.002598581
-0.1309	-0.143086216	-0.27330414	-0.001447619	-0.001185185	-0.010066667	0.004333333	-0.011814346	-0.014261603	-0.005740348
0.3927	0	0	0.010666667	0.007752137	0.015066667	0.0114	0.014261603	0.012658228	0.01196755
-0.1309	-0.143086216	-0.27330414	-0.008838095	-0.006017094	-0.002733333	-0.002066667	-0.005991561	-0.004725738	-0.005062081
-0.1309	0.028412595	0.981336815	-0.003352381	0.000575499	-0.004333333	-0.003666667	0.002447257	0.002194093	-0.001022589
-4.8088E-17	0	0	-0.006247619	0.001039886	-0.0012	0.001466667	-0.001350211	0.002362869	-0.000654735
-4.8088E-17	0	0	0.008914286	-0.000490028	0.0036	0.0046	-0.004725738	-0.009113924	0.000464099
0	-0.071439963	0.036010289	-0.00327619	0.000780627	-0.003266667	0.0004	0.00464135	0.006919831	0.001033159
-1.1782	0	0	-0.037714286	-0.014749288	-0.053066667	-0.032866667	-0.03907173	-0.049367089	-0.037805954
]
x1=cc(:,3); x2=cc(:, 4);x3=cc(:,5);y=b(:,10); %wi8i4
X = [ones(size(y)) x1 x2 x3];
[bxx,bint,r,rint,stats] = regress(y,X);bb= [bxx bint stats']
%-- 2020/11/19 12:45 --%
%一、replicator dynamics
syms x1 x2 x3 x4 y1 y2 y3 y4 real
payoff_matrix_1 = [1 -1 -1 -1; -1 -1 1 1; -1 1 -1 1; -1 1 1 -1]
payoff_matrix_2 = -payoff_matrix_1
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2 y3 y4]'
Payoff_vector_field_F_2 = payoff_matrix_2' *[x1 x2 x3 x4]'
mean_U_1 = [x1 x2 x3 x4] * Payoff_vector_field_F_1
mean_U_2 = [y1 y2 y3 y4] * Payoff_vector_field_F_2
V_Eq_1 = [x1 x2 x3 x4]'.*(Payoff_vector_field_F_1 - mean_U_1);
V_Eq_2 = [y1 y2 y3 y4]'.*(Payoff_vector_field_F_2 - mean_U_2);
V_Eq_0 = [V_Eq_1; V_Eq_2]
D_Eq_0 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ...
diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3') diff(V_Eq_0,'y4') ...
]
for i = 1:5
'column-i of the charactor matrix)'
D_Eq_0(:,i)
end
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]
%下面将纳什均衡填入
A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]);
x1=A(length(A(:,1)),1); x2=A(length(A(:,1)),2);
x3=A(length(A(:,1)),3); x4=A(length(A(:,1)),4);
y1=A(length(A(:,1)),5); y2=A(length(A(:,1)),6);
y3=A(length(A(:,1)),7); y4=A(length(A(:,1)),8);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
eigen_vector
eigen_vector = round(eigen_vector,4)
eigen_value
payoff_matrix_2
D_Eq_at_NE + rand(0.001)
rand(8)
rand(8)*0.0001
D_Eq_at_NE + rand(8)*0.0001
[v d] = eig(D_Eq_at_NE + rand(8)*0.0001)
[v d] = eig(D_Eq_at_NE + rand(8)*0.0000)
d(5,5)-d(7,7)
[v d] = eig(D_Eq_at_NE + rand(8)*0.0000)
d(5,5)-d(7,7)
[v d] = eig(D_Eq_at_NE + rand(8)*0.0001)
d(5,5)-d(7,7)
%-- 2020/11/19 15:14 --%
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end;evt=round(af,4);evl=round(diag(eigen_value)',4);
-0.0000 - 0.6124i	  -0.0000 + 0.6124i	  -0.7559 + 0.0000i	  -0.0000 + 0.0000i	  -0.0000 - 0.0000i	  -0.0000 + 0.0000i	  -0.0000 + 0.0000i	  -0.0000 - 0.0000i
0.0000 + 0.2041i	   0.0000 - 0.2041i	  -0.3780 + 0.0000i	   0.0000 + 0.0000i	  -0.0000 + 0.5353i	  -0.0000 - 0.5353i	   0.1557 - 0.0205i	   0.1557 + 0.0205i
0.0000 + 0.2041i	   0.0000 - 0.2041i	  -0.3780 + 0.0000i	  -0.0000 + 0.0000i	  -0.0425 - 0.0851i	  -0.0425 + 0.0851i	  -0.5589 + 0.0000i	  -0.5589 + 0.0000i
0.0000 + 0.2041i	   0.0000 - 0.2041i	  -0.3780 + 0.0000i	   0.0000 + 0.0000i	   0.0425 - 0.4501i	   0.0425 + 0.4501i	   0.4032 + 0.0205i	   0.4032 - 0.0205i
0.6124 + 0.0000i	   0.6124 + 0.0000i	  -0.0000 + 0.0000i	   0.7559 + 0.0000i	   0.0000 - 0.0000i	   0.0000 + 0.0000i	  -0.0000 + 0.0000i	  -0.0000 - 0.0000i
-0.2041 - 0.0000i	  -0.2041 + 0.0000i	  -0.0000 + 0.0000i	   0.3780 + 0.0000i	   0.5353 + 0.0000i	   0.5353 + 0.0000i	  -0.0205 - 0.1557i	  -0.0205 + 0.1557i
-0.2041 + 0.0000i	  -0.2041 - 0.0000i	  -0.0000 + 0.0000i	   0.3780 + 0.0000i	  -0.0851 + 0.0425i	  -0.0851 - 0.0425i	  -0.0000 + 0.5589i	  -0.0000 - 0.5589i
-0.2041 + 0.0000i	  -0.2041 - 0.0000i	   0.0000 + 0.0000i	   0.3780 + 0.0000i	  -0.4501 - 0.0425i	  -0.4501 + 0.0425i	   0.0205 - 0.4032i	   0.0205 + 0.4032i
gen_vector =	[
-0.0000 - 0.6124i	  -0.0000 + 0.6124i	  -0.7559 + 0.0000i	  -0.0000 + 0.0000i	  -0.0000 - 0.0000i	  -0.0000 + 0.0000i	  -0.0000 + 0.0000i	  -0.0000 - 0.0000i
0.0000 + 0.2041i	   0.0000 - 0.2041i	  -0.3780 + 0.0000i	   0.0000 + 0.0000i	  -0.0000 + 0.5353i	  -0.0000 - 0.5353i	   0.1557 - 0.0205i	   0.1557 + 0.0205i
0.0000 + 0.2041i	   0.0000 - 0.2041i	  -0.3780 + 0.0000i	  -0.0000 + 0.0000i	  -0.0425 - 0.0851i	  -0.0425 + 0.0851i	  -0.5589 + 0.0000i	  -0.5589 + 0.0000i
0.0000 + 0.2041i	   0.0000 - 0.2041i	  -0.3780 + 0.0000i	   0.0000 + 0.0000i	   0.0425 - 0.4501i	   0.0425 + 0.4501i	   0.4032 + 0.0205i	   0.4032 - 0.0205i
0.6124 + 0.0000i	   0.6124 + 0.0000i	  -0.0000 + 0.0000i	   0.7559 + 0.0000i	   0.0000 - 0.0000i	   0.0000 + 0.0000i	  -0.0000 + 0.0000i	  -0.0000 - 0.0000i
-0.2041 - 0.0000i	  -0.2041 + 0.0000i	  -0.0000 + 0.0000i	   0.3780 + 0.0000i	   0.5353 + 0.0000i	   0.5353 + 0.0000i	  -0.0205 - 0.1557i	  -0.0205 + 0.1557i
-0.2041 + 0.0000i	  -0.2041 - 0.0000i	  -0.0000 + 0.0000i	   0.3780 + 0.0000i	  -0.0851 + 0.0425i	  -0.0851 - 0.0425i	  -0.0000 + 0.5589i	  -0.0000 - 0.5589i
-0.2041 + 0.0000i	  -0.2041 - 0.0000i	   0.0000 + 0.0000i	   0.3780 + 0.0000i	  -0.4501 - 0.0425i	  -0.4501 + 0.0425i	   0.0205 - 0.4032i	   0.0205 + 0.4032i];
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end;evt=round(af,4);evl=round(diag(eigen_value)',4);
eigen_vector = gen_vector
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end;evt=round(af,4);evl=round(diag(eigen_value)',4);
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8; vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];end;end;af=[af ae];end;evt=round(af,4);
ev = [[-0.00000000000000 - 0.612400000000000i,0.00000000000000 + 0.612400000000000i,-0.755900000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.535300000000000i,-0.00000000000000 - 0.535300000000000i,0.155700000000000 - 0.0205000000000000i,0.155700000000000 + 0.0205000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,-0.0425000000000000 - 0.0851000000000000i,-0.0425000000000000 + 0.0851000000000000i,-0.558900000000000 + 0.00000000000000i,-0.558900000000000 + 0.00000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.0425000000000000 - 0.450100000000000i,0.0425000000000000 + 0.450100000000000i,0.403200000000000 + 0.0205000000000000i,0.403200000000000 - 0.0205000000000000i;0.612400000000000 + 0.00000000000000i,0.612400000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.755900000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,0.535300000000000 + 0.00000000000000i,0.535300000000000 + 0.00000000000000i,-0.0205000000000000 - 0.155700000000000i,-0.0205000000000000 + 0.155700000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,-0.0851000000000000 + 0.0425000000000000i,-0.0851000000000000 - 0.0425000000000000i,0.00000000000000 + 0.558900000000000i,-0.00000000000000 - 0.558900000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,-0.450100000000000 - 0.0425000000000000i,-0.450100000000000 + 0.0425000000000000i,0.0205000000000000 - 0.403200000000000i,0.0205000000000000 + 0.403200000000000i]];
lissa2Plot
angle(vm)
angle(vn)
lissa2Plot
plot(xm,xn,'-','color',linecolor,'linewidth',3);xlimt([-1 13])
plot(xm,xn,'-','color',linecolor,'linewidth',3);xlimit([-1 13])
plot(xm,xn,'-','color',linecolor,'linewidth',3);axis( [-1 1 -1 1] )
lissa2Plot
close all
lissa2Plot
set(gca, 'size',[30 30])
lissa2Plot
clear
ew=[-1.84688150281763e-16 - 0.325472277452060i
3.21300396058951e-17 + 0.130188910980824i
-1.33004764142749e-16 - 0.325472277452060i
-5.72627262363905e-17 + 0.520755643923295i
-0.325472277452060 + 3.01289563405940e-17i
0.130188910980824 - 5.14831112757637e-17i
-0.325472277452060 + 5.78845319562229e-17i
0.520755643923295 + 0.00000000000000i]
lissa2Plot
ew=[-1.84688150281763e-16 - 0.325472277452060i
3.21300396058951e-17 + 0.130188910980824i
-1.33004764142749e-16 - 0.325472277452060i
-5.72627262363905e-17 + 0.520755643923295i
-0.325472277452060 + 3.01289563405940e-17i
0.130188910980824 - 5.14831112757637e-17i
-0.325472277452060 + 5.78845319562229e-17i
0.520755643923295 + 0.00000000000000i]
vm=ew(m);vn=ew(n);  % 书洁修正
[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);
if area_i_j > 0; linecolor = 'r'; end
if abs(area_i_j) < 0.00001; linecolor = 'k'; end
if area_i_j < 0; linecolor = 'b'; end
subplot(2,3,figk)
t=0:0.1:10
xm = abs(vm) * cos((t + angle(vm)))
xn = abs(vn) * cos(t + angle(vn))
plot(xm,xn,'-','color',linecolor,'linewidth',3);axis( [-1 1 -1 1]);axis square;
%set(gca,'xtick',[]); set(gca,'ytick',[]);
xlabel(strcat('\it{x}_',num2str(m)));ylabel(strcat('\it{x}_',num2str(n)));
grid on
figk=1
for m=1:3
for n=m+1:4
vm=ew(m);vn=ew(n);  % 书洁修正
[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);
if area_i_j > 0; linecolor = 'r'; end
if abs(area_i_j) < 0.00001; linecolor = 'k'; end
if area_i_j < 0; linecolor = 'b'; end
subplot(4,7,figk)
t=0:0.1:10
xm = abs(vm) * cos((t + angle(vm)))
xn = abs(vn) * cos(t + angle(vn))
plot(xm,xn,'-','color',linecolor,'linewidth',3);axis( [-1 1 -1 1]);axis square;
%set(gca,'xtick',[]); set(gca,'ytick',[]);
xlabel(strcat('\it{x}_',num2str(m)));ylabel(strcat('\it{x}_',num2str(n)));
grid on
figk=figk+1
end
end
lissa2Plot
ew=[-1.84688150281763e-16 - 0.325472277452060i
3.21300396058951e-17 + 0.130188910980824i
-1.33004764142749e-16 - 0.325472277452060i
-5.72627262363905e-17 + 0.520755643923295i
-0.325472277452060 + 3.01289563405940e-17i
0.130188910980824 - 5.14831112757637e-17i
-0.325472277452060 + 5.78845319562229e-17i
0.520755643923295 + 0.00000000000000i]
lissa2Plot
ew=[-1.84688150281763e-16 - 0.325472277452060i
3.21300396058951e-17 + 0.130188910980824i
-1.33004764142749e-16 - 0.325472277452060i
-5.72627262363905e-17 + 0.520755643923295i
-0.325472277452060 + 3.01289563405940e-17i
0.130188910980824 - 5.14831112757637e-17i
-0.325472277452060 + 5.78845319562229e-17i
0.520755643923295 + 0.00000000000000i]
%-- 2020/11/19 23:39 --%
[1:4].*1.4
[1:4].*[1:4]
[1:4]*[1:4]
s16 = []; for I=1:4; for J=1:4; s16 = [s16; I J I*10+J max(I,J)*100+min(I,J)];end;end
unique(s16(:,4))
A=load('C:\Users\Think\Desktop\A25sj.csv');
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),JA(:,2))*10+max(A(:,1),A(:,2));
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),A(:,2))*10+max(A(:,1),A(:,2));
s16 = []; for I=1:4; for J=1:4; s16 = [s16; I J I*10+J min(I,J)*10+max(I,J)];end;end;s10=unique(s16(:,4));
s16 = []; for I=1:4; for J=1:4; s16 = [s16; I J I*10+J min(I,J)*10+max(I,J)];end;end;S10=unique(s16(:,4));S16=(s16(:,3));
M10=zeros(10); for K=1:length(A(:,1))-1; m=find(S10==A(K,4));m=find(S10==A(K+1,4));M10(m,n)=M10(m,n)+1;end
M10=zeros(10); for K=1:length(A(:,1))-1; m=find(S10==A(K,4));n=find(S10==A(K+1,4));M10(m,n)=M10(m,n)+1;end
sum(sum(M10))
Asymm = M10-M10';
A=load('C:\Users\Think\Desktop\A25sj.csv');
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),A(:,2))*10+max(A(:,1),A(:,2));
s16 = []; for I=1:4; for J=1:4; s16 = [s16; I J I*10+J min(I,J)*10+max(I,J)];end;end;S10=unique(s16(:,4));S16=(s16(:,3));
M10=zeros(10); for K=1:length(A(:,1))-1; m=find(S10==A(K,4));n=find(S10==A(K+1,4));M10(m,n)=M10(m,n)+1;end
sum(sum(M10))
Asymm = M10-M10';
clear
A=load('C:\Users\Think\Desktop\A44sj.csv');
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),A(:,2))*10+max(A(:,1),A(:,2));
s16 = []; for I=1:4; for J=1:4; s16 = [s16; I J I*10+J min(I,J)*10+max(I,J)];end;end;S10=unique(s16(:,4));S16=(s16(:,3));
M10=zeros(10); for K=1:length(A(:,1))-1; m=find(S10==A(K,4));n=find(S10==A(K+1,4));M10(m,n)=M10(m,n)+1;end
sum(sum(M10))
Asymm = M10-M10';
%-- 2020/11/20 4:33 --%
A=load('C:\Users\Think\Desktop\A25sj.csv');
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),A(:,2))*10+max(A(:,1),A(:,2));
S16=unique(A(:,3));
M16=zeros(16); for K=1:length(A(:,1))-1; m=find(S16==A(K,3));n=find(S16==A(K+1,3));M16(m,n)=M16(m,n)+1;end
M16a=M16-M16'
A=load('C:\Users\Think\Desktop\A44sj.csv');
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),A(:,2))*10+max(A(:,1),A(:,2));
S16=unique(A(:,3));
M16=zeros(16); for K=1:length(A(:,1))-1; m=find(S16==A(K,3));n=find(S16==A(K+1,3));M16(m,n)=M16(m,n)+1;end
M16a=M16-M16'
%-- 2020/11/20 12:49 --%
cc=[1	2	0	0	0
1	3	0	0	0
1	6	0.4284	0	0
1	7	0.4284	0	0
2	3	0	-0.172	-0.3518
2	5	0.4284	0	0
2	6	-0.1428	0.6841	0.5758
2	7	-0.1428	0.7596	-0.6702
3	5	0.4284	0	0
3	6	-0.1428	0.7596	-0.6702
3	7	-0.1428	0.8867	0.995
5	6	0	0	0
5	7	0	0	0
6	7	0	-0.172	-0.3518
1	5	-1.1781	0	0
]
b=[-4.8088E-17	0	0	0.003047619	0.001763533	-0.0078	-0.001333333	0.001687764	0.008016878	0.000897077
-4.8088E-17	0	0	-0.011047619	-0.00168661	0.0034	-0.000333333	-0.009113924	-0.001687764	-0.003411542
0.3927	0	0	0.009371429	0.010396011	0.017933333	0.018866667	0.007594937	0.020506329	0.014111451
0.3927	0	0	0.021104762	0.00182906	0.0204	0.006333333	0.015189873	0.014514768	0.013228633
0	-0.071439963	0.036010289	0.011047619	0.002541311	-0.009066667	0.000666667	0.008101266	0.006751055	0.003340208
0.3927	0	0	0.00952381	0.003461538	0.0206	0.004733333	0.014852321	0.017721519	0.01181542
-0.1309	0.900211091	0.077437117	-0.000838095	-0.000384615	0.003533333	-0.011066667	0.001687764	-0.008523207	-0.002598581
-0.1309	-0.143086216	-0.27330414	-0.001447619	-0.001185185	-0.010066667	0.004333333	-0.011814346	-0.014261603	-0.005740348
0.3927	0	0	0.010666667	0.007752137	0.015066667	0.0114	0.014261603	0.012658228	0.01196755
-0.1309	-0.143086216	-0.27330414	-0.008838095	-0.006017094	-0.002733333	-0.002066667	-0.005991561	-0.004725738	-0.005062081
-0.1309	0.028412595	0.981336815	-0.003352381	0.000575499	-0.004333333	-0.003666667	0.002447257	0.002194093	-0.001022589
-4.8088E-17	0	0	-0.006247619	0.001039886	-0.0012	0.001466667	-0.001350211	0.002362869	-0.000654735
-4.8088E-17	0	0	0.008914286	-0.000490028	0.0036	0.0046	-0.004725738	-0.009113924	0.000464099
0	-0.071439963	0.036010289	-0.00327619	0.000780627	-0.003266667	0.0004	0.00464135	0.006919831	0.001033159
-1.1782	0	0	-0.037714286	-0.014749288	-0.053066667	-0.032866667	-0.03907173	-0.049367089	-0.037805954
]
x1=cc(:,3); x2=cc(:, 4);x3=cc(:,5);y=b(:,10); %wi8i4
X = [ones(size(y)) x1 x2 x3];
[bxx,bint,r,rint,stats] = regress(y,X);bb= [bxx bint stats']
clear
yw=[  -0.0000 - 0.0000i	   0.0000 + 0.0000i
-0.0000 + 0.5353i	   0.2850 + 0.3201i
-0.0425 - 0.0851i	  -0.2850 + 0.1432i
0.0425 - 0.4501i	   0.0000 - 0.4633i
0.0000 - 0.0000i	   0.0000 + 0.0000i
0.5353 + 0.0000i	   0.3201 - 0.2850i
-0.0851 + 0.0425i	   0.1432 + 0.2850i
-0.4501 - 0.0425i	  -0.4633 + 0.0000i
]
abs(yw)
y=yw(:,1);w=yw(:,2)
abs(y)-abs(w)
abs(y-w))
abs(y-w)
abs(y-w*ezp(i*pi/4))
abs(y-w*exp(i*pi/4))
sum(abs(y-w))
sum(abs(y-w*exp(i*pi/4)))
t=1:16; sum(abs(y-w*exp(i*2*pi*t/16)))
t=1:16; sum(abs(y-w*exp(i*2*pi*t./16)))
err= []; for t=1:16; err= [err; t sum(abs(y-w*exp(i*2*pi*t/16)))];end
err= []; for t=1:32; err= [err; t sum(abs(y-w*exp(i*2*pi*t/16)))];end
plot(err(:,end))
err= []; for t=1:32; err= [err; t sum(abs(y(:,2:4)-w(:,2:4)*exp(i*2*pi*t/16)))];end
err= []; for t=1:32; err= [err; t sum(abs(y(2:4)-w(2:4)*exp(i*2*pi*t/16)))];end
plot(err(:,end))
err= []; for t=1:320; err= [err; t sum(abs(y(2:4)-w(2:4)*exp(i*2*pi*t/16)))];end
y(2:4)
polar(y(2:4))
plot(y(2:4))
hold on;plot(w(2:4))
grid on
w(2:4)
%-- 2020/11/21 0:24 --%
findBest_4i
%-- 2020/11/21 4:40 --%
findBest_4i
findBest_4i
findBest_4i
plot(rR(:,end))
findBest_4i
plot(rR(:,end))
findBest_4i
plot(rR(:,end))
findBest_4i
plot(rR(:,end))
findBest_4i
plot(rR(:,end))
findBest_4i
plot(rR(:,end))
%-- 2020/11/21 11:48 --%
DesignONeill
%-- 2020/11/21 13:26 --%
DesignONeill
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8;
vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];
end;end;af=[af ae];end;
evt=round(af,4); evl=round(diag(eigen_value),4);
unique(evt(:,1))
YQM20200830
DesignONeill
unique(evt(:,1))
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8;
vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];
end;end;af=[af ae];end;
evt=round(af,4); evl=round(diag(eigen_value),4);
DesignONeill
evt
evl'
DesignONeill
A43 = 2;
DesignONeill(A43)
A43 = 1/2;
DesignONeill(A43)
unique(evt(:,1))
DesignONeill
DesignONeill(A43)
unique(evt(:,1))
unique(evt(:,2))
unique(evt(:,3))
unique(evt(:,4))
unique(evt(:,5))
unique(evt(:,6))
unique(evt(:,7))
%-- 2020/11/22 5:20 --%
findBest_4i
lsline(x1,y)
lsline(ex1,y)
data = [ 12	-4.8088E-17	0	0	8	3.438888885	-23.4	-4	4	19	900	1	2	0	0	0
13	-4.8088E-17	0	0	-29	-3.288888888	10.2	-1	-21.6	-4	900	1	3	0	0	0
15	-1.1782	0	0	-99	-28.76111112	-159.2	-98.6	-92.6	-117	910	1	5	-1.2852	0	0
16	0.3927	0	0	24.6	20.27222221	53.8	56.6	18	48.6	910	1	6	0.4284	0	0
17	0.3927	0	0	55.4	3.566666669	61.2	19	36	34.4	910	1	7	0.4284	0	0
23	0	-0.071439963	0.036010289	29	4.955555547	-27.2	2	19.2	16	901	2	3	0	-0.172	-0.3518
25	0.3927	0	0	25	6.749999992	61.8	14.2	35.2	42	910	2	5	0.4284	0	0
26	-0.1309	0.900211091	0.077437117	-2.2	-0.749999997	10.6	-33.2	4	-20.2	911	2	6	-0.1428	0.6841	0.5758
27	-0.1309	-0.143086216	-0.27330414	-3.8	-2.31111111	-30.2	13	-28	-33.8	911	2	7	-0.1428	0.7596	-0.6702
35	0.3927	0	0	28	15.11666666	45.2	34.2	33.8	30	910	3	5	0.4284	0	0
36	-0.1309	-0.143086216	-0.27330414	-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2	911	3	6	-0.1428	0.7596	-0.6702
37	-0.1309	0.028412595	0.981336815	-8.8	1.12222222	-13	-11	5.8	5.2	911	3	7	-0.1428	0.8867	0.995
56	-4.8088E-17	0	0	-16.4	2.027777784	-3.6	4.4	-3.2	5.6	900	5	6	0	0	0
57	-4.8088E-17	0	0	23.4	-0.955555556	10.8	13.8	-11.2	-21.6	900	5	7	0	0	0
67	0	-0.071439963	0.036010289	-8.6	1.522222224	-9.8	1.2	11	16.4	901	6	7	0	-0.172	-0.3518
];
rr = [];
scatter(data(:,2),data(:,5),'ro')
fig_name_saved =   'testsave2Dfig.png'
figureFilePath = 'E:/TMP/';xlab ='Theo', ylab='ONeill1987'
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
scatter(data(:,2),data(:,5),'ro')
fig_name_saved =   'testsave2Dfig.png'
figureFilePath = 'E:/TMP/';xlab ='Theo', ylab='ONeill1987'
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
scatter(data(:,2),data(:,5),'ro');lsline(data(:,2),data(:,5))
fig_name_saved =   'testsave2Dfig'
figureFilePath = 'E:/TMP/';xlab ='Theo', ylab='ONeill1987'
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
ax1 = scatter(data(:,2),data(:,5),'ro');lsline(ax1)
fig_name_saved =   'testsave2Dfig'
figureFilePath = 'E:/TMP/';xlab ='Theo', ylab='ONeill1987'
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
ax1 = scatter(data(:,2),data(:,5),'ro');lsline
ax1 = scatter(data(:,2),data(:,5),'ro');kk = lsline
set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),'ro','fontSize',5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),'ro','markSize',5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),'ro','markerSize',5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),'ro','LineWidth',2.5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),'ro','MarkerFaceColor',[0 .7 .7],'LineWidth',2.5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),40,'MarkerFaceColor',[0 .7 .7],'LineWidth',2.5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),40,'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),40,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),40,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),40,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);kk = lsline; set(kk,'linewidth',4)
fig_name_saved =   'testsave2Dfig'
figureFilePath = 'E:/TMP/';xlab ='Theo', ylab='ONeill1987'
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
ax1 = scatter(data(:,2),data(:,5),'x','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);kk = lsline; set(kk,'linewidth',4)
fig_name_saved =   'testsave2Dfig'
figureFilePath = 'E:/TMP/';xlab ='Theo', ylab='ONeill1987'
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
ax1 = scatter(data(:,2),data(:,5),'X','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
%-- 2020/11/22 14:49 --%
data = [ 12	-4.8088E-17	0	0	8	3.438888885	-23.4	-4	4	19	900	1	2	0	0	0
13	-4.8088E-17	0	0	-29	-3.288888888	10.2	-1	-21.6	-4	900	1	3	0	0	0
15	-1.1782	0	0	-99	-28.76111112	-159.2	-98.6	-92.6	-117	910	1	5	-1.2852	0	0
16	0.3927	0	0	24.6	20.27222221	53.8	56.6	18	48.6	910	1	6	0.4284	0	0
17	0.3927	0	0	55.4	3.566666669	61.2	19	36	34.4	910	1	7	0.4284	0	0
23	0	-0.071439963	0.036010289	29	4.955555547	-27.2	2	19.2	16	901	2	3	0	-0.172	-0.3518
25	0.3927	0	0	25	6.749999992	61.8	14.2	35.2	42	910	2	5	0.4284	0	0
26	-0.1309	0.900211091	0.077437117	-2.2	-0.749999997	10.6	-33.2	4	-20.2	911	2	6	-0.1428	0.6841	0.5758
27	-0.1309	-0.143086216	-0.27330414	-3.8	-2.31111111	-30.2	13	-28	-33.8	911	2	7	-0.1428	0.7596	-0.6702
35	0.3927	0	0	28	15.11666666	45.2	34.2	33.8	30	910	3	5	0.4284	0	0
36	-0.1309	-0.143086216	-0.27330414	-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2	911	3	6	-0.1428	0.7596	-0.6702
37	-0.1309	0.028412595	0.981336815	-8.8	1.12222222	-13	-11	5.8	5.2	911	3	7	-0.1428	0.8867	0.995
56	-4.8088E-17	0	0	-16.4	2.027777784	-3.6	4.4	-3.2	5.6	900	5	6	0	0	0
57	-4.8088E-17	0	0	23.4	-0.955555556	10.8	13.8	-11.2	-21.6	900	5	7	0	0	0
67	0	-0.071439963	0.036010289	-8.6	1.522222224	-9.8	1.2	11	16.4	901	6	7	0	-0.172	-0.3518
];
fig_name_saved =   'testsave2DfigYY'
figureFilePath = 'E:/TMP/';xlab ='TheoryYY', ylab='ONeill1987'
ax1 = scatter(data(:,2),data(:,5),'X','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
%-- 2020/11/23 14:33 --%
data = [ 12	-4.8088E-17	0	0	8	3.438888885	-23.4	-4	4	19	900	1	2	0	0	0
13	-4.8088E-17	0	0	-29	-3.288888888	10.2	-1	-21.6	-4	900	1	3	0	0	0
15	-1.1782	0	0	-99	-28.76111112	-159.2	-98.6	-92.6	-117	910	1	5	-1.2852	0	0
16	0.3927	0	0	24.6	20.27222221	53.8	56.6	18	48.6	910	1	6	0.4284	0	0
17	0.3927	0	0	55.4	3.566666669	61.2	19	36	34.4	910	1	7	0.4284	0	0
23	0	-0.071439963	0.036010289	29	4.955555547	-27.2	2	19.2	16	901	2	3	0	-0.172	-0.3518
25	0.3927	0	0	25	6.749999992	61.8	14.2	35.2	42	910	2	5	0.4284	0	0
26	-0.1309	0.900211091	0.077437117	-2.2	-0.749999997	10.6	-33.2	4	-20.2	911	2	6	-0.1428	0.6841	0.5758
27	-0.1309	-0.143086216	-0.27330414	-3.8	-2.31111111	-30.2	13	-28	-33.8	911	2	7	-0.1428	0.7596	-0.6702
35	0.3927	0	0	28	15.11666666	45.2	34.2	33.8	30	910	3	5	0.4284	0	0
36	-0.1309	-0.143086216	-0.27330414	-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2	911	3	6	-0.1428	0.7596	-0.6702
37	-0.1309	0.028412595	0.981336815	-8.8	1.12222222	-13	-11	5.8	5.2	911	3	7	-0.1428	0.8867	0.995
56	-4.8088E-17	0	0	-16.4	2.027777784	-3.6	4.4	-3.2	5.6	900	5	6	0	0	0
57	-4.8088E-17	0	0	23.4	-0.955555556	10.8	13.8	-11.2	-21.6	900	5	7	0	0	0
67	0	-0.071439963	0.036010289	-8.6	1.522222224	-9.8	1.2	11	16.4	901	6	7	0	-0.172	-0.3518
];
fig_name_saved =   'testsave2DfigYY'
figureFilePath = 'E:/TMP/';xlab ='TheoryYY', ylab='ONeill1987'
class_color=[ .5 .5 .5 ; 1 0 .5 ;  0 0 1 ; 0 .5 1 ;]
ax1 = scatter(data(:,2),data(:,5),'X','MarkerEdgeColor',class_color(2),'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);
kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),'X','MarkerEdgeColor',class_color(2,:),'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);
kk = lsline; set(kk,'linewidth',4)
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
ax1 = scatter(data(:,2),data(:,5),'X','MarkerEdgeColor',class_color(3,:),'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
ax1 = scatter(data(:,2),data(:,5),'X','MarkerEdgeColor',class_color(4,:),'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);
kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_color(4,:),'MarkerFaceColor',[0 .7 .7],'LineWidth',4.5);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_color(1,:),'MarkerFaceColor',class_color(1,:),'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
class_color=[ .5 .5 .5 ; 1 0 0 ;  0 0 1 ; 1 .5 1 ;];
class_c = class_color(1,:)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
class_color=[ 1 0 0 ; .5 .5 .5 ;   0 0 1 ; 1 .5 1 ;];
class_c = class_color(1,:)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
class_color=[ 1 0 0 ; .5 .5 .5 ; .3 0 1 ; 1 .5 1 ;];
class_c = class_color(1,:)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
class_color=[ 1 0 0 ; .5 .5 .5 ; .3 0 1 ; 1 .5 1 ;];
class_c = class_color(2,:)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
class_color=[ 1 0 0 ; .13 .52 .13 ; .3 0 1 ; 1 .5 1 ;];
class_c = class_color(2,:)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
ax1 = scatter(data(:,2),data(:,5),10,'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'Markersize',6,'LineWidth',6);
ax1 = scatter(data(:,2),data(:,5), 'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'Markersize',6,'LineWidth',6);
ax1 = scatter(data(:,2),data(:,5), 'Markersize',6, 'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c );
ax1 = scatter(data(:,2),data(:,5),56, 'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c );
ax1 = scatter(data(:,2),data(:,5),126, 'MarkerFaceColor',class_c );
ax1 = scatter(data(:,2),data(:,5), 'Markersize',126, 'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
data = [  % Sam28	dimX	dimY	t8i(4)	t4i1	t4i2	Oneill(7)	Binmore	AIBT	ATBI(10)	AIBI	ATBT	total	E8i	E4i	Classid(15)
12	1	2	-4.81E-17	0	0	8	3.438888885	-23.4	-4	4	19	7.038888885	0	0	0
13	1	3	-4.81E-17	0	0	-29	-3.288888888	10.2	-1	-21.6	-4	-48.68888889	0	0	0
14	1	4	-4.81E-17	0	0	21	-0.149999998	13.2	5	17.6	-15	41.65	0	0	0
15	1	5	-1.1782	0	0	-99	-28.76111112	-159.2	-98.6	-92.6	-117	-595.1611111	1	0	2
16	1	6	0.3927	0	0	24.6	20.27222221	53.8	56.6	18	48.6	221.8722222	1	0	2
17	1	7	0.3927	0	0	55.4	3.566666669	61.2	19	36	34.4	209.5666667	1	0	2
18	1	8	0.3927	0	0	19	4.922222222	44.2	23	38.6	34	163.7222222	1	0	2
23	2	3	0	-0.071439963	0.036010289	29	4.955555547	-27.2	2	19.2	16	43.95555555	0	1	1
24	2	4	0	0.071549996	-0.036019622	-21	-1.516666671	3.8	-6	-15.2	3	-36.91666667	0	1	1
25	2	5	0.3927	0	0	25	6.749999992	61.8	14.2	35.2	42	184.95	1	0	2
26	2	6	-0.1309	0.900211091	0.077437117	-2.2	-0.749999997	10.6	-33.2	4	-20.2	-41.75	1	1	3
27	2	7	-0.1309	-0.143086216	-0.27330414	-3.8	-2.31111111	-30.2	13	-28	-33.8	-85.11111111	1	1	3
28	2	8	-0.1309	0.741983761	0.195831955	-19	-3.688888886	-42.2	6	-11.2	12	-58.08888889	1	1	3
34	3	4	0	-0.071441185	0.036059535	-1.33E-15	1.666666663	-17	1	-2.4	12	-4.733333337	0	1	1
35	3	5	0.3927	0	0	28	15.11666666	45.2	34.2	33.8	30	186.3166667	1	0	2
36	3	6	-0.1309	-0.143086216	-0.27330414	-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2	-74.73333332	1	1	3
37	3	7	-0.1309	0.028412595	0.981336815	-8.8	1.12222222	-13	-11	5.8	5.2	-20.67777778	1	1	3
38	3	8	-0.1309	-0.131098343	-0.707913246	4	-4.505555552	-24	-17	-25.4	-24	-90.90555555	1	1	3
45	4	5	0.3927	0	0	46	6.894444446	52.2	50.2	23.6	45	223.8944444	1	0	2
46	4	6	-0.1309	-0.756919893	0.195831955	0.8	-7.788888883	-56.2	-17.2	-7.8	-17.2	-105.3888889	1	1	3
47	4	7	-0.1309	0.114632307	-0.707913246	-42.8	-2.377777774	-18	-21	-13.8	-5.8	-103.7777778	1	1	3
48	4	8	-0.1309	-0.610696192	0.511996947	-4	3.27222222	22	-12	-2	-22	-14.72777778	1	1	3
56	5	6	-4.81E-17	0	0	-16.4	2.027777784	-3.6	4.4	-3.2	5.6	-11.17222222	0	0	0
57	5	7	-4.81E-17	0	0	23.4	-0.955555556	10.8	13.8	-11.2	-21.6	14.24444444	0	0	0
58	5	8	-4.81E-17	0	0	-7	-1.072222226	-7.2	-18.2	14.4	16	-3.072222226	0	0	0
67	6	7	0	-0.071439963	0.036010289	-8.6	1.522222224	-9.8	1.2	11	16.4	11.72222222	0	1	1
68	6	8	0	0.16585302	-0.036019622	-7.8	0.505555554	6.2	3.2	-14.2	-10.8	-22.89444445	0	1	1
78	7	8	0	0.032521274	0.036059535	14.8	0.56666667	1	15	-0.2	-5.2	25.96666667	0	1	1
];
field_id=[ 'Sam28'	'dimX'	'dimY'	't8i'	't4i1'	't4i2'	'ONeill'	'Binmore'	'AIBT'	'ATBI'	'AIBI'	'ATBT'	'total'	'E8i'	'E4i'	'Class id'];
for Yfig_ad=7:13
fig_name_saved =   strcat('testsave2DfigYY',num2str(Yfig_ad)
figureFilePath = 'E:/TMP/';xlab ='Theory(8i)', ylab = field_id(Yfig_ad)
class_color=[ 1 0 0 ; .13 .52 .13 ; .3 0 1 ; 1 .5 1 ;];
class_c = class_color(2,:)
ax1 = scatter(data(:,2),data(:,5), 126, 'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
end
for Yfig_id=7:13
fig_name_saved =   strcat('Yfig',num2str(Yfig_id)
figureFilePath = 'E:/TMP/';xlab ='Theory(8i)', ylab = field_id(Yfig_id)
class_color=[ 1 0 0 ; .13 .52 .13 ; .3 0 1 ; 1 .5 1 ;];
class_c = class_color(2,:)
ax1 = scatter(data(:,2),data(:,5), 126, 'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
end
for Yfig_id=7:13
fig_name_saved =   strcat('Yfig',num2str(Yfig_id))
figureFilePath = 'E:/TMP/';xlab ='Theory(8i)', ylab = field_id(Yfig_id)
class_color=[ 1 0 0 ; .13 .52 .13 ; .3 0 1 ; 1 .5 1 ;];
class_c = class_color(2,:)
ax1 = scatter(data(:,2),data(:,5), 126, 'MarkerEdgeColor',class_c,'MarkerFaceColor',class_c,'LineWidth',6);
kk = lsline; set(kk,'linewidth',4)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
end
testsave2Dfig
data(:,15)
class_c = class_color(data(:,16))
data(:,16)
class_c = class_color(data(:,16)+1)
class_c = class_color(data(:,16)+1,:)
testsave2Dfig
ax1 = scatter(data(:,4),data(:,Yfig_id), 126, 'MarkerEdgeColor',class_c,'LineWidth',6);
ax1 = scatter(data(:,4),data(:,Yfig_id), 126, 'MarkerEdgeColor',[class_c(:,1) class_c(:,2) class_c(:,3)], 'LineWidth',6);
testsave2Dfig
for Yfig_id=7:13
fig_name_saved = strcat('Yfig',num2str(Yfig_id))
figureFilePath = 'E:/TMP/';xlab ='Theory(8i)', ylab = field_id(Yfig_id)
class_color=[ 1 0 0 ; .13 .52 .13 ; .3 0 1 ; 1 .5 1 ;];
class_c = class_color(data(:,16)+1,:)
for n=1:28
ax1 = scatter(data(n,4),data(n,Yfig_id), 126, 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',6);hold on
end
ax2 = scatter(data(:,4),data(:,Yfig_id), 1, 'MarkerEdgeColor','w')
kk = lsline; set(kk,'linewidth',2)
plotfilesaveTo=save2Dfig(xlab, ylab, fig_name_saved, figureFilePath)
end
kk = lsline(ax2)
testsave2Dfig
save2Dfig
testsave2Dfig
line([1,2],[3,4])
testsave2Dfig
line([0,0], [-105,60],'LineWidth',2);hold on
testsave2Dfig
line([-1.5,1.5],[0,0], 'LineWidth',1,  'color',[.5 .5 .5] );hold on
line([0,0], [-105,60],'LineWidth',1,  'color',[.5 .5 .5]);hold on
testsave2Dfig
xlim([-1.5 0.5])
testsave2Dfig
kk = lsline ; set(kk,'linewidth',5,'linetype','dash');hold on;
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
testsave2Dfig
%-- 2020/11/23 20:50 --%
testsave2Dfig
save2Dfig
matlab.graphics.internal.PrintVertexChecker
testsave2Dfig
data(:,7)./TotalPeriod(1)
data(:,7:13)./TotalPeriod(1:7)
for kt=7:13;data(:,kt)=data(:,kt)./TotalPeriod(kt);end
for kt=7:13;data(:,kt)=data(:,kt)./TotalPeriod(kt-6);end
testsave2Dfig
for kt=7:13;data(:,kt)=data(:,kt)./TotalPeriod(kt-6);end
MinMax=[min(data) max(data) ]
MinMax=[min(data);max(data) ]
MinMax=[min(data);max(data) ]';
testsave2Dfig
save2Dfig
testsave2Dfig
uiopen('D:\MATLAB\R2016a\bin\html\Yfig8.fig',1)
xlabel(xlab,'FontSize',25,'interpreter','latex');
ylabel(ylab,'FontSize',25,'interpreter','latex');  axis square;    grid on;
set(gca,'FontSize',20,'linewidth',2);hold on;
box on;
grid off
set(gcf,'papersize',[15 15],'paperposition',[0,0,15,15]);hold on;
%           set(gca,'xlim',[100 400],'XTick',100:50:400)
%           set(gca,'ylim',[8 24],'YTick',8:4:24)
plotfilesaveTo =strcat( figureFilePath, fig_name_saved,'2.fig')
saveas(gcf, plotfilesaveTo);
plotfilesaveTo =strcat( figureFilePath, fig_name_saved,'2.pdf')
saveas(gcf, plotfilesaveTo);
plotfilesaveEPS =strcat( figureFilePath, fig_name_saved ,'2.eps')
saveas(gcf, plotfilesaveEPS,'psc2');
set(gca,'FontSize',20,'linewidth',2);hold on;
box on;
grid off
set(gcf,'papersize',[15 15],'paperposition',[0,0,15,15]);hold on;
%           set(gca,'xlim',[100 400],'XTick',100:50:400)
%           set(gca,'ylim',[8 24],'YTick',8:4:24)
plotfilesaveTo =strcat( figureFilePath, fig_name_saved,'2.fig')
saveas(gcf, plotfilesaveTo);
plotfilesaveTo =strcat( figureFilePath, fig_name_saved,'2.pdf')
saveas(gcf, plotfilesaveTo);
plotfilesaveEPS =strcat( figureFilePath, fig_name_saved ,'2.eps')
saveas(gcf, plotfilesaveEPS,'psc2');
close
testsave2Dfig
%-- 2020/11/24 13:30 --%
DesignONeill
A43 = 1/2;
DesignONeill(A43)
%-- 2020/11/24 22:30 --%
testsave2Dfig
a=data(find(data(:,16)==1),:))
a=data(find(data(:,16)==1),:)
boxplot(a(:,7:13))
J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J),:);boxplot(a(:,7:13)); end
for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J),:);boxplot(a(:,7:13)); end
for J=0:4; subplot(5,1,J+1);a=data(find(data(:,16)==J),:);boxplot(a(:,7:13)); end
for J=0:4; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,7:13)); end
for J=0:4; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); end
for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); end
for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); grid on;end
for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); ylim([-.2 .2]);grid on;end
for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); ylim([-.1 .1]);grid on;end
for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); ylim([-.1 .1]*.2);grid on;end
for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); ylim([-.1 .1]*.2);axis square;grid on;end
%-- 2020/11/25 16:21 --%
testsave2Dfig
s9=[0.077437117 -41.75
-0.27330414 -85.11111111
0.195831955 -58.08888889
-0.27330414 -74.73333332
0.981336815 -20.67777778
-0.707913246 -90.90555555
0.195831955 -105.3888889
-0.707913246 -103.7777778
0.511996947 -14.72777778];
clf
for n=1:28
ax1 = scatter(data(n,4),data(n,Yfig_id), 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);hold on
end
clf
class_c(n,1)
for n=1:9
ax1 = scatter(s9(n,1),s9(n,2), 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);hold on
end
for n=1:9
ax1 = scatter(s9(n,1),s9(n,2), 129, class_l(4), 'MarkerEdgeColor',[class_c(4,1) class_c(4,2) class_c(4,3)], 'LineWidth',2);hold on
end
class_l(4)
testsave2Dfig
s9 = data(find(data(:,16) == 3), :)
testsave2Dfig
s9 = data(find(data(:,16) == 3), :)
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2)
expO9=sum(da')'/dt
the4i2 = s9(:,6)
ax2 = scatter(data(:,4),data(:,Yfig_id), 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
ax2 = scatter(the4i2, expO9129, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
scatter(the4i2, expO9,129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
testsave2Dfig
clf
s9 = data_ori(find(data_ori(:,16) == n), :)
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2)
expO9 = sum(da')'/dt
the4i2 = s9(:,6)
ax2 = scatter(the4i2, expO9129, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
ax2 = scatter(the4i2, expO,9129, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
ax2 = scatter(the4i2, expO9, 129, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
s9 = data_ori(find(data_ori(:,16) == n), :)
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2)
expO9 = sum(da')'/dt
the4i2 = s9(:,6)
ax2 = scatter(the4i2, expO9, 129, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
s9 = data_ori(find(data_ori(:,16) == 3), :)
n=3;
s9 = data_ori(find(data_ori(:,16) == n), :)
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2)
expO9 = sum(da')'/dt
the4i2 = s9(:,6)
ax2 = scatter(the4i2, expO9, 129, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
ax2 = scatter(the4i2, expO9, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);
scatter(the4i2, expO9, 129, class_l(n+1), 'MarkerEdgeColor',[class_c(n+1,1) class_c(n+1,2) class_c(n+1,3)], 'LineWidth',2);
scatter(the4i2, expO9, 129, 's', 'MarkerEdgeColor',[1 0 0], 'LineWidth',2);
ax2 = scatter(the4i2, expO9, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, 's', 'MarkerEdgeColor',[1 0 0], 'LineWidth',2);
xlab ='Theory(4i(2))', ylab = '$L_{mn}^E$ (Pooled)',
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
testsave2Dfig
clf; n=3;
s9 = data_ori(find(data_ori(:,16) == n), :)
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2)
expO9 = sum(da')'/dt
the4i2 = s9(:,6)
ax2 = scatter(the4i2, expO9, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, 's', 'MarkerEdgeColor',[1 0 0], 'LineWidth',2);
xlab ='Theory ($\sigma 4i_{2}$)', ylab = '$L_{mn}^{P11}$',
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
xlab ='Theory ($\sigma_{4i_{2}}$)', ylab = '$L_{mn}^{P11}$',
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
xlab ='Theory ($\sigma_{4i_{2}}$)', ylab = 'Experiment $L_{mn}^{P11}$',
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
xlab ='Theory ($\sigma_{4i_{2}}$)', ylab = 'Experiment ($L_{mn}^{P11}$)',
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
testsave2Dfig
clf; n=3;
s9 = data_ori;
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2);
expO9 = sum(da')'/dt
the4i2 = s9(:,6)
testsave2Dfig
fig_name_saved = 'fig4i2';
text(-0.1, -0.001*1.5,'(b)','fontsize',20);
text(-0.1, -0.001*1.1,'(b)','fontsize',20);
testsave2Dfig
%-- 2020/11/26 19:56 --%
testsave2Dfig
s9 = data(find(data(:,16)==3),:)
ttest(s9(:,7))
rr=[];for j=7:12;rr=[rr;s9(:,j)];end
[h p]=ttest(rr)
mean(s9)
mean(s9(:,7:13))
%-- 2020/11/27 1:20 --%
testsave2Dfig
text(-1+0.02,-0.5,'$x10^3$')
text(-1+0.02,-0.5,'x10^3')
text(-1+0.02,-0.5,'x10^3','fontsize',20)
text(-1+0.03,-0.4,'x10^3','fontsize',20)
testsave2Dfig
a=3+5i
plot(a)
plot(a,'r')
plot(a,'ro');
quiver(a,'r');
quiver(0,0,real(a),imag(a),1)
r = DesignONeill(1)
DesignONeill
r = DesignONeill(1)
quiver(0,0,real(r),imag(r),1)
quiver(zeros(8,1),zeros(8,1),real(r),imag(r),1)
for k=1:8; subplot(2,4,k);for L=1:8; a=r(L,K);quiver(0,0,real(a),imag(a),1);hold on;end;end
for k=1:8; subplot(2,4,k);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);hold on;end;end
for k=5:8; subplot(2,2,k);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);hold on;end;end
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);hold on;end;end
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);hold on;xlim([-0.6 0.6]);end;end
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a),imag(a),num2str(L));hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
%-- 2020/11/27 7:31 --%
r = DesignONeill(1)
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
r = DesignONeill(1)
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
r = DesignONeill(1)
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
r = DesignONeill(1)
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
r = DesignONeill(1)
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
r = DesignONeill(1)
figure(2)
for k=5:8; subplot(2,2,k-4);for L=1:8; a=r(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;end
r = DesignONeill(1)
[r T]= DesignONeill(1)
T*r(:,1)
T*r(:,4)
T*r(:,5)
[0; 1; (-1+i*sqrt(3))/2; (-1+i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1+i*sqrt(3))/2*i;]
T*[0; 1; (-1+i*sqrt(3))/2; (-1+i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1+i*sqrt(3))/2*i;]
rw=[0; 1; (-1+i*sqrt(3))/2; (-1+i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1+i*sqrt(3))/2*i;];for L=1:8; a=rw(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);
end
figure(3);rw=[0; 1; (-1+i*sqrt(3))/2; (-1+i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1+i*sqrt(3))/2*i;];for L=1:8; a=rw(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);
end
figure(3);rw=[0; 1; (-1+i*sqrt(3))/2; (-1+i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1+i*sqrt(3))/2*i;]*0.5;for L=1:8; a=rw(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
figure(3);rw=[0; 1; (-1+i*sqrt(3))/2; (-1+i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1+i*sqrt(3))/2*i;]*0.5;for L=1:8; a=rw(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
close;figure(3);rw=[0; 1; (-1+i*sqrt(3))/2; (-1-i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1+-*sqrt(3))/2*i;]*0.5;for L=1:8; a=rw(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
figure(3);rw=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*0.5;for L=1:8; a=rw(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
rq=T*rw
rw
[rq rw*i]
figure(3);rw=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));for L=1:8; a=rw(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
[T*rw rw]
[T*rw rw*0.4i]
[T*rw-rw*0.4i]
[T*rw+rw*0.4i]
rw1=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]
rw1'
figure(4);rw=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));for L=1:8; a=conj(rw(L));quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
rw1=conj(rw)
[T*rw1 rw1*0.4i T*rw1-rw1*0.4i]
rwd2=[0; 1;  (-1 -i*sqrt(3))/2;(-1+i*sqrt(3))/2;
0; 1*i;  (-1-i*sqrt(3))/2*i;(-1+i*sqrt(3))/2 *i;]*(1/sqrt(6));
figure(5); ;for L=1:8; a=conj(rwd2(L));quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
[T*rwd2 rwd2*0.4i T*rwd2-rwd2*0.4i]
[T*rwd2 rwd2*0.4i T*rwd2+rwd2*0.4i]
T
rwd2=[0; 1;  (-1 -i*sqrt(3))/2;(-1+i*sqrt(3))/2;
0; 1*i;  (-1-i*sqrt(3))/2*i;(-1+i*sqrt(3))/2 *i;]*(1/sqrt(6));
[rq rw rw1 rw2]
R= [rq rw rw1 rwd2]
R= [rw rw1 rwd2]
figure(3);rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));for L=1:8; a=rwA(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
rwB=conj(rwA);
figure(4);rwB=conj(rwA);for L=1:8; a=rwA(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
figure(3);rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));for L=1:8; a=rwA(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
clear
figure(3);clf;rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));for L=1:8; a=rwA(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end
ID=4;figure(ID);clf;rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));for L=1:8; a=rwA(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(ID))
ID=4;figure(ID);clf;rwB=conj(rwA)
for L=1:8; a=rwB(L);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',20);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(ID))
rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6)); rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 =[rwA rwB rwC rwD];
for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
R4 = [rwD rwC rwB rwA ]
for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
clear
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
R4
sum(abs(R4))
2.4495*2.4495
abs(R4)
abs(rwA)
abs(rwA/sqrt(6))
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/6);rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
sum(abs(R4))
R4
real
factor(rwA)
clf;rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/sqrt(6));
rwA
[abs(rwA) angle(rwA)]
1/sqrt(6)
[abs(rwA) angle(rwA)/pi]
R4
abs(R4)
S = sym(rwA)
S = sym(R4)
latex(S)
%-- 2020/11/27 16:44 --%
ev = [[-0.00000000000000 - 0.612400000000000i,0.00000000000000 + 0.612400000000000i,-0.755900000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.535300000000000i,-0.00000000000000 - 0.535300000000000i,0.155700000000000 - 0.0205000000000000i,0.155700000000000 + 0.0205000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,-0.0425000000000000 - 0.0851000000000000i,-0.0425000000000000 + 0.0851000000000000i,-0.558900000000000 + 0.00000000000000i,-0.558900000000000 + 0.00000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.0425000000000000 - 0.450100000000000i,0.0425000000000000 + 0.450100000000000i,0.403200000000000 + 0.0205000000000000i,0.403200000000000 - 0.0205000000000000i;0.612400000000000 + 0.00000000000000i,0.612400000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.755900000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,0.535300000000000 + 0.00000000000000i,0.535300000000000 + 0.00000000000000i,-0.0205000000000000 - 0.155700000000000i,-0.0205000000000000 + 0.155700000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,-0.0851000000000000 + 0.0425000000000000i,-0.0851000000000000 - 0.0425000000000000i,0.00000000000000 + 0.558900000000000i,-0.00000000000000 - 0.558900000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,-0.450100000000000 - 0.0425000000000000i,-0.450100000000000 + 0.0425000000000000i,0.0205000000000000 - 0.403200000000000i,0.0205000000000000 + 0.403200000000000i]];
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/6);rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
R8=[ev(:,1:2)/sqrt(6)  ev(:,3)/sum(ev(:,3)) ev(:,4)/sum(ev(:,4))  R4]
sum(abs(R8))
lissa2Plot
xtick([])
set(gca,xtick([]))
set(gca,'xtick', [])
set(gca,'xtick', [],'ytick', [])
close
close all
lissa2Plot
close all
lissa2Plot
h=subplot(4,7,figk)
%-- 2020/11/27 18:54 --%
DesignONeill
A43=1;
[r D_Eq_at_NE]= DesignONeill(A43)
lissa2Plot
save('AreaTheo28x8_20201127',AreaTheo)
save('AreaTheo28x8_20201127.mat',AreaTheo)
save('AreaTheo28x8_20201127.mat','AreaTheo')
close all
R8s=round(R8,4)
R8s=round(R8,3)
R8s=round(R8,4)
latex(R8s)
matrix2latex(R8s,'1.tex')
matrix2latex(AreaTheo,'2.tex')
T28s=round(AreaTheo,4)
matrix2latex(T28s,'2.tex')
%-- 2020/11/27 20:33 --%
testsave2Dfig
[r D_Eq_at_NE]= DesignONeill(A43)
A43=1;
[r D_Eq_at_NE]= DesignONeill(A43)
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/6);rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
[D_Eq_at_NE*R8(:,1) R8(:,1)*4i]
R8=[ev(:,1:2)/sqrt(6)  ev(:,3)/sum(ev(:,3)) ev(:,4)/sum(ev(:,4))  R4]
ev = [[-0.00000000000000 - 0.612400000000000i,0.00000000000000 + 0.612400000000000i,-0.755900000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.535300000000000i,-0.00000000000000 - 0.535300000000000i,0.155700000000000 - 0.0205000000000000i,0.155700000000000 + 0.0205000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,-0.0425000000000000 - 0.0851000000000000i,-0.0425000000000000 + 0.0851000000000000i,-0.558900000000000 + 0.00000000000000i,-0.558900000000000 + 0.00000000000000i;0.00000000000000 + 0.204100000000000i,0.00000000000000 - 0.204100000000000i,-0.378000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.0425000000000000 - 0.450100000000000i,0.0425000000000000 + 0.450100000000000i,0.403200000000000 + 0.0205000000000000i,0.403200000000000 - 0.0205000000000000i;0.612400000000000 + 0.00000000000000i,0.612400000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.755900000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,0.535300000000000 + 0.00000000000000i,0.535300000000000 + 0.00000000000000i,-0.0205000000000000 - 0.155700000000000i,-0.0205000000000000 + 0.155700000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,-0.0851000000000000 + 0.0425000000000000i,-0.0851000000000000 - 0.0425000000000000i,0.00000000000000 + 0.558900000000000i,-0.00000000000000 - 0.558900000000000i;-0.204100000000000 + 0.00000000000000i,-0.204100000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.378000000000000 + 0.00000000000000i,-0.450100000000000 - 0.0425000000000000i,-0.450100000000000 + 0.0425000000000000i,0.0205000000000000 - 0.403200000000000i,0.0205000000000000 + 0.403200000000000i]];
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/6);rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
R8=[ev(:,1:2)/sqrt(6)  ev(:,3)/sum(ev(:,3)) ev(:,4)/sum(ev(:,4))  R4]
sum(abs(R8))
[D_Eq_at_NE*R8(:,5)-R8(:,5)*4i]
[D_Eq_at_NE*R8(:,5)-R8(:,5)*0.4i]
[D_Eq_at_NE*R8(:,7)-R8(:,7)*0.4i]
save('TheoEigCycle','R8','D_Eq_at_NE','AreaTheo')
%-- 2020/11/28 1:32 --%
DesignONeill
DesignONeill(1)
[S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]
DesignONeill(1)
%-- 2020/11/28 11:23 --%
load('TheoEigCycle.mat')
listID = [];for m=1:7;for n = m:8; listID = [listID;m*10+n m n];end;end
save('TheoEigCycle.mat','listID','-append')
Th=[AreaTheo listID];
for m=1:7;for n = m+1:8; listID = [listID;m*10+n m n];end;end
save('TheoEigCycle.mat','listID','-append')
load('TheoEigCycle.mat')
for m=1:7;for n = m+1:8; listID = [listID;m*10+n m n];end;end
listID = [];for m=1:7;for n = m+1:8; listID = [listID;m*10+n m n];end;end
save('TheoEigCycle.mat','listID','-append')
load('TheoEigCycle.mat')
c=round(AreaTheo,3)
c=round(AreaTheo,4)
c=[c listID]
f = latex2MxWithMxPrecision(c, 4)
f = latex2MxWithMxPrecision(c(:,[0 1:8]), 4)
f = latex2MxWithMxPrecision(c(:,[9 1:8]), 4)
%-- 2020/11/28 13:29 --%
load('TheoEigCycle.mat')
f = latex2MxWithMxPrecision(R8, 4)
f = latex2MxWithMxPrecision(Imag(R8), 4)
f = latex2MxWithMxPrecision(imag(R8), 4)
%-- 2020/11/28 20:18 --%
testsave2Dfig
data6expmean = data(:,7:12);
HP =[];for j=1:28; [h p] = ttest(data6expmean(j,:)); HP =[HP;j h p];end
HP(:.4) = data(:,1)
HP(:,4) = data(:,1)
%-- 2020/11/29 20:53 --%
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/6);rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
%-- 2020/11/29 23:21 --%
load('TheoEigCycle.mat')
testsave2Dfig
rwA=[0; 1; (-1+i*sqrt(3))/2; (-1 -i*sqrt(3))/2; 0; 1*i; (-1+i*sqrt(3))/2 *i; (-1-i*sqrt(3))/2*i;]*(1/6);rwB=conj(rwA);rwC=[rwA(1); rwA(2); rwA(4); rwA(3); rwA(5); rwA(6); rwA(8); rwA(7)];rwD=conj(rwC);
R4 = [rwD rwC rwB rwA ];for k=1:4; subplot(2,2,k);for L=1:8; a=R4(L,k);quiver(0,0,real(a),imag(a),1);text(real(a)*1.1,imag(a)*1.1,num2str(L),'fontsize',10);hold on;xlim([-0.6 0.6]);ylim([-0.6 0.6]);end;title(num2str(k));end
load('TheoEigCycle.mat')
sum(AreaTheo)
R8 = [0.00000000000000 - 0.250011253080070i,0.00000000000000 + 0.250011253080070i,0.399968252288481 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;0.00000000000000 + 0.0833234760836745i,0.00000000000000 - 0.0833234760836745i,0.200010582570506 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.166666666666667 + 0.00000000000000i,0.166666666666667 + 0.00000000000000i,0.166666666666667 + 0.00000000000000i,0.166666666666667 + 0.00000000000000i;0.00000000000000 + 0.0833234760836745i,0.00000000000000 - 0.0833234760836745i,0.200010582570506 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,-0.0833333333333333 + 0.144337567297406i,-0.0833333333333333 - 0.144337567297406i,-0.0833333333333333 - 0.144337567297406i,-0.0833333333333333 + 0.144337567297406i;0.00000000000000 + 0.0833234760836745i,0.00000000000000 - 0.0833234760836745i,0.200010582570506 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,-0.0833333333333333 - 0.144337567297406i,-0.0833333333333333 + 0.144337567297406i,-0.0833333333333333 + 0.144337567297406i,-0.0833333333333333 - 0.144337567297406i;0.250011253080070 + 0.00000000000000i,0.250011253080070 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.399968252288481 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i;-0.0833234760836745 + 0.00000000000000i,-0.0833234760836745 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.200010582570506 + 0.00000000000000i,0.00000000000000 - 0.166666666666667i,0.00000000000000 + 0.166666666666667i,0.00000000000000 - 0.166666666666667i,0.00000000000000 + 0.166666666666667i;-0.0833234760836745 + 0.00000000000000i,-0.0833234760836745 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.200010582570506 + 0.00000000000000i,0.144337567297406 + 0.0833333333333333i,0.144337567297406 - 0.0833333333333333i,-0.144337567297406 + 0.0833333333333333i,-0.144337567297406 - 0.0833333333333333i;-0.0833234760836745 + 0.00000000000000i,-0.0833234760836745 + 0.00000000000000i,0.00000000000000 + 0.00000000000000i,0.200010582570506 + 0.00000000000000i,-0.144337567297406 + 0.0833333333333333i,-0.144337567297406 - 0.0833333333333333i,0.144337567297406 + 0.0833333333333333i,0.144337567297406 - 0.0833333333333333i]% close all
N41 = (R8(:,5)+R8(:,7)*i)/2;N42 = (R8(:,5)-R8(:,7)*i)/2; N4 = [N41 N42];
%-- 2020/12/3 1:06 --%
testsave2Dfig
ax2 = scatter(data(:,TheoryCol),data(:,Yfig_id), 1, 'MarkerEdgeColor','b')
testsave2Dfig
line([0,0], [MinMax(Yfig_id,1),MinMax(Yfig_id,2)]*1.2,'LineWidth',2,  'color',[0.8 0.8 0.8]);hold on
testsave2Dfig
%-- 2020/12/3 1:48 --%
testsave2Dfig
r= []; for J=7:12; [h p] = ttest(data(:,J)); r=[r; J h p];end
r= []; for J=7:12; [h p] = ttest(data(:,J)); r=[r; J h p mean(data(:,J))];end
P11=mean(data(:,[7 9:12])')'
P11=mean(data(:,[7:12])')'
P11=mean(data(:,[7 9:12])')'
[h p] = ttest(P11)
P11=mean(data(:,[7:12])')'
[h p] = ttest(P11)
testsave2Dfig
scatter(getcolumn(data(:,[6,11]),1),getcolumn(data(:,[6,11]),2))
Oneill_replicator
length(A(:,1))
Oneill_replicator
D_Eq_at_NE
%-- 2020/12/6 19:28 --%
readYdata1206
r = readYdata1206
mean(r(:,4:11))
0.1865/0.2
readYdata1206
Lmn_t0t1 = [Lmn_t0t1; cross(v_t0,v_t1)]
Lmn_t0t1 = [Lmn_t0t1; cross([v_t0 0],[v_t1 0])]
readYdata1206
mean(Lmn_t0t1)'
mean_L_15327 = mean(Lmn_t0t1)';
plot(mean_L_15327)
stem(mean_L_15327)
save('F:\Spectrum2020\datasource\L_28.mat','L28_all15326','L28_P11')
readYdata1206
%-- 2020/12/6 22:28 --%
load('F:\Spectrum2020\datasource\L_28.mat')
unique(L28_all15326(:,[29]))
unique(L28_all15326(:,[29 30]))
unique(L28_all15326(:,[30]))
unique(L28_all15326(:,[29]))
[a b]=unique(L28_all15326(:,[29 30]))
[a b c]=unique(L28_all15326(:,[29 30]))
clear
load('F:\Spectrum2020\datasource\L_28.mat')
for I=1:6; tmpA=L28_all15326(find( L28_all15326(:,29)==I),:);L28_YTable4(:,I) = mean(tmpA)';end
mn28=[];for m=1:7;for n=m+1:8; mn28=[mn28;m*10+n];end
end
YTable4_expL = [mn28 L28_YTable4(1:28,:)]
corr(L28_YTable4(1:28,:))
[YTable_Spearman p]= corr(L28_YTable4(1:28,:),'type', 'spearman')
[YTable5_Spearman p]= corr(L28_YTable4(1:28,:),'type', 'spearman')
mn28=[];for m=1:7;for n=m+1:8; mn28=[mn28;m*10+n];end;end;for I=1:6; tmpA=L28_all15326(find( L28_all15326(:,29)==I),:);L28_YTable4(:,I) = mean(tmpA)';end; YTable4_expL = [mn28 L28_YTable4(1:28,:)]
%-- 2020/12/11 0:09 --%
Mix5_2000
N_cycle0315_Price6
[S.x1 S.x2 S.x3 S.x4 S.x5 S.x6]
x1=A(length(A(:,1)),1); x2=A(length(A(:,1)),2);
x3=A(length(A(:,1)),3);
x4=A(length(A(:,1)),4); x5=A(length(A(:,1)),5);
x6=A(length(A(:,1)),6);
[eigen_vector eigen_value] = eig(eval(D_V_F))
A=[S.x1 S.x2 S.x3 S.x4 S.x5 S.x6]
x1=A(length(A(:,1)),1); x2=A(length(A(:,1)),2);
x3=A(length(A(:,1)),3);
x4=A(length(A(:,1)),4); x5=A(length(A(:,1)),5);
x6=A(length(A(:,1)),6);
[eigen_vector eigen_value] = eig(eval(D_V_F))
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5 S.x6])
x1=A(length(A(:,1)),1); x2=A(length(A(:,1)),2);
x3=A(length(A(:,1)),3);
x4=A(length(A(:,1)),4); x5=A(length(A(:,1)),5);
x6=A(length(A(:,1)),6);
[eigen_vector eigen_value] = eig(eval(D_V_F))
%-- 2020/12/11 12:44 --%
Mix5_2000
%-- 2020/12/11 12:49 --%
Mix5_2000
sum(Area)
Area=[];
for J=1:6;
tmpArea=[];
B = eigen_vector(:,J);
for m=1:5;
for n=m+1:6;
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn) ;
tmpArea=[tmpArea ; area_m_n];
end;
end
Area=[Area tmpArea];
end
Area
exp15 =[17
4
-7
9
-6
-3
2
7
6
-6
-10
15
-2
-10
-6
]
C=[exp15 Area]
corr(C)
[h p] = corr(C(:,1),C(:,5),'type','spearman')
tmpArea=[];for m=1:5;
for n=m+1:6;;
m_n =  m*10+n;
tmpArea=[tmpArea ;  m_n];
end;
end
tmpArea
exp15=[17	4	-7	9	-6	-3	2	7	6	-6	-10	15	-2	-10	-6]';
data=[exp15 Area];
y=data(:,1); xx1=data(:,3);xx2=data(:,5); %wi8i4
X = [ones(size(y))  xx1 xx2];
[b,bint,r,rint,stats] = regress(y,X);bb= [b bint stats']; rr=[rr; bb];
exp15=[17	4	-7	9	-6	-3	2	7	6	-6	-10	15	-2	-10	-6]';
data=[exp15 Area];
y=data(:,1); xx1=data(:,3);xx2=data(:,5); %wi8i4
X = [ones(size(y))  xx1 xx2];
[b,bint,r,rint,stats] = regress(y,X);
[b,bint,r,rint,stats] = regress(y,X)
scatter3(xx1,xx2,y,'filled')
exp15 =  [17	4	-7	9	6	-3	2	7	-6	-6	-10	-15	-2	10	6]';  % dimension 6 oppose
data=[exp15 Area];
y=data(:,1); xx1=data(:,3);xx2=data(:,5); %wi8i4
X = [ones(size(y))  xx1 xx2];
[b,bint,r,rint,stats] = regress(y,X);
[b,bint,r,rint,stats] = regress(y,X)
%-- 2020/12/11 16:22 --%
Mix5_2000
[b,bint,r,rint,stats] = regress(y,X)
exp15 = [17	4	-7	9	-6	-3	2	7	6	-6	-10	15	-2	-10	-6]';
%exp15 = [17	4	-7	9	6	-3	2	7	-6	-6	-10	-15	-2	 10	 6]';  % dimension 6 oppose
data=[exp15 Area];
y=data(:,1); xx1=data(:,3);xx2=data(:,5); %wi8i4
X = [ones(size(y))  xx1 xx2];
[b,bint,r,rint,stats] = regress(y,X);
[b,bint,r,rint,stats] = regress(y,X)
scatter(getcolumn(data(:,[1,end]),1),getcolumn(data(:,[1,end]),2))
Mix5_2000
plot(eigen_vector(:,2))
plot(eigen_vector(:,4))
quiver(real(eigen_vector(:,2)),imag(eigen_vector(:,2)),1)
quiver(0,0,real(eigen_vector(:,2)),imag(eigen_vector(:,2)),1)
quiver(zero(6,1),zero(6,1),real(eigen_vector(:,2)),imag(eigen_vector(:,2)),1)
quiver(zeros(6,1),zeros(6,1),real(eigen_vector(:,2)),imag(eigen_vector(:,2)),1)
[abs(eigen_vector(:,2)) angle(eigen_vector(:,2))]
for tJ=1:6
quiver(0,0,real(eigen_vector(tJ,2)),imag(eigen_vector(tJ,2)),1);hold on
text(real(eigen_vector(tJ,2)),imag(eigen_vector(tJ,2)),num2str(eigen_vector(tJ,2)));hold on
end
for tI=1:6; subplot(2,3,tI);title(num2str(eigen_value(tI,tI)))
for tJ=1:6
quiver(0,0,real(eigen_vector(tJ,2)),imag(eigen_vector(tJ,2)),1);hold on
text(real(eigen_vector(tJ,2)),imag(eigen_vector(tJ,2)),num2str(eigen_vector(tJ,2)));hold on
end
end
for tI=1:6; subplot(2,3,tI);title(num2str(eigen_value(tI,tI)))
for tJ=1:6
quiver(0,0,real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),1);hold on
text(real(eigen_vector(tJ,tI)),imag(eigen_vector(tJtI2)),num2str(eigen_vector(tJ,tI)));hold on
end
end
for tI=1:6; subplot(2,3,tI);title(num2str(eigen_value(tI,tI)))
for tJ=1:6
quiver(0,0,real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),1);hold on
text(real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),num2str(eigen_vector(tJ,tI)));hold on
end
end
for tI=1:6; subplot(2,3,tI);title(num2str(eigen_value(tI,tI)))
for tJ=1:6
quiver(0,0,real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),1,'LineWidth',tJ);hold on
text(real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),num2str(eigen_vector(tJ,tI)));hold on
end
end
Mix5_2000
%% 观察特征向量的几何形态
for tI=1:6; subplot(2,3,tI);title(num2str(eigen_value(tI,tI)));axis square;hold on
for tJ=1:6
quiver(0,0,real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),1,'LineWidth',tJ);hold on
text(real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),num2str(eigen_vector(tJ,tI)));hold on
end
end
[abs(eigen_vector(:,2)) angle(eigen_vector(:,2))];
Mix5_2000
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price6\abed-1pop Strategy distributions (complete history)a.csv'
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22]
mean(ps1)
mean(ps1(:,2:7))
psd = [v.u1-v.u2 v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u22];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22];
mean(psd)
[h p]=ttest(psd(:,6),psd(:,5))
[h p]=ttest(psd(:,1),psd(:,5))
[h p]=ttest(psd(:,1),psd(:,6))
[h p]=ttest(psd(:,1),psd(:,2))
[h p]=ttest(psd(2001:2500,1),psd(2001:2500,2))
[h p]=ttest(psd(1001:1500,1),psd(1001:1500,2))
[h p]=ttest(psd(100:500,1),psd(100:500,2))
plot(psd(:,1:2),'DisplayName','psd(:,1:2)')
plot(psd,'DisplayName','psd')
bar(psd,'DisplayName','psd')
histogram(psd)
bar(psd,'DisplayName','psd')
angluar_momentum = from_N_colExp_out_am(psd)
data=[angluar_momentum Area];
y=data(:,1); xx1=data(:,3);xx2=data(:,5); %wi8i4
X = [ones(size(y))  xx1 xx2];
[b,bint,r,rint,stats] = regress(y,X)
Mix5_2000
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price6\abed-1pop Strategy distributions (complete history)a.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22];
abed_mean=mean(psd)
abed_angluar_momentum = from_N_colExp_out_am(psd);
data=[abed_angluar_momentum theo_eigcyc_Area];
y=data(:,1); xx1=data(:,3);xx2=data(:,5); %wi8i4
X = [ones(size(y))  xx1 xx2]
[b,bint,r,rint,stats] = regress(y,X)
Mix5_2000
f = latex2MxWithMxPrecision(payoff_matrix, 0)
diag(eigen_value)
facotr(diag(eigen_value))
factor(diag(eigen_value))
latex2MxWithMxPrecision(diag(eigen_value),4)
latex2MxWithMxPrecision(diag(eigen_value)',4)
latex(diag(eigen_value)')
latex(num2str(diag(eigen_value)'))
(num2str(diag(eigen_value)'))
-0.2886*3
0.28868*3
matrix2latex(diag(eigen_value)')
latex2MxWithMxPrecision(diag(eigen_vector)',4)
latex2MxWithMxPrecision((eigen_vector)',4)
roundn(theo_eigcyc_Area,4)
round(theo_eigcyc_Area,4)
latex2MxWithMxPrecision(theo_eigcyc_Area,4)
tA=[];for m = 1:5; for n = m+1:6; tA=[tA;m*10+n];end;end
latex2MxWithMxPrecision([tA theo_eigcyc_Area],4)
a=num2str(eigen_vector)
latex(a)
save('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212.mat')
abedvsTheo = [abed_angluar_momentum theo_eigcyc_Area]
scatter(getcolumn(abedvsTheo(:,[1,4]),1),getcolumn(abedvsTheo(:,[1,4]),2))
scatter(getcolumn(abedvsTheo(:,[1,6]),1),getcolumn(abedvsTheo(:,[1,6]),2))
scatter(getcolumn(abedvsTheo(:,[1,6]),1),getcolumn(abedvsTheo(:,[1,6]),2),'ro')
scatter(getcolumn(abedvsTheo(:,[1,6]),1),getcolumn(abedvsTheo(:,[1,6]),2),'ro');hold on
scatter(getcolumn(abedvsTheo(:,[1,4]),1),getcolumn(abedvsTheo(:,[1,4]),2),'bX')
plot(abedvsTheo(:,1))
%-- 2020/12/12 11:41 --%
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price6\abed-1pop Strategy distributions (complete history)b.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22];
abed_mean=mean(psd)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price6\abed-1pop Strategy distributions (complete history)b';
[num,txt,raw] = xlsread(abedfilesname);
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price6\abed-1pop Strategy distributions (complete history)b.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22];
abed_mean=mean(psd)
Mix5_2000
Mix5_2000
2.3199/0.7622
sqrt(3)/6
0.4534^2
save('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212c.mat')
Mix5_2000
clear
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212c.mat')
0.4534/pi
sqrt(3)/12
sqrt(3)*pi/12
drop u
sldrop u
u=[];
TE = [theo_eigcyc_Area abed_angluar_momentum]
close all
Mix5_2000
close all
Mix5_2000
abed_theo_mn=[abed_angluar_momentum theo_eigcyc_Area theo_subspace_id_mn]
y=abed_theo_mn(:,1); xx1=abed_theo_mn(:,3);xx2=abed_theo_mn(:,5); %wi8i4
X = [ones(size(y))  xx1 xx2]
[b,bint,r,rint,stats] = regress(y,X)
Mix5_2000
clear
close all
clear
Mix5_2000
y=theo_abed_mn(:,7); xx1=theo_abed_mn(:,2);xx2=theo_abed_mn(:,4); %wi8i4
X = [ones(size(y))  xx1 xx2]
[b,bint,r,rint,stats] = regress(y,X)
close all
clear
Mix5_2000
spectrumArr = [theo_abed_mn TE_multiLinerE];
spectrumArr = [theo_abed_mn TE_multiLinerE roundn(TE_multiLinerE*10)];
spectrumArr = [theo_abed_mn TE_multiLinerE roundn(TE_multiLinerE*10,0)];
spectrumArr = [theo_abed_mn TE_multiLinerE roundn(theo_abed_mn(:,2)*10,0) roundn(TE_multiLinerE*10,0)];
spectrumArr = [theo_abed_mn TE_multiLinerE roundn(theo_abed_mn(:,2)*1000,0) roundn(TE_multiLinerE*1000,0)];
tUni12 = unique(spectrumArr(:,12))
for tI=1:length(tUni12); spectrumArr(find(spectrumArr(:,12)==tUni12(tI)),8) end
spectrumArr = [theo_abed_mn TE_multiLinerE roundn(theo_abed_mn(:,2)*1000,0) roundn(TE_multiLinerE*1000,0)];
tUni12 = unique(spectrumArr(:,12))
for tI=1:length(tUni12);
spectrumArr(find(spectrumArr(:,12)==tUni12(tI)),8)
end
spectrumArr = [theo_abed_mn TE_multiLinerE roundn(theo_abed_mn(:,2)*1000,0) roundn(TE_multiLinerE*1000,0)];
tUni12 = unique(spectrumArr(:,12))
for tI=1:length(tUni12);
[ spectrumArr(find(spectrumArr(:,12)==tUni12(tI)),8)' tUni12(tI)]
end
tUni13 = unique(spectrumArr(:,13))
for tI=1:length(tUni13);
[ spectrumArr(find(spectrumArr(:,12)==tUni13(tI)),8)' tUni13(tI)]
end
tUni13 = unique(spectrumArr(:,13))
for tI=1:length(tUni13);
[ spectrumArr(find(spectrumArr(:,13)==tUni13(tI)),8)' tUni13(tI)]
end
Mix8_2000
[eigen_vector eigen_value] = eig(eval(D_V_F))
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5 S.x6 S.x7 S.x8]);
for k=1:length(A(:,1))
if prod(A(k,:)) > 0
x1=A(k,1); x2=A(k,2); x3=A(k,3); x4=A(k,4); x5=A(k,5); x6=A(k,6); x7=A(k,7); x8=A(k,8);
Rest_Ne = A(k,:);
else
'Not inner mixed strategy';
end
end
[eigen_vector eigen_value] = eig(eval(D_V_F))
for tI=1:8; subplot(2,4,tI);title(num2str(eigen_value(tI,tI)));axis square;hold on
for tJ=1:8
quiver(0,0,real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),1,'LineWidth',tJ);hold on
text(real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),num2str(eigen_vector(tJ,tI)));hold on
end
end
theo_eigcyc_Area=[];
for J=1:8;
theo_subspace_id_mn=zeros(28,3); tmpK=1;
tmpArea=[];
B = eigen_vector(:,J);
for m=1:7;
for n=m+1:8;
theo_subspace_id_mn(tmpK,:)=[m*10+n m n];tmpK=tmpK+1;
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn);
tmpArea=[tmpArea; area_m_n];
end;
end
theo_eigcyc_Area=[theo_eigcyc_Area tmpArea];
end
theo_eigcyc_Area
save('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price8-20201212c.mat')
f = latex2MxWithMxPrecision( theo_eigcyc_Area, precision)
f = latex2MxWithMxPrecision( theo_eigcyc_Area, 3)
f = latex2MxWithMxPrecision( theo_eigcyc_Area, 4)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price6\abed-1pop Strategy distributions (complete history)8a.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22 v.u26 v.u30];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30];
abed_mean=mean(psd)
NE = Rest_Ne
abed_angluar_momentum = from_N_colExp_out_am(psd,NE)/length(psd(:,1))-1;
abed_angluar_momentum = from_N_colExp_out_am(psd,NE)/(length(psd(:,1))-1);
theo_abed_mn=[ theo_eigcyc_Area abed_angluar_momentum theo_subspace_id_mn]
f = latex2MxWithMxPrecision( theo_eigcyc_Area, 4)
%-- 2020/12/13 18:13 --%
save('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212c.mat')
load
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212c.mat')
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212.mat')
medianPrice = median(psd')';
iqr6 = iqr(psd')';
med_iqr = [medianPrice iqr6];
plot(med_iqr,'DisplayName','med_iqr')
scatter(med_iqr(:,1),med_iqr(:,2))
medianPrice = ones(2588,1);
for tI=1:2588; for tJ = 6:-1:2;
tm = 7-tJ;
if ps1(tI,tJ) > 0.5; medianPrice(tI) =  7-tJ; end; end
end
medianPrice = ones(2588,1);
for tI=1:2588;
medianPrice(tI)=max(find(ps1(tI,:)> 0.5));
end
iqr6 = iqr(psd')';
med_iqr = [medianPrice iqr6];
scatter(med_iqr(:,1),med_iqr(:,2))
corr(med_iqr(:,1),med_iqr(:,2))
plot(med_iqr,'DisplayName','med_iqr')
clear
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price10-20201212.mat')
totalperiods = 5232;
medianPrice = ones(totalperiods,1);
for tI=1:totalperiods;
medianPrice(tI)=max(find(ps1(tI,:)> 0.5));
end
iqr6 = iqr(psd')';
med_iqr = [medianPrice iqr6];
scatter(med_iqr(:,1),med_iqr(:,2))
plot(med_iqr,'DisplayName','med_iqr')
scatter(med_iqr(:,1),med_iqr(:,2),'--')
scatter(med_iqr(:,1),med_iqr(:,2),'-')
scatter(med_iqr(:,1),med_iqr(:,2),'*')
line(med_iqr(:,1),med_iqr(:,2),'*')
line(med_iqr(:,1),med_iqr(:,2))
NC = mean(med_iqr)
med_iqr = [medianPrice iqr6];
NC = mean(med_iqr) ;
arr= [];
for tI = 1:totalperiods-1
arr=arr;cross([med_iqr(tI,:)-NC  0],[med_iqr(tI,:)-NC 0]);
end
arr= [];
for tI = 1:totalperiods-1
arr=[arr;cross([med_iqr(tI,:)-NC  0],[med_iqr(tI,:)-NC 0])];
end
arr= [];
for tI = 1:totalperiods-1
arr=[arr;cross([med_iqr(tI,:)-NC  0],[med_iqr(tI+1,:)-NC 0])];
end
mean(arr)
[h p]=ttest(arr(:,3))
arr= [];
for tI = 1:totalperiods-1
arr=[arr;cross([(med_iqr(tI,:)-NC)  0],[(med_iqr(tI+1,:)-NC) 0])];
end
[h p]=ttest(arr(:,3))
arr= [];
for tI = 1:totalperiods-1000
arr=[arr;cross([(med_iqr(tI,:)-NC)  0],[(med_iqr(tI+1,:)-NC) 0])];
end
[h p]=ttest(arr(:,3))
arr= [];
for tI = 1:totalperiods-3000
arr=[arr;cross([(med_iqr(tI,:)-NC)  0],[(med_iqr(tI+1,:)-NC) 0])];
end
[h p]=ttest(arr(:,3))
arr= [];
for tI = 1111:totalperiods-1
arr=[arr;cross([(med_iqr(tI,:)-NC)  0],[(med_iqr(tI+1,:)-NC) 0])];
end
[h p]=ttest(arr(:,3))
arr= [];
for tI = 2111:totalperiods-1
arr=[arr;cross([(med_iqr(tI,:)-NC)  0],[(med_iqr(tI+1,:)-NC) 0])];
end
[h p]=ttest(arr(:,3))
histc(arr(:,3))
histc(arr(:,3),10)
histc(arr(:,3),2)
arr(:,3)
clear
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212c.mat')
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212.mat')
eta1 = [3+7i, 2+5i];  area_i_j1  = lissa2Complex2areaPlot(eta1(1),eta1(2))
eta2 = [3+1.7i, 1.2+5i];  area_i_j2  = lissa2Complex2areaPlot(eta2(1),eta2(2))
eta3=eta1+eta2
area_i_j3 = lissa2Complex2areaPlot(eta3(1),eta3(2))
syms a1 b1 a2 b2 c1 c2 d1 d2
cross(a1+b1*i, a2+b2*i)
a1*a2 - b1*b2 + (a1*b2 + a2*b1)*i
c1*c2 - d1*d2 + (c1*d2 + c2*d1)*i
vt=eigen_vector
v2=0.5(vt(:,2)+vt(:,4))
v2 = 0.5*(vt(:,2) + vt(:,4))
tmpArea=[];
B = v2;
for m=1:6-1;
for n=m+1:10;
theo_subspace_id_mn(tmpK,:)=[m*100+n m n];tmpK=tmpK+1;
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn);
tmpArea=[tmpArea; area_m_n];
end;
end
tmpArea=[];
B = v2;
for m=1:6-1;
for n=m+1:10;
theo_subspace_id_mn(tmpK,:)=[m*100+n m n];
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn);
tmpArea=[tmpArea; area_m_n];
end;
end
tmpArea=[];
B = v2;
for m=1:6-1;
for n=m+1:10;
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn);
tmpArea=[tmpArea; area_m_n];
end;
end
tmpArea=[];
B = v2;
for m=1:6-1;
for n=m+1:6;
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn);
tmpArea=[tmpArea; area_m_n];
end;
end
%-- 2020/12/13 22:47 --%
Mix10_2000
payoffA
bimat(A,A')
bimat(payoffA,payoffA')
D_coff_ai_fi_price10
[payoffA mean_U V_F D_V_F Jacob_at_Ne eigen_vector eigen_value theo_eigcyc_Area] = D_coff_ai_fi_price10(e1,r)
eigen_value_d = diag(eigen_value)
1.5388*1.1
D_coff_ai_fi_price10
[payoffA mean_U V_F D_V_F Jacob_at_Ne eigen_vector eigen_value_d theo_eigcyc_Area] = D_coff_ai_fi_price10(e1,r)
D_coff_ai_fi_price10
[biMatA,biMatB,biMata,biMatb,biMatiterations,biMaterr,biMatms]=bimat(payoffA,payoffA')
payoffA
payoff_matrix-2
payoff_matrix-2*ones(10)
payoffA-2*ones(10)
B=payoffA-2*ones(10);  [biMatA,biMatB,biMata,biMatb,biMatiterations,biMaterr,biMatms]=bimat(B,B')
D_coff_ai_fi_price09
B=payoffA-2*ones(10);  [biMatA,biMatB,biMata,biMatb,biMatiterations,biMaterr,biMatms]=bimat(B,B')
B=payoffA-2*ones(9);  [biMatA,biMatB,biMata,biMatb,biMatiterations,biMaterr,biMatms]=bimat(B,B')
D_coff_ai_fi_price09
V_F
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5 S.x6 S.x7 S.x8 S.x9])
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5 S.x6 S.x7 S.x8 S.x9]);
for k=1:length(A(:,1))
if prod(A(k,:)) > 0
x1=A(k,1); x2=A(k,2); x3=A(k,3); x4=A(k,4); x5=A(k,5); x6=A(k,6); x7=A(k,7); x8=A(k,8); x9=A(k,9);
Rest_Ne = A(k,:);
else
'Not inner mixed strategy';
end
end
Rest_Ne
payoffA
D_coff_ai_fi_price09
D_coff_ai_fi_price07
clear
A=[0.128	0.334	0.481	0.153
0.124	0.316	0.497	0.156
0.116	0.303	0.503	0.165
0.104	0.292	0.509	0.173
0.1	0.283	0.513	0.179
0.1	0.271	0.513	0.191
0.1	0.255	0.515	0.205
0.096	0.244	0.513	0.219
0.096	0.23	0.517	0.229
0.112	0.23	0.508	0.234
0.112	0.228	0.503	0.241
0.112	0.214	0.503	0.255
0.104	0.21	0.498	0.266
0.108	0.197	0.498	0.278
0.116	0.185	0.483	0.303
0.128	0.173	0.475	0.32
0.132	0.17	0.461	0.336
0.14	0.17	0.452	0.343
0.14	0.172	0.44	0.353
0.152	0.168	0.429	0.365
0.16	0.158	0.423	0.379
0.164	0.154	0.417	0.388
0.192	0.143	0.405	0.404
0.2	0.141	0.401	0.408
0.212	0.137	0.39	0.42
0.232	0.132	0.379	0.431
0.232	0.132	0.366	0.444
0.232	0.127	0.358	0.457
0.252	0.125	0.341	0.471
0.268	0.127	0.333	0.473
0.288	0.13	0.324	0.474
0.284	0.131	0.317	0.481
0.296	0.124	0.303	0.499
0.312	0.122	0.291	0.509
0.344	0.122	0.284	0.508
0.332	0.127	0.279	0.511
0.344	0.133	0.276	0.505
0.36	0.13	0.272	0.508
0.4	0.131	0.26	0.509
0.396	0.138	0.253	0.51
0.416	0.132	0.251	0.513
0.432	0.133	0.238	0.521
0.48	0.14	0.234	0.506
0.476	0.143	0.228	0.51
0.504	0.149	0.224	0.501
0.52	0.155	0.212	0.503
0.54	0.163	0.21	0.492
0.548	0.171	0.206	0.486
0.56	0.183	0.203	0.474
0.568	0.192	0.199	0.467
0.576	0.206	0.195	0.455
0.608	0.207	0.192	0.449
0.636	0.223	0.187	0.431
0.636	0.239	0.189	0.413
0.644	0.254	0.183	0.402
0.656	0.272	0.184	0.38
0.668	0.286	0.182	0.365
0.688	0.292	0.185	0.351
0.684	0.314	0.179	0.336
0.664	0.331	0.179	0.324
0.676	0.338	0.179	0.314
0.64	0.363	0.174	0.303
0.612	0.385	0.18	0.282
0.604	0.408	0.18	0.261
0.6	0.424	0.181	0.245
0.572	0.444	0.182	0.231
0.528	0.459	0.192	0.217
0.508	0.472	0.192	0.209
0.488	0.482	0.195	0.201
0.504	0.485	0.202	0.187
0.496	0.487	0.206	0.183
0.484	0.503	0.202	0.174
0.46	0.51	0.206	0.169
0.432	0.517	0.217	0.158
0.404	0.529	0.222	0.148
0.392	0.535	0.228	0.139
0.388	0.523	0.243	0.137
0.372	0.525	0.247	0.135
0.368	0.521	0.254	0.133
0.348	0.527	0.263	0.123
0.328	0.525	0.273	0.12
0.304	0.537	0.265	0.122
0.3	0.539	0.27	0.116
0.296	0.529	0.281	0.116
0.288	0.523	0.29	0.115
0.272	0.519	0.301	0.112
0.268	0.516	0.307	0.11
0.252	0.508	0.321	0.108
0.232	0.507	0.324	0.111
0.208	0.495	0.34	0.113
0.192	0.483	0.353	0.116
0.184	0.469	0.367	0.118
0.168	0.452	0.392	0.114
0.152	0.436	0.413	0.113
0.16	0.423	0.421	0.116
0.152	0.404	0.438	0.12
0.144	0.397	0.442	0.125
0.14	0.379	0.451	0.135
0.14	0.368	0.454	0.143
0.14	0.341	0.468	0.156
0.128	0.322	0.48	0.166
0.136	0.318	0.476	0.172
0.132	0.302	0.483	0.182
0.112	0.291	0.493	0.188
0.104	0.272	0.501	0.201
0.104	0.265	0.496	0.213
0.1	0.256	0.494	0.225
0.096	0.248	0.495	0.233
0.088	0.239	0.491	0.248
0.088	0.231	0.487	0.26
0.088	0.215	0.49	0.273
0.084	0.206	0.482	0.291
0.084	0.201	0.475	0.303
0.084	0.194	0.477	0.308
0.076	0.188	0.476	0.317
0.076	0.185	0.468	0.328
0.076	0.177	0.457	0.347
0.08	0.172	0.454	0.354
0.08	0.165	0.446	0.369
0.084	0.16	0.435	0.384
0.092	0.154	0.422	0.401
0.088	0.146	0.412	0.42
0.08	0.141	0.397	0.442
0.092	0.127	0.387	0.463
0.092	0.127	0.37	0.48
0.096	0.122	0.361	0.493
0.112	0.119	0.348	0.505
0.116	0.115	0.34	0.516
0.116	0.112	0.327	0.532
0.124	0.103	0.324	0.542
0.128	0.097	0.31	0.561
0.16	0.093	0.297	0.57
0.168	0.093	0.289	0.576
0.192	0.091	0.278	0.583
0.208	0.093	0.273	0.582
0.22	0.093	0.27	0.582
0.24	0.088	0.257	0.595
0.264	0.087	0.247	0.6
0.284	0.087	0.245	0.597
0.324	0.089	0.236	0.594
0.356	0.088	0.23	0.593
0.384	0.085	0.221	0.598
0.404	0.085	0.22	0.594
0.456	0.087	0.213	0.586
0.46	0.092	0.205	0.588
0.496	0.096	0.2	0.58
0.528	0.102	0.19	0.576
0.552	0.11	0.182	0.57
0.584	0.112	0.177	0.565
0.612	0.119	0.169	0.559
0.66	0.127	0.161	0.547
0.684	0.137	0.158	0.534
0.724	0.148	0.154	0.517
0.788	0.154	0.146	0.503
0.832	0.156	0.138	0.498
0.884	0.166	0.135	0.478
0.872	0.188	0.129	0.465
0.876	0.213	0.124	0.444
0.9	0.235	0.123	0.417
0.904	0.254	0.123	0.397
0.868	0.294	0.115	0.374
0.852	0.329	0.114	0.344
0.844	0.354	0.116	0.319
0.832	0.379	0.118	0.295
0.816	0.396	0.12	0.28
0.796	0.426	0.125	0.25
0.748	0.447	0.129	0.237
0.708	0.477	0.13	0.216
0.664	0.488	0.135	0.211
0.612	0.517	0.136	0.194
0.588	0.537	0.141	0.175
0.564	0.552	0.143	0.164
0.52	0.563	0.15	0.157
0.496	0.58	0.152	0.144
0.472	0.602	0.145	0.135
0.464	0.612	0.148	0.124
0.416	0.621	0.152	0.123
0.376	0.621	0.164	0.121
0.364	0.618	0.178	0.113
0.344	0.62	0.181	0.113
0.34	0.616	0.188	0.111
0.324	0.611	0.196	0.112
0.308	0.603	0.211	0.109
0.272	0.602	0.225	0.105
0.272	0.591	0.237	0.104
0.256	0.582	0.257	0.097
0.24	0.573	0.272	0.095
0.228	0.553	0.296	0.094
0.22	0.535	0.313	0.097
0.216	0.517	0.337	0.092
0.204	0.501	0.35	0.098
0.204	0.494	0.358	0.097
0.188	0.482	0.371	0.1
0.184	0.471	0.384	0.099
0.18	0.458	0.397	0.1
0.176	0.451	0.403	0.102
0.172	0.441	0.414	0.102
0.172	0.432	0.422	0.103
0.172	0.424	0.427	0.106
0.148	0.416	0.437	0.11
0.148	0.403	0.446	0.114
0.144	0.39	0.457	0.117
0.144	0.37	0.471	0.123
0.144	0.353	0.484	0.127
0.14	0.345	0.487	0.133
0.136	0.325	0.494	0.147
0.132	0.31	0.503	0.154
0.128	0.293	0.513	0.162
0.116	0.277	0.52	0.174
0.112	0.265	0.521	0.186
0.112	0.265	0.514	0.193
0.116	0.255	0.514	0.202
0.12	0.255	0.508	0.207
0.116	0.246	0.504	0.221
0.12	0.243	0.499	0.228
0.124	0.235	0.493	0.241
0.116	0.225	0.496	0.25
0.112	0.219	0.489	0.264
0.108	0.213	0.473	0.287
0.096	0.21	0.467	0.299
0.1	0.205	0.464	0.306
0.108	0.201	0.457	0.315
0.112	0.199	0.449	0.324
0.12	0.191	0.446	0.333
0.136	0.19	0.444	0.332
0.132	0.18	0.445	0.342
0.136	0.174	0.436	0.356
0.14	0.169	0.433	0.363
0.144	0.17	0.426	0.368
0.164	0.158	0.424	0.377
0.172	0.15	0.419	0.388
0.172	0.147	0.407	0.403
0.168	0.148	0.39	0.42
0.172	0.144	0.38	0.433
0.172	0.138	0.371	0.448
0.18	0.133	0.355	0.467
0.192	0.133	0.342	0.477
0.208	0.134	0.333	0.481
0.228	0.134	0.322	0.487
0.248	0.132	0.308	0.498
0.256	0.131	0.306	0.499
0.264	0.127	0.302	0.505
0.288	0.124	0.295	0.509
0.304	0.131	0.289	0.504
0.324	0.127	0.289	0.503
0.34	0.126	0.281	0.508
0.364	0.127	0.275	0.507
0.392	0.126	0.265	0.511
0.404	0.131	0.263	0.505
0.416	0.134	0.26	0.502
0.464	0.142	0.251	0.491
0.476	0.143	0.243	0.495
0.496	0.145	0.241	0.49
0.496	0.15	0.233	0.493
0.52	0.157	0.229	0.484
0.544	0.165	0.221	0.478
0.576	0.165	0.207	0.484
0.628	0.167	0.201	0.475
0.648	0.177	0.195	0.466
0.664	0.187	0.191	0.456
0.704	0.186	0.186	0.452
0.728	0.204	0.178	0.436
0.752	0.206	0.174	0.432
0.724	0.232	0.17	0.417
0.732	0.238	0.173	0.406
0.748	0.256	0.166	0.391
0.72	0.279	0.156	0.385
0.736	0.297	0.156	0.363
0.7	0.328	0.15	0.347
0.716	0.343	0.145	0.333
0.712	0.361	0.143	0.318
0.652	0.379	0.152	0.306
0.628	0.398	0.151	0.294
0.584	0.416	0.155	0.283
0.548	0.436	0.161	0.266
0.536	0.45	0.16	0.256
0.508	0.468	0.162	0.243
0.472	0.477	0.168	0.237
0.452	0.485	0.178	0.224
0.416	0.499	0.181	0.216
0.416	0.502	0.183	0.211
0.42	0.501	0.191	0.203
0.42	0.496	0.203	0.196
0.392	0.511	0.205	0.186
0.372	0.522	0.212	0.173
0.352	0.524	0.219	0.169
0.336	0.521	0.229	0.166
0.32	0.524	0.238	0.158
0.308	0.516	0.25	0.157
0.296	0.516	0.254	0.156
0.272	0.516	0.264	0.152
0.268	0.514	0.274	0.145
0.264	0.503	0.287	0.144
0.268	0.483	0.308	0.142
0.256	0.481	0.317	0.138
0.244	0.473	0.328	0.138
0.24	0.462	0.338	0.14
0.22	0.465	0.341	0.139
0.212	0.455	0.354	0.138
0.192	0.449	0.365	0.138
0.188	0.437	0.373	0.143
0.188	0.431	0.384	0.138
0.184	0.413	0.395	0.146
0.172	0.409	0.397	0.151
0.176	0.4	0.405	0.151
0.168	0.39	0.411	0.157
0.16	0.384	0.414	0.162
0.152	0.368	0.429	0.165
0.148	0.36	0.435	0.168
0.144	0.353	0.441	0.17
0.144	0.354	0.436	0.174
0.148	0.343	0.435	0.185
0.144	0.334	0.435	0.195
0.14	0.325	0.435	0.205
0.14	0.312	0.44	0.213
0.136	0.307	0.443	0.216
0.128	0.294	0.451	0.223
0.12	0.29	0.459	0.221
0.12	0.284	0.457	0.229
0.124	0.273	0.467	0.229
0.128	0.264	0.455	0.249
0.12	0.256	0.455	0.259
0.116	0.251	0.456	0.264
0.116	0.244	0.454	0.273
0.12	0.228	0.456	0.286
0.116	0.219	0.454	0.298
0.128	0.213	0.447	0.308
0.124	0.206	0.444	0.319
0.128	0.199	0.443	0.326
0.144	0.189	0.428	0.347
0.148	0.181	0.426	0.356
0.152	0.182	0.417	0.363
0.148	0.167	0.405	0.391
0.144	0.165	0.395	0.404
0.14	0.161	0.39	0.414
0.144	0.163	0.38	0.421
0.16	0.17	0.369	0.421
0.16	0.17	0.357	0.433
0.172	0.175	0.344	0.438
0.172	0.17	0.34	0.447
0.188	0.167	0.328	0.458
0.2	0.163	0.321	0.466
0.224	0.16	0.31	0.474
0.244	0.154	0.304	0.481
0.264	0.15	0.298	0.486
0.292	0.151	0.285	0.491
0.308	0.147	0.274	0.502
0.336	0.148	0.266	0.502
0.332	0.15	0.263	0.504
0.364	0.158	0.256	0.495
0.38	0.157	0.259	0.489
0.42	0.148	0.256	0.491
0.432	0.155	0.245	0.492
0.424	0.164	0.239	0.491
0.444	0.177	0.229	0.483
0.464	0.183	0.219	0.482
0.456	0.193	0.216	0.477
0.464	0.195	0.217	0.472
0.488	0.199	0.215	0.464
0.496	0.21	0.208	0.458
0.504	0.219	0.201	0.454
0.536	0.235	0.193	0.438
0.544	0.246	0.19	0.428
0.592	0.253	0.187	0.412
0.64	0.262	0.18	0.398
0.632	0.281	0.179	0.382
0.62	0.295	0.181	0.369
0.624	0.316	0.173	0.355
0.612	0.343	0.176	0.328
0.584	0.362	0.18	0.312
0.584	0.369	0.183	0.302
0.56	0.394	0.181	0.285
0.536	0.414	0.178	0.274
0.56	0.429	0.178	0.253
0.54	0.444	0.187	0.234
0.528	0.461	0.191	0.216
0.516	0.46	0.197	0.214
0.508	0.468	0.197	0.208
0.464	0.481	0.206	0.197
0.46	0.481	0.211	0.193
0.452	0.481	0.219	0.187
0.46	0.48	0.224	0.181
0.456	0.486	0.229	0.171
0.432	0.5	0.227	0.165
0.42	0.508	0.228	0.159
0.384	0.524	0.22	0.16
0.372	0.525	0.226	0.156
0.36	0.531	0.229	0.15
0.336	0.527	0.244	0.145
0.332	0.509	0.269	0.139
0.316	0.51	0.279	0.132
0.28	0.507	0.291	0.132
0.268	0.501	0.298	0.134
0.264	0.502	0.303	0.129
0.268	0.492	0.317	0.124
0.244	0.48	0.338	0.121
0.236	0.47	0.35	0.121
0.22	0.463	0.361	0.121
0.212	0.451	0.372	0.124
0.192	0.433	0.392	0.127
0.188	0.422	0.403	0.128
0.192	0.411	0.41	0.131
0.184	0.405	0.413	0.136
0.18	0.392	0.422	0.141
0.176	0.383	0.425	0.148
0.164	0.363	0.439	0.157
0.156	0.359	0.44	0.162
0.164	0.353	0.444	0.162
0.16	0.348	0.444	0.168
0.16	0.341	0.442	0.177
0.184	0.329	0.443	0.182
0.184	0.323	0.442	0.189
0.184	0.316	0.449	0.189
0.18	0.312	0.453	0.19
0.164	0.312	0.452	0.195
0.164	0.306	0.451	0.202
0.172	0.292	0.449	0.216
0.18	0.281	0.452	0.222
0.184	0.271	0.445	0.238
0.188	0.26	0.447	0.246
0.188	0.255	0.439	0.259
0.2	0.247	0.427	0.276
0.208	0.241	0.422	0.285
0.212	0.246	0.408	0.293
0.208	0.25	0.392	0.306
0.212	0.256	0.379	0.312
0.216	0.248	0.374	0.324
0.228	0.25	0.367	0.326
0.22	0.246	0.361	0.338
0.232	0.243	0.357	0.342
0.24	0.237	0.363	0.34
0.228	0.237	0.356	0.35
0.248	0.236	0.351	0.351
0.256	0.232	0.346	0.358
0.272	0.224	0.347	0.361
0.256	0.227	0.342	0.367
0.268	0.222	0.335	0.376
0.268	0.227	0.319	0.387
0.272	0.228	0.314	0.39
0.296	0.235	0.305	0.386
0.292	0.23	0.3	0.397
0.292	0.229	0.292	0.406
0.312	0.225	0.289	0.408
0.32	0.234	0.286	0.4
0.316	0.235	0.283	0.403
0.292	0.238	0.282	0.407
0.28	0.246	0.287	0.397
0.284	0.25	0.29	0.389
0.284	0.248	0.291	0.39
0.288	0.256	0.284	0.388
0.304	0.256	0.286	0.382
0.3	0.257	0.287	0.381
0.308	0.26	0.283	0.38
0.312	0.27	0.275	0.377
0.312	0.262	0.272	0.388
0.3	0.252	0.284	0.389
0.288	0.259	0.289	0.38
0.276	0.27	0.284	0.377
0.276	0.272	0.29	0.369
0.264	0.284	0.286	0.364
0.268	0.277	0.289	0.367
0.276	0.278	0.289	0.364
0.284	0.276	0.289	0.364
0.276	0.277	0.282	0.372
0.284	0.277	0.281	0.371
0.3	0.275	0.286	0.364
0.308	0.291	0.28	0.352
0.296	0.292	0.282	0.352
0.316	0.297	0.276	0.348
0.328	0.287	0.285	0.346
0.32	0.289	0.282	0.349
0.34	0.286	0.278	0.351
0.336	0.284	0.276	0.356
0.336	0.284	0.278	0.354
0.336	0.289	0.276	0.351
0.352	0.291	0.28	0.341
0.344	0.288	0.281	0.345
0.348	0.291	0.276	0.346
0.352	0.296	0.276	0.34
0.372	0.3	0.272	0.335
0.38	0.307	0.264	0.334
0.392	0.302	0.264	0.336
0.404	0.305	0.267	0.327
0.412	0.306	0.272	0.319
0.396	0.313	0.276	0.312
0.4	0.316	0.285	0.299
0.384	0.32	0.292	0.292
0.372	0.326	0.29	0.291
0.368	0.326	0.294	0.288
0.368	0.333	0.298	0.277
0.372	0.335	0.295	0.277
0.376	0.329	0.307	0.27
0.368	0.324	0.308	0.276
0.372	0.328	0.297	0.282
0.4	0.333	0.288	0.279
0.388	0.331	0.286	0.286
0.392	0.33	0.286	0.286
0.396	0.323	0.285	0.293
0.388	0.329	0.283	0.291
0.372	0.333	0.288	0.286
0.368	0.342	0.285	0.281
0.356	0.347	0.286	0.278
0.332	0.361	0.282	0.274
0.308	0.357	0.296	0.27
0.284	0.359	0.304	0.266
0.272	0.358	0.306	0.268
0.272	0.355	0.312	0.265
0.268	0.354	0.313	0.266
0.264	0.359	0.315	0.26
0.272	0.342	0.32	0.27
0.264	0.341	0.323	0.27
0.26	0.335	0.329	0.271
0.28	0.33	0.332	0.268
0.276	0.319	0.334	0.278
0.28	0.317	0.332	0.281
0.272	0.314	0.33	0.288
0.272	0.312	0.323	0.297
0.26	0.309	0.322	0.304
0.248	0.306	0.32	0.312
0.244	0.309	0.312	0.318
0.244	0.312	0.311	0.316
0.252	0.315	0.308	0.314
0.256	0.314	0.308	0.314
0.236	0.318	0.304	0.319
0.244	0.307	0.31	0.322
0.268	0.304	0.309	0.32
0.284	0.298	0.314	0.317
0.28	0.295	0.323	0.312
0.288	0.29	0.331	0.307
0.28	0.287	0.326	0.317
0.272	0.277	0.327	0.328
0.272	0.276	0.322	0.334
0.264	0.287	0.322	0.325
0.276	0.282	0.312	0.337
0.264	0.275	0.32	0.339
0.264	0.271	0.32	0.343
0.272	0.271	0.313	0.348
0.28	0.275	0.315	0.34
0.268	0.285	0.313	0.335
0.264	0.283	0.316	0.335
0.26	0.289	0.311	0.335
0.26	0.288	0.307	0.34
0.272	0.295	0.294	0.343
0.26	0.291	0.299	0.345
0.26	0.285	0.307	0.343
0.252	0.278	0.305	0.354
0.248	0.278	0.303	0.357
0.26	0.27	0.304	0.361
0.26	0.261	0.308	0.366
0.264	0.262	0.31	0.362
0.28	0.26	0.305	0.365
0.28	0.26	0.302	0.368
0.268	0.26	0.303	0.37
0.288	0.259	0.301	0.368
0.308	0.264	0.295	0.364
0.304	0.26	0.298	0.366
0.304	0.274	0.289	0.361
0.288	0.275	0.289	0.364
0.288	0.267	0.289	0.372
0.288	0.275	0.281	0.372
0.292	0.279	0.276	0.372
0.3	0.276	0.275	0.374
0.292	0.287	0.271	0.369
0.284	0.286	0.281	0.362
0.3	0.287	0.277	0.361
0.304	0.279	0.281	0.364
0.3	0.278	0.282	0.365
0.316	0.275	0.271	0.375
0.32	0.279	0.265	0.376
0.3	0.279	0.264	0.382
0.312	0.283	0.262	0.377
0.32	0.288	0.266	0.366
0.324	0.285	0.269	0.365
0.344	0.287	0.268	0.359
0.348	0.282	0.269	0.362
0.34	0.273	0.272	0.37
0.336	0.279	0.265	0.372
0.376	0.281	0.263	0.362
0.384	0.28	0.268	0.356
0.376	0.28	0.271	0.355
0.38	0.285	0.273	0.347
0.368	0.291	0.272	0.345
0.36	0.289	0.272	0.349
0.36	0.289	0.271	0.35
0.368	0.292	0.268	0.348
0.372	0.3	0.265	0.342
0.356	0.299	0.263	0.349
0.38	0.291	0.263	0.351
0.4	0.286	0.269	0.345
0.388	0.29	0.269	0.344
0.392	0.301	0.267	0.334
0.408	0.311	0.261	0.326
0.396	0.304	0.273	0.324
0.392	0.31	0.269	0.323
0.38	0.325	0.263	0.317
0.38	0.326	0.263	0.316
0.384	0.332	0.263	0.309
0.4	0.329	0.261	0.31
0.388	0.336	0.262	0.305
0.372	0.341	0.265	0.301
0.384	0.346	0.266	0.292
0.372	0.348	0.27	0.289
0.352	0.351	0.269	0.292
0.348	0.354	0.27	0.289
0.356	0.359	0.27	0.282
0.356	0.357	0.272	0.282
0.336	0.353	0.279	0.284
0.324	0.343	0.29	0.286
0.332	0.345	0.284	0.288
0.296	0.359	0.283	0.284
0.308	0.349	0.289	0.285
0.312	0.336	0.296	0.29
0.316	0.335	0.296	0.29
0.312	0.335	0.295	0.292
0.308	0.329	0.304	0.29
0.292	0.332	0.304	0.291
0.292	0.329	0.307	0.291
0.288	0.336	0.309	0.283
0.288	0.336	0.307	0.285
0.296	0.339	0.311	0.276
0.308	0.333	0.312	0.278
0.316	0.326	0.316	0.279
0.316	0.327	0.316	0.278
0.324	0.318	0.323	0.278
0.3	0.327	0.32	0.278
0.304	0.322	0.324	0.278
0.308	0.329	0.325	0.269
0.296	0.325	0.335	0.266
0.296	0.319	0.339	0.268
0.3	0.332	0.326	0.267
0.296	0.325	0.328	0.273
0.292	0.319	0.343	0.265
0.284	0.316	0.352	0.261
0.276	0.315	0.348	0.268
0.264	0.321	0.346	0.267
0.276	0.321	0.342	0.268
0.288	0.317	0.336	0.275
0.28	0.312	0.335	0.283
0.276	0.306	0.337	0.288
0.268	0.301	0.339	0.293
0.264	0.289	0.337	0.308
0.272	0.285	0.331	0.316
0.276	0.285	0.335	0.311
0.292	0.286	0.332	0.309
0.288	0.293	0.322	0.313
0.284	0.3	0.314	0.315
0.292	0.295	0.31	0.322
0.3	0.291	0.308	0.326
0.308	0.288	0.308	0.327
0.3	0.285	0.308	0.332
0.316	0.288	0.304	0.329
0.316	0.291	0.308	0.322
0.332	0.292	0.302	0.323
0.328	0.294	0.303	0.321
0.336	0.294	0.306	0.316
0.34	0.286	0.31	0.319
0.336	0.29	0.309	0.317
0.344	0.29	0.303	0.321
0.336	0.3	0.296	0.32
0.332	0.303	0.293	0.321
0.336	0.304	0.291	0.321
0.348	0.309	0.284	0.32
0.356	0.302	0.282	0.327
0.348	0.313	0.28	0.32
0.348	0.31	0.287	0.316
0.38	0.309	0.288	0.308
0.396	0.32	0.279	0.302
0.38	0.319	0.277	0.309
0.396	0.328	0.274	0.299
0.404	0.337	0.27	0.292
0.4	0.35	0.259	0.291
0.4	0.356	0.265	0.279
0.376	0.365	0.266	0.275
0.4	0.363	0.264	0.273
0.404	0.37	0.262	0.267
0.396	0.369	0.268	0.264
0.38	0.387	0.258	0.26
0.4	0.388	0.257	0.255
0.384	0.389	0.258	0.257
0.4	0.386	0.257	0.257
0.4	0.4	0.252	0.248
0.388	0.405	0.251	0.247
0.404	0.413	0.248	0.238
0.396	0.414	0.257	0.23
0.38	0.41	0.268	0.227
0.38	0.408	0.268	0.229
0.364	0.411	0.267	0.231
0.36	0.42	0.264	0.226
0.344	0.421	0.266	0.227
0.356	0.422	0.269	0.22
0.348	0.415	0.279	0.219
0.324	0.412	0.286	0.221
0.328	0.417	0.286	0.215
0.34	0.413	0.292	0.21
0.34	0.416	0.297	0.202
0.332	0.414	0.31	0.193
0.324	0.413	0.31	0.196
0.304	0.41	0.326	0.188
0.304	0.411	0.325	0.188
0.288	0.413	0.329	0.186
0.284	0.414	0.327	0.188
0.26	0.414	0.329	0.192
0.248	0.413	0.331	0.194
0.248	0.402	0.341	0.195
0.224	0.405	0.337	0.202
0.216	0.392	0.343	0.211
0.224	0.379	0.35	0.215
0.212	0.38	0.352	0.215
0.208	0.389	0.336	0.223
0.208	0.378	0.342	0.228
0.208	0.372	0.35	0.226
0.208	0.365	0.358	0.225
0.212	0.353	0.367	0.227
0.212	0.349	0.361	0.237
0.22	0.343	0.361	0.241
0.224	0.341	0.367	0.236
0.224	0.331	0.377	0.236
0.212	0.324	0.378	0.245
0.216	0.31	0.382	0.254
0.22	0.299	0.392	0.254
0.216	0.285	0.4	0.261
0.216	0.287	0.399	0.26
0.232	0.29	0.391	0.261
0.24	0.291	0.383	0.266
0.248	0.296	0.372	0.27
0.252	0.292	0.37	0.275
0.26	0.28	0.377	0.278
0.26	0.279	0.373	0.283
0.256	0.273	0.377	0.286
0.256	0.272	0.38	0.284
0.244	0.271	0.376	0.292
0.224	0.276	0.373	0.295
0.24	0.278	0.368	0.294
0.24	0.283	0.365	0.292
0.24	0.281	0.36	0.299
0.22	0.268	0.361	0.316
0.232	0.264	0.358	0.32
0.228	0.258	0.358	0.327
0.22	0.258	0.353	0.334
0.212	0.252	0.351	0.344
0.24	0.248	0.348	0.344
0.26	0.238	0.348	0.349
0.268	0.246	0.341	0.346
0.276	0.241	0.335	0.355
0.264	0.246	0.331	0.357
0.272	0.242	0.328	0.362
0.272	0.246	0.327	0.359
0.272	0.239	0.321	0.372
0.288	0.232	0.318	0.378
0.288	0.223	0.32	0.385
0.3	0.226	0.312	0.387
0.32	0.22	0.309	0.391
0.324	0.223	0.308	0.388
0.328	0.22	0.304	0.394
0.332	0.223	0.298	0.396
0.344	0.22	0.293	0.401
0.36	0.218	0.28	0.412
0.38	0.213	0.279	0.413
0.372	0.218	0.276	0.413
0.376	0.221	0.272	0.413
0.372	0.215	0.27	0.422
0.38	0.217	0.269	0.419
0.404	0.219	0.265	0.415
0.404	0.228	0.259	0.412
0.44	0.234	0.253	0.403
0.424	0.251	0.248	0.395
0.452	0.259	0.251	0.377
0.448	0.265	0.249	0.374
0.452	0.276	0.241	0.37
0.456	0.3	0.232	0.354
0.444	0.307	0.234	0.348
0.452	0.304	0.236	0.347
0.436	0.313	0.237	0.341
0.42	0.317	0.237	0.341
0.436	0.327	0.235	0.329
0.444	0.342	0.222	0.325
0.412	0.349	0.227	0.321
0.388	0.377	0.22	0.306
0.4	0.386	0.219	0.295
0.388	0.398	0.215	0.29
0.38	0.407	0.215	0.283
0.372	0.414	0.217	0.276
0.372	0.412	0.222	0.273
0.364	0.415	0.229	0.265
0.36	0.415	0.237	0.258
0.344	0.426	0.239	0.249
0.34	0.428	0.245	0.242
0.336	0.433	0.25	0.233
0.316	0.435	0.259	0.227
0.312	0.437	0.261	0.224
0.292	0.434	0.271	0.222
0.264	0.442	0.275	0.217
0.264	0.437	0.279	0.218
0.244	0.449	0.277	0.213
0.248	0.439	0.289	0.21
0.24	0.437	0.293	0.21
0.248	0.434	0.296	0.208
0.252	0.43	0.3	0.207
0.26	0.424	0.309	0.202
0.248	0.423	0.314	0.201
0.256	0.411	0.315	0.21
0.248	0.406	0.321	0.211
0.244	0.4	0.335	0.204
0.244	0.39	0.346	0.203
0.248	0.378	0.35	0.21
0.256	0.364	0.357	0.215
0.248	0.365	0.357	0.216
0.252	0.359	0.355	0.223
0.236	0.351	0.363	0.227
0.232	0.349	0.365	0.228
0.224	0.341	0.369	0.234
0.204	0.332	0.376	0.241
0.204	0.325	0.374	0.25
0.204	0.334	0.361	0.254
0.2	0.322	0.361	0.267
0.196	0.309	0.361	0.281
0.196	0.305	0.357	0.289
0.196	0.301	0.361	0.289
0.196	0.294	0.367	0.29
0.204	0.291	0.364	0.294
0.212	0.286	0.364	0.297
0.22	0.282	0.364	0.299
0.224	0.266	0.37	0.308
0.228	0.259	0.372	0.312
0.224	0.264	0.355	0.325
0.228	0.258	0.35	0.335
0.236	0.253	0.341	0.347
0.232	0.248	0.338	0.356
0.252	0.244	0.338	0.355
0.26	0.242	0.332	0.361
0.248	0.236	0.329	0.373
0.252	0.231	0.326	0.38
0.256	0.224	0.319	0.393
0.26	0.227	0.31	0.398
0.284	0.221	0.316	0.392
0.292	0.223	0.315	0.389
0.3	0.217	0.314	0.394
0.32	0.218	0.306	0.396
0.32	0.214	0.312	0.394
0.336	0.212	0.305	0.399
0.348	0.204	0.3	0.409
0.368	0.206	0.296	0.406
0.384	0.203	0.287	0.414
0.4	0.203	0.281	0.416
0.416	0.213	0.268	0.415
0.412	0.218	0.259	0.42
0.424	0.222	0.255	0.417
0.42	0.227	0.248	0.42
0.408	0.234	0.246	0.418
0.452	0.241	0.232	0.414
0.46	0.248	0.23	0.407
0.448	0.262	0.232	0.394
0.456	0.269	0.232	0.385
0.444	0.28	0.233	0.376
0.444	0.289	0.231	0.369
0.44	0.3	0.227	0.363
0.448	0.303	0.226	0.359
0.428	0.308	0.229	0.356
0.436	0.318	0.225	0.348
0.436	0.325	0.232	0.334
0.444	0.331	0.235	0.323
0.46	0.335	0.24	0.31
0.464	0.344	0.249	0.291
0.448	0.364	0.243	0.281
0.428	0.375	0.239	0.279
0.46	0.375	0.234	0.276
0.456	0.379	0.234	0.273
0.448	0.387	0.232	0.269
0.452	0.388	0.24	0.259
0.468	0.389	0.248	0.246
0.476	0.387	0.252	0.242
0.488	0.394	0.246	0.238
0.464	0.396	0.25	0.238
0.456	0.409	0.247	0.23
0.452	0.405	0.255	0.227
0.432	0.411	0.259	0.222
0.44	0.421	0.255	0.214
0.416	0.422	0.266	0.208
0.392	0.432	0.267	0.203
0.36	0.447	0.267	0.196
0.36	0.437	0.275	0.198
0.376	0.434	0.281	0.191
0.352	0.443	0.287	0.182
0.32	0.443	0.303	0.174
0.312	0.423	0.32	0.179
0.292	0.419	0.327	0.181
0.288	0.406	0.341	0.181
0.284	0.4	0.344	0.185
0.268	0.399	0.349	0.185
0.268	0.399	0.352	0.182
0.232	0.407	0.357	0.178
0.236	0.399	0.365	0.177
0.228	0.392	0.375	0.176
0.22	0.387	0.38	0.178
0.216	0.38	0.383	0.183
0.216	0.375	0.388	0.183
0.204	0.366	0.394	0.189
0.196	0.363	0.394	0.194
0.188	0.361	0.397	0.195
0.176	0.358	0.399	0.199
0.168	0.346	0.407	0.205
0.184	0.329	0.416	0.209
0.196	0.318	0.42	0.213
0.196	0.309	0.423	0.219
0.188	0.298	0.423	0.232
0.176	0.293	0.425	0.238
0.196	0.293	0.426	0.232
0.18	0.29	0.424	0.241
0.176	0.296	0.416	0.244
0.18	0.291	0.417	0.247
0.172	0.283	0.418	0.256
0.176	0.277	0.413	0.266
0.168	0.273	0.409	0.276
0.168	0.263	0.415	0.28
0.16	0.26	0.418	0.282
0.16	0.25	0.418	0.292
0.16	0.247	0.414	0.299
0.16	0.242	0.41	0.308
0.164	0.237	0.398	0.324
0.172	0.232	0.389	0.336
0.164	0.233	0.38	0.346
0.168	0.225	0.373	0.36
0.18	0.222	0.365	0.368
0.184	0.227	0.362	0.365
0.2	0.215	0.368	0.367
0.188	0.219	0.365	0.369
0.2	0.226	0.357	0.367
0.208	0.222	0.351	0.375
0.224	0.217	0.345	0.382
0.224	0.217	0.339	0.388
0.228	0.211	0.333	0.399
0.224	0.208	0.326	0.41
0.236	0.205	0.321	0.415
0.252	0.205	0.315	0.417
0.248	0.21	0.31	0.418
0.268	0.204	0.31	0.419
0.268	0.194	0.31	0.429
0.284	0.188	0.306	0.435
0.316	0.186	0.299	0.436
0.324	0.189	0.294	0.436
0.348	0.194	0.285	0.434
0.348	0.191	0.287	0.435
0.376	0.191	0.288	0.427
0.4	0.192	0.286	0.422
0.42	0.196	0.283	0.416
0.452	0.195	0.277	0.415
0.46	0.204	0.27	0.411
0.48	0.202	0.265	0.413
0.472	0.204	0.266	0.412
0.492	0.213	0.252	0.412
0.52	0.226	0.243	0.401
0.488	0.24	0.248	0.39
0.508	0.244	0.241	0.388
0.528	0.248	0.234	0.386
0.508	0.26	0.233	0.38
0.532	0.264	0.226	0.377
0.52	0.28	0.22	0.37
0.504	0.278	0.228	0.368
0.488	0.285	0.227	0.366
0.492	0.291	0.228	0.358
0.516	0.3	0.223	0.348
0.504	0.313	0.22	0.341
0.516	0.32	0.214	0.337
0.52	0.327	0.213	0.33
0.508	0.342	0.207	0.324
0.48	0.357	0.205	0.318
0.508	0.357	0.211	0.305
0.496	0.371	0.203	0.302
0.504	0.383	0.204	0.287
0.5	0.398	0.206	0.271
0.492	0.408	0.208	0.261
0.484	0.42	0.204	0.255
0.484	0.434	0.204	0.241
0.484	0.438	0.208	0.233
0.46	0.45	0.214	0.221
0.448	0.452	0.216	0.22
0.448	0.453	0.214	0.221
0.428	0.459	0.211	0.223
0.412	0.47	0.215	0.212
0.396	0.473	0.223	0.205
0.384	0.472	0.233	0.199
0.38	0.473	0.237	0.195
0.388	0.459	0.246	0.198
0.34	0.46	0.26	0.195
0.324	0.455	0.268	0.196
0.304	0.45	0.28	0.194
0.276	0.446	0.294	0.191
0.272	0.433	0.304	0.195
0.264	0.423	0.31	0.201
0.26	0.418	0.314	0.203
0.256	0.404	0.327	0.205
0.232	0.407	0.335	0.2
0.24	0.389	0.353	0.198
0.24	0.374	0.36	0.206
0.24	0.366	0.362	0.212
0.224	0.365	0.362	0.217
0.228	0.366	0.36	0.217
0.216	0.362	0.36	0.224
0.208	0.356	0.362	0.23
0.208	0.349	0.369	0.23
0.2	0.354	0.365	0.231
0.2	0.354	0.363	0.233
0.204	0.34	0.372	0.237
0.212	0.334	0.373	0.24
0.208	0.321	0.38	0.247
0.196	0.317	0.387	0.247
0.188	0.313	0.39	0.25
0.172	0.299	0.399	0.259
0.16	0.289	0.403	0.268
0.16	0.277	0.416	0.267
0.156	0.272	0.417	0.272
0.152	0.266	0.417	0.279
0.16	0.257	0.417	0.286
0.156	0.247	0.417	0.297
0.152	0.25	0.417	0.295
0.16	0.246	0.414	0.3
0.188	0.237	0.406	0.31
0.188	0.242	0.398	0.313
0.188	0.238	0.389	0.326
0.192	0.236	0.38	0.336
0.2	0.226	0.379	0.345
0.22	0.219	0.377	0.349
0.208	0.217	0.376	0.355
0.208	0.215	0.372	0.361
0.22	0.211	0.362	0.372
0.236	0.212	0.351	0.378
0.244	0.211	0.343	0.385
0.236	0.215	0.339	0.387
0.24	0.22	0.33	0.39
0.24	0.224	0.325	0.391
0.24	0.228	0.323	0.389
0.248	0.228	0.322	0.388
0.256	0.232	0.303	0.401
0.248	0.23	0.308	0.4
0.244	0.232	0.31	0.397
0.24	0.231	0.305	0.404
0.252	0.236	0.299	0.402
0.268	0.23	0.297	0.406
0.272	0.221	0.295	0.416
0.288	0.221	0.297	0.41
0.3	0.217	0.292	0.416
0.288	0.227	0.284	0.417
0.28	0.231	0.287	0.412
0.292	0.225	0.289	0.413
0.288	0.233	0.282	0.413
0.284	0.229	0.277	0.423
0.304	0.23	0.274	0.42
0.324	0.227	0.269	0.423
0.332	0.227	0.266	0.424
0.356	0.238	0.257	0.416
0.368	0.235	0.257	0.416
0.36	0.244	0.253	0.413
0.372	0.239	0.255	0.413
0.388	0.236	0.249	0.418
0.404	0.236	0.254	0.409
0.4	0.256	0.245	0.399
0.404	0.264	0.248	0.387
0.412	0.266	0.251	0.38
0.424	0.265	0.254	0.375
0.412	0.271	0.255	0.371
0.44	0.277	0.242	0.371
0.452	0.284	0.236	0.367
0.448	0.298	0.231	0.359
0.46	0.306	0.226	0.353
0.436	0.314	0.228	0.349
0.456	0.323	0.219	0.344
0.456	0.319	0.221	0.346
0.46	0.322	0.224	0.339
0.472	0.333	0.223	0.326
0.472	0.343	0.22	0.319
0.488	0.353	0.216	0.309
0.456	0.36	0.221	0.305
0.432	0.368	0.229	0.295
0.432	0.386	0.222	0.284
0.428	0.394	0.232	0.267
0.428	0.402	0.233	0.258
0.408	0.421	0.23	0.247
0.4	0.417	0.244	0.239
0.38	0.425	0.25	0.23
0.352	0.43	0.256	0.226
0.348	0.432	0.256	0.225
0.332	0.424	0.259	0.234
0.328	0.41	0.266	0.242
0.32	0.412	0.273	0.235
0.308	0.42	0.274	0.229
0.304	0.41	0.285	0.229
0.332	0.4	0.292	0.225
0.344	0.389	0.301	0.224
0.36	0.388	0.301	0.221
0.352	0.39	0.304	0.218
0.34	0.388	0.315	0.212
0.312	0.389	0.323	0.21
0.268	0.39	0.331	0.212
0.272	0.382	0.34	0.21
0.26	0.388	0.341	0.206
0.248	0.38	0.349	0.209
0.24	0.379	0.342	0.219
0.236	0.369	0.356	0.216
0.236	0.36	0.361	0.22
0.216	0.357	0.365	0.224
0.228	0.343	0.383	0.217
0.232	0.344	0.382	0.216
0.232	0.34	0.382	0.22
0.24	0.328	0.392	0.22
0.236	0.325	0.394	0.222
0.216	0.329	0.391	0.226
0.212	0.329	0.4	0.218
0.212	0.33	0.405	0.212
0.2	0.327	0.413	0.21
0.184	0.328	0.409	0.217
0.196	0.32	0.409	0.222
0.2	0.313	0.411	0.226
0.188	0.307	0.414	0.232
0.196	0.298	0.413	0.24
0.2	0.292	0.403	0.255
0.196	0.285	0.401	0.265
0.192	0.282	0.398	0.272
0.184	0.282	0.402	0.27
0.18	0.274	0.411	0.27
0.18	0.268	0.405	0.282
0.168	0.274	0.394	0.29
0.172	0.286	0.382	0.289
0.164	0.28	0.381	0.298
0.172	0.278	0.377	0.302
0.16	0.269	0.378	0.313
0.168	0.266	0.375	0.317
0.168	0.263	0.374	0.321
0.168	0.26	0.384	0.314
0.184	0.256	0.384	0.314
0.196	0.248	0.381	0.322
0.196	0.237	0.387	0.327
0.184	0.233	0.384	0.337
0.192	0.226	0.383	0.343
0.2	0.21	0.383	0.357
0.22	0.201	0.375	0.369
0.22	0.199	0.376	0.37
0.216	0.204	0.366	0.376
0.22	0.205	0.357	0.383
0.224	0.202	0.348	0.394
0.22	0.202	0.339	0.404
0.24	0.201	0.344	0.395
0.232	0.2	0.343	0.399
0.224	0.195	0.342	0.407
0.224	0.191	0.345	0.408
0.24	0.19	0.343	0.407
0.248	0.192	0.333	0.413
0.264	0.192	0.325	0.417
0.288	0.191	0.319	0.418
0.3	0.19	0.321	0.414
0.312	0.194	0.318	0.41
0.316	0.192	0.316	0.413
0.336	0.198	0.296	0.422
0.368	0.204	0.288	0.416
0.38	0.212	0.286	0.407
0.404	0.203	0.282	0.414
0.428	0.206	0.278	0.409
0.44	0.214	0.276	0.4
0.46	0.228	0.266	0.391
0.468	0.237	0.253	0.393
0.468	0.235	0.255	0.393
0.484	0.24	0.251	0.388
0.528	0.241	0.244	0.383
0.576	0.243	0.236	0.377
0.564	0.255	0.227	0.377
0.568	0.268	0.222	0.368
0.556	0.287	0.213	0.361
0.548	0.302	0.201	0.36
0.536	0.324	0.196	0.346
0.568	0.329	0.19	0.339
0.564	0.337	0.192	0.33
0.564	0.35	0.188	0.321
0.596	0.357	0.178	0.316
0.592	0.367	0.187	0.298
0.596	0.369	0.192	0.29
0.6	0.382	0.192	0.276
0.568	0.406	0.191	0.261
0.54	0.426	0.186	0.253
0.536	0.426	0.196	0.244
0.52	0.434	0.195	0.241
0.504	0.447	0.193	0.234
0.472	0.459	0.192	0.231
0.44	0.466	0.198	0.226
0.432	0.463	0.205	0.224
0.412	0.464	0.22	0.213
0.416	0.467	0.229	0.2
0.388	0.473	0.238	0.192
0.372	0.471	0.251	0.185
0.364	0.47	0.261	0.178
0.344	0.479	0.26	0.175
0.332	0.482	0.263	0.172
0.308	0.471	0.275	0.177
0.3	0.48	0.274	0.171
0.296	0.475	0.276	0.175
0.288	0.469	0.285	0.174
0.28	0.47	0.293	0.167
0.268	0.466	0.302	0.165
0.256	0.454	0.317	0.165
0.26	0.443	0.325	0.167
0.272	0.436	0.332	0.164
0.264	0.439	0.336	0.159
0.252	0.439	0.337	0.161
0.224	0.429	0.351	0.164
0.224	0.422	0.353	0.169
0.208	0.432	0.354	0.162
0.196	0.426	0.362	0.163
0.188	0.426	0.363	0.164
0.184	0.412	0.376	0.166
0.184	0.41	0.383	0.161
0.176	0.394	0.397	0.165
0.176	0.387	0.398	0.171
0.168	0.386	0.406	0.166
0.144	0.371	0.416	0.177
0.136	0.364	0.418	0.184
0.128	0.355	0.422	0.191
0.132	0.339	0.434	0.194
0.12	0.328	0.438	0.204
0.128	0.318	0.434	0.216
0.132	0.304	0.437	0.226
0.132	0.286	0.442	0.239
0.14	0.274	0.438	0.253
0.144	0.26	0.446	0.258
0.148	0.253	0.45	0.26
0.152	0.243	0.445	0.274
0.16	0.235	0.442	0.283
0.16	0.225	0.435	0.3
0.148	0.221	0.435	0.307
0.14	0.217	0.434	0.314
0.144	0.211	0.434	0.319
0.14	0.208	0.424	0.333
0.14	0.202	0.418	0.345
0.132	0.192	0.415	0.36
0.14	0.193	0.406	0.366
0.14	0.181	0.405	0.379
0.14	0.171	0.398	0.396
0.148	0.159	0.389	0.415
0.14	0.151	0.383	0.431
0.144	0.148	0.369	0.447
0.16	0.144	0.356	0.46
0.168	0.143	0.351	0.464
0.176	0.138	0.346	0.472
0.184	0.134	0.332	0.488
0.204	0.138	0.327	0.484
0.208	0.131	0.326	0.491
0.212	0.13	0.316	0.501
0.236	0.135	0.3	0.506
0.256	0.13	0.299	0.507
0.264	0.132	0.287	0.515
0.292	0.135	0.277	0.515
0.304	0.133	0.264	0.527
0.336	0.132	0.254	0.53
0.356	0.126	0.258	0.527
0.356	0.133	0.254	0.524
0.408	0.136	0.246	0.516
0.428	0.139	0.24	0.514
0.424	0.139	0.24	0.515
0.44	0.134	0.238	0.518
0.448	0.14	0.231	0.517
0.492	0.141	0.22	0.516
0.512	0.149	0.215	0.508
0.516	0.159	0.214	0.498
0.516	0.158	0.214	0.499
0.556	0.16	0.202	0.499
0.596	0.169	0.2	0.482
0.604	0.178	0.202	0.469
0.596	0.197	0.202	0.452
0.612	0.207	0.2	0.44
0.644	0.208	0.198	0.433
0.672	0.212	0.194	0.426
0.676	0.222	0.186	0.423
0.7	0.227	0.19	0.408
0.704	0.242	0.187	0.395
0.72	0.251	0.183	0.386
0.732	0.27	0.177	0.37
0.74	0.293	0.173	0.349
0.696	0.316	0.17	0.34
0.664	0.343	0.164	0.327
0.672	0.363	0.164	0.305
0.652	0.386	0.161	0.29
0.644	0.401	0.161	0.277
0.592	0.432	0.17	0.25
0.564	0.456	0.165	0.238
0.544	0.466	0.168	0.23
0.512	0.484	0.17	0.218
0.48	0.502	0.164	0.214
0.444	0.515	0.163	0.211
0.448	0.522	0.16	0.206
0.436	0.536	0.16	0.195
0.408	0.537	0.167	0.194
0.404	0.546	0.164	0.189
0.396	0.548	0.165	0.188
0.376	0.56	0.168	0.178
0.372	0.563	0.171	0.173
0.328	0.575	0.183	0.16
0.304	0.571	0.197	0.156
0.296	0.558	0.207	0.161
0.296	0.544	0.221	0.161
0.296	0.532	0.23	0.164
0.264	0.53	0.241	0.163
0.248	0.524	0.253	0.161
0.228	0.525	0.259	0.159
0.228	0.509	0.273	0.161
0.236	0.496	0.282	0.163
0.224	0.492	0.29	0.162
0.208	0.478	0.303	0.167
0.188	0.465	0.319	0.169
0.212	0.458	0.325	0.164
0.192	0.447	0.342	0.163
0.188	0.436	0.353	0.164
0.192	0.42	0.366	0.166
0.192	0.413	0.366	0.173
0.18	0.404	0.372	0.179
0.188	0.401	0.378	0.174
0.172	0.394	0.381	0.182
0.176	0.373	0.391	0.192
0.172	0.355	0.396	0.206
0.16	0.346	0.398	0.216
0.16	0.342	0.4	0.218
0.152	0.328	0.405	0.229
0.156	0.321	0.407	0.233
0.144	0.307	0.415	0.242
0.16	0.301	0.419	0.24
0.156	0.286	0.428	0.247
0.16	0.278	0.428	0.254
0.176	0.27	0.432	0.254
0.168	0.272	0.43	0.256
0.164	0.26	0.437	0.262
0.172	0.255	0.438	0.264
0.184	0.252	0.439	0.263
0.176	0.251	0.434	0.271
0.172	0.249	0.427	0.281
0.164	0.238	0.426	0.295
0.168	0.236	0.425	0.297
0.168	0.227	0.422	0.309
0.168	0.225	0.417	0.316
0.164	0.225	0.406	0.328
0.172	0.222	0.401	0.334
0.18	0.221	0.39	0.344
0.176	0.22	0.385	0.351
0.172	0.214	0.382	0.361
0.172	0.212	0.372	0.373
0.184	0.207	0.368	0.379
0.188	0.207	0.358	0.388
0.192	0.2	0.359	0.393
0.2	0.206	0.348	0.396
0.2	0.203	0.346	0.401
0.204	0.194	0.348	0.407
0.212	0.189	0.341	0.417
0.224	0.182	0.343	0.419
0.24	0.181	0.337	0.422
0.256	0.182	0.328	0.426
0.264	0.176	0.33	0.428
0.268	0.171	0.326	0.436
0.276	0.167	0.322	0.442
0.284	0.156	0.327	0.446
0.288	0.155	0.323	0.45
0.288	0.154	0.321	0.453
0.3	0.151	0.312	0.462
0.312	0.15	0.307	0.465
0.332	0.148	0.3	0.469
0.348	0.148	0.293	0.472
0.38	0.155	0.28	0.47
0.396	0.162	0.264	0.475
0.392	0.159	0.261	0.482
0.404	0.16	0.255	0.484
0.42	0.162	0.246	0.487
0.432	0.169	0.237	0.486
0.464	0.184	0.226	0.474
0.512	0.196	0.216	0.46
0.552	0.208	0.208	0.446
0.572	0.233	0.2	0.424
0.58	0.244	0.193	0.418
0.572	0.249	0.193	0.415
0.572	0.254	0.191	0.412
0.588	0.27	0.186	0.397
0.612	0.278	0.184	0.385
0.612	0.283	0.183	0.381
0.624	0.297	0.176	0.371
0.604	0.305	0.181	0.363
0.628	0.309	0.182	0.352
0.624	0.323	0.183	0.338
0.636	0.34	0.177	0.324
0.636	0.362	0.172	0.307
0.632	0.386	0.16	0.296
0.62	0.406	0.15	0.289
0.576	0.422	0.158	0.276
0.552	0.444	0.163	0.255
0.504	0.463	0.161	0.25
0.48	0.466	0.17	0.244
0.46	0.468	0.176	0.241
0.452	0.48	0.178	0.229
0.452	0.49	0.181	0.216
0.432	0.498	0.19	0.204
0.444	0.499	0.192	0.198
0.436	0.502	0.198	0.191
0.416	0.5	0.208	0.188
0.396	0.501	0.212	0.188
0.392	0.509	0.212	0.181
0.36	0.518	0.217	0.175
0.336	0.521	0.221	0.174
0.312	0.514	0.243	0.165
0.3	0.506	0.258	0.161
0.28	0.512	0.262	0.156
0.272	0.502	0.276	0.154
0.24	0.506	0.286	0.148
0.228	0.495	0.295	0.153
0.228	0.477	0.311	0.155
0.216	0.467	0.319	0.16
0.2	0.452	0.338	0.16
0.196	0.444	0.352	0.155
0.192	0.441	0.358	0.153
0.184	0.427	0.377	0.15
0.184	0.421	0.385	0.148
0.176	0.415	0.392	0.149
0.168	0.399	0.406	0.153
0.176	0.392	0.412	0.152
0.184	0.375	0.419	0.16
0.176	0.363	0.428	0.165
0.18	0.353	0.437	0.165
0.18	0.349	0.439	0.167
0.18	0.341	0.441	0.173
0.176	0.332	0.445	0.179
0.144	0.324	0.455	0.185
0.14	0.311	0.458	0.196
0.148	0.299	0.462	0.202
0.152	0.299	0.454	0.209
0.168	0.294	0.454	0.21
0.168	0.28	0.459	0.219
0.18	0.278	0.449	0.228
0.18	0.265	0.446	0.244
0.18	0.262	0.439	0.254
0.184	0.255	0.44	0.259
0.2	0.254	0.432	0.264
0.2	0.251	0.43	0.269
0.2	0.243	0.433	0.274
0.196	0.239	0.43	0.282
0.208	0.235	0.431	0.282
0.196	0.235	0.432	0.284
0.184	0.231	0.43	0.293
0.172	0.223	0.428	0.306
0.164	0.215	0.43	0.314
0.152	0.217	0.426	0.319
0.148	0.202	0.422	0.339
0.148	0.194	0.422	0.347
0.148	0.188	0.415	0.36
0.16	0.18	0.41	0.37
0.16	0.171	0.405	0.384
0.152	0.169	0.401	0.392
0.152	0.17	0.396	0.396
0.168	0.167	0.391	0.4
0.172	0.164	0.378	0.415
0.168	0.165	0.368	0.425
0.18	0.163	0.364	0.428
0.188	0.166	0.35	0.437
0.184	0.166	0.344	0.444
0.188	0.15	0.349	0.454
0.204	0.152	0.347	0.45
0.212	0.151	0.34	0.456
0.208	0.155	0.331	0.462
0.22	0.153	0.325	0.467
0.236	0.148	0.319	0.474
0.24	0.147	0.312	0.481
0.26	0.146	0.304	0.485
0.28	0.151	0.3	0.479
0.284	0.156	0.298	0.475
0.284	0.154	0.302	0.473
0.3	0.154	0.293	0.478
0.324	0.164	0.284	0.471
0.316	0.168	0.275	0.478
0.328	0.173	0.268	0.477
0.352	0.178	0.258	0.476
0.38	0.177	0.255	0.473
0.388	0.177	0.256	0.47
0.404	0.181	0.249	0.469
0.42	0.187	0.243	0.465
0.428	0.187	0.244	0.462
0.444	0.184	0.238	0.467
0.472	0.192	0.228	0.462
0.496	0.198	0.22	0.458
0.504	0.214	0.222	0.438
0.528	0.226	0.223	0.419
0.524	0.245	0.22	0.404
0.52	0.251	0.221	0.398
0.508	0.264	0.221	0.388
0.516	0.28	0.215	0.376
0.52	0.288	0.205	0.377
0.532	0.3	0.199	0.368
0.512	0.308	0.196	0.368
0.52	0.317	0.198	0.355
0.528	0.332	0.193	0.343
0.524	0.346	0.189	0.334
0.516	0.364	0.186	0.321
0.516	0.367	0.188	0.316
0.488	0.387	0.178	0.313
0.496	0.404	0.178	0.294
0.476	0.426	0.174	0.281
0.44	0.439	0.18	0.271
0.42	0.452	0.183	0.26
0.408	0.454	0.187	0.257
0.36	0.473	0.182	0.255
0.34	0.474	0.193	0.248
0.352	0.479	0.198	0.235
0.348	0.484	0.201	0.228
0.34	0.484	0.21	0.221
0.34	0.482	0.217	0.216
0.32	0.489	0.22	0.211
0.32	0.488	0.228	0.204
0.304	0.48	0.24	0.204
0.304	0.477	0.247	0.2
0.3	0.478	0.255	0.192
0.292	0.477	0.264	0.186
0.296	0.465	0.273	0.188
0.292	0.463	0.281	0.183
0.284	0.453	0.296	0.18
0.264	0.46	0.302	0.172
0.248	0.456	0.31	0.172
0.252	0.436	0.329	0.172
0.224	0.44	0.331	0.173
0.216	0.425	0.34	0.181
0.204	0.418	0.344	0.187
0.2	0.4	0.361	0.189
0.2	0.399	0.361	0.19
0.196	0.384	0.377	0.19
0.18	0.385	0.379	0.191
0.176	0.371	0.391	0.194
0.188	0.352	0.4	0.201
0.18	0.335	0.418	0.202
0.188	0.328	0.411	0.214
0.192	0.317	0.414	0.221
0.188	0.31	0.415	0.228
0.184	0.302	0.416	0.236
0.176	0.293	0.418	0.245
0.188	0.283	0.416	0.254
0.192	0.271	0.42	0.261
0.184	0.279	0.409	0.266
0.188	0.271	0.414	0.268
0.188	0.261	0.421	0.271
0.184	0.26	0.423	0.271
0.212	0.249	0.417	0.281
0.22	0.245	0.412	0.288
0.22	0.238	0.411	0.296
0.22	0.239	0.405	0.301
0.22	0.235	0.403	0.307
0.216	0.235	0.404	0.307
0.24	0.226	0.403	0.311
0.236	0.228	0.401	0.312
0.252	0.238	0.385	0.314
0.248	0.234	0.379	0.325
0.264	0.236	0.375	0.323
0.264	0.241	0.366	0.327
0.268	0.237	0.367	0.329
0.268	0.23	0.367	0.336
0.26	0.236	0.361	0.338
0.26	0.227	0.358	0.35
0.272	0.233	0.356	0.343
0.276	0.237	0.348	0.346
0.292	0.233	0.347	0.347
0.308	0.231	0.334	0.358
0.324	0.23	0.329	0.36
0.324	0.24	0.321	0.358
0.328	0.235	0.326	0.357
0.332	0.239	0.328	0.35
0.34	0.242	0.32	0.353
0.34	0.245	0.32	0.35
0.344	0.247	0.318	0.349
0.348	0.248	0.309	0.356
0.356	0.258	0.306	0.347
0.392	0.257	0.307	0.338
0.404	0.25	0.31	0.339
0.392	0.257	0.31	0.335
0.416	0.262	0.3	0.334
0.4	0.27	0.301	0.329
0.364	0.279	0.297	0.333
0.368	0.276	0.294	0.338
0.4	0.286	0.278	0.336
0.392	0.298	0.267	0.337
0.404	0.298	0.264	0.337
0.404	0.304	0.263	0.332
0.412	0.312	0.26	0.325
0.436	0.311	0.257	0.323
0.444	0.308	0.257	0.324
0.44	0.321	0.257	0.312
0.456	0.323	0.255	0.308
0.428	0.339	0.257	0.297
0.416	0.341	0.262	0.293
0.4	0.351	0.265	0.284
0.408	0.355	0.263	0.28
0.404	0.362	0.261	0.276
0.396	0.364	0.263	0.274
0.38	0.368	0.265	0.272
0.376	0.362	0.27	0.274
0.392	0.369	0.269	0.264
0.384	0.378	0.271	0.255
0.356	0.376	0.279	0.256
0.352	0.386	0.274	0.252
0.344	0.386	0.275	0.253
0.308	0.382	0.284	0.257
0.296	0.376	0.297	0.253
0.288	0.375	0.302	0.251
0.296	0.369	0.31	0.247
0.28	0.364	0.319	0.247
0.276	0.359	0.325	0.247
0.26	0.348	0.333	0.254
0.264	0.337	0.339	0.258
0.256	0.335	0.344	0.257
0.236	0.34	0.338	0.263
0.224	0.344	0.342	0.258
0.22	0.334	0.35	0.261
0.224	0.33	0.35	0.264
0.216	0.323	0.353	0.27
0.216	0.315	0.357	0.274
0.224	0.314	0.357	0.273
0.22	0.308	0.358	0.279
0.216	0.303	0.361	0.282
0.224	0.309	0.359	0.276
0.22	0.301	0.363	0.281
0.208	0.304	0.361	0.283
0.208	0.299	0.366	0.283
0.196	0.296	0.367	0.288
0.188	0.293	0.37	0.29
0.184	0.288	0.372	0.294
0.172	0.29	0.373	0.294
0.184	0.277	0.38	0.297
0.18	0.278	0.373	0.304
0.184	0.263	0.375	0.316
0.184	0.256	0.381	0.317
0.18	0.257	0.371	0.327
0.18	0.249	0.383	0.323
0.184	0.248	0.377	0.329
0.192	0.24	0.377	0.335
0.204	0.241	0.375	0.333
0.216	0.229	0.379	0.338
0.22	0.217	0.381	0.347
0.232	0.208	0.371	0.363
0.228	0.199	0.375	0.369
0.232	0.193	0.362	0.387
0.252	0.189	0.349	0.399
0.256	0.197	0.342	0.397
0.256	0.199	0.342	0.395
0.264	0.199	0.332	0.403
0.264	0.199	0.324	0.411
0.276	0.2	0.317	0.414
0.3	0.2	0.305	0.42
0.312	0.197	0.304	0.421
0.328	0.201	0.295	0.422
0.328	0.202	0.295	0.421
0.348	0.199	0.286	0.428
0.352	0.197	0.284	0.431
0.372	0.194	0.277	0.436
0.396	0.191	0.278	0.432
0.408	0.19	0.275	0.433
0.428	0.188	0.277	0.428
0.44	0.196	0.272	0.422
0.46	0.201	0.262	0.422
0.448	0.214	0.247	0.427
0.464	0.226	0.237	0.421
0.484	0.231	0.238	0.41
0.488	0.244	0.235	0.399
0.532	0.252	0.229	0.386
0.536	0.263	0.224	0.379
0.544	0.268	0.222	0.374
0.548	0.281	0.214	0.368
0.552	0.289	0.216	0.357
0.572	0.306	0.211	0.34
0.584	0.322	0.211	0.321
0.572	0.337	0.213	0.307
0.568	0.352	0.21	0.296
0.56	0.355	0.217	0.288
0.532	0.366	0.214	0.287
0.544	0.376	0.215	0.273
0.548	0.382	0.213	0.268
0.58	0.391	0.21	0.254
0.556	0.405	0.219	0.237
0.54	0.421	0.215	0.229
0.548	0.42	0.221	0.222
0.504	0.448	0.213	0.213
0.492	0.466	0.214	0.197
0.464	0.474	0.219	0.191
0.424	0.491	0.219	0.184
0.4	0.51	0.22	0.17
0.364	0.515	0.228	0.166
0.376	0.507	0.243	0.156
0.356	0.517	0.241	0.153
0.356	0.512	0.25	0.149
0.348	0.51	0.256	0.147
0.34	0.505	0.276	0.134
0.324	0.515	0.273	0.131
0.304	0.511	0.28	0.133
0.28	0.499	0.298	0.133
0.272	0.492	0.303	0.137
0.264	0.483	0.309	0.142
0.244	0.482	0.317	0.14
0.228	0.464	0.341	0.138
0.208	0.449	0.361	0.138
0.2	0.445	0.368	0.137
0.188	0.442	0.376	0.135
0.18	0.419	0.4	0.136
0.172	0.404	0.419	0.134
0.164	0.395	0.424	0.14
0.16	0.379	0.435	0.146
0.164	0.366	0.443	0.15
0.164	0.356	0.445	0.158
0.164	0.348	0.451	0.16
0.168	0.342	0.452	0.164
0.16	0.335	0.458	0.167
0.14	0.324	0.468	0.173
0.148	0.314	0.474	0.175
0.16	0.309	0.478	0.173
0.152	0.303	0.484	0.175
0.152	0.307	0.469	0.186
0.148	0.298	0.473	0.192
0.144	0.289	0.477	0.198
0.14	0.283	0.475	0.207
0.14	0.272	0.48	0.213
0.14	0.264	0.478	0.223
0.132	0.256	0.481	0.23
0.136	0.25	0.478	0.238
0.136	0.243	0.473	0.25
0.156	0.231	0.465	0.265
0.148	0.228	0.465	0.27
0.144	0.217	0.467	0.28
0.152	0.206	0.465	0.291
0.148	0.205	0.449	0.309
0.144	0.197	0.447	0.32
0.144	0.193	0.44	0.331
0.14	0.188	0.436	0.341
0.152	0.183	0.423	0.356
0.152	0.184	0.404	0.374
0.156	0.173	0.405	0.383
0.16	0.166	0.4	0.394
0.164	0.159	0.395	0.405
0.152	0.159	0.391	0.412
0.16	0.155	0.387	0.418
0.168	0.152	0.382	0.424
0.168	0.149	0.376	0.433
0.192	0.153	0.361	0.438
0.2	0.146	0.362	0.442
0.212	0.142	0.358	0.447
0.22	0.139	0.348	0.458
0.244	0.137	0.342	0.46
0.248	0.137	0.338	0.463
0.248	0.139	0.331	0.468
0.252	0.141	0.32	0.476
0.268	0.149	0.31	0.474
0.268	0.151	0.298	0.484
0.276	0.148	0.291	0.492
0.288	0.151	0.284	0.493
0.304	0.147	0.279	0.498
0.3	0.156	0.27	0.499
0.332	0.158	0.265	0.494
0.344	0.162	0.252	0.5
0.348	0.158	0.249	0.506
0.352	0.164	0.245	0.503
0.344	0.166	0.241	0.507
0.348	0.167	0.232	0.514
0.36	0.171	0.229	0.51
0.392	0.174	0.224	0.504
0.388	0.185	0.224	0.494
0.424	0.187	0.217	0.49
0.42	0.187	0.208	0.5
0.424	0.194	0.204	0.496
0.432	0.198	0.201	0.493
0.464	0.202	0.195	0.487
0.488	0.208	0.189	0.481
0.492	0.212	0.189	0.476
0.52	0.215	0.195	0.46
0.54	0.222	0.19	0.453
0.588	0.228	0.19	0.435
0.62	0.231	0.191	0.423
0.648	0.244	0.193	0.401
0.628	0.266	0.192	0.385
0.632	0.28	0.182	0.38
0.64	0.296	0.179	0.365
0.62	0.322	0.179	0.344
0.608	0.342	0.184	0.322
0.6	0.358	0.184	0.308
0.588	0.374	0.183	0.296
0.596	0.389	0.179	0.283
0.584	0.408	0.181	0.265
0.596	0.412	0.185	0.254
0.56	0.429	0.19	0.241
0.524	0.461	0.184	0.224
0.512	0.467	0.187	0.218
0.5	0.48	0.189	0.206
0.472	0.496	0.194	0.192
0.476	0.497	0.199	0.185
0.432	0.507	0.21	0.175
0.42	0.512	0.216	0.167
0.4	0.528	0.212	0.16
0.38	0.526	0.224	0.155
0.38	0.527	0.224	0.154
0.36	0.519	0.239	0.152
0.352	0.513	0.247	0.152
0.332	0.51	0.259	0.148
0.316	0.508	0.267	0.146
0.288	0.504	0.279	0.145
0.28	0.491	0.294	0.145
0.272	0.488	0.3	0.144
0.248	0.482	0.31	0.146
0.24	0.476	0.326	0.138
0.224	0.463	0.338	0.143
0.2	0.443	0.36	0.147
0.188	0.427	0.379	0.147
0.184	0.406	0.4	0.148
0.176	0.409	0.406	0.141
0.168	0.401	0.414	0.143
0.156	0.385	0.427	0.149
0.144	0.368	0.439	0.157
0.136	0.358	0.445	0.163
0.128	0.344	0.455	0.169
0.112	0.335	0.463	0.174
0.12	0.325	0.468	0.177
0.12	0.321	0.466	0.183
0.116	0.306	0.471	0.194
0.112	0.293	0.475	0.204
0.1	0.278	0.481	0.216
0.104	0.267	0.49	0.217
0.096	0.251	0.492	0.233
0.1	0.244	0.49	0.241
0.1	0.233	0.487	0.255
0.108	0.223	0.482	0.268
0.112	0.213	0.476	0.283
0.108	0.2	0.47	0.303
0.112	0.194	0.463	0.315
0.1	0.184	0.462	0.329
0.1	0.18	0.458	0.337
0.092	0.174	0.451	0.352
0.088	0.175	0.438	0.365
0.096	0.163	0.439	0.374
0.104	0.152	0.43	0.392
0.108	0.15	0.418	0.405
0.104	0.148	0.41	0.416
0.096	0.147	0.406	0.423
0.096	0.144	0.392	0.44
0.092	0.143	0.377	0.457
0.088	0.138	0.378	0.462
0.084	0.138	0.371	0.47
0.084	0.139	0.361	0.479
0.092	0.14	0.356	0.481
0.096	0.133	0.35	0.493
0.108	0.135	0.334	0.504
0.124	0.127	0.326	0.516
0.128	0.128	0.317	0.523
0.14	0.125	0.303	0.537
0.136	0.122	0.299	0.545
0.152	0.118	0.293	0.551
0.168	0.113	0.287	0.558
0.176	0.115	0.279	0.562
0.2	0.108	0.272	0.57
0.208	0.105	0.268	0.575
0.22	0.106	0.264	0.575
0.228	0.107	0.264	0.572
0.244	0.104	0.254	0.581
0.264	0.1	0.253	0.581
0.28	0.102	0.247	0.581
0.316	0.1	0.234	0.587
0.328	0.097	0.231	0.59
0.336	0.104	0.224	0.588
0.356	0.106	0.22	0.585
0.384	0.105	0.215	0.584
0.388	0.1	0.215	0.588
0.416	0.111	0.214	0.571
0.452	0.116	0.204	0.567
0.484	0.116	0.199	0.564
0.508	0.12	0.195	0.558
0.552	0.127	0.193	0.542
0.608	0.135	0.186	0.527
0.632	0.148	0.176	0.518
0.668	0.156	0.172	0.505
0.684	0.166	0.168	0.495
];
a=[0.128	0.334	0.481	0.153
0.124	0.316	0.497	0.156
0.116	0.303	0.503	0.165
0.104	0.292	0.509	0.173
0.1	0.283	0.513	0.179
0.1	0.271	0.513	0.191
0.1	0.255	0.515	0.205
0.096	0.244	0.513	0.219
0.096	0.23	0.517	0.229
0.112	0.23	0.508	0.234
0.112	0.228	0.503	0.241
0.112	0.214	0.503	0.255
0.104	0.21	0.498	0.266
0.108	0.197	0.498	0.278
0.116	0.185	0.483	0.303
0.128	0.173	0.475	0.32
0.132	0.17	0.461	0.336
0.14	0.17	0.452	0.343
0.14	0.172	0.44	0.353
0.152	0.168	0.429	0.365
0.16	0.158	0.423	0.379
0.164	0.154	0.417	0.388
0.192	0.143	0.405	0.404
0.2	0.141	0.401	0.408
0.212	0.137	0.39	0.42
0.232	0.132	0.379	0.431
0.232	0.132	0.366	0.444
0.232	0.127	0.358	0.457
0.252	0.125	0.341	0.471
0.268	0.127	0.333	0.473
0.288	0.13	0.324	0.474
0.284	0.131	0.317	0.481
0.296	0.124	0.303	0.499
0.312	0.122	0.291	0.509
0.344	0.122	0.284	0.508
0.332	0.127	0.279	0.511
0.344	0.133	0.276	0.505
0.36	0.13	0.272	0.508
0.4	0.131	0.26	0.509
0.396	0.138	0.253	0.51
0.416	0.132	0.251	0.513
0.432	0.133	0.238	0.521
0.48	0.14	0.234	0.506
0.476	0.143	0.228	0.51
0.504	0.149	0.224	0.501
0.52	0.155	0.212	0.503
0.54	0.163	0.21	0.492
0.548	0.171	0.206	0.486
0.56	0.183	0.203	0.474
0.568	0.192	0.199	0.467
0.576	0.206	0.195	0.455
0.608	0.207	0.192	0.449
0.636	0.223	0.187	0.431
0.636	0.239	0.189	0.413
0.644	0.254	0.183	0.402
0.656	0.272	0.184	0.38
0.668	0.286	0.182	0.365
0.688	0.292	0.185	0.351
0.684	0.314	0.179	0.336
0.664	0.331	0.179	0.324
0.676	0.338	0.179	0.314
0.64	0.363	0.174	0.303
0.612	0.385	0.18	0.282
0.604	0.408	0.18	0.261
0.6	0.424	0.181	0.245
0.572	0.444	0.182	0.231
0.528	0.459	0.192	0.217
0.508	0.472	0.192	0.209
0.488	0.482	0.195	0.201
0.504	0.485	0.202	0.187
0.496	0.487	0.206	0.183
0.484	0.503	0.202	0.174
0.46	0.51	0.206	0.169
0.432	0.517	0.217	0.158
0.404	0.529	0.222	0.148
0.392	0.535	0.228	0.139
0.388	0.523	0.243	0.137
0.372	0.525	0.247	0.135
0.368	0.521	0.254	0.133
0.348	0.527	0.263	0.123
0.328	0.525	0.273	0.12
0.304	0.537	0.265	0.122
0.3	0.539	0.27	0.116
0.296	0.529	0.281	0.116
0.288	0.523	0.29	0.115
0.272	0.519	0.301	0.112
0.268	0.516	0.307	0.11
0.252	0.508	0.321	0.108
0.232	0.507	0.324	0.111
0.208	0.495	0.34	0.113
0.192	0.483	0.353	0.116
0.184	0.469	0.367	0.118
0.168	0.452	0.392	0.114
0.152	0.436	0.413	0.113
0.16	0.423	0.421	0.116
0.152	0.404	0.438	0.12
0.144	0.397	0.442	0.125
0.14	0.379	0.451	0.135
0.14	0.368	0.454	0.143
0.14	0.341	0.468	0.156
0.128	0.322	0.48	0.166
0.136	0.318	0.476	0.172
0.132	0.302	0.483	0.182
0.112	0.291	0.493	0.188
0.104	0.272	0.501	0.201
0.104	0.265	0.496	0.213
0.1	0.256	0.494	0.225
0.096	0.248	0.495	0.233
0.088	0.239	0.491	0.248
0.088	0.231	0.487	0.26
0.088	0.215	0.49	0.273
0.084	0.206	0.482	0.291
0.084	0.201	0.475	0.303
0.084	0.194	0.477	0.308
0.076	0.188	0.476	0.317
0.076	0.185	0.468	0.328
0.076	0.177	0.457	0.347
0.08	0.172	0.454	0.354
0.08	0.165	0.446	0.369
0.084	0.16	0.435	0.384
0.092	0.154	0.422	0.401
0.088	0.146	0.412	0.42
0.08	0.141	0.397	0.442
0.092	0.127	0.387	0.463
0.092	0.127	0.37	0.48
0.096	0.122	0.361	0.493
0.112	0.119	0.348	0.505
0.116	0.115	0.34	0.516
0.116	0.112	0.327	0.532
0.124	0.103	0.324	0.542
0.128	0.097	0.31	0.561
0.16	0.093	0.297	0.57
0.168	0.093	0.289	0.576
0.192	0.091	0.278	0.583
0.208	0.093	0.273	0.582
0.22	0.093	0.27	0.582
0.24	0.088	0.257	0.595
0.264	0.087	0.247	0.6
0.284	0.087	0.245	0.597
0.324	0.089	0.236	0.594
0.356	0.088	0.23	0.593
0.384	0.085	0.221	0.598
0.404	0.085	0.22	0.594
0.456	0.087	0.213	0.586
0.46	0.092	0.205	0.588
0.496	0.096	0.2	0.58
0.528	0.102	0.19	0.576
0.552	0.11	0.182	0.57
0.584	0.112	0.177	0.565
0.612	0.119	0.169	0.559
0.66	0.127	0.161	0.547
0.684	0.137	0.158	0.534
0.724	0.148	0.154	0.517
0.788	0.154	0.146	0.503
0.832	0.156	0.138	0.498
0.884	0.166	0.135	0.478
0.872	0.188	0.129	0.465
0.876	0.213	0.124	0.444
0.9	0.235	0.123	0.417
0.904	0.254	0.123	0.397
0.868	0.294	0.115	0.374
0.852	0.329	0.114	0.344
0.844	0.354	0.116	0.319
0.832	0.379	0.118	0.295
0.816	0.396	0.12	0.28
0.796	0.426	0.125	0.25
0.748	0.447	0.129	0.237
0.708	0.477	0.13	0.216
0.664	0.488	0.135	0.211
0.612	0.517	0.136	0.194
0.588	0.537	0.141	0.175
0.564	0.552	0.143	0.164
0.52	0.563	0.15	0.157
0.496	0.58	0.152	0.144
0.472	0.602	0.145	0.135
0.464	0.612	0.148	0.124
0.416	0.621	0.152	0.123
0.376	0.621	0.164	0.121
0.364	0.618	0.178	0.113
0.344	0.62	0.181	0.113
0.34	0.616	0.188	0.111
0.324	0.611	0.196	0.112
0.308	0.603	0.211	0.109
0.272	0.602	0.225	0.105
0.272	0.591	0.237	0.104
0.256	0.582	0.257	0.097
0.24	0.573	0.272	0.095
0.228	0.553	0.296	0.094
0.22	0.535	0.313	0.097
0.216	0.517	0.337	0.092
0.204	0.501	0.35	0.098
0.204	0.494	0.358	0.097
0.188	0.482	0.371	0.1
0.184	0.471	0.384	0.099
0.18	0.458	0.397	0.1
0.176	0.451	0.403	0.102
0.172	0.441	0.414	0.102
0.172	0.432	0.422	0.103
0.172	0.424	0.427	0.106
0.148	0.416	0.437	0.11
0.148	0.403	0.446	0.114
0.144	0.39	0.457	0.117
0.144	0.37	0.471	0.123
0.144	0.353	0.484	0.127
0.14	0.345	0.487	0.133
0.136	0.325	0.494	0.147
0.132	0.31	0.503	0.154
0.128	0.293	0.513	0.162
0.116	0.277	0.52	0.174
0.112	0.265	0.521	0.186
0.112	0.265	0.514	0.193
0.116	0.255	0.514	0.202
0.12	0.255	0.508	0.207
0.116	0.246	0.504	0.221
0.12	0.243	0.499	0.228
0.124	0.235	0.493	0.241
0.116	0.225	0.496	0.25
0.112	0.219	0.489	0.264
0.108	0.213	0.473	0.287
0.096	0.21	0.467	0.299
0.1	0.205	0.464	0.306
0.108	0.201	0.457	0.315
0.112	0.199	0.449	0.324
0.12	0.191	0.446	0.333
0.136	0.19	0.444	0.332
0.132	0.18	0.445	0.342
0.136	0.174	0.436	0.356
0.14	0.169	0.433	0.363
0.144	0.17	0.426	0.368
0.164	0.158	0.424	0.377
0.172	0.15	0.419	0.388
0.172	0.147	0.407	0.403
0.168	0.148	0.39	0.42
0.172	0.144	0.38	0.433
0.172	0.138	0.371	0.448
0.18	0.133	0.355	0.467
0.192	0.133	0.342	0.477
0.208	0.134	0.333	0.481
0.228	0.134	0.322	0.487
0.248	0.132	0.308	0.498
0.256	0.131	0.306	0.499
0.264	0.127	0.302	0.505
0.288	0.124	0.295	0.509
0.304	0.131	0.289	0.504
0.324	0.127	0.289	0.503
0.34	0.126	0.281	0.508
0.364	0.127	0.275	0.507
0.392	0.126	0.265	0.511
0.404	0.131	0.263	0.505
0.416	0.134	0.26	0.502
0.464	0.142	0.251	0.491
0.476	0.143	0.243	0.495
0.496	0.145	0.241	0.49
0.496	0.15	0.233	0.493
0.52	0.157	0.229	0.484
0.544	0.165	0.221	0.478
0.576	0.165	0.207	0.484
0.628	0.167	0.201	0.475
0.648	0.177	0.195	0.466
0.664	0.187	0.191	0.456
0.704	0.186	0.186	0.452
0.728	0.204	0.178	0.436
0.752	0.206	0.174	0.432
0.724	0.232	0.17	0.417
0.732	0.238	0.173	0.406
0.748	0.256	0.166	0.391
0.72	0.279	0.156	0.385
0.736	0.297	0.156	0.363
0.7	0.328	0.15	0.347
0.716	0.343	0.145	0.333
0.712	0.361	0.143	0.318
0.652	0.379	0.152	0.306
0.628	0.398	0.151	0.294
0.584	0.416	0.155	0.283
0.548	0.436	0.161	0.266
0.536	0.45	0.16	0.256
0.508	0.468	0.162	0.243
0.472	0.477	0.168	0.237
0.452	0.485	0.178	0.224
0.416	0.499	0.181	0.216
0.416	0.502	0.183	0.211
0.42	0.501	0.191	0.203
0.42	0.496	0.203	0.196
0.392	0.511	0.205	0.186
0.372	0.522	0.212	0.173
0.352	0.524	0.219	0.169
0.336	0.521	0.229	0.166
0.32	0.524	0.238	0.158
0.308	0.516	0.25	0.157
0.296	0.516	0.254	0.156
0.272	0.516	0.264	0.152
0.268	0.514	0.274	0.145
0.264	0.503	0.287	0.144
0.268	0.483	0.308	0.142
0.256	0.481	0.317	0.138
0.244	0.473	0.328	0.138
0.24	0.462	0.338	0.14
0.22	0.465	0.341	0.139
0.212	0.455	0.354	0.138
0.192	0.449	0.365	0.138
0.188	0.437	0.373	0.143
0.188	0.431	0.384	0.138
0.184	0.413	0.395	0.146
0.172	0.409	0.397	0.151
0.176	0.4	0.405	0.151
0.168	0.39	0.411	0.157
0.16	0.384	0.414	0.162
0.152	0.368	0.429	0.165
0.148	0.36	0.435	0.168
0.144	0.353	0.441	0.17
0.144	0.354	0.436	0.174
0.148	0.343	0.435	0.185
0.144	0.334	0.435	0.195
0.14	0.325	0.435	0.205
0.14	0.312	0.44	0.213
0.136	0.307	0.443	0.216
0.128	0.294	0.451	0.223
0.12	0.29	0.459	0.221
0.12	0.284	0.457	0.229
0.124	0.273	0.467	0.229
0.128	0.264	0.455	0.249
0.12	0.256	0.455	0.259
0.116	0.251	0.456	0.264
0.116	0.244	0.454	0.273
0.12	0.228	0.456	0.286
0.116	0.219	0.454	0.298
0.128	0.213	0.447	0.308
0.124	0.206	0.444	0.319
0.128	0.199	0.443	0.326
0.144	0.189	0.428	0.347
0.148	0.181	0.426	0.356
0.152	0.182	0.417	0.363
0.148	0.167	0.405	0.391
0.144	0.165	0.395	0.404
0.14	0.161	0.39	0.414
0.144	0.163	0.38	0.421
0.16	0.17	0.369	0.421
0.16	0.17	0.357	0.433
0.172	0.175	0.344	0.438
0.172	0.17	0.34	0.447
0.188	0.167	0.328	0.458
0.2	0.163	0.321	0.466
0.224	0.16	0.31	0.474
0.244	0.154	0.304	0.481
0.264	0.15	0.298	0.486
0.292	0.151	0.285	0.491
0.308	0.147	0.274	0.502
0.336	0.148	0.266	0.502
0.332	0.15	0.263	0.504
0.364	0.158	0.256	0.495
0.38	0.157	0.259	0.489
0.42	0.148	0.256	0.491
0.432	0.155	0.245	0.492
0.424	0.164	0.239	0.491
0.444	0.177	0.229	0.483
0.464	0.183	0.219	0.482
0.456	0.193	0.216	0.477
0.464	0.195	0.217	0.472
0.488	0.199	0.215	0.464
0.496	0.21	0.208	0.458
0.504	0.219	0.201	0.454
0.536	0.235	0.193	0.438
0.544	0.246	0.19	0.428
0.592	0.253	0.187	0.412
0.64	0.262	0.18	0.398
0.632	0.281	0.179	0.382
0.62	0.295	0.181	0.369
0.624	0.316	0.173	0.355
0.612	0.343	0.176	0.328
0.584	0.362	0.18	0.312
0.584	0.369	0.183	0.302
0.56	0.394	0.181	0.285
0.536	0.414	0.178	0.274
0.56	0.429	0.178	0.253
0.54	0.444	0.187	0.234
0.528	0.461	0.191	0.216
0.516	0.46	0.197	0.214
0.508	0.468	0.197	0.208
0.464	0.481	0.206	0.197
0.46	0.481	0.211	0.193
0.452	0.481	0.219	0.187
0.46	0.48	0.224	0.181
0.456	0.486	0.229	0.171
0.432	0.5	0.227	0.165
0.42	0.508	0.228	0.159
0.384	0.524	0.22	0.16
0.372	0.525	0.226	0.156
0.36	0.531	0.229	0.15
0.336	0.527	0.244	0.145
0.332	0.509	0.269	0.139
0.316	0.51	0.279	0.132
0.28	0.507	0.291	0.132
0.268	0.501	0.298	0.134
0.264	0.502	0.303	0.129
0.268	0.492	0.317	0.124
0.244	0.48	0.338	0.121
0.236	0.47	0.35	0.121
0.22	0.463	0.361	0.121
0.212	0.451	0.372	0.124
0.192	0.433	0.392	0.127
0.188	0.422	0.403	0.128
0.192	0.411	0.41	0.131
0.184	0.405	0.413	0.136
0.18	0.392	0.422	0.141
0.176	0.383	0.425	0.148
0.164	0.363	0.439	0.157
0.156	0.359	0.44	0.162
0.164	0.353	0.444	0.162
0.16	0.348	0.444	0.168
0.16	0.341	0.442	0.177
0.184	0.329	0.443	0.182
0.184	0.323	0.442	0.189
0.184	0.316	0.449	0.189
0.18	0.312	0.453	0.19
0.164	0.312	0.452	0.195
0.164	0.306	0.451	0.202
0.172	0.292	0.449	0.216
0.18	0.281	0.452	0.222
0.184	0.271	0.445	0.238
0.188	0.26	0.447	0.246
0.188	0.255	0.439	0.259
0.2	0.247	0.427	0.276
0.208	0.241	0.422	0.285
0.212	0.246	0.408	0.293
0.208	0.25	0.392	0.306
0.212	0.256	0.379	0.312
0.216	0.248	0.374	0.324
0.228	0.25	0.367	0.326
0.22	0.246	0.361	0.338
0.232	0.243	0.357	0.342
0.24	0.237	0.363	0.34
0.228	0.237	0.356	0.35
0.248	0.236	0.351	0.351
0.256	0.232	0.346	0.358
0.272	0.224	0.347	0.361
0.256	0.227	0.342	0.367
0.268	0.222	0.335	0.376
0.268	0.227	0.319	0.387
0.272	0.228	0.314	0.39
0.296	0.235	0.305	0.386
0.292	0.23	0.3	0.397
0.292	0.229	0.292	0.406
0.312	0.225	0.289	0.408
0.32	0.234	0.286	0.4
0.316	0.235	0.283	0.403
0.292	0.238	0.282	0.407
0.28	0.246	0.287	0.397
0.284	0.25	0.29	0.389
0.284	0.248	0.291	0.39
0.288	0.256	0.284	0.388
0.304	0.256	0.286	0.382
0.3	0.257	0.287	0.381
0.308	0.26	0.283	0.38
0.312	0.27	0.275	0.377
0.312	0.262	0.272	0.388
0.3	0.252	0.284	0.389
0.288	0.259	0.289	0.38
0.276	0.27	0.284	0.377
0.276	0.272	0.29	0.369
0.264	0.284	0.286	0.364
0.268	0.277	0.289	0.367
0.276	0.278	0.289	0.364
0.284	0.276	0.289	0.364
0.276	0.277	0.282	0.372
0.284	0.277	0.281	0.371
0.3	0.275	0.286	0.364
0.308	0.291	0.28	0.352
0.296	0.292	0.282	0.352
0.316	0.297	0.276	0.348
0.328	0.287	0.285	0.346
0.32	0.289	0.282	0.349
0.34	0.286	0.278	0.351
0.336	0.284	0.276	0.356
0.336	0.284	0.278	0.354
0.336	0.289	0.276	0.351
0.352	0.291	0.28	0.341
0.344	0.288	0.281	0.345
0.348	0.291	0.276	0.346
0.352	0.296	0.276	0.34
0.372	0.3	0.272	0.335
0.38	0.307	0.264	0.334
0.392	0.302	0.264	0.336
0.404	0.305	0.267	0.327
0.412	0.306	0.272	0.319
0.396	0.313	0.276	0.312
0.4	0.316	0.285	0.299
0.384	0.32	0.292	0.292
0.372	0.326	0.29	0.291
0.368	0.326	0.294	0.288
0.368	0.333	0.298	0.277
0.372	0.335	0.295	0.277
0.376	0.329	0.307	0.27
0.368	0.324	0.308	0.276
0.372	0.328	0.297	0.282
0.4	0.333	0.288	0.279
0.388	0.331	0.286	0.286
0.392	0.33	0.286	0.286
0.396	0.323	0.285	0.293
0.388	0.329	0.283	0.291
0.372	0.333	0.288	0.286
0.368	0.342	0.285	0.281
0.356	0.347	0.286	0.278
0.332	0.361	0.282	0.274
0.308	0.357	0.296	0.27
0.284	0.359	0.304	0.266
0.272	0.358	0.306	0.268
0.272	0.355	0.312	0.265
0.268	0.354	0.313	0.266
0.264	0.359	0.315	0.26
0.272	0.342	0.32	0.27
0.264	0.341	0.323	0.27
0.26	0.335	0.329	0.271
0.28	0.33	0.332	0.268
0.276	0.319	0.334	0.278
0.28	0.317	0.332	0.281
0.272	0.314	0.33	0.288
0.272	0.312	0.323	0.297
0.26	0.309	0.322	0.304
0.248	0.306	0.32	0.312
0.244	0.309	0.312	0.318
0.244	0.312	0.311	0.316
0.252	0.315	0.308	0.314
0.256	0.314	0.308	0.314
0.236	0.318	0.304	0.319
0.244	0.307	0.31	0.322
0.268	0.304	0.309	0.32
0.284	0.298	0.314	0.317
0.28	0.295	0.323	0.312
0.288	0.29	0.331	0.307
0.28	0.287	0.326	0.317
0.272	0.277	0.327	0.328
0.272	0.276	0.322	0.334
0.264	0.287	0.322	0.325
0.276	0.282	0.312	0.337
0.264	0.275	0.32	0.339
0.264	0.271	0.32	0.343
0.272	0.271	0.313	0.348
0.28	0.275	0.315	0.34
0.268	0.285	0.313	0.335
0.264	0.283	0.316	0.335
0.26	0.289	0.311	0.335
0.26	0.288	0.307	0.34
0.272	0.295	0.294	0.343
0.26	0.291	0.299	0.345
0.26	0.285	0.307	0.343
0.252	0.278	0.305	0.354
0.248	0.278	0.303	0.357
0.26	0.27	0.304	0.361
0.26	0.261	0.308	0.366
0.264	0.262	0.31	0.362
0.28	0.26	0.305	0.365
0.28	0.26	0.302	0.368
0.268	0.26	0.303	0.37
0.288	0.259	0.301	0.368
0.308	0.264	0.295	0.364
0.304	0.26	0.298	0.366
0.304	0.274	0.289	0.361
0.288	0.275	0.289	0.364
0.288	0.267	0.289	0.372
0.288	0.275	0.281	0.372
0.292	0.279	0.276	0.372
0.3	0.276	0.275	0.374
0.292	0.287	0.271	0.369
0.284	0.286	0.281	0.362
0.3	0.287	0.277	0.361
0.304	0.279	0.281	0.364
0.3	0.278	0.282	0.365
0.316	0.275	0.271	0.375
0.32	0.279	0.265	0.376
0.3	0.279	0.264	0.382
0.312	0.283	0.262	0.377
0.32	0.288	0.266	0.366
0.324	0.285	0.269	0.365
0.344	0.287	0.268	0.359
0.348	0.282	0.269	0.362
0.34	0.273	0.272	0.37
0.336	0.279	0.265	0.372
0.376	0.281	0.263	0.362
0.384	0.28	0.268	0.356
0.376	0.28	0.271	0.355
0.38	0.285	0.273	0.347
0.368	0.291	0.272	0.345
0.36	0.289	0.272	0.349
0.36	0.289	0.271	0.35
0.368	0.292	0.268	0.348
0.372	0.3	0.265	0.342
0.356	0.299	0.263	0.349
0.38	0.291	0.263	0.351
0.4	0.286	0.269	0.345
0.388	0.29	0.269	0.344
0.392	0.301	0.267	0.334
0.408	0.311	0.261	0.326
0.396	0.304	0.273	0.324
0.392	0.31	0.269	0.323
0.38	0.325	0.263	0.317
0.38	0.326	0.263	0.316
0.384	0.332	0.263	0.309
0.4	0.329	0.261	0.31
0.388	0.336	0.262	0.305
0.372	0.341	0.265	0.301
0.384	0.346	0.266	0.292
0.372	0.348	0.27	0.289
0.352	0.351	0.269	0.292
0.348	0.354	0.27	0.289
0.356	0.359	0.27	0.282
0.356	0.357	0.272	0.282
0.336	0.353	0.279	0.284
0.324	0.343	0.29	0.286
0.332	0.345	0.284	0.288
0.296	0.359	0.283	0.284
0.308	0.349	0.289	0.285
0.312	0.336	0.296	0.29
0.316	0.335	0.296	0.29
0.312	0.335	0.295	0.292
0.308	0.329	0.304	0.29
0.292	0.332	0.304	0.291
0.292	0.329	0.307	0.291
0.288	0.336	0.309	0.283
0.288	0.336	0.307	0.285
0.296	0.339	0.311	0.276
0.308	0.333	0.312	0.278
0.316	0.326	0.316	0.279
0.316	0.327	0.316	0.278
0.324	0.318	0.323	0.278
0.3	0.327	0.32	0.278
0.304	0.322	0.324	0.278
0.308	0.329	0.325	0.269
0.296	0.325	0.335	0.266
0.296	0.319	0.339	0.268
0.3	0.332	0.326	0.267
0.296	0.325	0.328	0.273
0.292	0.319	0.343	0.265
0.284	0.316	0.352	0.261
0.276	0.315	0.348	0.268
0.264	0.321	0.346	0.267
0.276	0.321	0.342	0.268
0.288	0.317	0.336	0.275
0.28	0.312	0.335	0.283
0.276	0.306	0.337	0.288
0.268	0.301	0.339	0.293
0.264	0.289	0.337	0.308
0.272	0.285	0.331	0.316
0.276	0.285	0.335	0.311
0.292	0.286	0.332	0.309
0.288	0.293	0.322	0.313
0.284	0.3	0.314	0.315
0.292	0.295	0.31	0.322
0.3	0.291	0.308	0.326
0.308	0.288	0.308	0.327
0.3	0.285	0.308	0.332
0.316	0.288	0.304	0.329
0.316	0.291	0.308	0.322
0.332	0.292	0.302	0.323
0.328	0.294	0.303	0.321
0.336	0.294	0.306	0.316
0.34	0.286	0.31	0.319
0.336	0.29	0.309	0.317
0.344	0.29	0.303	0.321
0.336	0.3	0.296	0.32
0.332	0.303	0.293	0.321
0.336	0.304	0.291	0.321
0.348	0.309	0.284	0.32
0.356	0.302	0.282	0.327
0.348	0.313	0.28	0.32
0.348	0.31	0.287	0.316
0.38	0.309	0.288	0.308
0.396	0.32	0.279	0.302
0.38	0.319	0.277	0.309
0.396	0.328	0.274	0.299
0.404	0.337	0.27	0.292
0.4	0.35	0.259	0.291
0.4	0.356	0.265	0.279
0.376	0.365	0.266	0.275
0.4	0.363	0.264	0.273
0.404	0.37	0.262	0.267
0.396	0.369	0.268	0.264
0.38	0.387	0.258	0.26
0.4	0.388	0.257	0.255
0.384	0.389	0.258	0.257
0.4	0.386	0.257	0.257
0.4	0.4	0.252	0.248
0.388	0.405	0.251	0.247
0.404	0.413	0.248	0.238
0.396	0.414	0.257	0.23
0.38	0.41	0.268	0.227
0.38	0.408	0.268	0.229
0.364	0.411	0.267	0.231
0.36	0.42	0.264	0.226
0.344	0.421	0.266	0.227
0.356	0.422	0.269	0.22
0.348	0.415	0.279	0.219
0.324	0.412	0.286	0.221
0.328	0.417	0.286	0.215
0.34	0.413	0.292	0.21
0.34	0.416	0.297	0.202
0.332	0.414	0.31	0.193
0.324	0.413	0.31	0.196
0.304	0.41	0.326	0.188
0.304	0.411	0.325	0.188
0.288	0.413	0.329	0.186
0.284	0.414	0.327	0.188
0.26	0.414	0.329	0.192
0.248	0.413	0.331	0.194
0.248	0.402	0.341	0.195
0.224	0.405	0.337	0.202
0.216	0.392	0.343	0.211
0.224	0.379	0.35	0.215
0.212	0.38	0.352	0.215
0.208	0.389	0.336	0.223
0.208	0.378	0.342	0.228
0.208	0.372	0.35	0.226
0.208	0.365	0.358	0.225
0.212	0.353	0.367	0.227
0.212	0.349	0.361	0.237
0.22	0.343	0.361	0.241
0.224	0.341	0.367	0.236
0.224	0.331	0.377	0.236
0.212	0.324	0.378	0.245
0.216	0.31	0.382	0.254
0.22	0.299	0.392	0.254
0.216	0.285	0.4	0.261
0.216	0.287	0.399	0.26
0.232	0.29	0.391	0.261
0.24	0.291	0.383	0.266
0.248	0.296	0.372	0.27
0.252	0.292	0.37	0.275
0.26	0.28	0.377	0.278
0.26	0.279	0.373	0.283
0.256	0.273	0.377	0.286
0.256	0.272	0.38	0.284
0.244	0.271	0.376	0.292
0.224	0.276	0.373	0.295
0.24	0.278	0.368	0.294
0.24	0.283	0.365	0.292
0.24	0.281	0.36	0.299
0.22	0.268	0.361	0.316
0.232	0.264	0.358	0.32
0.228	0.258	0.358	0.327
0.22	0.258	0.353	0.334
0.212	0.252	0.351	0.344
0.24	0.248	0.348	0.344
0.26	0.238	0.348	0.349
0.268	0.246	0.341	0.346
0.276	0.241	0.335	0.355
0.264	0.246	0.331	0.357
0.272	0.242	0.328	0.362
0.272	0.246	0.327	0.359
0.272	0.239	0.321	0.372
0.288	0.232	0.318	0.378
0.288	0.223	0.32	0.385
0.3	0.226	0.312	0.387
0.32	0.22	0.309	0.391
0.324	0.223	0.308	0.388
0.328	0.22	0.304	0.394
0.332	0.223	0.298	0.396
0.344	0.22	0.293	0.401
0.36	0.218	0.28	0.412
0.38	0.213	0.279	0.413
0.372	0.218	0.276	0.413
0.376	0.221	0.272	0.413
0.372	0.215	0.27	0.422
0.38	0.217	0.269	0.419
0.404	0.219	0.265	0.415
0.404	0.228	0.259	0.412
0.44	0.234	0.253	0.403
0.424	0.251	0.248	0.395
0.452	0.259	0.251	0.377
0.448	0.265	0.249	0.374
0.452	0.276	0.241	0.37
0.456	0.3	0.232	0.354
0.444	0.307	0.234	0.348
0.452	0.304	0.236	0.347
0.436	0.313	0.237	0.341
0.42	0.317	0.237	0.341
0.436	0.327	0.235	0.329
0.444	0.342	0.222	0.325
0.412	0.349	0.227	0.321
0.388	0.377	0.22	0.306
0.4	0.386	0.219	0.295
0.388	0.398	0.215	0.29
0.38	0.407	0.215	0.283
0.372	0.414	0.217	0.276
0.372	0.412	0.222	0.273
0.364	0.415	0.229	0.265
0.36	0.415	0.237	0.258
0.344	0.426	0.239	0.249
0.34	0.428	0.245	0.242
0.336	0.433	0.25	0.233
0.316	0.435	0.259	0.227
0.312	0.437	0.261	0.224
0.292	0.434	0.271	0.222
0.264	0.442	0.275	0.217
0.264	0.437	0.279	0.218
0.244	0.449	0.277	0.213
0.248	0.439	0.289	0.21
0.24	0.437	0.293	0.21
0.248	0.434	0.296	0.208
0.252	0.43	0.3	0.207
0.26	0.424	0.309	0.202
0.248	0.423	0.314	0.201
0.256	0.411	0.315	0.21
0.248	0.406	0.321	0.211
0.244	0.4	0.335	0.204
0.244	0.39	0.346	0.203
0.248	0.378	0.35	0.21
0.256	0.364	0.357	0.215
0.248	0.365	0.357	0.216
0.252	0.359	0.355	0.223
0.236	0.351	0.363	0.227
0.232	0.349	0.365	0.228
0.224	0.341	0.369	0.234
0.204	0.332	0.376	0.241
0.204	0.325	0.374	0.25
0.204	0.334	0.361	0.254
0.2	0.322	0.361	0.267
0.196	0.309	0.361	0.281
0.196	0.305	0.357	0.289
0.196	0.301	0.361	0.289
0.196	0.294	0.367	0.29
0.204	0.291	0.364	0.294
0.212	0.286	0.364	0.297
0.22	0.282	0.364	0.299
0.224	0.266	0.37	0.308
0.228	0.259	0.372	0.312
0.224	0.264	0.355	0.325
0.228	0.258	0.35	0.335
0.236	0.253	0.341	0.347
0.232	0.248	0.338	0.356
0.252	0.244	0.338	0.355
0.26	0.242	0.332	0.361
0.248	0.236	0.329	0.373
0.252	0.231	0.326	0.38
0.256	0.224	0.319	0.393
0.26	0.227	0.31	0.398
0.284	0.221	0.316	0.392
0.292	0.223	0.315	0.389
0.3	0.217	0.314	0.394
0.32	0.218	0.306	0.396
0.32	0.214	0.312	0.394
0.336	0.212	0.305	0.399
0.348	0.204	0.3	0.409
0.368	0.206	0.296	0.406
0.384	0.203	0.287	0.414
0.4	0.203	0.281	0.416
0.416	0.213	0.268	0.415
0.412	0.218	0.259	0.42
0.424	0.222	0.255	0.417
0.42	0.227	0.248	0.42
0.408	0.234	0.246	0.418
0.452	0.241	0.232	0.414
0.46	0.248	0.23	0.407
0.448	0.262	0.232	0.394
0.456	0.269	0.232	0.385
0.444	0.28	0.233	0.376
0.444	0.289	0.231	0.369
0.44	0.3	0.227	0.363
0.448	0.303	0.226	0.359
0.428	0.308	0.229	0.356
0.436	0.318	0.225	0.348
0.436	0.325	0.232	0.334
0.444	0.331	0.235	0.323
0.46	0.335	0.24	0.31
0.464	0.344	0.249	0.291
0.448	0.364	0.243	0.281
0.428	0.375	0.239	0.279
0.46	0.375	0.234	0.276
0.456	0.379	0.234	0.273
0.448	0.387	0.232	0.269
0.452	0.388	0.24	0.259
0.468	0.389	0.248	0.246
0.476	0.387	0.252	0.242
0.488	0.394	0.246	0.238
0.464	0.396	0.25	0.238
0.456	0.409	0.247	0.23
0.452	0.405	0.255	0.227
0.432	0.411	0.259	0.222
0.44	0.421	0.255	0.214
0.416	0.422	0.266	0.208
0.392	0.432	0.267	0.203
0.36	0.447	0.267	0.196
0.36	0.437	0.275	0.198
0.376	0.434	0.281	0.191
0.352	0.443	0.287	0.182
0.32	0.443	0.303	0.174
0.312	0.423	0.32	0.179
0.292	0.419	0.327	0.181
0.288	0.406	0.341	0.181
0.284	0.4	0.344	0.185
0.268	0.399	0.349	0.185
0.268	0.399	0.352	0.182
0.232	0.407	0.357	0.178
0.236	0.399	0.365	0.177
0.228	0.392	0.375	0.176
0.22	0.387	0.38	0.178
0.216	0.38	0.383	0.183
0.216	0.375	0.388	0.183
0.204	0.366	0.394	0.189
0.196	0.363	0.394	0.194
0.188	0.361	0.397	0.195
0.176	0.358	0.399	0.199
0.168	0.346	0.407	0.205
0.184	0.329	0.416	0.209
0.196	0.318	0.42	0.213
0.196	0.309	0.423	0.219
0.188	0.298	0.423	0.232
0.176	0.293	0.425	0.238
0.196	0.293	0.426	0.232
0.18	0.29	0.424	0.241
0.176	0.296	0.416	0.244
0.18	0.291	0.417	0.247
0.172	0.283	0.418	0.256
0.176	0.277	0.413	0.266
0.168	0.273	0.409	0.276
0.168	0.263	0.415	0.28
0.16	0.26	0.418	0.282
0.16	0.25	0.418	0.292
0.16	0.247	0.414	0.299
0.16	0.242	0.41	0.308
0.164	0.237	0.398	0.324
0.172	0.232	0.389	0.336
0.164	0.233	0.38	0.346
0.168	0.225	0.373	0.36
0.18	0.222	0.365	0.368
0.184	0.227	0.362	0.365
0.2	0.215	0.368	0.367
0.188	0.219	0.365	0.369
0.2	0.226	0.357	0.367
0.208	0.222	0.351	0.375
0.224	0.217	0.345	0.382
0.224	0.217	0.339	0.388
0.228	0.211	0.333	0.399
0.224	0.208	0.326	0.41
0.236	0.205	0.321	0.415
0.252	0.205	0.315	0.417
0.248	0.21	0.31	0.418
0.268	0.204	0.31	0.419
0.268	0.194	0.31	0.429
0.284	0.188	0.306	0.435
0.316	0.186	0.299	0.436
0.324	0.189	0.294	0.436
0.348	0.194	0.285	0.434
0.348	0.191	0.287	0.435
0.376	0.191	0.288	0.427
0.4	0.192	0.286	0.422
0.42	0.196	0.283	0.416
0.452	0.195	0.277	0.415
0.46	0.204	0.27	0.411
0.48	0.202	0.265	0.413
0.472	0.204	0.266	0.412
0.492	0.213	0.252	0.412
0.52	0.226	0.243	0.401
0.488	0.24	0.248	0.39
0.508	0.244	0.241	0.388
0.528	0.248	0.234	0.386
0.508	0.26	0.233	0.38
0.532	0.264	0.226	0.377
0.52	0.28	0.22	0.37
0.504	0.278	0.228	0.368
0.488	0.285	0.227	0.366
0.492	0.291	0.228	0.358
0.516	0.3	0.223	0.348
0.504	0.313	0.22	0.341
0.516	0.32	0.214	0.337
0.52	0.327	0.213	0.33
0.508	0.342	0.207	0.324
0.48	0.357	0.205	0.318
0.508	0.357	0.211	0.305
0.496	0.371	0.203	0.302
0.504	0.383	0.204	0.287
0.5	0.398	0.206	0.271
0.492	0.408	0.208	0.261
0.484	0.42	0.204	0.255
0.484	0.434	0.204	0.241
0.484	0.438	0.208	0.233
0.46	0.45	0.214	0.221
0.448	0.452	0.216	0.22
0.448	0.453	0.214	0.221
0.428	0.459	0.211	0.223
0.412	0.47	0.215	0.212
0.396	0.473	0.223	0.205
0.384	0.472	0.233	0.199
0.38	0.473	0.237	0.195
0.388	0.459	0.246	0.198
0.34	0.46	0.26	0.195
0.324	0.455	0.268	0.196
0.304	0.45	0.28	0.194
0.276	0.446	0.294	0.191
0.272	0.433	0.304	0.195
0.264	0.423	0.31	0.201
0.26	0.418	0.314	0.203
0.256	0.404	0.327	0.205
0.232	0.407	0.335	0.2
0.24	0.389	0.353	0.198
0.24	0.374	0.36	0.206
0.24	0.366	0.362	0.212
0.224	0.365	0.362	0.217
0.228	0.366	0.36	0.217
0.216	0.362	0.36	0.224
0.208	0.356	0.362	0.23
0.208	0.349	0.369	0.23
0.2	0.354	0.365	0.231
0.2	0.354	0.363	0.233
0.204	0.34	0.372	0.237
0.212	0.334	0.373	0.24
0.208	0.321	0.38	0.247
0.196	0.317	0.387	0.247
0.188	0.313	0.39	0.25
0.172	0.299	0.399	0.259
0.16	0.289	0.403	0.268
0.16	0.277	0.416	0.267
0.156	0.272	0.417	0.272
0.152	0.266	0.417	0.279
0.16	0.257	0.417	0.286
0.156	0.247	0.417	0.297
0.152	0.25	0.417	0.295
0.16	0.246	0.414	0.3
0.188	0.237	0.406	0.31
0.188	0.242	0.398	0.313
0.188	0.238	0.389	0.326
0.192	0.236	0.38	0.336
0.2	0.226	0.379	0.345
0.22	0.219	0.377	0.349
0.208	0.217	0.376	0.355
0.208	0.215	0.372	0.361
0.22	0.211	0.362	0.372
0.236	0.212	0.351	0.378
0.244	0.211	0.343	0.385
0.236	0.215	0.339	0.387
0.24	0.22	0.33	0.39
0.24	0.224	0.325	0.391
0.24	0.228	0.323	0.389
0.248	0.228	0.322	0.388
0.256	0.232	0.303	0.401
0.248	0.23	0.308	0.4
0.244	0.232	0.31	0.397
0.24	0.231	0.305	0.404
0.252	0.236	0.299	0.402
0.268	0.23	0.297	0.406
0.272	0.221	0.295	0.416
0.288	0.221	0.297	0.41
0.3	0.217	0.292	0.416
0.288	0.227	0.284	0.417
0.28	0.231	0.287	0.412
0.292	0.225	0.289	0.413
0.288	0.233	0.282	0.413
0.284	0.229	0.277	0.423
0.304	0.23	0.274	0.42
0.324	0.227	0.269	0.423
0.332	0.227	0.266	0.424
0.356	0.238	0.257	0.416
0.368	0.235	0.257	0.416
0.36	0.244	0.253	0.413
0.372	0.239	0.255	0.413
0.388	0.236	0.249	0.418
0.404	0.236	0.254	0.409
0.4	0.256	0.245	0.399
0.404	0.264	0.248	0.387
0.412	0.266	0.251	0.38
0.424	0.265	0.254	0.375
0.412	0.271	0.255	0.371
0.44	0.277	0.242	0.371
0.452	0.284	0.236	0.367
0.448	0.298	0.231	0.359
0.46	0.306	0.226	0.353
0.436	0.314	0.228	0.349
0.456	0.323	0.219	0.344
0.456	0.319	0.221	0.346
0.46	0.322	0.224	0.339
0.472	0.333	0.223	0.326
0.472	0.343	0.22	0.319
0.488	0.353	0.216	0.309
0.456	0.36	0.221	0.305
0.432	0.368	0.229	0.295
0.432	0.386	0.222	0.284
0.428	0.394	0.232	0.267
0.428	0.402	0.233	0.258
0.408	0.421	0.23	0.247
0.4	0.417	0.244	0.239
0.38	0.425	0.25	0.23
0.352	0.43	0.256	0.226
0.348	0.432	0.256	0.225
0.332	0.424	0.259	0.234
0.328	0.41	0.266	0.242
0.32	0.412	0.273	0.235
0.308	0.42	0.274	0.229
0.304	0.41	0.285	0.229
0.332	0.4	0.292	0.225
0.344	0.389	0.301	0.224
0.36	0.388	0.301	0.221
0.352	0.39	0.304	0.218
0.34	0.388	0.315	0.212
0.312	0.389	0.323	0.21
0.268	0.39	0.331	0.212
0.272	0.382	0.34	0.21
0.26	0.388	0.341	0.206
0.248	0.38	0.349	0.209
0.24	0.379	0.342	0.219
0.236	0.369	0.356	0.216
0.236	0.36	0.361	0.22
0.216	0.357	0.365	0.224
0.228	0.343	0.383	0.217
0.232	0.344	0.382	0.216
0.232	0.34	0.382	0.22
0.24	0.328	0.392	0.22
0.236	0.325	0.394	0.222
0.216	0.329	0.391	0.226
0.212	0.329	0.4	0.218
0.212	0.33	0.405	0.212
0.2	0.327	0.413	0.21
0.184	0.328	0.409	0.217
0.196	0.32	0.409	0.222
0.2	0.313	0.411	0.226
0.188	0.307	0.414	0.232
0.196	0.298	0.413	0.24
0.2	0.292	0.403	0.255
0.196	0.285	0.401	0.265
0.192	0.282	0.398	0.272
0.184	0.282	0.402	0.27
0.18	0.274	0.411	0.27
0.18	0.268	0.405	0.282
0.168	0.274	0.394	0.29
0.172	0.286	0.382	0.289
0.164	0.28	0.381	0.298
0.172	0.278	0.377	0.302
0.16	0.269	0.378	0.313
0.168	0.266	0.375	0.317
0.168	0.263	0.374	0.321
0.168	0.26	0.384	0.314
0.184	0.256	0.384	0.314
0.196	0.248	0.381	0.322
0.196	0.237	0.387	0.327
0.184	0.233	0.384	0.337
0.192	0.226	0.383	0.343
0.2	0.21	0.383	0.357
0.22	0.201	0.375	0.369
0.22	0.199	0.376	0.37
0.216	0.204	0.366	0.376
0.22	0.205	0.357	0.383
0.224	0.202	0.348	0.394
0.22	0.202	0.339	0.404
0.24	0.201	0.344	0.395
0.232	0.2	0.343	0.399
0.224	0.195	0.342	0.407
0.224	0.191	0.345	0.408
0.24	0.19	0.343	0.407
0.248	0.192	0.333	0.413
0.264	0.192	0.325	0.417
0.288	0.191	0.319	0.418
0.3	0.19	0.321	0.414
0.312	0.194	0.318	0.41
0.316	0.192	0.316	0.413
0.336	0.198	0.296	0.422
0.368	0.204	0.288	0.416
0.38	0.212	0.286	0.407
0.404	0.203	0.282	0.414
0.428	0.206	0.278	0.409
0.44	0.214	0.276	0.4
0.46	0.228	0.266	0.391
0.468	0.237	0.253	0.393
0.468	0.235	0.255	0.393
0.484	0.24	0.251	0.388
0.528	0.241	0.244	0.383
0.576	0.243	0.236	0.377
0.564	0.255	0.227	0.377
0.568	0.268	0.222	0.368
0.556	0.287	0.213	0.361
0.548	0.302	0.201	0.36
0.536	0.324	0.196	0.346
0.568	0.329	0.19	0.339
0.564	0.337	0.192	0.33
0.564	0.35	0.188	0.321
0.596	0.357	0.178	0.316
0.592	0.367	0.187	0.298
0.596	0.369	0.192	0.29
0.6	0.382	0.192	0.276
0.568	0.406	0.191	0.261
0.54	0.426	0.186	0.253
0.536	0.426	0.196	0.244
0.52	0.434	0.195	0.241
0.504	0.447	0.193	0.234
0.472	0.459	0.192	0.231
0.44	0.466	0.198	0.226
0.432	0.463	0.205	0.224
0.412	0.464	0.22	0.213
0.416	0.467	0.229	0.2
0.388	0.473	0.238	0.192
0.372	0.471	0.251	0.185
0.364	0.47	0.261	0.178
0.344	0.479	0.26	0.175
0.332	0.482	0.263	0.172
0.308	0.471	0.275	0.177
0.3	0.48	0.274	0.171
0.296	0.475	0.276	0.175
0.288	0.469	0.285	0.174
0.28	0.47	0.293	0.167
0.268	0.466	0.302	0.165
0.256	0.454	0.317	0.165
0.26	0.443	0.325	0.167
0.272	0.436	0.332	0.164
0.264	0.439	0.336	0.159
0.252	0.439	0.337	0.161
0.224	0.429	0.351	0.164
0.224	0.422	0.353	0.169
0.208	0.432	0.354	0.162
0.196	0.426	0.362	0.163
0.188	0.426	0.363	0.164
0.184	0.412	0.376	0.166
0.184	0.41	0.383	0.161
0.176	0.394	0.397	0.165
0.176	0.387	0.398	0.171
0.168	0.386	0.406	0.166
0.144	0.371	0.416	0.177
0.136	0.364	0.418	0.184
0.128	0.355	0.422	0.191
0.132	0.339	0.434	0.194
0.12	0.328	0.438	0.204
0.128	0.318	0.434	0.216
0.132	0.304	0.437	0.226
0.132	0.286	0.442	0.239
0.14	0.274	0.438	0.253
0.144	0.26	0.446	0.258
0.148	0.253	0.45	0.26
0.152	0.243	0.445	0.274
0.16	0.235	0.442	0.283
0.16	0.225	0.435	0.3
0.148	0.221	0.435	0.307
0.14	0.217	0.434	0.314
0.144	0.211	0.434	0.319
0.14	0.208	0.424	0.333
0.14	0.202	0.418	0.345
0.132	0.192	0.415	0.36
0.14	0.193	0.406	0.366
0.14	0.181	0.405	0.379
0.14	0.171	0.398	0.396
0.148	0.159	0.389	0.415
0.14	0.151	0.383	0.431
0.144	0.148	0.369	0.447
0.16	0.144	0.356	0.46
0.168	0.143	0.351	0.464
0.176	0.138	0.346	0.472
0.184	0.134	0.332	0.488
0.204	0.138	0.327	0.484
0.208	0.131	0.326	0.491
0.212	0.13	0.316	0.501
0.236	0.135	0.3	0.506
0.256	0.13	0.299	0.507
0.264	0.132	0.287	0.515
0.292	0.135	0.277	0.515
0.304	0.133	0.264	0.527
0.336	0.132	0.254	0.53
0.356	0.126	0.258	0.527
0.356	0.133	0.254	0.524
0.408	0.136	0.246	0.516
0.428	0.139	0.24	0.514
0.424	0.139	0.24	0.515
0.44	0.134	0.238	0.518
0.448	0.14	0.231	0.517
0.492	0.141	0.22	0.516
0.512	0.149	0.215	0.508
0.516	0.159	0.214	0.498
0.516	0.158	0.214	0.499
0.556	0.16	0.202	0.499
0.596	0.169	0.2	0.482
0.604	0.178	0.202	0.469
0.596	0.197	0.202	0.452
0.612	0.207	0.2	0.44
0.644	0.208	0.198	0.433
0.672	0.212	0.194	0.426
0.676	0.222	0.186	0.423
0.7	0.227	0.19	0.408
0.704	0.242	0.187	0.395
0.72	0.251	0.183	0.386
0.732	0.27	0.177	0.37
0.74	0.293	0.173	0.349
0.696	0.316	0.17	0.34
0.664	0.343	0.164	0.327
0.672	0.363	0.164	0.305
0.652	0.386	0.161	0.29
0.644	0.401	0.161	0.277
0.592	0.432	0.17	0.25
0.564	0.456	0.165	0.238
0.544	0.466	0.168	0.23
0.512	0.484	0.17	0.218
0.48	0.502	0.164	0.214
0.444	0.515	0.163	0.211
0.448	0.522	0.16	0.206
0.436	0.536	0.16	0.195
0.408	0.537	0.167	0.194
0.404	0.546	0.164	0.189
0.396	0.548	0.165	0.188
0.376	0.56	0.168	0.178
0.372	0.563	0.171	0.173
0.328	0.575	0.183	0.16
0.304	0.571	0.197	0.156
0.296	0.558	0.207	0.161
0.296	0.544	0.221	0.161
0.296	0.532	0.23	0.164
0.264	0.53	0.241	0.163
0.248	0.524	0.253	0.161
0.228	0.525	0.259	0.159
0.228	0.509	0.273	0.161
0.236	0.496	0.282	0.163
0.224	0.492	0.29	0.162
0.208	0.478	0.303	0.167
0.188	0.465	0.319	0.169
0.212	0.458	0.325	0.164
0.192	0.447	0.342	0.163
0.188	0.436	0.353	0.164
0.192	0.42	0.366	0.166
0.192	0.413	0.366	0.173
0.18	0.404	0.372	0.179
0.188	0.401	0.378	0.174
0.172	0.394	0.381	0.182
0.176	0.373	0.391	0.192
0.172	0.355	0.396	0.206
0.16	0.346	0.398	0.216
0.16	0.342	0.4	0.218
0.152	0.328	0.405	0.229
0.156	0.321	0.407	0.233
0.144	0.307	0.415	0.242
0.16	0.301	0.419	0.24
0.156	0.286	0.428	0.247
0.16	0.278	0.428	0.254
0.176	0.27	0.432	0.254
0.168	0.272	0.43	0.256
0.164	0.26	0.437	0.262
0.172	0.255	0.438	0.264
0.184	0.252	0.439	0.263
0.176	0.251	0.434	0.271
0.172	0.249	0.427	0.281
0.164	0.238	0.426	0.295
0.168	0.236	0.425	0.297
0.168	0.227	0.422	0.309
0.168	0.225	0.417	0.316
0.164	0.225	0.406	0.328
0.172	0.222	0.401	0.334
0.18	0.221	0.39	0.344
0.176	0.22	0.385	0.351
0.172	0.214	0.382	0.361
0.172	0.212	0.372	0.373
0.184	0.207	0.368	0.379
0.188	0.207	0.358	0.388
0.192	0.2	0.359	0.393
0.2	0.206	0.348	0.396
0.2	0.203	0.346	0.401
0.204	0.194	0.348	0.407
0.212	0.189	0.341	0.417
0.224	0.182	0.343	0.419
0.24	0.181	0.337	0.422
0.256	0.182	0.328	0.426
0.264	0.176	0.33	0.428
0.268	0.171	0.326	0.436
0.276	0.167	0.322	0.442
0.284	0.156	0.327	0.446
0.288	0.155	0.323	0.45
0.288	0.154	0.321	0.453
0.3	0.151	0.312	0.462
0.312	0.15	0.307	0.465
0.332	0.148	0.3	0.469
0.348	0.148	0.293	0.472
0.38	0.155	0.28	0.47
0.396	0.162	0.264	0.475
0.392	0.159	0.261	0.482
0.404	0.16	0.255	0.484
0.42	0.162	0.246	0.487
0.432	0.169	0.237	0.486
0.464	0.184	0.226	0.474
0.512	0.196	0.216	0.46
0.552	0.208	0.208	0.446
0.572	0.233	0.2	0.424
0.58	0.244	0.193	0.418
0.572	0.249	0.193	0.415
0.572	0.254	0.191	0.412
0.588	0.27	0.186	0.397
0.612	0.278	0.184	0.385
0.612	0.283	0.183	0.381
0.624	0.297	0.176	0.371
0.604	0.305	0.181	0.363
0.628	0.309	0.182	0.352
0.624	0.323	0.183	0.338
0.636	0.34	0.177	0.324
0.636	0.362	0.172	0.307
0.632	0.386	0.16	0.296
0.62	0.406	0.15	0.289
0.576	0.422	0.158	0.276
0.552	0.444	0.163	0.255
0.504	0.463	0.161	0.25
0.48	0.466	0.17	0.244
0.46	0.468	0.176	0.241
0.452	0.48	0.178	0.229
0.452	0.49	0.181	0.216
0.432	0.498	0.19	0.204
0.444	0.499	0.192	0.198
0.436	0.502	0.198	0.191
0.416	0.5	0.208	0.188
0.396	0.501	0.212	0.188
0.392	0.509	0.212	0.181
0.36	0.518	0.217	0.175
0.336	0.521	0.221	0.174
0.312	0.514	0.243	0.165
0.3	0.506	0.258	0.161
0.28	0.512	0.262	0.156
0.272	0.502	0.276	0.154
0.24	0.506	0.286	0.148
0.228	0.495	0.295	0.153
0.228	0.477	0.311	0.155
0.216	0.467	0.319	0.16
0.2	0.452	0.338	0.16
0.196	0.444	0.352	0.155
0.192	0.441	0.358	0.153
0.184	0.427	0.377	0.15
0.184	0.421	0.385	0.148
0.176	0.415	0.392	0.149
0.168	0.399	0.406	0.153
0.176	0.392	0.412	0.152
0.184	0.375	0.419	0.16
0.176	0.363	0.428	0.165
0.18	0.353	0.437	0.165
0.18	0.349	0.439	0.167
0.18	0.341	0.441	0.173
0.176	0.332	0.445	0.179
0.144	0.324	0.455	0.185
0.14	0.311	0.458	0.196
0.148	0.299	0.462	0.202
0.152	0.299	0.454	0.209
0.168	0.294	0.454	0.21
0.168	0.28	0.459	0.219
0.18	0.278	0.449	0.228
0.18	0.265	0.446	0.244
0.18	0.262	0.439	0.254
0.184	0.255	0.44	0.259
0.2	0.254	0.432	0.264
0.2	0.251	0.43	0.269
0.2	0.243	0.433	0.274
0.196	0.239	0.43	0.282
0.208	0.235	0.431	0.282
0.196	0.235	0.432	0.284
0.184	0.231	0.43	0.293
0.172	0.223	0.428	0.306
0.164	0.215	0.43	0.314
0.152	0.217	0.426	0.319
0.148	0.202	0.422	0.339
0.148	0.194	0.422	0.347
0.148	0.188	0.415	0.36
0.16	0.18	0.41	0.37
0.16	0.171	0.405	0.384
0.152	0.169	0.401	0.392
0.152	0.17	0.396	0.396
0.168	0.167	0.391	0.4
0.172	0.164	0.378	0.415
0.168	0.165	0.368	0.425
0.18	0.163	0.364	0.428
0.188	0.166	0.35	0.437
0.184	0.166	0.344	0.444
0.188	0.15	0.349	0.454
0.204	0.152	0.347	0.45
0.212	0.151	0.34	0.456
0.208	0.155	0.331	0.462
0.22	0.153	0.325	0.467
0.236	0.148	0.319	0.474
0.24	0.147	0.312	0.481
0.26	0.146	0.304	0.485
0.28	0.151	0.3	0.479
0.284	0.156	0.298	0.475
0.284	0.154	0.302	0.473
0.3	0.154	0.293	0.478
0.324	0.164	0.284	0.471
0.316	0.168	0.275	0.478
0.328	0.173	0.268	0.477
0.352	0.178	0.258	0.476
0.38	0.177	0.255	0.473
0.388	0.177	0.256	0.47
0.404	0.181	0.249	0.469
0.42	0.187	0.243	0.465
0.428	0.187	0.244	0.462
0.444	0.184	0.238	0.467
0.472	0.192	0.228	0.462
0.496	0.198	0.22	0.458
0.504	0.214	0.222	0.438
0.528	0.226	0.223	0.419
0.524	0.245	0.22	0.404
0.52	0.251	0.221	0.398
0.508	0.264	0.221	0.388
0.516	0.28	0.215	0.376
0.52	0.288	0.205	0.377
0.532	0.3	0.199	0.368
0.512	0.308	0.196	0.368
0.52	0.317	0.198	0.355
0.528	0.332	0.193	0.343
0.524	0.346	0.189	0.334
0.516	0.364	0.186	0.321
0.516	0.367	0.188	0.316
0.488	0.387	0.178	0.313
0.496	0.404	0.178	0.294
0.476	0.426	0.174	0.281
0.44	0.439	0.18	0.271
0.42	0.452	0.183	0.26
0.408	0.454	0.187	0.257
0.36	0.473	0.182	0.255
0.34	0.474	0.193	0.248
0.352	0.479	0.198	0.235
0.348	0.484	0.201	0.228
0.34	0.484	0.21	0.221
0.34	0.482	0.217	0.216
0.32	0.489	0.22	0.211
0.32	0.488	0.228	0.204
0.304	0.48	0.24	0.204
0.304	0.477	0.247	0.2
0.3	0.478	0.255	0.192
0.292	0.477	0.264	0.186
0.296	0.465	0.273	0.188
0.292	0.463	0.281	0.183
0.284	0.453	0.296	0.18
0.264	0.46	0.302	0.172
0.248	0.456	0.31	0.172
0.252	0.436	0.329	0.172
0.224	0.44	0.331	0.173
0.216	0.425	0.34	0.181
0.204	0.418	0.344	0.187
0.2	0.4	0.361	0.189
0.2	0.399	0.361	0.19
0.196	0.384	0.377	0.19
0.18	0.385	0.379	0.191
0.176	0.371	0.391	0.194
0.188	0.352	0.4	0.201
0.18	0.335	0.418	0.202
0.188	0.328	0.411	0.214
0.192	0.317	0.414	0.221
0.188	0.31	0.415	0.228
0.184	0.302	0.416	0.236
0.176	0.293	0.418	0.245
0.188	0.283	0.416	0.254
0.192	0.271	0.42	0.261
0.184	0.279	0.409	0.266
0.188	0.271	0.414	0.268
0.188	0.261	0.421	0.271
0.184	0.26	0.423	0.271
0.212	0.249	0.417	0.281
0.22	0.245	0.412	0.288
0.22	0.238	0.411	0.296
0.22	0.239	0.405	0.301
0.22	0.235	0.403	0.307
0.216	0.235	0.404	0.307
0.24	0.226	0.403	0.311
0.236	0.228	0.401	0.312
0.252	0.238	0.385	0.314
0.248	0.234	0.379	0.325
0.264	0.236	0.375	0.323
0.264	0.241	0.366	0.327
0.268	0.237	0.367	0.329
0.268	0.23	0.367	0.336
0.26	0.236	0.361	0.338
0.26	0.227	0.358	0.35
0.272	0.233	0.356	0.343
0.276	0.237	0.348	0.346
0.292	0.233	0.347	0.347
0.308	0.231	0.334	0.358
0.324	0.23	0.329	0.36
0.324	0.24	0.321	0.358
0.328	0.235	0.326	0.357
0.332	0.239	0.328	0.35
0.34	0.242	0.32	0.353
0.34	0.245	0.32	0.35
0.344	0.247	0.318	0.349
0.348	0.248	0.309	0.356
0.356	0.258	0.306	0.347
0.392	0.257	0.307	0.338
0.404	0.25	0.31	0.339
0.392	0.257	0.31	0.335
0.416	0.262	0.3	0.334
0.4	0.27	0.301	0.329
0.364	0.279	0.297	0.333
0.368	0.276	0.294	0.338
0.4	0.286	0.278	0.336
0.392	0.298	0.267	0.337
0.404	0.298	0.264	0.337
0.404	0.304	0.263	0.332
0.412	0.312	0.26	0.325
0.436	0.311	0.257	0.323
0.444	0.308	0.257	0.324
0.44	0.321	0.257	0.312
0.456	0.323	0.255	0.308
0.428	0.339	0.257	0.297
0.416	0.341	0.262	0.293
0.4	0.351	0.265	0.284
0.408	0.355	0.263	0.28
0.404	0.362	0.261	0.276
0.396	0.364	0.263	0.274
0.38	0.368	0.265	0.272
0.376	0.362	0.27	0.274
0.392	0.369	0.269	0.264
0.384	0.378	0.271	0.255
0.356	0.376	0.279	0.256
0.352	0.386	0.274	0.252
0.344	0.386	0.275	0.253
0.308	0.382	0.284	0.257
0.296	0.376	0.297	0.253
0.288	0.375	0.302	0.251
0.296	0.369	0.31	0.247
0.28	0.364	0.319	0.247
0.276	0.359	0.325	0.247
0.26	0.348	0.333	0.254
0.264	0.337	0.339	0.258
0.256	0.335	0.344	0.257
0.236	0.34	0.338	0.263
0.224	0.344	0.342	0.258
0.22	0.334	0.35	0.261
0.224	0.33	0.35	0.264
0.216	0.323	0.353	0.27
0.216	0.315	0.357	0.274
0.224	0.314	0.357	0.273
0.22	0.308	0.358	0.279
0.216	0.303	0.361	0.282
0.224	0.309	0.359	0.276
0.22	0.301	0.363	0.281
0.208	0.304	0.361	0.283
0.208	0.299	0.366	0.283
0.196	0.296	0.367	0.288
0.188	0.293	0.37	0.29
0.184	0.288	0.372	0.294
0.172	0.29	0.373	0.294
0.184	0.277	0.38	0.297
0.18	0.278	0.373	0.304
0.184	0.263	0.375	0.316
0.184	0.256	0.381	0.317
0.18	0.257	0.371	0.327
0.18	0.249	0.383	0.323
0.184	0.248	0.377	0.329
0.192	0.24	0.377	0.335
0.204	0.241	0.375	0.333
0.216	0.229	0.379	0.338
0.22	0.217	0.381	0.347
0.232	0.208	0.371	0.363
0.228	0.199	0.375	0.369
0.232	0.193	0.362	0.387
0.252	0.189	0.349	0.399
0.256	0.197	0.342	0.397
0.256	0.199	0.342	0.395
0.264	0.199	0.332	0.403
0.264	0.199	0.324	0.411
0.276	0.2	0.317	0.414
0.3	0.2	0.305	0.42
0.312	0.197	0.304	0.421
0.328	0.201	0.295	0.422
0.328	0.202	0.295	0.421
0.348	0.199	0.286	0.428
0.352	0.197	0.284	0.431
0.372	0.194	0.277	0.436
0.396	0.191	0.278	0.432
0.408	0.19	0.275	0.433
0.428	0.188	0.277	0.428
0.44	0.196	0.272	0.422
0.46	0.201	0.262	0.422
0.448	0.214	0.247	0.427
0.464	0.226	0.237	0.421
0.484	0.231	0.238	0.41
0.488	0.244	0.235	0.399
0.532	0.252	0.229	0.386
0.536	0.263	0.224	0.379
0.544	0.268	0.222	0.374
0.548	0.281	0.214	0.368
0.552	0.289	0.216	0.357
0.572	0.306	0.211	0.34
0.584	0.322	0.211	0.321
0.572	0.337	0.213	0.307
0.568	0.352	0.21	0.296
0.56	0.355	0.217	0.288
0.532	0.366	0.214	0.287
0.544	0.376	0.215	0.273
0.548	0.382	0.213	0.268
0.58	0.391	0.21	0.254
0.556	0.405	0.219	0.237
0.54	0.421	0.215	0.229
0.548	0.42	0.221	0.222
0.504	0.448	0.213	0.213
0.492	0.466	0.214	0.197
0.464	0.474	0.219	0.191
0.424	0.491	0.219	0.184
0.4	0.51	0.22	0.17
0.364	0.515	0.228	0.166
0.376	0.507	0.243	0.156
0.356	0.517	0.241	0.153
0.356	0.512	0.25	0.149
0.348	0.51	0.256	0.147
0.34	0.505	0.276	0.134
0.324	0.515	0.273	0.131
0.304	0.511	0.28	0.133
0.28	0.499	0.298	0.133
0.272	0.492	0.303	0.137
0.264	0.483	0.309	0.142
0.244	0.482	0.317	0.14
0.228	0.464	0.341	0.138
0.208	0.449	0.361	0.138
0.2	0.445	0.368	0.137
0.188	0.442	0.376	0.135
0.18	0.419	0.4	0.136
0.172	0.404	0.419	0.134
0.164	0.395	0.424	0.14
0.16	0.379	0.435	0.146
0.164	0.366	0.443	0.15
0.164	0.356	0.445	0.158
0.164	0.348	0.451	0.16
0.168	0.342	0.452	0.164
0.16	0.335	0.458	0.167
0.14	0.324	0.468	0.173
0.148	0.314	0.474	0.175
0.16	0.309	0.478	0.173
0.152	0.303	0.484	0.175
0.152	0.307	0.469	0.186
0.148	0.298	0.473	0.192
0.144	0.289	0.477	0.198
0.14	0.283	0.475	0.207
0.14	0.272	0.48	0.213
0.14	0.264	0.478	0.223
0.132	0.256	0.481	0.23
0.136	0.25	0.478	0.238
0.136	0.243	0.473	0.25
0.156	0.231	0.465	0.265
0.148	0.228	0.465	0.27
0.144	0.217	0.467	0.28
0.152	0.206	0.465	0.291
0.148	0.205	0.449	0.309
0.144	0.197	0.447	0.32
0.144	0.193	0.44	0.331
0.14	0.188	0.436	0.341
0.152	0.183	0.423	0.356
0.152	0.184	0.404	0.374
0.156	0.173	0.405	0.383
0.16	0.166	0.4	0.394
0.164	0.159	0.395	0.405
0.152	0.159	0.391	0.412
0.16	0.155	0.387	0.418
0.168	0.152	0.382	0.424
0.168	0.149	0.376	0.433
0.192	0.153	0.361	0.438
0.2	0.146	0.362	0.442
0.212	0.142	0.358	0.447
0.22	0.139	0.348	0.458
0.244	0.137	0.342	0.46
0.248	0.137	0.338	0.463
0.248	0.139	0.331	0.468
0.252	0.141	0.32	0.476
0.268	0.149	0.31	0.474
0.268	0.151	0.298	0.484
0.276	0.148	0.291	0.492
0.288	0.151	0.284	0.493
0.304	0.147	0.279	0.498
0.3	0.156	0.27	0.499
0.332	0.158	0.265	0.494
0.344	0.162	0.252	0.5
0.348	0.158	0.249	0.506
0.352	0.164	0.245	0.503
0.344	0.166	0.241	0.507
0.348	0.167	0.232	0.514
0.36	0.171	0.229	0.51
0.392	0.174	0.224	0.504
0.388	0.185	0.224	0.494
0.424	0.187	0.217	0.49
0.42	0.187	0.208	0.5
0.424	0.194	0.204	0.496
0.432	0.198	0.201	0.493
0.464	0.202	0.195	0.487
0.488	0.208	0.189	0.481
0.492	0.212	0.189	0.476
0.52	0.215	0.195	0.46
0.54	0.222	0.19	0.453
0.588	0.228	0.19	0.435
0.62	0.231	0.191	0.423
0.648	0.244	0.193	0.401
0.628	0.266	0.192	0.385
0.632	0.28	0.182	0.38
0.64	0.296	0.179	0.365
0.62	0.322	0.179	0.344
0.608	0.342	0.184	0.322
0.6	0.358	0.184	0.308
0.588	0.374	0.183	0.296
0.596	0.389	0.179	0.283
0.584	0.408	0.181	0.265
0.596	0.412	0.185	0.254
0.56	0.429	0.19	0.241
0.524	0.461	0.184	0.224
0.512	0.467	0.187	0.218
0.5	0.48	0.189	0.206
0.472	0.496	0.194	0.192
0.476	0.497	0.199	0.185
0.432	0.507	0.21	0.175
0.42	0.512	0.216	0.167
0.4	0.528	0.212	0.16
0.38	0.526	0.224	0.155
0.38	0.527	0.224	0.154
0.36	0.519	0.239	0.152
0.352	0.513	0.247	0.152
0.332	0.51	0.259	0.148
0.316	0.508	0.267	0.146
0.288	0.504	0.279	0.145
0.28	0.491	0.294	0.145
0.272	0.488	0.3	0.144
0.248	0.482	0.31	0.146
0.24	0.476	0.326	0.138
0.224	0.463	0.338	0.143
0.2	0.443	0.36	0.147
0.188	0.427	0.379	0.147
0.184	0.406	0.4	0.148
0.176	0.409	0.406	0.141
0.168	0.401	0.414	0.143
0.156	0.385	0.427	0.149
0.144	0.368	0.439	0.157
0.136	0.358	0.445	0.163
0.128	0.344	0.455	0.169
0.112	0.335	0.463	0.174
0.12	0.325	0.468	0.177
0.12	0.321	0.466	0.183
0.116	0.306	0.471	0.194
0.112	0.293	0.475	0.204
0.1	0.278	0.481	0.216
0.104	0.267	0.49	0.217
0.096	0.251	0.492	0.233
0.1	0.244	0.49	0.241
0.1	0.233	0.487	0.255
0.108	0.223	0.482	0.268
0.112	0.213	0.476	0.283
0.108	0.2	0.47	0.303
0.112	0.194	0.463	0.315
0.1	0.184	0.462	0.329
0.1	0.18	0.458	0.337
0.092	0.174	0.451	0.352
0.088	0.175	0.438	0.365
0.096	0.163	0.439	0.374
0.104	0.152	0.43	0.392
0.108	0.15	0.418	0.405
0.104	0.148	0.41	0.416
0.096	0.147	0.406	0.423
0.096	0.144	0.392	0.44
0.092	0.143	0.377	0.457
0.088	0.138	0.378	0.462
0.084	0.138	0.371	0.47
0.084	0.139	0.361	0.479
0.092	0.14	0.356	0.481
0.096	0.133	0.35	0.493
0.108	0.135	0.334	0.504
0.124	0.127	0.326	0.516
0.128	0.128	0.317	0.523
0.14	0.125	0.303	0.537
0.136	0.122	0.299	0.545
0.152	0.118	0.293	0.551
0.168	0.113	0.287	0.558
0.176	0.115	0.279	0.562
0.2	0.108	0.272	0.57
0.208	0.105	0.268	0.575
0.22	0.106	0.264	0.575
0.228	0.107	0.264	0.572
0.244	0.104	0.254	0.581
0.264	0.1	0.253	0.581
0.28	0.102	0.247	0.581
0.316	0.1	0.234	0.587
0.328	0.097	0.231	0.59
0.336	0.104	0.224	0.588
0.356	0.106	0.22	0.585
0.384	0.105	0.215	0.584
0.388	0.1	0.215	0.588
0.416	0.111	0.214	0.571
0.452	0.116	0.204	0.567
0.484	0.116	0.199	0.564
0.508	0.12	0.195	0.558
0.552	0.127	0.193	0.542
0.608	0.135	0.186	0.527
0.632	0.148	0.176	0.518
0.668	0.156	0.172	0.505
0.684	0.166	0.168	0.495
]
[h p]=ttest(a(:,1),a(:,2))
%-- 2020/12/14 16:18 --%
getData
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\YsimB.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
psd2 = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14]; %-v.u18 ]; %v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
%% 仿真abed值的测量
% load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price10-20201212.mat')
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\YsimA.csv';
%         abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\YsimB.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [v.u1 v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\YsimB.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1B = [v.u1 v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
ps_8col = [ps1A(:,1:4) ps1B(:,1:4)];
ps_8col = [ps1A(:,2:5) ps1B(:,2:5)];
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJa41.csv.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [v.u1 v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
psdA = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14]; %-v.u18 v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJa41.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [v.u1 v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
psdA = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14]; %-v.u18 v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJa41.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
psdA = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14]; %-v.u18 v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
mean(ps1A)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJa251.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
psdA = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14]; %-v.u18 v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
mean(ps1A)
[h p]=ttest(ps1A(:,1),ps1A(:,3) )
[h p]=ttest(ps1A(:,1),ps1A(:,2) )
[h p]=ttest(ps1A(:,3),ps1A(:,2) )
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ5_a251.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 ]; %v.u18 v.u22 v.u26 v.u30 v.u34 v.u38];
psdA = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18]; %-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
[[ 0 0 0 0 0.25 ]
[ 1 0 0 0 0 ]
[ 0 1 0 0 0 ]
[ 0 0 1 0 0 ]
[ 0 0 0 1 0 ]]
mean(ps1A)
Mix5_2000
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_1.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18  v.u22 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
[0.190476  0.190476  0.190476  0.190476  0.190476  0.047619
0.2251    0.1656    0.1646    0.1699    0.1640    0.1108]'
[0.111111  0.111111  0.111111  0.111111  0.111111  0.444444]'
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_2.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18  v.u22 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
[0.2231    0.1678    0.1671    0.1678    0.1608    0.1133;0.2251    0.1656    0.1646    0.1699    0.1640    0.1108]'
[0.190476  0.190476  0.190476  0.190476  0.190476  0.047619;0.2231    0.1678    0.1671    0.1678    0.1608    0.1133]'
[0.1245    0.1528    0.1726    0.1619    0.1688    0.2194]'
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_3.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18  v.u22 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
%-- 2020/12/14 22:01 --%
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18  v.u22 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
t=length(text(:,1));abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4';
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4';
[num,txt,raw] = xlsread(abedfilesname);
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18  v.u22 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18  v.u22 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18  v.u22 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14   ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
mean(ps1A(1:1000,:))
mean(ps1A(1:500,:))
mean(ps1A(1000:2100,:))
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14   ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
mean(ps1A(8000:end,:))
mean(ps1A(5000:end,:))
mean(ps1A(10000:end,:))
mean(ps1A(100:2000,:))
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14   ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
%-- 2020/12/15 18:51 --%
createTetrahedron4RPSD
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
[xState oS e123]=createTetrahedron4RPSD(6)
createTetrahedron4RPSD
[xState oS e123]=createTetrahedron4RPSD(6)
zlabel('Z')
ylabel('Y')
xlabel('X')
createTetrahedron4RPSD
[xState oS e123]=createTetrahedron4RPSD(6)
Axis_Exp_ZSJ
fileName = strcat('6A',DataFileIndex,'.csv');
%
Axis_Exp_ZSJ
DataFileIndex(tF)
strcat('6A',DataFileIndex(tF),'.csv')
expdata_fileName = strcat('6A',DataFileIndex(tF),'.csv');
Col4 = load(expdata_fileName);
cell2str(expdata_fileName)
cellstr(expdata_fileName)
Col4 = load(cellstr(expdata_fileName));
Axis_Exp_ZSJ
Col_simplex = Axis_Exp_ZSJ()
Axis_Exp_ZSJ
mean(AnM3_simplex(:,:,3))
Axis_Exp_ZSJ
%-- 2020/12/15 22:58 --%
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1A = [ v.u2 v.u6 v.u10  v.u14   ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
mean(ps1A)
psdB = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14]
mean(psdB)
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
psdB = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14]
mean(psdB)
%-- 2020/12/17 0:31 --%
1.2071/0.5
0.5/2.071
-0.3536 * sqrt(2)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price6\abed-1pop Strategy distributions (complete history)8b.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22 v.u26 v.u30];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30];
abed_mean=mean(psd)
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\Price8\abed-1pop Strategy distributions (complete history)8b.csv';
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22 v.u26 v.u30];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30];
abed_mean=mean(psd)
plot(psd(2:57,:),'DisplayName','psd(2:57,:)')
plot(psd(1:100,:),'DisplayName','psd(1:100,:)')
plot(psd,'DisplayName','psd')
freq_subspace_price8
AngleFunction2D3
freq_subspace_price8
plot(angleRetaRR)
freq_subspace_price8
mean(angleRetaRR)
mean(AngleMonentumRR)
ttest(angleRetaRR)
[h p]=ttest(angleRetaRR)
ttest(AngleMonentumRR)
[h p]=ttest(AngleMonentumRR)
freq_subspace_price8
m_angle = mean(angleRetaRR28col)
m_angle = mean(angleRetaRR28col)'
plot(m_angle)
plot(real(m_angle))
Omega_label=[]; for m=1:7;for n=m+1:8; Omega_label=[Omega_label;m*10+n]; end;end;
plot(1:28',real(m_angle),'ro')
plot(1:28',real(m_angle),'ro-')
text(1:28',real(m_angle),num2str(Omega_label))
m_angMo = mean(AngleMonentumRR28col)
figure;plot(1:28',real(m_angle),'bo-');
figure;plot(1:28',real(m_angMo),'bo-');
text(1:28',real(m_angMo),num2str(Omega_label))
grid on
title('Angular Momentum','FontSize',20)
close all
clear
figure;
scatter(real(m_angle),real(m_angMo),'rs');
text(real(m_angle),real(m_angMo),num2str(Omega_label));
grid on
title('AngMo vs AngTr','FontSize',20);
%-- 2020/12/17 11:10 --%
0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.2500 - 0.2500i  -0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 + 0.2500i   0.2500 - 0.2500i  -0.3536 - 0.0000i  -0.3536 + 0.0000i  -0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.3536 + 0.0000i   0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.0000 - 0.3536i   0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i   0.2500 + 0.2500i   0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 - 0.2500i  -0.2500 + 0.2500i  -0.3536 + 0.0000i  -0.3536 - 0.0000i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 - 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 - 0.0000i   0.3536 + 0.0000i   0.3536 - 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i
[   0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.2500 - 0.2500i  -0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 + 0.2500i   0.2500 - 0.2500i  -0.3536 - 0.0000i  -0.3536 + 0.0000i  -0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.3536 + 0.0000i   0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.0000 - 0.3536i   0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i   0.2500 + 0.2500i   0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 - 0.2500i  -0.2500 + 0.2500i  -0.3536 + 0.0000i  -0.3536 - 0.0000i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 - 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 - 0.0000i   0.3536 + 0.0000i   0.3536 - 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i]
[-7.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i  -0.5000 + 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.00]
[  -7.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i  -0.5000 + 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i ]
[  -7.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i  -0.5000 + 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 0.2071i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.0000i]
%-- 2020/12/17 14:44 --%
a=[4.60045E-05
3.36165E-05
9.54047E-06
4.59597E-06
-2.01164E-05
-3.44212E-05
-3.92458E-05
4.77402E-05
3.12672E-05
1.18073E-05
2.20473E-06
-1.64196E-05
-3.06305E-05
5.30745E-05
2.23431E-05
1.76658E-05
-1.13404E-06
-1.04876E-05
5.16958E-05
3.25767E-05
1.41392E-05
-4.47613E-06
4.77235E-05
2.83632E-05
1.42938E-05
4.29365E-05
3.70194E-05
3.34005E-05
]
histogram(a)
histogram(a,8)
histogram(a,20)
histogram(a,120)
barh(a)
spy(a)
bar(a)
pie(a)
area(a)
plot(a)
bar(a)
semilogx(a)
plot(a)
comet(a)
polar(a)
rose(a)
compass(a)
b=[  -7.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i  -0.5000 + 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 1.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 0.5000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.2071i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 - 0.2071i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.5000 + 0.0000i]
compass(a)
compass(a)
compass(b)
c=[[   0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.2500 - 0.2500i  -0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 + 0.2500i   0.2500 - 0.2500i  -0.3536 - 0.0000i  -0.3536 + 0.0000i  -0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.3536 + 0.0000i   0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.0000 - 0.3536i   0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i   0.2500 + 0.2500i   0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 - 0.2500i  -0.2500 + 0.2500i  -0.3536 + 0.0000i  -0.3536 - 0.0000i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 - 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 - 0.0000i   0.3536 + 0.0000i   0.3536 - 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i]]
compass(c)
for k=1:8; subplot(2,4,k);compass(c(:,k));end
for k=1:8; subplot(2,4,k);compass(c(:,k),'LineWidth',k);end
for k=1:8; subplot(2,4,k);compass(c(:,k),'LineWidth',1:8);end
for k=1:8; subplot(2,4,k);compass(c(:,k),'LineWidth',[1:8]');end
%-- 2020/12/17 23:44 --%
classification
c=k1*A(:,3) + k2*A(:,5) + k3*A(:,7)
d = unique(c)
[A(:,1) c]
c=k1*A(:,3) + k2*A(:,5) + k3*A(:,7)
c =
a*k1 + a*k3 + b*k2
b*k1 - b*k3
a*k1 + a*k3 - b*k2
0
b*k2 - a*k3 - a*k1
b*k3 - b*k1
- a*k1 - a*k3 - b*k2
a*k1 + a*k3 + b*k2
b*k1 - b*k3
a*k1 + a*k3 - b*k2
0
b*k2 - a*k3 - a*k1
b*k3 - b*k1
a*k1 + a*k3 + b*k2
b*k1 - b*k3
a*k1 + a*k3 - b*k2
0
b*k2 - a*k3 - a*k1
a*k1 + a*k3 + b*k2
b*k1 - b*k3
a*k1 + a*k3 - b*k2
0
a*k1 + a*k3 + b*k2
b*k1 - b*k3
a*k1 + a*k3 - b*k2
a*k1 + a*k3 + b*k2
b*k1 - b*k3
a*k1 + a*k3 + b*k2
A=[   12         0    0.2777   -0.2777    0.3928   -0.3928    0.2777   -0.2777         0
13         0    0.3928   -0.3928         0         0   -0.3927    0.3927         0
14         0    0.2777   -0.2777   -0.3928    0.3928    0.2777   -0.2777         0
15         0         0         0         0         0         0         0         0
16         0   -0.2777    0.2777    0.3928   -0.3928   -0.2777    0.2777         0
17         0   -0.3928    0.3928         0         0    0.3927   -0.3927         0
18         0   -0.2777    0.2777   -0.3928    0.3928   -0.2777    0.2777         0
23         0    0.2777   -0.2777    0.3928   -0.3928    0.2777   -0.2777         0
24         0    0.3927   -0.3927         0         0   -0.3928    0.3928         0
25         0    0.2777   -0.2777   -0.3928    0.3928    0.2777   -0.2777         0
26         0         0         0         0         0         0         0         0
27         0   -0.2777    0.2777    0.3928   -0.3928   -0.2777    0.2777         0
28         0   -0.3927    0.3927         0         0    0.3928   -0.3928         0
34         0    0.2777   -0.2777    0.3928   -0.3928    0.2777   -0.2777         0
35         0    0.3928   -0.3928         0         0   -0.3927    0.3927         0
36         0    0.2777   -0.2777   -0.3928    0.3928    0.2777   -0.2777         0
37         0         0         0         0         0         0         0         0
38         0   -0.2777    0.2777    0.3928   -0.3928   -0.2777    0.2777         0
45         0    0.2777   -0.2777    0.3928   -0.3928    0.2777   -0.2777         0
46         0    0.3927   -0.3927         0         0   -0.3928    0.3928         0
47         0    0.2777   -0.2777   -0.3928    0.3928    0.2777   -0.2777         0
48         0         0         0         0         0         0         0         0
56         0    0.2777   -0.2777    0.3928   -0.3928    0.2777   -0.2777         0
57         0    0.3928   -0.3928         0         0   -0.3927    0.3927         0
58         0    0.2777   -0.2777   -0.3928    0.3928    0.2777   -0.2777         0
67         0    0.2777   -0.2777    0.3928   -0.3928    0.2777   -0.2777         0
68         0    0.3927   -0.3927         0         0   -0.3928    0.3928         0
78         0    0.2777   -0.2777    0.3928   -0.3928    0.2777   -0.2777         0]
A(:,3).*A(:,5)
sum(A(:,3).*A(:,5))
sum(A(:,3).*A(:,7))
sum(A(:,4).*A(:,7))
C=[   0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.2500 - 0.2500i  -0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 + 0.2500i   0.2500 - 0.2500i  -0.3536 - 0.0000i  -0.3536 + 0.0000i  -0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.3536 + 0.0000i   0.3536 + 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i   0.2500 - 0.2500i   0.2500 + 0.2500i   0.3536 + 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 + 0.0000i
0.3536 + 0.0000i   0.0000 - 0.3536i   0.0000 + 0.3536i  -0.0000 - 0.3536i  -0.0000 + 0.3536i   0.2500 + 0.2500i   0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 - 0.2500i  -0.2500 + 0.2500i  -0.3536 + 0.0000i  -0.3536 - 0.0000i  -0.0000 - 0.3536i  -0.0000 + 0.3536i  -0.3536 + 0.0000i
0.3536 + 0.0000i  -0.3536 + 0.0000i  -0.3536 - 0.0000i   0.0000 + 0.3536i   0.0000 - 0.3536i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 + 0.0000i
0.3536 + 0.0000i  -0.2500 + 0.2500i  -0.2500 - 0.2500i   0.3536 - 0.0000i   0.3536 + 0.0000i   0.3536 - 0.0000i   0.3536 + 0.0000i  -0.3536 + 0.0000i]
sum(A(:,3).*A(:,5))
sum(A(:,3).*A(:,7))
sum(A(:,5).*A(:,7))
sum(A(:,6).*A(:,7))
clear
run('F:\CFH2014\dan2005\pictures.m')
load('F:\CFH2014\dan2005\snapshot.mat')
run('F:\CFH2014\dan2005\pictures2.m')
run('F:\CFH2014\dan2005\pictures_dan.m')
multpic
pictures
run('F:\CFH2014\dan2005\run199_3.m')
run('F:\CFH2014\dan2005\wrapper.m')
run('F:\CFH2014\dan2005\update_repl.m')
run('F:\CFH2014\dan2005\test.m')
run('F:\CFH2014\dan2005\pictures_dan.m')
run('F:\CFH2014\dan2005\edgew.m')
edgew
run2
rundata = wrapper(m_exp,table_entry,c,q,dyntype,a,step)
multpicsupp
run('F:\CFH2014\dan2005\markov_logit.m')
minimal
dyntype = 1
alpha = 0.10
tmax = 2000
matlab = 1
npic = 20
per = 36:85
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
disp(d)
min3 = av2
save -ascii minimal3 min3
dyntype = 1
alpha = 0.10
tmax = 2000
matlab = 1
npic = 20
per = 36:85
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
minimal
%-- 2020/12/18 20:01 --%
minimal
run('F:\CFH2014\dan2005\pic_dan.m')
run('F:\CFH2014\dan2005\pictures.m')
%-- 2020/12/18 22:05 --%
edgew
minimal
edgew
minimal
plot(p,f)
plot(p,f,p,av_f)
minimal
m_exp = [0.73 0.11 0.09 0.07; 0.24 0.59 0.13 0.04; 0.08 0.25 0.58 0.09; 0.03 0.11 0.14 0.72]
minimal
figure(1);
plot(p,f)
gset xrange [6:30]
gset yrange [0:0.175]
gset nokey
replot
figure(1)
plot(p,f)
axis([lp rp 0 0.16])
pause(0.001)
fine
pic_dan
pictures_dan
pic_dan(0.20,[10:step:18],dyntype)
snapshot
for i=1:3
imagesc(magic(i))
end
load('F:\CFH2014\dan2005\snapshot.mat')
pictures_dan
pic_dan(0.20,[10:step:18],dyntype)
load snapshot
pictures_dan
plot(ft,'DisplayName','ft')
ffft = ft'
plot(ffft,'DisplayName','ffft')
minimal
%-- 2020/12/19 13:52 --%
minimal
dyntype = 1
alpha = 0.10
tmax = 2000
matlab = 1
npic = 20
per = 36:450
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
minimal
dyntype = 1
alpha = 0.10
tmax = 2000
matlab = 1
npic = 20
per = 336:450
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
disp(d)
min1 = av5
dyntype = 1
alpha = 0.10
tmax = 2000
matlab = 1
npic = 20
per = 1:10 %36:45
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
disp(d)
min1 = av5
dyntype = 1
alpha = 0.010
tmax = 2000
matlab = 1
npic = 20
per = 1:10 %36:45
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
dyntype = 1
alpha = 0.010
tmax = 2000
matlab = 1
npic = 2
per = 1:10 %36:45
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
dyntype = 1
alpha = 0.010
tmax = 2000
matlab = 1
npic = 200
per = 1:10 %36:45
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
dyntype = 1
alpha = 0.010
tmax = 2000
matlab = 1
npic = 200
per = 1:100 %36:45
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
m_exp = [0.73 0.11 0.09 0.07; 0.24 0.59 0.13 0.04; 0.08 0.25 0.58 0.09; 0.03 0.11 0.14 0.72];
c = 20;
q = 1/3;
dyntype = 1
alpha = 0.010
tmax = 2000
matlab = 1
npic = 200
per = 1:100 %36:45
[d, av1, av2, av3, av4, av5, av6, av7, av8, av9, av0] = ...
edgew(alpha,m_exp,c,q,per,tmax,dyntype,matlab,npic);
WZJreadDan2005_20201219
ylim([0 0.4])
WZJreadDan2005_20201219
Wrr = [Wrr; count*ones(14,1) p',f',p',av_f']
%
WZJreadDan2005_20201219
save('Wrr14','Wrr')
WZJreadDan2005_20201219
%-- 2020/12/19 23:41 --%
read312data
size(num)
size(num,1)
col8=[];
for tI = 1:size(num,1)
tCol=zeros(1,8);
B= num(find(num(:,4))==tI,10);
for tJ = 1:6
tCol(B(tJ)+1) = tCol(B(tJ)+1) +1/6;
end
col8=[col8; tCol];
end
col8
read312data
B(tJ)
B =  num(find(num(:,4))==tI,11);
B =  num(find(num(:,4)==tI),11)
col8=[];
for tI = 1:size(num,1)
tCol=zeros(1,8);
B =  num(find(num(:,4)==tI),11);
for tJ = 1:6
tCol(B(tJ)+1) = tCol(B(tJ)+1) +1/6;
end
col8=[col8; tCol];
end
read312data
plot(col8(:,end))
read312data
plot(col8_01219W12(:,end))
plot(col8_01219W11(:,end))
e3 = [col8_01219W11;col8_01219W12]
mean(e3)
e3m=mean(e3)'
p=[];m=1:7;n=m+1:8; p= [p; e3m(m)*e3m(n)];end;end
p=[];m=1:7;n=m+1:8; p= [p; e3m(m)*e3m(n)] ;end;end
p=[];m=1:7;n=m+1:8; p= [p; e3m(m)*e3m(n)];end
p=[];for m=1:7; for n=m+1:8; p= [p; e3m(m)*e3m(n)] ;end;end
e3m=mean(e3(1:1000))'
e3m=mean(e3(1:1000,:))'
p=[];for m=1:7; for n=m+1:8; p= [p; e3m(m)*e3m(n)] ;end;end
e3m=mean(e3(1001:end,:))'
p=[];for m=1:7; for n=m+1:8; p= [p; e3m(m)*e3m(n)] ;end;end
C=[2.112333332	0.465888886
0.778999998	0.709555553
1.025333334	1.491222224
0.555555557	0.367444445
-0.333333331	-0.638888889
-1.444444444	-0.333333333
-2.69444444	-2.061888891
0.925083335	-1.027194447
1.436833335	0.131194446
1.125916666	0.332916667
-0.898027778	0.003333333
0.836694446	0.530166667
-1.31416667	0.49547222
0.981444445	-1.319194447
0.125916668	-0.028027778
0.935305557	0.34
1.031138887	0.115888889
-1.369722222	0.573694446
0.637055556	0.410222223
0.9905	0.212222224
1.32661111	0.076166667
0.489444446	-0.395388889
0.361111113	0.191111109
0.555555556	0.108722222
1.527777779	0.782722223
0.277777779	-0.111111113
0.777777777	0.218888889
2.583333333	0.386500001
]
rr(c)>
rr(C)
crr(C)
corr(C)
reg(ones(28,1),C(:,1),C(:,2))
regress(ones(28,1),C(:,1),C(:,2))
[b rb d s]=regress(ones(28,1),C(:,1),C(:,2))
X = [ones(size(y)) ,C(:,1)]; y = C(:,2)
y = C(:,2);X = [ones(size(y)) ,C(:,1)];
[b,bint,r,rint,stats] = regress(y,X)
corr(C(:,1),C(:,2),'style','spearman')
corr(C(:,1),C(:,2),'type' ,'Spearman')
[h p]=corr(C(:,1),C(:,2),'type' ,'Spearman')
[h p]=corrcoef(C(:,1),C(:,2))
help corrcoef
[r,lags] = xcorr(col8_01219W11(:,1:2) )
plot(r(:,1))
plot(r(:,end))
plot(r(:,3))
plot(r(:,2))
plot(r(:,1))
scatter(getcolumn(r(:,2:3),1),getcolumn(r(:,2:3),2))
scatter(getcolumn(r(:,3:end),1),getcolumn(r(:,3:end),2))
scatter(getcolumn(r(:,2:3),1),getcolumn(r(:,2:3),2))
semilogy(r(:,2:3),'DisplayName','r(:,2:3)')
scatter(getcolumn(r(:,[1,end]),1),getcolumn(r(:,[1,end]),2))
plot(r(:,[1,end]),'DisplayName','r(:,[1,end])')
bar(r(:,[1,end]),'DisplayName','r(:,[1,end])')
[h p]=corrcoef(C(:,1),C(:,2))
from_N_colExp_out_am
read312data
[ Yret11] = from_N_colExp_out_am(col8_01219W11,mean(col8_01219W11))
[ Yret12] = from_N_colExp_out_am(col8_01219W12,mean(col8_01219W12))
C = [ Yret11 Yret12];
y = C(:,2);X = [ones(size(y)) ,C(:,1)];
[b,bint,r,rint,stats] = regress(y,X)
[h p]=corrcoef(C(:,1),C(:,2))  %
[h p]=corr(C(:,1),C(:,2),'type' ,'Spearman')
C=[0.2777	0.0232	0.025617444
0.3928	0.014413	0.017840111
0.2777	0.019546	0.016032222
0	0.013543	0.007982
-0.2777	0.009541	0.004093333
-0.3928	0.025491	0.002933556
-0.2777	0.03799	0.088279556
0.2777	0.011044444	0.010910361
0.3927	0.014977778	0.009804722
0.2777	0.010377778	0.0048815
0	0.007311111	0.002503333
-0.2777	0.019533333	0.001794056
-0.3927	0.029111111	0.053988556
0.2777	0.009304944	0.006828056
0.3928	0.006447194	0.0033995
0.2777	0.004542028	0.001743333
0	0.012135083	0.001249389
-0.2777	0.018085278	0.037597889
0.2777	0.008743278	0.003055
0.3927	0.006159611	0.001566667
0.2777	0.016456833	0.001122778
0	0.024526111	0.033787778
0.2777	0.004267861	0.00078
0.3928	0.011402583	0.000559
0.2777	0.016993611	0.016822
0.2777	0.008033083	0.000286667
0.3927	0.011971944	0.008626667
0.2777	0.031985833	0.006182444
]
[h p]=corr(C(:,1),C(:,2),'type' ,'Spearman')
[h p]=corr(C(:,1),C(:,3),'type' ,'Spearman')
[h p]=corrcoef(C(:,1),C(:,2))
[h p]=corrcoef(C(:,1),C(:,3))
y = C(:,2);X = [ones(size(y)) ,C(:,1)];
[b,bint,r,rint,stats] = regress(y,X)
y = C(:,3);X = [ones(size(y)) ,C(:,1)];
[b,bint,r,rint,stats] = regress(y,X)
%-- 2020/12/20 12:38 --%
pic_dan
dyntype=1; step = 40;
pic_dan(0.20,[90:step:280],dyntype)
plot(ft,'DisplayName','ft')
pic_dan
dyntype=1; step = 40;
pic_dan(0.20,[90:step:280],dyntype)
v=ft'
plot(v,'DisplayName','v')
dyntype=1; step = 40;
pic_dan(0.20,[90:step:1280],dyntype)
v=ft'
plot(v,'DisplayName','v')
dyntype=1; step = 40;
pic_dan(0.20,[90:step:1280],dyntype)
v=ft'
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,8:16),'DisplayName','v(:,8:16)')
plot(v(:,17:23),'DisplayName','v(:,17:23)')
plot(v(:,24:end),'DisplayName','v(:,24:end)')
dyntype=2; step = 40;
pic_dan(0.20,[90:step:1280],dyntype)
v=ft'
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,8:16),'DisplayName','v(:,8:16)')
plot(v(:,17:23),'DisplayName','v(:,17:23)')
plot(v(:,24:end),'DisplayName','v(:,24:end)')
dyntype=2; step = 40;
pic_dan(1.20,[90:step:1280],dyntype)
v=ft'
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,8:16),'DisplayName','v(:,8:16)')
dyntype=2; step = 40;
pic_dan(0.0020,[90:step:1280],dyntype)
v=ft'
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,8:16),'DisplayName','v(:,8:16)')
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,17:23),'DisplayName','v(:,17:23)')
plot(v(:,24:end),'DisplayName','v(:,24:end)')
plot(v(:,17:23),'DisplayName','v(:,17:23)')
plot(v(:,8:16),'DisplayName','v(:,8:16)')
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,24:end),'DisplayName','v(:,24:end)')
dyntype=2; step = 400;
pic_dan(0.0020,[90:step:12800],dyntype)
v=ft'
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,8:16),'DisplayName','v(:,8:16)')
plot(v(:,1:7),'DisplayName','v(:,1:7)')
plot(v(:,17:23),'DisplayName','v(:,17:23)')
plot(v(:,1:end),'DisplayName','v(:,1:7)')
dyntype=2; step = 2000;
pic_dan(0.0020,[90:step:128000],dyntype)
v=ft'
plot(v(:,1:end),'DisplayName','v(:,1:7)')
plot(ft(:,1:3),'DisplayName','ft(:,1:3)')
plot(ft(:,20:end),'DisplayName','ft(:,20:end)')
plot(ft(:,13:end),'DisplayName','ft(:,13:end)')
plot(ft(:,9:end),'DisplayName','ft(:,9:end)')
plot(ft,'DisplayName','ft')
pic_dan_WZJ20201220
dyntype=1; step = 40; pic_dan(0.20,[90:step:1280],dyntype)
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[90:step:1280],dyntype)
23*0.0417
24*0.0417
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[90:step:1280],dyntype)
plot(av_f)
pic_dan_WZJ20201220
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[90:step:1280],dyntype)
%-- 2020/12/20 16:14 --%
read312data
from_N_colExp_out_am
read312data
mean(col8_01219W12)
mean(col8_01219W13)
[ Yret13] = from_N_colExp_out_am(col8_01219W13,mean(col8_01219W13))
C = [   Yret13];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28()
y = C(:,1);X = [ones(size(y)) ,r3col_TheoAuction1Yuan28(:,1)];
[b,bint,r,rint,stats] = regress(y,X)
r3col_TheoAuction1Yuan28(:,1)
y = C(:,1);X = [ones(size(y)) ,r3col_TheoAuction1Yuan28(:,1) ,r3col_TheoAuction1Yuan28(:,2) ,r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
read312data
mean(col8_01219W11)
[mean(col8_01219W11); mean(col8_01219W12); mean(ezp01220W13)]'
c = [mean(col8_01219W11); mean(col8_01219W12); mean(ezp01220W13)]'
corr(c)
dyntype=2; step = 2000;
pic_dan(0.0020,[90:step:128000],dyntype)
v=ft'
dyntype=2; step = 2000;
pic_dan(0.0020,[90:step:128000],dyntype)
F:\CFH2014\dan2005\pic_dan.m
pic_dan
dyntype=2; step = 2000;
pic_dan(0.0020,[90:step:128000],dyntype)
load snapshot
plot(ft,'DisplayName','ft')
plot(ft(:,1:5),'DisplayName','ft(:,1:5)')
plot(ft(:,13:18),'DisplayName','ft(:,13:18)')
plot(ft(:,1:4:24),'DisplayName','ft(:,13:18)')
plot(ft(:,2:6:24),'DisplayName','ft(:,13:18)')
v=ft'
plot(v(:,17:23),'DisplayName','v(:,17:23)')
plot(v(:,1:end),'DisplayName','v(:,1:7)')
plot(v(:,1:4:end),'DisplayName','v(:,1:7)')
plot(ft(:,8:12),'DisplayName','ft(:,8:12)')
plot(ft(:,20:end),'DisplayName','ft(:,20:end)')
plot(ft(:,1:7),'DisplayName','ft(:,1:7)')
plot(ft(:,10:14),'DisplayName','ft(:,10:14)')
plot(ft(:,18:end),'DisplayName','ft(:,18:end)')
plot(ft(:,9:15),'DisplayName','ft(:,9:15)')
plot(ft,'DisplayName','ft')
abedfilesname = 'F:\CFH2014\Auction1Yuan\data\ezp01220W14.csv';  %
col8_01220W14 = get312data(abedfilesname)
read312data
C = [ Yret14 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
abedfilesname = 'F:\CFH2014\Auction1Yuan\data\ezp01220W15.csv';  % ezp01220W15
col8_01220W15 = get312data(abedfilesname)
read312data
abedfilesname = 'F:\CFH2014\Auction1Yuan\data\ezp01220W15.csv';  % ezp01220W15
col8_01220W15 = get312data(abedfilesname)
[ Yret15] = from_N_colExp_out_am(col8_01220W14,mean(col8_01220W15))
[ Yret15] = from_N_colExp_out_am(col8_01220W15,mean(col8_01220W15))
C = [ Yret15 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
mean(col8_01220W15)
mean(col8_01220W14)
%-- 2020/12/20 22:46 --%
read312data
C = [ Yret15+Yret14 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [ Yret15+Yret14+Yret13 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [ Yret15+Yret13 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [ Yret15 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [  Yret14+Yret13 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [  Yret14+Yret15 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [  Yret11+Yret12+Yret13+Yret14+Yret15 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [  Yret11+Yret13+Yret14+Yret15 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
C = [Yret14+Yret15 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
scatter(C,r3col_TheoAuction1Yuan28)
scatter(C,r3col_TheoAuction1Yuan28(:,1))
scatter(r3col_TheoAuction1Yuan28(:,1),C)
scatter(r3col_TheoAuction1Yuan28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
C = [Yret14+Yret15+Yret13 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
scatter(r3col_TheoAuction1Yuan28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
[b,bint,r,rint,stats] = regress(y,X)
read312data
scatter(r3col_TheoAuction1Yuan28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
C = [Yret14+Yret15 ];
r3col_TheoAuction1Yuan28=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_TheoAuction1Yuan28(:,1) , ...
r3col_TheoAuction1Yuan28(:,2) , ...
r3col_TheoAuction1Yuan28(:,3)];
scatter(r3col_TheoAuction1Yuan28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
[b,bint,r,rint,stats] = regress(y,X)
find(abs(r3col_TheoAuction1Yuan28(:,1)-0.277)<-0.1)
find(abs(r3col_TheoAuction1Yuan28(:,1)-0.277)<0.1)
read312data
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
[h p]=corrcoef(Yret14, Yret15)  %
[h p]=corr(C(:,1),C(:,2), 'type', 'Spearman')
[h p]=corr(Yret14, Yret15, 'type', 'Spearman')
%-- 2020/12/21 21:17 --%
read312data
[ Yret17] = from_N_colExp_out_am(col8_01220W15,mean(col8_01220W17))
[ Yret17] = from_N_colExp_out_am(col8_01220W15,mean(col8_01221W17))
[ Yret17] = from_N_colExp_out_am(col8_01221W17,mean(col8_01221W17))
C = [Yret17 ]; %[Yret14+Yret15 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
C = [Yret17+Yret14+Yret15 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
mean(col8_01221W17)
mean(col8_01221W17)'
[ Yret16] = from_N_colExp_out_am(col8_01221W16,mean(col8_01221W16))
[ Yret17] = from_N_colExp_out_am(col8_01221W17,mean(col8_01221W17))
C = [Yret16]; %7+Yret14+Yret15 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
abedfilesname = 'F:\CFH2014\Auction1Yuan\data\ezp01221W16.csv';  %
col8_01221W16 = get312data(abedfilesname)
[ Yret16] = from_N_colExp_out_am(col8_01221W16,mean(col8_01221W16))
[ Yret17] = from_N_colExp_out_am(col8_01221W17,mean(col8_01221W17))
C = [Yret16]; %7+Yret14+Yret15 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
mean(col8_01221W16)'
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
C = [Yret16+Yret17+Yret14+Yret15 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
C = [Yret16+Yret17 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
read312data
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1)]% , ...
%                         r3col_Theo28(:,2) , ...
%                         r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
scatter(r3col_Theo28(:,1),C,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
C = [ Yret16+Yret17 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1)]  , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1)  , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
mean(col8_01221W17)
mean(col8_01221W16)
C = [ Yret14+Yret15+Yret16+Yret17 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
y = C(:,1);X = [ones(size(y)) ,...
r3col_Theo28(:,1)  , ...
r3col_Theo28(:,2) , ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
corr([ Yret14 Yret15 Yret16 Yret17 ])
corr([Yret13 Yret14 Yret15 Yret16 Yret17 ])
read312data
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
rra = [];
for tI = 1:size(C)
y = C(:,tI);X = [ones(size(y)) ,...
r3col_Theo28(:,1), ...
r3col_Theo28(:,2), ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X);
rra = [rra; ones(4,1)*tI b bint stats];
end
C = [ Yret14 Yret15 Yret16 Yret17 ];
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
rra = [];
for tI = 1:size(C)
y = C(:,tI);X = [ones(size(y)) ,...
r3col_Theo28(:,1), ...
r3col_Theo28(:,2), ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X);
rra = [rra; ones(4,1)*tI b bint stats];
end
rra = [];
for tI = 1:size(C)
y = C(:,tI);X = [ones(size(y)) ,...
r3col_Theo28(:,1), ...
r3col_Theo28(:,2), ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X);
rra = [rra; ones(4,1)*tI b bint stats'];
end
size(C)
rra = [];
for tI = 1:size(C,2)
y = C(:,tI);X = [ones(size(y)) ,...
r3col_Theo28(:,1), ...
r3col_Theo28(:,2), ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X);
rra = [rra; ones(4,1)*tI b bint stats'];
end
find(class_id == 1)
read312data
size(Omn)
class_id =zeros(size(Omn,1),1);
class_id(find(abs(tA(:,1)-0.3927)<0.1)) = 1
class_id(find(abs(tA(:,2)-0.3927)<0.1)) = 1
class_id(find(abs(tA(:,2)-0.3927)<0.1)) = 1
class_id(find(abs(tA(:,2)-0.2777)<0.1)) = 2
class_id(find(abs(tA(:,2) )<0.1)) = 3
class_id(find(abs(tA(:,2)+0.2777)<0.1)) = 4
class_id(find(abs(tA(:,2)+0.3927)<0.1)) = 5
find(class_id == 1)
read312data
D=C(find(class_id == 1),:)
D1=C(find(class_id == 1),:)
D2=C(find(class_id == 2),:)
ranksum(D1,D2)
ranksum(D1(:,1),D2(:,2))
[p h]=ranksum(D1(:,1),D2(:,2))
[p h]=ranksum(D1(:,1),D2(:,1))
[p h]=ranksum(D1(:,2),D2(:,2))
[p h]=ranksum(D1(:,3),D2(:,3))
[p h]=ranksum(D1(:,4),D2(:,4))
mean(D1)
mean(D2)
DD1=[D1(:,1);D1(:,2);D1(:,3);D1(:,4); ]
DD2=[D2(:,1);D2(:,2);D2(:,3);D2(:,4); ]
[p h]=ranksum(DD1,DD2)
[p h]=ttest(DD1)
[p h]=ttest(DD3)
[p h]=ttest(DD2)
%-- 2020/12/22 16:11 --%
read312data
scatter(r3col_Theo28(:,1),Yret19,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
read312data
CC = mean(C')
CC = mean(C')'
y=CC
X = [ones(size(y)) ,...
r3col_Theo28(:,1), ...
r3col_Theo28(:,2), ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X);
mean(col8_01222W18)
mean(col8_01222W19)
rint
scatter(r3col_Theo28(:,1),CC,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
scatter(r3col_Theo28(:,1),Yret15,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
figure(2);scatter(r3col_Theo28(:,1),Yret15,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
figure(18);scatter(r3col_Theo28(:,1),Yret18,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
figure(19);scatter(r3col_Theo28(:,1),Yret19,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
figure(13);scatter(r3col_Theo28(:,1),Yret13,129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
DDD = [r3col_Theo28(:,1),C]
[class_id DDD]
EE = [class_id DDD]
%-- 2020/12/22 18:18 --%
read312data
size(C,2)
D1=[]; for tK = 1:size(C,2)-1; D1=[D1; tD1(:,tK)];end
tD1=C(find(class_id == tI),:);
D1=[]; for tK = 1:size(C,2)-1; D1=[D1; tD1(:,tK)];end
tD2=C(find(class_id == tJ),:)
D2=[]; for tK = 1:size(C,2)-1; D2=[D2; tD2(:,tK)];end
[p h]=ranksum(D1(:,1),D2(:,1))
ranksumRR = [];
for tI = 1:size(C,2)-1
for tJ = tI+1:size(C,2)
tD1=C(find(class_id == tI),:);
D1=[]; for tK = 1:size(C,2)-1; D1=[D1; tD1(:,tK)];end
tD2=C(find(class_id == tJ),:)
D2=[]; for tK = 1:size(C,2)-1; D2=[D2; tD2(:,tK)];end
[p h]=ranksum(D1(:,1),D2(:,1))
ranksumRR = [ranksumRR; tI tJ h p size(D1) size(D2)]
end
end
ranksumRR = [];
for tI = 1:size(C,2)-1
for tJ = tI+1:size(C,2)
tD1=C(find(class_id == tI),:);
D1=[]; for tK = 1:size(C,2)-1; D1=[D1; tD1(:,tK)];end
tD2=C(find(class_id == tJ),:)
D2=[]; for tK = 1:size(C,2)-1; D2=[D2; tD2(:,tK)];end
[p h]=ranksum(D1,D2)
ranksumRR = [ranksumRR; tI tJ h p size(D1) size(D2)]
end
end
ranksumRR = [];
for tI = 1:5
for tJ = tI+1:5
tD1=C(find(class_id == tI),:);
D1=[]; for tK = 1:size(C,2)-1; D1=[D1; tD1(:,tK)];end
tD2=C(find(class_id == tJ),:)
D2=[]; for tK = 1:size(C,2)-1; D2=[D2; tD2(:,tK)];end
[p h]=ranksum(D1,D2)
ranksumRR = [ranksumRR; tI tJ h p size(D1) size(D2)]
end
end
ranksumRR = [];
for tI = 1:5
for tJ = tI+1:5
tD1=C(find(class_id == tI),:);
D1=[]; for tK = 1:size(C,2)-1; D1=[D1; tD1(:,tK)];end
tD2=C(find(class_id == tJ),:)
D2=[]; for tK = 1:size(C,2)-1; D2=[D2; tD2(:,tK)];end
[p h]=ranksum(D1,D2)
ranksumRR = [ranksumRR; tI tJ h p size(D1,1) size(D2,1)  mean(D1,1) mean(D2,1) ]
end
end
ranksumRR = [];
for tI = 1:5
for tJ = tI+1:5
tD1=C(find(class_id == tI),:);
D1=[]; for tK = 1:size(C,2); D1=[D1; tD1(:,tK)];end
tD2=C(find(class_id == tJ),:)
D2=[]; for tK = 1:size(C,2); D2=[D2; tD2(:,tK)];end
[p h]=ranksum(D1,D2)
ranksumRR = [ranksumRR; tI tJ h p size(D1,1) size(D2,1)  mean(D1,1) mean(D2,1) ]
end
end
read312data
ranksumRR = []; ranksumABS = [];
for tI = 1:5
for tJ = tI+1:5
tD1=C(find(class_id == tI),:);
D1=[]; for tK = 1:size(C,2); D1=[D1; tD1(:,tK)];end
tD2=C(find(class_id == tJ),:)
D2=[]; for tK = 1:size(C,2); D2=[D2; tD2(:,tK)];end
[p1 h1]=ranksum(D1,D2)
ranksumRR = [ranksumRR; tI tJ h1 p1 size(D1,1) size(D2,1)  mean(D1,1) mean(D2,1) ]
[p2 h2]=ranksum(abs(D1),abs(D2))
ranksumABS = [ranksumABS; tI tJ h2 p2 size(D1,1) size(D2,1)  mean(D1,1) mean(D2,1) ]
end
end
read312data
scatter(r3col_Theo28(:,1),C(:,8),129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
scatter(r3col_Theo28(:,1),C(:,7),129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),C,num2str(Omn));
scatter(r3col_Theo28(:,1),C(:,7),129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
text(r3col_Theo28(:,1),,C(:,7),num2str(Omn));
read312data
mean(col8_01222W21)
read312data
mean(col8_01222W21)
mean(col8_01222W20)
read312data
mean(col8_01222W21)
%-- 2020/12/23 0:24 --%
pic_dan_WZJ20201220
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[90:step:1280],dyntype)
load snapshot
plot(f)
plot(ft,'DisplayName','ft')
plot(ft','DisplayName','ft')
plot(ft(1:5:end)','DisplayName','ft')
plot((ft(1:5:end))','DisplayName','ft')
plot(ft','DisplayName','ft')
tft = ft'
plot(tft(:,1:25),'DisplayName','tft(:,1:25)')
pic_dan_WZJ20201220
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[90:step:1280],dyntype)
sum(f)
p
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[90:step:1280],dyntype)
plot(ft(:,1:5),'DisplayName','ft(:,1:5)')
plot(ft(1,:))
plot(ft(2,:))
plot(ft(3,:))
plot(ft(4,:))
plot(ft(5,:))
plot(ft(6,:))
plot(ft(7,:))
plot(ft(8,:))
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[10:step:1280],dyntype)
dyntype=1; step = 2; pic_dan_WZJ20201220(0.20,[10:step:1280],dyntype)
dyntype=1; step = 10; pic_dan_WZJ20201220(0.20,[10:step:10000],dyntype)
plot(ft(1:644,:),'DisplayName','ft(1:644,:)')
WZJ8eigencycle
plot(ft,'DisplayName','ft')
WZJ8eigencycle
plot(ft,'DisplayName','ft')
plot(ft','DisplayName','ft')
mean(ft)
plot(ans)
m_ft = mean(ft)
ft. - m_ft
d_ft = ft .-  m_ft
d_ft = ft - ones(1000,24).* m_ft
d_ft = ft - ones(1000,24)* m_ft.
d_ft = ft - ones(1000,1).* m_ft
d_ft = ft -   m_ft .* ones(1000,1)
tm = []; for I=1:1000; tm=[tm; m_ft];end
d_ft = ft -  tm;
plot(d_ft,'DisplayName','d_ft')
plot(d_ft','DisplayName','d_ft')
%-- 2020/12/23 10:04 --%
WZJ8eigencycle
[ Yret24] = from_N_colExp_out_am(ft,mean(ft))
plot(Yret24)
Omega_mn=[];
for m=1:split_price-1;for n=m+1:split_price;
Omega_mn = [Omega_mn; m n];
end;end
Omega_mn=[];
for m=1:split_price-1;for n=m+1:split_price;
Omega_mn = [Omega_mn; m n];
end;end
L_mn = [Omega_mn Yret24]
grid on
Omega_mn=[]; tK=1; L_matrix = zeros(split_price);
for m=1:split_price-1;for n=m+1:split_price;
Omega_mn = [Omega_mn; m n];
L_matrix(m,n) = Yret24(tK);
tK=tK+1;
end;end
L_mn = [Omega_mn Yret24]
Omega_mn=[]; tK=1; L_matrix = zeros(split_price);
for m=1:split_price-1;for n=m+1:split_price;
Omega_mn = [Omega_mn; m n];
L_matrix(m,n) = Yret24(tK);
L_matrix(n,m) = -Yret24(tK);
tK=tK+1;
end;end
L_mn = [Omega_mn Yret24]
x=1:24;
y=1:24;%注意这里有20个点
[xx,yy]=meshgrid(x,y);
figure;
hold on;
h=pcolor(xx,yy,L_matrix);
axis tight;
h=pcolor(L_matrix);
colormap(jet(8))
h=pcolor(L_matrix);
colormap(jet(8))
WZJ8eigencycle
bar off
colorbar off
colorbar on
%-- 2020/12/23 14:49 --%
WZJ8eigencycle
WZJreadDan2005_20201219
%-- 2020/12/23 22:52 --%
WZJ8eigencycle
r=[];for tI=1:size(ft,1);
acc=0; bcc=[];
for tJ = 1:30;
bcc=[bcc sum(ft(tI,1:tJ)) ];
end
m25 = min(find(bcc > 0.25)
m50 = min(find(bcc > 0.50)
m75 = min(find(bcc > 0.75)
r = [r;m75-m25 m50]
end
r=[];for tI=1:size(ft,1);
acc=0; bcc=[];
for tJ = 1:30;
bcc=[bcc sum(ft(tI,1:tJ)) ];
end
m25 = min(find(bcc > 0.25))
m50 = min(find(bcc > 0.50))
m75 = min(find(bcc > 0.75))
r = [r;m75-m25 m50];
end
plot(r,'DisplayName','r')
scatter(r(:,1),r(:,2))
r=[];for tI=1:size(ft,1);
acc=0; bcc=[];
for tJ = 1:30;
bcc=[bcc sum(ft(tI,1:tJ)) ];
end
m25 = min(find(bcc > 0.25))
m50 = min(find(bcc > 0.50))
m75 = min(find(bcc > 0.75))
r = [r;m75-m25 m50 m75 m25];
end
scatter(getcolumn(r(289:419,1:2),1),getcolumn(r(289:419,1:2),2))
scatter(getcolumn(r(289:419,1:2),1),getcolumn(r(289:419,1:2),2),'-')
scatter(r(289:419,1:2),'-')
plot(r(289:419,1:2),'-')
plot(r(289:end,1:2),'-')
line(r(289:end,1:2),'-')
line(r(289:end,1:2),'-.')
%-- 2020/12/25 17:12 --%
OneillYMQ20201118
A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]);
% A=eval([S.x1 S.x2 S.x3  S.y1 S.y2 S.y3 ]);
INs = length(A(:,1))-4
x1=A(INs,1); x2=A(INs,2);
x3=A(INs,3);  x4=A(length(A(:,1))-2,4);
y1=A(INs,5-1); y2=A(INs,6-1);
y3=A(INs,7-1);  y4=A(length(A(:,1))-2,8);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]);
% A=eval([S.x1 S.x2 S.x3  S.y1 S.y2 S.y3 ]);
INs = length(A(:,1))-4
x1=A(INs,1); x2=A(INs,2);
x3=A(INs,3);  x4=A(INs,4);
y1=A(INs,5); y2=A(INs,6);
y3=A(INs,7);  y4=A(INs,8);
% D_Eq_at_NE = eval(D_Eq_0)
% A6= D_Eq_0([1:3 5:7],[1:3 5:7])
D_Eq_at_NE = eval(D_Eq_0)
%          D_Eq_at_NE = eval(A6)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
%-- 2020/12/25 22:12 --%
read312data
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
rra = [];
for tI = 1:size(C,2)
y = C(:,tI);X = [ones(size(y)) ,...
r3col_Theo28(:,1)] , ...
r3col_Theo28(:,2), ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X);
rra = [rra; ones(2,1)*tI b bint stats'];
end
[r3col_Theo28 Omn class_id]=TheoAuction1Yuan28();
rra = [];
for tI = 1:size(C,2)
y = C(:,tI);X = [ones(size(y)) ,...
r3col_Theo28(:,1) , ...
r3col_Theo28(:,2), ...
r3col_Theo28(:,3)];
[b,bint,r,rint,stats] = regress(y,X);
rra = [rra; ones(4,1)*tI b bint stats'];
end
DD=mean(C')'
DD=mean(C(7:10)')'
DD=mean(C(:,7:10)')'
corr(C)
DE = corr(C)
%-- 2020/12/26 7:30 --%
tmp20201226
num = sortrows(num,[1 3])
expid = unique(num(:,1))
tmp20201226
sess = num(find(num(:,1)==expid(eid)),:)
sess = num(find(num(:,1)==expid(eid)),:)
sess = sortrows(sess,[1 3])
for tI = 1:size(sess,1)/6
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:6
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/6;
end
col8=[col8; tCol];
end
sess = num(find(num(:,1)==expid(eid)),:)
sess = sortrows(sess,[1 3])
for tI = 1:size(sess,1)/6
tCol=zeros(1,5);
B =  sess(find(sess(:,3)==tI),6);
for tJ = 1:6
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/5;
end
col8=[col8; tCol];
end
for eid = 1:4
sess = num(find(num(:,1)==expid(eid)),:)
sess = sortrows(sess,[1 3])
for tI = 1:size(sess,1)/6
tCol=zeros(1,5);
B =  sess(find(sess(:,3)==tI),6);
for tJ = 1:6
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/5;
end
col8=[col8; tCol];
end
end
[ Yret3 ] = from_N_colExp_out_am(col8,mean(col8))
r3col_Theo10 = [0.3693	0.5976
-0.5976	0.3693
-0.3693	-0.5976
0.5976	-0.3693
0.3693	0.5976
0.5976	-0.3693
-0.5976	0.3693
-0.5976	0.3693
0.3693	0.5976
-0.3693	-0.5976]
y = Yret3;X = [ones(size(y)) ,...
r3col_Theo10(:,1) , ...
r3col_Theo10(:,2) ];
[b,bint,r,rint,stats] = regress(y,X)
y = Yret3;X = [ones(size(y)) ,...
r3col_Theo10(:,1)  ];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
A=[0.3693	0.2
-0.5976	-1.1492
-0.3693	-1
0.5976	1.9492
0.3693	0.4882
0.5976	-0.28
-0.5976	-0.0082
-0.5976	-2.6591
0.3693	1.9981
-0.3693	-3.9391
]
corr(A)
corr(A,'type','spearman')
[r h]=corr(A,'type','spearman')
tmp20201226
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
psdA = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18]; %-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
tmp20201226
ps1A = [ v.u2 v.u6 v.u10  v.u14 v.u18 ]; %v.u22 v.u26 v.u30 v.u34 v.u38];
psdA = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18]; %-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
col8 = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18]; %-v.u22 v.u22-v.u26 v.u26-v.u30 v.u30-v.u34 v.u34-v.u38 v.u38];
[ Yret3 ] = from_N_colExp_out_am(col8,mean(col8))
r3col_Theo10 = [0.3693	0.5976
-0.5976	0.3693
-0.3693	-0.5976
0.5976	-0.3693
0.3693	0.5976
0.5976	-0.3693
-0.5976	0.3693
-0.5976	0.3693
0.3693	0.5976
-0.3693	-0.5976]
y = Yret3;X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
exp10=[0.2
-1.1492
-1
1.9492
0.4882
-0.28
-0.0082
-2.6591
1.9981
-3.9391
]
mean(col8)
tmp20201226
%-- 2020/12/26 14:00 --%
tmp20201226
read312data
tmp20201226
size(sess,1)
tmp20201226
expid = unique(num(:,7))
for eid = 1:3
sess = num(find(num(:,7)==expid(eid)),:)
sess = sortrows(sess,[4 8])
for tI = 1:800
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:6
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/5;
end
col8=[col8; tCol];
end
end
col8=[];
expid = unique(num(:,7))
for eid = 1:3
sess = num(find(num(:,7)==expid(eid)),:)
sess = sortrows(sess,[4 8])
for tI = 1:800/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:6
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/5;
end
col8=[col8; tCol];
end
end
col8=[];
expid = unique(num(:,7))
for eid = 1:3
sess = num(find(num(:,7)==expid(eid)),:)
sess = sortrows(sess,[4 8])
for tI = 1:800/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/5;
end
col8=[col8; tCol];
end
end
[ Yret3 ] = from_N_colExp_out_am(col8,mean(col8))
r3col_Theo10 = [0.3693	0.5976
-0.5976	0.3693
-0.3693	-0.5976
0.5976	-0.3693
0.3693	0.5976
0.5976	-0.3693
-0.5976	0.3693
-0.5976	0.3693
0.3693	0.5976
-0.3693	-0.5976]
y = Yret3;X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
[ Yret3 ] = from_N_colExp_out_am(col8,mean(col8))
r3col_Theo10 = [0.3693	0.5976
-0.5976	0.3693
-0.3693	-0.5976
0.5976	-0.3693
0.3693	0.5976
0.5976	-0.3693
-0.5976	0.3693
-0.5976	0.3693
0.3693	0.5976
-0.3693	-0.5976]
y = Yret3;X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
[ Yret3 ] = from_N_colExp_out_am(col8,mean(col8))
r3col_Theo10 = [0.3693	0.5976
-0.5976	0.3693
-0.3693	-0.5976
0.5976	-0.3693
0.3693	0.5976
0.5976	-0.3693
-0.5976	0.3693
-0.5976	0.3693
0.3693	0.5976
-0.3693	-0.5976]
y = Yret3;X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
r3col_Theo10 = [0.3693	0.5976
-0.5976	0.3693
-0.3693	-0.5976
0.5976	-0.3693
0.3693	0.5976
0.5976	-0.3693
-0.5976	0.3693
-0.5976	0.3693
0.3693	0.5976
-0.3693	-0.5976]
y = Yret3;X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
mean(col8)
tmp20201226
Yret3
tmp20201226
mean(col8)
tmp20201226
mean(col8)
A=[16.90283333	1.622277778	4.043611111	5.665888889
0.328	-2.083333333	-3.333611111	-5.416944444
-16.97033334	-1.694444444	-0.145555556	-1.84
-0.2605	2.1555	-0.564444444	1.591055556
32.88583333	0.722944444	5.813611111	6.536555556
6.517	1.437388889	0.223611111	1.661
-22.5	-0.538055556	-1.993611111	-2.531666667
9.141	-1.833333333	0.128611111	-1.704722222
24.07283334	0.472944444	2.351388889	2.824333333
-1.312333333	-2.090388889	0.206666667	-1.883722222
]
corr(A)
tmp20201226
clear
load snapshot
corr(Yret_pop,Yret_fix)
corr(Yret3_pop,Yret3_fix)
y = Yret3_fix;X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
y = Yret3_pop;
X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
for eid = 1:2
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
tmp20201226
for eid = 1:3
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
y = Yret3_fix;
X = [ones(size(y)) ,...
r3col_Theo10(:,1),...
r3col_Theo10(:,2)];
[b,bint,r,rint,stats] = regress(y,X)
tmp20201226
col8_fix=[];
expid = unique(num(:,7))
for eid = 1:3
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
clear
tmp20201226
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
mean(col8_fix)
col8_fix=[];
expid = unique(num(:,7))
for eid = 2:2
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
col8_fix=[];
expid = unique(num(:,7))
for eid = 3:3
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
tmp20201226
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
mean(col8_fix)
col8_fix=[];
expid = unique(num(:,7))
for eid = 1:1
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
col8_fix=[];
expid = unique(num(:,7))
for eid = 3:3
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
clear
tmp20201226
col8_fix=[];
expid = unique(num(:,7))
for eid = 1:3
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:2
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/2;
end
col8_fix=[col8_fix; tCol];
end
end
[ Yret3_fix ] = from_N_colExp_out_am(col8_fix,mean(col8_fix))
tmp20201226
col8_fix=[]; groupsize = 2;    %%***********************
expid = unique(num(:,7))
for eid = 1:size(expid)
sess = num(find(num(:,7)==expid(eid) & num(:,4) < 501),:)
sess = sortrows(sess,[4 8])
for tI = 1:1000/2
tCol=zeros(1,5);
B =  sess(find(sess(:,4)==tI),11);
for tJ = 1:groupsize
tCol(B(tJ)+1) = tCol(B(tJ)+1) + 1/groupsize;
end
col8_fix=[col8_fix; tCol];
end
end
tmp20201226
clear
tmp20201226
mean(col8_fix)
tmp20201226
%-- 2020/12/27 0:40 --%
test1226result
[b,bint,r,rint,stats] = regress(y,X)
rr=[];
y = a(:,5)
X = [ones(size(y)), a(:,2),  a(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
rr = [rr; 5 b' stats]
rr=[];
for sessionid = 5:18
y = a(:,sessionid)
X = [ones(size(y)), a(:,2),  a(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
rr = [rr; sessionid b' stats]
end
sum(a(:,16:18)')'
b = sum(a(:,16:18)')'
y=b;
X = [ones(size(y)), a(:,2),  a(:,3)];
[b,bint,r,rint,stats] = regress(y,X)
rr = [rr; sessionid b' stats]
scatter(a(:,2),b,'S')
bb = sum(a(:,16:18)')'
scatter(a(:,2),bb,'S')
scatter(a(:,2)+1/3*a(:,3),bb,'S')
scatter(a(:,2),bb,'S')
scatter(a(:,2)+1/3*a(:,3),bb,'S')
scatter(a(:,2)+1/3.5*a(:,3),bb,'S')
scatter(a(:,2)+a(:,3),bb,'S')
co=corr(a(:,5:18))
mesh(co)
surf(co)
pie(co)
scatter(getcolumn(rr(:,3:4),1),getcolumn(rr(:,3:4),2))
r1=[];
for sessionid = 5:18
y = a(:,sessionid)
X = [ones(size(y)), a(:,2) ];
[b,bint,r,rint,stats] = regress(y,X)
r1 = [r1; sessionid b' stats]
end
%-- 2020/12/27 12:54 --%
N_cycle0315_Price6
N_cycle0315
figure(1)
plot(real(t2(2:end,[1 2 3 4])),'linewidth',4)
plot3(t2(:,[1]),t2(:,[2]),t2(:,[3]));hold on
plot3(1/4,1/4,1/4,'r*')
plot3(t2(1:100:1000,[1]),t2(1:100:1000,[2]),t2(1:100:1000,[3]),'r.-')
text(t2(100,[1]),t2(100,[2]),t2(100,[3]),'100')
xlabel('x');ylabel('y');zlabel('z');
ZSJ1227TimeSeries
N_cycle0315
ZSJ1227TimeSeries
t2=[[[x1 x2 x3 x4] + 0.1*c2] 0]
for k=0.01:0.02:50
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*0.01; %时间步长取0.01
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*0.01;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*0.01;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*0.01;
t2=[t2; x1 x2 x3 x4 k];
end
figure(1)
plot(real(t2(:,[1 2 3 4])),'linewidth',4)
ZSJ1227TimeSeries
t2=[[[x1 x2 x3 x4] + 0.1*c2] 0]
tstep = 0.02; %时间步长取0.02
for k=0:tstep:4*pi
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
t2=[t2; x1 x2 x3 x4 k];
end
figure(1)
plot(real(t2(:,[1 2 3 4])),'linewidth',4)
ZSJ1227TimeSeries
figure(2)
plot3(t2(:,[1]),t2(:,[2]),t2(:,[3]));hold on
plot3(1/4,1/4,1/4,'r*')
plot3(t2(1:100:1000,[1]),t2(1:100:1000,[2]),t2(1:100:1000,[3]),'r.-')
text(t2(100,[1]),t2(100,[2]),t2(100,[3]),'100')
xlabel('x');ylabel('y');zlabel('z');
ZSJ1227TimeSeries
figure(3)
e1=[1 0 0];
e2=[1/2 sqrt(3/4) 0];
phi3=(pi/2-acos(...
((sqrt(3)/2)^2 *2 - 1)/ ((sqrt(3)/2)^2 *2) ...
));
e3=1/3*(e2+e1) + [0 0 cos(phi3)*sqrt(3/4)];
t3=[];
for m=1:size(t2,1)
tmp = real(t2(m,1:3))* [e1; e2; e3]' ;
t3=[t3; tmp 1-sum(tmp) t2(m,5)];
end
plot3(t3(:,[1]),t3(:,[2]),t3(:,[3]));hold on
Nsah = [1/4,1/4,1/4]*[e1; e2; e3]';
plot3(Nsah(1),Nsah(2),Nsah(3),'r*')
plot3(t3(1:100:1000,[1]),t3(1:100:1000,[2]),t3(1:100:1000,[3]),'r.')
text(t3(100,[1]),t3(100,[2]),t3(100,[3]),'100')
xlabel('x');ylabel('y');zlabel('z');
figure(3)
e1=[1 0 0];
e2=[1/2 sqrt(3/4) 0];
phi3=(pi/2-acos(...
((sqrt(3)/2)^2 *2 - 1)/ ((sqrt(3)/2)^2 *2) ...
));
e3=1/3*(e2+e1) + [0 0 cos(phi3)*sqrt(3/4)];
t3=[];
sampleSize = size(t2,1)
for m=1:size(t2,1)
tmp = real(t2(m,1:3))* [e1; e2; e3]' ;
t3=[t3; tmp 1-sum(tmp) t2(m,5)];
end
plot3(t3(:,[1]),t3(:,[2]),t3(:,[3]));hold on
Nsah = [1/4,1/4,1/4]*[e1; e2; e3]';
plot3(Nsah(1),Nsah(2),Nsah(3),'r*')
plot3(t3(1:2:sampleSize,[1]),t3(1:2:sampleSize,[2]),t3(1:2:sampleSize,[3]),'r.')
%     text(t3(100,[1]),t3(100,[2]),t3(100,[3]),'100')
xlabel('x');ylabel('y');zlabel('z');
xlim([0,1]); ylim([0,1]); zlim([0,1]);
ZSJ1227TimeSeries
NsahT = Ne(1:3)*[e1; e2; e3]';
plot3(NsahT(1),NsahT(2),NsahT(3),'r*')
plot3(t3(1:2:sampleSize,[1]),t3(1:2:sampleSize,[2]),t3(1:2:sampleSize,[3]),'r.')
%     text(t3(100,[1]),t3(100,[2]),t3(100,[3]),'100')
xlabel('x');ylabel('y');zlabel('z');
xlim([0,1]); ylim([0,1]); zlim([0,1]);
ZSJ1227TimeSeries
figure(2)
t3=[];
sampleSize = size(t2,1)
for m=1:size(t2,1)
tmp = real(t2(m,1:3));
t3=[t3; tmp 1-sum(tmp) t2(m,5)];
end
plot3(t3(:,[1]),t3(:,[2]),t3(:,[3]));hold on
NsahT = Ne(1:3) ;
plot3(NsahT(1),NsahT(2),NsahT(3),'r*')
plot3(t3(1:2:sampleSize,[1]),t3(1:2:sampleSize,[2]),t3(1:2:sampleSize,[3]),'r.')
%     text(t3(100,[1]),t3(100,[2]),t3(100,[3]),'100')
xlabel('x');ylabel('y');zlabel('z');
xlim([0,1]); ylim([0,1]); zlim([0,1]);
ZSJ1227TimeSeries
grid on
plot3([0 1],[0 1],[0 1] )
ZSJ1227TimeSeries
%-- 2020/12/27 22:06 --%
ZSJ1227TimeSeries
text(t3(1:10:end,[1]),t3(1:10:end,[2]),t3(1:10:end,[3]),'100')
text(t3(1:10:end,[1]),t3(1:10:end,[2]),t3(1:10:end,[3]),num2str(1:10:end))
plot3(NsahT(1),NsahT(2),NsahT(3),'r*')
plot3(t3(1:2:sampleSize,[1]),t3(1:2:sampleSize,[2]),t3(1:2:sampleSize,[3]),'r.')
text(t3(1:10:end,[1]),t3(1:10:end,[2]),t3(1:10:end,[3]),num2str(1:10:end))
xlabel('x');ylabel('y');zlabel('z');
xlim([0,1]); ylim([0,1]); zlim([0,1]);
ZSJ1227TimeSeries
text(t3(1:10:end,1),t3(1:10:end,2),t3(1:10:end,3),num2str(1:10:end))
text(t3(1:10:end,1),t3(1:10:end,2),t3(1:10:end,3),num2str(1:10:170))
text(t3(1:10:end,1),t3(1:10:end,2),t3(1:10:end,3),num2str([1:10:170]))
ZSJ1227TimeSeries
view(0,0)
view(0,90)
view(90,90)
view(90,0)
xlabel('x');ylabel('y');zlabel('z','fontsize',20);
xlabel('x','fontsize',30);ylabel('y','fontsize',30);zlabel('z','fontsize',30);
xlim([0,.5]); ylim([0,.5]); zlim([0,.5]);
xlim([0,.5]); ylim([0,.5]); zlim([0,.5]); grid on;box on
ZSJ1227TimeSeries
legend on
legend([1:4])
K>> legend([1:4],'location','NorthEast')
legend([1:4],'location','NorthEast')
legend(1:4,'location','NorthEast')
legend('X1','X2','X3','X4','location','NorthEast')
grid on
ZSJ1227TimeSeries
legend('X1','X2','X3','X4','location','NorthEast')
mean(t2)
mean(real(t2))
ShujieLandscape1228
eval(V_F)
X = t2(length(t2(:,1)),:) + eval(V_F)*tstep;
x1 = t2(length(t2(:,1)),1); x2 = t2(length(t2(:,1)),2); x3 = t2(length(t2(:,1)),3); x4 = t2(length(t2(:,1)),4);
X = t2(length(t2(:,1)),:) + eval(V_F)*tstep;
t2(length(t2(:,1)),:)
x1 = t2(length(t2(:,1)),1); x2 = t2(length(t2(:,1)),2); x3 = t2(length(t2(:,1)),3); x4 = t2(length(t2(:,1)),4);
X = t2(length(t2(:,1)),:)' + eval(V_F)*tstep;
eval(V_F)*tstep
t2(length(t2(:,1)),:)'
x1 = t2(length(t2(:,1)),1); x2 = t2(length(t2(:,1)),2); x3 = t2(length(t2(:,1)),3); x4 = t2(length(t2(:,1)),4);
X = t2(length(t2(:,1)),1:4)' + eval(V_F)*tstep;
ShujieLandscape1228
t2=[t2; X' k];
ShujieLandscape1228
sum(abs(t2'))
sss =sum(abs(t2'))
sss(1:4)
ShujieLandscape1228
%-- 2020/12/28 11:06 --%
ShujieLandscape1228
mean(real(t2))
ShujieLandscape1228
mean(real(t2))
c2=eigen_vector(:,2)'
period_para=abs(imag(eigen_value(2,2)));
% 这是举例，沿着第二个本征方向上，偏离了 0.01倍 系统会怎样。
t2=[[[x1 x2 x3 x4] + 0.2*c2] 0]
tstep = 0.5; %时间步长取0.02
Anoise = 0.3;
for k=0:tstep:4*pi/period_para
%      x1 = t2(length(t2(:,1)),1); x2 = t2(length(t2(:,1)),2); x3 = t2(length(t2(:,1)),3); x4 = t2(length(t2(:,1)),4);
%      X = [x1 x2 x3 x4]' + eval(V_F)*tstep;
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
%     X=X/sum(abs(X));
%     t2=[t2; X' k];
tseed = rand(1,4); tNoise0 = (tseed - mean(tseed))*Anoise;
AddNoise0 = [x1 x2 x3 x4] + tNoise0;
AddNoise1 = AddNoise0 / sum(abs(AddNoise0));
t2=[t2; AddNoise0 k];
end
figure
ShujieLandscape1228
plot(t2(:,4))
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
mean(real(t2))
mean(abs(t2))
mean(abs(t2(:,1:4)))
sum(mean(abs(t2(:,1:4))))
mean(abs(t2(:,1:4)))
sum(mean(abs(t2(:,1:4))))
mean(abs(t2))
mean(abs(t2(:,1:4)))
%-- 2020/12/28 17:40 --%
ShujieLandscape1228
%-- 2020/12/28 17:56 --%
ZSJ1227TimeSeries
mean(abs(t2(:,1:4)))
mean((t2(:,1:4)))
mean(real(t2(:,1:4)))
view(0,90)
view(90,0)
view(0,0)
ZSJ1227TimeSeries
view(0,90)
view(90,0)
view(0,0)
ZSJ1227TimeSeries
view(0,0)
view(90,0)
view(0,90)
view(0,0)
view(90,0)
view(0,0)
xlim([0.1 0.2])
ylim([0.1 0.2])
zlim([0.1 0.2])
xlim([0.13 0.16]);   zlim([0.13 0.16]);
xlim([0.245 0.265]);   zlim([0.245 0.265]);
xlim([0.235 0.265]);   zlim([0.235 0.265]);  % a=1
xlim([0.285 0.315]);   zlim([0.285 0.315]);  % a=4
xlim([0.285 0.325]);   zlim([0.285 0.325]);  % a=4
ZSJ1227TimeSeries
view(0,0)
xlim([0.13 0.16]);   zlim([0.13 0.16]);  % a=1/4
xlim([0.235 0.265]);   zlim([0.235 0.265]);  % a=1
xlim([0.285 0.325]);   zlim([0.285 0.325]);  % a=4
xlim([0.235 0.265]);   zlim([0.235 0.265]);  % a=1
xlim([0.13 0.16]);   zlim([0.13 0.16]);  % a=1/4
%-- 2020/12/30 13:38 --%
ZSJ1227TimeSeries
mean(t2)
mean(t2(:, 1:4))
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
ZSJ1227TimeSeries
size(kp,1)
ZSJ1227TimeSeries
plotmatrix(rrYret3)
plotmatrix(rrYret3)
t2r = real(t2);
t2r = real(t2(:,1:4));
tt=sum(t2r)'
tt=sum(t2r')'
tt=sum(abs(t2r)')'
tt=sum(abs(t2'))'
tt=sum(abs(t2(:,1:4)'))'
tt=sum((t2(:,1:4)'))'
t2r = real(t2(:,1:4));
tt=sum(t2r')'
plot(tt)
plotmatrix(rrYret3)
plotmatrix(rrYret3)
ZSJ1227TimeSeries
plotmatrix(rrYret3)
plotmatrix(rrYret3)
plotmatrix(rrYret3)
plotmatrix(rrYret3(:,6:end))
ZSJ1227TimeSeries
mean(t2)
4/13
1/13
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
grid on
%-- 2021/1/1 13:05 --%
ZSJ1227TimeSeries
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
ZSJ1227TimeSeries
plot(t2(1:41,1:4),'DisplayName','t2(1:41,1:4)')
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
mean(t2)
8*pi/period_para
mean(t2)
ZSJ1227TimeSeries
for k=0:tstep:18*pi/period_para
%                 seed = 0.01*rand(4,1);
%                 seed = seed - mean(seed);
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
x_addnoise = addNoise2timeseries([x1 x2 x3 x4]', noise_amplit);
t2=[t2; x_addnoise' k];
end
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
51/34
5192/3461
5192/8653
3461/8653
mean(t2)
mean(t2(:,1:4))
mean(t2(1:end/2,1:4))
mean(t2(1:end/3,1:4))
mean(t2(end/3:end,1:4))
ZSJ1227TimeSeries
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
plot(t2(1:41,1:4),'DisplayName','t2(1:41,1:4)')
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
mean(t2(:,1:4))
ZSJ1227TimeSeries
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
18*pi/period_para
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
mean(t2(:,1:4))
ZSJ1227TimeSeries
mean(t2(:,1:4))
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
for k=0:tstep:18*pi/period_para
%                 seed = 0.01*rand(4,1);
%                 seed = seed - mean(seed);
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
x_addnoise = addNoise2timeseries([x1 x2 x3 x4]', noise_amplit);
t2=[t2; x_addnoise' k];
end
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
mean(t2(:,1:4))
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
for k=0:tstep:18*pi/period_para
%                 seed = 0.01*rand(4,1);
%                 seed = seed - mean(seed);
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
x_addnoise = addNoise2timeseries([x1 x2 x3 x4]', noise_amplit);
t2=[t2; x_addnoise' k];
end
mean(t2(:,1:4))
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
ZSJ1227TimeSeries
mean(t2(:,1:4))
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
noise_amplit = 0.05;
for k=0:tstep:18*pi/period_para
%                 seed = 0.01*rand(4,1);
%                 seed = seed - mean(seed);
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
x_addnoise = addNoise2timeseries([x1 x2 x3 x4]', noise_amplit);
t2=[t2; x_addnoise' k];
end
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
rrYret3 = [rrYret3 Yret3];
mean(t2(:,1:4))
for k=0:tstep:18*pi/period_para
%                 seed = 0.01*rand(4,1);
%                 seed = seed - mean(seed);
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
x_addnoise = addNoise2timeseries([x1 x2 x3 x4]', noise_amplit);
t2=[t2; x_addnoise' k];
end
[ Yret3 ] = from_N_colExp_out_am(real(t2(:, 1:4)),mean(real(t2(:, 1:4))))
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
ZSJ1227TimeSeries
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
mean(t2(:,1:4))
for k=0:tstep:18*pi/period_para
%                 seed = 0.01*rand(4,1);
%                 seed = seed - mean(seed);
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
x_addnoise = addNoise2timeseries([x1 x2 x3 x4]', noise_amplit);
t2=[t2; x_addnoise' k];
end
mean(t2(:,1:4))
0.7454/4659
2795/4659
4650/7454
4659/7454
0.7365^2
0.4603^2
0.4603^(1/2)
0.7365^(1/2)
%-- 2021/1/2 13:00 --%
ZSJ1227TimeSeries
size(t3,1)/4
ZSJ1227TimeSeries
view(0,0) % x-z
ZSJ1227TimeSeries
view(0,0) % x-z
xlim([0.1 0.4])
ylim([0.1 0.4])
xlim([0.1 0.4])
ylim([0.1 0.4])
figure(2)
ylim([0.1 0.4])
2lim([0.1 0.4])
zlim([0.1 0.4])
dyntype=1; step = 10; split_price = 8; pic_dan_WZJ20201220(0.20,[10:step:10000],dyntype,split_price)
pic_dan_WZJ20201220
dyntype=1; step = 10; split_price = 8; pic_dan_WZJ20201220(0.20,[10:step:10000],dyntype,split_price)
clear
dyntype=1; step = 40; pic_dan_WZJ20201220(0.20,[90:step:1280],dyntype)
dyntype=1; step = 10; split_price = 8; pic_dan_WZJ20201220(0.20,[10:step:10000],dyntype,split_price)
load snap
load snapshot
WZJreadDan2005_20201219
aa=Wrr(:,3)-Wrr(:,5);
bb=[]; for tI = 1:24;length(aa); bb=[bb;aa(tI:tI+23)'];end
bb=[]; for tI = 1:24:length(aa); bb=[bb;aa(tI:tI+23)'];end
mean(bb)
plot(bb(:,1:3),'DisplayName','bb(:,1:3)')
plot(bb,'DisplayName','bb')
plot(bb(1:100,:),'DisplayName','bb')
plot(bb(1:200,:),'DisplayName','bb')
plot(bb(1:200,1:3:end),'DisplayName','bb')
plot(bb(1:200,1:3:end),'DisplayName','bb','linewidth',4)
plot(bb(1:200,1:3:end),'DisplayName','bb','linewidth',3)
clear
WZJreadDan2005_20201219
aa=Wrr(:,3)-Wrr(:,5);
bb=[]; for tI = 1:24:length(aa); bb=[bb;aa(tI:tI+23)'];end
plot(bb(1:200,1:3:end),'DisplayName','bb','linewidth',3)
c=xcorr(bb(:.1:2),1000)
c=xcorr(bb(:,1:2),1000)
[r lags]=xcorr(bb(:,1:2),1000)
plot(lags)
stem(lags,c)
plot(bb(1:200,1:3:end),'DisplayName','bb','linewidth',3)
[r lags]=xcorr(bb(:,1:2),100)
stem(lags,c)
[r lags]=xcorr(bb(:,1:2),100)
stem(lags,c)
c
stem(lags,r)
[r lags]=xcorr(bb(:,[1 5]),100)
stem(lags,r)
[r lags]=xcorr(bb(:,[1 12]),100)
stem(lags,r)
plot(bb(:,[1 12]))
plot(r(:,1))
grid on
%-- 2021/1/2 17:39 --%
WZJreadDan2005_20201219
aa=Wrr(:,3)-Wrr(:,5);
bb=[]; for tI = 1:24:length(aa); bb=[bb;aa(tI:tI+23)'];end
plot(bb,'DisplayName','bb')
bb1=bb(:,1);
for ti=2:length(bb(:,1))-1; bb1max=[];
if bb1(ti) > bb1(ti-1)  & bb1(ti) > bb1(ti+1)
bb1max=[bb1max bb1(ti)];
end ; end
length(bb(:,1))-1
plot(bb1)
bb1=bb(:,1);bb1max=[];for ti=2:length(bb(:,1))-1;
if bb1(ti) > bb1(ti-1)  & bb1(ti) > bb1(ti+1)
bb1max=[bb1max bb1(ti)];
end ; end
bbmax=[];for tj=1:24; bb1=bb(:,tj);bb1max=[];for ti=2:length(bb(:,1))-1;
if bb1(ti) > bb1(ti-1)  & bb1(ti) > bb1(ti+1)
bb1max=[bb1max bb1(ti)];
end ; end; bbmax=[bbmax;bb1max];end
bbmax=[];for tj=1:24; bb1=bb(:,tj);bb1max=[];for ti=2:length(bb(:,1))-1;
if bb1(ti) > bb1(ti-1)  & bb1(ti) > bb1(ti+1)
bb1max=[bb1max bb1(ti)];
end ; end; bbmax=[bbmax;bb1max(1:2)];end
bbmax=[];for tj=1:24; bb1=bb(:,tj);bb1max=[];for ti=2:length(bb(:,1))-1;
if bb1(ti) > bb1(ti-1)  & bb1(ti) > bb1(ti+1)
bb1max=[bb1max bb1(ti)];
end ; end; bbmax=[bbmax;bb1max(1:3)];end
%-- 2021/1/2 19:22 --%
WZJreadDan2005_20201219
aa=Wrr(:,3)-Wrr(:,5);
bb=[]; for tI = 1:24:length(aa); bb=[bb;aa(tI:tI+23)'];end
plot(bb,'DisplayName','bb')
bb=[]; for tI = 1:24:length(aa); bb=[bb;aa(tI:tI+23)'];end
[ Yret3 ] = from_N_colExp_out_am(bb,mean(bb))
plot(Yret3)
[ Yret3 mn] = from_N_colExp_out_am(bb,mean(bb))
A=[mn Yret3];
bar3(A)
contourf(A)
surf(A)
surf(A)
plot3(A)
plot3(A(:,1),A(:,2),A(:,3),'.')
dd=Wrr(:,3);cc=[]; for tI = 1:24:length(dd); cc=[cc;dd(tI:tI+23)'];end
sum(cc')
m24x24=zeros(24); for it=1:2000-1; for m=1:24;for n=1:24;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
plot(m24x24,'DisplayName','m24x24')
sum(sum(m24x24))
sum(m24x24)
tmp1=[sum(m24x24)' sum(m24x24')']
tmp2= m24x24 - m24x24
tmp2= m24x24 - m24x24'
plot(m24x24)
imagesc(m24x24);
colorbar
imagesc(cc);
colorbar
imagesc(tmp2);
colorbar
surf(tmp2);
hold on
imagesc(tmp2);
T_a= (m24x24 - m24x24')/2
T_s= (m24x24 + m24x24')/2
m24x24./sum(m24x24)
m24x24 ./ sum(m24x24)
m24x24 ./ sum(m24x24)'
m24x24. \ sum(m24x24)'
m24x24. * sum(m24x24)'
m24x24.* sum(m24x24)'
m24x24. + sum(m24x24)'
m24x24. + sum(m24x24).'
m24x24 / sum(m24x24).'
m24x24.* sum(m24x24).'
m24x24.* sum(m24x24)
n24=[]; for it=1:24; n24=[n24; m24x24(it,:)/sum(m24x24(it,:))]; end
[v d]=eig(n24)
[v d]=eig(n24')
diag(d)
dd=diag(d)
plot(v(:,1))
v2=v(:,1)/sum(v(:,1))
scatter(getcolumn(m24x24(:,[2,20]),1),getcolumn(m24x24(:,[2,20]),2))
scatter(getcolumn(m24x24(:,[1,end]),1),getcolumn(m24x24(:,[1,end]),2))
scatter(getcolumn(m24x24(:,[1,13]),1),getcolumn(m24x24(:,[1,13]),2))
scatter(getcolumn(m24x24(:,[1,8]),1),getcolumn(m24x24(:,[1,8]),2))
mesh(m24x24(:,[1,8]))
sum(n24)
sum(n24')
WZJreadDan2005_20201219
dd=diag(d)
plot(v(:,1))
v2=v(:,1)/sum(v(:,1))
plot(v2)
plot(v(:,7))
plot(cc,'DisplayName','cc')
grid on
WZJreadDan2005_20201219
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
tmp2= m24x24 - m24x24'
%eigen of markbian
n24=[]; for it=1:split_price; n24=[n24; m24x24(it,:)/sum(m24x24(it,:))]; end;
[v d]=eig(n24')
plot(v(:,6))
WZJreadDan2005_20201219
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
plot(cc,'DisplayName','cc')
plot(bb,'DisplayName','bb')
WZJreadDan2005_20201219
plot(bb(:,1:7),'DisplayName','bb(:,1:7)')
WZJreadDan2005_20201219
plot(bb,'DisplayName','bb')
WZJreadDan2005_20201219
grid on
plot(Yret3)
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
tmp2= m24x24 - m24x24'
sum(sum(m24x24))
surf(tmp2); hold on; imagesc(tmp2);
sum(sum(abs(tmp2)))
WZJreadDan2005_20201219
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
%plot asym trans
tmp2= m24x24 - m24x24'
surf(tmp2); hold on; imagesc(tmp2);
sum(sum(abs(tmp2)))
%eigen of markbian
n24=[]; for it=1:split_price; n24=[n24; m24x24(it,:)/sum(m24x24(it,:))]; end;
[v d]=eig(n24')
WZJreadDan2005_20201219
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
%plot asym trans
tmp2= m24x24 - m24x24'
surf(tmp2); hold on; imagesc(tmp2);
sum(sum(abs(tmp2)))
sum(sum(abs( m24x24)))
WZJreadDan2005_20201219
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
%plot asym trans
tmp2= m24x24 - m24x24'
surf(tmp2); hold on; imagesc(tmp2);
sum(sum(abs(tmp2)))
n24=[]; for it=1:split_price; n24=[n24; m24x24(it,:)/sum(m24x24(it,:))]; end;
[v d]=eig(n24')
diag(d)
plot(v(:,2))
plot(v(:,2),'o')
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
%plot asym trans
tmp2= m24x24 - m24x24'
surf(tmp2); hold on; imagesc(tmp2);
%eigen of markbian
n24=[]; for it=1:split_price; n24=[n24; m24x24(it,:)/sum(m24x24(it,:))]; end;
[v d]=eig(n24')
diag(d)
plot(v(:,4),'o')
plot(v(:,6),'o')
WZJreadDan2005_20201219
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
%plot asym trans
tmp2= m24x24 - m24x24'
surf(tmp2); hold on; imagesc(tmp2);
sum(sum(abs(tmp2)))
WZJreadDan2005_20201219
aa=Wrr(:,3)-Wrr(:,5)
bb=[]; for tI = 1:split_price:length(aa); bb=[bb;aa(tI:tI+split_price-1)'];end
plot(bb(:,1:6),'DisplayName','bb(:,1:6)')
WZJreadDan2005_20201219
dd=Wrr(:,3);cc=[]; for tI = 1:split_price:length(dd); cc=[cc;dd(tI:tI+split_price-1)'];end
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + cc(it,m)*cc(it+1,n); end;end;end
tmp2= m24x24 - m24x24'
surf(tmp2); hold on; imagesc(tmp2);
sum(sum(abs(tmp2)))
plot(tmp2(1,:))
sum(tmp2(1,:))
sum(abs(tmp2(1,:)))
sum(abs(tmp2(:,1)))
sum((tmp2(:,1)))
plot(tmp2(2,:))
plot(tmp2(3,:))
plot(tmp2(4,:))
plot(tmp2(5,:))
plot(tmp2(6,:))
plot(Yret3)
%-- 2021/1/3 8:20 --%
WZJ8eigencycle
axis ij
axis square
WZJ8eigencycle
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top')
h=pcolor(L_matrix);
colormap(jet(8));
axis ij
axis square
hold on
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
'XAxisLocation','top')
h=pcolor(L_matrix);
colormap(jet(8));
axis ij
axis square
set(xtick, 'XtickLocation','top')
set(gca, 'XAxisLocation','top')
h=pcolor(L_matrix);
colormap(jet(8));
axis ij
axis square
h=pcolor(L_matrix,, 'XAxisLocation','top');
colormap(jet(8));
axis ij
axis square
colorbar;
L_mn = [Omega_mn Yret24]
h=pcolor(L_matrix,  'XAxisLocation','top');
colormap(jet(8));
axis ij
axis square
colorbar;
set(gcf, 'XAxisLocation','top');
hold on;set(gca, 'XAxisLocation','top');
h=pcolor(L_matrix);
colormap(jet(8));
axis ij
axis square
colorbar;
hold on;set(gca, 'XAxisLocation','top');
WZJ8eigencycle
[lon,lat]=meshgrid(1:200,1:200);
contourf(lon,lat,L_matrix);
axis ij
axis square
colorbar;
hold on;set(gca, 'XAxisLocation','top');
colormap(jet(8));
colormap(parula(5))
colormap(hot(5))
colormap(hot(10))
colormap(hot(20))
colormap(hot(10))
colormap(hot(0))
colormap(jet(1))
colormap(jet(10))
colormap(jet(20))
colormap(jet(7))
colormap(jet(12))
h=pcolor(L_matrix);
colormap(jet(8));
axis ij
axis square
colorbar;
hold on;set(gca, 'XAxisLocation','top');
grid off
WZJ8eigencycle
plot(ft,'DisplayName','ft')
xlim([0 100])
WZJ8eigencycle
w2021pricecycle
plot(ft,'DisplayName','ft')
m24x24=zeros(split_price); for it=1:2000-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + ft(it,m)*ft(it+1,n); end;end;end
m24x24=zeros(split_price); for it=1:length(ft(:,1))-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + ft(it,m)*ft(it+1,n); end;end;end
tmp2= m24x24 - m24x24'
surf(tmp2); hold on; imagesc(tmp2);
tmp2= m24x24 - m24x24'
h=pcolor(tmp2);
colormap(jet(8)); axis ij; axis square; colorbar;
hold on;set(gca, 'XAxisLocation','top');
figure
h=pcolor(L_matrix);
colormap(jet(8)); axis ij; axis square; colorbar;
hold on;set(gca, 'XAxisLocation','top');
figure
h=pcolor(m24x24);
colormap(jet(8)); axis ij; axis square; colorbar;
hold on;set(gca, 'XAxisLocation','top');
colormap(jet(48));
h=pcolor(log(m24x24));
clf
h=pcolor(log(m24x24));
colormap(jet(8)); axis ij; axis square; colorbar;
hold on;set(gca, 'XAxisLocation','top');
w2021pricecycle
load snapshot
clear
close all
w2021pricecycle
v
plot(v(:,1))
plot(v(:,2))
plot(v(:,4))
plot(v(:,6))
plot(v(:,8))
plot(v(:,9))
plot(v(:,10))
plot(v(:,13))
dd=diag(d)
plot(dd)
plot(dd,ro'');
plot(dd,'ro');
plot(dd,'ro');grid on;xlabel('real');ylabel('image')
quiver(0,0,real(dd),imag(dd))
quiver(zeros(size(dd,2)),zeros(size(dd,2)),real(dd),imag(dd))
zeros(size(dd,2))
zeros(size(dd,1))
(size(dd,1))
quiver(zeros(size(dd,1)),zeros(size(dd,1)),real(dd),imag(dd))
real(dd)
quiver(zeros(size(dd,1),1),zeros(size(dd,1),1),real(dd),imag(dd))
quiver(zeros(size(dd,1),1),zeros(size(dd,1),1),real(dd),imag(dd),1)
quiver(zeros(size(dd,1),1),zeros(size(dd,1),1),real(dd),imag(dd),1)
plot(dd,'ro');grid on;xlabel('real');ylabel('image')
%gen subspace cycle-value matric
subplot(2,3,1)
h=pcolor(L_matrix);title('Subspace cycle value')
colormap(jet(8)); axis ij; axis square; colorbar;
hold on;set(gca, 'XAxisLocation','top');
%gen trans matric
subplot(2,3,2)
m24x24=zeros(split_price); for it=1:length(ft(:,1))-1; for m=1:split_price;for n=1:split_price;
m24x24(m,n) = m24x24(m,n) + ft(it,m)*ft(it+1,n); end;end;end
h=pcolor(log(m24x24));title('Transit matrix full')
colormap(jet(8)); axis ij; axis square; colorbar;
hold on;set(gca, 'XAxisLocation','top');
%plot asym trans
subplot(2,3,3)
Tasym= m24x24 - m24x24'
h=pcolor(Tasym); title('Transit matrix asym')
colormap(jet(8)); axis ij; axis square; colorbar;
hold on;set(gca, 'XAxisLocation','top');
%eigen of markbian
subplot(2,3,4)
n24=[]; for it=1:split_price; n24=[n24; m24x24(it,:)/sum(m24x24(it,:))]; end;
[v d]=eig(n24')
quiver(zeros(size(dd,1),1),zeros(size(dd,1),1),real(dd),imag(dd),1);hold on
plot(dd,'ro');grid on;xlabel('real');ylabel('image')
quiver(zeros(size(dd,1),1),zeros(size(dd,1),1),real(dd),imag(dd),1);hold on
plot(dd,'ro');grid on;xlabel('real');ylabel('image');axis square;
subplot(2,3,5)
distr=v(:,1)/sum(v(:,1));
histogram(distr)
bar(distr)
distr=v(:,1)/sum(v(:,1)); bar(distr);xlabel('price');ylabel('density');axis square;
subplot(2,3,5)
distr=v(:,1)/sum(v(:,1)); bar(distr);xlabel('price');ylabel('density');axis square;
plot(v(:,2))
plot(v(:,2),'r.')
plot(v(:,2),'r.-')
w2021pricecycle
%-- 2021/1/3 14:08 --%
w2021pricecycle
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=theory_figure202101(dos_x_t)
w2021pricecycle
theory_figure202101(dos_x_t)
w2021pricecycle
axis off
w2021pricecycle
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set] ...
=theory_figure202101(f_x_t,parameter_c_q) ;
w2021pricecycle
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set] ...
=theory_figure202101(f_x_t,parameter_c_q) ;
w2021pricecycle
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set] ...
=theory_figure202101(f_x_t,parameter_c_q) ;
end
w2021pricecycle
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set] ...
=theory_figure202101(f_x_t,parameter_c_q) ;
w2021pricecycle
%-- 2021/1/3 15:52 --%
w2021pricecycle
set(gca, max)
%-- 2021/1/3 19:48 --%
w2021pricecycle
title('time series');hold on;
plot(dos_x_t);
xlabel('time');ylabel('density');axis square;grid on;
xlim([0 100])
m24x24=zeros(split_price);
sum(sum( m24x24 ))
w2021pricecycle
L_v2 = from_eigenvector_out_am(eigenvector_set(:,2))
L_v2 = from_eigenvector_out_am(eigenvector_set(:,2))
L_v2trix = zeros(split_price);tK=1;
for m=1:split_price-1;for n=m+1:split_price;
L_v2trix(m,n) = L_v2(tK);L_v2trix(n,m) = -L_v2(tK);
tK=tK+1;
end;end;
h=pcolor(L_v2trix);
colormap(jet(8)); axis ij; axis square; colorbar;;xlim([1 split_price]);ylim([1 split_price])
hold on;set(gca, 'XAxisLocation','top');
subplot(2,4,8)
%        title('parameter [q c] ');hold on;
%        text(0.1,0.8,strcat('[',num2str(parameter_c_q'),']'));axis off
L_v2 = from_eigenvector_out_am(eigenvector_set(:,2))
L_v2trix = zeros(split_price);tK=1;
for m=1:split_price-1;for n=m+1:split_price;
L_v2trix(m,n) = L_v2(tK);L_v2trix(n,m) = -L_v2(tK);
tK=tK+1;
end;end;
h=pcolor(L_v2trix);
colormap(jet(8)); axis ij; axis square; colorbar;;xlim([1 split_price]);ylim([1 split_price])
hold on;set(gca, 'XAxisLocation','top');
w2021pricecycle
subplot(2,4,8)
%        title('parameter [q c] ');hold on;
%        text(0.1,0.8,strcat('[',num2str(parameter_c_q'),']'));axis off
title('Eigencycle v2 / M');hold on;
w2021pricecycle
%-- 2021/1/4 10:47 --%
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price10-20201212.mat')
totalperiods = 5232;
medianPrice = ones(totalperiods,1);
for tI=1:totalperiods;
medianPrice(tI)=max(find(ps1(tI,:)> 0.5));
end
iqr6 = iqr(psd')';
med_iqr = [medianPrice iqr6];
scatter(med_iqr(:,1),med_iqr(:,2))
plot(med_iqr,'DisplayName','med_iqr')
scatter(med_iqr(:,1),med_iqr(:,2))
line(med_iqr(:,1),med_iqr(:,2))
clear
w2021pricecycle
totalperiods = size(f_x_t,1);
split_price  = size(f_x_t,2);
medianPrice = ones(totalperiods,1);
accumulatedistribution = f_x_t * 0
accumulatedistribution(:,1) = f_x_t(:,1)
for tI=2:split_price;
accumulatedistribution(:,tI)=sum(f_x_t(:,1:f_x_t)')';
end
accumulatedistribution(:,1) = f_x_t(:,1)
for tI=2:split_price;
accumulatedistribution(:,tI)=sum(f_x_t(:,1:tI)')';
end
median_25_50_75 = zeros(totalperiods,3);
for tI=1:totalperiods;
median_25_50_75(tI,1)=max(find(accumulatedistribution(tI,:) < 0.25));
median_25_50_75(tI,2)=max(find(accumulatedistribution(tI,:) < 0.5));
median_25_50_75(tI,2)=max(find(accumulatedistribution(tI,:) < 0.75));
end
w2021pricecycle
totalperiods = size(f_x_t,1);
split_price  = size(f_x_t,2);
accumulatedistribution = f_x_t * 0;
accumulatedistribution(:,1) = f_x_t(:,1)
for tI=2:split_price;
accumulatedistribution(:,tI)=sum(f_x_t(:,1:tI)')';
end
median_25_50_75 = zeros(totalperiods,3);
for tI=1:totalperiods;
median_25_50_75(tI,1)=max(find(accumulatedistribution(tI,:) < 0.25));
median_25_50_75(tI,2)=max(find(accumulatedistribution(tI,:) < 0.5));
median_25_50_75(tI,3)=max(find(accumulatedistribution(tI,:) < 0.75));
end
plot(accumulatedistribution(:,1))
median_25_50_75 = zeros(totalperiods,3);
for tI=1:totalperiods;
for q=1:3
if accumulatedistribution(tI,1) > 0.25*q
median_25_50_75(tI,q) = 1
else
median_25_50_75(tI,q)=max(find(accumulatedistribution(tI,:) < 0.25*q));
end
end
end
plot(median_25_50_75(:,1))
median_25_50_75(:,4) = median_25_50_75(:,3) - median_25_50_75(:,2)
plot(median_25_50_75(:,[1,end]),'DisplayName','median_25_50_75(:,[1,end])')
scatter(getcolumn(median_25_50_75(:,[1,end]),1),getcolumn(median_25_50_75(:,[1,end]),2))
scatter(getcolumn(median_25_50_75(:,[1,end]),1),getcolumn(median_25_50_75(:,[1,end]),2),'r.--')
scatter(getcolumn(median_25_50_75(:,[1,end]),1),getcolumn(median_25_50_75(:,[1,end]),2),'r.-')
andrewsplot(median_25_50_75(:,[1,end]))
compare(iddata(median_25_50_75(:,[1,end])),arx(iddata(median_25_50_75(:,[1,end])),2*eye(2)),1);
compare(iddata(median_25_50_75(:,[1,end])),arx(iddata(median_25_50_75(:,[1,end])),2*eye(2)),1);
compare(iddata(median_25_50_75(:,[1,end])(:,1),median_25_50_75(:,[1,end])(:,2:end)),n4sid(iddata(median_25_50_75(:,[1,end])(:,1),median_25_50_75(:,[1,end])(:,2:end))));
scatter(getcolumn(median_25_50_75(:,[1,end]),1),getcolumn(median_25_50_75(:,[1,end]),2))
scatter(getcolumn(median_25_50_75(:,[1,end]),1),getcolumn(median_25_50_75(:,[1,end]),2))
bar(median_25_50_75(:,[1,end]),'DisplayName','median_25_50_75(:,[1,end])')
w2021pricecycle
dan_median_quantile(f_x_t,parameter_c_q)
totalperiods = size(f_x_t,1);
split_price  = size(f_x_t,2);
accumulatedistribution = f_x_t * 0;
accumulatedistribution(:,1) = f_x_t(:,1)
for tI=2:split_price;
accumulatedistribution(:,tI)=sum(f_x_t(:,1:tI)')';
end
median_25_50_75 = zeros(totalperiods,3);
for tI=1:totalperiods;
for q=1:3
if accumulatedistribution(tI,1) > 0.25*q
median_25_50_75(tI,q) = 1
else
median_25_50_75(tI,q)=max(find(accumulatedistribution(tI,:) < 0.25*q));
end
end
end
median_25_50_75(:,4) = median_25_50_75(:,3) - median_25_50_75(:,2)
bar(median_25_50_75(:,[1,end]),'DisplayName','median_25_50_75(:,[1,end])')
plot(median_25_50_75(:,1),median_25_50_75(:,4),'r-')
bar(median_25_50_75(:,[1,end]),'DisplayName','median_25_50_75(:,[1,end])')
bar(median_25_50_75(:,1),median_25_50_75(:,4))
bar(1:1000;median_25_50_75(:,1),median_25_50_75(:,4))
bar(1:1000,median_25_50_75(:,1),median_25_50_75(:,4))
bar(median_25_50_75(:,[1,end]),'DisplayName','median_25_50_75(:,[1,end])')
bar(median_25_50_75(:,[1,4]),'DisplayName','median_25_50_75(:,[1,end])')
plot(median_25_50_75(:,1),median_25_50_75(:,4))
scatter(median_25_50_75(:,1),median_25_50_75(:,4))
w2021pricecycle
totalperiods = size(f_x_t,1);
split_price  = size(f_x_t,2);
accumulatedistribution = f_x_t * 0;
accumulatedistribution(:,1) = f_x_t(:,1)
for tI=2:split_price;
accumulatedistribution(:,tI)=sum(f_x_t(:,1:tI)')';
end
median_25_50_75 = zeros(totalperiods,3);
for tI=1:totalperiods;
for q=1:3
if accumulatedistribution(tI,1) > 0.25*q
median_25_50_75(tI,q) = 1
else
median_25_50_75(tI,q)=max(find(accumulatedistribution(tI,:) < 0.25*q));
end
end
end
median_25_50_75(:,4) = median_25_50_75(:,3) - median_25_50_75(:,2)
bar(median_25_50_75(:,[1,4]),'DisplayName','median_25_50_75(:,[1,end])')
ylim([0 10])
plot(median_25_50_75(:,1),median_25_50_75(:,4))
plot(median_25_50_75(:,1),median_25_50_75(:,4)*10)
plot(median_25_50_75(:,1)*10,median_25_50_75(:,4))
plot(median_25_50_75(:,1),median_25_50_75(:,4),'ro')
plot(median_25_50_75(:,2),median_25_50_75(:,4),'ro')
totalperiods = size(f_x_t,1);
split_price  = size(f_x_t,2);
accumulatedistribution = f_x_t * 0;
accumulatedistribution(:,1) = f_x_t(:,1)
for tI=2:split_price;
accumulatedistribution(:,tI)=sum(f_x_t(:,1:tI)')';
end
median_25_50_75 = zeros(totalperiods,3);
for tI=1:totalperiods;
for q=1:3
if accumulatedistribution(tI,1) > 0.25*q
median_25_50_75(tI,q) = 1;
else
median_25_50_75(tI,q)=max(find(accumulatedistribution(tI,:) < 0.25*q));
end
end
end
median_25_50_75(:,4) = median_25_50_75(:,3) - median_25_50_75(:,1)
plot(median_25_50_75(:,2),median_25_50_75(:,4),'ro')
w2021pricecycle
plot(median_25_50_75(:,2),median_25_50_75(:,4),'ro')
plot(dan_median_quantile,'ro')
plot(dan_median_quantile)
plot(log(dan_median_quantile))
w2021pricecycle
plot(log(dan_median_quantile))
w2021pricecycle
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]  = w2021pricecycle()
w2021pricecycle
close all
dan_median_quantile=dan_median_quantile(f_x_t,parameter_c_q)
w2021pricecycle
line(dan_median_quantile)
line(dan_median_quantile(:,1), dan_median_quantile(:,2))
line(dan_median_quantile(1:20,1), dan_median_quantile(1:20,2))
line(dan_median_quantile(10:20,1), dan_median_quantile(10:20,2))
text(dan_median_quantile(10:20,1), dan_median_quantile(10:20,2),num2str(10:20))
for j=10:20
text(dan_median_quantile(j,1), dan_median_quantile(j,2),num2str(j));end
w2021pricecycle
for j=10:20
text(dan_median_quantil(j,1), dan_median_quantil(j,2),num2str(j));end
for j=10:20
text(dan_median_quantil(j,1), dan_median_quantil(j,2),num2str(j));hold on;
end
dan_median_quantil(j,1), dan_median_quantil(j,2)
for j=10:20
plot(dan_median_quantil(j,1), dan_median_quantil(j,2),'ro');hold on;
text(dan_median_quantil(j,1), dan_median_quantil(j,2),num2str(j));hold on;
end
for j=10:15
plot(dan_median_quantil(j,1), dan_median_quantil(j,2),'ro');hold on;
text(dan_median_quantil(j,1), dan_median_quantil(j,2),num2str(j));hold on;
end
clf
x=10
for j=10+x:15+x
plot(dan_median_quantil(j,1), dan_median_quantil(j,2),'ro');hold on;
text(dan_median_quantil(j,1), dan_median_quantil(j,2),num2str(j));hold on;
end
clf
x=100
for j=10+x:18+x
plot(dan_median_quantil(j,1), dan_median_quantil(j,2),'ro');hold on;
text(dan_median_quantil(j,1), dan_median_quantil(j,2),num2str(j));hold on;
end
%-- 2021/1/4 14:02 --%
w2021pricecycle
scatter(getcolumn(median_25_50_75(:,[2,end]),1),getcolumn(median_25_50_75(:,[2,end]),2))
clf
x=100
for j=10+x:18+x
plot(dan_median_quantil(j,1), dan_median_quantil(j,2),'ro');hold on;
text(dan_median_quantil(j,1), dan_median_quantil(j,2),num2str(j));hold on;
end
plot(dan_median_quantil(1:100,1), dan_median_quantil(1:100,2),'ro');
plot(dan_median_quantil,'DisplayName','dan_median_quantil')
plot(dan_median_quantil(1:20,:),'DisplayName','dan_median_quantil')
plot(dan_median_quantil(1:120,:),'DisplayName','dan_median_quantil')
grid on
clf
load('F:\Spectrum2020\datasource\L_28.mat')
plot(dan_median_quantil(1:20,:),'DisplayName','dan_median_quantil')
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( L28_all15326)
load('F:\Spectrum2020\datasource\L_28.mat')
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( L28_all15326)
load('F:\Spectrum2020\datasource\L_28.mat')
clear
load('F:\Spectrum2020\datasource\L_28.mat')
load('F:\Spectrum2020\datasource\oneill15327.mat')
d=data15327(:,4:11);
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( d)
sum(distr(1:4))
z=sum(distr(1:4))
z=sum(distr(5:8))
close all
timesies2figure8
w2021pricecycle
load('F:\Spectrum2020\datasource\oneill15327.mat')
d=data15327(:,4:11);
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( d)
w2021pricecycle
timesies2figure8
d=data15327(:,4:11);timesies2figure8
d=data15327(:,4:11);timesies2figure8(d)
pcolor(1:8,1:8,Tasym)
pcolor(1:9,1:9,Tasym)
pcolor(0.5:8.5,0.5:8.5,Tasym)
[X,Y]=meshgrid(0.5:1:8.5,0.5:1:8.5);
pcolor(X,Y,Tasym);
colorbar;
colormap(flipud(hot));
timesies2figure8
d=data15327(:,4:11);timesies2figure8(d)
timesies2figure8
d=data15327(:,4:11);timesies2figure8(d)
timesies2figure8
d=data15327(:,4:11);timesies2figure8(d)
timesies2figure8
d=data15327(:,4:11);timesies2figure8(d)
scatter(r3col_Theo28(:,1),C(:,7),129,'O'); kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');axis square;grid on;hold on;
testsave2Dfig
load('F:\Spectrum2020\datasource\oneill15327.mat')
d=data15327(:,4:11);
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( d)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( d)
suptitle('总标题')
suptitle('总标题','fontsize',25)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( d)
suptitle('总标题')
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8( d)
suptitle('总标题')
%-- 2021/1/5 15:00 --%
N_cycle0315
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
ShujieLandscape1228
[ Yret3 mn] = from_N_colExp_out_am(t2(:,1:4),mean(t2(:,1:4)))
Lmn = from_eigenvector_out_am(eigen_vector(:,2))
scatter(Lmn,Yret3)
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
abs(eigen_vector(2))
abs(eigen_vector(:,2))
eigen_vector
ShujieLandscape1228
t2=[[[x1 x2 x3 x4] + 0.02*real(c2)] 0]
tstep = 0.02; %时间步长取0.02
Anoise = 0;
for k=0:tstep:2*pi/period_para
%      x1 = t2(length(t2(:,1)),1); x2 = t2(length(t2(:,1)),2); x3 = t2(length(t2(:,1)),3); x4 = t2(length(t2(:,1)),4);
%      X = [x1 x2 x3 x4]' + eval(V_F)*tstep;
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
%     X=X/sum(abs(X));
%     t2=[t2; X' k];
tseed = rand(1,4); tNoise0 = (tseed - mean(tseed))*Anoise;
AddNoise0 = [x1 x2 x3 x4] + tNoise0;
AddNoise1 = AddNoise0 / sum(abs(AddNoise0))
t2=[t2; AddNoise1 k];
end
[ Yret3 mn] = from_N_colExp_out_am(t2(:,1:4),mean(t2(:,1:4)))
Lmn = from_eigenvector_out_am(eigen_vector(:,2))
scatter(Lmn,Yret3)
abs(eigen_vector(:,2))
mean_U
ut=[];
for ampt=0.01:0.01:0.03
utt=[];
t2=[[[x1 x2 x3 x4] + ampt*real(c2)] 0]
tstep = 0.02; %时间步长取0.02
Anoise = 0;
for k=0:tstep:2*pi/period_para
%      x1 = t2(length(t2(:,1)),1); x2 = t2(length(t2(:,1)),2); x3 = t2(length(t2(:,1)),3); x4 = t2(length(t2(:,1)),4);
%      X = [x1 x2 x3 x4]' + eval(V_F)*tstep;
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
%     X=X/sum(abs(X));
%     t2=[t2; X' k];
tseed = rand(1,4); tNoise0 = (tseed - mean(tseed))*Anoise;
AddNoise0 = [x1 x2 x3 x4] + tNoise0;
AddNoise1 = AddNoise0 / sum(abs(AddNoise0))
t2=[t2; AddNoise1 k];
utt=[utt;eval(mean_U)];
end
ut=[ut; ampt mean(utt)];
[ Yret3 mn] = from_N_colExp_out_am(t2(:,1:4),mean(t2(:,1:4)))
Lmn = from_eigenvector_out_am(eigen_vector(:,2))
figure; scatter(Lmn,Yret3);title(ampt)
end
ShujieLandscape1228
0.25/1.75
4/13
ShujieLandscape1228
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
2*pi/period_para
ShujieLandscape1228
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
ShujieLandscape1228
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
ShujieLandscape1228
plot(utt)
mean(utt)
ShujieLandscape1228
clear
ShujieLandscape1228
clear
ShujieLandscape1228
mean_U
tmp20210105
mean(u)
tmp20210105
mean(u)
tmp20210105
mu = mean(u)'
tmp20210105
mu = mean(u)'
tmp20210105
mu = mean(u)'
%-- 2021/1/6 0:42 --%
tmp20210105
expand(Um)
D_Um = [diff(Um,'x1') diff(Um,'x2') diff(Um,'x3') diff(Um,'x4')]
D_Um = [diff(Um,'x1') diff(Um,'x2') diff(Um,'x3') diff(Um,'x4')]'
xx=[ 0 + a*(0 + 1/(3*a + 1)) + a/(3*a + 1)
0 + 0 + (2*a)/(3*a + 1)
0 + 0 + a/(3*a + 1) + 1/(3*a + 1)
0 + a/(3*a + 1) + a*(0 + a/(3*a + 1))]
xx*(3*a+1)
simplify(xx*(3*a+1))
N_cycle0315
D_Um = [diff(mean_U,'x1') diff(mean_U,'x2') diff(mean_U,'x3') diff(mean_U,'x4')]'
W = subs(subs(subs(subs(D_Um,'x3','x2'),'x1','x2'),'x4','x2/a'),'x2','a/(3*a + 1)')
a=0.25;eval(W)
a=4;eval(W)
a=1;eval(W)
%-- 2021/1/6 7:26 --%
N_cycle0315
A.x1
S.x1
S.x2
S.x3
S.x4
N_cycle0315
N_cycle0315_Price6
Mix5_2000
eigen_vector
Mix5_2000
Payoff_vector_field_F.^3
Mix5_2000
eigen_vector
Mix5_2000
%-- 2021/1/6 13:29 --%
OneillYMQ20201118
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
OneillYMQ20201118
find(prod(A')>0)
INs = find(prod(A')>0);
x1=A(INs,1);  x2=A(INs,2);
x3=A(INs,3);  x4=A(INs,4);
y1=A(INs,5);  y2=A(INs,6);
y3=A(INs,7);  y4=A(INs,8);
D_Eq_at_NE = eval(D_Eq_0);
% A6= D_Eq_0([1:3 5:7],[1:3 5:7])
%         D_Eq_at_NE = eval(D_Eq_0)
%          D_Eq_at_NE = eval(A6)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn from_eigenvector_out_am(v)];
end;
Em =diag(eigen_value)
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn from_eigenvector_out_am(eigen_vector(:,k))];
end;
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn round(from_eigenvector_out_am(eigen_vector(:,k)),4)];
end;
OneillYMQ20201118
clear
OneillYMQ20201118
%-- 2021/1/6 13:57 --%
OneillYMQ20201118
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_0);
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn round(from_eigenvector_out_am(eigen_vector(:,k)),4)];
end;
Em =diag(eigen_value)
OneillYMQ20201118
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_0);
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn round(from_eigenvector_out_am(eigen_vector(:,k)),4)];
end;
Em =diag(eigen_value)
x1=1/5;x2=2/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_0);
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn round(from_eigenvector_out_am(eigen_vector(:,k)),4)];
end;
Em =diag(eigen_value)
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_0);
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn round(from_eigenvector_out_am(eigen_vector(:,k)),4)];
end;
Em =diag(eigen_value)
OneillYMQ20201118
Em =diag(eigen_value)
OneillYMQ20201118
Em =diag(eigen_value)
OneillYMQ20201118
w2021pricecycle
load('F:\Spectrum2020\datasource\oneill15327.mat'); d=data15327(:,4:11);
timesies2figure8(d);
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
plot(ft(1:200,:),'DisplayName','ft')
plot(ft(1:600,:),'DisplayName','ft')
price20210115
c = 0.6; q1 = 2/3; q2=1/3; q3=0;
highp = c/(1-q1)                 %**************** 1
lowp = highp*q1/(q1+2*q2+3*q3)   %**************** 2
%S:
c = 0.36; q1 = 0.8; q2=0.1; q3=0.1;
highp = c/(1-q1)                 %**************** 1
lowp = highp*q1/(q1+2*q2+3*q3)   %**************** 2
%-- 2021/1/7 0:12 --%
price20210115
plot(ft,'DisplayName','ft')
plot(ft(1:61,:),'DisplayName','ft(1:61,:)')
price20210115
plot(ft,'DisplayName','ft')
n=5
(1-(1:n)/n)
(1-(1:n)/n).^2
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
p
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
p
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
plot(ft,'DisplayName','ft')
price20210115
p(12)
price20210115
p(1/n)
p(6)
price20210115
plot(ft(end,:))
f.*p
f*p'
price20210115
ft = price20210115(1)
ft = price20210115(2)
ft = price20210115(1)
f*p'
ft = price20210115(1)
price20210115
ft = price20210115(1)
ft = price20210115(2)
ft = price20210115(1)
%-- 2021/1/7 10:36 --%
ft = price20210115(1)
ft = price20210115(2)
%-- 2021/1/7 16:00 --%
seed=rand(1,10)
price20210115addnoise
ft = price20210115(2)
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft(1000:end,:));
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft(500:end,:));
price20210115addnoise
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(2)
ft = price20210115addnoise(1)
price20210115addnoise
ft = price20210115addnoise(1)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft(500:2000,:));
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft(500:1000,:));
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft(500:100:100000,:));
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft(500:10:10000,:));
price20210115addnoise
ft = price20210115addnoise(1)
price20210115addnoise
ft = price20210115addnoise(1)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft);
price20210115addnoise
ft = price20210115addnoise(1)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft);
ft = price20210115addnoise(2)
ft = price20210115addnoise(1)
ft = price20210115addnoise(2)
price20210115addnoise
ft = price20210115addnoise(1)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft);
ft = price20210115addnoise(2)
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft);
%-- 2021/1/9 10:12 --%
Mix3_2000
for tI=1:3; subplot(1,3,tI);title(num2str(eigen_value(tI,tI)));axis square;hold on
for tJ=1:3
quiver(0,0,real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),1,'LineWidth',tJ);hold on
text(real(eigen_vector(tJ,tI)),imag(eigen_vector(tJ,tI)),num2str(eigen_vector(tJ,tI)));hold on
end
end
theo_eigcyc_Area=[];
for J=1:3;
theo_subspace_id_mn=zeros(3,3); tmpK=1;
tmpArea=[];
B = eigen_vector(:,J);
for m=1:7;
for n=m+1:8;
theo_subspace_id_mn(tmpK,:)=[m*10+n m n];tmpK=tmpK+1;
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn);
tmpArea=[tmpArea; area_m_n];
end;
end
theo_eigcyc_Area=[theo_eigcyc_Area tmpArea];
end
%% 生成 eigencycle set
theo_eigcyc_Area=[];
for J=1:3;
theo_subspace_id_mn=zeros(3,3); tmpK=1;
tmpArea=[];
B = eigen_vector(:,J);
for m=1:2;
for n=m+1:3;
theo_subspace_id_mn(tmpK,:)=[m*10+n m n];tmpK=tmpK+1;
vm=B(m);
vn=B(n);
[area_m_n, Em, En] = lissa2Complex2areaPlot(vm,vn);
tmpArea=[tmpArea; area_m_n];
end;
end
theo_eigcyc_Area=[theo_eigcyc_Area tmpArea];
end
%-- 2021/1/9 14:54 --%
addpath 'F:\Spectrum2020\DynEquation'
type lotkavolterra
sol1 = ode23(@lotkavolterra,[0 20],[ 20 20 ]);
sol2 = ode23(@lotkavolterra,[0 20],[ 100 100 ]);
xint = linspace(0,20,500);
yint1 = deval(sol1,xint); yint2 = deval(sol2,xint);
subplot(121)
plot(xint,yint1,'LineW',1,'Color','k'); hold on
plot(xint,yint2,'LineW',2,'Color',[0.6 0.6 0.6])
axis square
xlabel('time'); ylabel('population');
subplot(122)
plot(yint1(1,:),yint1(2,:),'LineW',1,'Color','k'); hold on
plot(yint2(1,:),yint2(2,:),'LineW',2,'Color',[0.6 0.6 0.6])
axis square
xlabel('number of preys'); ylabel('number of predators');
plot(sol1.x)
run('F:\Spectrum2020\DynEquation\multiscaleDemo.m')
multiscaleDemo
close all
plot(t2,y2,'LineW',1,'Color',[0.6 0.6 0.6])
axis([0 20 0 350])
axis square
plot(t1,y1,'LineW',4,'Color','k'); hold on
figure(5)
subplot(121)
plot(t1,y1,'LineW',4,'Color','k'); hold on
plot(t2,y2,'LineW',1,'Color',[0.6 0.6 0.6])
axis([0 20 0 350])
axis square
axis square
xlabel('number of preys'); ylabel('number of predators');
%-- 2021/1/9 20:16 --%
multiscaleDemo
addpath 'F:\Spectrum2020\DynEquation'
type lotkavolterra
multiscaleDemo
set(gca,'FontSize',14)
multiscaleDemo
type sdeEuler
clear
h=.01;
sh=sqrt(h);
t=tspan(1):h:tspan(2);
y=zeros(length(ic),length(t));
y(:,1)=ic;
stopcond=false;
for i=1:length(t)-1,
y(:,i+1)=y(:,i)+h*feval(f,t,y(:,i),varargin{:})+sh*sdnoise.*randn(size(sdnoise));
if stopcond && any(y(:,i+1)<0)
yz =find(y(:,i+1)<0);
y(yz,i+1)=0;
sdnoise(yz)=0;
end
end
clear
multiscaleDemo
y
%-- 2021/1/10 13:37 --%
sqrt{3}
sqrt(3)
a=[0	0	0	0
1/6 	 1/6 	 1/6 	 1/6
(-1+1.7321 i)/12 	  (-1-1.7321i)/12 	  (-1-1.7321i)/12 	  (-1+1.7321i)/12
(-1 -1.7321 i)/12 	  (-1+1.7321i)/12 	  (-1+1.7321i)/12 	  (-1-1.7321i)/12
0	0	0	0
- i/6 	 i/6 	 - i/6 	 i/6
(1.7321 + i)/12 	  (1.7321-i)/12 	  (-1.7321+i)/12 	  (-1.7321-i)/12
(-1.7321+i)/12 	  (-1.7321-i)/12 	  (1.7321+i)/12 	  (1.7321-i)/12
]
[0	0	0	0
1/6 	 1/6 	 1/6 	 1/6
(-1+1.7321i)/12 	  (-1-1.7321i)/12 	  (-1-1.7321i)/12 	  (-1+1.7321i)/12
(-1 -1.7321i)/12 	  (-1+1.7321i)/12 	  (-1+1.7321i)/12 	  (-1-1.7321i)/12
0	0	0	0
- i/6 	 i/6 	 - i/6 	 i/6
(1.7321 + i)/12 	  (1.7321-i)/12 	  (-1.7321+i)/12 	  (-1.7321-i)/12
(-1.7321+i)/12 	  (-1.7321-i)/12 	  (1.7321+i)/12 	  (1.7321-i)/12
]>
a=[0	0	0	0
1/6 	 1/6 	 1/6 	 1/6
(-1+1.7321i)/12 	  (-1-1.7321i)/12 	  (-1-1.7321i)/12 	  (-1+1.7321i)/12
(-1 -1.7321i)/12 	  (-1+1.7321i)/12 	  (-1+1.7321i)/12 	  (-1-1.7321i)/12
0	0	0	0
- i/6 	 i/6 	 - i/6 	 i/6
(1.7321 + i)/12 	  (1.7321-i)/12 	  (-1.7321+i)/12 	  (-1.7321-i)/12
(-1.7321+i)/12 	  (-1.7321-i)/12 	  (1.7321+i)/12 	  (1.7321-i)/12
]
a=[ 1/6 	 1/6 	 1/6 	 1/6
(-1+1.7321i)/12 	  (-1-1.7321i)/12 	  (-1-1.7321i)/12 	  (-1+1.7321i)/12
(-1 -1.7321i)/12 	  (-1+1.7321i)/12 	  (-1+1.7321i)/12 	  (-1-1.7321i)/12
]
b=[ -i/6 	i/6 	 -i/6 	i/6
(1.7321 +i)/12 	  (1.7321-i)/12 	  (-1.7321+i)/12 	  (-1.7321-i)/12
(-1.7321+i)/12 	  (-1.7321-i)/12 	  (1.7321+i)/12 	  (1.7321-i)/12
]
a=[0	0	0	0
1/6 	 1/6 	 1/6 	 1/6
(-1+1.7321i)/12 	  (-1-1.7321i)/12 	  (-1-1.7321i)/12 	  (-1+1.7321i)/12
(-1 -1.7321i)/12 	  (-1+1.7321i)/12 	  (-1+1.7321i)/12 	  (-1-1.7321i)/12
0	0	0	0
-i/6 	i/6 	 -i/6 	i/6
(1.7321 +i)/12 	  (1.7321-i)/12 	  (-1.7321+i)/12 	  (-1.7321-i)/12
(-1.7321+i)/12 	  (-1.7321-i)/12 	  (1.7321+i)/12 	  (1.7321-i)/12
]
alphabeta0110
close all
alphabeta0110
close all
alphabeta0110
xlabel(num2str(m),'fontsize',6,'location','top');
xlabel(num2str(m),'fontsize',6,'position','top');
alphabeta0110
%-- 2021/1/11 0:25 --%
alphabeta0110
invariantEVec
a3=round(vec_rr,3)
invariantEVec
vec_rr
invariantEVec
vec_rr
invariantEVec
vec_rr
%-- 2021/1/11 11:37 --%
invariantEVec
d=real(0.02*eigen_vector(:,8));x1=x1+d(1);
x2=x2+d(2);
x3=x3+d(3);
x4=x4+d(4);y1=y1+d(5);
y2=y2+d(6);
y3=y3+d(7);
y4=y4+d(8);
eval([mean_U_1 mean_U_2])
invariantEVec
d=real(0.02*eigen_vector(:,8));x1=x1+d(1);
x2=x2+d(2);
x3=x3+d(3);
x4=x4+d(4);y1=y1+d(5);
y2=y2+d(6);
y3=y3+d(7);
y4=y4+d(8);
eval([mean_U_1 mean_U_2])
mean_U_1
d=real(0*eigen_vector(:,8));x1=x1+d(1);
x2=x2+d(2);
x3=x3+d(3);
x4=x4+d(4);y1=y1+d(5);
y2=y2+d(6);
y3=y3+d(7);
y4=y4+d(8);
eval([mean_U_1 mean_U_2])
mean_U_1 = [x1 x2 x3 x4] * Payoff_vector_field_F_1;
mean_U_2 = [y1 y2 y3 y4] * Payoff_vector_field_F_2;
d=real(0*eigen_vector(:,8));x1=x1+d(1);
x2=x2+d(2);
x3=x3+d(3);
x4=x4+d(4);y1=y1+d(5);
y2=y2+d(6);
y3=y3+d(7);
y4=y4+d(8);
eval([mean_U_1 mean_U_2])
d=real(0*eigen_vector(:,1));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
invariantEVec
eval([mean_U_1 mean_U_2])
d=real(0*eigen_vector(:,1));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
d=real(0.01*eigen_vector(:,1));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
invariantEVec
d=real(0.1*eigen_vector(:,1));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
invariantEVec
eval([mean_U_1 mean_U_2])
d=real(0.1*eigen_vector(:,8));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.1*eigen_vector(:,2));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.1*eigen_vector(:,3));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.1*eigen_vector(:,4));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.1*eigen_vector(:,5));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.1*eigen_vector(:,6));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.1*eigen_vector(:,7));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.01*eigen_vector(:,7));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
d=real(0.01*eigen_vector(:,1));
x1=x1+d(1);x2=x2+d(2);x3=x3+d(3);x4=x4+d(4);
y1=y1+d(5);y2=y2+d(6);y3=y3+d(7);y4=y4+d(8);
eval([mean_U_1 mean_U_2])
eval([mean_U_1 mean_U_2 Payoff_vector_field_F_1 Payoff_vector_field_F_2])
eval([Payoff_vector_field_F_1 Payoff_vector_field_F_2])
eval([V_Eq_0])
invariantEVec
abs(eigen_vector(:,2)).^2
abs(eigen_vector(:,2))
eval(mean_U_1)
eigen_vector(:,2)'*eigen_vector(:,2)
eigen_vector(:,2)*eigen_vector(:,2)'
abs(eigen_vector(:,2)*eigen_vector(:,2)')
abs(eigen_vector(:,2)*eigen_vector(:,2)').^2
eval([mean_U_1  mean_U_2])
invariantEVec
clear
invariantEVec
abs(eigen_vector(:,5)).^2
abs(eigen_vector(:,7)).^2
invariantEVec
D_Eq_at_NE = eval(D_Eq_0) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
%-- 2021/1/11 13:20 --%
invariantEVec
D_Eq_1 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ...
diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3') diff(V_Eq_0,'y4') ...
];
invariantEVec
D_Eq_at_NE = eval(D_Eq_0) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
D_Eq_at_NE = eval(D_Eq_2) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
invariantEVec
D_Eq_at_NE = eval(D_Eq_2) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_1) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
invariantEVec
D_Eq_1 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ...
diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3') diff(V_Eq_0,'y4') ...
];
D_Eq_at_NE = eval(D_Eq_1) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
invariantEVec
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_1) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
invariantEVec
D_Eq_1 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ...
diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3') diff(V_Eq_0,'y4') ...
]
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_1) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
invariantEVec
abs(eigen_vector)
x1=2/5/2;x2=1/5;x3=1/5;x4=1/5/2;y1=2/5/2;y2=1/5/2;y3=1/5/2;y4=1/5/2;
D_Eq_at_NE = eval(D_Eq_1) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
invariantEVec
D_Eq_at_NE = eval(D_Eq_1) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
vec_rr = [vec_rr; eigen_vector];
x1=2/5;x2=1/5;x3=1/5;x4=1/5;y1=2/5;y2=1/5;y3=1/5;y4=1/5;
D_Eq_at_NE = eval(D_Eq_1) ;
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
invariantEVec
eval(mean_U_1)
invariantEVec
%-- 2021/1/11 16:43 --%
invariantEVec
x5_YQM
plot(eigen_vector(:,3))
angle(eigen_vector(:,3))
quiver(0,0,real(eigen_vector(:,3)),image(eigen_vector(:,3)),1)
quiver(0,0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
x5_YQM
x1=444/2987;x2=22/103;x3=641/2987;x4=288/2987;x5=976/2987
x5_YQM
Payoff_vector_field_F
%-- 2021/1/11 17:37 --%
x5_YQM
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
x5_YQM
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
x5_YQM
invariantEVec
x5_YQM
invariantEVec
x5_YQM
eigen_vector(:,1) / sum(eigen_vector(:,1))
eigen_vector(:,1) / sum(eigen_vector(:,1)) - [x1 x2 x3 x4 x5]'
x5_YQM
S=solve(V_Eq_0)
S=solve(V_Eq_0(1)=V_Eq_0(2), V_Eq_0(1)=V_Eq_0(3) , V_Eq_0(1)=V_Eq_0(4) , V_Eq_0(1)=V_Eq_0(5), 'x1+x2+x3+x4+x5=1' );
S=solve(V_Eq_0(1)-V_Eq_0(2), V_Eq_0(1)-V_Eq_0(3) , V_Eq_0(1)-V_Eq_0(4) , V_Eq_0(1)-V_Eq_0(5), 'x1+x2+x3+x4+x5=1' );
[S.x1 S.x2 S.x3 S.x4 S.x5]
V_Eq_0(1)
V_Eq_0(1)-V_Eq_0(2)
S=solve(V_Eq_0(1)-V_Eq_0(2), V_Eq_0(1)-V_Eq_0(3) , V_Eq_0(1)-V_Eq_0(4) , V_Eq_0(1)-V_Eq_0(5), 'x1+x2+x3+x4+x5-1' );
[S.x1 S.x2 S.x3 S.x4 S.x5]
V_Eq_0(1:4)-V_Eq_0(2:5)
S=solve(V_Eq_0(1:4)-V_Eq_0(2:5),x1+x2+x3+x4+x5-1)
S=solve([V_Eq_0(1:4)-V_Eq_0(2:5)]',x1+x2+x3+x4+x5-1)
S=solve( V_Eq_0(1:4)-V_Eq_0(2:5) )
[S.x1 S.x2 S.x3 S.x4 S.x5]
[S.x1 S.x2 S.x3 S.x4]
x5_YQM
x1=444/2987;x2=22/103;x3=641/2987;x4=298/2987;x5=966/2987;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
x1=444/2987;x2=22/103;x3=641/2987;x4=288/2987;x5=976/2987;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
%-- 2021/1/12 13:10 --%
timesies2figure8
load('F:\Spectrum2020\datasource\oneill15327.mat')
d=data15327(:,4:11)
timesies2figure8(d)
timesies2figure8_arxiv2
load('F:\Spectrum2020\datasource\oneill15327.mat'); d=data15327(:,4:11);
timesies2figure8(d);
timesies2figure8_arxiv2(d)  % 3 figures
timesies2figure8_arxiv2
timesies2figure8_arxiv2(d)
% timesies2figure8_arxiv2(d)  % 3 figures
timesies2figure8_arxiv2(d)  % 3 figures
w2021pricecycle
ft = price20210115(1)
ft = price20210115(2)
ft = price20210115addnoise(2)
ft = price20210115addnoise(2)
ft = price20210115(1)
timesies2figure8_arxiv2(ft)  % 3 figures
ft = price20210115(1)
timesies2figure8_arxiv2(ft)  % 3 figures
price20210115
ft = price20210115(1)
timesies2figure8_arxiv2(ft)  % 3 figures
ft = price20210115(2)
timesies2figure8_arxiv2(ft)  % 3 figures
%-- 2021/1/13 17:34 --%
x5_YQM
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
x5_YQM
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
abs(eigen_vector(:,3))
S=solve(V_Eq_0)
x5_YQM
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3 S.x4 S.x5]
D_Eq_at_NE - D_Eq_at_NE'
payoff_matrix
payoff_matrix - payoff_matrix'
x5_YQM
eval(V_Eq_0)
eval(Payoff_vector_field_F)
D_Eq_at_NE
jordan(D_Eq_at_NE)
[V,J]=jordan(D_Eq_at_NE)
A = [ 1 -3 -2;
-1  1 -1;
2  4  5]; [vA dA] = eig(A)
alphabeta0110
ev(:,1)'*ev(:,2)
for m=1:8;for n=1:8;test(m,n) = ev(:,m)'*ev(:,n);end;end
alphabeta0110
0.6124 * 0.378
ev=[0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
-0.1964	0.1964	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
]
for m=1:8;for n=1:8;test(m,n) = ev(:,m)'*ev(n,:);end;end
for m=1:8;for n=1:8;test(m,n) = ev(:,m)'*ev(:,n);end;end
ev=[0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
-0.1964	0.1964	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
];
for m=1:8;for n=1:8;test(m,n) = ev(:,m)'*ev(:,n);end;end
ev=[0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
-0.1964	0.1964	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
]
alphabeta0110
for m=1:8;for n=1:8;test(m,n) = ev(:,m)'*ev(:,n);end;end
OneillYMQ20201118
[V,J]=jordan(D_Eq_at_NE)
ev=V; for m=1:8;for n=1:8;test(m,n) = ev(:,m)'*ev(:,n);end;end
OneillYMQ20201118
Po=eye(8)
Po=eye(8)-ones(8)/8
D_Eq_at_NE = Po*D_Eq_at_NE
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
[V J]=jordan(D_Eq_at_NE)
[Ve Je]=eig(D_Eq_at_NE)
OneillYMQ20201118
[V,J]=jordan(D_Eq_at_NE)
cond = J == V\D_Eq_at_NE*V;
isAlways(cond)
eigen_vector(:,3)^2
eigen_vector(:,3).^2
eigen_vector(:,3).^(1/2)
clear
OneillYMQ20201118
[V,J]=jordan(D_Eq_at_NE)
OneillYMQ20201118
[V,J]=jordan(D_Eq_at_NE)
for m=1:8;for n=1:8;test(m,n) = V(:,m)'*V(:,n);end;end
%-- 2021/1/15 11:50 --%
ft = price20210115addnoise(1)
figure'
figure
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(ft(4000:end,:));
[v d] = eig(L_matrix)
plot(v(:,1))
plot(v(:,3))
clf
[Lmn Tmn]= from_eigenvector_out_am(v(:,1))
[v d] = jordan(L_matrix)
[U,H] = hermiteForm(L_matrix )
[v d] = eig(L_matrix)
%-- 2021/1/17 8:46 --%
price20210115
ft = price20210115(1)
v=eigenvector_set(:,maxImag_eigenvalueId)
[RPSvalue RPSgroup]= from_eigenvector_out_RPSvalue(v)
v_angle = angle(v)
RPS_index(1) = find(v_angle >= 0 & v_angle < 2*pi/3)
RPS_index1 = find(v_angle >= 0 & v_angle < 2*pi/3)
RPS_index3 = find(v_angle < 0 & v_angle > -2*pi/3)
RPS_index2 = find(ismember(1:Ndim,[RPS_index(1) RPS_index(3)])==0);
RPS_index2 = find(ismember(1:Ndim,[RPS_index1 RPS_index3])==0);
RPS_index2 = find(ismember(1:Ndim,[RPS_index1' RPS_index3'])==0);
Ndim=length(v);
RPS_index2 = find(ismember(1:Ndim,[RPS_index1' RPS_index3'])==0);
Ndim=length(v);
RPSgroup = zeros(Ndim,3);
v_angle = angle(v);
RPS_index1 = find(v_angle >= 0 & v_angle < 2*pi/3)
RPS_index3 = find(v_angle < 0 & v_angle > -2*pi/3)
RPS_index2 = find(ismember(1:Ndim,[RPS_index1' RPS_index3'])==0);
RPS_index2 = find(ismember(1:Ndim,[RPS_index1' RPS_index3'])==0)';
R = sum(v(RPS_index1))
P = sum(v(RPS_index2))
S = sum(v(RPS_index3))
RPSv(1) = sum(v(RPS_index1))
RPSv(2) = sum(v(RPS_index2))
RPSv(3) = sum(v(RPS_index3))
area_m_n = zeros(3,3)
for m=1:3-1
for n=m+1:3
vm=RPSv(m);
vn=RPSv(n);
area_m_n = [m n abs(vm)*abs(vn)*pi*sin(angle(vm)-angle(vn))];
end
end
area_m_n = zeros(3,3)
for m=1:3-1
for n=m+1:3
vm=RPSv(m);
vn=RPSv(n);
area_m_n(mn,1) = m;
area_m_n(mn,2) =n;
area_m_n(mn,3) = abs(vm)*abs(vn)*pi*sin(angle(vm)-angle(vn)) ;
end
end
area_m_n = zeros(3,3)
mn=1
for m=1:3-1
for n=m+1:3
vm=RPSv(m);
vn=RPSv(n);
area_m_n(mn,1) = m;
area_m_n(mn,2) =n;
area_m_n(mn,3) = abs(vm)*abs(vn)*pi*sin(angle(vm)-angle(vn)) ;
mn=mn+1;
end
end
RPSgroup(1,RPS_index1) = 1;
RPSgroup = zeros(Ndim,3);
RPSgroup(RPS_index1,1) = 1;
Ndim=length(v);
RPSgroup = zeros(Ndim,3);
v_angle = angle(v);
RPS_index1 = find(v_angle >= 0 & v_angle < 2*pi/3);
RPS_index3 = find(v_angle < 0 & v_angle > -2*pi/3);
RPS_index2 = find(ismember(1:Ndim,[RPS_index1' RPS_index3'])==0)';
RPSgroup(RPS_index1,1) = 1;
RPSgroup(RPS_index2,2) = 1;
RPSgroup(RPS_index3,3) = 1;
[RPSvalue RPSgroup]= from_eigenvector_out_RPSvalue(v)
%-- 2021/1/18 15:43 --%
payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ;6 10 10 0 3 ; 10 10 10 4 0 ]
p_1=444/2987;p_2=22/103;p_3=641/2987;p_4=288/2987;p_5=976/2987;
Pn=[];
p=[p_1 p_2 p_3 p_4 p_5];
for t1=1:5
payoff_matrix(:,t1)*p(t1)
Pn=[Pn payoff_matrix(:,t1)*p(t1)];
end
g = [Pn(:,1)./p_1 Pn(:,4)./p_4  sum(Pn(:,[2 3 5])')']
delta=p_2 + p_3 + p_5;
sigma=[0 p_2  p_3 0 p_5]/(p_2 + p_3 + p_5);
t2 = sum(sigma(2)*g(2,:) + sigma(3)*g(3,:) + sigma(5)*g(5,:))
f = [g(1,:); g(4,:); t2]
delta=p_2 + p_3 + p_5;
sigma=[0 p_2  p_3 0 p_5]/(p_2 + p_3 + p_5);
t2 =  (sigma(2)*g(2,:) + sigma(3)*g(3,:) + sigma(5)*g(5,:))
f = [g(1,:); g(4,:); t2]
delta=p_2 + p_3 + p_5;
sigma=[0 p_2  p_3 0 p_5]/(p_2 + p_3 + p_5);
t2 = sum(Pn(:,[2 3 5])')'./delta
g = [Pn(:,1)./p_1 Pn(:,4)./p_4  t2]
t3 =  (sigma(2)*g(2,:) + sigma(3)*g(3,:) + sigma(5)*g(5,:))
f = [g(1,:); g(4,:); t3]
x5_YQM
syms x1 x2 x3 x4 x5 real
payoff_matrix = [ 0, 11, 6.7468
6, 0, 6.9703
6.3113, 7.4018, 5.9637
]
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 ]'
mean_U = [x1 x2 x3 ] * Payoff_vector_field_F
V_Eq_0 = [x1 x2 x3 ]'.*(Payoff_vector_field_F - mean_U);
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3')  ]
x1=444/2987;x2=288/2987;x3=1-x1-x2;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
[eigen_vector1 eigen_value1] = jordan(D_Eq_at_NE)
x5_YQM
%-- 2021/1/19 7:44 --%
x5_pool3
x5_YQM
-0.6674/0.8259
x5_pool3
B(end)
g = [Pn(:,A(1))./p(A(1)) Pn(:,A(2))./p(A(4))  t2];
t2 = sum(Pn(:,B)')'./delta;
x5_pool3
g = [Pn(:,A(1))./p(A(1)) Pn(:,A(2))./p(A(4))  t2];
g = [Pn(:,A(1))./p(A(1)) Pn(:,A(2))./p(A(2))  t2];
x5_pool3
sigma=[0 p(B(1))  p(B(2)) 0 p(B(3))]/delta;
t3 = sigma(2)*g(2,:) + sigma(3)*g(3,:) + sigma(5)*g(5,:);
eigen_vector
x5_pool3
M3
Ne1=[0.148642  0.213596  0.637762];
syms x1 x2 x3 x4 x5 real
payoff_matrix = M3;
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 ]';
mean_U = [x1 x2 x3 ] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 ]'.*(Payoff_vector_field_F - mean_U);
D_Eq_0 = [ diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3')  ]
% x1=444/2987;x2=288/2987;x3=1-x1-x2;
x1=Ne1(1);x2=Ne1(2);x3=1-x1-x2;
D_Eq_at_NE = eval(D_Eq_0);
[eigen_vector eigen_value] = eig(D_Eq_at_NE);
x5_pool3
x5_dynamic3
x5_YQM
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
v=eigen_vector(:,3);
text(real(v),imag(v), [1:5]')
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
v=eigen_vector(:,3);hold on
text(real(v),imag(v), [1:5]')
text(real(v),imag(v), num2str([1:5]'))
text(real(v),imag(v), num2str([1:5]'),'fontsize',15);axis square;
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
v=eigen_vector(:,3);hold on
text(real(v),imag(v), num2str([1:5]'),'fontsize',15);axis square;
clf
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
v=eigen_vector(:,3);hold on
text(real(v)*1.05,imag(v)*1.05, num2str([1:5]'),'fontsize',15);axis square;
abs(v))
abs(v)
%-- 2021/1/21 14:14 --%
read312data
u=raw(1:length(raw(:,1)),:);
read312data
t=length(txt(:,1));
u=raw(1:length(raw(:,1)),:);
v=cell2table(u);
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
aa=find(raw(:,3)=='T001')
aa=find(raw(:,3)='T001')
aa=find(strcmp(raw(:,3),'T001')
aa=find(strcmp(raw(:,3),'T001'))
abedfilesname = 'C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-01\data20201230\202012pool.csv';  %
[num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
%-- 2021/1/21 19:35 --%
Mix5_2000
[A]=bimat(payoffA)
[A]=bimat(payoffA,payoffA')
%-- 2021/2/1 20:35 --%
load('C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-02\pt01.mat')
seizures_data
Arr= seizures_data(1,1);
B = corr(Arr)
Arr= seizures_data(1,1).data;
Arr= seizures_data(1,1).
Arr= transpose(str2num(cell2mat(seizures_data(1,1))))
Arr= transpose((cell2mat(seizures_data(1,1))))
Arr= transpose((cell2mat(seizures_data(1,1))));
B = corr(Arr);
B = corr(Arr');
plot(B)
h=pcolor(data_add1c1l(B));
h=pcolor(B);
colormap(jet(8));
load('F:\Spectrum2020\datasource\oneill15327.mat'); d=data15327(:,4:11);
clear
load('F:\Spectrum2020\datasource\oneill15327.mat'); d=data15327(:,4:11);
timesies2figure8(d);
load('C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-02\pt01.mat')
clear
load('C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-02\pt01.mat')
Arr= transpose((cell2mat(seizures_data(1,1))));
timesies2figure8(Arr');
Arr= transpose((cell2mat(seizures_data(1,1))))';
Cm=mean(Arr);
f =ones(50000,1) * Cm;
Brr=Arr-f;
timesies2figure8(Brr);
t1=mean(Brr)
sum(t1)
Brt = Brr'; Cm2 = mean(Brt);
f2 = ones(75,1) * Cm2;
Crr = (Brt-f2)';
timesies2figure8(Crr);
title('time series');hold on;
plot(dos_x_t);
xlabel('time');ylabel('density');axis square;grid on;
mmax = max(max(Crr))
mmin = min(min(Crr))
Crt = Crr - mmin;
figure; plot(Crt);
Cm3 = sum(Crt')';
f3=[]; for i=1:75; f3 = [f3 Cm3];end
Drr = Crt ./ f3
figure; plot(Drr);
timesies2figure8(Drr);
save('C:\Users\Think\Desktop\2021毕设\EightFig.mat')
%-- 2021/2/2 20:31 --%
load('C:\Users\Think\Desktop\2021毕设\EightFig.mat')
A = m24x24-m24x24';
[v d] = eig(A)
plot(v(:,1))
plot(v(:,2))
plot(v(:,3))
plot(v(:,4))
plot(v(:,1))
plot(v(:,3))
plot(v(:,5))
plot(v(:,7))
diag(d)
absv=abs(v);
plot(absv(:,1))
plot(absv,'DisplayName','absv')
plot(absv(:,1:2:10),'DisplayName','absv')
plot(absv(:,1:2:10),'DisplayName','absv','linewidth',2)
load('C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-02\pt01.mat')
%-- 2021/2/7 7:08 --%
OneillYMQ20201118
[jigen_vector jigen_value] = eig(D_Eq_at_NE)
[jigen_vector jigen_value] = jordan(D_Eq_at_NE)
v=jigen_vector(:,5);[Lmn Tmn]= from_eigenvector_out_am(v)
%-- 2021/2/9 14:48 --%
run('C:\Users\Think\Desktop\2021毕设\HAVOK\FIG01_LORENZ_GEN.m')
run('C:\Users\Think\Desktop\2021毕设\HAVOK\FIG01_LORENZ_PLOT.m')
FIG01_LORENZ_GEN
%-- 2021/2/9 17:45 --%
FIG01_LORENZ_GEN
length(V)
plot(x(:,1))
plot(x(:,2))
plot(x(:,3))
plot(x(:,5))
plot(x(:,10))
plot(x(:,end))
plot(x(:,11))
plot(x(:,end))
plot(x(:,7))
plot(x(:,1))
plot(x(:,2))
plot(x(:,1))
figure
plot(x(:,end))
figure
plot(x(:,2))
figure
plot(x(:,3))
figure; plot(x(:,1:3));
figure; plot(x(1:4000,1:3));
figure; plot(x(1:40000,1:3));
FIG01_LORENZ_GEN
U
plot(U(:,1))
plot(U(:,2))
plot(U(:,3))
plot(U(:,4))
plot(U(:,5))
plot(U(:,6))
plot(U(:,7))
plot(U(:,8))
plot(U(:,9))
plot(U(:,10))
plot(U(:,end))
plot(U(:,99))
plot(U(:,87))
plot(U(:,81))
FIG01_LORENZ_GEN
plot(V(:,1))
plot(V(:,2))
plot(V(:,end))
plot(U,'DisplayName','U')
plot(U(:,1))
plot(U(:,1:2),'DisplayName','U(:,1:2)')
plot(U(:,1:3),'DisplayName','U(:,1:3)')
plot(V,'DisplayName','V')
plot(V(:,end))
plot(V(:,1))
plot(V(:,end))
plot(V(:,end));hold  on; plot(V(:,1))
load('C:\Users\Think\Desktop\2021毕设\博弈实验时间序列-群体.xlsx')
xlsread('C:\Users\Think\Desktop\2021毕设\博弈实验时间序列-群体.xlsx')
a=ans;
xdat = a;
FIG01_LORENZ_GEN
a=xlsread('C:\Users\Think\Desktop\2021毕设\博弈实验时间序列-群体.xlsx');
xdat = a(:,5:8);
xdat = a(:,6:8);
plot(V(:,end));hold  on; plot(V(:,1))
hold  on; plot(xdat(:,1))
[y,t] = lsim(sys,x(L,r),dt*(L-1),x(1,1:r-1));
FIG01_LORENZ_PLOT
load('DATA/FIG01_LORENZ.mat'); %WZJ
FIG01_LORENZ_GEN
save ./DATA/FIG01_LORENZ.mat
clear
FIG01_LORENZ_PLOT
load('C:\Users\Think\Desktop\2021毕设\HAVOK\DATA\Lorenz.mat')
clear
load('C:\Users\Think\Desktop\2021毕设\HAVOK\DATA\sleepeeg.mat')
plot(sleepeeg(1,:))
plot(sleepeeg(2,:))
plot(sleepeeg(1:2,:))
%-- 2021/2/17 9:12 --%
help neteork
net = network(1, 2, [1; 1], [1; 0], [0 0; 1 0], [0 1]);
net.inputs{1}.size = 122; % input size
net.layers{1}.size = 25; % hidden layer size
net.layers{2}.size = 1; % output layer size
net.view;
%-- 2021/2/19 15:49 --%
A=[-21,19,-20;19,-21,20;40,-40,-40];
B=[0,1,2]’;
C=[1,0,2];
D=[0];
sf=ss(A,B,C,D);
g=tf(sf)
A=[-21,19,-20;19,-21,20;40,-40,-40];
B=[0,1,2]';
C=[1,0,2];
D=[0];
sf=ss(A,B,C,D);
g=tf(sf)
h = tf([1 0],[1 2 10])
zeta = 0.25;
w0 = 3;
H = tf(w0^2,[1,2*zeta*w0,w0^2])
stepplot(H)
clear;
Fs = 5;
dt = 1/Fs;
N = 50;
t = dt*(0:N-1);
A = [cos(dt) sin(dt);-sin(dt) cos(dt)];
B = [1-cos(dt);sin(dt)];
C = [-1 0];
D = 1;
H = tf(zeros(1,1,10));
s = tf('s')
for k=1:10,
H(:,:,k) = k/(s^2+s+k);
end
filter
Fofx=inline('x.^2.*cos(a.*x)-b','x','a','b');g = Fofx([pi/3 pi/3.5],4,1)
clear
https://blog.csdn.net/EliminatedAcmer/article/details/80487859
PEd
clear;
f=sym('y+2*x/y^2');a=0;b=2;
h=0.4;
n=(b-a)/h+1;
x=0;
y=1;
szj=[x,y];%数值解
for i=1:n-1
y=y+h*subs(f,{'x','y'},{x,y});%subs，替换函数
x=x+h;
szj=[szj;x,y];
end;
szj;
plot(szj(:,1),szj(:,2))
num_1=[2 14 24];  den_1=[1 5 8 4];
sys_1=tf(num_1,den_1);
[a,b,c,d]=ssdata(sys_1);
[v,j]=jordan(a);
sys_2=ss(va*v,v,c*v,d)
sys_2=ss(v,a*v,v,c*v,d)
sys_2=ss(a*v,v,c*v,d)
num=[0 0 10 10];
den=[1 6 6 10];
gs = tf(num,den)
[a,b,c,d]=ssdata(sys_1);
[v,j]=jordan(a);
[a,b,c,d]=ssdata(sys_1)
gain
%-- 2021/2/21 15:30 --%
run('F:\8x8code\Matlab4Student2018\granger_cause.m')
granger_cause
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1:2000,:);
[F,c_v] = granger_cause(data2(:,1),data2(:,2),0.05,20)
[F,c_v] = granger_cause(data2(:,2),data2(:,1),0.05,20)
plot(data2(:,1:2),'DisplayName','data2(:,1:2)')
[F,c_v] = granger_cause(data2(:,2),data2(:,1),0.005,20)
[F,c_v] = granger_cause(data2(:,2),data2(:,1),0.00005,20)
[F,c_v] = granger_cause(data2(:,1),data2(:,2),0.00005,20)
[F,c_v] = granger_cause(data2(:,1),data2(:,2),0.0000000000005,20)
data2=data1(1:200,:);
[F,c_v] = granger_cause(data2(:,1),data2(:,2),0.0000000000005,20)
[F,c_v] = granger_cause(data2(:,2),data2(:,1),0.0000000000005,20)
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,2),data2(:,1),0.0000000000005,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0000000000005,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=data1(1:2000,:);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0000000000005,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
ZSJ1227TimeSeries
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
data2= data2 - ones(size(data2,1),1)*mean(data2)
data2=t2(1:20:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:10:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:10:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,5);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:1:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,5);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:1:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,50);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:10:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,5);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
plot(t2(1:442,1:4),'DisplayName','t2(1:442,1:4)')
%           t2=[[[x1 x2 x3 x4] + kp(t02i)*c2] 0]
t2=[[(Ne) + kp(t02i)*real(c2)] 0];
tstep = 0.2; %时间步长取0.02
noise_amplit = 0;
for k=0:tstep:4*pi/period_para*10
%                 seed = 0.01*rand(4,1);
%                 seed = seed - mean(seed);
x1 = t2(length(t2(:,1)),1) + eval(V_F(1))*tstep;
x2 = t2(length(t2(:,1)),2) + eval(V_F(2))*tstep;
x3 = t2(length(t2(:,1)),3) + eval(V_F(3))*tstep;
x4 = t2(length(t2(:,1)),4) + eval(V_F(4))*tstep;
x_addnoise = addNoise2timeseries([x1 x2 x3 x4]', noise_amplit);
t2=[t2; x_addnoise' k];  % x_addnoise' replaced by [x1 x2 x3 x4]
end
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
data2=t2(1:1:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,5);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:1:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,25);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
plot(data2(:,1:4),'DisplayName','data2(:,1:4)')
data2=t2(1:10:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,25);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:20:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,25);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data2=t2(1:5:end,:);
data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999; for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
clear
ZSJ1227TimeSeries
data1=load('a_1repdone.csv');  %书洁的仿真结果
data2=data1(1:2000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0000000000005,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1:2000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0000000000005,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
[F,c_v] = granger_cause(data2(:,j),data2(:,j),0.0000000000005,20)
[F,c_v] = granger_cause(data2(:,j),data2(:,j),0.00005,20)
[F,c_v] = granger_cause(data2(:,j),data2(:,j),0.05,20)
[F,c_v] = granger_cause(data2(:,j),data2(:,j),0.05,2)
[F,c_v] = granger_cause(data2(:,j),data2(:,j),0.05,200)
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1000:3000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0005,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1000:3000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0005,2);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1000:3000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0005,10);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1000:3:9000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0005,10);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1:2000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.0005,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1:2000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
V=MF-Mc_v
F
MF
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1:2000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1:2000,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,20);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,10);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,5);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_4repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,5);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_025repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,5);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_025repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_4repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1800:2001,:);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_025repdone.csv');
data2=data1(1800:2001,:);data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_1repdone.csv');
data2=data1(1800:2001,:);data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
data1=load('F:\8x8code\ZSJ20210221\data_zsj\a_4repdone.csv');
data2=data1(1800:2001,:);data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
%-- 2021/2/24 12:14 --%
v_4strategy_logit
ZSJ1227TimeSeries
data2=t2(1800:2001,:);data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
granger_cause
ZSJ1227TimeSeries
data2=t2(1800:2001,:);data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
plot(data2(:,1:4),'DisplayName','data2(:,1:4)')
data2=t2(1:2001,:);data2= data2 - ones(size(data2,1),1)*mean(data2);
MF=ones(4)*999;Mc_v=ones(4)*999;
for j=1:4;for k=1:4;
[F,c_v] = granger_cause(data2(:,j),data2(:,k),0.05,1);
MF(j,k)=F;Mc_v(j,k)=c_v; end;end
MF
V=MF-Mc_v
plot(data2(:,1:4),'DisplayName','data2(:,1:4)')
data2=t2(1800:2001,:);data3= data2 - ones(size(data2,1),1)*mean(data2);
data4 = data3. / ones(size(data3,1),1)*max(data3)
data4 = data3. \ ones(size(data3,1),1)*max(data3)
data4 = data3 .\ ones(size(data3,1),1)*max(data3)
data4 = data3 ./ ones(size(data3,1),1)*max(data3)
data4 = data3 ./ max(data3)
data4 = data3./ max(data3)
data4 = data3. / max(data3)
data4 = data3. / max(data3).
data4 = data3 ./ max(data3).
data4 = data3 ./ max(data3)
1/max(data3).
1/max(data3)
1\max(data3)
max(data3)
1/max(data3)
data4 = data3. / max(data3(:,1:4))
data4 = data3(:,1:4) ./ max(data3(:,1:4))
data4 = max(data3(:,1:4)) .\ data3(:,1:4)
data4 = max(data3(:,1:4)) ./ data3(:,1:4)
max(data3(:,1:4))
1\max(data3(:,1:4))
1/max(data3(:,1:4))
1/max(data3(:,1:4)).
data4 = data3 ./ (ones(size(data3,1),1)*max(data3))
data2=t2(1:2001,:);data3= data2 - ones(size(data2,1),1)*mean(data2);
data4 = data3 ./ (ones(size(data3,1),1)*max(data3))
plot(data4(:,1:4),'DisplayName','data4(:,1:4)')
data2=t2(1:100,:);data3= data2 - ones(size(data2,1),1)*mean(data2);
data4 = data3 ./ (ones(size(data3,1),1)*max(data3))
plot(data4(:,1:4),'DisplayName','data4(:,1:4)')
data2=t2(1:1000,:);data3= data2 - ones(size(data2,1),1)*mean(data2);
data4 = data3 ./ (ones(size(data3,1),1)*max(data3))
plot(data4(:,1:4),'DisplayName','data4(:,1:4)')
(ones(size(data3,1),1)*max(data3))
max(data3)
plot(data3(:,1:4),'DisplayName','data3(:,1:4)')
figure(1)
plot(real(t2(:,[1 2 3 4])),'linewidth',4)
data2=t2(1:1000,:);data3= data2 - ones(size(data2,1),1)*mean(data2);;
plot(data3(1:1000,1:4),'linewidth',4);
ZSJ1227TimeSeries
%-- 2021/2/26 16:31 --%
timesies2figure8
load('F:\cDownload\PokerAnan20200716\8x8abedata\realization_linear_4096.mat')
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(data')
mean(data')
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(data')
%-- 2021/2/27 15:20 --%
syms x1 x2 x3 x4 x5 x6 real
%% 动力学方程
dy=zeros(3,1);  %建立 一个 三维列 向量
dy(1)=10*(-y(1)+y(2));
dy(2)=28*y(1)-y(2)-Y(1)*y(3);
dy(3)=Y(1)*y(2)-8*y(3)/3;
V_F = dy;
syms y
dy=zeros(3,1);  %建立 一个 三维列 向量
dy(1)=10*(-y(1)+y(2));
dy(2)=28*y(1)-y(2)-Y(1)*y(3);
dy(3)=Y(1)*y(2)-8*y(3)/3;
V_F = dy;
syms y
dy=zeros(3,1);  %建立 一个 三维列 向量
dy1=10*(-y(1)+y(2));
dy2=28*y(1)-y(2)-Y(1)*y(3);
dy3=y(1)*y(2)-8*y(3)/3;
syms y
% dy=zeros(3,1);  %建立 一个 三维列 向量
dy1=10*(-y(1)+y(2));
dy2=28*y(1)-y(2)-Y(1)*y(3);
dy3=y(1)*y(2)-8*y(3)/3;
syms y1 y2 y3
% dy=zeros(3,1);  %建立 一个 三维列 向量
dy1=10*(-y1+y2);
dy2=28*y1-y2-y1*y3;
dy3=y1*y2-8*y3/3;
V_F = [dy1;dy2; dy3]
S=solve(dy1,dy2,dy3)
S.y1
Restpoint = [S.y1 S.y2 S.y3]
syms y1 y2 y3
% dy=zeros(3,1);  %建立 一个 三维列 向量
dy1=10*(-y1+y2);
dy2=28*y1-y2-y1*y3;
dy3=y1*y2-8*y3/3;
% V_F = [dy1;dy2; dy3];
S=solve(dy1,dy2,dy3);
Restpoint = [S.y1 S.y2 S.y3]
syms y1 y2 y3
% dy=zeros(3,1);  %建立 一个 三维列 向量
dy1=10*(-y1+y2);
dy2=28*y1-y2-y1*y3;
dy3=y1*y2-8*y3/3;
V_Eq_0 = [dy1;dy2; dy3];
S=solve(dy1,dy2,dy3);
Restpoint = [S.y1 S.y2 S.y3]
D_Eq_0 = [ diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3')]
subs(D_Eq_0, 'y3', Restpoint(2,3))
subs(D_Eq_0, 'y1', Restpoint(2,1), 'y3', Restpoint(2,3))
subs(D_Eq_0, '[y1 y3]', [Restpoint(2,1),  Restpoint(2,3)])
subs(D_Eq_0, '[y1 y3]', [Restpoint(2,1) Restpoint(2,3)])
subs(D_Eq_0, 'y1 y3', [Restpoint(2,1) Restpoint(2,3)])
subs(subs(D_Eq_0, 'y1', Restpoint(2,1)),'y3', Restpoint(2,3)])
subs(subs(D_Eq_0, 'y1', Restpoint(2,1)),'y3', Restpoint(2,3))
subs(subs(subs(D_Eq_0, 'y1', Restpoint(2,1)),'y2', Restpoint(2,2)),'y3', Restpoint(2,3))
A = subs(subs(subs(D_Eq_0, 'y1', Restpoint(2,1)),'y2', Restpoint(2,2)),'y3', Restpoint(2,3))
[V D]=eig(A)
[V D]=eig(eval(A))
plot(V(:,2))
eval(Restpoint)
s2=1; A = subs(subs(subs(D_Eq_0, 'y1', Restpoint(s2,1)),'y2', Restpoint(s2,2)),'y3', Restpoint(s2,3))
[V D]=eig(eval(A))
s2=2; A = subs(subs(subs(D_Eq_0, 'y1', Restpoint(s2,1)),'y2', Restpoint(s2,2)),'y3', Restpoint(s2,3))
s2=2; A = subs(subs(subs(D_Eq_0, 'y1', Restpoint(s2,1)),'y2', Restpoint(s2,2)),'y3', Restpoint(s2,3));[V D]=eig(eval(A))
s2=3; A = subs(subs(subs(D_Eq_0, 'y1', Restpoint(s2,1)),'y2', Restpoint(s2,2)),'y3', Restpoint(s2,3));[V D]=eig(eval(A))
%-- 2021/2/28 15:22 --%
x5_dynamic3
uiopen('F:\cDownload\PokerAnan20200716\8x8abedata\test1_LorezSim.slx',1)
upIndex =-1:1
a=4^upIndex;
syms x1 x2 x3 x4  real
%% 支付矩阵
% payoff_matrix = [0 -1 0 1; 1 0 -1 0; 0 1 0 -1; -1 0 1 0]
payoff_matrix = [0 0 0 a;
1 0 0 0;
0 1 0 0;
0 0 1 0];
%% 空间各个点的支付向量
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4]'
%% 各点的支付均值
mean_U = [x1 x2 x3 x4] * Payoff_vector_field_F
%% 动力学方程
V_F = [x1 x2 x3 x4]'.*(Payoff_vector_field_F - mean_U);
upIndex =-1
a=4^upIndex;
syms x1 x2 x3 x4  real
%% 支付矩阵
% payoff_matrix = [0 -1 0 1; 1 0 -1 0; 0 1 0 -1; -1 0 1 0]
payoff_matrix = [0 0 0 a;
1 0 0 0;
0 1 0 0;
0 0 1 0];
%% 空间各个点的支付向量
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4]'
%% 各点的支付均值
mean_U = [x1 x2 x3 x4] * Payoff_vector_field_F
%% 动力学方程
V_F = [x1 x2 x3 x4]'.*(Payoff_vector_field_F - mean_U)
%-- 2021/3/1 15:10 --%
ZSJ0224TimeSeries
%-- 2021/3/3 13:05 --%
N_cycle0315_Price6
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price6-20201212.mat')
diag(eigen_value)
clear
load('F:\cDownload\PokerAnan20200716\8x8abedata\Price6\Price8-20201212.mat')
%-- 2021/3/10 19:26 --%
A=load('C:\Users\Think\Desktop\prob2.csv')
histogram(A)
histfit(A)
histfit(A)
histfit(A)
%-- 2021/3/26 1:17 --%
OneillYMQ20201118
filenameYQM='E:\AbedData\01\Y202103\8colA.csv'; [Yret1 Yret3] = in_8colExp_out_am(filenameYQM);
filenameYQM='E:\AbedData\01\Y202103\8colB.csv'; [Yret1 Yret3] = in_8colExp_out_am(filenameYQM);
%-- 2021/3/28 14:41 --%
OneillYMQ20201118
%-- 2021/4/6 13:46 --%
OneillYMQ20201118
%-- 2021/4/7 19:36 --%
ft = price20210115(1)
price20210115(2)
price20210115
price20210115(2)
ft = price20210115(2)
%-- 2021/4/11 17:16 --%
a=load('F:\CFH2014\dan2021\dan2021Fastall.csv')
a=xlsread('F:\CFH2014\dan2021\dan2021Fastall.csv')
a1=a(find(a(:,14)==1),:);
a2=a(find(a(:,14)==2),:);
unique(a1(:,10)
unique(a1(:,10))
sessionida2 = unique(a2(:,10))
sessionida1 = unique(a1(:,10))
tmp=[sessionida1 sessionida2]
tmp(:,1) - tmp(:,2)
hist(:,4:7)
hist(a1(:,4:7))
figure; hist(a2(:,4:7))
%-- 2021/4/11 19:59 --%
a=xlsread('F:\CFH2014\dan2021\dan2021Fastall.csv')
a1=a(find(a(:,14)==1),:);
a2=a(find(a(:,14)==2),:);
b1 = []; for j=1:6;b1=[b1;a1(:,3+j)];end
b2 = []; for j=1:6;b2=[b2;a2(:,3+j)];end
b1 = []; for j=1:6;b1=[b1;round(a1(:,3+j),1)];end
b2 = []; for j=1:6;b2=[b2;round(a2(:,3+j),1)];end
for k=2:length(b1(:,1))-1; b1(k,2)=b1(k,1)-b1(k-1,1);b1(k,3)=b1(k+1,1)-b1(k,1);end
b1(:,4) = (b1(:,2) + b1(:,3))/2;
mean(b1(:,4))
hist(b1(:,4),40)
for k=2:length(b2(:,1))-1; b2(k,2)=b2(k,1)-b2(k-1,1);b2(k,3)=b2(k+1,1)-b2(k,1);end
b2(:,4) = (b2(:,2) + b2(:,3))/2; figure; hist(b2(:,4),40)
v1=[]; for j=1:20; t1=b1(find(b1(:,1)==0.1*j),:); v1=[v1; 0.1*j mean(t1)];end
v1=[]; for j=1:20; t1=b1(find(b1(:,1)==0.1*j),:); v1=[v1; 0.1*j mean(t1(:,4))];end
v1=[]; for j=1:20; t1=b1(find(b1(:,1)==0.1*j),:); v1=[v1; 0.1*j count(t1(:,4))];end
v1=[]; for j=1:20; t1=b1(find(b1(:,1)==0.1*j),:); v1=[v1; 0.1*j length(t1(:,4))];end
sum(v1(:,2))
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*j)<0.001),:); v1=[v1; 0.1*j length(t1(:,4))];end
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*j)<0.001),:); v1=[v1; 0.1*j mean(t1(:,4))];end
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*j)<0.001),:); v1=[v1; 0.1*j length(t1(:,4))];end
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*j)<0.001),:); v1=[v1; 0.1*j length(t1(:,4)) sum(t1(:,4))];end
sum(v1(:,2))
v1=[]; for j=1:120; t1=b1(find(abs(b1(:,1)-0.1*j)<0.001),:); v1=[v1; 0.1*j length(t1(:,4)) sum(t1(:,4))];end
sum(v1(:,2))
v1=[]; for j=1:200; t1=b1(find(abs(b1(:,1)-0.01*j)<0.001),:); v1=[v1; 0.1*j length(t1(:,4)) sum(t1(:,4))];end
v1=[]; for j=1:200; t1=b1(find(abs(b1(:,1)-0.01*j)<0.001),:); v1=[v1; 0.01*j length(t1(:,4)) sum(t1(:,4))];end
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*j)<0.001),:); v1=[v1; 0.1*j length(t1(:,4)) sum(t1(:,4))];end
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*j)<0.01),:); v1=[v1; 0.1*j length(t1(:,4)) sum(t1(:,4))];end
sum(v1(:,2))
max(b1(:,1))
min(b1(:,1))
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*(j-1))<0.01),:); v1=[v1; 0.1*(j-1) length(t1(:,4)) sum(t1(:,4))];end
sum(v1(:,2))
mean(b(:,1))
mean(b1(:,1))
find(empty(a(:,4)))
find(isnull(a(:,4)))
find(~empty(a(:,4)))
v1=[]; for j=1:20; t1=b1(find(abs(b1(:,1)-0.1*(j-1))<0.01),:); v1=[v1; 0.1*(j-1) length(t1(:,4)) sum(t1(:,4))];end
sum(v1(:,2))
%-- 2021/4/12 13:14 --%
ZSJ0224TimeSeries
fffffffffffffffffff
eval(D_V_F)
%-- 2021/4/12 20:11 --%
ft = price20210115(1)
price20210115
ft = price20210115(1)
price20210115
ft = price20210115(1)
%-- 2021/4/13 14:00 --%
ShujieLandscape1228
latex(eigen_vector)
latex(eigen_value)
latex(simplify(eigen_vector))
%-- 2021/4/16 11:17 --%
ShujieLandscape1228
latex(Ne)
%-- 2021/4/17 23:24 --%
ShujieLandscape1228
Ne*eigen_vector(:,1)
simplify(Ne*eigen_vector(:,1))
simplify(Ne*eigen_vector(:,4))
Ne*eigen_vector(:,1)
Ne*eigen_vector(:,1)'
eigen_vector(:,1)'
eigen_vector(:,1)
eigen_vector(:,2)'*eigen_vector(:,1)
eigen_vector(:,3)'*eigen_vector(:,1)
[V,J] = jordan(eval(D_V_F))
simplify(Ne*V(:,4))
simplify(Ne*V(:,3))
simplify(Ne*V(:,2))
simplify(Ne*V(:,1))
V(:,1)
simplify(V(:,1))
a=4;simplify(Ne*eigen_vector(:,4))
a=4;eval(simplify(Ne*eigen_vector(:,4)))
a=4;eval(simplify(V(:,4)))
a=1;eval(simplify(V(:,4)))
a=1/4;eval(simplify(V(:,4)))
a=1/4;eval(abs(V(:,4)))
a=4;eval(abs(V(:,4)))
a=4;eval(abs(eigen_vector(:,4)))
a=1/4;eval(abs(eigen_vector(:,4)))
%-- 2021/4/18 15:32 --%
ShujieLandscape1228
a=1/4;eval(abs(eigen_vector(:,4)))
a=4;eval(abs(eigen_vector(:,4)))
a=4;p0=eval(abs(eigen_vector(:,4)));p1=p0/sum(p0)
eval(Ne)
eval(Ne)'
a=4;p0=eval(abs(eigen_vector(:,4)));p1=[p0/sum(p0) eval(Ne)']
mean_U
payoff_matrix
payoff_matrix * p1(:,1)
eval(payoff_matrix * p1(:,1))
p1
;
a=4;p0=eval(abs(eigen_vector(:,4)));p1=[p0/sum(p0) eval(Ne)']
eval(payoff_matrix * p1(:,1)) eval(payoff_matrix * p1(:,2))
[eval(payoff_matrix * p1(:,1)) eval(payoff_matrix * p1(:,2))]
c=[eval(payoff_matrix * p1(:,1)) eval(payoff_matrix * p1(:,2))]
mean(c)
a=1/4;p0=eval(abs(eigen_vector(:,4)));p1=p0/sum(p0)
[eval(payoff_matrix * p1(:,1)) eval(payoff_matrix * p1(:,2))]
a=1/4;p0=eval(abs(eigen_vector(:,4)));p1=[p0/sum(p0) eval(Ne)']
c=[eval(payoff_matrix * p1(:,1)) eval(payoff_matrix * p1(:,2))]
%-- 2021/4/19 13:43 --%
ShujieLandscape1228
eval(D_V_F)
a=3; matrix = eval(D_V_F)
a=3; matrix = vap(D_V_F)
a=3; matrix = eval(D_V_F)
a=3; matrix = vpa(D_V_F)
a=3; matrix = vpa(eval(D_V_F))
a=3; matrix = (eval(D_V_F))
a=3; matrix = (eval(D_V_F)); vpa(matrix)
a=3.000; matrix = (eval(D_V_F)); vpa(matrix)
a=3.000; matrix = (eval(D_V_F)); eval(matrix)
a=3; matrix = (eval(D_V_F)); eval(matrix)
matrix
a=3; matrix = (eval(D_V_F)); eval(matrix)
%-- 2021/4/20 14:19 --%
ShujieLandscape1228
a=1/4;p0=eval(abs(eigen_vector(:,4)));p1=[p0/sum(p0) eval(Ne)']
c=[eval(payoff_matrix * p1(:,1)) eval(payoff_matrix * p1(:,2))]
ShujieLandscape1228
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
ShujieLandscape1228
plot(t2(:,end))
area(t2(:,end))
plot(t2(:,[1,end]),'DisplayName','t2(:,[1,end])')
plot(t2(:,[3,end]),'DisplayName','t2(:,[3,end])')
ylim([0.13 0.15])
plot(t2(:,[1:3,end]),'DisplayName','t2(:,[1:3,end])')
ylim([0.13 0.15])
mean(t2(:,[1:3,end]))
ev=mean(t2(:,[1:3,end]))
ev-ev(4)
ev=mean(t2(:,[1:end]))
ev-1/7
mean(t2(:,[1:3,end])) -1/7
ShujieLandscape1228
a=1/4;p0=eval(abs(eigen_vector(:,4)));p1=[p0/sum(p0) eval(Ne)']
p0=eval(abs(eigen_vector(:,4)));p1=[p0/sum(p0) eval(Ne)']
p0= (abs(eigen_vector(:,4)));p1=[p0/sum(p0)  (Ne)']
c=[ (payoff_matrix * p1(:,1))  (payoff_matrix * p1(:,2))]
ev=mean(t2(:,[1:3,end]))
ev-ev(4)
ev=mean(t2(:,[1:3,end]))
ev-ev(4)
ShujieLandscape1228
plot(t2(:,1:3),'DisplayName','t2(:,1:3)')
abs(c2)
v = abs(c2)
v/sum(v)
Ne
%-- 2021/4/21 13:40 --%
eigenvector2cycle
%-- 2021/4/24 9:23 --%
ShujieLandscape1228
run('F:\cDownload\PokerAnan20200716\8x8abedata\from_N_colExp_out_am.m')
v=eigen_vector(:,2);[Lmn Tmn]= from_eigenvector_out_am(v)
ShujieLandscape1228
v=eigen_vector(:,2);[Lmn Tmn]= from_eigenvector_out_am(v)
ShujieLandscape1228
v=eigen_vector(:,2);[Lmn Tmn]= from_eigenvector_out_am(v)
ShujieLandscape1228
v=eigen_vector(:,2);[Lmn Tmn]= from_eigenvector_out_am(v)
ShujieLandscape1228
v=eigen_vector(:,2);[Lmn Tmn]= from_eigenvector_out_am(v)
ShujieLandscape1228
v=eigen_vector(:,2);[Lmn Tmn]= from_eigenvector_out_am(v)
%-- 2021/4/30 16:35 --%
SchodingerEquation20210430
plot(D)
figure
plot(x_t,V(:,2),'k--','linewidth',2)
plot(x_t,U0,'b','linewidth',2);
plot(x_t,V(:,2),'k--','linewidth',2)
hold on;plot(x_t,V(:,1),'r--','linewidth',2)
figure
plot(x_t,U0,'b','linewidth',2);
%-- 2021/5/1 8:17 --%
YQM20200830
MultiLRYQM
see5strategyPatternYQM
abedfilesname = 'F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\SJ6_a44_4.csv';
run('F:\cDownload\PokerAnan20200716\8x8abedata\YQM-ONeill-Simulation\getData.m')
clear
getData
YQM20200830
[eigen_vector eigen_value] = eig(D_Eq_0)
eigen_vector
YQM20200830
af=[]; for k=1:8;v=eigen_vector(:,k);ae=[];for m=1:7;for n=m+1:8;
vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];
end;end;af=[af ae];end;
evt=round(af,4); evl=round(diag(eigen_value),4);
af=[]; for k=1:6;v=eigen_vector(:,k);ae=[];for m=1:5;for n=m+1:7;
vm=v(m);vn=v(n);[area_i_j, Ei, Ej] = lissa2Complex2areaPlot(vm,vn);ae=[ae; area_i_j];
end;end;af=[af ae];end;
evt=round(af,4); evl=round(diag(eigen_value),4);
YQM20200830
[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
see5strategyPatternYQM_a
plotmatrix(a)
plotmatrix(A)
Strategy5
uiopen('F:\cDownload\PokerAnan20200716\8x8abedata\Lissajous(1,1)_YQM20201127_evtID=1.fig',1)
run('F:\cDownload\PokerAnan20200716\8x8abedata\Yret3_fix.m')
path='F:\cDownload\PokerAnan20200716\8x8abedata'
setpath
cd('F:\cDownload\PokerAnan20200716\8x8abedata')
plotfilesaveTo =strcat( figureFilePath, fig_name_saved,'1112.fig')
DesignONeill(1)
Oneill_replicator
load('eigenVV8x8.mat')
Oneill_replicator
NEatK = 0
for k=1:length(A(:,1))
if prod(A(k,1:8)) > 0
NEatK = k
end
end
if NEatK==0
'Wrong'
end
Oneill_replicator
clear
Oneill_replicator
L28x8=[]
for k=1:8
v = eigen_vector(:,k)
[Lmn Tmn]= from_eigenvector_out_am(v)
L28x8(:,k)=Lmn;
end
A(INs,:)'
A(INs,:)
NEat = A(INs,:)
Oneill_replicator
for k=1:8
v = eigen_vector(:,k)
[Lmn Tmn]= from_eigenvector_out_am(v)
L28x8(:,k)=Lmn;
end
Oneill_replicator
for k=1:8
v = eigen_vector(:,k)
[Lmn Tmn]= from_eigenvector_out_am(v)
L28x8(:,k)=Lmn;
end
plot(eigen_vector(:,3))
plot(eigen_vector(:,5))
plot(eigen_vector(:,7))
plot(eigen_vector(:,end))
clear
Oneill_replicator
V_Eq_0
V_Eq_0
Oneill_replicator
for k=1:8
v = eigen_vector(:,k)
[Lmn Tmn]= from_eigenvector_out_am(v)
L28x8(:,k)=Lmn;
end
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
for k=1:8
v = eigen_vector(:,k)
[Lmn Tmn]= from_eigenvector_out_am(v)
L28x8(:,k)=Lmn;
end
Oneill_replicator
for k=1:8
v = eigen_vector(:,k)
[Lmn Tmn]= from_eigenvector_out_am(v)
L28x8(:,k)=Lmn;
end
v_4strategy_two_by_two
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
clear
syms x1 x2 x3 x4 y1 y2 y3 y4 real
payoff_matrix_1 = [1 -1 -1 -1; -1 -1 1 1; -1 1 -1 1; -1 1 1 -1]
payoff_matrix_2 = -payoff_matrix_1
Payoff_vector_field_F_1 = payoff_matrix_1 *[y1 y2 y3 y4]'
Payoff_vector_field_F_2 = payoff_matrix_2 *[x1 x2 x3 x4]'
lamda=0.05
V_Eq_1a =(exp(lamda*Payoff_vector_field_F_1(1))/(exp(lamda*Payoff_vector_field_F_1(1))+exp(lamda*Payoff_vector_field_F_1(2))+exp(lamda*Payoff_vector_field_F_1(3))+exp(lamda*Payoff_vector_field_F_1(4)))-x1);
V_Eq_1b =(exp(lamda*Payoff_vector_field_F_1(2))/(exp(lamda*Payoff_vector_field_F_1(1))+exp(lamda*Payoff_vector_field_F_1(2))+exp(lamda*Payoff_vector_field_F_1(3))+exp(lamda*Payoff_vector_field_F_1(4)))-x2);
V_Eq_1c =(exp(lamda*Payoff_vector_field_F_1(3))/(exp(lamda*Payoff_vector_field_F_1(1))+exp(lamda*Payoff_vector_field_F_1(2))+exp(lamda*Payoff_vector_field_F_1(3))+exp(lamda*Payoff_vector_field_F_1(4)))-x3);
V_Eq_1d =(exp(lamda*Payoff_vector_field_F_1(4))/(exp(lamda*Payoff_vector_field_F_1(1))+exp(lamda*Payoff_vector_field_F_1(2))+exp(lamda*Payoff_vector_field_F_1(3))+exp(lamda*Payoff_vector_field_F_1(4)))-x4);
V_Eq_2a =(exp(lamda*Payoff_vector_field_F_2(1))/(exp(lamda*Payoff_vector_field_F_2(1))+exp(lamda*Payoff_vector_field_F_2(2))+exp(lamda*Payoff_vector_field_F_2(3))+exp(lamda*Payoff_vector_field_F_2(4)))-y1);
V_Eq_2b =(exp(lamda*Payoff_vector_field_F_2(2))/(exp(lamda*Payoff_vector_field_F_2(1))+exp(lamda*Payoff_vector_field_F_2(2))+exp(lamda*Payoff_vector_field_F_2(3))+exp(lamda*Payoff_vector_field_F_2(4)))-y2);
V_Eq_2c =(exp(lamda*Payoff_vector_field_F_2(3))/(exp(lamda*Payoff_vector_field_F_2(1))+exp(lamda*Payoff_vector_field_F_2(2))+exp(lamda*Payoff_vector_field_F_2(3))+exp(lamda*Payoff_vector_field_F_2(4)))-y3);
V_Eq_2d =(exp(lamda*Payoff_vector_field_F_2(4))/(exp(lamda*Payoff_vector_field_F_2(1))+exp(lamda*Payoff_vector_field_F_2(2))+exp(lamda*Payoff_vector_field_F_2(3))+exp(lamda*Payoff_vector_field_F_2(4)))-y4);
V_Eq_1=[V_Eq_1a;V_Eq_1b;V_Eq_1c;V_Eq_1d]
V_Eq_2=[V_Eq_2a;V_Eq_2b;V_Eq_2c;V_Eq_2d]
V_Eq_0 = [V_Eq_1; V_Eq_2]
D_Eq_0 = [diff(V_Eq_0,'x1') diff(V_Eq_0,'x2') diff(V_Eq_0,'x3') diff(V_Eq_0,'x4') ...
diff(V_Eq_0,'y1') diff(V_Eq_0,'y2') diff(V_Eq_0,'y3') diff(V_Eq_0,'y4') ...
]
for i = 1:5
'column-i of the charactor matrix)'
D_Eq_0(:,i)
end
S=solve(V_Eq_0);
[S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]
%下面将纳什均衡填入
A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4]);
x1=A(length(A(:,1)),1); x2=A(length(A(:,1)),2);
x3=A(length(A(:,1)),3); x4=A(length(A(:,1)),4);
y1=A(length(A(:,1)),5); y2=A(length(A(:,1)),6);
y3=A(length(A(:,1)),7); y4=A(length(A(:,1)),8);
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
for k=1:8
v = eigen_vector(:,k)
[Lmn Tmn]= from_eigenvector_out_am(v)
L28x8(:,k)=Lmn;
end
clear
v_4strategy_2x2_logit
L=[];for k=0:10; r=0.001*2^k;L =[L; v_4strategy_2x2_logit(r)];end
%-- 2021/5/1 23:32 --%
L=[];for k=0:2:10; r=0.001*2^k;L =[L; v_4strategy_2x2_logit(r)];end
L=[];for k=0:5:40; r=0.001*2^k;L =[L; v_4strategy_2x2_logit(r)];end
L1=[];L2=[];L3=[];L4=[];for k=0:4:20; r=0.001*2^k;[L28x8 D_Eq_at_NE eigen_vector eigen_value] = v_4strategy_2x2_logit(r);
L1 =[L1; L28x8];L2 =[L2; D_Eq_at_NE];L3 =[L3; eigen_vector];L4 =[L4; eigen_value];end
clean
clear
L1=[];L2=[];L3=[];L4=[];L5=[];for k=0:4:20; r=0.001*2^k;[L28x8 D_Eq_at_NE eigen_vector eigen_value NEat] = v_4strategy_2x2_logit(r);
L1 =[L1; L28x8];L2 =[L2; D_Eq_at_NE];L3 =[L3; eigen_vector];L4 =[L4; eigen_value];L5 =[L5; NEat];end
plot(L5,'DisplayName','L5')
plot(L1(:,1))
plot(L3(:,3))
plot(L3(:,4))
plot(L3(:,1))
plot(L3(:,2))
plot(L3(:,3))
plot(L3(:,4))
plot(L3(:,end))
plot(L3(:,6))
plot(L1(:,3))
plot(L1(:,4))
plot(L1(:,5))
plot(L1(:,6))
plot(L1(:,7))
plot(L1(:,end))
plot(L1(:,3))
plot(L1(:,6))
%-- 2021/5/2 12:00 --%
L1=[];L2=[];L3=[];L4=[];L5=[];for k=0:4:20; r=0.001*2^k;[L28x8 D_Eq_at_NE eigen_vector eigen_value NEat] = v_4strategy_2x2_logit(r);
L1 =[L1 L28x8];L2 =[L2; D_Eq_at_NE];L3 =[L3; eigen_vector];L4 =[L4; eigen_value];L5 =[L5; NEat];end
tL1=L1(:,1:8:end);
L1=[];L2=[];L3=[];L4=[];L5=[];for k=0:1:20; r=0.001*2^k;[L28x8 D_Eq_at_NE eigen_vector eigen_value NEat] = v_4strategy_2x2_logit(r);
L1 =[L1 L28x8];L2 =[L2; D_Eq_at_NE];L3 =[L3; eigen_vector];L4 =[L4; eigen_value];L5 =[L5; NEat];end
plot(L5,'DisplayName','L5')
plot(L5,'DisplayName','L5','linewidth',3)
tL1=L1(:,1:8:end);
plot(tL1(4,:))
plot(L4(:,1))
plot(L4(:,3))
plot(L4(:,4))
plot(L4(:,5))
plot(L4(:,6))
plot(L4(:,7))
plot(L4(:,end))
plot(L4(:,5:end),'DisplayName','L4(:,5:end)')
plot(L4(:,1:2),'DisplayName','L4(:,1:2)')
plot(L4(:,3:4),'DisplayName','L4(:,3:4)')
0.001*2^20
Oneill_replicator
[L28x8 D_Eq_at_NE eigen_vector eigen_value NEat] = v_4strategy_2x2_logit(1048)
[L28x8 D_Eq_at_NE eigen_vector eigen_value NEat] = v_4strategy_2x2_logit(10)
[L28x8 D_Eq_at_NE eigen_vector eigen_value NEat] = v_4strategy_2x2_logit(1)
ShujieLandscape1228
Oneill_replicator
eval(mean_U_1)
eval(mean_U_2)
Oneill_replicator
eval(mean_U_1)
shujie4logit
S=solve(V_Eq_1);
[S.x1 S.x2 S.x3 S.x4 ];
A=eval([S.x1 S.x2 S.x3 S.x4 ])
clear
shujie4logit
A=eval([S.x1 S.x2 S.x3 S.x4 ])
shujie4logit
eval(mean_U)
shujie4logit
size(A,1)
shujie4logit
eval(mean_U)
shujie4logit
uiopen('F:\cDownload\DT20191219\2019讲义\eigenvectorphase1.fig',1)
clear
%-- 2021/5/2 22:18 --%
shujie4logit
A(INs,:) = 1/4
D_Eq_at_NE = eval(D_V_F)
[eigen_vector eigen_value] = jordan(D_Eq_at_NE)
x1=a/(3*a + 1);x2=a/(3*a + 1);;x3=a/(3*a + 1);1/4;x4=1/(3*a + 1);
Ne = [x1 x2 x3 x4];
[eigen_vector eigen_value] = eig(eval(D_V_F))
%-- 2021/5/5 15:03 --%
diag(eigen_value)
angle(diag(eigen_value))
angle(diag(eigen_value))/pi
tmp=diag(eigen_value); [ abs(tmp) angle(tmp)/pi]
L1=[];L2=[];L3=[];L4=[];L5=[];L6=[];for k=-10:2:10; r=2^k;[L28x8 D_Eq_at_NE eigen_vector eigen_value NEat] = v_4strategy_2x2_logit(r);
L1 =[L1 L28x8];L2 =[L2; D_Eq_at_NE];L3 =[L3; eigen_vector];L4 =[L4; eigen_value];L5 =[L5; NEat];end
plot(L4(:,1))
plot(L4(:,3))
t=round(L3,4)
%-- 2021/5/5 22:08 --%
run('D:\MATLAB\R2016a\bin\V2_5005_Exp.m')
%-- 2021/5/6 14:23 --%
v_5strategy
eigen_vector(:,1)'*eigen_vector(:,1)
eigen_vector(:,1)'* payoff_matrix_1*eigen_vector(:,1)
payoff_matrix_1
payoff_matrix_1*eigen_vector(:,1)
ta=payoff_matrix_1*eigen_vector(:,1)
tb=eigen_vector(:,1)'*ta
sum(eigen_vector(:,1))
tc = sum(eigen_vector(:,1))
td = (eigen_vector(:,1)'* payoff_matrix_1*eigen_vector(:,1)) / tc^2
v_5strategy
x1=0.2; x2=x1;x3=x1;x4=x1;x5=x1;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
tc = sum(eigen_vector(:,1))
td = (eigen_vector(:,1)'* payoff_matrix_1*eigen_vector(:,1)) / tc^2
eigen_vector(:,1)
eigen_vector(:,1)/tc
v_5strategy
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
tc = sum(eigen_vector(:,3))
td = (eigen_vector(:,3)'* payoff_matrix_1*eigen_vector(:,3))
td = (eigen_vector(:,4)'* payoff_matrix_1*eigen_vector(:,3))
td = (eigen_vector(:,3)'* D_Eq_at_NE *eigen_vector(:,3))
td = (eigen_vector(:,3)'* eigen_vector(:,3))
td = (eigen_vector(:,4)'* eigen_vector(:,3))
td = (eigen_vector(:,4)'* eigen_vector(:,4))
td = (eigen_vector(:,1)'* D_Eq_at_NE *eigen_vector(:,1))
td = (eigen_vector(:,1)'* payoff_matrix_1*eigen_vector(:,1))
td = (eigen_vector(:,1)'* payoff_matrix_1*eigen_vector(:,1))/tc
tc = sum(eigen_vector(:,1))
td = (eigen_vector(:,1)'* payoff_matrix_1*eigen_vector(:,1))/tc
td = (eigen_vector(:,1)'* payoff_matrix_1*eigen_vector(:,1))/tc^2
tc = sum(abs(eigen_vector(:,3)))
td = (eigen_vector(:,3)'* payoff_matrix_1*eigen_vector(:,3))/tc^2
td = (eigen_vector(:,3)'* payoff_matrix_1*eigen_vector(:,3))/tc
td = abs((eigen_vector(:,3))'* payoff_matrix_1*abs(eigen_vector(:,3))
abs(eigen_vector(:,3))
abs((eigen_vector(:,3))'
td = abs(eigen_vector(:,3))'* payoff_matrix_1*abs(eigen_vector(:,3))
tc = sum(eigen_vector(:,1))
td/tc^2
td = abs(eigen_vector(:,3)')* payoff_matrix_1*abs(eigen_vector(:,3))
td = abs(eigen_vector(:,3)')*D_Eq_at_NE *abs(eigen_vector(:,3))
td =  (eigen_vector(:,3)')*D_Eq_at_NE * (eigen_vector(:,3))
td =  (eigen_vector(:,3)')*payoff_matrix_1 * (eigen_vector(:,3))
3.6551/0.6674
2.3514/0.4294
td1 =  (eigen_vector(:,3)')*D_Eq_at_NE * (eigen_vector(:,3))
td2 =  (eigen_vector(:,3)')*payoff_matrix_1 * (eigen_vector(:,3))
td2/td1
zz=td2/td1
[eigen_vectorP eigen_valueP] = eig(payoff_matrix_1)
[diag(eigen_value); diag(eigen_value)P]
[diag(eigen_value); diag(eigen_valueP)]
[diag(eigen_value) diag(eigen_valueP)]
abs([diag(eigen_value) diag(eigen_valueP)])
[eigen_vector eigen_vectorP]
vd=eigen_vector; vp=eigen_vectorP;
vd(3)'*D_Eq_at_NE*vd(3) - eigen_value(3,3)
vd(:,3)'*D_Eq_at_NE*vd(:,3) - eigen_value(3,3)
vd(:,3)'*payoff_matrix_1*vd(:,3) / eigen_value(3,3)
vd(:,1)'*payoff_matrix_1*vd(:,1) / eigen_value(1,1)
vd(:,2)'*payoff_matrix_1*vd(:,2) / eigen_value(2,2)
vd(:,4)'*payoff_matrix_1*vd(:,4) / eigen_value(4,4)
vd(:,5)'*payoff_matrix_1*vd(:,5) / eigen_value(5,5)
mean(payoff_matrix_1)
mean(payoff_matrix_1')
vp(:,3)'*D_Eq_at_NE *vp(:,3) / eigen_value(3,3)
vp(:,3)'*D_Eq_at_NE *vp(:,3) / eigen_valueP(3,3)
vp(:,3)'*D_Eq_at_NE *vp(:,3) / eigen_valueP(4,4)
vp(:,4)'*D_Eq_at_NE *vp(:,4) / eigen_valueP(4,4)
v_5strategy
[abs(eigen_vectorP(:,3)) angle(eigen_vectorP(:,3))]
[abs(eigen_vector(:,3)) angle(eigen_vector(:,3))]
[angle(eigen_vectorP(:,3)) angle(eigen_vector(:,3))]
[angle(eigen_vectorP(:,3))-angle(eigen_vector(:,3))]
%-- 2021/5/6 19:11 --%
v_5strategy
v = eigen_vector(:,3);[Lmn Tmn]= from_eigenvector_out_am(v)
[eigen_vector eigen_value] = eig(payoff_matrix_1)
v = eigen_vector(:,3);[Lmn Tmn]= from_eigenvector_out_am(v)
%-- 2021/5/6 22:45 --%
ShujieLandscape1228
[eigen_vectorP eigen_valueP] = eig(eval(payoff_matrix))
[eigen_vectorP eigen_valueP] = eig((payoff_matrix))
vd=eigen_vector; vp=eigen_vectorP;
vd(:,3)'*D_Eq_at_NE*vd(:,3) - eigen_value(3,3)
vd(:,3)'*eval(D_V_F)*vd(:,3) - eigen_value(3,3)
vd(:,3)'*eval(D_V_F)*vd(:,3)
vd(:,3)'*payoff_matrix*vd(:,3)
vp(:,3)'*payoff_matrix*vp(:,3)
0.730296743340221/0.365148371670111
0.730296743340221/0.516397779494322
ans^2
0.730296743340221/0.235282268592295
ans^2
0.730296743340221/0.258198889747161
(0.730296743340221/0.258198889747161)^2
v = eigen_vector(:,3);[Lmn Tmn]= from_eigenvector_out_am(v)
v = eigen_vectorP(:,3);[Lmn Tmn]= from_eigenvector_out_am(v)
v = eigen_vector(:,2);[Lmn Tmn]= from_eigenvector_out_am(v)
v = eigen_vectorP(:,3);[LmnP]= from_eigenvector_out_am(v)
v = eigen_vector(:,2);[LmnR]= from_eigenvector_out_am(v)
[LmnP LmnR]
eigen_valueP
payoff_matrix
exp(0.0001*payoff_matrix)
exp(payoff_matrix)
[v d]=eig(payoff_matrix)
[v d]=eig(exp(payoff_matrix))
[v d]=eig(payoff_matrix')
[v d]=eig(payoff_matrix)
[v d]=jordan(payoff_matrix)
1/1.75
1/1.75/4
s=1/1.75/4
payoff_matrix-eye(4)*s
t=payoff_matrix-eye(4)*s
[v d] = eig(t)
s=1/1.75
t=payoff_matrix-eye(4)*s
[v d] = eig(t)
eigen_vectorP
(1/payoff_matrix)*[1;1;1;1]
(1\payoff_matrix)*[1;1;1;1]
(payoff_matrix^(-1))*[1;1;1;1]
t=(payoff_matrix^(-1))*[1;1;1;1]
t/sum(t)
u=t/sum(t)
u'*payoff_matrix*u
%-- 2021/5/7 19:01 --%
ShujieLandscape1228
[U,S,V] = svd(payoff_matrix )
[V,D,W] = eig(payoff_matrix)
[V,D,W] = eig(payoff_matrix')
%-- 2021/5/7 20:02 --%
ShujieLandscape1228
[V,D,W] = eig(eval(D_V_F))
ShujieLandscape1228
[V,D,W] = eig(eval(D_V_F))
V(:,1)'*W(:,1)
V(:,1).*W(:,1)
V(:,1).*W(:,2)
V(:,1).*W(:,3)
V(:,1).*W(:,4)
[V,D,W] = eig(payoff_matrix)
t=(payoff_matrix^(-1))*[1;1;1;1]
t=(payoff_matrix^(-1))
a=4
t=(payoff_matrix^(-1))
payoff_matrix = [0 0 0 a;
1 0 0 0;
0 1 0 0;
0 0 1 0];
t=(payoff_matrix^(-1))
t=(payoff_matrix^(-1))*[1;1;1;1]
t=eig(payoff_matrix^(-1))
[v d]=eig(payoff_matrix^(-1))
[v d]=eig(payoff_matrix)
clear
ShujieLandscape1228
V=eigen_vector;D=eigen_value;A=eval(D_V_F);
A*V - V*D
A*V
W*A
w*A
A*w
w
w(:,4)
w*V
%-- 2021/5/8 13:34 --%
GenTheoMarkoW20190408
plotK20190408
plot(PS(:,16))
plot(PS(:,32))
plotK20190408
payoff88(Para(Pa,1),Para(Pa,2))
payoff_matrix =  [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0]
Strategy5
payoff_matrix =  [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0]
[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
Oneill_replicator
plotK20190408
%-- 2021/5/8 19:06 --%
uiopen('C:\Users\Think\Desktop\2021课程\硕士\liu\eigen3.fig',1)
uiopen('C:\Users\Think\Desktop\2021课程\硕士\liu\test\exp.fig',1)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\msrep7.csv';           [num,txt,raw] = xlsread(abedfilesname);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
d=num(15:62596,1:16);
d2=num(15:62596,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent.csv
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent.csv';[num,txt,raw] = xlsread(abedfilesname);t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
d=num(15:end,1:16)
d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent2.csv';[num,txt,raw] = xlsread(abedfilesname);t=length(txt(:,1));
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent2';[num,txt,raw] = xlsread(abedfilesname);t=length(txt(:,1));
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent2.csv';[num,txt,raw] = xlsread(abedfilesname);t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
mean(d3(500:end,:))
mean(d3(1500:end,:))
mean(d3(2500:end,:))
mean(d3(2:2500,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent3.csv';[num,txt,raw] = xlsread(abedfilesname);t=length(txt(:,1));
d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
mean(d3(1500:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent4.csv';[num,txt,raw] = xlsread(abedfilesname);t=length(txt(:,1));
d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
mean(d3(1500:end,:))
mean(d3(500:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent5.csv';[num,txt,raw] = xlsread(abedfilesname);t=length(txt(:,1));d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
mean(d3(1500:end,:))
mean(d3(500:end,:))
mean(d3(2500:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent5.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent6.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
mean(d3(200:end,:))
sum(ans)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
mean(d3(200:end,:))
mean(d3(500:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent8.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3)
mean(d3(500:end,:))
mean(d3(200:end,:))
mean(d3(500:end,:))
mean(d3(800:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agent9.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(500:end,:))
mean(d3(800:end,:))
mean(d3(1200:end,:))
mean(d3(1500:end,:))
mean(d3(1900:end,:))
mean(d3(2200:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;sum(d3(1:k,:));end
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;sum(d3(1:k,:))];end
plot(e3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agenta.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(2200:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;sum(d3(1:k,:))];end
plot(e3)
s3=load('C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-05\a_4msrepdone.csv');
d3=s3;
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;sum(d3(1:k,:))];end
plot(e3)
[h p]=ttest(d3(:,1),d3(:,3))
[h p]=ttest(d3(2000:end,1),d3(2000:end,3))
[h p]=ttest(d3(20000:end,1),d3(20000:end,3))
[h p]=ttest(d3(50000:end,1),d3(50000:end,3))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3)
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agentb.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(2200:end,:))
mean(d3(1:end,:))
mean(d3(201:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
plot(e3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4-200agentc.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch1.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
plot(e3(:,1)-e3(:,3))
mean(d3(1:end,:))
mean(d3(201:end,:))
mean(d3(2001:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3)
mean(d3(201:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch2.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch3.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
plot(e3(:,1)-e3(:,3))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
%-- 2021/5/9 15:05 --%
ShujieLandscape1228
uiopen('C:\Users\Think\Desktop\2021课程\硕士\liu\eigen4.fig',1)
run('C:\Users\Think\Desktop\2021课程\硕士\liu\eigen_cy.m')
%-- 2021/5/9 15:17 --%
EigenCycle_Copy
uiopen('C:\Users\Think\Desktop\2021课程\硕士\liu\test\exp.fig',1)
equi
a = log([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10;]);
origin = a./(3.*a+1);
ms_rep = [0.140506,0.182685,0.200329,0.213793,0.22199,0.229062,0.234378,0.243280,0.244692,0.250915,0.278179,0.292941,0.302367,0.30952,0.317782,0.321202,0.324545,0.327989,0.332006; ...
0.175026,0.210753,0.227401,0.236551,0.239636,0.242012,0.241739,0.244918,0.250675,0.251416,0.256285,0.260204,0.265230,0.266439,0.267704,0.270041,0.270804,0.271630,0.274573; ...
0.123803,0.166392,0.191142,0.207515,0.217650,0.226747,0.237176,0.239038,0.244176,0.248952,0.280419,0.292804,0.302553,0.309025,0.314203,0.317619,0.321913,0.322804,0.323487;];
plot(a,origin,'r.',a,ms_rep(1,:),'b.',a,ms_rep(2,:),'m.',a,ms_rep(3,:),'c.');
legend('origin','ρ1','ρ2','ρ3');
a = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10;];
origin = a./(3.*a+1);
ms_rep = [0.140506,0.182685,0.200329,0.213793,0.22199,0.229062,0.234378,0.243280,0.244692,0.250915,0.278179,0.292941,0.302367,0.30952,0.317782,0.321202,0.324545,0.327989,0.332006; ...
0.175026,0.210753,0.227401,0.236551,0.239636,0.242012,0.241739,0.244918,0.250675,0.251416,0.256285,0.260204,0.265230,0.266439,0.267704,0.270041,0.270804,0.271630,0.274573; ...
0.123803,0.166392,0.191142,0.207515,0.217650,0.226747,0.237176,0.239038,0.244176,0.248952,0.280419,0.292804,0.302553,0.309025,0.314203,0.317619,0.321913,0.322804,0.323487;];
plot(log(a),origin,'r.',log(a),ms_rep(1,:),'b.',log(a),ms_rep(2,:),'m.',log(a),ms_rep(3,:),'c.');
legend('origin','ρ1','ρ2','ρ3');
plot(log(a),origin,'r.--',log(a),ms_rep(1,:),'b.',log(a),ms_rep(2,:),'m.',log(a),ms_rep(3,:),'c.');
legend('origin','ρ1','ρ2','ρ3');
plot(log(a),origin,'r.--',log(a),ms_rep(1,:),'b*--',log(a),ms_rep(2,:),'m.',log(a),ms_rep(3,:),'c.');
plot(log(a),origin,'r.--',log(a),ms_rep(1,:),'b*-',log(a),ms_rep(2,:),'m.',log(a),ms_rep(3,:),'c.');
plot(log(a),origin,'r.--',log(a),ms_rep(1,:),'b*-',log(a),ms_rep(2,:),'m.',log(a),ms_rep(3,:),'c.=');
plot(log(a),origin,'r.--',log(a),ms_rep(1,:),'b*-',log(a),ms_rep(2,:),'m.',log(a),ms_rep(3,:),'c.-');
log(0.25)
plot(log(a),origin,'r.--',log(a),ms_rep(1,:),'b*-',log(a),ms_rep(2,:),'m*-',log(a),ms_rep(3,:),'c.-');
log(4)
ShujieLandscape1228
syms x1 x2 x3 x4 a  real
%% 支付矩阵
% payoff_matrix = [0 -1 0 1; 1 0 -1 0; 0 1 0 -1; -1 0 1 0]
payoff_matrix = [0 0 0 a;
1 0 0 0;
0 1 0 0;
0 0 1 0];
%% 空间各个点的支付向量
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4]'
%% 各点的支付均值
mean_U = [x1 x2 x3 x4] * Payoff_vector_field_F
%% 动力学方程
V_F = [x1 x2 x3 x4]'.*(Payoff_vector_field_F - mean_U);
%% 全微分矩阵
D_V_F = [diff(V_F,'x1') diff(V_F,'x2') diff(V_F,'x3') diff(V_F,'x4')]
%% 求0点
% S=solve(V_F);
%
x1=a/(3*a + 1);x2=a/(3*a + 1);;x3=a/(3*a + 1);1/4;x4=1/(3*a + 1);
Ne = [x1 x2 x3 x4];
[eigen_vector eigen_value w] = eig(eval(D_V_F))
EigenCycle_Copy
m=eigen_vector(1,1); n=eigen_vector(2,1); eigen_cycle = pi*sin(atan2(imag(m),real(m))-atan2(imag(n),real(n)))*norm(m)*norm(n);
eigen_cycle
clear
syms x1 x2 x3 x4 a  real
%% 支付矩阵
% payoff_matrix = [0 -1 0 1; 1 0 -1 0; 0 1 0 -1; -1 0 1 0]
payoff_matrix = [0 0 0 a;
1 0 0 0;
0 1 0 0;
0 0 1 0];
%% 空间各个点的支付向量
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4]'
%% 各点的支付均值
mean_U = [x1 x2 x3 x4] * Payoff_vector_field_F
%% 动力学方程
V_F = [x1 x2 x3 x4]'.*(Payoff_vector_field_F - mean_U);
%% 全微分矩阵
D_V_F = [diff(V_F,'x1') diff(V_F,'x2') diff(V_F,'x3') diff(V_F,'x4')]
%% 求0点
% S=solve(V_F);
%
x1=a/(3*a + 1);x2=a/(3*a + 1);;x3=a/(3*a + 1);1/4;x4=1/(3*a + 1);
Ne = [x1 x2 x3 x4];
[eigen_vector eigen_value w] = eig(eval(D_V_F))
x2
m=eigen_vector(1,1); n=eigen_vector(2,1); eigen_cycle = pi*sin(atan2(imag(m),real(m))-atan2(imag(n),real(n)))*norm(m)*norm(n);
eigen_cycle
ShujieLandscape1228
m=eigen_vector(1,3); n=eigen_vector(2,3); eigen_cycle = pi*sin(atan2(imag(m),real(m))-atan2(imag(n),real(n)))*norm(m)*norm(n);
eigen_cycle
eigen_cycle = []; for J=1:3; for K=J+1:4; m=eigen_vector(1,3); n=eigen_vector(2,3); tmp = pi*sin(atan2(imag(m),real(m))-atan2(imag(n),real(n)))*norm(m)*norm(n); eigen_cycle = [ eigen_cycle; tmp]; end;end
latex(eigen_cycle)
eigen_cycle = []; for J=1:3; for K=J+1:4; m=eigen_vector(J,3); n=eigen_vector(K,3); tmp = pi*sin(atan2(imag(m),real(m))-atan2(imag(n),real(n)))*norm(m)*norm(n); eigen_cycle = [ eigen_cycle; tmp]; end;end
latex(eigen_cycle)
t2=simplify(eigen_cycle)
latex(t2)
abs(eigen_cycle)
abs(eigen_vector)
t3=abs(eigen_vector(:,3))
a=0.1:0.1:10; plot(a,t3)
a=0.1:0.1:10; plot(a,eval(t3))
a=0.1:0.1:10; plot(a,eval(t3(:,1)))
a=0.1:0.1:10; plot(a,eval(t3(1)))
a=0.1:0.1:10; plot(a,eval(t3))
a=0.1:0.1:10; plot(a,eval(t3(:,4)))
a=0.1:0.1:10; plot(a,eval(t3(1:4)))
a=0.1:0.1:10; plot(a,eval(t3(1))); plot(a,eval(t3(2))); plot(a,eval(t3(3))); plot(a,eval(t3(4)));hold on
a=0.1:0.1:10; plot(a,eval(t3(1)));hold on; plot(a,eval(t3(2)));hold on; plot(a,eval(t3(3)));hold on; plot(a,eval(t3(4)));hold on
b=eval(t3)
t3
b=eval(t3')
a=0.1:0.1:10; plot(a,eval(t3(1)));hold on; plot(a,eval(t3(2)));hold on; plot(a,eval(t3(3)));hold on; plot(a,eval(t3(4)));hold on
t3
s3=sum(abs(t3))
t4=t3/s3
t3=t4
a=0.1:0.1:10; plot(a,eval(t3(1)));hold on; plot(a,eval(t3(2)));hold on; plot(a,eval(t3(3)));hold on; plot(a,eval(t3(4)));hold on
legend('1','2','3','4');
grid on
dot(eigen_vector(:,1),eigen_vector(:,2),)
dot(eigen_vector(:,1),eigen_vector(:,2))
dot(eigen_vector(:,1),eigen_vector(:,3))
dot(eigen_vector(:,1),eigen_vector(:,4))
eigen_vector
dot(eigen_vector(:,2),eigen_vector(:,4))
dot(eigen_vector(:,2),eigen_vector(:,3))
dot(a*eigen_vector(:,1)+eigen_vector(:,2),eigen_vector(:,3))
dot(eigen_vector(:,1).*a+eigen_vector(:,2),eigen_vector(:,3))
dot(eigen_vector(:,1)*a + eigen_vector(:,2),eigen_vector(:,3))
eigen_vector(:,1)*a + eigen_vector(:,2)
eigen_vector(:,1)*a
syms a
dot(eigen_vector(:,1).*a+eigen_vector(:,2),eigen_vector(:,3))
simplify(ans)
eigen_vector(:,1)
eigen_vector(:,1).*a
eigen_vector(:,1).*a+eigen_vector(:,2)
eigen_vector(:,3)
simplify(dot(eigen_vector(:,1).*a+eigen_vector(:,2),eigen_vector(:,3)+[ 0 0 0 (0)]))
simplify(dot(eigen_vector(:,1).*a+eigen_vector(:,2),eigen_vector(:,3)+[ 0 0 0 (0)]'))
simplify(dot(eigen_vector(:,1).*a+eigen_vector(:,2),eigen_vector(:,3)+[ 0 0 0 (-1)]'))
simplify(dot(eigen_vector(:,1).*a+eigen_vector(:,2),eigen_vector(:,3)+[ 0 0 0 (-1+a)]'))
eigen_vector(:,3)+[ 0 0 0 (-1+a)]'
abs(eigen_vector(:,3)+[ 0 0 0 (-1+a)]' )
a=1;eval(abs(eigen_vector(:,3)+[ 0 0 0 (-1+a)]' ))
a=4;eval(abs(eigen_vector(:,3)+[ 0 0 0 (-1+a)]' ))
a=4;eval(abs(eigen_vector(:,3)+[ 0 (-1+a)/a 0 0]' ))
a=4;eval(abs(eigen_vector(:,3)+[ 0 0 (-1+a)/a 0]' ))
a=4;eval(abs(eigen_vector(:,4)+[ 0 0 (-1+a)/a 0]' ))
figure
a=0.1:0.1:10; plot(a,eval(t3(1)));hold on; plot(a,eval(t3(2)));hold on; plot(a,eval(t3(3)));hold on; plot(a,eval(t3(4)));hold on
a=0.1:0.1:10; plot(a,eval(t3(1)),'*');hold on; plot(a,eval(t3(2)));hold on; plot(a,eval(t3(3)));hold on; plot(a,eval(t3(4)));hold on
a=0.1:0.1:10; plot(a,eval(t3(1)),'*');hold on; plot(a,eval(t3(2)),'d');hold on; plot(a,eval(t3(3)),'s');hold on; plot(a,eval(t3(4)));hold on
a=0.1:0.2:5; plot(a,eval(t3(1)),'*');hold on; plot(a,eval(t3(2)),'d');hold on; plot(a,eval(t3(3)),'s');hold on; plot(a,eval(t3(4)));hold on
a=0.1:0.2:5; plot(a,eval(t3(1)),'*');hold on; plot(a,eval(t3(2)),'d');hold on; plot(a,eval(t3(3)),'s');hold on; plot(a,eval(t3(4)),'S');hold on
a=0.1:0.2:5; plot(log(a),eval(t3(1)),'*');hold on; plot(log(a),eval(t3(2)),'d');hold on; plot(log(a),eval(t3(3)),'s');hold on; plot(log(a),eval(t3(4)),'S');hold on
grid on
a=0.1:0.05:5; plot(log(a),eval(t3(1)),'*');hold on; plot(log(a),eval(t3(2)),'d');hold on; plot(log(a),eval(t3(3)),'s');hold on; plot(log(a),eval(t3(4)),'S');hold on
%-- 2021/5/10 22:51 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch4.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
plot(e3)
plot(e3(:,1)-e3(:,3))
plot(e3)
plot(e3(:,1)-e3(:,3))
plot(e3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch4.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch5.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
mean(d3(2000:end,:))
mean(d3(200:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch6.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
mean(d3(1000:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
mean(d3(1000:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
mean(d3(1000:end,:))
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch8.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch8.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch8.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatch9.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
plot(e3(:,1)-e3(:,3))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatcha.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
plot(e3(:,1)-e3(:,3))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatchb.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4agentRanMatchc.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
%-- 2021/5/11 0:00 --%
getEigen
res=getEigen(2,1)
getEigen
getEigen(2,1)
D_Eq_0
%-- 2021/5/11 9:57 --%
getEigen(2,1)
diag(ans(:,17:end))
ev = diag(ans(:,17:end))
getEigen(2,1)
cleaar
getEigen(2,1)
save(strc('ZMX',num2str([m n]),'mat'))
save(strcat('ZMX',num2str([m n]),'mat'))
clear
getEigen(2,1)
clear
load('ZMX1  2mat.mat')
clear
getEigen(2,1)
getEigen(3,2)
getEigen(2,1)
getEigen(4, 2)
getEigen(3, 2)
load('ZMX1  2mat.mat')
load('ZMX2  3mat.mat')
getEigen(2,1)
load('ZMX1  2mat.mat')
abs(eigen_vector)
getEigen(2,1)
plot(A)
bar(A)
histogram(A)
bar(A)
barh(A)
getEigen(2,1)
diff(V_Eq_0,'x1')
sum(NEat)
sum(NEat(1:8))
sum(NEat(9:16))
getEigen(2,1)
t1=image(eigen_vector)
t1=imag(eigen_vector)
sum(t2)
sum(t1)
t0=real(eigen_vector)
getEigen(2,1)
[eigen_vector,eigen_value] = eig(D_Eq_at_NE,'nobalance');
v4=round(eigen_vector,4)
expdata20
clear
load('PE26.mat')
a = sortrows(PS,[1 2 3])
a = sortrows(PS,[1 21 3]);
expdata20
a = sortrows(PS,[1 21 2 3]);
e1 = a(1:12000,[4:11 13 20]);
e1 = a(1:12000,[4:11 13:20]);
[ Yret3 mn] = from_N_colExp_out_am(e1,mean(e1))
e1 = a(12001:24000,[4:11 13:20]);
e1 = a(1:12000,[4:11 13 20]);
e2 = a(12001:24000,[4:11 13:20]);
[ Yret3 mn] = from_N_colExp_out_am(e2,mean(e2))
e3 = a(24001:36000,[4:11 13:20]);
[ Yret3 mn] = from_N_colExp_out_am(e3,mean(e3))
[ Yret2 mn] = from_N_colExp_out_am(e2,mean(e2))
mean(e1)
e1 = a(12001:24000,[4:11 13:20]);
e1 = a(1:12000,[4:11 13:20]);
mean(e1)
[ Yret1 mn] = from_N_colExp_out_am(e1,mean(e1))
getEigen(3, 2)
%-- 2021/5/11 17:06 --%
price20210115
price20210115(1)
load('PE26.mat')
a = sortrows(PS,[1 21 2 3]);
e1 = a(1:12000,[4:11 13:20]);
mean(e1)
[ Yret1 mn] = from_N_colExp_out_am(e1,mean(e1))
dos_x_t=e1;  [L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(dos_x_t)
e2 = a(12001:24000,[4:11 13:20]);
dos_x_t=e2;  [L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(dos_x_t)
e3 = a(24001:36000,[4:11 13:20]);
dos_x_t=e3;  [L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(dos_x_t)
figure
dos_x_t=e1;  [L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(dos_x_t)
figure
dos_x_t=e2;  [L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(dos_x_t)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a1.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a1.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a2.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(s3(:,1)-d3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
[h p]=ttest(d3(:,1)-d3(:,2))
[h p]=ttest(d3(:,3)-d3(:,2))
[h p]=ttest(d3(:,3)-d3(:,4))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a2.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a2.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(s3(:,1)-d3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a2.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a2.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
%-- 2021/5/11 22:12 --%
abedfilesname='C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-05\0511.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
d2=load('C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4_1-1data.csv');
d3=d2(:,2:5);
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
d3=d2(:,[2:5]+5);
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
d3=d2(:,[2:5]+5*2);
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
d2=load('C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4_1-2data.csv');
d3=d2(:,2:5);
d2=load('C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4_1-2data.csv');
d3=d2(:,2:5);
mean(d3(1:end,:))
d2=load('C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4_1-3data.csv');
d3=d2(:,2:5);
mean(d3(1:end,:))
zz=mean(d2)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a2.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a3.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
%-- 2021/5/12 0:02 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a4.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
(d3(:,1)-d3(:,3))
histogram(ans)
plot(ans)
histogram(d3(:,1))
histogram(d3)
bar(d3,'DisplayName','d3')
barh(d3,'DisplayName','d3')
plotmatrix(d3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a4.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a5.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
%-- 2021/5/12 8:05 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a6.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
plot(e3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a6.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
plot(e3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a7.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d2=num(3000:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a8.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a9.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(15:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a4a9.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(215:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a25a.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(215:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
%-- 2021/5/12 12:05 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a25b.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(215:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
plot(e3)
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a25b.csv';[num,txt,raw] = xlsread(abedfilesname); d2=num(215:end,2:4:15);
d3=[d2(:,1)-d2(:,2) d2(:,2)-d2(:,3) d2(:,3)-d2(:,4) d2(:,4)];
mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021a.csv';[num,txt,raw] = xlsread(abedfilesname);
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\a25b.csv'; d3=num(35:end,2:5);
mean(d3(1:end,:))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021a.csv';[num,txt,raw] = xlsread(abedfilesname);
d3=num(35:end,3:6);
mean(d3(1:end,:))
d3=num(35:end,[3:6]+5); mean(d3(1:end,:))
e3=d3(1,:); for k=2:length(d3(:,1)); e3=[e3;(d3(k,:)+e3(k-1,:))];end
plot(e3(:,1)-e3(:,3))
[h p]=ttest(d3(:,1)-d3(:,3))
d3=num(35:end,[3:6]+5*2); mean(d3(1:end,:))
d3=num(35:end,[3:6]+5*3); mean(d3(1:end,:))
%-- 2021/5/12 13:32 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021b.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'DisplayName','rr')
rr=[]; for k=0:19; d3=num(25000:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'DisplayName','rr')
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'O')
plot(rr,'O-')
%-- 2021/5/12 14:21 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021b2.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'DisplayName','rr')
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'O-')
r=rr;
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021b.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'DisplayName','rr')
r=[r;rr];
plot(r,'O-')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021c100.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'DisplayName','rr')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021c100m001.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021cR001.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021cR001_1.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021a.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021b.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:19; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021c.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:23; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
legend('x1','x2','x3','x4',)
legend('x1','x2','x3','x4')
legend('x1','x2','x3','x4','fontsize',20)
legend('x1','x2','x3','x4','font',20)
title('ShujieTwist2021c')
%-- 2021/5/12 22:23 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021d2000.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:23; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
plot(rr,'DisplayName','rr')
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021c.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:23; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
plotmatrix(d3)
plotmatrix(d3)
plotmatrix(d3(1:20,:))
plotmatrix(d3(1:40,:))
plotmatrix(d3(1:60,:))
plotmatrix(d3(61:80,:),'ro')
plotmatrix(d3(1:60,:),'b.')
hold on
plotmatrix(d3(61:80,:),'ro')
hold on
plotmatrix(d3(1:60,:),'b.')
plotmatrix(d3(1:60,:),'b.'d3(61:80,:),'ro')
plotmatrix(d3(1:60,:),'b.',d3(61:80,:),'ro')
clf; plotmatrix(d3(1:60,:),'b.'); hold on ; plotmatrix(d3(61:80,:),'ro')
clf; plotmatrix(d3(1:60,:),'b.'); hold on ; plotmatrix(d3(61:80,:),'ro'); hold on ;
clf; plotmatrix(d3(1:60,:),'b.'); hold on ; figure;plotmatrix(d3(61:80,:),'ro'); hold on ;
clf; plotmatrix(d3(1:60,:),'b.'); hold on ; figure;plotmatrix(d3(1:80,:),'ro'); hold on ;
figure;plotmatrix(d3(1:800,:),'.'); hold on ;
figure;plotmatrix(d3(1:400,:),'r.'); hold on ;
figure;plotmatrix(d3(401:1200,:),'r.'); hold on ;
%-- 2021/5/12 23:53 --%
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021e.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:23; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
grid on
%-- 2021/5/13 0:40 --%
load('C:\Users\Think\Desktop\2021毕设\周名铉\zmx20210513\2_1(6)\exp\exp_6.mat')
eval(x10)
eval([x1 x2 x3 x4 x5 x6 x7 x8])
eval(data)
w1=eval(data);
%-- 2021/5/13 4:19 --%
timesies2figure8
dos_x_t=load('C:\Users\Think\Desktop\12.csv');
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=timesies2figure8(dos_x_t)
timesies2figure8
timesies2figure8(dos_x_t)
dos_x_t=load('C:\Users\Think\Desktop\12.csv');
timesies2figure8(dos_x_t)
dos_x_t=load('C:\Users\Think\Desktop\12.csv');
timesies2figure8(dos_x_t)
tp1 = abs(eigenvector_set(:,2))
plot(tp1)
dos_x_t=load('C:\Users\Think\Desktop\13.csv');
timesies2figure8(dos_x_t)
dos_x_t=load('C:\Users\Think\Desktop\13.csv');
timesies2figure8(dos_x_t)
dos_x_t=load('C:\Users\Think\Desktop\12.csv');
timesies2figure8(dos_x_t)
%-- 2021/5/13 5:03 --%
dos_x_t=load('C:\Users\Think\Desktop\6.csv');
timesies2figure8(dos_x_t)
mean(dos_x_t)
tp1 =  (eigenvector_set(:,1))-mean(dos_x_t)'
sum(mean(dos_x_t)')
sum((eigenvector_set(:,1))')
sum((eigenvector_set(:,1)).^2)
((eigenvector_set(:,1)).^2)
((eigenvector_set(:,1)).^(1/2))
((eigenvector_set(:,1)))
readc
%-- 2021/5/13 10:23 --%
load('F:\cDownload\PokerAnan20200716\8x8abedata\PS26.mat')
uiopen('C:\Users\Think\AppData\Local\Temp\Rar$DRa0.936\8.fig',1)
load('F:\8x8code\8x8abedata\PS25.mat')
load('F:\8x8code\8x8abedata\PS24.mat')
load('F:\8x8code\8x8abedata\eigenVV8x8.mat')
readc
Ctime
load('PE26.mat');
dos_x_t=PS(:,[4:11 13:20]);
timesies2figure8(dos_x_t)
load('PE26.mat')
a = sortrows(PS,[1 21 2 3]);
e1 = a(1:12000,[4:11 13:20]);
timesies2figure8(e1)
a400=a(find(a(:,3)>400),:)
e1 = a400(1:7200,[4:11 13:20]);
timesies2figure8(e1)
e1 = a400(1+7200:7200+7200,[4:11 13:20]);
figure; timesies2figure8(e1)
e1 = a400(1+7200+7200:7200+7200+7200,[4:11 13:20]);
figure; timesies2figure8(e1)
figure; timesies2figure8(e1)
eigenvector_set(:,maxImag_eigenvalueId)
abs(eigenvector_set(:,maxImag_eigenvalueId))
abs(eigenvector_set(:,13))
svd(dos_x_t)
[U,S,V] = svd(dos_x_t)
[U,S,V] = svd(dos_x_t);
e1 = a400(1:7200,[4:11 13:20]);
[U,S,V] = svd(e1);
plot(U(1,:))
plot(V(:,1))
plot(V(:,2))
plot(V(:,1:2),'DisplayName','V(:,1:2)')
plot(V(:,1:3),'DisplayName','V(:,1:3)')
e1 = a400(1+7200:7200+7200,[4:11 13:20]);
[U,S,V] = svd(e1);
plot(U(1,:))
plot(U(2,:))
plot(V(:,1:2),'DisplayName','V(:,1:2)')
plot(V(:,1))
e1 = a400(1:7200,[4:11 13:20]);
[U,S,V] = svd(e1);
e1 = a400(1+7200+7200:7200+7200+7200,[4:11 13:20]);
[U,S,V] = svd(e1);
plot(V(:,1))
plot(U(:,265))
e1 = a400(1:7200,[4:11 13:20]);
[U,S,V] = svd(e1);
plot(V(:,1))
plot(V(:,2))
%-- 2021/5/13 14:34 --%
load('PS26.mat')
load('F:\8x8code\8x8abedata\PS26.mat')
plot(PS(:,4:11),'DisplayName','PS(:,4:11)')
plot(PS(:,13:end),'DisplayName','PS(:,13:end)')
subplot(1,2,1); plot(PS(:,4:11),'DisplayName','PS(:,4:11)');
subplot(1,2,2); plot(PS(:,13:end),'DisplayName','PS(:,13:end)');
subplot(2,1,1); plot(PS(:,4:11),'DisplayName','PS(:,4:11)');subplot(2,1,2); plot(PS(:,13:end),'DisplayName','PS(:,13:end)');
a = sortrows(PS,[1 21 2 3]);
a = sortrows(PS,[1 2 3]);
b=PS(find(PS(:,1)==6),:);
a = sortrows(b,[1 2 3]);
timesies2figure8(a)
e1 = a(:,[4:11 13:20]);
timesies2figure8(e1)
a = sortrows(b,[1 2 3]);
e1 = a(:,[4:11 13:20]);
timesies2figure8(e1)
load('PS26.mat')
clear
load('PS26.mat')
load('F:\8x8code\8x8abedata\PS26.mat')
b=PS(find(PS(:,1)==6),:);
a = sortrows(b,[1 2 3]);
e1 = a(:,[4:11 13:20]);
timesies2figure8(e1)
b=PS(find(PS(:,1)==8),:);
a = sortrows(b,[1 2 3]);
e1 = a(:,[4:11 13:20]);
timesies2figure8(e1)
b=PS(find(PS(:,1)==10),:);
a = sortrows(b,[1 2 3]);
e1 = a(:,[4:11 13:20]);
timesies2figure8(e1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\scatter6.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\scatter8.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\scatter6_abed.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\scatter10_abed.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\scatter8_abed.fig',1)
load('C:\Users\Think\Desktop\2021毕设\周名铉\EigenLogic50.mat')
clear
load('C:\Users\Think\Desktop\2021毕设\周名铉\EigenLogic50.mat')
eval(EigenVector_10)
EigenVector_10d = eval(EigenVector_10)
eval(EigenValue_10)
EigenValue_10d = eval(EigenValue_10)
diag(EigenValue_10d)
EigenValue_6d = eval(EigenValue_6)
EigenVector_6d = eval(EigenVector_6)
diag(EigenValue_6d)
absvec6=abs(EigenVector_6)
absvec6=abs(EigenVector_6d)
plot(absvec6,'DisplayName','absvec6')
clear
[a b c]=xlsread('C:\Users\Think\Desktop\2021毕设\周名铉\allData.xlsx')
plotmatrix(a(:,1:9))
plotmatrix(a(:,1:3:9))
figure(plotmatrix(a(:,1:3:9)))
figure;(plotmatrix(a(:,1:3:9)))
figure;plotmatrix(a(:,2:3:9),'r')
figure;plotmatrix(a(:,2:3:9),'r.')
figure;plotmatrix(a(:,3:3:9),'k.')
%-- 2021/5/13 23:50 --%
abedfilesname='C:\Users\Think\Desktop\abed-Pocker-Experiment40.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:35; d3=num(38:end,[3:18]+17*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
figure; rr=[]; for k=0:35; d3=num(38:end,[3:18]+17*k); rr=[rr; mean(d3(501:end,:))];end
plot(rr,'S-')
abedfilesname='F:\cDownload\phase_exp Strategy distributions_0.1.csv';[num,txt,raw] = xlsread(abedfilesname);
d3=num(16:end,[2:4:14]);
d4 = [d3(:,1)-d3(:,2) d3(:,2)-d3(:,3)  d3(:,3)-d3(:,4)  d3(:,4)];
plot(d4,'DisplayName','d4')
mean(d4)
mean(d4(200:end,:))
mean(d4(600:end,:))
mean(d4(400:end,:))
ShujieLandscape1228
plot(t2(:,1:4),'DisplayName','t2(:,1:4)')
mean(t2(:,1:4))
plot(ans(1,1:3))
ShujieLandscape1228
mean(t2(:,1:4))
plot(ans(1,1:3))
mean(t2(:,1:4))
plot(ans(1,1:3))
plot(t2(:,end))
plot(t2(:,2:3),'DisplayName','t2(:,2:3)')
plot(t2(:,1:2),'DisplayName','t2(:,1:2)')
plot(t2(:,[1,3]),'DisplayName','t2(:,[1,3])')
scatter(getcolumn(t2(:,[1,3]),1),getcolumn(t2(:,[1,3]),2))
scatter(getcolumn(t2(1:100,[1,3]),1),getcolumn(t2(:,[1,3]),2))
scatter(getcolumn(t2(1:100,[1,3]),1),getcolumn(t2(1:100,[1,3]),2))
scatter(getcolumn(t2(1:30,[1,3]),1),getcolumn(t2(1:30,[1,3]),2))
scatter(getcolumn(t2(1:31,[1,3]),1),getcolumn(t2(1:31,[1,3]),2))
mean(t2(2000:endm1:4))
mean(t2(2000:end,1:4))
plot(ans(1,1:3))
ans(1)-ans(3)
ans(2)-ans(3)
mean(t2(2000:end,1:4))
ans(2)-ans(3)
figure; scatter(t2(2001:end,2), t2(2001:end,4))
figure; scatter(t2(2001:end,2), t2(2001:end,6))
figure; scatter(t2(2001:end,6), t2(2001:end,4))
figure;plot(t2(5000:end,1:3))
abs(c2)
abs(c2).*abs(c2)
(c2)*(c2')
(c2).*(c2)
(conj(c2)).*(c2)
figure; scatter(t2(10:3000,3), t2(10:3000,6))
figure; scatter(t2(10:3000,1), t2(10:3000,6))
figure; scatter(t2(10:3000,2), t2(10:3000,6))
ShujieLandscape1228
plot(t2(:,3:4),'DisplayName','t2(:,3:4)')
scatter(getcolumn(t2(:,2:3),1),getcolumn(t2(:,2:3),2))
scatter(getcolumn(t2(:,1:2),1),getcolumn(t2(:,1:2),2))
scatter(getcolumn(t2(:,[1,3]),1),getcolumn(t2(:,[1,3]),2))
scatter(getcolumn(t2(:,[2,4]),1),getcolumn(t2(:,[2,4]),2))
scatter(getcolumn(t2(:,[4,end]),1),getcolumn(t2(:,[4,end]),2))
scatter(getcolumn(t2(10:end,[4,end]),1),getcolumn(t2(10:end,[4,end]),2))
%-- 2021/5/14 10:54 --%
abedfilesname='F:\cDownload\phase_exp Strategy distributions_0.1.csv';[num,txt,raw] = xlsread(abedfilesname);
d3=num(16:end,[2:4:14]);
[ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
[ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021d2000.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=23:23; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
rr=[]; for k=3:3; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
[ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
plot(am_eigencycleDim_t(:,3))
rr=[]; for k=3:3; d3=num(35000:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
[ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
sum(am_eigencycleDim_t)
clear
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021d2000.csv';[num,txt,raw] = xlsread(abedfilesname);
[ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
rr=[]; for k=3:3; d3=num(35000:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
[ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
sum(am_eigencycleDim_t)
%-- 2021/5/14 13:25 --%
load('F:\8x8code\8x8abedata\PS26.mat')
load('PE26.mat')
a = sortrows(PS,[1 21 2 3]);
e1 = a(1:12000,[4:11 13:20]);
d3=e1; [ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
plot(am_eigencycleDim_t,'DisplayName','am_eigencycleDim_t')
histogram(am_eigencycleDim_t(:,1:2))
plot(am_eigencycleDim_t(:,23))
histogram(am_eigencycleDim_t(:,23))
plot(Yret3)
plot(am_eigencycleDim_t(:,23))
bar(am_eigencycleDim_t(:,23))
histogram(am_eigencycleDim_t(:,23))
histogram(am_eigencycleDim_t(:,25))
d3=e1; [ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3)*0)
histogram(am_eigencycleDim_t(:,25))
histogram(am_eigencycleDim_t(:,25),5)
histogram(Yret3)
plot(Yret3)
e1 = a(1+12000:12000+12000,[4:11 13:20]);
d3=e1; [ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3))
d3=e1; [ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(d3,mean(d3)*0);
plot(Yret3)
semilogy(Yret3)
semilogx(Yret3)
loglog(Yret3)
comet(Yret3)
stem(Yret3)
stairs(Yret3)
barh(Yret3)
bar(Yret3)
tmp=[0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
-0.1964	0.1964	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0.0654	-0.0654	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
0.0654	-0.0654	0	0	0	0	0	0
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	-0.0436	0.0436	-0.0436	0.0436
-0.0218	0.0218	0	0	0.0873	-0.0873	0.0873	-0.0873
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
0	0	0	0	0.0756	-0.0756	-0.0756	0.0756
0	0	0	0	-0.0756	0.0756	0.0756	-0.0756
]
tmp(:,1)'*tmp(:,5)
tmp(:,1)'*tmp(:,7)
tmp(:,5)'*tmp(:,7)
Strategy5
clear
payoff_matrix = -0.001 * [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0];
[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
tmp= eigen_vector; tmp(:,5)'*tmp(:,3)
tmp= eigen_vector; tmp(:,5)'*tmp(:,1)
tmp= eigen_vector; tmp(:,3)'*tmp(:,4)
tmp= eigen_vector; tmp(:,3).*tmp(:,4)
tmp= eigen_vector; tmp(:,3).*tmp(:,3)
tmp= eigen_vector; tmp(:,1).*tmp(:,2)
tmp= eigen_vector; sum(tmp(:,1).*tmp(:,2))
tmp= eigen_vector; sum(tmp(:,4).*tmp(:,2))
[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
[eigen_vector eigen_value] = jordan(eval(D_V_F))
tmp= eigen_vector; sum(tmp(:,1).*tmp(:,2))
tmp= eigen_vector; tmp(:,5)'*tmp(:,1)
Lmn=[]; for k=1:size(eigen_vector,2)
Lmn = [Lmn round(from_eigenvector_out_am(eigen_vector(:,k)),4)];
end;
Em =diag(eigen_value)
sum(Lmn)
7365/4603
eigen_vector(:,1)
sum(real(eigen_vector(:,1)))
maxVnash = real(eigen_vector(:,3))'/sum(real(eigen_vector(:,3))) - Nash;
clear
ShujieLandscape1228
p=[1 1 1 4]/7; v=eigen_vector(:,2)
p
p'
p*v
del=0.02; p(1)=p(1)+delta; p(3)=p(3)-delta
del=0.02; p(1)=p(1)+del; p(3)=p(3)-del
del=0.02; p(2)=p(2)+2*del;
p*v
p(4)=p(4) - 3*del;
p*v
ShujieLandscape1228
eigen_vector(:,3)
syms x1 x2 x3 x4 a real
r= (a+x1)*(a*(1/4 + 3i/4) - 1/4 + 1i/4) ...
+ (a+x2)*(- a/2 - 1/2) ...
+ (a+x3)*(a*(1/4 - 3i/4) - 1/4 - 1i/4) ...
+ (a+x4)*(1)
s=subs(r,x4,-x1-x2-x3)
real(s)
image(s)
imag(s)
clear
syms x1 x2 x3 x4 a real
r= (a+x1)*(a*(1/4 + 3i/4) - 1/4 + 1i/4) ...
+ (a+x2)*(- a/2 - 1/2) ...
+ (a+x3)*(a*(1/4 - 3i/4) - 1/4 - 1i/4) ...
+ (1+x4)*(1)
s=subs(r,x4,-x1-x2-x3)
imag(s)
real(s)
expand(imag(s))
factor(imag(s))
reals=real(s)
subs(reals,x1,x3)
realt = subs(reals,x1,x3)
factor(realt)
S=solve('2*a + 3*x2 + 5*x3 + a*x2 - a*x3 - 2',x2)
syms x1 x2 x3 x4 a real
r = (a+x1)*(a*(1/4 + 3i/4) - 1/4 + 1i/4) ...
+ (a+x2)*(- a/2 - 1/2) ...
+ (a+x3)*(a*(1/4 - 3i/4) - 1/4 - 1i/4) ...
+ (1+x4)*(1)
solve(real(r),imag(r),'x1+x2+x3+x4','x1')
real(r)
imag(r)
ShujieLandscape1228
v=eigen_vector;r= (v(:,1)*x1+v(:,2)*x2)'*v(:,3)
v=eigen_vector;syms x1 x3 real; r= (v(:,1)*x1+v(:,2)*x2)'*v(:,3)
real(r)
imag(r)
x2
v=eigen_vector;syms x1 x2 real; r= (v(:,1)*x1+v(:,2)*x2)'*v(:,3)
real(r)
imag(r)
simplify(r)
v=eigen_vector;syms x1 x2 real; r= (v(:,1)*x1+v(:,2)*x2).*v(:,3)
abs(r)
v=eigen_vector;syms x1 x2 a real; r= (v(:,1)*x1+v(:,2)*x2).*v(:,3)
abs(r)
sum(r)
%-- 2021/5/14 22:52 --%
ShujieLandscape1228
v=eigen_vector;syms x1 x2 a real; r= (v(:,1)*x1+v(:,2)*x2).*v(:,3)
sum(r)
r2=sum(r)
r3=real(r2)
r4=imag(r2)
v
a=4; w=eval(v)
w3 = abs(w).*cos(angle(w))
w3 = abs(w(:,3)).*cos(angle(w(:,3)))
syms x1 x2 a real; r= (v(:,1)*x1+v(:,2)*x2).*w3
clear;abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021d2000.csv';[num,txt,raw] = xlsread(abedfilesname);
clear;abedfilesname='C:\Users\Think\Desktop\a4.csv';[num,txt,raw] = xlsread(abedfilesname);
d3=num(16:end,[2:4:14]);
d4 = [d3(:,1)-d3(:,2) d3(:,2)-d3(:,3)  d3(:,3)-d3(:,4)  d3(:,4)];
plot(d4,'DisplayName','d4')
plot(d4(:,[2 4]))
grid on
mean(d4)
clear;abedfilesname='C:\Users\Think\Desktop\a4.csv';[num,txt,raw] = xlsread(abedfilesname);
d3=num(16:end,[2:4:14]);
d4 = [d3(:,1)-d3(:,2) d3(:,2)-d3(:,3)  d3(:,3)-d3(:,4)  d3(:,4)];
mean(d4)
plot(d4(:,[2 4]))
clear
abedfilesname='C:\Users\Think\Desktop\2021课程\硕士\liu\wzj0508\ShujieTwist2021c.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:23; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
mean(d3)
plot(d3(:,[2,end]),'DisplayName','d3(:,[2,end])')
plot(d3(4:109,:),'DisplayName','d3(4:109,:)')
plot(d3(21590:22044,:),'DisplayName','d3(21590:22044,:)')
plot(d3(21090:22044,:),'DisplayName','d3(21590:22044,:)')
plot(d3(21090:26044,:),'DisplayName','d3(21590:22044,:)')
plot(d3(21090:26044,:),'DisplayName','d3(21590:22044,:)');grid on
clear
abedfilesname='C:\Users\Think\Desktop\a4.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:3; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
plot(d3,'DisplayName','d3')
mean(d3(1:end,:))
plot(d4(:,[2 4]))
plot(d3(:,[2 4]))
grid on
plot(d3(:,[1 3 4]))
plot(d3(:,[4]))
plot(d3(:,[2]))
plot(d3(:,[3]))
plot(d3(:,[1]))
sactter(d3(:,[2 4]))
scatter(d3(:,[2 4]))
scatter(d3(:,[2 ]), d3(:,[4]))
comet(d3(:,2))
comet(d3(:,[2 ]), d3(:,[4]))
animatedline(d3(:,[2 ]), d3(:,[4]))
comet(d3(:,[2 ]), d3(:,[4]))
for k=1:10; comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[4]));end
for k=1:50; comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[4]));end
for k=1:50; comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[3]));end
xlim([0,1])
for k=1:50; comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[3]));end
for k=1:50; xlim([0,1]);hold on;comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[3]));end
for k=1:50; xlim([0,1]);ylim([0,1]);hold on;comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[3]));end
for k=1:50; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);hold on;comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[3]));end;end
for k=1:50; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square; hold on;comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[3]));end;end
for k=1:50; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:100]+100*k,[2 ]), d3([1:100]+100*k,[3]));end;end
o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:50; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:100]+100*k,[o(fid,1)]), d3([1:100]+100*k,[o(fid,2)]));end;end
clf;o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:50; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:100]+100*k,[o(fid,1)],0.3), d3([1:100]+100*k,[o(fid,2)]),0.3);end;end
clf;o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:50; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:100]+100*k,[o(fid,1)]), d3([1:100]+100*k,[o(fid,2)]),0.3);end;end
clf;o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:5; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
uiopen('C:\Users\Think\Desktop\a4.fig',1)
clf;o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:25:1000]+1000*k,[o(fid,1)]), d3([1:25:1000]+1000*k,[o(fid,2)]),0.3);end;end
clf;o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:15:1000]+1000*k,[o(fid,1)]), d3([1:15:1000]+1000*k,[o(fid,2)]),0.3);end;end
clf;o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
%-- 2021/5/15 5:33 --%
abedfilesname='C:\Users\Think\Desktop\a4.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:3; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
clf;o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
mean(d3)
clf;a=4; NE = [a a a 1]/(3*a+1); o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;
plot(NE(o(fid,1)),NE(o(fid,1)),'X');
comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
clf;a=4; NE = [a a a 1]/(3*a+1); o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'X');
comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
clf;a=4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'X');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
abedfilesname='C:\Users\Think\Desktop\a4.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:3; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];D(:,:,k+1)=d3;end
d3=D(:,:,1);
clf;a=1/4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'X');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
d3=D(:,:,4);
figure;a=4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'X');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
d3=D(:,:,1);
clf;a=1/4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'X');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
clf;a=1/4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'rX');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
%-- 2021/5/15 10:56 --%
abedfilesname='C:\Users\Think\Desktop\a4-rmf.csv';[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:3; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];end
plot(rr,'S-')
figure;a=4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'X');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
rr=[]; for k=0:3; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];D(:,:,k+1)=d3;end; d3=D(:,:,1);
clf;a=1/4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'rX');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
rr=[]; for k=0:3; d3=num(35:end,[3:6]+5*k); rr=[rr; mean(d3(1:end,:))];D(:,:,k+1)=d3;end; d3=D(:,:,2);
clf;a=1/4; NE = [a a a 1]/(3*a+1); ME=mean(d3);o=[1 2;1 3;1 4; 2 3; 2 4; 3 4]; for k=1:2; for fid=1:6; subplot(2,3,fid); xlim([0,1]);ylim([0,1]);axis square;box on; hold on;plot(NE(o(fid,1)),NE(o(fid,2)),'rX');plot(ME(o(fid,1)),ME(o(fid,2)),'s');comet(d3([1:5:1000]+1000*k,[o(fid,1)]), d3([1:5:1000]+1000*k,[o(fid,2)]),0.3);end;end
clear
ShujieMovie
plot(NE(o(fid,1)),NE(o(fid,2)),'rX','markersize',20);
plot(NE(o(fid,1)),NE(o(fid,2)),'rX','markersize',20,'linewidth',12);
plot(NE(o(fid,1)),NE(o(fid,2)),'rX','markersize',20,'linewidth',4);
ShujieMovie
xlabel(num2str(o(fid,1)),'fontsize',12);
ShujieMovie
%-- 2021/5/15 13:07 --%
ShujieMovie
uiopen('C:\Users\Think\Desktop\新建文件夹 (2)\ShujieTwist2021e.fig',1)
A=load('C:\Users\Think\Desktop\A25sj.csv');
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),A(:,2))*10+max(A(:,1),A(:,2));
s16 = []; for I=1:4; for J=1:4; s16 = [s16; I J I*10+J min(I,J)*10+max(I,J)];end;end;S10=unique(s16(:,4));S16=(s16(:,3));
M10=zeros(10); for K=1:length(A(:,1))-1; m=find(S10==A(K,4));n=find(S10==A(K+1,4));M10(m,n)=M10(m,n)+1;end
sum(sum(M10))
Asymm = M10-M10';
A=load('C:\Users\Think\Desktop\书洁\A25sj.csv');
A(:,3)= A(:,1)*10+A(:,2);A(:,4)= min(A(:,1),A(:,2))*10+max(A(:,1),A(:,2));
s16 = []; for I=1:4; for J=1:4; s16 = [s16; I J I*10+J min(I,J)*10+max(I,J)];end;end;S10=unique(s16(:,4));S16=(s16(:,3));
M10=zeros(10); for K=1:length(A(:,1))-1; m=find(S10==A(K,4));n=find(S10==A(K+1,4));M10(m,n)=M10(m,n)+1;end
sum(sum(M10))
Asymm = M10-M10'
clear
ShujieMovie
%-- 2021/5/16 10:29 --%
a=[1 2 3 4]; b=[a; -1*a;4*a;6*a]
[v d]=eig(b)
v(:,1)'*v(:,2)
[U S V]=svd(b)
U(:,1)*V(:,1)'
U(:,1)'*V(:,1)
U(:,1)*V(:,1)'*40.24
U(:,1)*V(:,1)'*S(1,1)
%-- 2021/5/16 13:04 --%
ShujieLandscape1228
eval(D_V_F)
eval(V_F)
simplify(eval(V_F))
latex(eigen_vector)
latex(factor(eigen_vector))
latex(simplify(eigen_vector))
latex(simplify(eigen_value))
latex(simplify(diag(eigen_value))')
latex([simplify(diag(eigen_value))';simplify(eigen_vector)])
angle(eigen_vector(1,3),eigen_vector(2,3))
angle([1 0], [1/2 1/2])
eigen_cy
eigen_cy(eigen_vector(1,3)), eigen_vector(1,4)))
eigen_cy(eigen_vector(1,3)), eigen_vector(1,4))
eigen_cy(eigen_vector(1,3), eigen_vector(1,4))
latex(eigen_cy(eigen_vector(1,3), eigen_vector(1,4)))
latex(simplify(eigen_cy(eigen_vector(1,3), eigen_vector(1,4))))
eigen_vector(1,4)
latex(simplify(eigen_cy(eigen_vector(1,3), eigen_vector(4,3))))
latex(simplify(eigen_cy(eigen_vector(3,3), eigen_vector(4,3))))
latex(simplify(eigen_cy(eigen_vector(2,3), eigen_vector(3,3))))
latex(simplify(eigen_cy(eigen_vector(1,3), eigen_vector(3,3))))
latex(simplify(eigen_cy(eigen_vector(1,3), eigen_vector(3,3))) / (eigen_cy(eigen_vector(1,3), eigen_vector(4,3)))
latex(simplify(eigen_cy(eigen_vector(1,3), eigen_vector(3,3))) / (eigen_cy(eigen_vector(1,3), eigen_vector(4,3))) ))
latex(simplify(eigen_cy(eigen_vector(1,3), eigen_vector(3,3))  / (eigen_cy(eigen_vector(1,3), eigen_vector(4,3))))
latex(simplify(eigen_cy(eigen_vector(1,3), eigen_vector(3,3))  / (eigen_cy(eigen_vector(1,3), eigen_vector(4,3)))))
-0.7406/-0.1746
0.9177/4
%-- 2021/5/16 22:08 --%
load('C:\Users\Think\Desktop\2021毕设\周名铉\eigen.mat')
ShujieLandscape1228
id=[1 2; 1 3; 1 4;2 3; 2 4;3 4]; eygencycleset_rr=[]; for k=1:6; eygencycleset_rr=[eygencycleset_rr; eigen_cy(eigen_vector(id(k,1),3), eigen_vector(id(k,2),3)) ;end
id=[1 2; 1 3; 1 4;2 3; 2 4;3 4]; eygencycleset_rr=[]; for k=1:6; eygencycleset_rr=[eygencycleset_rr; eigen_cy(eigen_vector(id(k,1),3), eigen_vector(id(k,2),3))] ;end
ShujieLandscape1228
id=[1 2; 1 3; 1 4;2 3; 2 4;3 4]; eygencycleset_rr=[]; for k=1:6; eygencycleset_rr=[eygencycleset_rr; eigen_cy(eigen_vector(id(k,1),3), eigen_vector(id(k,2),3))] ;end
a=0.1:0.1:10; plot(eygencycleset_rr)
a=0.1:0.1:10; plot(a,eygencycleset_rr)
a=0.1:0.1:10; plot(a,eval(eygencycleset_rr))
a=0.1:0.1:10; plot(log(a),eval(eygencycleset_rr))
norm([1 2 3])
a=0.1:0.1:10; plot(log(a),eval(eygencycleset_rr)/norm(eval(eygencycleset_rr)))
norm(eval(eygencycleset_rr))
a=0.1;norm(eval(eygencycleset_rr))
a=0.1:0.1:10; plot(log(a),eval(eygencycleset_rr)/1)
for a=0.1:0.1:10;  rr=[tt;eval(eygencycleset_rr)/norm(eval(eygencycleset_rr))];end
tt=[];for a=0.1:0.1:10;  tt=[tt;eval(eygencycleset_rr)/norm(eval(eygencycleset_rr))];end
tt=[];for a=0.1:0.1:10;  tt=[tt;eval(eygencycleset_rr)/norm(eval(eygencycleset_rr))'];end
tt=[];for a=0.1:0.1:10;  tt=[tt;eval(eygencycleset_rr')/norm(eval(eygencycleset_rr))];end
plot(tt,'DisplayName','tt')
tt=[];for a=0.1:0.1:10;  tt=[tt;log(a) eval(eygencycleset_rr')/norm(eval(eygencycleset_rr))];end
plot(tt,'DisplayName','tt')
plot(tt(:,1),tt(:,[2:7]) )
tt=[];for a=0.05:0.1:10;  tt=[tt;log(a) eval(eygencycleset_rr')/norm(eval(eygencycleset_rr))];end
plot(tt(:,1),tt(:,[2:7]),'o-' )
tt=[];for a=0.05:0.1:10;  tt=[tt;log(a) eval(eygencycleset_rr')/norm(eval(eygencycleset_rr)) a];end
plot(tt(:,1),tt(:,[2:7]),'-','linewidth',4);xlim([-2.2 2.2]);ylim([-0.8 0.8]);box on;grid on;
set(gca,'fontsize',15)
line([log(0.25) log(0.25)],[-0.8 0.8],'LineStyle','--','linewidth',2)
line([log(4) log(4)],[-0.8 0.8],'LineStyle','--','linewidth',2)
line([log(1) log(1)],[-0.8 0.8],'LineStyle','--','linewidth',2)
legend('\sigma(1 2)','\sigma(1 3)','\sigma(1 4)','\sigma(2 3)','\sigma(2 4)','\sigma(3 4)','Location','NorthOutside','Orientation','horizontal');legend('boxoff')
xlabel('ln($a$)','interpreter','Latex'),ylabel('eigencycle value $\sigma$', 'interpreter','Latex')
tool_figure_eigencycle_value_a
clear
tool_figure_eigencycle_value_a
S=solve(V_F)
S.x1
S.x2
S.x3
S.x4
tool_figure_eigencycle_value_a
eval(S)
cell2str(S)
cellstr(S)
cellstr(S.x2)
x1(S.x1)
x1=S.x1
x1=S.x1(2)
tool_figure_eigencycle_value_a
real(m)
imag(m)
real(n)
cross([real(m) imag(m) 0],[real(n) imag(n) 0]
cross([real(m) imag(m) 0],[real(n) imag(n) 0])
tool_figure_eigencycle_value_a
%-- 2021/5/17 8:46 --%
tool_figure_eigencycle_value_a
norm(eval(eygencycleset_rr)
norm(eval(eygencycleset_rr))
tool_figure_eigencycle_value_a
norm(eval(eygencycleset_rr))
eygencycleset_rr.*((3*a)/4 + 1/4)
eygencycleset_rr./((3*a)/4 + 1/4)
eygencycleset_rr./((3*a)/4 + 1/4)/pi
norm(eygencycleset_rr./((3*a)/4 + 1/4)/pi)
((a/2 - 1/2)^2 + 2(a/2 + 1/2)^2 + 2)^(1/2)
((a/2 - 1/2)^2 + 2*(a/2 + 1/2)^2 + 2)^(1/2)
((a/2 - 1/2)^2 + 2*(a/2 + 1/2)^2 + 2)
simplify((a/2 - 1/2)^2 + 2*(a/2 + 1/2)^2 + 2)
eygencycleset_rr./((3*a)/4 + 1/4)/pi
x = eygencycleset_rr./((3*a)/4 + 1/4)/pi
y=x'*x
y=sqrt(x'*x)
latex(((a/2 - 1/2)^2 + 2*(a/2 + 1/2)^2 + 2)^(1/2))
%-- 2021/5/17 22:36 --%
function r=testsave2Dfig()
% data_ori = [  % Sam28	dimX	dimY	t8i(4)	t4i1	t4i2	Oneill(7)	Binmore	AIBT	ATBI(10)	AIBI	ATBT	total	E8i	E4i	Classid(15)
% 12	1	2	-4.81E-17	0	0	8	3.438888885	-23.4	-4	4	19	7.038888885	0	0	0
% 13	1	3	-4.81E-17	0	0	-29	-3.288888888	10.2	-1	-21.6	-4	-48.68888889	0	0	0
% 14	1	4	-4.81E-17	0	0	21	-0.149999998	13.2	5	17.6	-15	41.65	0	0	0
% 15	1	5	-1.1782	0	0	-99	-28.76111112	-159.2	-98.6	-92.6	-117	-595.1611111	1	0	2
% 16	1	6	0.3927	0	0	24.6	20.27222221	53.8	56.6	18	48.6	221.8722222	1	0	2
% 17	1	7	0.3927	0	0	55.4	3.566666669	61.2	19	36	34.4	209.5666667	1	0	2
% 18	1	8	0.3927	0	0	19	4.922222222	44.2	23	38.6	34	163.7222222	1	0	2
% 23	2	3	0	-0.071439963	0.036010289	29	4.955555547	-27.2	2	19.2	16	43.95555555	0	1	1
% 24	2	4	0	0.071549996	-0.036019622	-21	-1.516666671	3.8	-6	-15.2	3	-36.91666667	0	1	1
% 25	2	5	0.3927	0	0	25	6.749999992	61.8	14.2	35.2	42	184.95	1	0	2
% 26	2	6	-0.1309	0.900211091	0.077437117	-2.2	-0.749999997	10.6	-33.2	4	-20.2	-41.75	1	1	3
% 27	2	7	-0.1309	-0.143086216	-0.27330414	-3.8	-2.31111111	-30.2	13	-28	-33.8	-85.11111111	1	1	3
% 28	2	8	-0.1309	0.741983761	0.195831955	-19	-3.688888886	-42.2	6	-11.2	12	-58.08888889	1	1	3
% 34	3	4	0	-0.071441185	0.036059535	-1.33E-15	1.666666663	-17	1	-2.4	12	-4.733333337	0	1	1
% 35	3	5	0.3927	0	0	28	15.11666666	45.2	34.2	33.8	30	186.3166667	1	0	2
% 36	3	6	-0.1309	-0.143086216	-0.27330414	-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2	-74.73333332	1	1	3
% 37	3	7	-0.1309	0.028412595	0.981336815	-8.8	1.12222222	-13	-11	5.8	5.2	-20.67777778	1	1	3
% 38	3	8	-0.1309	-0.131098343	-0.707913246	4	-4.505555552	-24	-17	-25.4	-24	-90.90555555	1	1	3
% 45	4	5	0.3927	0	0	46	6.894444446	52.2	50.2	23.6	45	223.8944444	1	0	2
% 46	4	6	-0.1309	-0.756919893	0.195831955	0.8	-7.788888883	-56.2	-17.2	-7.8	-17.2	-105.3888889	1	1	3
% 47	4	7	-0.1309	0.114632307	-0.707913246	-42.8	-2.377777774	-18	-21	-13.8	-5.8	-103.7777778	1	1	3
% 48	4	8	-0.1309	-0.610696192	0.511996947	-4	3.27222222	22	-12	-2	-22	-14.72777778	1	1	3
% 56	5	6	-4.81E-17	0	0	-16.4	2.027777784	-3.6	4.4	-3.2	5.6	-11.17222222	0	0	0
% 57	5	7	-4.81E-17	0	0	23.4	-0.955555556	10.8	13.8	-11.2	-21.6	14.24444444	0	0	0
% 58	5	8	-4.81E-17	0	0	-7	-1.072222226	-7.2	-18.2	14.4	16	-3.072222226	0	0	0
% 67	6	7	0	-0.071439963	0.036010289	-8.6	1.522222224	-9.8	1.2	11	16.4	11.72222222	0	1	1
% 68	6	8	0	0.16585302	-0.036019622	-7.8	0.505555554	6.2	3.2	-14.2	-10.8	-22.89444445	0	1	1
% 78	7	8	0	0.032521274	0.036059535	14.8	0.56666667	1	15	-0.2	-5.2	25.96666667	0	1	1
% ];
data_ori =[%Sam28	dimX	dimY	t8i(4)	t4i1	t4i2	Oneill(7)	Binmore	AIBT	ATBI(10)	AIBI	ATBT	LP11?¨total)	E8i	E4i	Classid(15)
12	1	2	0.00E+00	0	0	8	3.438888885	-23.4	-4	4	19	3.6	0	0	0
13	1	3	0.00E+00	0	0	-29	-3.288888888	10.2	-1	-21.6	-4	-45.4	0	0	0
14	1	4	0.00E+00	0	0	21	-0.149999998	13.2	5	17.6	-15	41.8	0	0	0
15	1	5	-0.7503335	0	0	-99	-28.76111112	-159.2	-98.6	-92.6	-117	-566.4	1	0	2
16	1	6	0.2498565	0	0	24.6	20.27222221	53.8	56.6	18	48.6	201.6	1	0	2
17	1	7	0.2498565	0	0	55.4	3.566666669	61.2	19	36	34.4	206	1	0	2
18	1	8	0.2498565	0	0	19	4.922222222	44.2	23	38.6	34	158.8	1	0	2
23	2	3	0	-0.2887246	0.2887246	29	4.955555547	-27.2	2	19.2	16	39	0	1	1
24	2	4	0	0.2887246	-0.2887246	-21	-1.516666671	3.8	-6	-15.2	3	-35.4	0	1	1
25	2	5	0.2498565	0	0	25	6.749999992	61.8	14.2	35.2	42	178.2	1	0	2
26	2	6	-0.0832855	0.3334082	0.3334082	-2.2	-0.749999997	10.6	-33.2	4	-20.2	-41	1	1	3
27	2	7	-0.0832855	-0.1665132	-0.1665132	-3.8	-2.31111111	-30.2	13	-28	-33.8	-82.8	1	1	3
28	2	8	-0.0832855	-0.1665132	-0.1665132	-19	-3.688888886	-42.2	6	-11.2	12	-54.4	1	1	3
34	3	4	0	-0.2887246	0.2887246	-1.33E-15	1.666666663	-17	1	-2.4	12	-6.4	0	1	1
35	3	5	0.2498565	0	0	28	15.11666666	45.2	34.2	33.8	30	171.2	1	0	2
36	3	6	-0.0832855	-0.1665132	-0.1665132	-23.2	-11.73333332	-8.2	-6.2	-14.2	-11.2	-63	1	1	3
37	3	7	-0.0832855	0.3334082	0.3334082	-8.8	1.12222222	-13	-11	5.8	5.2	-21.8	1	1	3
38	3	8	-0.0832855	-0.1665132	-0.1665132	4	-4.505555552	-24	-17	-25.4	-24	-86.4	1	1	3
45	4	5	0.2498565	0	0	46	6.894444446	52.2	50.2	23.6	45	217	1	0	2
46	4	6	-0.0832855	-0.1665132	-0.1665132	0.8	-7.788888883	-56.2	-17.2	-7.8	-17.2	-97.6	1	1	3
47	4	7	-0.0832855	-0.1665132	-0.1665132	-42.8	-2.377777774	-18	-21	-13.8	-5.8	-101.4	1	1	3
48	4	8	-0.0832855	0.3334082	0.3334082	-4	3.27222222	22	-12	-2	-22	-18	1	1	3
56	5	6	0.00E+00	0	0	-16.4	2.027777784	-3.6	4.4	-3.2	5.6	-13.2	0	0	0
57	5	7	0.00E+00	0	0	23.4	-0.955555556	10.8	13.8	-11.2	-21.6	15.2	0	0	0
58	5	8	0.00E+00	0	0	-7	-1.072222226	-7.2	-18.2	14.4	16	-2	0	0	0
67	6	7	0	-0.2887246	0.2887246	-8.6	1.522222224	-9.8	1.2	11	16.4	10.2	0	1	1
68	6	8	0	0.2887246	-0.2887246	-7.8	0.505555554	6.2	3.2	-14.2	-10.8	-23.4	0	1	1
78	7	8	0	-0.2887246	0.2887246	14.8	0.56666667	1	15	-0.2	-5.2	25.4	0	1	1
];
TotalPeriod = [2625 1950 3000 3000 2375 2375 sum([2625 1950 3000 3000 2375 2375])];
data=data_ori; for kt=7:13;data(:,kt)=data_ori(:,kt)./TotalPeriod(kt-6);end
%%%%% for J=0:3; subplot(2,2,J+1);a=data(find(data(:,16)==J & data(4,4)<1),:);boxplot(a(:,[7 9:13])); ylim([-.1 .1]*.2);axis square;grid on;end
MinMax=[13.65	1.05	3.15	-1.23711	-0.794765888	-0.743308908	-103.95	-30.19916668	-167.16	-103.53	-97.23	-122.85	-624.9191667
81.9	7.35	8.4	0.412335	0.945221646	1.030403656	58.17	21.28583332	64.89	59.43	40.53	51.03	235.0891666
]';
MinMax=[min(data);max(data) ]';
% field_id={ 'Sam28'	'dimX'	'dimY'	't8i'	't4i1'	't4i2'	'ONeill'	'Binmore'	'AIBT'	'ATBI'	'AIBI'	'ATBT'	'Total'	'E8i'	'E4i'	'Class id'};
field_id={ 'Sam28'	'dimX'	'dimY'	't8i'	't4i1'	't4i2'	'Experiment ($L_O$)'	'Experiment ($L_B$)'	'Experiment ($L_{IT}$)'	...
'Experiment ($L_{TI}$)'	'Experiment ($L_{II}$)'	'Experiment ($L_{TT}$)'	'Experiment ($L_{P11}$)'	'E8i'	'E4i'	'Class id'};
titel_id={ '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(a)'};
for Yfig_id=7:13
if Yfig_id == 13; data(:,Yfig_id) = LP11(data_ori,TotalPeriod);end
fig_name_saved = strcat('Yfig',num2str(Yfig_id))
figureFilePath = 'D:/M/grd_co/twobytwo-2Dig/';xlab ='Theory ($\sigma_{8i}$)', ylab = field_id(Yfig_id)
class_color=[  1 .5 1 ; .13 .52 .13 ; .3 0 1 ;1 0 0 ;];
class_label=[ 'd' 'X' 'o' 's'];
class_c = class_color(data(:,16)+1,:)
class_l = class_label(data(:,16)+1)
ax2 = scatter(data(:,4),data(:,Yfig_id), 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
line([-1.5,1.5],[0,0], 'LineWidth',2,  'color',[0.8 0.8 0.8]);hold on
line([0,0], [MinMax(Yfig_id,1),MinMax(Yfig_id,2)]*1.2,'LineWidth',2,  'color',[0.8 0.8 0.8]);hold on
for n=1:28
ax1 = scatter(data(n,4),data(n,Yfig_id), 129, class_l(n), 'MarkerEdgeColor',[class_c(n,1) class_c(n,2) class_c(n,3)], 'LineWidth',2);hold on
end
xlim([-1 0.5]);ylim([MinMax(Yfig_id,1),MinMax(Yfig_id,2)]*1.2);text(-0.6, MinMax(Yfig_id,2), titel_id(Yfig_id - 6),'fontsize',20);
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
end
%% R3 figure 4i dependence alpha
clf; n=3;
s9 = data_ori(find(data_ori(:,16) == n), :)
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2)
expO9 = sum(da')'/dt
the4i2 = s9(:,6)
ax2 = scatter(the4i2, expO9, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, 's', 'MarkerEdgeColor',[1 0 0], 'LineWidth',2);
line([-0.2 0.4],[0,0], 'LineWidth',2,  'color',[0.8 0.8 0.8]);hold on
line([0,0], [-0.008,0.001],'LineWidth',2,  'color',[0.8 0.8 0.8]);hold on
xlab ='Theory ($\sigma_{\alpha}$)'; ylab = 'Experiment ($L_{P11}$)';
fig_name_saved = 'figalpha';
text(-0.01, 0.001*0.6,'(b)','fontsize',20);
%                 text(-0.1, -0.001*1.1,'(b)','fontsize',20);
xlim([-0.2 0.4]);ylim([-0.008,0.001]);
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
%% beta
clf; n=1;
s9 = data_ori(find(data_ori(:,16) == n), :)
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2)
expO9 = sum(da')'/dt
the4i2 = s9(:,6)
ax2 = scatter(the4i2, expO9, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(the4i2, expO9, 129, 'd', 'MarkerEdgeColor',[1 .5 1], 'LineWidth',2);
line([-0.4 0.4],[0,0], 'LineWidth',2,  'color',[0.8 0.8 0.8]);hold on
line([0,0], [-0.003,0.003],'LineWidth',2,  'color',[0.8 0.8 0.8]);hold on
xlab ='Theory ($\sigma_{\beta}$)'; ylab = 'Experiment ($L_{P11}$)';
fig_name_saved = 'figbeta';
text(-0.1, 0.001*2.7,'(c)','fontsize',20);
%                 text(-0.4, -0.001*1.1,'(b)','fontsize',20);
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
%%    multiple
alpha=(data_ori(:,6)+data_ori(:,5))./2;
beta=(data_ori(:,6)-data_ori(:,5))./2;
mul8i_alpha_beta=0.0565190.*data_ori(:,4)+0.0057181.*alpha+0.0045221.*beta
all_lp11= data_ori(:,13)./dt
ax2 = scatter(mul8i_alpha_beta, all_lp11, 1, 'MarkerEdgeColor','w')
kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
scatter(mul8i_alpha_beta, all_lp11, 129, 's', 'MarkerEdgeColor',[1 0 0], 'LineWidth',2);
xlab ='Theory(Multiple regression) '; ylab = 'Experiment ($L_{P11}$)';
fig_name_saved = 'figMultiple';
text(-0.1, -0.001*1.3,'(b)','fontsize',20);
%                 text(-0.1, -0.001*1.1,'(b)','fontsize',20);
plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
% %%
%
%  alpha=(data_ori(:,6)+data_ori(:,5))./2;
% beta=(data_ori(:,6)-data_ori(:,5))./2;
%  data_ori(:,17)=alpha;
%         data_ori(:,18)=beta;
% data_ori11_9=data_ori(find(abs(data_ori(:,17))>0.0001),:);
%             ax2 = scatter(data_ori11_9(:,17), data_ori11_9(:,12), 1, 'MarkerEdgeColor','w')
%                 kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
%             scatter(data_ori11_9(:,17),  data_ori11_9(:,12), 129, 's', 'MarkerEdgeColor',[1 0 0], 'LineWidth',2);
%                 xlab ='Theory ($\sigma_{alpha}$)'; ylab = 'Experiment ($L_{P11}$)';
%                 fig_name_saved = 'figalpha';
%                 text(-0.1, -0.001*1.3,'(b)','fontsize',20);
% %                 text(-0.1, -0.001*1.1,'(b)','fontsize',20);
%             plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
end
function rLP11 = LP11(data_ori11,TotalPeriod)
clf; n=3;
s9 = data_ori11;
da=s9(:,[7 9:12])
dt = sum(TotalPeriod(1:6))-TotalPeriod(2);
rLP11 = sum(da')'/dt
%         the8i = s9(:,4)
%                 ax2 = scatter(the8i, expO9, 1, 'MarkerEdgeColor','w')
%                 kk = lsline ; set(kk,'linewidth',5,'LineStyle',':');hold on;
%             scatter(the8i, expO9, 129, 's', 'MarkerEdgeColor',[1 0 0], 'LineWidth',2);
%                 xlab ='Theory ($\sigma_{8i}$)'; ylab = 'Experiment ($L_{P11}$)';
%                 fig_name_saved = 'fig8i';
%             plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath);
end
function plotfilesaveTo=save2Dfig20201123(xlab, ylab, fig_name_saved, figureFilePath)
xlabel(xlab,'FontSize',25,'interpreter','latex');
ylabel(ylab,'FontSize',25,'interpreter','latex');  axis square;
set(gca,'FontSize',20,'linewidth',2);hold on;
box on;
grid off
set(gcf,'papersize',[15 15],'paperposition',[0,0,15,15]);hold on;
%           set(gca,'xlim',[100 400],'XTick',100:50:400)
%           set(gca,'ylim',[8 24],'YTick',8:4:24)
plotfilesaveTo =strcat( figureFilePath, fig_name_saved,'.png')
saveas(gcf, plotfilesaveTo);
plotfilesaveTo =strcat( figureFilePath, fig_name_saved,'.pdf')
saveas(gcf, plotfilesaveTo);
plotfilesaveEPS =strcat( figureFilePath, fig_name_saved ,'.eps')
saveas(gcf, plotfilesaveEPS,'psc2');
close
end
WangYao2105_hypofine
set(gca,'tickdir','out');
set(gca,'tickdir','in');
set(gca,'ticklabeldir','in');
set(gca,'labeldir','in');
set(gca,'tickdir','in');
%-- 2021/5/18 14:28 --%
[v d]=eig([4 9 -13
-5 -9 14
1 0 -1])
load('C:\Users\Think\Desktop\2021毕设\周名铉\eigen.mat')
v = [     -0.0015190517410030831822535279648902i
0.49213660213971016988350157594816i
-0.00000041580361041727194178560523229358i
0.00011945162452490785029171263968735i
-0.031520296874740631362638818176905i
-0.45909583424060977839689323358007i
-0.0000077303670440890614331393569463273i
-0.00011272473722707845863278390380825i
-0.001130913
-0.521868881
-0.000169506
0.523152192
-3.85839E-07
-0.000170198
-5.50584E-08
0.000187747
];
plot(v)
labda= - 1.0 - 5.5788083301784925639030018689538i
t=0:0.01:2; w=v.*exp(labda*t)*1
t=0:0.01:2; w=v*exp(labda*t)*1
plot(w(1,:))
xa=real(w(2,:))
xb=real(w(10,:))
plot(xa,xb)
clf
plot(xa,xb)
%-- 2021/5/19 12:00 --%
lissa2Plot
%-- 2021/5/19 12:33 --%
lissa2Plot
%-- 2021/5/20 15:08 --%
[c*(am + i*bm)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))
]
sysm am bm an bn a_lambda b_lambda  t c real
[c*(am + i*bm)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))
]
syms am bm an bn a_lambda b_lambda  t c real
[c*(am + i*bm)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))
]
real(ans)
[c*(an + i*bn)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))
]
syms am bm an bn a_lambda b_lambda  t c real
x_m = [c*(am + i*bm)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))]
x_n = [c*(an + i*bn)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))]
r_x_m = real(x_m)
r_x_n = real(x_n)
ddt_r_x_m = diff(r_x_m,'t')
ddt_r_x_n = diff(r_x_n,'t')
cross([r_x_m  r_x_n  0], [ddt_r_x_m  ddt_r_x_n  0])
instanL = cross([r_x_m  r_x_n  0], [ddt_r_x_m  ddt_r_x_n  0])
simplify(instanL(3))
clear
instanL
simplify(r_x_m)
instanL
simplify(ddt_r_x_m)
%-- 2021/5/20 18:21 --%
x_m = [c*(am + i*bm)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))]
syms am bm an bn a_lambda b_lambda  t c real
x_m = [c*(am + i*bm)*(cos((a_lambda + i*b_lambda)*t) + i*sin((a_lambda + i*b_lambda)*t))]
latex(x_m)
%-- 2021/5/20 18:50 --%
instanL
%-- 2021/5/20 21:24 --%
toms t
p = tomPhase('p', t, 0, 0.2, 20);
setPhase(p);
tomStates x1 x2 x3 x4 x5 x6 x7
tomControls u1 u2 u3 u4
x = [x1; x2; x3; x4; x5; x6; x7];
u = [u1; u2; u3; u4];
x0i = [0.1883;0.2507;0.0467;0.0899;0.1804;0.1394;0.1046];
x0 = icollocate({x1==x0i(1),x2==x0i(2),x3==x0i(3),x4==x0i(4),x5==x0i(5),x6==x0i(6),x7==x0i(7)});
% Box constraints and boundary
uL = zeros(4,1); uU = [20;6;4;20];
cbb = {collocate(uL <= u <= uU)
initial(x == x0i)};
% ODEs and path constraints
q = u(1)+u(2)+u(4);
ceq = collocate({
dot(x1) == u4-q.*x1-17.6*x1.*x2-23*x1.*x6.*u3;
dot(x2) == u1-q.*x2-17.6*x1.*x2-146*x2.*x3;
dot(x3) == u2-q.*x3-73*x2.*x3;
dot(x4) == -q.*x4+35.2*x1.*x2-51.3*x4.*x5;
dot(x5) == -q.*x5+219*x2.*x3-51.3*x4.*x5;
dot(x6) == -q.*x6+102.6*x4.*x5-23*x1.*x6.*u3;
dot(x7) == -q.*x7+46*x1.*x6.*u3});
% Objective
objective = integrate(-(5.8*(q.*x1-u4)-3.7*u1-4.1*u2+...
q.*(23*x4+11*x5+28*x6+35*x7)-5.0*u3.^2-0.099));
cd c:\tomlab\
startup
vanDerPol
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_2.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_10.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_2.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_imag2_10.fig',1)
abs([      -0.15494774856142694532083465380928i
-0.25734160533365550950802612982817i
0.000027317005470930115908124106269321i
-0.0000050568892148258081028807735565442i
0.3912257035863691569268080378674i
0.021025687404089776042930036777495i
0.000014804893778527784752587858957484i
0.00000089789458888976656487780088575144i
-0.127499826
-0.534511527
-0.078745556
0.62377761
-0.024262227
-0.07878076
-0.008328933
0.228351219
])
abs([ 0.606051539
-0.40620472
-0.000242122
-3.86071E-05
-0.151385831
-0.048153343
-2.35437E-05
-3.37201E-06
0.051357551702254107007554985325012i
-0.31786992125280511494700471894338i
0.13216040987084861188840127453643i
-0.43518118974249019295518327998378i
0.033583000145554480979481442348376i
0.13209286814164252156558587059153i
0.075817261105617834279249694287044i
0.32804002002937775218191473183876i
])
abs([        0.00042678400803120157910732645138997i
0.30183258943764705159114112169191i
-0.00000037165299980714718052727347513137i
0.000000056193815335362514509025941837735i
-0.29074440734485148556486792756242i
-0.011514547046458234918945249980749i
-0.000000099015704579470645038307596445364i
-0.0000000045794794814311242140450003121288i
-0.270444509
-0.531092115
-0.003522008
0.665822298
-0.02578894
-0.003521952
0.010722104
0.157825123
])
abs([ 0.589780368
-0.376979682
-2.96544E-06
-3.67031E-07
-0.18034799
-0.032449162
-1.81363E-07
-2.00144E-08
0.17158616616253713682012265028865i
-0.54403545845985707616661175654067i
0.16642709566674112294818577150941i
-0.25000123413009188460430106000159i
0.072946882480634329827369630581566i
0.16642647017352215499578909349551i
0.061610703538245025098327695862537i
0.15503937456826919108111797480459i
])
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_2.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\8_exp2_imag1_12.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\8_exp2_real2_10.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_imag2_12.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\8_exp2_imag2_10.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_imag2_12.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\10_exp11_real2_12.fig',1)
uiopen('C:\Users\Think\Desktop\2021毕设\周名铉\20210520ZMX\fig\6_exp13_imag2_10.fig',1)
%-- 2021/5/21 9:27 --%
instanL
L2 = cross([real(r_x_m)  real(r_x_n)  0], [real(ddt_r_x_m)  real(ddt_r_x_n)  0])
simplify(L2)
L3 = cross([(r_x_m)  (r_x_n)  0], [(ddt_r_x_m)  (ddt_r_x_n)  0])
clear
instanL
x_m = c*(eta_m_real + i*eta_m_imag)* ...
(cos((lambda_real + i*lambda_imag)*t) ...
+ i*sin((lambda_real + i*lambda_imag)*t))
x_n = c*(eta_n_real + i*eta_n_imag)*(cos((lambda_real + ...
i*lambda_imag)*t) + i*sin((lambda_real + ...
i*lambda_imag)*t))
r_x_m = real(x_m)
r_x_n = real(x_n)
% r_x_m =  (x_m)
% r_x_n =  (x_n)
ddt_r_x_m = diff(r_x_m,'t')
ddt_r_x_n = diff(r_x_n,'t')
instantaneousL = cross([r_x_m  r_x_n  0], [ddt_r_x_m  ddt_r_x_n  0])
L = simplify(instantaneousL(3))
instanL
x_m = c*(eta_m_real + i*eta_m_imag)* ...
(cos((lambda_real + i*lambda_imag)/i*t) ...
+ i*sin((lambda_real + i*lambda_imag)/i*t))
x_n = c*(eta_n_real + i*eta_n_imag)*(cos((lambda_real + ...
i*lambda_imag)/i*t) + i*sin((lambda_real + ...
i*lambda_imag)/i*t))
r_x_m = real(x_m)
r_x_n = real(x_n)
% r_x_m =  (x_m)
% r_x_n =  (x_n)
ddt_r_x_m = diff(r_x_m,'t')
ddt_r_x_n = diff(r_x_n,'t')
instantaneousL = cross([r_x_m  r_x_n  0], [ddt_r_x_m  ddt_r_x_n  0])
L = simplify(instantaneousL(3))
instanL
replace(L, 'lambda_imag','d')
replace(sym2str(L), 'lambda_imag','d')
replace(sim2str(L), 'lambda_imag','d')
replace(string(L), 'lambda_imag','d')
findsym
display(L)
a=display(L)
h=display(L)
display(L)
sym2cell(L)
a=sym2cell(L)
a(1,1)
a(1,1),tex
a(1,1).text
a.text
a
cellstr(a)
symstr()
syms(L)
instanL
%-- 2021/5/21 23:40 --%
int(x)
syms x
int(x)
int(x,[2 4])
syms x r t real
syms x r t T real; int(exp(2*r*t),'t',[0 T])
%-- 2021/5/22 16:26 --%
ShujieLandscape1228
eigen_vector
syms c1 c2 c3 a real
X = c1*[1; 0; 1; 0] + ...
c2*[0; a; 0; 1] + ...
c3*
ShujieLagrange
eigen_value
syms c1 c2 c3 a t real
eigen_vector = [
[ 1, 0, a*(1/4 + 3i/4) - 1/4 + 1i/4, a*(1/4 - 3i/4) - 1/4 - 1i/4]
[ 0, a,                 - a/2 - 1/2,                 - a/2 - 1/2]
[ 1, 0, a*(1/4 - 3i/4) - 1/4 - 1i/4, a*(1/4 + 3i/4) - 1/4 + 1i/4]
[ 0, 1,                           1,                           1]]
eigen_value = [
[ -1,  0,   0,  0]
[  0, -1,   0,  0]
[  0,  0, -1i,  0]
[  0,  0,   0, 1i]]
X = c1*eigen_vector(:,1) * exp(eigen_value(1,1)*t) + ...
c2*eigen_vector(:,2) * exp(eigen_value(2,2)*t) + ...
c3*eigen_vector(:,3) * exp(eigen_value(3,3)*t)
ShujieLagrange
sum(X)
simplify(sum(X))
simplify(X/sum(X))
real(simplify(X/sum(X)))
simplify(real(simplify(X/sum(X))))
Xn = simplify(real(simplify(X/sum(X))))
Xn(3)-Xn(1)
simplify(Xn(3)-Xn(1))
ShujieLagrange
mean_U = Xn' * Payoff_vector_field_F
simplify(mean_U)
simplify(simplify(mean_U))
V_F = mean_U
D_V_F = [diff(V_F,'x1') diff(V_F,'x2') diff(V_F,'x3') diff(V_F,'x4')]
D_V_F = [diff(V_F,'c1') diff(V_F,'c2') diff(V_F,'c3')]
simplify(D_V_F)
simplify(simplify(D_V_F)')
ShujieLagrange
V_F = mean_U
D_V_F = [diff(V_F,'c1') diff(V_F,'c2') diff(V_F,'c3') diff(V_F,'c4')]
simplify(D_V_F)
simplify(simplify(simplify(D_V_F))))
simplify(simplify(simplify(D_V_F)))
simplify(simplify(simplify(D_V_F))')
Eq=simplify(simplify(simplify(D_V_F))')
Eq(3)-Eq(4)
a=4; t=0.001; eval(Eq)
Solve(eval(Eq))
solve(eval(Eq))
S = solve(eval(Eq))
S.c1
S.c2
S.c3
S.c4
Eq2=Eq(1:3)
Eq2=subs(Eq(1:3),'c4','c3')
S = solve(eval(Eq2))
S.c1
S.c2
S.c3
eval(Eq2)
double(eval(Eq2))
format short
eval(Eq2)
for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1; rr=[rr; c1 c2 c3
eval(mean_U)
for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1; rr=[rr; c1 c2 c3 eval(mean_U))];end
for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1; rr=[rr; c1 c2 c3 eval(mean_U))];end;end;end
for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1; rr=[rr; c1 c2 c3 eval(mean_U)];end;end;end
rr=[];for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1; rr=[rr; c1 c2 c3 eval(mean_U)];end;end;end
vpa(mean_U)
eval(mean_U)
vpa(eval(mean_U))
rr=[];for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1; rr=[rr; c1 c2 c3 vpa(eval(mean_U))];end;end;end
c4=c3
rr=[];for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1;c4=c3; rr=[rr; c1 c2 c3 vpa(eval(mean_U))];end;end;end
vpa(eval(mean_U))
eval(mean_U)
rr=[];for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1;c4=c3; rr=[rr; c1 c2 c3 vpa(eval(mean_U)];end;end;end
rr=[];for c1=0:0.1:1; for c2=0:0.1:1; for c3=0:0.1:1;c4=c3; rr=[rr; c1 c2 c3 eval(mean_U)];end;end;end
plot(rr(:,end))
rr=[];for c1=0:0.1:1; for c2=0:0.1:1; for c3=(1-c1-c2)/2;c4=c3; rr=[rr; c1 c2 c3 eval(mean_U)];end;end
rr=[];for c1=0:0.1:1; for c2=0:0.1:1;  c3=(1-c1-c2)/2;c4=c3; rr=[rr; c1 c2 c3 eval(mean_U)];end;end
plot(rr(:,end))
rr=[];for c1=0:0.1:1; for c2=0:0.1:1;  c3=sqrt(1-c1^2-c2^2)/2;c4=c3; rr=[rr; c1 c2 c3 eval(mean_U)];end;end
plot(rr(:,end))
rr=[];for c1=0:0.1:1; for c2=0:0.1:sqrt(1-c1^2);  c3=sqrt(1-c1^2-c2^2)/2;c4=c3; rr=[rr; c1 c2 c3 eval(mean_U)];end;end
plot(rr(:,end))
rr=[];for c1=0:0.1:1; for c2=0:0.1:sqrt(1-c1^2);  c3=sqrt(1-c1^2-c2^2)/2;c4=c3; rr=[rr; c1 c2 c3 real(eval(mean_U))];end;end
plot(rr(:,end))
rr=[];for c1=0:0.01:1; for c2=0:0.01:sqrt(1-c1^2);  c3=sqrt(1-c1^2-c2^2)/2;c4=c3; rr=[rr; c1 c2 c3 real(eval(mean_U))];end;end
plot(rr(:,end))
histogram(rr(:,end))
histogram(rr(:,end))
histogram(rr(:,end))
plot(rr(:,end))
a=0.25
rr=[];for c1=0:0.01:1; for c2=0:0.01:sqrt(1-c1^2);  c3=sqrt(1-c1^2-c2^2)/2;c4=c3; rr=[rr; c1 c2 c3 real(eval(mean_U))];end;end
plot(rr(:,end))
plot3(rr(1:2,end))
plot3(rr([1 2 4],:))
plot3(rr([1],:),rr([ 2 ],:),rr([ 4],:))
hist3(rr([1],:),rr([ 2 ],:),rr([ 4],:))
pcolor(rr)
ShujieLagrange
clear
ShujieLagrange
eigencycle
ShujieLagrange
eigencycle
ShujieLagrange
ShujieLandscape1228
eigencycle = [];
for m=1:4-1
for n=m+1:4
L0= cross([real(eigen_vector(3,m)) imag(eigen_vector(3,m)) 0 ], ...
[real(eigen_vector(3,n)) imag(eigen_vector(3,n)) 0 ])
eigencycle = [eigencycle; L0(3)];
end
end
eigencycle
ShujieLagrange
syms c1 c2 c3 c4 a t real
eigen_vector = [
[ 1, 0, a*(1/4 + 3i/4) - 1/4 + 1i/4, a*(1/4 - 3i/4) - 1/4 - 1i/4]
[ 0, a,                 - a/2 - 1/2,                 - a/2 - 1/2]
[ 1, 0, a*(1/4 - 3i/4) - 1/4 - 1i/4, a*(1/4 + 3i/4) - 1/4 + 1i/4]
[ 0, 1,                           1,                           1]]
eigen_value = [
[ -1,  0,   0,  0]
[  0, -1,   0,  0]
[  0,  0, -1i,  0]
[  0,  0,   0, 1i]]
X = c1*eigen_vector(:,1) * exp(eigen_value(1,1)*t) + ...
c2*eigen_vector(:,2) * exp(eigen_value(2,2)*t) + ...
c3*eigen_vector(:,3) * exp(eigen_value(3,3)*t)
Xn = simplify(real(simplify(X/sum(X))))
dXn13 = simplify(Xn(3)-Xn(1))
ShujieLagrange
mean_U
simplify(mean_U)
int(simplify(mean_U),'t',[0,tau])
syms tau; int(simplify(mean_U),'t',[0,tau])
limit(int(simplify(mean_U),'t',[0,tau]),'tau',0)
limit(int(simplify(mean_U),'t',[0,tau])/tau,'tau',0)
simplify(ans)
meanUtau = limit(int(simplify(mean_U),'t',[0,tau])/tau,'tau',0)
simplify(meanUtau)
meanUtau = simplify(limit(int(simplify(mean_U),'t',[0,tau])/tau,'tau',0))
ShujieLagrange
L_motion = eigencycle*c3^2*imag(eigen_value(3,3))
L1= sqrt(L_motion'*L_motion)
L1= simplify(sqrt(L_motion'*L_motion))
L1= sqrt(L_motion'*L_motion)
latex(L1)
latex(L1/c3^2)
latex(simplify(L1/c3^2))
simplify(L1/c3^2)
ShujieLagrange
imag(eigen_value(3,3)
imag(eigen_value(3,3))
ShujieLagrange
L1/c3
aa = (c3^4*((3*a)/4 + 1/4)^2 + c3^4*(a/2 + 1/2)^2*((3*a)/4 + 1/4)^2 + 2*c3^4*(a/4 - 1/4)^2*((3*a)/4 + 1/4)^2)^(1/2) /c3
simplify(aa)
aa = (c3^4*((3*a)/4 + 1/4)^2 + c3^4*(a/2 + 1/2)^2*((3*a)/4 + 1/4)^2 + 2*c3^4*(a/4 - 1/4)^2*((3*a)/4 + 1/4)^2)
simplify(aa)
solve('3*a^2 + 2*a + 11')
solve('3*a^2 + 2*a + 11=0',a)
solve('3*a^2 + 2*a + 11=0','a')
a=0.1:0.1:2;plot(3*a^2 + 2*a + 11)
a=0.1:0.1:2;plot(a,3*a^2 + 2*a + 11)
a=0.1:0.1:2;y= 3*a.^2 + 2*a + 11
plot(a,y)
a=-2:0.1:2;y= 3*a.^2 + 2*a + 11
plot(a,y)
syms a
solve('3*a^2 + 2*a + 11=0','a')
syms a imag
solve('3*a^2 + 2*a + 11=0','a')
imag(eigen_value(3,3))
L1 = sqrt(L_motion'*L_motion)
L1 = sqrt(L_motion.*L_motion)
L1 = sqrt(sum(L_motion.*L_motion))
ShujieLagrange
clear
ShujieLagrange
simplify(W)
simplify(W,'c3 ~= 0')
latex(W)
F = -(c3*exp(t)*sin(t)*(3*a + 1))/(4*c1 + 2*c2 + 2*a*c2)
simplify(limit(int(F,'t',[0,tau])/tau,'tau',0));
simplify(limit(int(F,'t',[0,tau])/tau,'tau',0))
simplify(limit(int(F,'t',[0,tau])/tau,'tau',1))
simplify(limit(int(F,'t',[0,tau])/tau,'tau',pi))
%-- 2021/5/22 22:11 --%
ShujieLagrange
V_F = W
D_V_F = [diff(V_F,'c1') diff(V_F,'c2') diff(V_F,'c3') ]
simplify(D_V_F)
simplify(D_V_F)'
[(c2 - 2*c1 + c3 + a*c2 - a*c3));(4*c1 - c3 + a*c3)*(c2 - 2*c1 + a*c2);(c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4 + (c2*(3*a + 1)*(a - 1))]
(c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4 + (c2*(3*a + 1)*(a - 1))
[(c2 - 2*c1 + c3 + a*c2 - a*c3));(4*c1 - c3 + a*c3)*(c2 - 2*c1 + a*c2);  (c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4 + (c2*(3*a + 1)*(a - 1))]
(4*c1 - c3 + a*c3)*(c2 - 2*c1 + a*c2);
(c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4 + (c2*(3*a + 1)*(a - 1))
[ (c2 - 2*c1 + c3 + a*c2 - a*c3))   ;    ; c2*(3*a + 1)*(a - 1) + (c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4]
[ (c2 - 2*c1 + c3 + a*c2 - a*c3)   ;    ; c2*(3*a + 1)*(a - 1) + (c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4]
[ (c2 - 2*c1 + c3 + a*c2 - a*c3)   ;   (4*c1 - c3 + a*c3)*(c2 - 2*c1 + a*c2); ; c2*(3*a + 1)*(a - 1) + (c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4]
Q = [ (c2 - 2*c1 + c3 + a*c2 - a*c3)   ;   (4*c1 - c3 + a*c3)*(c2 - 2*c1 + a*c2); ; c2*(3*a + 1)*(a - 1) + (c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4]
solve(Q)
Q.c1
S = solve(Q)
S.c1
a=4; Q1 = eval(Q)
subs(Q1,'c1', '(5*c2 - 3*c3)/2')
subs(Q1,'c2', '(3*c3)/10')
subs(Q1,'c1', '(5*c2 - 3*c3)/2','c2', '(3*c3)/10')
subs(subs(Q1,'c1', '(5*c2 - 3*c3)/2',)'c2', '(3*c3)/10')
subs(subs(Q1,'c1', '(5*c2 - 3*c3)/2'), 'c2', '(3*c3)/10')
Q
subs(subs(Q1, 'c2', '1')
subs(Q1, 'c2', '1')
W
simplify(W)
latex(simplify(W))
W30 = (c1*c2*(3*a + 1))/(2*c1 + c2 + a*c2)^2
diff(W30,c2)
c2=1;
eval(diff(W30,c2))
S=solve(eval(diff(W30,c2)))
c2=0.1;
S=solve(eval(diff(W30,c1)))
c2=1;
S=solve(eval(diff(W30,c1)))
c2=0.02;
S=solve(eval(diff(W30,c1)))
W30
[x y]=0:0.1:1
[x]=0:0.1:1
[y]=0:0.1:1
c1=0:0.1:1
mesh(c1,c2,V30)
mesh(c1,c2,W30)
[c1 c2]=meshgrid(0:0.1:1)
mesh(W30)
W30
mesh((13*c1.*c2)/(2*c1. + 5*c2.)^2)
mesh((13*c1.*c2)/(2*c1 + 5*c2)^2)
pcolor((13*c1.*c2)/(2*c1 + 5*c2)^2)
[c1 c2]=meshgrid(0:0.01:1)
pcolor((13*c1.*c2)/(2*c1 + 5*c2)^2)
mesh((13*c1.*c2)/(2*c1 + 5*c2)^2)
Q = [ (c2 - 2*c1 + c3 + a*c2 - a*c3)   ;   (4*c1 - c3 + a*c3)*(c2 - 2*c1 + a*c2); ; c2*(3*a + 1)*(a - 1) + (c3*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/4]
simplify(D_V_F)'
W
clear
ShujieLagrange
W = L2 - meanUtau
W30 = subs(W, c3, 0)
W30a4 = subs(subs(W, c3, 0),a,4)
W30a4 = subs( subs(subs(W, c3, 0),a,4), c1,1-c2)
W30a4 = subs( subs(subs(W, c3, 0),a,4), sqrt(c1,1-c2^2))
W30a4 = subs( subs(subs(W, c3, 0),a,4), c1,sqrt(1-c2^2))
W30a4 = subs( subs(subs(W, c3, 0),a,4), c2,sqrt(1-c1^2))
diff(W30a4)
c1=0:0.1:1; eval(diff(W30a4))
c1=0:0.01:1; z = eval(diff(W30a4)); plot(z)
W
fmincon
simplify(W)
-(c1*c2*(3*a + 1))/(2*c1 + c2 + a*c2)^2
W30 = -(c1*c2*(3*a + 1))/(2*c1 + c2 + a*c2)^2
ShujieLagrange
W30 = -(c1*c2*(3*a + 1))/(2*c1 + c2 + a*c2)^2
diff(W30,c1)
diff(W30,c1) * (2*c1 + c2 + a*c2)^3
simplify(diff(W30,c1) * (2*c1 + c2 + a*c2)^3)
diff(W30,c1)
S=solve((4*c1*c2*(3*a + 1))/(2*c1 + c2 + a*c2)^3 - (c2*(3*a + 1))/(2*c1 + c2 + a*c2)^2,c1)
S = sovle(diff(W30,c1),c1)
S = solve(diff(W30,c1),c1)
S = solve(diff(W30,c2),c2)
S = solve(diff(W30,c1),c1)
ShujieLagrange
eigen_vector(:,1) * S1 + c2*eigen_vector(:,2)
dist = eigen_vector(:,1) * S1 + c2*eigen_vector(:,2)
p = dist/sum(dist)
p = simplify(dist/sum(dist))
dist = eigen_vector(:,1)*c1 + S2*eigen_vector(:,2)
p = simplify(dist/sum(dist))
W31=(c3^2*abs(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/8 - (c2*(3*a + 1)*(4*c1 - c3 + a*c3))/(4*(2*c1 + c2 + a*c2)^2)
W31=(c3^2*(3*a + 1)*(3*a^2 + 2*a + 11)^(1/2))/8 - (c2*(3*a + 1)*(4*c1 - c3 + a*c3))/(4*(2*c1 + c2 + a*c2)^2)
simplify(W31)
S = solve(diff(W31,c1),c1)
solve(S^2+c2^2+c3^2=1,c3)
solve(S^2+c2^2+c3^2-1,c3)
C3 = solve(S^2+c2^2+c3^2-1,c3))
C3 = solve(S^2+c2^2+c3^2-1,c3)
simplify(C3)
S
C1 = subs(S,c3,C3)
W
W2 = subs(subs(W,c1,C1),c3,C3)
W
C1
S
S = solve(diff(W31,c1),c1)
W31
subs(W31,c1,sqrt(1-c2^2-c3^2))
subs(W31,c3,sqrt(1-c1^2-c2^2))
W312 = subs(W31,c3,sqrt(1-c1^2-c2^2))
C2=solve(W312,c2)
latex(W312)
W31
W312 = subs(W31,c3,sqrt(1-c1^2-c2^2))
subs(W312,a,4)
solve(subs(W312,a,4), c2)
XX = solve(subs(W312,a,4), c2)
XX = solve(subs(W312,a,4), c2,'ReturnConditions', true)
XX.c2
XX = solve(subs(W312,a,10), c2,'ReturnConditions', true)
XX.c2
x
X
XX = solve(subs(W312,a,10.01), c2,'ReturnConditions', true)
c2
W312
subs(W312,a,10.01)
vap
vpa
vpa(W312,a,10.01)
vpa(subs(W312,a,10.01))
solve(vpa(subs(W312,a,10.01)))
vpasolve(vpa(subs(W312,a,10.01)))
c1
vpasolve(vpa(subs(W312,a,10.01),c2))
vpasolve(vpa(subs(W312,a,10.01)),c2)
vpa(subs(W312,a,10.01))
vpasolve(vpa(subs(W312,a,10.01)),'c2')
vpasolve(vpa(subs(W312,a,10.01))=0)
vpasolve((subs(W312,a,4)))
subs(W312,a,4)
solve(subs(W312,a,4),c2)
XX = solve(subs(W312,a,4), c1,'ReturnConditions', true)
XX = solve(subs(W312,c1,0.24), c2,'ReturnConditions', true)
c2
XX = vpasolve(subs(subs(W312,c1,0.24),a,4) c2,'ReturnConditions', true)
XX = vpasolve(subs(subs(W312,c1,0.24),a,4), c2,'ReturnConditions', true)
XX = vpasolve(subs(subs(W312,c1,0.24),a,4), c2)
c1^2+c2^2
0.24^2+XX^2
XX = vpasolve(subs(subs(W312,c1,0.8),a,4), c2)
c1
0.8^2+XX^2
XX = vpasolve(subs(subs(W312,c1,0.9),a,4), c2)
XX = vpasolve(subs(subs(W312,c1,0.99),a,4), c2)
XX = vpasolve(subs(subs(W312,c1,0.99),a,0.2), c2)
XX = vpasolve(subs(subs(W312,c1,0.5),a,0.2), c2)
XX = vpasolve(subs(subs(W312,c1,0.5),a,4), c2)
XX = vpasolve(subs(subs(W312,c1,0.707),a,4), c2)
XX = vpasolve(subs(subs(W312,c1,0.707),a,1/4), c2)
ShujieLagrange
W
subs(W,c3,0)
solve(subs(W,c3,0),c1)
solve(subs(W,c3,0),c2)
solve(-((3*a + 1))/(2*c1 + c2 + a*c2)^2,c2)
solve(diff(subs(W,c3,0),c2),c2)
solve(diff(subs(W,c3,0.2),c2),c2)
vpasolve(diff(subs(W,c3,0.2),c2),c2)
vpasolve(diff(subs(W,c3,0.2),c1),c1)
vpasolve(diff(subs(W,c3,0),c1),c1)
vpasolve(diff(subs(W,c3,0.5),c1),c1)
vpasolve(diff(subs(W,c3,0.8),c1),c1)
vpasolve(diff(subs(W,c3,0.18),c1),c1)
vpasolve(diff(subs(W,c3,0.8),c1),c1)
vpasolve(diff(subs(W,c3,0.1),c1),c1)
W
solve(diff(subs(W),c1),c1)
solve(c2/2 + c3/2 + (a*c2)/2 - (a*c3)/2, c3)
solve(diff(subs(W),c2),c2)
solve(diff(subs(W),c3),c3)
simplify(ans)
W
solve(diff(subs(W),c2),c2)
solve(diff(subs(W),c1),c1)
S2 = solve(diff(W30,c2),c2)
%-- 2021/5/23 9:02 --%
mycon
clc
clear
close all
fun=@(x)x(1)^2+x(2)^2+8;
x0=rand(2,1);
A=[];
b=[];
Aeq=[];
beq=[];
vlb=[0,0];
vub=[];
exitflag=1;
[x,fval,exitflag]=fmincon(fun,x0,A,b,Aeq,beq,vlb,vub,'mycon')
Shujie_fmincon
4/13
a_value=1;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value)
clear
a_value=1;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value)
a_value=1;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);
c1c2c3
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,fval];end
plot(rr(:,1:2),'DisplayName','rr(:,1:2)')
plot(rr(:,end))
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,fval/(a/(3*a+1))];end
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,fval/(a_value/(3*a_value+1))];end
plot(rr(1:39,end))
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,fval, fval/(a_value/(3*a_value+1))];end
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
plot(rr(:,4))
plot(rr(:,end))
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
ShujieLagrange
Wplus
subs(Wplus,c3,'x(3)')
subs(subs(subs(Wplus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)')
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
subs(subs(subs(subs(Wplus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)',)a,'a_value')
subs(subs(subs(subs(Wplus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'),a,'a_value')
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
Shujie_fmincon
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
clear
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
[g,ceq]=mycon(1)
[g,ceq]=mycon(x)
syms x real
[g,ceq]=mycon(x)
syms x(1:3) real
syms x(1) x(2) x(3) real
[g,ceq]=mycon(1)
x=sym('x',[3,1])
[g,ceq]=mycon(x)
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
clear
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
%-- 2021/5/23 13:04 --%
ShujieLagrange
clear; rr=[]; for a_value=0.1:0.1:4;
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];
end
plot(rr(:,1:2),'DisplayName','rr(:,1:2)')
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
4/13
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
ShujieLagrange
eigen_vector(:,2)/norm(eigen_vector(:,2))
[eigen_vector(:,1)/norm(eigen_vector(:,1)) eigen_vector(:,2)/norm(eigen_vector(:,2)) eigen_vector(:,3)/norm(eigen_vector(:,3)) eigen_vector(:,4)/norm(eigen_vector(:,4))]
ShujieLagrange
meanUtau
subs(subs(subs(subs(meanUtau,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'),a,'a_value')
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
ShujieLagrange
meanUtau
subs(subs(subs(subs(meanUtau,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'),a,'a_value')
clear
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1))];end
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
ShujieLagrange
dXn13
simplify(limit(int(simplify(dXn13),'t',[0,tau])/tau,'tau',0));
ShujieLagrange
meanUtau
Xn
meanXntau = simplify(limit(int(simplify(Xn),'t',[0,tau])/tau,'tau',0));
meanXntau = (limit(int((Xn),'t',[0,tau])/tau,'tau',0));
clear
ShujieLagrange
meanUtau
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value)
a_value=4
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value)
subs(subs(subs(subs(meanUtau,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'),a,'a_value')
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
ShujieLagrange
tmp
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
sum(Xn)
simplifiy(sum(Xn))
simplify(sum(Xn))
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
tmpW
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
tmpW
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
ShujieLagrange
norm(eigencycle)
norm(eigencycle))
clear
ShujieLagrange
clear
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
plot(rr(:,3))
ShujieLagrange
kl=Xn./[a/(3*a+1) a/(3*a+1) a/(3*a+1) 1/(3*a+1)]
kl=Xn./[a/(3*a+1) a/(3*a+1) a/(3*a+1) 1/(3*a+1)]'
Entropy = sum(Xn*log(kl))
log(kl)
Entropy = sum(Xn.*log(kl))
ShujieLagrange
a=4
kl
subs(kl, a, 4)
subs(kl, 'a', 4)
Entropy_D = sum(Xn.*log(kl))
Entropy_D = subs(Entropy_D,'a',4)
meanUtau = simplify(limit(int(simplify(mean_U - Entropy_D),'t',[0,tau])/tau,'tau',0))
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
ShujieLagrange
tmpWplus
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
clear
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) exitflag];end
rr=[]; for a_value=0.1:0.2:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
rr=[]; for a_value=0.1:0.3:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
ShujieLagrange
rr=[]; for a_value=0.1:0.3:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
ShujieLagrange
rr=[]; for a_value=0.1:0.3:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
ShujieLagrange
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
ShujieLagrange
%-- 2021/5/23 20:10 --%
rr=[]; for a_value=0.1:0.1:4;[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(a_value);rr=[rr; c1c2c3,-fval, -fval/(a_value/(3*a_value+1)) a_value exitflag];end
ShujieLagrange
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(tmp_W_minus)
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(@tmp_W_minus)
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(tmp_W_minus)
tmp_W_minus
sym3str( tmp_W_minus)
symstr( tmp_W_minus)
syms( tmp_W_minus)
g = matlabFunction(tmp_W_minus )
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(r)
g = matlabFunction(tmp_W_minus )
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(g)
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(g)
g = matlabFunction(tmp_W_minus )
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(g)
g = matlabFunction(tmp_W_minus,'Vars',{x} )
g = matlabFunction(tmp_W_minus,'Vars',{[x]} )
g = matlabFunction(tmp_W_minus,'Vars','x' )
ShujieLagrange
[c1c2c3,fval,exitflag]=Shujie_fmincon_in_a_out_c1c2c3(g)
fun
fun2
ShujieLagrange
x0= [0.707 0.707 0]'+0.1*rand(3,1)
A=[];
b=[];
Aeq=[];
beq=[];
% vlb=[0,0];
vlb=[-2 -2 -2];
vub=[2 2 2];
[x_result,fval,exitflag]=fmincon(g,x0,A,b,Aeq,beq,vlb,vub,'mycon') ;
ShujieLagrange
clear
ShujieLagrange
rr=[];
for a=0.1:0.1:2
% tmp_W_minus = subs(subs(subs(subs(W_minus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'),a,'a_value')
tmp_W_minus = subs(subs(subs(subs(W_minus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'))
g = matlabFunction(tmp_W_minus,'Vars','x' )
x0= [0.5 0.5 0]'+ 0.1*rand(3,1)
A=[];
b=[];
Aeq=[];
beq=[];
% vlb=[0,0];
vlb=[-2 -2 -2];
vub=[2 2 2];
[x_result,fval,exitflag]=fmincon(g,x0,A,b,Aeq,beq,vlb,vub,'mycon') ;
c1c2c3=x_result';
rr=[rr; c1c2c3,-fval, (a/(3*a+1)) -fval/(a/(3*a+1)) a exitflag];
end
rr=[];
for a=0.0:0.25:4
% tmp_W_minus = subs(subs(subs(subs(W_minus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'),a,'a_value')
tmp_W_minus = subs(subs(subs(subs(W_minus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'))
g = matlabFunction(tmp_W_minus,'Vars','x' )
x0= [0.5 0.5 0]'+ 0.1*rand(3,1)
A=[];
b=[];
Aeq=[];
beq=[];
% vlb=[0,0];
vlb=[-2 -2 -2];
vub=[2 2 2];
[x_result,fval,exitflag]=fmincon(g,x0,A,b,Aeq,beq,vlb,vub,'mycon') ;
c1c2c3=x_result';
rr=[rr; c1c2c3,-fval, (a/(3*a+1)) -fval/(a/(3*a+1)) a exitflag];
end
ShujieLagrange
clear
ShujieLagrange
clear
ShujieLagrange
clear
ShujieLagrange
clear
ShujieLagrange
clear
clear
ShujieLagrange
clear
ShujieLagrange
rr
ShujieLagrange
%-- 2021/5/24 0:45 --%
ShujieLagrange
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);eval(Xn)
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);vpa(Xn)
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);eval(Xn)
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);vpa(eval(Xn))
format short
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);vpa(eval(Xn))
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);simplify(vpa(eval(Xn)))
limit( (simplify(real(simplify(vpa(eval(Xn))))) ) ,'t',0)
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);p=round(limit( (simplify(real(simplify(vpa(eval(Xn))))) ) ,'t',0)',4)
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);limit( (simplify(real(simplify(vpa(eval(Xn))))) ) ,'t',0)'
round(ans,4)
roundn(ans,4)
roundn(double(ans),4)
double(ans)
ans
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);limit( (simplify(real(simplify(vpa(eval(Xn))))) ) ,'t',0)'
double(ans)
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);double(limit( (simplify(real(simplify(vpa(eval(Xn))))) ) ,'t',0)')
ShujieLagrange
clear
ShujieLagrange
plot(rr(:,9:end),'DisplayName','rr(:,9:end)')
ShujieLagrange
limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0)
xxxx = limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0)
eval(xxxx)
xxxx = eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
ShujieLagrange
plot(rr(:,13:end),'DisplayName','rr(:,13:end)')
ShujieLagrange
clear
format long
ShujieLagrange
clear
syms c1 c2 c3 c4 a t tau real
eigen_vector = [
[ 1, 0, a*(1/4 + 3i/4) - 1/4 + 1i/4, a*(1/4 - 3i/4) - 1/4 - 1i/4]
[ 0, a,                 - a/2 - 1/2,                 - a/2 - 1/2]
[ 1, 0, a*(1/4 - 3i/4) - 1/4 - 1i/4, a*(1/4 + 3i/4) - 1/4 + 1i/4]
[ 0, 1,                           1,                           1]]
% eigen_vector = [eigen_vector(:,1)/norm(eigen_vector(:,1)) eigen_vector(:,2)/norm(eigen_vector(:,2)) eigen_vector(:,3)/norm(eigen_vector(:,3)) eigen_vector(:,4)/norm(eigen_vector(:,4))]
eigen_value = [
[ -1,  0,   0,  0]
[  0, -1,   0,  0]
[  0,  0, -1i,  0]
[  0,  0,   0, 1i]] %* a/(3*a+1)
X0 = c1*eigen_vector(:,1) * exp(eigen_value(1,1)*t) + ...
c2*eigen_vector(:,2) * exp(eigen_value(2,2)*t) + ...
c3*eigen_vector(:,3) * exp(eigen_value(3,3)*t)
X1 =   real(X0)
Xn =  X1/sum(X1)
payoff_matrix = [0 0 0 a;
1 0 0 0;
0 1 0 0;
0 0 1 0];
eigen_vector =
[ 1, 0, a*(1/4 + 3i/4) - 1/4 + 1i/4, a*(1/4 - 3i/4) - 1/4 - 1i/4]
[ 0, a,                 - a/2 - 1/2,                 - a/2 - 1/2]
[ 1, 0, a*(1/4 - 3i/4) - 1/4 - 1i/4, a*(1/4 + 3i/4) - 1/4 + 1i/4]
[ 0, 1,                           1,                           1]
eigen_value =
-1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i  -1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 - 1.0000i   0.0000 + 0.0000i
0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 1.0000i
X0 =
c1*exp(-t) + c3*exp(-t*1i)*(a*(1/4 + 3i/4) - 1/4 + 1i/4)
a*c2*exp(-t) - c3*exp(-t*1i)*(a/2 + 1/2)
c1*exp(-t) + c3*exp(-t*1i)*(a*(1/4 - 3i/4) - 1/4 - 1i/4)
c2*exp(-t) + c3*exp(-t*1i)
X1 =
c1*exp(-t) + c3*cos(t)*(a/4 - 1/4) + c3*sin(t)*((3*a)/4 + 1/4)
a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2)
c1*exp(-t) + c3*cos(t)*(a/4 - 1/4) - c3*sin(t)*((3*a)/4 + 1/4)
c3*cos(t) + c2*exp(-t)
Xn =
(c1*exp(-t) + c3*cos(t)*(a/4 - 1/4) + c3*sin(t)*((3*a)/4 + 1/4))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
(a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
(c1*exp(-t) + c3*cos(t)*(a/4 - 1/4) - c3*sin(t)*((3*a)/4 + 1/4))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
(c3*cos(t) + c2*exp(-t))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
Payoff_vector_field_F = payoff_matrix * Xn
Payoff_vector_field_F =
(a*(c3*cos(t) + c2*exp(-t)))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
(c1*exp(-t) + c3*cos(t)*(a/4 - 1/4) + c3*sin(t)*((3*a)/4 + 1/4))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
(a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
(c1*exp(-t) + c3*cos(t)*(a/4 - 1/4) - c3*sin(t)*((3*a)/4 + 1/4))/(c3*cos(t) + 2*c1*exp(-t) + c2*exp(-t) + a*c2*exp(-t) - c3*cos(t)*(a/2 + 1/2) + 2*c3*cos(t)*(a/4 - 1/4))
mean_U = Xn' * Payoff_vector_field_F
meanUtau = simplify(limit(int(simplify(real(mean_U)),'t',[0,tau])/tau,'tau',0))
tmp_meanUtau = subs(subs(subs(subs(meanUtau,c1,'x(1)'),c2,'x(2)'), ...
c3,'x(3)'),a,'a_value')
eigencycle = [];
for m=1:4-1
for n=m+1:4
L0= cross([real(eigen_vector(m,3)) imag(eigen_vector(m,3)) 0 ], ...
[real(eigen_vector(n,3)) imag(eigen_vector(n,3)) 0 ])
eigencycle = [eigencycle; L0(3)];
end
end  %eigen_cycle = pi*real(cross([real(m) imag(m) 0],[real(n) imag(n) 0]));
L_motion = eigencycle * c3^2 * imag(eigen_value(3,3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
L1 = sqrt(L_motion'*L_motion)
L2 = abs(L1*imag(eigen_value(3,3)))/2
W_minus = L2 - meanUtau
rr=[];
for a=0.0:0.25:4
tmp_W_minus = subs(subs(subs(subs(W_minus,c1,'x(1)'),c2,'x(2)'),c3,'x(3)'))
g = matlabFunction(tmp_W_minus,'Vars','x' )
x0= [0.5 0.5 0]'+ 0.1*rand(3,1)
A=[];
b=[];
Aeq=[];
beq=[];
% vlb=[0,0];
vlb=[-2 -2 -2];
vub=[2 2 2];
[x_result,fval,exitflag]=fmincon(g,x0,A,b,Aeq,beq,vlb,vub,'mycon') ;
c1c2c3=x_result';
rr=[rr; c1c2c3, -fval, (a/(3*a+1)) -fval/(a/(3*a+1)) a exitflag];end
ShujieLagrange
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0.1)')
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0.1)')
p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',1/(2*pi)))'
p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',(2*pi)))'
p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',.3))'
plot(rr(:,2:4),'DisplayName','rr(:,2:4)')
sum(Xn)
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);t=0;
eval(sum(Xn))
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);t=0.3;
eval(sum(Xn))
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);t=3.3;
eval(sum(Xn))
eval((Xn))
c1=c1c2c3(1);c2=c1c2c3(2);c3=c1c2c3(3);t=.3;
eval((Xn))
%-- 2021/5/24 10:58 --%
ShujieLagrange
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0.1)')
p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',1/(2*pi)))'
pp=[pp; rr(k,2),   p1234 p4tau];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0.1)')
p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',1/(2*pi)))'
pp=[pp; rr(k,1),   p1234 p4tau];
end
ShujieLagrange
pp=[]; Ne=[a a a 1]/(3*a+1);
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4); a = rr(k,1);
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp; a p1234  Ne];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4); a = ;
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) ];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) ];
end
clear
ShujieLagrange
x0=   rand(3,1).*[1 1 -1]'
ShujieLagrange
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + c3*real(eigen_vector_n(:,3))
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + c3*eval(real(subs(eigen_vector_n(:,3),a,rr(k,1)))
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + ...
c3*eval(real(subs(eigen_vector_n(:,3),a,rr(k,1))))
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + ...
c3*eval(real(subs(eigen_vector_n(:,3),a,0.00001+rr(k,1))))
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);  a=rr(k,1);
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + ...
c3*eval(real(eigen_vector_n(:,3)))
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
eigen_vector_n(:,3)
vpa(eigen_vector_n(:,3))
eval(eigen_vector_n(:,3))
ShujieLagrange
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + ...
c3*eval(real(eigen_vector_n(:,3)))
eval( (eigen_vector_n(:,3)))
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);  a=rr(k,1);
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + ...
c3*real(eval(eigen_vector_n(:,3)))
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
eval( (eigen_vector_n(:,3)))
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);  a=rr(k,1);
tmp=eval(eigen_vector_n(:,3))
p1234_C =  c1*eigen_vector_n(:,1)  + c2*eigen_vector_n(:,2) + ...
c3*real(tmp)
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);  a=rr(k,1);
tmp=eval(eigen_vector_n(:,3))
p1234_C =  c1*eigen_vector_n(:,1)  + eval(c2*eigen_vector_n(:,2)) + ...
c3*real(tmp)
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
pp=[];
for k=1:size(rr,1)
c1=rr(k,2);c2=rr(k,3);c3=rr(k,4);  a=rr(k,1);
tmp=eval(eigen_vector_n(:,3))
p1234_C =  eval(c1*eigen_vector_n(:,1))  + eval(c2*eigen_vector_n(:,2)) + ...
c3*real(tmp)
p1234 = double(limit(real(simplify(vpa(eval(Xn)))) ,'t',0)')
%       p4tau= eval(limit(int(simplify(real(Xn)),'t',[0,tau])/tau,'tau',0))'
pp=[pp;   rr(k,1) p1234  [rr(k,1) rr(k,1) rr(k,1) 1]/(3*rr(k,1)+1) p1234_C'];
end
ShujieLagrange
Xn
t=0
eval(Xn)
%-- 2021/5/24 14:19 --%
ShujieLagrange
[V D]=eig(rand(5))
sum(V)
Strategy5
payoff_matrix = -0.001 * [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0];
Strategy5
Strategy5([0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0])
eval(D_V_F)
sum(eigen_vector)
maxVnash = real(eigen_vector(:,3))'/sum(real(eigen_vector(:,3))) - Nash
x1=A(length(A(:,1)),1)-0.02; x2=A(length(A(:,1)),2)-0.02;
x3=A(length(A(:,1)),3); x4=A(length(A(:,1)),4); x5=A(length(A(:,1)),5)+0.02;
x1=A(length(A(:,1)),1)-0.02; x2=A(length(A(:,1)),2)-0.02;
x3=A(length(A(:,1)),3); x4=A(length(A(:,1)),4); x5=A(length(A(:,1)),5)+0.04;
[eigen_vector eigen_value] = eig(eval(D_V_F))
sum(eigen_vector)
x1=A(length(A(:,1)),1)-0.02; x2=A(length(A(:,1)),2)-0.06;
x3=A(length(A(:,1)),3)+0.04; x4=A(length(A(:,1)),4); x5=A(length(A(:,1)),5)+0.04;
[eigen_vector eigen_value] = eig(eval(D_V_F))
sum(eigen_vector)
eval(D_V_F)
sum(D_V_F)
sum(eval(D_V_F))
sum(eval(D_V_F)')
clear
Strategy5([0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0])
eval(D_V_F)
D_V_F
sum(D_V_F)
simplify(sum(D_V_F))
a=  simplify(sum(D_V_F))
a(1)-a(2)
simplify(a(1)-a(2))
solve(a(1)-a(2),a(1)-a(3),a(1)-a(4),a(1)-a(5), 1-x1-x2-x3-x4 )
S = solve(a(1)-a(2),a(1)-a(3),a(1)-a(4),a(1)-a(5), 1-x1-x2-x3-x4 )
S.x1
simplify(sum(D_V_F))
latex(D_V_f)
latex(D_V_F)
latex(D_V_F')
latex(simplify(D_V_F)')
a=rand(5)
a-sum(a)*eyes(5)
b=a-sum(a)*eye(5)
sum(a)
b=a-eye(5)*sum(a)
b=a-eye(5).*sum(a)
b=a-eye(sum(a))
b=a.-sum(a)
b=a(1,:) - sum(a)
a(1,:)=a(1,:) - sum(a)
[V D]=eig(a)
sum(V)
eval(D_V_F)
mean(mean( eval(D_V_F) ))
eval(D_V_F) - mean(mean( eval(D_V_F) ))
sum(diag(eval(D_V_F) - mean(mean( eval(D_V_F) ))))
X1 = eval(D_V_F) - mean(mean( eval(D_V_F) ))
diag(X1)/-6.1540
[V D] = eig(X1)
sum(sum(X1))
sum(V)
X0 = eval(D_V_F)
sum(X0)
Oneill_replicator
sum(D_Eq_at_NE)
ShujieLandscape1228
X0 = eval(D_V_F)
sum(X0)
a=sum(X0)
a(1)-a(2)
a(1)-a(3)
simplify(a(1)-a(3))
simplify(a(1)-a(4))
X0 = eval(D_V_F)
simplify(D_V_F)
[eigen_vector eigen_value w] = eig(eval(D_V_F))
ShujieLandscape1228
%% Left Eigenvectors
%%
% Create a 3-by-3 matrix.
A = [1 7 3; 2 9 12; 5 22 7];
%%
% Calculate the right eigenvectors, |V|, the eigenvalues, |D|, and the left
% eigenvectors, |W|.
[V,D,W] = eig(A)
%%
% Verify that the results satisfy |W'*A = D*W'|.
W'*A - D*W'
ShujieLandscape1228
eval(D_V_F)
sum( eval(D_V_F) )
sum( eval(D_V_F)' )
sum( eval(D_V_F))
sum( eval(D_V_F))'
simplify(sum( eval(D_V_F))')
%-- 2021/5/24 21:30 --%
ShujieLandscape1228
eval(D_V_F)
sum( eval(D_V_F) )
sum( eval(D_V_F)' )
simplify(sum( eval(D_V_F))')
simplify(sum( eval(D_V_F)))
ShujieLandscape1228
DD_V_F = [diff(D_V_F(1,:)','x1') diff(D_V_F(2,:)','x2') diff(D_V_F(3,:)','x3') diff(D_V_F(4,:)','x4')]
sum(DD_V_F)
sum(DD_V_F)'
mean_U
v_F = [mean_U,'x1'  diff(mean_U,'x2'  diff(mean_U,'x3'  diff(mean_U,'x4' ]'
v_F = [diff(mean_U,'x1')  diff(mean_U,'x2')  diff(mean_U,'x3')  diff(mean_U,'x4') ]'
V_F
d_v_F = [diff(v_F,'x1')  diff(v_F,'x2')  diff(v_F,'x3')  diff(v_F,'x4') ]
[v d]=eig(d_v_F)
sum(v)
ShujieLandscape1228
x1=a/(3*a + 1);x2=a/(3*a + 1);;x3=a/(3*a + 1);1/4;x4=1/(3*a + 1);
eval(D_V_F)
A = eval(D_V_F)
T = A-A'
simplify(T)
latex(simplify(T))
sum(T)
simplify(sum(T))
simplify(sum(T*(3*a + 1)^2))
A = eval(D_V_F)
sum(A)
simplify(sum(A))
simplify(sum(A))*eye(4)
diag(A)
A
eigenvector=eig(A)
[eigenvector d]=eig(A)
ShujieLagrange
eigencycle
latex(eigencycle/((3*a)/4 + 1/4))
ShujieLandscape1228
(eval(D_V_F))
eval(DD_V_F)
x1=a/(3*a + 1);x2=a/(3*a + 1);;x3=a/(3*a + 1);1/4;x4=1/(3*a + 1);
Ne = [x1 x2 x3 x4];
eval(DD_V_F)
eval(D_V_F)
eval(DD_V_F) - eval(DD_V_F)'
simplify(eval(DD_V_F) - eval(DD_V_F)')
simplify(eval(DD_V_F) - eval(DD_V_F)')*(3*a+1)
simplify(eval(D_V_F) - eval(D_V_F)')*(3*a+1)
simplify(eval(D_V_F) - eval(D_V_F)')*(3*a+1)^2
simplify(eval(D_V_F) - eval(D_V_F)')*(3*a+1)^2/a
sum(eval(D_V_F))
simplify(sum(eval(D_V_F)))
ShujieLandscape1228
Jacobian = eval(D_V_F)
[eigen_vector eigen_value w] = eig(Jacobian)
simplify(Jacobian)
Jacobian = simplify(eval(D_V_F))
[eigen_vector eigen_value w] = eig(Jacobian)
sum(Jacobian)
simplify(sum(Jacobian))
sum(Jacobian')
simplify(sum(Jacobian'))
simplify(sum(Jacobian'))'
sum(simplify(sum(Jacobian'))')
simplify(sum(simplify(sum(Jacobian'))'))
a=1;eval(simplify(sum(Jacobian'))')
a=4;eval(simplify(sum(Jacobian'))')
a=1/4;eval(simplify(sum(Jacobian'))')
Jacobian
Jacobian'
[v1 d1]=eig(Jacobian')
[v1 d1 w]=eig(Jacobian')
[V,D,W] = eig(Jacobian')
[V,D,W] = eig(Jacobian)
Jacobian
Jacobian/(-(3*a + 1))
Jacobian*(-(3*a + 1))
syms a
Jacobian*(-(3*a + 1))
ShujieLandscape1228
simplify(sum(Jacobian))
Jacobian'.*Ne'
Ne'
Jt = Jacobian';
for i=1:4;Jt(i,:)= Jt(i,:)*Ne(i);end
Jt
Asym = Jt-Jt'
Asym = simplify(Jt-Jt')
Asym = simplify(Jt-Jt')*(3*a + 1)^2
sum([     -(2*a^2)/(3*a + 1)^2,  (a*(a + 1))/(3*a + 1)^2,     -(2*a^2)/(3*a + 1)^2,   -(2*a)/(3*a + 1)^2])
sim[lify(sum([     -(2*a^2)/(3*a + 1)^2,  (a*(a + 1))/(3*a + 1)^2,     -(2*a^2)/(3*a + 1)^2,   -(2*a)/(3*a + 1)^2]))
simplify(sum([     -(2*a^2)/(3*a + 1)^2,  (a*(a + 1))/(3*a + 1)^2,     -(2*a^2)/(3*a + 1)^2,   -(2*a)/(3*a + 1)^2]))
sum(Jt)
simplify(sum(Jt))
simplify(sum(Jt'))
Asym = simplify(Jt-Jt')*(3*a + 1)^2/a
%-- 2021/5/25 23:04 --%
NoDyamics
v_F = [diff(D_V ,'x1'); diff(D_V ,'x2'); diff(D_V ,'x3'); diff(D_V ,'x4')]
D_F = [diff(V_F,'x1') diff(V_F,'x2') diff(V_F,'x3') diff(V_F,'x4')]
V_F = [diff(D_V ,'x1'); diff(D_V ,'x2'); diff(D_V ,'x3'); diff(D_V ,'x4')]
D_F = [diff(V_F,'x1') diff(V_F,'x2') diff(V_F,'x3') diff(V_F,'x4')]
[v d] = eig(D_F)
NoDyamics
diff(diff(D_V ,'x1'),'x3')
diff(diff(D_V ,'x1'),'x2')
[v d] = eig(D_F)
simplify(s)
simplify(d)
mean_U.*[x1 x2 x3 x4]'
NoDyamics
simplify(V_F)
[v d] = eig(DV_F)
[v d] = eig(V_F)
[v d] = eig(V_F')
NoDyamics
V_F = [diff(D_V ,'x1'); diff(D_V ,'x2'); diff(D_V ,'x3'); diff(D_V ,'x4')].*[x1 x2 x3 x4]'
D_F = [diff(V_F,'x1') diff(V_F,'x2') diff(V_F,'x3') diff(V_F,'x4')]
[v d] = eig(D_F)
simplify(d)
x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);
eval(d)
simplify(eval(d))
subs(d,a,4)
x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);
simplify(eval(d))
vpa(subs(d,a,4))
v
v(:,1)
v(:,2)
v(:,3)
v(:,4)
x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);a=4;eval(v)
a=4;x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);simplify(eval(v))
a=4;x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);(eval(v))
a=0.25;x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);(eval(v))
a=0.25;x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);(eval(d))
NoDyamics
x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);a=4;eval( D_V)
a=4;x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);eval( D_V)
[v d] = eig( eval( D_V) )
NoDyamics
a=4;x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);eval( D_F)
[v d] = eig( eval( D_F) )
[v d] = eig( eval( D_F') )
NoDyamics
a=4;x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1);x4=1/(3*a + 1);eval( D_F)
[v d] = eig( eval( D_F') )
ShujieLandscape1228
simplify(eigen_value)
simplify(simplify(eigen_value))
a=4;eval(d)
a=4;eval(eigen_value)
ShujieLandscape1228
V_F = [x1 x2 x3 x4]'.*(Payoff_vector_field_F.^2 - mean_U^2)  %/mean_U;
ShujieLandscape1228
simplify(eigen_vector)
simplify(eigen_value)
ShujieLandscape1228
NoDyamics
simplify(eigen_vector)
NoDyamics
J2 = simplify(eval(DD_V_F))
x1=a/(3*a + 1); x2=a/(3*a + 1); x3=a/(3*a + 1); x4=1/(3*a + 1);
J2 = simplify(eval(DD_V_F))
k=1
J2 = simplify(eval(DD_V_F))
J2 = simplify(eval(DD_V_F))
[v d] = eig( eval( J2) )
a=4
J2 = simplify(eval(DD_V_F))
subs(J2,a,4)
subs(J2,'a',4)
[v d] = eig( eval( subs(J2,'a',4)) )
NoDyamics
eval(D_V_F)
simplify(eval(D_V_F))
eigen_valueeigen_value
eigen_value
eigen_vector
NoDyamics
x1
NoDyamics
%-- 2021/5/26 9:08 --%
YQM20200830
eval(D_Eq_0)
sum(eval(D_Eq_0))
NoDyamics
w'*Jacobian - eigen_value*w'
w*Jacobian - eigen_value*w
%-- 2021/5/26 13:46 --%
WangYao2105_hypofine
[h p]=ttest([0.009		0.01		0.018		0.019		0.008		0.02
0.021		0.002		0.02		0.006		0.015		0.014
0.01		0.003		0.021		0.005		0.015		0.018
0.011		0.008		0.015		0.011		0.014		0.013
])
[h p]=ttest([0.009		0.01		0.018		0.019		0.008		0.02  ...
0.021		0.002		0.02		0.006		0.015		0.014 ...
0.01		0.003		0.021		0.005		0.015		0.018 ...
0.011		0.008		0.015		0.011		0.014		0.013 ...
])
%-- 2021/5/26 20:49 --%
Oneill_replicator
simplifify(D_Eq_0)
simplify(D_Eq_0)
simplify(V_Eq_0)
simplify(simplify(V_Eq_0))
expand(V_Eq_0)
simplify(V_Eq_0)
factor(D_Eq_at_NE)
national(D_Eq_at_NE)
rat(D_Eq_at_NE)
D_Eq_at_NE
format rat
(D_Eq_at_NE)
latex(D_Eq_at_NE)
D_Eq_at_NE*25
a=D_Eq_at_NE*25
latex(a)
latex(num2str(a))
num2str(a)
matrix2latex(a)
latex2MxWithMxPrecision(A, precision)
latex2MxWithMxPrecision(a, 0)
%-- 2021/5/27 10:36 --%
ShujieLandscape1228
Jacobian
Jacobian*(3*a+1)
latex(Jacobian*(3*a+1))
ShujieLandscape1228
latex((3*a + 1)^2)
v_4strategy_2x2_logit
sum(D_Eq_at_NE)
%-- 2021/5/27 15:50 --%
syms a t b;int(exp(a*t),'t',[0 b])
instanL_2_components
simplify(L)
int(L,'t',[0, T])/T
syms T;int(L,'t',[0, T])/T
syms T;lim(L,t,0)
syms T;limit(L,t,0)
syms T;int(L,'t',[0, t])/t
%-- 2021/5/27 23:52 --%
syms a b t T real; R = int(exp(a*t)*sin(b*t),'t',[0, T])/T
R/exp(T*a)
R
limit(R,'T', inf)
limit(R,'T', Inf)
limit(R,'T', 0)
limit(R,'T', -Inf)
%-- 2021/5/28 16:08 --%
v_4strategy_2x2_logit
S=vpasolve(V_Eq_0)
A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4])
v_4strategy_2x2_logit
S=vpasolve(V_Eq_0)
A=eval([S.x1 S.x2 S.x3 S.x4 S.y1 S.y2 S.y3 S.y4])
shujie4logit
S=solve(V_F);
S=vpasolve(V_Eq_0)
S=vpasolve(V_Eq_1)
eval(S)
[S.x1 S.x2 S.x3 S.x4]
[S.lambda S.x1 S.x2 S.x3 S.x4]
[S.lamda S.x1 S.x2 S.x3 S.x4]
S=solve(V_F,x1,x2,x3,x4);
%-- 2021/5/29 20:44 --%
syms a x real
diff(cos(a*x),'x')
syms a b x real
int(cos(a*x)*sin(b*x),'x')
latex(int(cos(a*x)*sin(b*x),'x'))
syms a b c d x real
int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x')
R= int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x')
latex(R)
R= int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x',[0,T])
R= int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x',[0,x])
latex(R)
lim(R/x,'x',Inf)
limit(R/x,'x',Inf)
Q = limit(R/x,'x',Inf)
latex(Q)
limit(R/x,'x',0)
latex(limit(R/x,'x',0))
latex(simplify(limit(R/x,'x',0)))
%-- 2021/5/29 23:41 --%
syms a b c d x real
R= int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x')
R0x = int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x',[0,x])
Q = limit(R0x/x,'x',Inf)
P = limit(R0x/x,'x',0)
clear
syms a b c d x real
R= int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x')
R0x = int(-cos(a*x+c)*sin(b*x+d)*b + sin(a*x+c)*cos(b*x+d)*a,'x',[0,x])
Q = limit(R0x/x,'x',Inf)
P = limit(R0x/x,'x',0)
latex([R0x; Q; P])
latex(simplify(P)
latex(simplify(P))
latex([R0x; Q; simplify(P)])
syms w_1 w_2 phi_1 phi_2 r_1 r_2 t real
R0x = int(-cos(w_1*t+phi_1)*sin(w_2*t+phi_2)*w_2 + sin(w_1*t+phi_1)*cos(w_2*t+phi_2)*w_1,'t',[0 t])
R0x = int( exp(r_1*t+r_2*t)*(-cos(w_1*t+phi_1)*sin(w_2*t+phi_2)*w_2 + sin(w_1*t+phi_1)*cos(w_2*t+phi_2)*w_1),'t',[0 t])
P = limit(R0x/x,'x',0)
P = limit(R0x/t,'t',0)
real(P)
Q = limit(R0x/t,'t',Inf)
R0x = int( exp(r_1*t)*(-cos(w_1*t+phi_1)*sin(w_2*t+phi_2)*w_2 + sin(w_1*t+phi_1)*cos(w_2*t+phi_2)*w_1),'t',[0 t])
P = limit(R0x/t,'t',0)
clear
R0x = int( exp(r_1*t)*(-cos(w_1*t+phi_1)*sin(w_2*t+phi_2)*w_2 + sin(w_1*t+phi_1)*cos(w_2*t+phi_2)*w_1),'t',[0 t])
syms w_1 w_2 phi_1 phi_2 r_1 r_2 t real
clear
syms w_1 w_2 phi_1 phi_2 r_1 r_2 t real
R0x = int( exp(r_1*t)*(-cos(w_1*t+phi_1)*sin(w_2*t+phi_2)*w_2 + sin(w_1*t+phi_1)*cos(w_2*t+phi_2)*w_1),'t',[0 t])
%-- 2021/5/31 14:05 --%
A=[12	0.003313	0.001703	0.00355	0.005107	0.001595	0.000945
13	-0.00154	-0.00089	-0.0031	-0.00354	-0.00097	-0.0012
14	-0.00139	-0.00106	-0.0013	-0.00277	-0.0015	-0.00092
15	-0.00038	0.00025	0.000846	0.001202	0.000882	0.001174
23	0.003499	0.003364	0.004895	0.00684	0.001811	0.002117
24	-0.0009	-0.0003	-0.00018	-0.00136	0.000597	-0.00128
25	0.000719	-0.00136	-0.00116	-0.00037	-0.00081	0.000104
34	0.002043	0.00114	0.002469	0.004032	0.000983	0.002386
35	-8.8E-05	0.001335	-0.00067	-0.00073	-0.00015	-0.00147
45	-0.00025	-0.00022	0.000985	-0.0001	7.66E-05	0.00019
]
[A(:,1) A(:,2:end)*1000]
[A(:,1) A(:,2:end)*100]
[A(:,1) A(:,2:end)*1000]
latex2MxWithMxPrecision(A, 3)
latex2MxWithMxPrecision([A(:,1) A(:,2:end)*1000], 3)
B = [L_E6	系数	系数t值	常数项	常数项t值
Replicator	0.0036 ***	6.48	0.00019	0.77
msReplicator	0.0036 ***	6.48	0.00019	0.77
logit2000	0.0033***	5.23	0.00026	0.87
logit00625	0.0035***	6.5	0.0002	0.79
logit002	0.0035 ***	6.49	0.0002	0.78
]
B = [
Replicator	0.0036 ***	6.48	0.00019	0.77
msReplicator	0.0036 ***	6.48	0.00019	0.77
logit2000	0.0033***	5.23	0.00026	0.87
logit00625	0.0035***	6.5	0.0002	0.79
logit002	0.0035 ***	6.49	0.0002	0.78
]
payoff_matrix =  [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0];[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
Strategy5([0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0])
payoff_matrix =  [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0];[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
r = [A(length(A(:,1)),:); diag(eigen_value)]
payoff_matrix =  [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0];[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
A(length(A(:,1)),:)
payoff_matrix =  [0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0];[eigen_vector, eigen_value, V_F, D_V_F, Nash, maxVnash]=Strategy5(payoff_matrix)
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);
r = [Nash; diag(eigen_value)']
latex(r)
latex2MxWithMxPrecision(r, 3)
latex2MxWithMxPrecision([r; eigen_vector], 3)
diag(eigen_value)]
diag(eigen_value)
diag(eigen_value)'
eva = diag(eigen_value)'
eve = eigen_vector
eve = round(eigen_vector,3)
eve = round(eigen_vector,4)
b=[0.1486	0.161083
0.2136	0.2245
0.2146	0.210556
0.0964	0.115583
0.3268	0.288278
]
scatter(b(:,1),b(:,2))
bar(b,'DisplayName','b')
%-- 2021/6/1 19:20 --%
a=load('F:\cDownload\PokerAnan20200716\8x8abedata\ZSJ20201215\6A25.csv');
[ sum_rr,inst_rr] = from_N_colExp_out_linear_regress(a)
%-- 2021/6/2 23:48 --%
read_abed_z19data
d3=num(15:end,2:4:14);
clear
read_abed_z19data
abedfilesname=strcat('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\data\1pop4-100-', ...
num2str(a(j_a)),'-',num2str(k_session),'.csv');
[num] = xlsread(abedfilesname);
d3=num(15:end,2:4:14);
save(num2str(j_a),d3);
save(strcat(num2str(j_a),'.mat'),d3);
strcat(num2str(j_a),'.mat')
save(strcat(num2str(j_a),'.mat'),'d3');
a=[0.25 1 4]
for j_a = 1:3
rr=[];
for k_session=1:5
abedfilesname=strcat('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\data\1pop4-100-', ...
num2str(a(j_a)),'-',num2str(k_session),'.csv');
[num] = xlsread(abedfilesname);
d3=num(15:end,[2:4:14 1]);
rr=[rr; d3];
end
save(strcat(num2str(j_a),'.mat'),'d3');
end
load('3.mat')
clear
load('3.mat')
read_abed_z19data
a=unique(rr,[1:4])
[C,ia,ic] = unique(rr,'rows');
[C,ia,ic] = unique(rr(:,1:4),'rows');
read_abed_z19data
clear
read_abed_z19data
for j_a = 1:3
rr=[];
for k_session=1:5
abedfilesname=strcat(pathd, 'data\1pop4-100-', ...
num2str(a(j_a)),'-',num2str(k_session),'.csv');
[num] = xlsread(abedfilesname);
d3=num(15:end,[2:4:14 1]);
d4=[d3(:,1)-d3(:,2) d3(:,2)-d3(:,3) d3(:,3)-d3(:,4) d3(:,4)];
rr=[rr; d3];
end
save(strcat(pathd, 'dataMatlab\',num2str(j_a),'.mat'),'rr');
end
for j_a = 1:3
rr=[];
for k_session=1:5
abedfilesname=strcat(pathd, 'data\1pop4-100-', ...
num2str(a(j_a)),'-',num2str(k_session),'.csv');
[num] = xlsread(abedfilesname);
d3=num(15:end,[2:4:14 1]);
d4=[d3(:,1)-d3(:,2) d3(:,2)-d3(:,3) d3(:,3)-d3(:,4) d3(:,4)];
rr=[rr; d4];
end
save(strcat(pathd, 'dataMatlab\',num2str(j_a),'.mat'),'rr');
end
read_abed_z19data
save(strcat(pathd, 'dataMatlab\mean_4.mat'),'mean_4');
%-- 2021/6/3 9:48 --%
popsize=6
S28a=[]; id=1; %k=0; a11=1; a12=1/2; a21=0; a22=sqrt(3)/2;
for i=0:1:popsize % 6 100
for j=0:1:popsize-i
for k=0:1:popsize-i-j
id=id+1;
S28a=[S28a; i j k popsize-i-j-k  ...
nchoosek(popsize,i) * nchoosek(popsize-i,j) ...
* nchoosek(popsize-i-j,k) id];
end
end
end
popsize=2
S28a=[]; id=1; %k=0; a11=1; a12=1/2; a21=0; a22=sqrt(3)/2;
for i=0:1:popsize % 6 100
for j=0:1:popsize-i
for k=0:1:popsize-i-j
S28a=[S28a; i j k popsize-i-j-k  ...
nchoosek(popsize,i) * nchoosek(popsize-i,j) ...
* nchoosek(popsize-i-j,k) id];
id=id+1;
end
end
end
sum(S28)
sum(S28a)
read_abed_z19data
sum(S28a)
clear
read_abed_z19data
sum(S28a)
1024*1024
read_abed_z19data
clear
read_abed_z19data
m0=mean_4(2,1:4);
m=round(m0*100);
clear
read_abed_z19data
nchoosek(popsize,i)
nchoosek(popsize-i,j)
nchoosek(popsize-i-j,k)
nchoosek(popsize,i) * nchoosek(popsize-i,j) ...
* nchoosek(popsize-i-j,k)
read_abed_z19data
m(3)-10
read_abed_z19data
[v c]=max(S28a(:,5))
semilogy(S28a(:,5))
read_abed_z19data
semilogy(S28a(1:7512,5))
semilogy(S28a(:,5))
read_abed_z19data
plot(S28a(:,end))
sum(S28a(:,7))
read_abed_z19data
sum(S28a(:,7))
read_abed_z19data
sum(S28a(:,7))
read_abed_z19data
sum(S28a(:,7))
read_abed_z19data
sum(S28a(:,7))
plot(p_mean(:,end))
%-- 2021/6/3 13:18 --%
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\p_mean.mat')
sum(S28a(:,7))
sum(p_mean(:,7))
read_abed_z19data
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\p_mean.mat')
read_abed_z19data
strcat(pathd, 'dataMatlab\p_mean' & num2str(a(parameter_a)) &
'.mat')
strcat(pathd, 'dataMatlab\p_mean' & num2str(a(parameter_a)) & '.mat')
save(strcat(pathd, 'dataMatlab\p_mean' & num2str(a(parameter_a)) , '.mat'),'p_mean');
save(strcat(pathd, 'dataMatlab\p_mean' , num2str(a(parameter_a)) , '.mat'),'p_mean');
read_abed_z19data
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\p_mean4.mat')
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\p_mean4.mat')
plot(p_mean(:,end))
clear
read_abed_z19data
cleat
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1.mat')
[C,ia,ic] = unique(rr(:,1:4),'rows')
histc(C)
histc(ic,13637)
histc(ic,13636)
histc(ic,1:end)
histc(ic,1:13637)
sum(histc(ic,1:13637))
sum(ia)
C100=round(C*100)
b = histc(ic,1:13637)
d=[C100 b];
histogram(d(:,end))
plot(d(:,end))
find(C100(C100(:,1)==13 & C100(:,2)==20 & C100(:,1)==18 & C100(:,1)==49),: ) % [13,20,18,49,239]
C100(find(C100(:,1)==13 & C100(:,2)==20 & C100(:,1)==18 & C100(:,1)==49),: ) % [13,20,18,49,239]
find(C100(:,1)==13
find(C100(:,1)==13)
C100(find(C100(:,1)==13 & C100(:,2)==20 ),: ) % [13,20,18,49,239]
C100(find(C100(:,1)==13 & C100(:,2)==20 & C100(:,3)==18 & C100(:,4)==49),: ) % [13,20,18,49,239]
C100=round(rr*100); C100(find(C100(:,1)==13 & C100(:,2)==20 & C100(:,3)==18 & C100(:,4)==49),: ) % [13,20,18,49,239]
sum(d(:,5))
C100=round(C*100);b = histc(ic,size(C,1))
C100=round(C*100);b = histc(ic,1:size(C,1))
[C,ia,ic] = unique(rr(:,1:4),'rows');
C100=round(C*100);b = histc(ic,1:size(C,1)); state_count=[C b];
sum(state_count(:,5))
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo4.mat')
d=abed_count_theo(find(abed_count_theo(:,6)==2),:);
d2=[d round(d(:,9)*500005) ];
d2(:,11) = d2(:,8) - d2(:,10);
e=d2(find(abs(d2(:,11))>100),:);
for k=1:size(e,1); if e(k,11)>0; colo=[1 0 0];else; colo=[0 0 1];end;
scatter3(e(k,1), e(k,2), e(k,3),abs(e(k,11))-100, 'o','MarkerFaceColor',colo);hold on;end
%-- 2021/6/3 18:22 --%
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo4.mat')
d=abed_count_theo(find(abed_count_theo(:,6)==2),:);
d2=[d round(d(:,9)*500005) ];
d2(:,11) = d2(:,8) - d2(:,10);
e=d2(find(abs(d2(:,11))>100),:);
for k=1:size(e,1); if e(k,11)>0; colo=[1 0 0];else; colo=[0 0 1];end;
scatter3(e(k,1), e(k,2), e(k,3),abs(e(k,11))-100, 'o','MarkerFaceColor',colo);hold on;end
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo4.mat')
d=abed_count_theo(find(abed_count_theo(:,6)==2),:);
d2=[d round(d(:,9)*500005) ];
d2(:,11) = d2(:,8) - d2(:,10);
e=d2(find(abs(d2(:,11))>50),:);
for k=1:size(e,1); if e(k,11)>0; colo=[1 0 0];else; colo=[0 0 1];end;
scatter3(e(k,1), e(k,2), e(k,3),abs(e(k,11))-100, 'o','MarkerFaceColor',colo);hold on;end
clear;load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo4.mat')
d=abed_count_theo(find(abed_count_theo(:,6)==2),:);
d2=[d round(d(:,9)*500005) ];
d2(:,11) = d2(:,8) - d2(:,10);
e=d2(find(abs(d2(:,11))>50),:);
for k=1:size(e,1); if e(k,11)>0; colo=[1 0 0];else; colo=[0 0 1];end;
scatter3(e(k,1), e(k,2), e(k,3),abs(e(k,11))-49, 'o','MarkerFaceColor',colo);hold on;end
saveas(gca,'C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\a11_50.fig')
clear
read_abed_z19data
uiopen('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_theo0.25_50.fig',1)
figure
uiopen('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_theo4_50.fig',1)
figure
uiopen('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_theo1_50.fig',1)
%-- 2021/6/3 22:38 --%
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo4.mat')
read_abed_z19data
for i in range(posr.shape[0]):
xyz = (posr-nash_equilibrium)[i, 0:3]
J_posr  = J_posr + probr[i] * (xyz.dot(xyz) * E - np.transpose([xyz]).dot([xyz]))
w, v = np.linalg.eig(J_posr)
read_abed_z19data
load(strcat(pathd, 'dataMatlab\mean_4.mat'),'mean_4');
sum(d(:,8))
read_abed_z19data
[v d] = eig(J_posr)
[v d] = eig(eye(4)-J_posr)
read_abed_z19data
xyzw(1:3).*xyzw(1:3)
eye(3)*xyzw(1:3).*xyzw(1:3)
eye(3)*(xyzw(1:3).*xyzw(1:3))
(xyzw(1:3).*xyzw(1:3))
eye(3).*(xyzw(1:3).*xyzw(1:3))
eye(3).*(xyzw(1:3).*xyzw(1:3))'
(xyzw(1:3).*xyzw(1:3))'
eye(3).*(xyzw(1:3).*xyzw(1:3))'
eye(3) *. (xyzw(1:3).*xyzw(1:3))'
(xyzw(1:3)'*xyzw(1:3))
(eye(3) *  (xyzw(1:3)'*xyzw(1:3))
eye(3) *  (xyzw(1:3)'*xyzw(1:3))
read_abed_z19data
(xyzw(1:3)'*xyzw(1:3))* eye(3)
(xyzw(1:3)'*xyzw(1:3)).* eye(3)
read_abed_z19data
(xyzw(1:3).*xyzw(1:3)')
xyzw(1:3)'  * eye(3)*xyzw(1:3)
xyzw(1:3)'
read_abed_z19data
(xyzw(1:3)'*xyzw(1:3)).* eye(3)
read_abed_z19data
d(k,1:3)/100
d(k,8)
read_abed_z19data
27626*3
clear
read_abed_z19data
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo4.mat')
for parameter_a = 1:3
d=abed_count_theo(find(abed_count_theo(:,6)==parameter_a),:);
d2=[d round(d(:,9)*500005) ];
d2(:,11) = d2(:,8) - d2(:,10);
e=d2(find(abs(d2(:,11))>50),:);
clf;
for k=1:size(e,1); if e(k,11)>0; colo=[1 0 0];else; colo=[0 0 1];end;
scatter3(e(k,1), e(k,2), e(k,3),abs(e(k,11))-49, 'o','MarkerFaceColor',colo);hold on;
end
saveas(gca,strcat( pathd ,'dataMatlab\abed_theo',num2str(a(parameter_a))  ,'_50.fig'))
end
read_abed_z19data
%-- 2021/6/4 9:05 --%
read_abed_z19data
xx=[0.447109341
-0.759107701
0.473126553
]
]
xx'*xx
read_abed_z19data
a(parameter_a)/(3* a(parameter_a) +1)
mean_4(parameter_a ,1:3)
read_abed_z19data
%-- 2021/6/4 23:02 --%
read_abed_z19data
%-- 2021/6/5 3:37 --%
0.057*0.083
%-- 2021/6/5 21:01 --%
read_abed_z19data
d3(k,8)
sum(d3(:,8))
sum(d2(:,8))
sum(d2(:,10))
read_abed_z19data
a=[0.25 1 4];
for j_a = 1:3
rr=[];
for k_session=1:5
abedfilesname=strcat(pathd, 'data\1pop4-100-', ...
num2str(a(j_a)),'-',num2str(k_session),'.csv');
[num] = xlsread(abedfilesname);
d3=num(15:end,[2:4:14 1]);
d4=[d3(:,1)-d3(:,2) d3(:,2)-d3(:,3) d3(:,3)-d3(:,4) d3(:,4)];
rr=[rr; d4 d3(:,5)];
end
save(strcat(pathd, 'dataMatlab\',num2str(j_a),'_',num2str(a(parameter_a)),'_abed.mat'),'rr');
end
mean(rr)
mean(rr(:,1:4))
%-- 2021/6/6 0:12 --%
A=[1	0	0
0	1	0
-1	0	0
0	-1	0
]
xyzw=A
J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3])
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3)) * 1 /4;
end
[V D]=eig(J_posr)
A=[1	0	0
0	2	0
-1	0	0
0	-1	0
]
J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3])
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3)) * 1 /4;
end
[V D]=eig(J_posr)
A=[0	1	0
0	0	2
0	-1	0
0	0	-1
]
J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3])
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3)) * 1 /4;
end
[V D]=eig(J_posr)
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\s_count_1_abed.mat')
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo_ne_0.25.mat')
a=abed_count_theo
sum(a(:,9))
sum(a(:,9).*a(:,8))
clear
pathd='C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\'
save(strcat(pathd, 'dataMatlab\abed_count_theo_ne2_' , num2str(a(parameter_a)) , '.mat'),'abed_count_theo');
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo_ne2_0.25.mat')
a=abed_count_theo
sum(a(:,9))
sum(a(:,9).*a(:,8))
clear
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo_ne2_0.25.mat')
sum(a(:,9))
a=abed_count_theo
sum(a(:,9))
read_abed_z19data
sum( s_count(:,5))
%-- 2021/6/6 1:10 --%
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed_count_theo_ne2_0.25.mat')
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\s_count_1_abed.mat')
sum(s_count(5))
sum(s_count(:,5))
A=s_count; J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3])
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3)) * A(k,5) /500005;
end
[V D]=eig(J_posr)
clear;load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\s_count_4_abed.mat')
sum(s_count(:,5))
A=s_count; J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3])
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3)) * A(k,5) /500005;
end
[V D]=eig(J_posr)
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\p_mean_1.mat')
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
A=rr; J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3]);
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3))  ;
end
J_posr
[V D]=eig(J_posr)
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\2_abed.mat')
A=rr; J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3]);
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3))  ;
end
J_posr
[V D]=eig(J_posr)
clear;load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\3_abed.mat')
A=rr; J_posr=0; for k =1:size(A,1)
xyzw = A(k,[1 2 3]);
J_posr = J_posr + ((xyzw(1:3)*xyzw(1:3)')*eye(3)  ...
- xyzw(1:3)'*xyzw(1:3))  ;
end
J_posr
[V D]=eig(J_posr)
%-- 2021/6/6 9:38 --%
read_abed_z19data
a=[];for i1=0:10;for i2=0:10-i1;for i3=0:10-i1-i2; i4=10-i1-i2-i3;
a=[a;i1 i2 i3 i4]; end;end;end
s_count=a;
for L=1:size(s_count,1)
x=s_count(L,1:4)
nck = nchoosek(popsize,x(1)) * nchoosek(popsize-x(1),x(2)) ...
* nchoosek(popsize-x(1)-x(2),x(3));
abed_count_theo=[abed_count_theo; x  nck  parameter_a a(parameter_a)  s_count(L,5) ...
nck * m0(1)^x(1) * m0(2)^x(2) * m0(3)^x(3) * m0(4)^x(4)];
end
for L=1:size(s_count,1)
x=s_count(L,1:4)
nck = nchoosek(popsize,x(1)) * nchoosek(popsize-x(1),x(2)) ...
* nchoosek(popsize-x(1)-x(2),x(3));
abed_count_theo=[abed_count_theo; x  nck  parameter_a a(parameter_a)  ...
nck * m0(1)^x(1) * m0(2)^x(2) * m0(3)^x(3) * m0(4)^x(4)];
end
aaa=[];for i1=0:10;for i2=0:10-i1;for i3=0:10-i1-i2; i4=10-i1-i2-i3;
aaa=[aaa;i1 i2 i3 i4]; end;end;end
s_count=aaa;
a=[0.25 1 4]; popsize=100;
for L=1:size(s_count,1)
x=s_count(L,1:4)
nck = nchoosek(popsize,x(1)) * nchoosek(popsize-x(1),x(2)) ...
* nchoosek(popsize-x(1)-x(2),x(3));
abed_count_theo=[abed_count_theo; x  nck  parameter_a a(parameter_a)  ...
nck * m0(1)^x(1) * m0(2)^x(2) * m0(3)^x(3) * m0(4)^x(4)];
end
abed_count_theo=[];
for L=1:size(s_count,1)
x=s_count(L,1:4)
nck = nchoosek(popsize,x(1)) * nchoosek(popsize-x(1),x(2)) ...
* nchoosek(popsize-x(1)-x(2),x(3));
abed_count_theo=[abed_count_theo; x  nck  parameter_a a(parameter_a)  ...
nck * m0(1)^x(1) * m0(2)^x(2) * m0(3)^x(3) * m0(4)^x(4)];
end
popsize=10
for L=1:size(s_count,1)
x=s_count(L,1:4)
nck = nchoosek(popsize,x(1)) * nchoosek(popsize-x(1),x(2)) ...
* nchoosek(popsize-x(1)-x(2),x(3));
abed_count_theo=[abed_count_theo; x  nck  parameter_a a(parameter_a)  ...
nck * m0(1)^x(1) * m0(2)^x(2) * m0(3)^x(3) * m0(4)^x(4)];
end
abed_count_theo=[];
for L=1:size(s_count,1)
x=s_count(L,1:4)
nck = nchoosek(popsize,x(1)) * nchoosek(popsize-x(1),x(2)) ...
* nchoosek(popsize-x(1)-x(2),x(3));
abed_count_theo=[abed_count_theo; x  nck  parameter_a a(parameter_a)  ...
nck * m0(1)^x(1) * m0(2)^x(2) * m0(3)^x(3) * m0(4)^x(4)];
end
sum(abed_count_theo(:,8))
plot(abed_count_theo(:,end))
clear
read_abed_z19data
%-- 2021/6/6 12:33 --%
a=[0.211479852014798 0.216014498550143 0.179746325367461 0.392759324067599
0.268911908809119 0.255113888611138 0.228454754524550 0.247519448055196
0.310324067593239 0.290094290570943 0.311819418058194 0.0877622237776218]
plot(a,'s')
plot(a,'s-')
xlabel('0.25','1','4')
xlabel(['0.25','1','4'])
xtick(['0.25','1','4'])
set(gca,'xtick', ['0.25','1','4'])
log(0.25)
log(0.25,2)
log(0.25)/log(2)
log(4)/log(2)
log(a)/log(2) + 1
set(gca,'xticklabel', ['0.25','1','4'])
set(gca,'xticklabel', {'0.25','1','4'})
clf
plot(a,'s-')
set(gca,'xticklabel', {'0.25','1','4'})
plot(a,'s-','xticklabel', {'0.25','1','4'})
ShujieMovie
plot(rr,'S-','linewidth',2);
plot(rr,'S-','linewidth',2,'fontsize',6);
plot(rr,'S-','linewidth',2,'markersize',6);
plot(rr,'S-','linewidth',2,'markersize',10);
set(gca,'xticklabel', {'','',''})
ShujieMovie
set(gca,'fontsize', 15)
'xticklabel', {'0.25','','1','','4'}
set(gca,'xticklabel',  {'0.25','','1','','4'})
set(gca,'xticklabel',  {'0.25','','1','','4'},'xlabel','a')
xlabel('a')
ylabel('p')
ledgen on
legend on
ylabel('\it{p}')
xlabel('\it{a}')
plot(rr,'S-','linewidth',2,'markersize',10); legend on
K>> ylabel('\it{p}')
K>> xlabel('\it{a}')
set(gca,'xticklabel',  {'0.25','','1','','4'})
plot(rr,'S-','linewidth',2,'markersize',10); legend on
ylabel('\it{p}'); xlabel('\it{a}')
set(gca,'xticklabel',  {'0.25','','1','','4'})
plot(rr,'S-','linewidth',2,'markersize',10);
ylabel('\it{p}'); xlabel('\it{a}')
set(gca,'xticklabel',  {'0.25','','1','','4'})
legend('x_1','x_1','x_1','x_1')
plot(rr,'S-','linewidth',2,'markersize',10);
ylabel('\it{p}'); xlabel('\it{a}')
set(gca,'xticklabel',  {'0.25','','1','','4'})
legend('x_1','x_2','x_3','x_4')
set(gca,'xtick', 1:1:3)
abedfilesname='C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-06\a4_0522simple.csv';
[num,txt,raw] = xlsread(abedfilesname);
rr=[]; for k=0:2; d=num(35:end,[3:6]+5*k);
rr=[rr; mean(d(1:end,:))];
%     D(:,:,k+1)=d;
end;
plot(rr,'S-','linewidth',2,'markersize',10);
ylabel('\it{p}'); xlabel('\it{a}')
legend('x_1','x_2','x_3','x_4')
set(gca,'fontsize', 15)
plot(rr,'S-','linewidth',2,'markersize',10);
ylabel('\it{p}'); xlabel('\it{a}')
legend('x_1','x_2','x_3','x_4')
set(gca,'fontsize', 15)
set(gca,'xticklabel',  {'0.25','','1','','4'})
%-- 2021/6/6 17:14 --%
ShujieMovie
plot(rr,'S-','linewidth',2,'markersize',10);
ylabel('\it{p}','fontsize', 15); xlabel('\it{a}')
legend('x_1','x_2','x_3','x_4','Orientation','horizontal')
set(gca,'fontsize', 15)
set(gca,'xticklabel',  {'0.25','','1','','4'})
plot(rr,'S-','linewidth',2,'markersize',10);
ylabel('\it{p}','fontsize', 20); xlabel('\it{a}','fontsize', 20)
legend('x_1','x_2','x_3','x_4','Orientation','horizontal')
set(gca,'fontsize', 15)
set(gca,'xticklabel',  {'0.25','','1','','4'})
plot(rr,'S-','linewidth',2,'markersize',10);
legend('x_1','x_2','x_3','x_4','Orientation','horizontal')
set(gca,'fontsize', 15)
ylabel('\it{p}','fontsize', 20); xlabel('\it{a}','fontsize', 20)
set(gca,'xticklabel',  {'0.25','','1','','4'})
plot(rr,'S-','linewidth',2,'markersize',10);
legend('1','2','3','4','Orientation','horizontal')
set(gca,'fontsize', 15)
ylabel('\it{x}','fontsize', 20); xlabel('\it{a}','fontsize', 20)
set(gca,'xticklabel',  {'0.25','','1','','4'})
ylim([0.1 0.5])
ShujieMovie
ylim([0 0.5])
plot(rr,'S-','linewidth',2,'markersize',10);
legend('1','2','3','4','Orientation','horizontal')
set(gca,'fontsize', 15)
ylabel('\it{x}','fontsize', 20); xlabel('\it{a}','fontsize', 20)
set(gca,'xticklabel',  {'1/4','','1','','4'})
ylim([0 0.5])
xlim([0.8 2.2])
plot(rr,'S-','linewidth',2,'markersize',10);
legend('1','2','3','4','Orientation','horizontal')
set(gca,'fontsize', 15)
ylabel('\it{x}','fontsize', 20); xlabel('\it{a}','fontsize', 20)
set(gca,'xticklabel',  {'1/4','','1','','4'})
ylim([0 0.5])
%-- 2021/6/7 15:45 --%
read_abed_z19data
xyzw(1:3)*xyzw(1:3)')*eye(3)
(xyzw(1:3)*xyzw(1:3)')*eye(3)
(xyzw(1:3).*xyzw(1:3))*eye(3)
eye(3).*(xyzw(1:3)'.*xyzw(1:3))'
(xyzw(1:3)'.*xyzw(1:3))
(xyzw(1:3).*xyzw(1:3))'
eye(3).*(xyzw(1:3).*xyzw(1:3))'
eye(3).*(xyzw(1:3)'*xyzw(1:3))
read_abed_z19data
xyzw(1:3)'*xyzw(1:3)
read_abed_z19data
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
tmp0607
[a(parameter_a)/(3* a(parameter_a) +1)*[1 1 1 1/a(parameter_a)]
]
mean(d3)
tmp0607
mean(d3)
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\3_abed.mat')
tmp0607
mean(d3)
clear
tmp0607
clear
tmp0607
clear
tmp0607
clear
tmp0607
plot(d3(1:48,:),'DisplayName','d3(1:48,:)')
plot(d3,'DisplayName','d3')
tmp0607
clear
tmp0607
eig_matrix=[eig_matrix; J_posr ve de  (Me4(1:3))' Ne4(1:3)']
tmp0607
clear
tmp0607
%-- 2021/6/8 1:29 --%
tmp0607
t1
(  (xyzw([1 2 3])*xyzw([1 2 3])')*eye(3)  ...
%                             - xyzw([1 2 3])'*xyzw([1 2 3]))
(  (xyzw([1 2 3])*xyzw([1 2 3])')*eye(3)  ...
- xyzw([1 2 3])'*xyzw([1 2 3]))
tmp0607
Oneill_replicator
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr)
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\3_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3]))
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[2 3 4]))
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\2_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\3_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\3_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
%-- 2021/6/8 9:09 --%
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\2_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\2_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\3_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\1_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
inertia_tensor
sqre(inertia_tensor)
sqrt(inertia_tensor)
eig(sqrt(inertia_tensor))
[v d]=eig(sqrt(inertia_tensor))
plot(v(:,3))
zz=sqrt(inertia_tensor)+sqrt(inertia_tensor)';[v d]=eig(zz)
v(:,4)/sum(v(:,4))
sqrt(v(:,4))/sqrt(sum(v(:,4)))
sqrt(v(:,4).*v(:,4))/(sum(v(:,4))^2)
(v(:,4).*v(:,4))/(sum(v(:,4))^2)
(v(:,4).*v(:,4))/sum(v(:,4).*v(:,4))
zz=sqrt(inertia_tensor)-sqrt(inertia_tensor)';[v d]=eig(zz)
zz
zz=sqrt(inertia_tensor)-sqrt(inertia_tensor)';[v d]=eig(zz(1:3,1:3))
zz=sqrt(inertia_tensor)+sqrt(inertia_tensor)';[v d]=eig(zz(1:3,1:3))
zz=sqrt(inertia_tensor)+sqrt(inertia_tensor)';[v d]=eig(zz(1:4,1:4))
zz= (inertia_tensor);[v d]=eig(zz(1:4,1:4))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\3_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
sum(v)
sum(ve)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
sum(ve)
ve.*ve
sum(ans)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\2_abed.mat')
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(rr(:,[1 2 3 4]))
tmp0607
%-- 2021/6/9 18:26 --%
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\s_count_1_abed.mat')
ckear
clear
read_abed_z19data
tmp0609_netcurrent
sum(ia)
sum(ic)
tmp0609_netcurrent
rr(2:lrr-1,1:4)
tmp0609_netcurrent
C(k,:)
tmp0609_netcurrent
C(k,:)
C(k,:) length(stateindex1) state1_velocity state1_jumpout state1_jumpin parameter_a
length(stateindex1)
C(k,:) length(stateindex1) state1_velocity state1_jumpout state1_jumpin parameter_a
length(stateindex1)
C(k,:) length(stateindex1)
length(stateindex1)
tmp0609_netcurrent
clear
tmp0609_netcurrent
plot(c18(:,5))
pathd='C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\';
a=[0.25 1 4];
save(strcat(pathd, 'dataMatlab\abed__velocity_count_' , num2str(a(parameter_a)) , '.mat'),'c18');
tmp0609_netcurrent
(state1_jumpout + state1_jumpin)/2
tmp0609_netcurrent
(rr(stateindex1, 1:dim) - rr(stateindex0, 1:dim))
tmp0609_netcurrent
clear
tmp0609_netcurrent
save(strcat(pathd, 'dataMatlab\abed__velocity_count_' , num2str(a(parameter_a)) , '.mat'),'c18','rs');
tmp0609_netcurrent
clear
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed__velocity_count_4.mat')
rs=[];
plot(c18(:,5))
c300=c18(find(c19(:,5)>300) ,:)
c300=c18(find(c18(:,5)>300) ,:);
c200=c18(find(c18(:,5)>200) ,:);
c400=c18(find(c18(:,5)>400) ,:);
%-- 2021/6/9 22:46 --%
tmp0609_netcurrent
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed__velocity_count_4.mat')
rs=[];
x=c18(:,1:3);
v=c18(:,6:8);
s=[0.3 0.3 0.3];
streamline(x,v,s)
streamline(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),s(:,1),s(:,2),s(:,3))
clear
load wind
[sx,sy,sz] = meshgrid(80,20:10:50,0:5:15);
streamline(stream3(x,y,z,u,v,w,sx,sy,sz))
view(3);
[X,Y,Z]=meshgrid(0.25:0.01:0.35,0.25:0.01:0.35,0.25:0.01:0.35);
quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3))
clear
pathd='C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\';
a=[0.25 1 4];
dim=4;
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed__velocity_count_4.mat')
rs=[];
x=c18(:,1:3);
v=c18(:,6:8);
quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3))
quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),2)
clear
pathd='C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\';
a=[0.25 1 4];
dim=4;
load('C:\Users\Think\Desktop\2021毕设\曾弈久\NASH_DIFFERENCE\NASH_DIFFERENCE\dataMatlab\abed__velocity_count_4.mat')
rs=[];
x=c18(:,1:3);
v=c18(:,6:8);
quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),5)
plot(c18(:,5))
c182=c18(find(c18(:,5)>400) ,:);
x=c182(:,1:3);
v=c182(:,6:8);
quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),5)
c182=c18(find(c18(:,5)>400) ,:);
x=c182(:,1:3);
v=c182(:,6:8);
quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),2)
c182=c18(find(c18(:,5)>600) ,:);
x=c182(:,1:3);
v=c182(:,6:8);
figure(6);quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),2)
rs=[];
c182=c18(find(c18(:,5)>400) ,:);
x=c182(:,1:3);
v=c182(:,6:8);
figure(6);quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),1)
c182=c18(find(c18(:,5)>200) ,:);
x=c182(:,1:3);
v=c182(:,6:8);
figure(6);quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),1)
% s=[0.3 0.3 0.3];
[X,Y,Z]=meshgrid(0.25:0.01:0.35,0.25:0.01:0.35,0.25:0.01:0.35);
streamline(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),s(:,1),s(:,2),s(:,3))
%-- 2021/6/10 3:51 --%
[X,Y,Z]=meshgrid(0.25:0.01:0.35,0.25:0.01:0.35,0.25:0.01:0.35);
tmp0609_netcurrent
(0.25:0.01:0.35 - 0.25)*100 + 1
0.25:0.01:0.35 - 0.25
0.25:0.01:0.35
[0.25:0.01:0.35] - 0.25
([0.25:0.01:0.35] - 0.25)*100 + 1
tmp0609_netcurrent
size(t,1)
clear
tmp0609_netcurrent
clear
tmp0609_netcurrent
([0.25:0.01:0.35] - 0.25)*100 + 1
tmp0609_netcurrent
startx=0.3;
starty=0.25:0.01:0.35;
startz=0.25:0.01:0.35;
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
[startx,starty,startz]= meshgrid(0.27,0.25:0.01:0.35,0.25:0.01:0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
[startx,starty,startz]= meshgrid(0.27:0.03:0.35,0.25:0.01:0.35,0.25:0.01:0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
[startx,starty,startz]= meshgrid(0.27,0.3,0.25:0.01:0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
tmp0609_netcurrent
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
[startx,starty,startz]= meshgrid(0.27,0.3,0.25:0.02:0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz,'r')
[startx,starty,startz]= meshgrid(0.27,0.3,0.25:0.02:0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
tmp0609_netcurrent
clear
tmp0609_netcurrent
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.27,0.27,0.25:0.02:0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
[startx,starty,startz]= meshgrid(0.27,0.33,0.25:0.02:0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
tmp0609_netcurrent
clear
tmp0609_netcurrent
%-- 2021/6/10 8:01 --%
tmp0609_netcurrent
clear
tmp0609_netcurrent
round((e_z - s_z)*100 + 1)
tmp0609_netcurrent
clear
tmp0609_netcurrent
round((intv_y - s_y)*100 + 1)
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.35, 0.35);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
[startx,starty,startz]= meshgrid(0.25,0.35, 0.25);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
[startx,starty,startz]= meshgrid(0.25,0.35, 0.3);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
tmp0609_netcurrent
v=c182(:,11:14);
[startx,starty,startz]= meshgrid(0.25,0.35, 0.3);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
plot3([startx,starty,startz],'ro')
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.3);
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.4); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.3); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.27); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.24); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.32); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.25); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.26); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
clf
[startx,starty,startz]= meshgrid(0.25,0.35, 0.26); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.35, 0.26); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.28); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.30); ;hold on
streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.30); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.35, 0.30); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.28); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.30); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.26); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.2:0.04:0.32); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro')
[startx,starty,startz]= meshgrid(0.25,0.35, 0.2:0.04:0.32); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
[startx,starty,startz]= meshgrid(0.25,0.35, 0.2:0.04:0.32); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
[startx,starty,startz]= meshgrid(0.25,0.35,  0.32); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
[startx,starty,startz]= meshgrid(0.25,0.35,  0.31); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.35,  0.31); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
figure(6);quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),1)
plot(v,'DisplayName','v')
c182=c18(find(c18(:,5)>2) ,:);
x=c182(:,1:3);
v=c182(:,6:8);
figure(6);quiver3(x(:,1),x(:,2),x(:,3),v(:,1),v(:,2),v(:,3),1)
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.35,  0.31); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
[startx,starty,startz]= meshgrid(0.25,0.30,  0.31); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.30,  0.31); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.30,  0.31); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
[startx,starty,startz]= meshgrid(0.25,0.30,  0.27); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
round((intv_y - s_y)*100 + 1)
round((intv_x - s_x)*100 + 1)
tmp0609_netcurrent
[startx,starty,startz]= meshgrid(0.25,0.30,  0.27); ;hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
plot3( startx,starty,startz ,'ro');hold on
[startx,starty,startz]= meshgrid(0.01:0.03:0.30,0.01:0.03:0.30,0.01:0.03:0.30); hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
[startx,starty,startz]= meshgrid(0.01:0.03:0.50,0.01:0.03:0.50,0.01:0.03:0.30); hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
[startx,starty,startz]= meshgrid(0.01:0.03:0.50,0.01:0.03:0.01,0.01:0.03:0.30); hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
[startx,starty,startz]= meshgrid(0.01:0.03:0.50,0.01:0.03:0.51,0.01:0.03:0.30); hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
[startx,starty,startz]= meshgrid(0.01:0.1:0.50,0.01:0.1:0.51,0.01:0.03:0.30); hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
[startx,starty,startz]= meshgrid(0.01:0.1:0.50,0.01:0.1:0.51,0.01:0.03:0.30); hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz,'linecolor','red');hold on
[startx,starty,startz]= meshgrid(0.01:0.1:0.50,0.01:0.1:0.51,0.01:0.03:0.30); hold on
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz);hold on
set('h','linecolor','red')
set(gca,'linecolor','red')
h
h.line(2)
h.Line
s=h;for i = 1:length(startx)
s(i).Color = [0.8 0 0];
s(i).LineStyle = '--';
s(i).LineWidth = 1;
end
[startx,starty,startz]= meshgrid(0.01:0.1:0.50,0.01:0.1:0.51,0.01:0.03:0.30);
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
s=h;for i = 1:length(startx)
s(i).Color = [0.8 0 0];
s(i).LineStyle = '--';
s(i).LineWidth = 1;
end
[startx,starty,startz]= meshgrid(0.01:0.1:0.50,0.2:0.1:0.51,0.01:0.03:0.30);
h=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
s=h;for i = 1:length(startx)
s(i).Color = [0.8 0 0];
s(i).LineStyle = '--';
s(i).LineWidth = 1;
end
[startx,starty,startz]= meshgrid(0.01:0.1:0.50,0.2:0.1:0.51,0.01:0.03:0.30);
s=streamline(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
for i = 1:length(startx)
s(i).Color = [0.8 0 0];
s(i).LineStyle = '--';
s(i).LineWidth = 1;
end
[startx,starty,startz]= meshgrid(0.01:0.1:0.50,0.2:0.1:0.51,0.01:0.03:0.30);
s=streamslice(Xx,Xy,Xz,Vx,Vy,Vz,startx,starty,startz)
for i = 1:length(startx)
s(i).Color = [0.8 0 0];
s(i).LineStyle = '--';
s(i).LineWidth = 1;
end
%-- 2021/6/11 18:46 --%
DV=[0 -1 0; 0 0 -1; -1 0 0]
DV=[0 -1 0; 0 0 -1; -1 0 0]*1/2
[v d] = eig(DV)
%-- 2021/6/12 13:46 --%
tmp202106122
plot(T,Y(:,1),'o')
hold on
plot(T,Y(:,2),'o')
plot(T,Y(:,1),'o')
figure; plot(T,Y(:,1),'*')
figure; plot(T,Y(:,2),'o')
figure; plot(T,Y(:,1),'*');figure; plot(T,Y(:,2),'o');
close all
tmp202106122
figure; plot(T,Y(:,1),'*');figure; plot(T,Y(:,2),'o');
tmp202106122
figure; plot(T,Y(:,1),'*');figure; plot(T,Y(:,2),'o');
tmp202106122
close all
figure; plot(T,Y(:,1),'*');figure; plot(T,Y(:,2),'o');
tmp202106122
figure; plot(T,Y(:,1),'*');figure; plot(T,Y(:,2),'o');
clear
a4_theo_num_replicator
V_F
tmp202106122
plot(T,Y(:,1),'.')
mean(Y)
plot(T,Y(1:3000,1),'.')
v=1:300;plot(T(v),Y(v,1),'.')
v=[1:300];plot(T(v),Y(v,1),'.')
v=[1:300];plot(T([1:300]),Y([1:300],1),'.')
plot(T([1:300]),Y([1:300],1),'.')
plot(T([1:300]),Y([1:300],1),'.-')
plot(T([1:3000]),Y([1:3000],1),'.-')
plot(T([2001:3000]),Y([2001:3000],1),'.-')
plot3(Y([2001:3000],1), Y([2001:3000],2), Y([2001:3000],3),   '.-')
tmp202106122
plot3(Y([2001:3000],1), Y([2001:3000],2), Y([2001:3000],3),   '.-')
figure;plot3(Y([1:3000],1), Y([1:3000],2), Y([1:3000],3),   '.-')
tmp202106122
plot3(Y([2001:3000],1), Y([2001:3000],2), Y([2001:3000],3),   '.-')
figure;plot3(Y([1:3000],1), Y([1:3000],2), Y([1:3000],3),   '.-')
figure;plot3(Y([2001:6000],1), Y([2001:6000],2), Y([2001:6000],3),   '.-')
plot(Y(:,3:end),'DisplayName','Y(:,3:end)')
scatter(getcolumn(Y(:,3:end),1),getcolumn(Y(:,3:end),2))
scatter(getcolumn(Y(:,2:3),1),getcolumn(Y(:,2:3),2))
plotmatrix(Y)
mean(Y)
mean(Y(3000:end,:))
mean(Y(2000:end,:))
tmp202106122
mean(Y(3000:end,:))
mean(Y(30000:end,:))
tmp202106122
figure;plot3(Y([2001:6000],1), Y([2001:6000],2), Y([2001:6000],3),   '.-')
mean(Y(3000:end,:))
mean(Y(30000:end,:))
figure;plot3(Y([20001:60000],1), Y([20001:60000],2), Y([20001:60000],3),   '.-')
tmp202106122
plotmatrix(Y)
clear
tmp202106122
rr = [rr; T Y ones(size(T),1)]
rr = [rr; T Y ones(size(T,1),1)]
tmp202106122
plot(T,Y(:,1),'.-')
tmp202106122
close all
tmp202106122
YQM20200830
uiopen('F:\cDownload\PokerAnan20200716\8x8abedata\6.fig',1)
x5_YQM
syms x(5) real
a = sym('x',[5,1])
sym('x',[5,1])
sym('x',[5,1]) real
payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ;6 10 10 0 3 ; 10 10 10 4 0 ]
%payoff_matrix_1 = [0 -1 1 1 -1; 1 0 -1 -1 1;-1 1 0 1 -1;-1 1 -1 0 1;1 -1 1 -1 0 ]
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
sym('x',[5,1])
payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ;6 10 10 0 3 ; 10 10 10 4 0 ]
%payoff_matrix_1 = [0 -1 1 1 -1; 1 0 -1 -1 1;-1 1 0 1 -1;-1 1 -1 0 1;1 -1 1 -1 0 ]
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
V_Eq_0
syms x
payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ;6 10 10 0 3 ; 10 10 10 4 0 ]
%payoff_matrix_1 = [0 -1 1 1 -1; 1 0 -1 -1 1;-1 1 0 1 -1;-1 1 -1 0 1;1 -1 1 -1 0 ]
Payoff_vector_field_F = payoff_matrix *[x(1) x(2) x(3) x(4) x(5)]';
syms x(1:5)
payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ;6 10 10 0 3 ; 10 10 10 4 0 ]
%payoff_matrix_1 = [0 -1 1 1 -1; 1 0 -1 -1 1;-1 1 0 1 -1;-1 1 -1 0 1;1 -1 1 -1 0 ]
Payoff_vector_field_F = payoff_matrix *[x(1) x(2) x(3) x(4) x(5)]';
syms x([1:5])
YQM5
sym('x',[1:5])
sym('x', 1:5 )
syms x(1) x(2) x(3) x(4) x(5)
syms x %(1) x(2) x(3) x(4) x(5)
payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ;6 10 10 0 3 ; 10 10 10 4 0 ]
%payoff_matrix_1 = [0 -1 1 1 -1; 1 0 -1 -1 1;-1 1 0 1 -1;-1 1 -1 0 1;1 -1 1 -1 0 ]
Payoff_vector_field_F = payoff_matrix *[x(1) x(2) x(3) x(4) x(5)]';
mean_U = [x(1) x(2) x(3) x(4) x(5)] * Payoff_vector_field_F;
V_Eq_0 = [x(1) x(2) x(3) x(4) x(5)]'.*(Payoff_vector_field_F - mean_U);
YQM5
x(1)
x(2)
YQM5
sym('x', [5,1])
Payoff_vector_field_F = payoff_matrix *x;
Payoff_vector_field_F
sym('x', [5,1])
Payoff_vector_field_F = payoff_matrix *[x(1) x(2) x(3) x(4) x(5)]'
syms x [1 5] x
syms x [1 5]
YQM5
tmp202106122
mean(Y)
mean(Y(2000:end,:))
tmp202106122
mean(Y(2000:end,:))
tmp202106122
plot(rr(:,5))
plot(rr(:,1:5),'DisplayName','rr(:,1:5)')
m-Y.-mean(Y);
m=Y.-mean(Y);
m=Y. - mean(Y);
m=Y - mean(Y);
clear
[T,Y] = call_osc()
tmp202106122
[T,Y] = call_osc()
[T,Y] = tmp202106122()
%-- 2021/6/12 19:54 --%
tmp202106122
m=Y. - mean(Y);
m=Y - mean(Y);
[T,Y] = tmp202106122()
m=Y. - mean(Y);
m=Y - mean(Y);
m = Y. - mean(Y);
m = Y(:,1:4). - mean(Y(:,1:4));
m = Y(:,1:4) .- mean(Y(:,1:4));
m = Y(:,1:4) - ones(73477,4).*mean(Y(:,1:4));
m=[]; for k=1:73477; m = [m; Y(k,1:4)-mean(Y(:,1:4))];end
m=[]; dd= mean(Y(:,1:4)) ; for k=1:73477; m = [m; Y(k,1:4)-dd];end
plotmatrix(m(10000:end,:))
plotmatrix(m(10000:end,:),'o')
plotmatrix(m(100:end,:),'o')
m=[]; dd= mean(Y(1:23477,1:5)) ; for k=1:23477; m = [m; Y(k,1:5)-dd];end
plot(m(5559:5623,:),'DisplayName','m(5559:5623,:)')
tmp202106122
OneillYMQ20201118
V_Eq_0
tmp202106122
figure;plotmatrix(Y(400,:));
plot(Y(785:end,:),'DisplayName','Y(785:end,:)')
clear
tmp202106122
plot(Y(1785:end,:),'DisplayName','Y(1785:end,:)')
tmp202106122
plot(Y,'DisplayName','Y')
tmp202106122
plot(Y,'DisplayName','Y')
plotmatrix(Y,'.-')
clear
tmp202106122
%-- 2021/6/12 23:18 --%
[T,Y] = tmp202106122()
tmp202106122
plot(Y,'DisplayName','Y')
tmp202106122
plotmatrix(Y,'.-')
plot(Y,'DisplayName','Y')
tmp202106122
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
dos_x_t=Y;
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
close all
tmp202106122
[L_matrix m24x24 Tasym eigenvalue_set distr eigenvector_set]=from_timesies_figure8(Y)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(Y)
sum(ve)
tmp202106122
close all
tmp202106122
plotmatrix(Yr)
close all
tmp202106122
plotmatrix(Yr)
corr(Yr)
tmp202106122
plotmatrix(Yr)
plotmatrix(Y,'.-')
tmp202106122
plotmatrix(Yr)
corr(Yr)
tmp202106122
plotmatrix(Yr)
mean(Y)
tmp202106122
plotmatrix(Yr)
tmp202106122
mean(Y')
mean(Yr')
s = mean(Yr')
hist(s)
hist(s,10)
hist(s,50)
plot(s)
tmp202106122
s = mean(Yr')
hist(s,50)
tmp202106122
s = mean(Yr')
hist(s,50)
tmp202106122
s = mean(Yr')
hist(s,50)
tmp202106122
hist(s,50)
s = mean(Yr')
hist(s,50)
plot(Y,'DisplayName','Y')
plot(Y(:,[1,5]),'DisplayName','Y(:,[1,5])')
%-- 2021/6/13 20:49 --%
load('F:\Spectrum2020\datasource\oneill15327.mat'); d=data15327(:,4:11);
from_timesies_figure8(d);
from_timesies_figure2arxiv(dos_x_t)
from_timesies_figure2arxiv(d);
%-- 2021/6/14 0:16 --%
tmp202106122
s = mean(Yr')
tmp202106122
scatter(getcolumn(Yr(:,99:end),1),getcolumn(Yr(:,99:end),2))
scatter(getcolumn(Yr(:,86:87),1),getcolumn(Yr(:,86:87),2))
scatter(getcolumn(Yr(:,83:84),1),getcolumn(Yr(:,83:84),2))
scatter(getcolumn(Yr(:,80:81),1),getcolumn(Yr(:,80:81),2))
mean(Y')
tmp202106122
mean(Y')
mean(Y)
plot(Y,'DisplayName','Y')
plot(Y(:,6:7),'DisplayName','Y(:,6:7)')
plot(Y(:,6:end),'DisplayName','Y(:,6:end)')
plot(Y(:,2:4),'DisplayName','Y(:,2:4)')
scatter(getcolumn(Y(:,1:2),1),getcolumn(Y(:,1:2),2))
scatter(getcolumn(Y(:,[1,5]),1),getcolumn(Y(:,[1,5]),2))
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
tmp202106124
mean(Ym)
plot(Y(:,1:2),'DisplayName','Y(:,1:2)')
plot(Y(:,3:end),'DisplayName','Y(:,3:end)')
scatter(getcolumn(Y(:,3:end),1),getcolumn(Y(:,3:end),2))
plot(Y,'DisplayName','Y')
s = mean(Yr')
%-- 2021/6/16 15:08 --%
tmp0607
%-- 2021/6/24 11:06 --%
c1dot := -alpha*c1 + p1*(c1+beta1*c2)
p1dot := p1*(1-c1-beta1*c2)
c2dot := -alpha*c2 + p2*(c2+beta2*c1)
p2dot := p2*(1-c2-beta2*c1)
syms alpha c1 p1 beta1 beta2 c2 p2
c1dot  = -alpha*c1 + p1*(c1+beta1*c2) ;
p1dot  = p1*(1-c1-beta1*c2) ;
c2dot  = -alpha*c2 + p2*(c2+beta2*c1) ;
p2dot  = p2*(1-c2-beta2*c1) ;
eqpts  =
solve({c1dot=0,p1dot=0,c2dot=0,p2dot=0},{c1,p1,c2,
p2});
solve(c1dot=0,p1dot=0,c2dot=0,p2dot=0, c1,p1,c2,p2);,
solve([c1dot ,p1dot ,c2dot ,p2dot ])
A=solve([c1dot ,p1dot ,c2dot ,p2dot ])
A.p1
[ A.c1 A.c2 A.p1 A.p2]
Jacobian([c1dot,p1dot,c2dot,p2dot],[c1,p1,c2,p2])
jacobian([c1dot,p1dot,c2dot,p2dot],[c1,p1,c2,p2])
clear
syms x_1 x_2 x_3 x_4 a a_1 a_2
x_1  =  a/(3a+1) + a_1 * (a - 1)/4 + a_2
x_2  =   a/(3a+1) + a_1 * (- a - 1)/2
x_3  =   a/(3a+1) + a_1 * (a - 1)/4 - a_2
x_4  =   1/(3a+1) + a_1
x_1  =  a/(3*a+1) + a_1 * (a - 1)/4 + a_2
x_2  =   a/(3*a+1) + a_1 * (- a - 1)/2
x_3  =   a/(3*a+1) + a_1 * (a - 1)/4 - a_2
x_4  =   1/(3*a+1) + a_1
subs(x_2,'a_1','x_4 - 1/(3*a+1)')
x_1_2 = subs(x_1,'a_1','x_4 - 1/(3*a+1)')
x_1_3 = subs(x_1_2,'a_2','(x_1-x_3)/2')
x_2_2 = subs(x_2,'a_1','x_4 - 1/(3*a+1)')
%-- 2021/6/24 22:42 --%
a4_theo_num_replicator
latex(V_F)
x1=2*(- x3/2 + a/(3*a + 1) + ((x4 - 1/(3*a + 1))*(a - 1))/4)
vv=[-((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*(((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*(a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2)  - a*x4 + (a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2) *x3 + x3*x4 + a*((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*x4)
-(a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2) *(((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*(a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2)  - ((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2) + (a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2) *x3 + x3*x4 + a*((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*x4)
-x3*(((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*(a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2)  - (a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2)  + (a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2) *x3 + x3*x4 + a*((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*x4)
-x4*(((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*(a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2)  - x3 + (a/(3*a + 1) - ((x4 - 1/(3*a + 1))*(a + 1))/2) *x3 + x3*x4 + a*((2*a)/(3*a + 1) - x3 + ((x4 - 1/(3*a + 1))*(a - 1))/2)*x4)]
simplify(vv)
simplify(simplify(vv))
(- a^2*x4^2 + 2*a*x4 + x4^2 - 1)/(x4*(a*x4 - x4 + 1))
simplify((- a^2*x4^2 + 2*a*x4 + x4^2 - 1)/(x4*(a*x4 - x4 + 1)))
factor(ans)
factor((- a^2*x4^2 + 2*a*x4 + x4^2 - 1))
V_F
x1=a/(3*a + 1)+y1; x2=a/(3*a + 1)+y2; x3=a/(3*a + 1)+y3; x4=1/(3*a + 1)+y4;
clear
syms y1 y2 y3 y4
x1=a/(3*a + 1)+y1; x2=a/(3*a + 1)+y2; x3=a/(3*a + 1)+y3; x4=1/(3*a + 1)+y4;
V_F = [ -x1*(x1*x2 - a*x4 + x2*x3 + x3*x4 + a*x1*x4)
-x2*(x1*x2 - x1 + x2*x3 + x3*x4 + a*x1*x4)
-x3*(x1*x2 - x2 + x2*x3 + x3*x4 + a*x1*x4)
-x4*(x1*x2 - x3 + x2*x3 + x3*x4 + a*x1*x4)]
syms y1 y2 y3 y4 a
x1=a/(3*a + 1)+y1; x2=a/(3*a + 1)+y2; x3=a/(3*a + 1)+y3; x4=1/(3*a + 1)+y4;
V_F = [ -x1*(x1*x2 - a*x4 + x2*x3 + x3*x4 + a*x1*x4)
-x2*(x1*x2 - x1 + x2*x3 + x3*x4 + a*x1*x4)
-x3*(x1*x2 - x2 + x2*x3 + x3*x4 + a*x1*x4)
-x4*(x1*x2 - x3 + x2*x3 + x3*x4 + a*x1*x4)]
simplify(V_F)
%-- 2021/6/27 8:45 --%
tmp202106124
tmp0607
x5_YQM
jiacobian(V_Eq_0)
jacobian(V_Eq_0)
x1=444/2987;x2=22/103;x3=641/2987;x4=288/2987;x5=976/2987;
J = jacobian(V_Eq_0)
eval(J)
sum(J)
sum(eval(J))
sum(eval(J'))
%-- 2021/6/29 23:19 --%
x5_YQM
jacobian(V_Eq_0) - D_Eq_0
S =solve(V_Eq_0)
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5])
tmp202106129
find(s>0)
s=prod(A,2)
id = find(s>0)
x1=A(id,1); x2=A(id,2); x3=A(id,3); x4=A(id,4); x5=A(id,5);
D_Eq_at_NE = eval(Jac)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
tmp202106129
D_Eq_at_NE
sum(D_Eq_at_NE)
tmp202106129
plot(V(:,end))
tmp202106129
VI = imag(V)
payoff_matrix = [1	1	4	8	7
6	1	2	5	9
1	10	2	5	7
5	6	4	2	10
5	4	4	9	10
]
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U)/mean_U;
Jac = jacobian(V_Eq_0);
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);
s=prod(A,2);
syms x1 x2 x3 x4 x5 real
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U)/mean_U;
Jac = jacobian(V_Eq_0);
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);
s=prod(A,2);
id = find(s>0);
x1=A(id,1); x2=A(id,2); x3=A(id,3); x4=A(id,4); x5=A(id,5);
payoff_matrix
D_Eq_at_NE = eval(Jac);
[eigen_vector, eigen_value] = eig(D_Eq_at_NE);
R=[R; payoff_matrix];
V=[V; diag(eigen_value)'];
tmp202106129
VI = imag(V)
tmp202106129
VI = imag(V)
VI2 = sum(abs(sign(imag(V)))')'
find(VI2>3)
tmp202106129
diag(eigen_value)'
A(id,1:5)
tmp202106129
clear
tmp202106129
%-- 2021/6/30 10:54 --%
tmp202106129
bimat(payoff_matrix)
bimat(payoff_matrix,payoff_matrix')
tmp202106129
plot(eigen_vector(:,2))
tmp202106129
save('manifoldtry.mat', 'R','Ne')
load('manifoldtry.mat')
clear
load('manifoldtry.mat')
%-- 2021/6/30 19:06 --%
tmp202106129
clear
load('manifoldtry.mat')
%-- 2021/7/1 1:19 --%
tmp202106129
save('manifoldtry1.mat', 'R','Ne','-append')
save('manifoldtry1.mat', 'R','Ne' )
tmp202106129
save('manifoldtry6.mat', 'R','Ne' )
clear
payoff_matrix =[1,5,1,5,1;2,2,5,5,1;1,5,2,2,3;1,5,4,1,2;4,3,3,3,1]
syms x1 x2 x3 x4 x5 real
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U)/mean_U;
Jac = jacobian(V_Eq_0);
S=solve(V_Eq_0);
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);
s=prod(A,2);
id = find(s>0);
if size(id,1)==1
x1=A(id,1); x2=A(id,2); x3=A(id,3); x4=A(id,4); x5=A(id,5);
D_Eq_at_NE = eval(Jac);
[eigen_vector, eigen_value] = eig(D_Eq_at_NE);
end
load('manifoldtry6.mat', 'R','Ne' )
payoff_matrix =[4,1,3,4,4;5,3,3,2,2;5,4,1,2,3;4,5,3,2,1;2,5,4,4,1]
syms x1 x2 x3 x4 x5 real
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U)/mean_U;
Jac = jacobian(V_Eq_0);
S=solve(V_Eq_0);
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5]);
s=prod(A,2);
id = find(s>0);
if size(id,1)==1
x1=A(id,1); x2=A(id,2); x3=A(id,3); x4=A(id,4); x5=A(id,5);
D_Eq_at_NE = eval(Jac);
[eigen_vector, eigen_value] = eig(D_Eq_at_NE);
end
%-- 2021/7/2 5:16 --%
syms x1 x2 x3 x4 x5 real
payoff_matrix = [ 0 3 4 11 11 ; 5 0 2 11 12 ; 2 5 0 9 12 ; 6 10 10 0 3 ; 10 10 10 4 0 ]
%payoff_matrix_1 = [0 -1 1 1 -1; 1 0 -1 -1 1;-1 1 0 1 -1;-1 1 -1 0 1;1 -1 1 -1 0 ]
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
%         V_Eq_0 = [V_Eq]
%         V_Eq_0 = [[x1 x2 x3 x4 x5]'.*Payoff_vector_field_F - sum([x1 x2 x3 x4 x5]'.*Payoff_vector_field_F)]
%         V_Eq_0 = [ x1 x2 x3 x4 x5]'.*Payoff_vector_field_F ]
%         V_Eq_0 = mean_U
D_Eq_0 = jacobian(V_Eq_0)
x1=444/2987;x2=22/103;x3=641/2987;x4=288/2987;x5=976/2987;
D_Eq_at_NE = eval(D_Eq_0)
[eigen_vector eigen_value] = eig(D_Eq_at_NE)
clf
quiver([1:5]'*0,[1:5]'*0,real(eigen_vector(:,3)),imag(eigen_vector(:,3)),1)
v=eigen_vector(:,3);hold on
text(real(v)*1.05,imag(v)*1.05, num2str([1:5]'),'fontsize',15);axis square;
tmp202107021
clear
tmp202107021
s = bimat(payoff_matrix,payoff_matrix');
tmp202107021
id = find(s~= 0)
size(id,1)
tmp202107021
Jac_at_NE = eval(Jac)
tmp202107021
simplify(eigen_vector)
vap(eigen_vector)
vpa(eigen_vector)
vpa(eigen_value)
clear
tmp202107021
%-- 2021/7/2 19:02 --%
perm(1:5)
prem(1:5)
perms(1:5)
s=perms(1:5)
rand(5,[1 120])
randi(5,[1 120])
randi([1 120],5)
randi([1 120],[5 1])
t=randi([1 120],[5 1])
s(t,:)
payoff_matrix = s(t,:)
x = bimat(payoff_matrix,payoff_matrix');
tmp202107021
payoff_matrix = s(t,:)
clear
s=perms(1:5)
t=randi([1 120],[5 1])
payoff_matrix = s(t,:)
syms x1 x2 x3 x4 x5 a real
tmp202107021
x = bimat(payoff_matrix,payoff_matrix');
tmp202107021
clear
tmp202107021
vv=imag(V);
vv=[real(V)  imag(V)];
save('manifold5try.mat', 'R','V','Ne' )
%-- 2021/7/2 22:44 --%
load('manifold5try.mat')
vv=[real(V)  imag(V)];
r=R(76:80,:)
v4=(:,7:10);
v4=vv(:,7:10);
prod(v4)
prod(v4')
prod(v4')'
vv(:,11)=prod(v4')'
346*5+1
r=R(1731:1731+4,:)
run('D:\MATLAB\R2016a\bin\YQM5.m')
%-- 2021/7/3 15:26 --%
tmp0607
clear
fn = 'F:\cDownload\PokerAnan20200716\8x8abedata\ZSJ20201215\6A44.csv';
dos_x_t=load(fn)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
dos_x_t=dos_x_t4(:,1:3)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
dos_x_t=dos_x_t4(:,1:3)
dos_x_t4=load(fn)
dos_x_t=dos_x_t4(:,1:3)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
fn = 'F:\cDownload\PokerAnan20200716\8x8abedata\ZSJ20201215\6A11.csv';
dos_x_t4=load(fn);
dos_x_t=dos_x_t4(:,1:3);
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
tmp0607
fn = 'F:\cDownload\PokerAnan20200716\8x8abedata\ZSJ20201215\6A11.csv'
dos_x_t4=load(fn);
dos_x_t=dos_x_t4(:,1:3)
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
fn = 'F:\cDownload\PokerAnan20200716\8x8abedata\ZSJ20201215\6A25.csv';
dos_x_t4=load(fn);
dos_x_t=dos_x_t4(:,1:3);
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
fn = 'F:\cDownload\PokerAnan20200716\8x8abedata\ZSJ20201215\6A25.csv';
dos_x_t4=load(fn); dos_x_t=dos_x_t4(:,1:3);
[inertia_tensor, ve, de, Me] = from_timesies_out_inertia(dos_x_t)
Me=[0.1429 0.1429 0.1429]
tmp202107021
figure
v=diag(eigen_vector(:,1));hold on
quiver([1:5]'*0,[1:5]'*0,real(v),imag(v),1)
text(real(v)*1.05,imag(v)*1.05, num2str([1:5]'),'fontsize',15);axis square;
abs(v)
figure
v= (eigen_vector(:,1));hold on
quiver([1:5]'*0,[1:5]'*0,real(v),imag(v),1)
text(real(v)*1.05,imag(v)*1.05, num2str([1:5]'),'fontsize',15);axis square;
abs(v)
figure
v= (eigen_vector(:,3));hold on
quiver([1:5]'*0,[1:5]'*0,real(v),imag(v),1)
text(real(v)*1.05,imag(v)*1.05, num2str([1:5]'),'fontsize',15);axis square;
abs(v)
sin(2*pi/5)
sin(4*pi/5)
%-- 2021/7/13 14:53 --%
A = [.3, .1, .05, .2; .1, .2, .3, .3; .3, .5, .2, .3; .3, .2, .45, .2]
[V, D] = eig(A);
[V, D] = eig(A)
T = [277, 444, 14, 735, 1123, 35, 51, 1209, 1944; ...
587, 11148, 1884, 13619, 8174, 4497, 3934, 16605, 30224; ...
236, 2915, 1572, 4723, 11657, 430, 1452, 13539, 18262; ...
1100, 14507, 3470, 19077, 20954, 4962, 5437, 31353, 50430; ...
133, 2844, 676, 3653, 1770, 250, 273, 2293, 5946; ...
3, 134, 42, 179, -90, -177, 88, -179, 0; ...
-246, 499, 442, 695, 2675, 100, 17, 2792, 3487; ...
954, 12240, 13632, 26826, 0, 0, 0, 0, 26826; ...
844, 15717, 14792, 31353, 4355, 173, 378, 4906, 36259; ...
1944, 30224, 18262, 50430, 25309, 5135, 5815, 36259, 86689];
%-- 2021/7/14 16:09 --%
read312data
ff='C:\Users\Think\Documents\WeChat Files\wxid_b1np1m8smrk612\FileStorage\File\2021-06\abed-1pop Strategy distributions (recent history).csv'
[num,txt,raw] = xlsread(ff);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18 v.u22];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18-v.u22 v.u22];
abed_mean=mean(psd)
[num,txt,raw] = xlsread(ff);
t=length(txt(:,1));
u=raw(t+1:length(raw(:,1)),:);
v=cell2table(u);
ps1 = [v.u1 v.u2 v.u6 v.u10  v.u14 v.u18];
psd = [v.u2-v.u6 v.u6-v.u10  v.u10-v.u14 v.u14-v.u18 v.u18];
abed_mean=mean(psd)
%-- 2021/7/14 20:24 --%
a=[0.014204	0.028104	0.017376	0.023405	0.017253 0.018632 0.0192 0.118979 VI
0.025277	0.049922	0.030861	0.041568	0.030642 0.033091 0.0341 0.211311 V
0.263961	0.522355	0.322904	0.434939	0.32062 0.346239 0.3568 2.211018 I
0.136423	0.270108	0.166973	0.224906	0.165792 0.179039 0.1875 1.14331 III
0.055189	0.109214	0.067513	0.090937	0.067036 0.072392 0.0746 0.462281 IV
0.244726	0.484291	0.299374	0.403245	0.297257 0.321008 0.3308 2.049901 II
]
A=[0.014204	0.028104	0.017376	0.023405	0.017253 0.018632 0.0192 0.118979 VI
0.025277	0.049922	0.030861	0.041568	0.030642 0.033091 0.0341 0.211311 V
0.263961	0.522355	0.322904	0.434939	0.32062 0.346239 0.3568 2.211018 I
0.136423	0.270108	0.166973	0.224906	0.165792 0.179039 0.1875 1.14331 III
0.055189	0.109214	0.067513	0.090937	0.067036 0.072392 0.0746 0.462281 IV
0.244726	0.484291	0.299374	0.403245	0.297257 0.321008 0.3308 2.049901 II
]
A=[0.014204	0.028104	0.017376	0.023405	0.017253	0.018632
0.025277	0.049922	0.030861	0.041568	0.030642	0.033091
0.263961	0.522355	0.322904	0.434939	0.32062	0.346239
0.136423	0.270108	0.166973	0.224906	0.165792	0.179039
0.055189	0.109214	0.067513	0.090937	0.067036	0.072392
0.244726	0.484291	0.299374	0.403245	0.297257	0.321008
]
[v d] = eig(A)
v(:,1)/sum(v(:,1))
%-- 2021/7/15 13:19 --%
tmp202107021
v
t=0:0.1:10
x=v*exp(0.4294i *t);
xr=real(x);
xr=real(x)';
plot(xr,'DisplayName','xr')
t=0:0.1:30
x=v*exp(0.4294i *t);
xr=real(x)';
plot(xr,'DisplayName','xr')
scatter(getcolumn(xr(:,1:2),1),getcolumn(xr(:,1:2),2))
scatter(xr(:,1),xr(:,2),'-')
scatter(xr(:,1),xr(:,2),'.')
plot(xr(:,1),xr(:,2),'.')
plot(xr(:,1),xr(:,2),'-')
clear
load('manifold5try.mat')
%-- 2021/7/16 13:25 --%
load('manifold5try.mat')
riV = [real(V) imag(V)]
for k = 1:size(riV,1); riv(k,11)=riv(k,7)*riv(k,8)*riv(k,9)*riv(k,10) ;end;
for k = 1:size(riV,1); riv(k,11)=riV(k,7)*riV(k,8)*riV(k,9)*riV(k,10) ;end;
for k = 1:size(riV,1); riV(k,11)=riV(k,7)*riV(k,8)*riV(k,9)*riV(k,10) ;end;
for k = 1:size(riV,1);riV(k,12)=k; riV(k,12)=riV(k,7)*riV(k,8)*riV(k,9)*riV(k,10) ;end;
for k = 1:size(riV,1);riV(k,11)=k; riV(k,12)=riV(k,7)*riV(k,8)*riV(k,9)*riV(k,10) ;end
riV = round([real(V) imag(V)],4);
for k = 1:size(riV,1);riV(k,11)=k; riV(k,12)=riV(k,7)*riV(k,8)*riV(k,9)*riV(k,10) ;end
tmp202107021
clear
load('F:\cDownload\PokerAnan20200716\8x8abedata\5d_500000ticks_psd.mat')
[ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(psd,mean(psd))
clear
load('manifold5try.mat')
riV = [real(V) imag(V)]
for k = 1:size(riV,1);riV(k,11)=k; riV(k,12)=riV(k,7)*riV(k,8)*riV(k,9)*riV(k,10) ;end
%-- 2021/7/16 19:12 --%
load('F:\yqm代码+数据\5策略仿真simula\MSafter.csv')
dos_x_t = MSafter; [ Yret3,mn,am_eigencycleDim_t] = from_N_colExp_out_am(dos_x_t,mean(MSafter))
a2=[0.6945   -0.0031
-0.5687  -0.0114
-0.3742  0.0007
0.2483   0.0138
0.6926   -0.0056
-0.2306  0.0041
0.2326              -0.0016
0.562               -0.0025
-0.4381  -0.0144
-0.0428  0.0022]
scatter(Yret3, a2(:,1))
scatter(Yret3, a2(:,2))
scatter( a2(:,1), Yret3)
%-- 2021/7/18 9:17 --%
load('manifold5try.mat')
riV = [real(V) imag(V)]
for k = 1:size(riV,1);riV(k,11)=k; riV(k,12)=riV(k,7)*riV(k,8)*riV(k,9)*riV(k,10) ;end
find(riV(:,12)>0)
class0=riV(find(riV(:,7)==0 & riV(:,8)==0 & riV(:,9)==0 & riV(:,10)==0  ),:)
tmp1 = sum(sign(class0(:,2:5))')'
class0(:,13) = sum(sign(class0(:,2:5))')'
bimat([[1,5,3,2,4;1,2,3,5,4;3,2,4,5,1;1,3,4,2,5;2,3,1,5,4]])
A=([[1,5,3,2,4;1,2,3,5,4;3,2,4,5,1;1,3,4,2,5;2,3,1,5,4]])
bimat(A,A')
A=[[1,2,4,3,5;3,1,5,4,2;1,5,2,3,4;3,5,4,2,1;3,4,5,1,2]]
bimat(A,A')
class0(find(:,13)>0)
class0(find(class0(:,13)>0))
class0(find(class0(:,13)>0),:)
class0(find(class0(:,13)<0),:)
R(227*5-4:227*5,:)
R(53*5-4:53*5,:)
tmp1 = sum(sign(class0(:,2:5))')'
tmp2 = sum(sign(riV(:,2:5))')'
riV(:,13) = sum(sign(riV(:,2:5))')'
find(abs(riV(:,13))==4)
class4=riV(find(abs(riV(:,13))==4),:)
class4pn=riV(find(abs(riV(:,13))==0 & riV(:,12)>0),:)
tmp202106129
clear
load('manifold5try.mat')
bimat(;[5,2,3,4,1;4,3,1,2,5;4,2,3,5,1;3,5,2,4,1;3,2,5,4,1]')
bimat([5,2,3,4,1;4,3,1,2,5;4,2,3,5,1;3,5,2,4,1;3,2,5,4,1])
A=[[5,2,3,4,1;4,3,1,2,5;4,2,3,5,1;3,5,2,4,1;3,2,5,4,1]];
w=bimat(A,A')
w=bimat(A',A)
w=bimat(A,A')
A
payoff_matrix = A
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
syms x1 x2 x3 x4 x5   real
Payoff_vector_field_F = payoff_matrix *[x1 x2 x3 x4 x5]';
mean_U = [x1 x2 x3 x4 x5] * Payoff_vector_field_F;
V_Eq_0 = [x1 x2 x3 x4 x5]'.*(Payoff_vector_field_F - mean_U);
S=solve(V_Eq_0);
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5])
t20 =[]; for i=1:491; A=R((i-1)*5+1: (i-1)*5+1 ,:); w=bimat(A,A'); t20 =[t20;w];end
t20 =[]; for i=1:491; A=R((i-1)*5+1: (i-1)*5+5 ,:); w=bimat(A,A'); t20 =[t20;w];end
%-- 2021/7/18 21:49 --%
str5_exp_line
replace(tmp1,'M','j')
aa=['a' 'b' 'c']
aa={'a' 'b' 'c'}
strcat('v',aa(1))
cell2str(ans)
cellstr(ans)
bb = cellstr(ans)
strcat('v',bb)
dd = {ezp01205A551
ezp01206A552
ezp01206A553
ezp01229A554
ezp01230A555
ezp01230A556}
dd = {'ezp01205A551'
'ezp01206A552'
'ezp01206A553'
'ezp01229A554'
'ezp01230A555'
'ezp01230A556'}
f=strcat('f:/',dd(1))
f=strcat('D:\M\5_s_exp\',dd(1),'.csv')
[R meanR meanAll A]=f_Readxls5col(f)
[R meanR meanAll A]=f_Readxls5col(cellstr(f))
cellstr(f)
dd = ['ezp01205A551'
'ezp01206A552'
'ezp01206A553'
'ezp01229A554'
'ezp01230A555'
'ezp01230A556']
f=strcat('D:\M\5_s_exp\',dd(1),'.csv')
[R meanR meanAll A]=f_Readxls5col(cellstr(f))
[R meanR meanAll A]=f_Readxls5col(f)
f
dd(1)
dd(1,1:end)
f=strcat('D:\M\5_s_exp\',dd(1,1:end),'.csv')
clear
str5_exp_line
dd = ['ezp01205A551'
'ezp01206A552'
'ezp01206A553'
'ezp01229A554'
'ezp01230A555'
'ezp01230A556']; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据',dd(1,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1)
str5_exp_line_w
tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据\',dd(1,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1)
m6=[];for k=1:6; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据\',dd(k,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1); m6=[m5; meanAll];end
cleat;m6=[];for k=1:6; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据\',dd(k,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1); m6=[m6; meanAll];end
clear;m6=[];for k=1:6; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据\',dd(k,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1); m6=[m6; meanAll];end
dd = ['ezp01205A551'
'ezp01206A552'
'ezp01206A553'
'ezp01229A554'
'ezp01230A555'
'ezp01230A556']; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据',dd(1,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1)
dd = ['ezp01205A551'
'ezp01206A552'
'ezp01206A553'
'ezp01229A554'
'ezp01230A555'
'ezp01230A556'];
clear;m6=[];for k=1:6; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据\',dd(k,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1); m6=[m6; meanAll];end
m6=[];for k=1:6; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据\',dd(k,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1); m6=[m6; meanAll];end
dd = ['ezp01205A551'
'ezp01206A552'
'ezp01206A553'
'ezp01229A554'
'ezp01230A555'
'ezp01230A556'];
m6=[];for k=1:6; tmp1=strcat('F:\yqm代码+数据\5_s_exp5策略实验数据\',dd(k,1:end),'.csv'); meanAll=str5_exp_line_w(tmp1); m6=[m6; meanAll];end
boxplot(m6)
md'
m6'
x=[0.1486	0.1632	0.1258	0.176	0.1965	0.1573	0.1477
0.2136	0.2383	0.2532	0.2397	0.2212	0.2055	0.1892
0.2146	0.2102	0.1637	0.211	0.2057	0.227	0.2458
0.0964	0.0803	0.1377	0.1317	0.1138	0.1177	0.1123
0.3268	0.308	0.3197	0.2417	0.2628	0.2925	0.305
]
catter(x(:,1),x(:,2:6))>
catter(x(:,1),x(:,2:6))
scatter(x(:,1),x(:,2:6))
scatter(x(:,1),x(:,2))
for k=1:6; scatter(x(:,1),x(:,k));hold on;end
latex(x)
f = latex2MxWithMxPrecision(x, 3)
%-- 2021/7/20 13:04 --%
x1 =[ 2 3 1 4 5 ];
x0 = [2 2 2 2 2];
x = x1 - x0;
r = norm(x);
dqrt(15)
sqrt(15)
%-- 2021/7/21 0:30 --%
tmp202106129
eval(mean_U)
tmp202106129
eval(mean_U)
eval(Payoff_vector_field_F)
perms(1:5)
randi([1 120],[5 1])
tmp202107A21
csvwrite
save('manifold0718try2.mat', 'R','V','Ne', '-append')
save('manifold0718try2.mat', 'R','Ne' )
save('manifold0718try2.mat', 'R','V','Ne', '-append')
load('manifold0718try2.mat')
clear
load('manifold0718try2.mat')
clear
tmp202107A21
load('manifold0718try2.mat')
clear
load('manifold0718try2.mat')
tmp202107A21
clear
load('manifold0718try2.mat')
tmp202107A21
riV = [real(V) imag(V)]
plot(riV(:,7))
riV4 = round(riV,4);
riV4 = [round(riV,4) [1:size(riV,1)]'];
riV4 = [round(riV,4) [1:size(riV,1)]'  prod(riV')'];
riV4 = [round(riV,4) [1:size(riV,1)]'  prod(riV(:,1:5)')'];
riV4 = [round(riV,4) [1:size(riV,1)]'  prod(riV4(:,1:5)')'];
clear
load('manifold0718try2.mat')
riV4 = round(riV,4);
riV4 = [round(riV,4) [1:size(riV,1)]'  prod(riV4(:,1:5)')'];
riV = [real(V) imag(V)]
riV4 = round(riV,4);
riV4 = [round(riV,4) [1:size(riV,1)]'  prod(riV4(:,1:5)')'];
centermanifoldpara=riV4(find(riV4(:,13) ==0),:);
centermanifoldpara=riV4(find(riV4(:,12) ==0),:);
retM5x5 = get_payoffmatrix_from_index(R,438)
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(payoff_matrix)
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
A=[12	1.97E-07	5.59E-05	-5.1E-06	4.502E-05	4.514E-05
13	-1.8E-07	-4.5E-05	-8.8E-06	-3.5E-05	-3.6E-05
14	-1.5E-07	-4.3E-05	1.27E-05	-3.55E-05	-3.41E-05
15	1.37E-07	3.23E-05	1.2E-06	2.551E-05	2.5E-05
23	2.19E-07	5.76E-05	-1.9E-06	4.417E-05	4.435E-05
24	-4.6E-08	5.12E-06	2.66E-06	2.191E-06	2.452E-06
25	2.37E-08	-6.8E-06	-5.9E-06	-1.34E-06	-1.66E-06
34	2.01E-07	3.99E-05	-1.7E-05	3.289E-05	3.128E-05
35	-1.6E-07	-2.8E-05	6.71E-06	-2.37E-05	-2.29E-05
45	6.53E-10	2.07E-06	-2E-06	-4.27E-07	-4.17E-07
]
B=[A(:,1) A(:,2:5)*10^6]
f = latex2MxWithMxPrecision(B, 3)
B=[A(:,1) A(:,2:end)*10^6]
f = latex2MxWithMxPrecision(B, 3)
%-- 2021/7/21 20:09 --%
load('manifold0718try2.mat')
retM5x5 = get_payoffmatrix_from_index(R,438)
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
retM5x5=[
1 2 3 4 5
5 1 2 3 4
4 5 1 2 3
3 4 5 1 2
2 3 4 5 1]
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
plot(eigen_vector(:,2))
scatter(eigen_vector(:,2))
angle(eigen_vector(:,2))
x = angle(eigen_vector(:,2))
x(3)-x(4)
(x(3)-x(4))*5/pi
tM = retM5x5;
tM(3,:) = 3
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(tM)
x = angle(eigen_vector(:,2))
plot(eigen_vector(:,2))
%-- 2021/7/22 16:29 --%
retM5x5=[
1 2 3 4 5
5 1 2 3 4
4 5 1 2 3
3 4 5 1 2
2 3 4 5 1]
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
retM5x5=[
5 1 2 3 4
1 2 3 4 5
4 5 1 2 3
3 4 5 1 2
2 3 4 5 1]
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
retM5x5=[
5 1 2 3 4
1 2 3 4 5
4 5 1 2 3
3 4 5 1 2
2 3 4 5 1]
retM5x5=[
5 1 2 3 4
4 5 1 2 3
3 4 5 1 2
2 3 4 5 1
1 2 3 4 5]
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
retM5x5=[
2 1 2 3 4
4 2 1 2 3
3 4 2 1 2
2 3 4 2 1
1 2 3 4 2]
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
bimat(retM5x5,retM5x5')
retM5x5=[
2 1 2 3 4
4 2 1 2 3
3 4 2 1 2
2 3 4 2 1
1 2 3 4 2.0001]
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
retM5x5=[
2 1 2 3 4
4 2 1 2 3
3 4 2 1 2
2 3 4 2 1
1 2 3 4 2]
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(retM5x5)
syms a; rM5x5=[
a 1 2 3 4
4 a 1 2 3
3 4 a 1 2
2 3 4 a 1
1 2 3 4 a]; eva = []; for k=1:5; a = k; [eigen_vector, d_eigen_value, rest_point_no]
= get_eigensystem_from_5x5(eval(rM5x5)); eva = [eva d_eigen_value]; end
syms a; rM5x5=[
a 1 2 3 4
4 a 1 2 3
3 4 a 1 2
2 3 4 a 1
1 2 3 4 a]; eva = []; for k=1:5; a = k; [eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(eval(rM5x5)); eva = [eva d_eigen_value]; end
syms a; rM5x5=[
a 1 2 3 4
4 a 1 2 3
3 4 a 1 2
2 3 4 a 1
1 2 3 4 a]; eva = []; for k=1:0.5:5; a = k; [eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(eval(rM5x5)); eva = [eva d_eigen_value]; end
compass(eva(:,2))
syms a; rM5x5=[
a 1 2 3 4
4 a 1 2 3
3 4 a 1 2
2 3 4 a 1
1 2 3 4 a]; eva = []; for k=2:0.1:3; a = k; [eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(eval(rM5x5)); eva = [eva d_eigen_value]; end
r_eva = real(eva)
r_eva = imag(eva)
r_eva(2,:)./r_eva(4,:)
syms a; rM5x5=[
a 1 2 3 4
4 a 1 2 3
3 4 a 1 2
2 3 4 a 1
1 2 3 4 a]; eva = [];eve = []; for k=2:0.5:3; a = k; [eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(eval(rM5x5)); eva = [eva d_eigen_value];eve = [eve; eigen_vector]; end
abseve=abs(eve);
angeve=angle(eve);
rM5x5
rM5x5 - 2
rM5x5=[
a 1 2 3 4+b
4 a 1 2+b 3
3 4 a+b 1 2
2 3+b 4 a 1
1+b 2 3 4 a]
syms a b;
rM5ab=[
a 1 2 3 4+b
4 a 1 2+b 3
3 4 a+b 1 2
2 3+b 4 a 1
1+b 2 3 4 a]
rM5ab=[
a 1 2 3 4+b
4 a 1 2+b 3
3 4 a+b 1 2
2 3+b 4 a 1
1+b 2 3 4 a]; eva = [];eve = []; for k=2:0.5:3; for b=-1:1; a = k; [eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(eval(rM5ab)); eva = [eva d_eigen_value];eve = [eve; eigen_vector]; end; end
Parameter20210722
eve4=round(eve,4);
eva4=round(eva,4);
plot(eve4(21:25,2))
plot(eve4(6:10,2))
plot(eve4(11:15,4))
plot(eve4(1:5,4))
figure;
plot(eve4(6:10,4))
figure;
plot(eve4(6:10,2))
figure;
plot(eve4(11:15,4))
vv=[-0.1022 + 0.154i
0.5417 + 0.0952i
0.4370 - 0.0952i
-0.2716 - 0.154i
-0.6048 + 0.000i
]
plot(vv)
figure;
vw=[0.6048 + 0.000i
0.2716 - 0.154i
-0.4370 - 0.0952i
-0.5417 + 0.0952i
0.1022 + 0.154i
]
plot(vw)
abs(vv)
abs(vw)
angle(vw)
angle(vv)
from_eigenvector_out_am(vv)
from_eigenvector_out_am([vv vw])
from_eigenvector_out_am(vw)
v0=[0.4472 + 0.000i
0.1382 + 0.4253i
-0.3618 + 0.2629i
-0.3618 - 0.2629i
0.1382 - 0.4253i
]
from_eigenvector_out_am(v0)
%-- 2021/7/23 0:01 --%
syms a b;
rM5ab=[ a 1 2 3 4+b
4 a 1 2+b 3
3 4 a+b 1 2
2 3+b 4 a 1
1+b 2 3 4 a];
eva = [];eve = [];
for k=2:0.5:3;
for b=-3:3:3;
a = k;
eval_M5ab = eval(rM5ab)
[eigen_vector, d_eigen_value, rest_point_no] = ...
get_eigensystem_from_5x5(eval_M5ab);
eva = [eva d_eigen_value];
eve = [eve; eigen_vector];
end;
end
syms a b;
rM5ab=[ a 1 2 3 4+b
4 a 1 2+b 3
3 4 a+b 1 2
2 3+b 4 a 1
1+b 2 3 4 a];
eva = [];eve = []; Mabed = [];
for k=2:0.5:3;
for b=-3:3:3;
a = k;
eval_M5ab = eval(rM5ab)
Mabed = [Mabed; 98 eval_M5ab(1,:) eval_M5ab(2,:) eval_M5ab(3,:) eval_M5ab(4,:) eval_M5ab(5,:) 99 ];
end;
end
syms a b;
rM5ab=[ a 1 2 3 4+b
4 a 1 2+b 3
3 4 a+b 1 2
2 3+b 4 a 1
1+b 2 3 4 a];
eva = [];eve = []; Mabed = [];
for k=2:0.5:3;
for b=-3:3:3;
a = k;
eval_M5ab = eval(rM5ab)
Mabed = [Mabed; 98 eval_M5ab(1,:) 99 98 eval_M5ab(2,:)  99 98 eval_M5ab(3,:)  99 98 eval_M5ab(4,:)  99 98 eval_M5ab(5,:) 99 ];
end;
end
syms a b;
rM5ab=[ a 1 2 3 4+b
4 a 1 2+b 3
3 4 a+b 1 2
2 3+b 4 a 1
1+b 2 3 4 a];
eva = [];eve = []; Mabed = [];
for k=2:0.5:3;
for b=-3:3:3;
a = k;
eval_M5ab = eval(rM5ab)
Mabed = [Mabed; 98 98 eval_M5ab(1,:) 99 98 eval_M5ab(2,:)  99 98 eval_M5ab(3,:)  99 98 eval_M5ab(4,:)  99 98 eval_M5ab(5,:) 99 99 ];
end;
end
%-- 2021/7/23 15:46 --%
fileFolder=fullfile('F:\tmpdata\result_exp202107');
dirOutput=dir(fullfile(fileFolder,'*.csv'))
fileNames={dirOutput.name}
fileNames(2)
cellstr(fileNames(2))
cellstr(fileNames)
fileNames={dirOutput.name}'
fileNames{1}
x = fileNames{1}
x = fileNames{1}; y= strrep(strrep(strrep(x,'[',''),']',''),'.csv','')
str2num(y)
z=str2num(y)
rNatrix=[]; for k=0:4; rNatrix=[rNatrix;z(1+k*5:5+k*5) ]; end
clear
file_path='F:\tmpdata\result_exp202107';
fileFolder=fullfile();file_path
dirOutput=dir(fullfile(fileFolder,'*.csv'))
fileNames={dirOutput.name}'
file_path='F:\tmpdata\result_exp202107';
fileFolder=fullfile(file_path);
dirOutput=dir(fullfile(fileFolder,'*.csv'))
fileNames={dirOutput.name}'
size(fileNames)
size(fileNames,1)
for fj = 1:size(fileNames,1);
x = fileNames{fj}; y= strrep(strrep(strrep(x,'[',''),']',''),'.csv','');z=str2num(y)
rNatrix=[]; for k=0:4; rNatrix=[rNatrix;z(1+k*5:5+k*5) ]; end;
end
clear
%
file_path='F:\tmpdata\result_exp202107';
fileFolder=fullfile(file_path);
dirOutput=dir(fullfile(fileFolder,'*.csv'));
fileNames={dirOutput.name}';
M=zeros(5,5,:)
for fj = 1:size(fileNames,1);
x = fileNames{fj};
y= strrep(strrep(strrep(x,'[',''),']',''),'.csv','');
z=str2num(y)
rMatrix=[];
for k=0:4;
rMatrix=[rMatrix;z(1+k*5:5+k*5) ];
end;
M=zeros(:,:,fj)=rMatrix;
end
%
file_path='F:\tmpdata\result_exp202107';
fileFolder=fullfile(file_path);
dirOutput=dir(fullfile(fileFolder,'*.csv'));
fileNames={dirOutput.name}';
M=zeros(5,5,:)
for fj = 1:size(fileNames,1);
x = fileNames{fj};
y= strrep(strrep(strrep(x,'[',''),']',''),'.csv','');
z=str2num(y)
rMatrix=[];
for k=0:4;
rMatrix=[rMatrix;z(1+k*5:5+k*5) ];
end;
M(:,:,fj)=rMatrix;
end
zeros(5)
%
file_path='F:\tmpdata\result_exp202107';
fileFolder=fullfile(file_path);
dirOutput=dir(fullfile(fileFolder,'*.csv'));
fileNames={dirOutput.name}';
% M=zeros(5,5,:)
for fj = 1:size(fileNames,1);
x = fileNames{fj};
y= strrep(strrep(strrep(x,'[',''),']',''),'.csv','');
z=str2num(y)
rMatrix=[];
for k=0:4;
rMatrix=[rMatrix;z(1+k*5:5+k*5) ];
end;
M(:,:,fj)=rMatrix;
end
file_path='F:\tmpdata\result_exp202107';[M5x5, fileNames]=get_5x5_from_abedExp(file_path)
%-- 2021/7/24 9:29 --%
gdt_integrate_1
dos_x_t = get_dos_x_t_from_abeddata(abedfilesname)
[Yret3,mn,am_eigencycleDim_t, mean_x] = from_N_colExp_out_am(dos_x_t)
gdt_integrate_1
rr_dos_x_t(:,:,fj) = dos_x_t
gdt_integrate_1
rr_abed_csv(:,:,fj) = abed_csv;
gdt_integrate_1
save('manifold20210724try.mat', 'rr_abed_csv', ...
'rr_dos_x_t', ...
'rr_E_mean_std', ...
'rr_E_angularmom', ...
'rr_eigen_vector', ...
'rr_d_eigen_value', ...
'rr_T_eigencycle', ...
'rr_rest_no', ...
'rr_solu_list', ...
'rr_coef_coefint3col', ...
'rr_resid_residint3col', ...
'rr_R2_F_p_evar');
clear
load('manifold20210724try.mat')
cellstr(rr_abed_csv)
cellstr(rr_abed_csv.value)
{rr_abed_csv}
rr_abed_csv(:,:,1)
cellstr(rr_abed_csv(:,:,1))
cell2mat(rr_solu_list(:,:,3))
rr_abed_csv{:,:,1}
load('C:\Users\Think\AppData\Local\Temp\Rar$DRa0.217\result_ana20210724\manifold20210724try.mat')
cell2mat(rr_solu_list(:,:,3))
clear
Parameter20210722
gen_abed_payoff
%-- 2021/7/24 19:00 --%
Parameter20210722
clear
gdt_integrate_1('F:\tmpdata\result_exp202107\','manifold20210724t9.mat')
load('manifold20210724t9.mat')
plot(rr_dos_x_t(:,:,5))
figure; plot(rr_dos_x_t(:,:,2));
figure; plot(rr_dos_x_t(:,:,8));
rr_E_angularmom(:,:,2:3:8)
[rr_E_angularmom(:,:,2) rr_E_angularmom(:,:,5) rr_E_angularmom(:,:,8) ]
plotmatrix(ans)
zz=[rr_T_eigencycle(:,:,2) rr_T_eigencycle(:,:,5) rr_T_eigencycle(:,:,8) ]
zval=[rr_d_eigen_value(:,:,2) rr_d_eigen_value(:,:,5) rr_d_eigen_value(:,:,8) ]
figure; plot(rr_dos_x_t(:,:,7));
cell2mat(rr_solu_list(:,:,7))
%-- 2021/7/24 23:54 --%
load('manifold0718try2 - 0721.mat')
clear
load('manifold5try.mat')
tmp202107A24
%-- 2021/7/25 9:06 --%
tmp202107A24
%-- 2021/7/25 14:40 --%
tmp202107A24
%gdt_integrate_1('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200.mat')
gdt_integrate_1('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200.mat')
gdt_integrate_1
gdt_integrate_1('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200.mat')
zeros(5)
gdt_integrate_1('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200.mat')
%-- 2021/7/25 19:27 --%
gdt_integrate_1('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200.mat')
%-- 2021/7/25 19:53 --%
gdt_integrate_onlyMean
gdt_integrate_onlyMean('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200.mat')
gdt_integrate_onlyMean('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200_onlyMean.mat')
gdt_integrate_onlyMean
gdt_integrate_onlyMean('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725t001_200_onlyMean.mat')
load('F:\tmpdata\result_ana202107\InvMani0725t001_200_onlyMean.mat')
%-- 2021/7/25 20:54 --%
load('F:\tmpdata\result_ana202107\InvMani0725t001_200_onlyMean.mat')
a=[];for k=1:200; a=[a;rr_E_mean_std(:,1,k)'];end
a(:,6:7) =[min(a')' max(a')'];
plot(a(:,6))
plot(a(:,end))
%-- 2021/7/25 21:58 --%
load('F:\tmpdata\result_ana202107\InvMani0725t001_200_onlyMean.mat')
a=[];for k=1:200; a=[a;rr_E_mean_std(:,1,k)'];end
a(:,6:7) =[min(a')' max(a')'];
find(a(:,6) > 0.15
find(a(:,6) > 0.15)
find(a(:,6) > 0.18)
xx=find(a(:,6) > 0.18)
for kk=1:size(xx)
fj=xx(kk)
E_angularmom = rr_E_angularmom(:,:,fj)
% theoretical result
[eigen_vector, d_eigen_value, T_eigencycle, n0, solu_list] = ...
get_eigensystem_from_5x5(rr_M5x5(:,:,fj));
% statics analysis
[coef_coefint3col,resid_residint3col,R2_F_p_evar] = ...
multi_reg_Eam_Tec(E_angularmom,T_eigencycle);
coef_coefint3col'
R2_F_p_evar
end
xx=find(a(:,6) > 0.18)
for kk=1:size(xx)
fj=xx(kk)
E_angularmom = rr_E_angularmom(:,:,fj);
% theoretical result
[eigen_vector, d_eigen_value, T_eigencycle, n0, solu_list] = ...
get_eigensystem_from_5x5(rr_M5x5(:,:,fj));
% statics analysis
[coef_coefint3col,resid_residint3col,R2_F_p_evar] = ...
multi_reg_Eam_Tec(E_angularmom,T_eigencycle);
bbb=coef_coefint3col'
R2_F_p_evar
end
clear
tmp202107A25c
gdt_integrate_onlyMean('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725_onlyMean.mat','F:\tmpdata\result_ana202107\InvMani0725_dosXt.mat')
%-- 2021/7/26 1:16 --%
gdt_integrate_onlyMean('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725_onlyMean.mat','F:\tmpdata\result_ana202107\InvMani0725_dosXt.mat')
gdt_integrate_onlyMean('F:\tmpdata\result_exp202107\','F:\tmpdata\result_ana202107\InvMani0725_onlyMean.mat','F:\tmpdata\result_ana202107\InvMani0725_dosXt.mat')
save(file_save_mat, 'rr_abed_csv', ...
'rr_E_mean_std', ...
'rr_E_angularmom', ... %%%
'rr_M5x5');
save(file_dos_xt_mat, 'rr_dos_x_t', 'rr_M5x5')
%-- 2021/7/26 11:11 --%
load('F:\tmpdata\result_ana202107\InvMani0725_onlyMean.mat')
mean_std10 = []; for k=1:489; mean_std10 = [mean_std10; rr_E_mean_std(:,1,k)']; end
mean_std10(:,6:7) =[min(mean_std10')' max(mean_std10')'];
bbb=[]; xx=find(mean_std10(:,6) > 0.18)
for kk=1:size(xx)
fj=xx(kk)
E_angularmom = rr_E_angularmom(:,:,fj);
% theoretical result
[eigen_vector, d_eigen_value, T_eigencycle, n0, solu_list] = ...
get_eigensystem_from_5x5(rr_M5x5(:,:,fj));
% statics analysis
[coef_coefint3col,resid_residint3col,R2_F_p_evar] = ...
multi_reg_Eam_Tec(E_angularmom,T_eigencycle);
bbb=[bbb; mean_std10(fj,:) 89 coef_coefint3col(:,1)' 99 R2_F_p_evar]
end
tmp202107A26a
rr_M5x5(:,:,fj)
tmp202107A26a
save('F:\tmpdata\result_ana202107\InvMani0725_onlyMeanbbb.mat','bbb' )
save('F:\tmpdata\result_ana202107\InvMani0725_onlyMeanbbb.mat','bbb','singularity' )
tmp202107A26a
norm(E_angularmom)
tmp202107A26a
m=[9
37
48
56
61
70
73
77
98
101
116
117
120
124
126
129
130
149
194
195
198
290
324
362
423
];
s=rr_M5x5(:,:,m)
save('F:\tmpdata\result_ana202107\InvMani0725_onlyMeanbbb18.mat','bbb','singularity' )
tmp202107A26a
plot(bbb(:,32))
bar(bbb(:,32))
loglog(bbb(:,[18,32]),'DisplayName','bbb(:,[18,32])')
scatter(getcolumn(bbb(:,[18,32]),1),getcolumn(bbb(:,[18,32]),2))
loglog(bbb(:,[18,32]),'DisplayName','bbb(:,[18,32])')
zz=log(bbb(:,[18,32]))
scatter(zz(:,1),zz(:,2))
xxx = bbb(:,22:25)
xx1 = xxx+abs(xxx);
xx2=xx1/2; xx3=sqrt(xx2(:,1).*xx2(:,2))+ sqrt(xx2(:,3).*xx2(:,4))
scatter(xx3, bbb(:,32))
find(xxx == xx3)
in=[];for j=1:size(xx3); in=[in; find(aba(xxx(j,:)-xx3(j))<0.00001)];end
in=[];for j=1:size(xx3); in=[in; find(abs(xxx(j,:)-xx3(j))<0.00001)];end
in=[];for j=1:size(xx3); tmp=find(abs(xxx(j,:)-xx3(j))<0.00001); if size()>0; tmp1=tmp;else;tmp1=0;end;in=[in; tmp1]; end
in=[];for j=1:size(xx3); tmp=find(abs(xxx(j,:)-xx3(j))<0.00001); if size(tmp)>0; tmp1=tmp;else;tmp1=0;end;in=[in; tmp1]; end
in=[];for j=1:size(xx3); tmp=find(abs(xxx(j,:)-xx3(j))<0.00001); if size(tmp)>0; tmp1=tmp(1);else;tmp1=0;end;in=[in; tmp1]; end
yyy = bbb(:,27:30);
yy1=zeros(size(yyy,1),1);
for j=1:size(yyy,1); if in(j,1)>0; yy1(j)=yyy(j,in(j,1));end;end
rrr=[20.4762979	0.120858471	0.167782073
0.175362294	0	0
80.50659933	0.1	0.387298335
39.07223058	0.05720757	0.465050311
7.684059297	0.033228298	0.394447893
127.0293276	0.179315485	0.494028359
0.144335572	0	0
62.38837923	0.060090325	0.558178024
0.183986507	0	0
0.071309745	0	0
0.300427754	0	0
0.043562268	0	0
111.8268084	0.197073533	0.415100498
0.375447576	0	0
34.90245173	0.048287601	0.380807831
0.121712013	0	0
0.206471011	0	0
5.850993969	0.015091108	0.30455888
0.108554836	0	0
0.109461501	0	0
0.069682137	0	0
1.067655942	0	0
0.108072778	0	0
0	0	0
75.7331202	0.155407809	0.496257623
0.060342041	0	0
0.17556643	0	0
0.14164861	0	0
47.9026081	0.105919744	0.413263975
104.4076565	0.20306276	0.39298761
233.7232364	0.315235163	0.814228175
10.56535191	0.054241393	0.188290444
178.0089282	0.206551196	0.673174888
33.64832357	0.067359306	0.416667527
0.215473908	0	0
178.9153096	0.300720395	0.596426562
0.224845176	0	0
0.226513217	0	0
]
y = rrr(:,1);X = [ones(size(y)) ,...
rrr(:,2),...
rrr(:,3) ];
[b,bint,r,rint,stats] = regress(y,X)
stats(3)
y = rrr(:,1);X = [ones(size(y)) ,...
rrr(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
stats(3)
y = rrr(:,1);X = [ones(size(y)) ,...
rrr(:,3)  ];
[b,bint,r,rint,stats] = regress(y,X)
stats(3)
scatter(y, rrr(:,3))
scatter(y, rrr(:,2))
rrr=[20.4762979	0.120858471	0.167782073
80.50659933	0.1	0.387298335
39.07223058	0.05720757	0.465050311
7.684059297	0.033228298	0.394447893
127.0293276	0.179315485	0.494028359
62.38837923	0.060090325	0.558178024
111.8268084	0.197073533	0.415100498
34.90245173	0.048287601	0.380807831
5.850993969	0.015091108	0.30455888
75.7331202	0.155407809	0.496257623
47.9026081	0.105919744	0.413263975
104.4076565	0.20306276	0.39298761
233.7232364	0.315235163	0.814228175
10.56535191	0.054241393	0.188290444
178.0089282	0.206551196	0.673174888
33.64832357	0.067359306	0.416667527
178.9153096	0.300720395	0.596426562
]
y = rrr(:,1);X = [ones(size(y)) ,...
rrr(:,2),...
rrr(:,3) ];
[b,bint,r,rint,stats] = regress(y,X)
stats(3)
y = rrr(:,1);X = [ones(size(y)) ,...
rrr(:,2)  ];
[b,bint,r,rint,stats] = regress(y,X)
stats(3)
y = rrr(:,1);X = [ones(size(y)) ,...
rrr(:,3)  ];
[b,bint,r,rint,stats] = regress(y,X)
stats(3)
stats(1)
clear
Parameter50
clear
Parameter50
%-- 2021/7/26 21:19 --%
syms a b;
rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
eva = [];eve = []; Mabed=[];
gen_abed_payoff = [];
for k=2:0.5:3; for b=-1:1:1;
a = k;
eval_M5ab = eval(rM5ab)
[eigen_vector, d_eigen_value, rest_point_no] = ...
get_eigensystem_from_5x5(eval_M5ab);
eva = [eva d_eigen_value];
eve = [eve; eigen_vector];
str1 = [98 98 eval_M5ab(1,:) 99 98 eval_M5ab(2,:)  99 98 eval_M5ab(3,:)  99 98 eval_M5ab(4,:)  99 98 eval_M5ab(5,:) 99 99];
Mabed = [Mabed; str1];
str2 = strcat(' "[[', num2str(eval_M5ab(1,:)), '] [', num2str(eval_M5ab(2,:)), '] [', num2str(eval_M5ab(3,:)), '] [',num2str(eval_M5ab(4,:)), '] [',num2str(eval_M5ab(5,:)), ']]" ');
gen_abed_payoff = [gen_abed_payoff str2];
end;end
%-- 2021/7/27 2:43 --%
syms a b;
rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
get_eigensystem_from_5x5
b=1;a=2.3;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
b=1;a=2.2;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
b=1;a=2.16;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
b=1;a=2.18;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
b=1;a=2.19;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
b=1;a=2.195;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
b=1;a=2.191;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
b=1;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
clear
b=1;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
syms a; b=1;rM5ab=[ a 1+b 2 3 4
4 a 1+b 2 3
3 4 a 1+b 2
2 3 4 a 1+b
1+b 2 3 4 a];
[eigen_vector, d_eigen_value, rest_point_no] = get_eigensystem_from_5x5(rM5ab)
x1=0.2;x2=0.2;x3=0.2;x4=0.2;x5=0.2;
D_Eq_at_NE = eval(Jac)
[eigen_vector, eigen_value] = eig(Jac)
Jac
eig([1 a;-a,1])
syms a
eig([1 a;-a,1])
Parameter5ab
D_Eq_at_NE =  subs(Jac,x1,0.2)
D_Eq_at_NE =  subs(Jac,'x1',0.2)
D_Eq_at_NE =  subs(subs(subs(subs(subs(Jac,'x1',0.2),'x2',0.2),'x3',0.2),'x4',0.2),'x5',0.2);
[eigen_vector, eigen_value] = eig(D_Eq_at_NE)
syms a
eig([1 a;-a,1])
clear
Parameter5ab
payoff_matrix
S=solve(V_Eq_0)
A=eval([S.x1 S.x2 S.x3 S.x4 S.x5])
Jac
D_Eq_at_NE
Parameter5ab
[eigen_vector, eigen_value] = eig(D_Eq_at_NE)
vap(eigen_value)
vpa(eigen_value)
format short
help vpa
xx=vpa(eigen_value(:,:),4)
xx=vpa(eigen_value(:,:),2)