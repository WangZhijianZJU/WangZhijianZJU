% D:\MATLAB\R2016a\bin\PlatomicSolid.m
function r = PlatomicSolid()   %郭经纬取数据用 12 24 2019
%    PlatomicSolid_Original  姚芹梅出图
   r=1;
%    for K = 1:5
   for K = 6:6
     [V,F,E]=platonic_solid(K,r);  %郭经纬取数据用
%       subplot(2,3,K);
     patch('Faces',F,'Vertices',V,'FaceColor','g','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;
     BlueForceArr = [1 3 5 7];
        for j=1:length(V(:,1))
            props.Marker='o';
    %         props.MarkerEdgeColor='r';
    %         props.MarkerFaceColor='g';
            if find(BlueForceArr==j) > 0; props.MarkerFaceColor='b'; else; props.MarkerFaceColor='r'; end
            props.LineStyle='none';
            drawPoint3d(gca,V(j,:),props,'MarkerSize',10)
            text(V(j,1),V(j,2),V(j,3),num2str(j),'fontsize',20)
        end
        for k=1:length(E(:,1))
            A=V(E(k,1),:);B=V(E(k,2),:);
            plot3([A(:,1),B(:,1)],[A(:,2),B(:,2)],[A(:,3),B(:,3)],'k');
        end
   end
end


function r = PlatomicSolid_Original()
clear all; close all; clc;

r=1;
figure;fig=gcf; clf(fig); colordef (fig, 'white'); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units); set(fig,'Color',[1 1 1]);
hold on; 

[V,F]=platonic_solid(1,r);
subplot(2,3,1);
patch('Faces',F,'Vertices',V,'FaceColor','w','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;

[V,F]=platonic_solid(2,r);
subplot(2,3,2);
patch('Faces',F,'Vertices',V,'FaceColor','w','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;

[V,F]=platonic_solid(3,r);
subplot(2,3,3);
patch('Faces',F,'Vertices',V,'FaceColor','w','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;

[V,F]=platonic_solid(4,r);
subplot(2,3,4);
patch('Faces',F,'Vertices',V,'FaceColor','w','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;

[V,F]=platonic_solid(5,r);
subplot(2,3,5);
patch('Faces',F,'Vertices',V,'FaceColor','w','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;

end



function [V,F,E]=platonic_solid(n,r)
% function [V,F]=platonic_solid(n,r)
% ------------------------------------------------------------------------
% PLATONIC_SOLID Creates the PATCH data, the vertices (V) and faces (F) for
% a given platonic solid (according to "n" see below) with radius (r)
%
% n=1 -> Tetrahedron
% n=2 -> Cube
% n=3 -> Octahedron
% n=4 -> Icosahedron
% n=5 -> Dodecahedron
%
% %%%Example
% clear all; close all; clc;
% 
% r=1;
% figure;fig=gcf; clf(fig); colordef (fig, 'white'); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units); set(fig,'Color',[1 1 1]);
% hold on; 
% 
% [V,F]=platonic_solid(1,r);
% subplot(2,3,1);
% patch('Faces',F,'Vertices',V,'FaceColor','b','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;
% 
% [V,F]=platonic_solid(2,r);
% subplot(2,3,2);
% patch('Faces',F,'Vertices',V,'FaceColor','b','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;
% 
% [V,F]=platonic_solid(3,r);
% subplot(2,3,3);
% patch('Faces',F,'Vertices',V,'FaceColor','b','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;
% 
% [V,F]=platonic_solid(4,r);
% subplot(2,3,4);
% patch('Faces',F,'Vertices',V,'FaceColor','b','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;
% 
% [V,F]=platonic_solid(5,r);
% subplot(2,3,5);
% patch('Faces',F,'Vertices',V,'FaceColor','b','FaceAlpha',0.6,'EdgeColor','k','LineWidth',2); axis equal; grid on; hold on; view(3); axis off;
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 13/11/2009
% ------------------------------------------------------------------------
phi=(1+sqrt(5))/2;
switch n
    case 1 % Tetrahedron
        V1=[-0.5;0.5;0;0;];
        V2=[-sqrt(3)/6;  -sqrt(3)/6; sqrt(3)/3; 0];
        V3=[-0.25.*sqrt(2/3); -0.25.*sqrt(2/3); -0.25.*sqrt(2/3);  0.75.*sqrt(2/3)];
        F= [1 2 3; 1 2 4; 2 3 4; 1 3 4;];
    case 2 % Cube
        V1=[-1;  1; 1; -1; -1;  1; 1; -1;];
        V2=[-1; -1; 1;  1; -1; -1; 1;  1;];
        V3=[-1; -1;-1; -1;  1;  1; 1;  1;];
        F= [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8;];
    case 3 % Octahedron
        V1=[-1;  1; 1; -1;  0;   0;];
        V2=[-1; -1; 1;  1;  0;   0;];
        V3=[ 0;   0; 0;  0; -1;  1;];
        F= [1 2 5; 2 3 5; 3 4 5; 4 1 5; 1 2 6; 2 3 6; 3 4 6; 4 1 6;];
    case 4 % Icosahedron
        V1=[0;0;0;0;-1;-1;1;1;-phi;phi;phi;-phi;];
        V2=[-1;-1;1;1;-phi;phi;phi;-phi;0;0;0;0;];
        V3=[-phi;phi;phi;-phi;0;0;0;0;-1;-1;1;1;];
        F= [1 4 9;1 5 9;1 8 5;1 8 10;1 10 4; 12 2 5; 12 2 3; 12 3 6; 12 6 9; 12 9 5; 11 7 10; 11 10 8; 11 8 2; 11 2 3; 11 3 7; 2 5 8; 10 4 7; 3 6 7; 6 7 4; 6 4 9; ];
    case 5 % Dodecahedron
        V1=[1;(1/phi);-phi;phi;-1;0;-phi;1;-1;-1;1;(1/phi);-1;0;0;-(1/phi);phi;-(1/phi);1;0;];
        V2=[1;0;-(1/phi);(1/phi);1;-phi;(1/phi);-1;1;-1;-1;0;-1;-phi;phi;0;-(1/phi);0;1;phi;];
        V3=[[1;phi;0;0;-1;-(1/phi);0;1;1;1;-1;-phi;-1;(1/phi);-(1/phi);phi;0;-phi;-1;(1/phi);]];
        F=[1,2,16,9,20;2,16,10,14,8;16,9,7,3,10;7 9 20 15 5;5,7,3,13,18;3,13,6,14,10;6,13,18,12,11;6,11,17,8,14;11,12,19,4,17;1,2,8,17,4;1,4,19,15,20;12,18,5,15,19];
    case 6 % C60
        [V12,F12,E12]=platonic_solid(4,r)
        V60=[]; component=[];
        for L = 1:length(E12(:,1))
            delta = (V12(E12(L,2),:) - V12(E12(L,1),:) )% / ((sqrt(5)+1)/2)
            V60 = [V60; V12(E12(L,1),:)+delta/3; V12(E12(L,1),:)+2*delta/3 ]
            component=[component; E12(L,1) E12(L,2) ; E12(L,1) E12(L,2) ];
        end
        V1 = V60(:,1);V2 = V60(:,2);V3 = V60(:,3);
%         F=[1 7 5 3 9;19 25 21 23 27; 2 15 11 13 17; 35 29 31 33 20; 52 60 56 10 36;16 26 50 57 54;22 37 30 39 41;14 47 6 49 51; ...
%             8 34 46 59 42; 4 38 12 45 43; 48 32 53 55 24; 40 58 44 28 18];
        F=[1 7 3 5 9; 12 27 21 23 25;11 17 15 13 19;2 33 29 31 35; 4 37 14 41 39;30 43 22 47 45; ...
             44 24 51 49 32;16 38 6 53 55; 18 26 52 60 56; 57 46 34 8 40; 59 50 36 10 54;28 20 42 58 48];
    otherwise
        warning('False input for n')
end
    [THETA,PHI,R]=cart2sph(V1,V2,V3);
    R=r.*ones(size(V1(:,1)));
%     [V1,V2,V3]=sph2cart(THETA,PHI,R);
    V=[V1 V2 V3];
    [E V]= getEdges(V);
end


function [rE v1]= getEdges(V)
    if length(V(:,1))==12;theta=atan(1/((sqrt(5)+1)/2));else; theta=0;end
    Rot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    v1=[];
        for k=1:length(V(:,1))
            vtmp = Rot*V(k,:)';
            v1=[v1; vtmp'];
        end
        v=v1;
%     v=V;
    L=length(V(:,1))
    D=100*ones(L-1,1);
    for I=2:L
        dv = (v(1,:) - v(I,:));
        D(I) = sqrt(dv * dv')
    end
    d=min(D);
E = [];
    for m = 1:length(v(:,1)) 
        for n = m+1:length(v(:,1))
                du = v(m,:) - v(n,:);
                dmn = sqrt(du * du');
                if dmn < d*1.01 && dmn > d*0.99;
                    E=[E; m n];
                end
        end
    end
rE = E;

end