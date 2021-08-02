% function r=payoff88()
%   r = fpayoff88(2,1) %6
%       payoff88(3,2) %8
%       payoff88(4,2) %10 
% end


function P=payoff88(n,m)

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

for i=1:8%for Ann %À„Ananµƒ ’“Ê
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

