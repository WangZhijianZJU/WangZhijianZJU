
function [n8x8list,stateList,RPara] = Get8x8Data(Para)
n8x8list = load('C:\Users\Think\Downloads\n8x8.csv');
n8x8list(:,1) = n8x8list(:,1) * 100 - 8000000;

%构造 State_t0,State_tp1(+1,后一轮),State_tm1(-1,前一轮)
n8x8list(:,10) = n8x8list(:,6)*10 + n8x8list(:,7);
stateList=round(unique(n8x8list(:,10)));

n8x8list(:,11:12) = 999;
for L = 1:72000-1; n8x8list(L,11) = n8x8list(L+1,10);end
for L = 1+1:72000; n8x8list(L,12) = n8x8list(L-1,10);end

 
%分参数，每个参数有24000记录
RPara = n8x8list(find(n8x8list(:,9) == Para & n8x8list(:,5) == 1),:); 
% ------------------------------------------------- ^^ Ann data enough
end 