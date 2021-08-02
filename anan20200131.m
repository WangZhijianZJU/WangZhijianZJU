close all
a=load('C:\Users\Think\Documents\WeChat Files\WangZhijianZJU\FileStorage\File\2020-01\von_poker_6_cfr_cr.csv')
%a=load('C:\Users\Think\Documents\WeChat Files\WangZhijianZJU\FileStorage\File\2020-01\von_poker_8_cfr_cr.csv')
%a=load('C:\Users\Think\Documents\WeChat Files\WangZhijianZJU\FileStorage\File\2020-01\von_poker_10_cfr_cr.csv')

a1=a(1:2:400,:)
a2=a(2:2:400,:)

figure(1)
plot(a1); title('Ann cr')
figure(2)
plot(a2); title('Bob cr')

acc1=a1(1,:);
acc2=a2(1,:);
for k=2:40
    atmp1=sum(a1(1:k,:));
    atmp2=sum(a2(1:k,:));
    acc1=[acc1;atmp1];
    acc2=[acc2;atmp2];
end

figure(3)
plot(acc1); title('Alg X')
figure(4)
plot(acc2); title('Alg Y')

