clear
clc
epsilong=2;
tau=0.2;
k=0.2;
ddex1dez = @(t,y,Z) [y(2);epsilong*(1-y(1)^2)*y(2)-y(1)+epsilong*k*Z(1,1)];
sol = dde23(ddex1dez,[tau,tau],[0 1],[0, 200]);;
tint=[0:0.01:200];
Sint=deval(sol,tint);
data=Sint(1,10001:end)';
plot(sol.y(1,:),sol.y(2,:))
title('时滞微分方程组');
a=[tint' Sint(1,:)'];
%plot(sol.x,sol.y(1,:))
plot(sol.y(1,:),sol.y(2,:))
A=[sol.x',sol.y(1,:)',sol.y(2,:)'];