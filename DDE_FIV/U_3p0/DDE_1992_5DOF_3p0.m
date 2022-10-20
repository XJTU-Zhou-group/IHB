clear
clc
tic
beta=0.24;
lamdal=zeros(5,1);
b=[4.730041 7.853205 10.995608 14.137166 17.278759];
v=zeros(5,1);
a=zeros(4,5);
gama=zeros(5,1);
phi=zeros(5,1);
U=3.0;
tau=2*pi/U;
k=1000;
for i=1:1:5
    lamdal(i)=(b(i))^2;
    v(i)=lamdal(i)/lamdal(1);
    gama(i)=(sinh(b(i))+sin(b(i)))/(cos(b(i))-cosh(b(i)));
    phi(i)=cos(b(i)*0.5)-cosh(b(i)*0.5)+(sin(b(i)*0.5)-sinh(b(i)*0.5))*(gama(i));
end
for j=1:1:5
   a(1,j)=0.0145*v(j);
   a(2,j)=0.00524;
   a(3,j)=0.76*v(j)^2;
   a(4,j)=0.026;
end

ddex1dez = @(t,y,Z) [y(2);-(a(1,1)+a(2,1)*U)*y(2)-a(3,1)*y(1)-a(4,1)*U^2*Z(1,1)-(1-beta)*k*phi(1)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(4);-(a(1,2)+a(2,2)*U)*y(4)-a(3,2)*y(3)-a(4,2)*U^2*Z(3,1)-(1-beta)*k*phi(2)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(6);-(a(1,3)+a(2,3)*U)*y(6)-a(3,3)*y(5)-a(4,3)*U^2*Z(5,1)-(1-beta)*k*phi(3)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(8);-(a(1,4)+a(2,4)*U)*y(8)-a(3,4)*y(7)-a(4,4)*U^2*Z(7,1)-(1-beta)*k*phi(4)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(10);-(a(1,5)+a(2,5)*U)*y(10)-a(3,5)*y(9)-a(4,5)*U^2*Z(9,1)-(1-beta)*k*phi(5)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3];
                 opts = ddeset('RelTol',1e-3,'AbsTol',1e-6);
                 sol = dde23(ddex1dez,[tau,tau],[0 0.01 0 0.01 0 0.01 0 0.01 0 0.01],[0, 800],opts);
tint=[0:0.05:800];
Sint=deval(sol,tint);
S=Sint(:,end-1000:end)';
eta1=phi(1)*Sint(1,:)+phi(2)*Sint(3,:)+phi(3)*Sint(5,:)+phi(4)*Sint(7,:)+phi(5)*Sint(9,:);
%A=eta1(1,end-2000:end)';
eta2=phi(1)*Sint(2,:)+phi(2)*Sint(4,:)+phi(3)*Sint(6,:)+phi(4)*Sint(8,:)+phi(5)*Sint(10,:);
eta3=phi(1)*sol.y(1,:)+phi(2)*sol.y(3,:)+phi(3)*sol.y(5,:)+phi(4)*sol.y(7,:)+phi(5)*sol.y(9,:);
eta4=phi(1)*sol.y(2,:)+phi(2)*sol.y(4,:)+phi(3)*sol.y(6,:)+phi(4)*sol.y(8,:)+phi(5)*sol.y(10,:);

figure;
%plot(sol.x,sol.y(1,:));
plot(sol.x,eta3);
%hold on
%plot(tint,eta1,'r'); 
%grid on
figure;
%plot(eta1,eta2);
%hold on
%plot(eta3,eta4,'r');
plot(eta1(end-1000:end),eta2(end-1000:end),'r');
%grid on
result=[double(eta1(end-1000:end)') double(eta2(end-1000:end)')];
%s=xlswrite('C:\Users\sunpan\Desktop\SDOF1992\4P3.xls',result);
%xlabel('displacement');
%ylabel('velocity');
%H=eta1';
toc