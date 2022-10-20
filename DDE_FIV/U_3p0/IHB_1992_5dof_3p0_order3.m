clear all 
clc
tic
%%%定义各参数
syms t
epsR=0.0001;
w0=1.4; %初始角频率
beta=0.24;
U=3.0; % 流速
tau=2*pi/U; % 时间延迟
k1=1000; % 考虑间隙简化的非线性立方弹簧刚度
lamdal=[4.730041 7.853205 10.995608 14.137166 17.278759]'; %两端固支梁固有频率
v=lamdal.^2/lamdal(1)^2; % 某一阶段固有频率与第一阶固有频率平方比
a1=0.0145*v;
a2=0.00524;
a3=0.76*v.^2;
a4=0.026;
m=eye(5);
c=(a1+a2*U).*eye(5);
k=a3.*eye(5);
e=a4*U^2*eye(5);

gamma=(sinh(lamdal)+sin(lamdal))./(cos(lamdal)-cosh(lamdal));%以下两式为中点除的振型
phi=cos(lamdal*0.5)-cosh(lamdal*0.5)+(sin(lamdal*0.5)-sinh(lamdal*0.5)).*(gamma);

f=(1-beta)*phi*k1;
q=3;
for l=1:1:q
    Cs(l)=cos((2*l-1)*t);
    Cs(q+l)=sin((2*l-1)*t);
end
S=blkdiag(Cs,Cs,Cs,Cs,Cs);
A1=[0.0176 -0.002 0.0001 0.0142 0.0004 0.0001]';
A2=[0 0 0 0 0 0]';
A3=[0.0011 0.0007 0.0001 0.0009 0.0015 0.0001 ]';
A4=[0 0 0 0 0 0]';
A5=[0.0002 0.0001 0.0001 0.0001 0.0001 0.0001]';
A0=[A1;A2;A3;A4;A5];
S1=diff(S,t,1);
S2=diff(S,t,2);
%fm=inline(S'*m*S2);
fm= matlabFunction(S'*m*S2);
%M=quadv(fm,0,2*pi);
M = integral(fm,0,2*pi,'ArrayValued',true);
%fc=inline(S'*c*S1);
fc=matlabFunction(S'*c*S1);
%C=quadv(fc,0,2*pi);
C=integral(fc,0,2*pi,'ArrayValued',true);
%fk=inline(S'*k*S);
fk=matlabFunction(S'*k*S);
%K=quadv(fk,0,2*pi);
K=integral(fk,0,2*pi,'ArrayValued',true);
n=0;
tol=100;
while (tol>epsR)||(w0<0) 
    k3=f*diag(phi'*S*A0).^2*phi';
    %fk3=inline(S'*k3*S);
    fk3=matlabFunction(S'*k3*S);
    %K3=quadv(fk3,0,2*pi);
    K3=integral(fk3,0,2*pi,'ArrayValued',true);
    taud=w0*tau;
    for r=1:1:q
        v1(r)=-sin((2*r-1)*taud);
        v2(r)=cos((2*r-1)*taud);
        v2(q+r)=cos((2*r-1)*taud);
        v3(r)=sin((2*r-1)*taud);
    end
    gamm=diag(v2)+diag(v1,q)+diag(v3,-q);
    Gamma=blkdiag(gamm,gamm,gamm,gamm,gamm);
    %fe=inline(S'*e*S*Gamma);
    fe=matlabFunction(S'*e*S*Gamma);
    %E=quadv(fe,0,2*pi);
    E=integral(fe,0,2*pi,'ArrayValued',true);
    Kmc=w0^2*M+w0*C+K+E+3*K3;
    R=-(w0^2*M+w0*C+K+E+K3)*A0;
    Rmc=-(2*w0*M+C)*A0;
    tol=norm(R)
    if(n>1000)
        disp('迭代步数太多，可能不收敛')
        return;
    end
    Kmc11=-Rmc(:,1);
    Kmcr=[Kmc11 Kmc(:,2:size(Kmc,2))];
    AA=inv(Kmcr)*R;
    ww=AA(1);
    A01=A0+[0.0;AA(2:length(A0),1)];
    A0=A01;
    w01=w0+ww;
    w0=w01;
    n=n+1
end
A0
w0
ww

X0=S*A0;
dX0=S1*A0*w0;
tt=0:.01:10;
xo1=subs(X0(1),tt);
xo2=subs(X0(2),tt);
xo3=subs(X0(3),tt);
xo4=subs(X0(4),tt);
xo5=subs(X0(5),tt);
dxo1=subs(dX0(1),tt);
dxo2=subs(dX0(2),tt);
dxo3=subs(dX0(3),tt);
dxo4=subs(dX0(4),tt);
dxo5=subs(dX0(5),tt);

eta=phi(1)*xo1+phi(2)*xo2+phi(3)*xo3+phi(4)*xo4+phi(5)*xo5;
eta2=phi(1)*dxo1+phi(2)*dxo2+phi(3)*dxo3+phi(4)*dxo4+phi(5)*dxo5;
IHB=[double((tt/w0)') double(eta')];
ddex1dez = @(x,y,Z) [y(2);-(a1(1)+a2*U)*y(2)-a3(1)*y(1)-a4*U^2*Z(1,1)-(1-beta)*k1*phi(1)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(4);-(a1(2)+a2*U)*y(4)-a3(2)*y(3)-a4*U^2*Z(3,1)-(1-beta)*k1*phi(2)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(6);-(a1(3)+a2*U)*y(6)-a3(3)*y(5)-a4*U^2*Z(5,1)-(1-beta)*k1*phi(3)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(8);-(a1(4)+a2*U)*y(8)-a3(4)*y(7)-a4*U^2*Z(7,1)-(1-beta)*k1*phi(4)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3;
                     y(10);-(a1(5)+a2*U)*y(10)-a3(5)*y(9)-a4*U^2*Z(9,1)-(1-beta)*k1*phi(5)*(phi(1)*y(1)+phi(2)*y(3)+phi(3)*y(5)+phi(4)*y(7)+phi(5)*y(9))^3];
                 opts = ddeset('RelTol',1e-3,'AbsTol',1e-6);
                 sol = dde23(ddex1dez,[tau,tau],[0.001 0 0.001 0 0.001 0 0.001 0 0.001 0],[0, 300],opts);
eta11=phi(1)*sol.y(1,:)+phi(2)*sol.y(3,:)+phi(3)*sol.y(5,:)+phi(4)*sol.y(7,:)+phi(5)*sol.y(9,:);
eta12=phi(1)*sol.y(2,:)+phi(2)*sol.y(4,:)+phi(3)*sol.y(6,:)+phi(4)*sol.y(8,:)+phi(5)*sol.y(10,:);

figure(1)
plot(eta,eta2,'b','linewidth',1)
hold on
plot(eta11(5000:end),eta12(5000:end),'r','linewidth',1);
grid on
%title('范德波极限环')
legend('IHB-3','dde23');
xlabel('x0')
ylabel('dx0')
result=[double(eta') double(eta2')];
s=xlswrite('C:\Users\sunpan\Desktop\IHB\comparison\IHB_1992_5dof\order_3\IHB_3P0_order3.xls',result);
toc