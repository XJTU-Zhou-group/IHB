clear
clc
tic
%%%定义各参数
syms t
epsR=0.00002;
w0=1.35; %初始角频率
beta=0.24;
U=3.0; % 流速
tau=2*pi/U; % 时间延迟
k1=1000; % 考虑间隙简化的非线性立方弹簧刚度
%%
load('mm2.mat');
m=mm2;
load('cc2.mat');
c=cc2;
load('kk2.mat');
k=kk2;
load('ff2.mat');
e=-ff2;

f=zeros(18,18);
f(9,9)=(1-beta)*k1;
%%
q=3;
for l=1:1:q
    Cs(l)=cos((2*l-1)*t);
    Cs(q+l)=sin((2*l-1)*t);
end
S=blkdiag(Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs);
A0=zeros(108,1);%初值
% load('C.mat');
% A0(1:6:end,1)=0.8*C;
% A0(4:6:end,1)=0.6*C;
%A0(1:3:end,1)=[0.00404431325681056;-0.00427308491917423;0.0712370859462629;-0.0753548107488114;0.0123973947108985;-0.0131316117164316;0.0873438026875697;-0.0929101570772176;0.0201043200318101;-0.0213698973462907;0.0616495217524868;-0.0663860724024350;0.0242021618791352;-0.0258351067412279;0.0205320662812020;-0.0229958035201260;0.0248929651998753;-0.0266446380814026;-1.45036320274085e-15;4.74410891601538e-17;0.0242021618791350;-0.0258351067412280;-0.0205320662812042;0.0229958035201265;0.0201043200318097;-0.0213698973462906;-0.0616495217524870;0.0663860724024350;0.0123973947108982;-0.0131316117164315;-0.0873438026875681;0.0929101570772174;0.00404431325681045;-0.00427308491917422;-0.0712370859462610;0.0753548107488111];
%load('A0.mat');
load('B_5.mat');
A0=B_5;
A0(55:60)=[0;0;0;0;0;0];


S1=diff(S,t,1);
S2=diff(S,t,2);
fm= matlabFunction(S'*m*S2);
M = integral(fm,0,2*pi,'ArrayValued',true);
fc=matlabFunction(S'*c*S1);
C=integral(fc,0,2*pi,'ArrayValued',true);
fk=matlabFunction(S'*k*S);
K=integral(fk,0,2*pi,'ArrayValued',true);
n=0;
tol=100;
%%
while (tol>epsR)||(w0<0) 
    k3=f*diag(S*A0).^2;
    fk3=matlabFunction(S'*k3*S);
    K3=integral(fk3,0,2*pi,'ArrayValued',true);
    taud=w0*tau;
    for r=1:1:q
        v1(r)=-sin((2*r-1)*taud);
        v2(r)=cos((2*r-1)*taud);
        v2(q+r)=cos((2*r-1)*taud);
        v3(r)=sin((2*r-1)*taud);
    end
    gamm=diag(v2)+diag(v1,q)+diag(v3,-q);
    Gamma=blkdiag(gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm,gamm);
    fe=matlabFunction(S'*e*S*Gamma);
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
    w0=w01
    n=n+1
end
A0
w0
ww

%%
X0=S*A0;
dX0=S1*A0*w0;
tt=0:.01:10;
xo9=subs(X0(9),tt);
dxo9=subs(dX0(9),tt);

load('data3.mat');

figure(1)
plot(xo9,dxo9,'b','linewidth',2)
hold on
plot(data(end-10000:end,2),data(end-10000:end,3),'-.r','linewidth',2);
grid on
legend('IHB-3','dde23');
xlabel('x0')
ylabel('dx0')
result1=[double(xo9') double(dxo9')];
result2=data(end-2000:end,2:3);
%s=xlswrite('IHB_3P0_order3.xls',result);
toc