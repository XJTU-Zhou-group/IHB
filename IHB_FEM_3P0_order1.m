clear
clc
tic
%%%定义各参数
syms t
epsR=0.00007;
w0=1.4; %初始角频率
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
q=1;
for l=1:1:q
    Cs(l)=cos((2*l-1)*t);
    Cs(q+l)=sin((2*l-1)*t);
end
S=blkdiag(Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs,Cs);
A0=zeros(36,1);%初值
load('C.mat');
A0=[0.00404431325681056;-0.00269232252568129;0.0716601499426500;-0.0477263126866847;0.0125555522767557;-0.00836643089661432;0.0903472287409401;-0.0602990220450916;0.0207204012904830;-0.0138253539344120;0.0676444522916330;-0.0453394540045132;0.0254686369994998;-0.0170197178528385;0.0267783079959394;-0.0181476314052987;0.0265413879976585;-0.0177534948598000;-1.34537722637446e-15;6.09833410378476e-16;0.0254686369994996;-0.0170197178528384;-0.0267783079959416;0.0181476314052996;0.0207204012904827;-0.0138253539344119;-0.0676444522916331;0.0453394540045133;0.0125555522767554;-0.00836643089661423;-0.0903472287409385;0.0602990220450911;0.00404431325681045;-0.00269232252568126;-0.0716601499426482;0.0477263126866832];
A0(19:20)=[0;0];



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
    if(n>400)
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
plot(xo9,dxo9,'o','linewidth',1)
hold on
plot(data(end-10000:end,2),data(end-10000:end,3),'-.r','linewidth',2);
grid on
legend('IHB-1','dde23');
xlabel('x0')
ylabel('dx0')
result1=[double(xo9') double(dxo9')]
toc