clear all 
clc
tic
%%%定义各参数
syms t
w0=3;
epsilon=2.0;
epsR=0.00001;
m=[1];
k=[1];
Cs=[cos(t) cos(3*t) cos(5*t) sin(t) sin(3*t) sin(5*t)];
%Cs=[cos(t) sin(t)];
Cs1=diff(Cs,t,1);
S=[Cs];
%A1=[2 1 1 1 1 1 1 1 1 1]';
A1=[1.5 1 1 1 1 1]'
A0=[A1];
S2=diff(S,t,2);
fm=matlabFunction(S'*m*S2);
M=quadv(fm,0,2*pi);
fk=matlabFunction(S'*k*S);
K=quadv(fk,0,2*pi);
S1=diff(S,t,1);
c=[1];
fc=matlabFunction(S'*c*S1);
C=quadv(fc,0,2*pi);
 n=0;
 tol=100;  
 while tol>epsR
         c3=diag(S*A0).^2;
         fc3=matlabFunction(S'*c3*S1);
         C3=quadv(fc3,0,2*pi);
         k2=diag(S*A0).*diag(S1*A0);
         fk2=matlabFunction(S'*k2*S);
         K2=quadv(fk2,0,2*pi);
        td=w0*0.2;
        T=[cos(td)  0          0         -sin(td)  0           0;
           0        cos(3*td)  0         0         -sin(3*td)  0;
           0        0         cos(5*td)  0         0           -sin(5*td); 
           sin(td)  0         0          cos(td)   0           0;
           0       sin(3*td)  0          0         cos(3*td)   0;
           0        0         sin(5*td)  0         0           cos(5*td)];    
         fkd=matlabFunction(S'*S*T);
         kd=quadv(fkd,0,2*pi);
         Kmc=w0^2*M+epsilon*w0*(C3-C)+K+2*epsilon*w0*K2-0.2*epsilon*kd;
         R=-(w0^2*M+epsilon*w0*(C3-C)+K-0.2*epsilon*kd)*A0;
         Rmc=-(2*w0*M+epsilon*(C3-C))*A0;
         %%%%%
          tol=norm(R)
          if(n>50)
              disp('迭代步数太多，可能不收敛')
             return;
          end
          Kmc11=-Rmc(:,1);
          Kmcr=[Kmc11 Kmc(:,2:size(Kmc,2))]; 
          AA=inv(Kmcr)*R;
          ww=AA(1); 
          A01=A0+[0.0;AA(2:length(A0),1)];
          w01=w0+ww;
          A0=A01;
          w0=w01;
          n=n+1;
 end
 A0
 w0
  X0=S*A0;
  dX0=S1*A0*w0; 
tt=0:.01:50;
xo1=subs(X0(1),tt);
dxo1=subs(dX0(1),tt);
figure(1)
plot(xo1,dxo1,'b*')
hold on

tau=0.2;
ddex1dez = @(x,y,Z) [y(2);epsilon*(1-y(1)^2)*y(2)-y(1)+epsilon*0.2*Z(1,1)];
sol = dde23(ddex1dez,[tau,tau],[0 1],[0, 50]);
tint=[0:0.01:50];
sint=deval(sol,tint);

%figure;
plot(sint(1,:),sint(2,:) )
title('limit cycle of single degree van der pol equation')
legend('IHB-1','dde23');
xlabel('x0')
ylabel('dx0')
dde=[double(sol.y(1,:)') double(sol.y(2,:)')];
IHB=[double(xo1') double(dxo1')];
IHB2=[tint' sint(1,:)'];
IHB3=[double(tt') double(xo1')];
t=xlswrite('C:\Users\sunpan\Desktop\IHB\comparison\IHB_delay_VDP\IHB_3.xls',IHB);
toc