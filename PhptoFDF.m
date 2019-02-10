A=[0.0001 -0.8027;-1.0135 0.1000];
Ad=[0.1993 -0.1000;0.0999 -0.1021];
B=[0.1200 -0.5950]';
C1=[0.0100 -0.0100];
C2=0.7827;
G=[2.0000 1.0007]';
D=[0.1 0.1]';
K=[0.1 0.2];
S=0.0023;
u=0.5*3^(-3);
rufa=0.8;
bata=0.9;
rou1=[rufa*(1-rufa)]^(0.5);
rou2=[bata*(1-bata)]^(0.5);
L=[-0.7045 -9.1678]';
% 系统矩阵及所求得矩阵
f1=zeros(201,1); % 故障信号
a1=zeros(201,1); % 丢包随机变量 0.8
b1=zeros(201,1); % 丢包随机变量 0.9
d1=zeros(201,1); % 外部扰动
r1=zeros(201,1); % 残差值
r2=zeros(201,1); % 残差值没有故障case
J1=zeros(201,1);  % 评估函数
J2=zeros(201,1);  % 评估函数没有故障case
J1q=zeros(201,1);  % 评估函数开方
J2q=zeros(201,1);  % 评估函数没有故障case开方
e1=zeros(2,201); % 状态差
e2=zeros(2,201); % 状态差没有故障case
x1=zeros(2,201);  % 状态变量
x2=zeros(2,201);  % 状态变量没有故障case

Temx1=zeros(2,1); % 临时变量
Temx2=zeros(2,1); % 临时变量
Teme1=zeros(2,1); % 临时变量
Teme2=zeros(2,1); % 临时变量
Tem3=0;
Tem4=0;
Asx1=zeros(2,1);
Asx2=zeros(2,1);
Bse1=zeros(2,1);
Bse2=zeros(2,1);
Cs1=0;
Cs2=0;
% 故障信号和丢包随机变量的生成
for i=1:1:201
    if(i>=31 & i<=61)
        f1(i)=4*sin(i-1);
    else
        f1(i)=0;
    end
    d1(i)=2*exp(-0.01*(i-1))*randsrc(1,1,[-0.01 0.01;1/2 1/2]);
   a1(i)=binornd(1,0.8);
   b1(i)=binornd(1,0.9); 
end
x1(:,[1])=[0 0]';
x2(:,[1])=[0 0]';
e1(:,[1])=[0 0]';
e2(:,[1])=[0 0]';
for k=1:1:200
    if(k>1)
      for d=1:1:(k-1)
          Temx1=3^(-3-d)*Ad*x1(:,[k-d]);
          Temx2=3^(-3-d)*Ad*x2(:,[k-d]);
          Teme1=2^(-3-d)*Ad*e1(:,[k-d]);
          Teme2=2^(-3-d)*Ad*e2(:,[k-d]);
          Asx1=Asx1+Temx1;
          Asx2=Asx2+Temx2;
          Bse1=Bse1+Teme1;
          Bse2=Bse2+Teme2;
      end
    else 
          Asx1=[0 0]';
          Asx2=[0 0]';
          Bse1=[0 0]';
          Bse2=[0 0]'; 
    end
       x1(:,[k+1])=(A+bata*D*K)*x1(:,[k])-bata*D*K*e1(:,[k])+[bata-b1(k)]*(-1*D*K*x1(:,[k])+D*K*e1(:,[k]))+Asx1+B*d1(k)+G*f1(k);
       x2(:,[k+1])=(A+bata*D*K)*x2(:,[k])-bata*D*K*e2(:,[k])+[bata-b1(k)]*(-1*D*K*x2(:,[k])+D*K*e2(:,[k]))+Asx2+B*d1(k);
       e1(:,[k+1])=(A-rufa*L*C1)*e1(:,[k])+[rufa-a1(k)]*L*C1*x1(:,[k])+[bata-b1(k)]*[(-D)*K*x1(:,[k])+D*K*e1(:,[k])]+Bse1+(B-L*C2)*d1(k)+G*f1(k);
       e2(:,[k+1])=(A-rufa*L*C1)*e2(:,[k])+[rufa-a1(k)]*L*C1*x2(:,[k])+[bata-b1(k)]*[(-D)*K*x2(:,[k])+D*K*e2(:,[k])]+Bse2+(B-L*C2)*d1(k);
       r1(k)=rufa*S*C1*e1(:,[k])+(rufa-a1(k))*(-S*C1*x1(:,[k]))+S*C2*d1(k);
       r2(k)=rufa*S*C1*e2(:,[k])+(rufa-a1(k))*(-S*C1*x2(:,[k]))+S*C2*d1(k);
       Tem3=r1(k)*r1(k);
       Tem4=r2(k)*r2(k);
       Cs1=Cs1+Tem3;
       Cs2=Cs2+Tem4;
       J1(k)=Cs1;
       J2(k)=Cs2;
          Asx1=0;
          Asx2=0;
          Bse1=0;
          Bse2=0;
end

r1(201)=rufa*S*C1*e1(:,[201])+(rufa-a1(201))*(-S*C1*x1(:,[201]))+S*C2*d1(201);
r2(201)=rufa*S*C1*e2(:,[201])+(rufa-a1(201))*(-S*C1*x2(:,[201]))+S*C2*d1(201);
J1(201)=J1(200)+r1(201)*r1(201);
J2(201)=J2(200)+r2(201)*r2(201);
 for k=1:1:201   
      J1q(k)=J1(k)^(0.5);
      J2q(k)=J2(k)^(0.5);
 end
%   J1q(37)
 % 开始画图
 j=1:1:201;
 plot(j-1,J1q);
  hold on;
j=1:1:201;
 plot(j-1,J2q);





