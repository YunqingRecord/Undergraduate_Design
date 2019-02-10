Gama=1.3500;
emc=0.05;
afaba=0.95;
Aemc=[0.8 0.07 0.015;0 -0.5 -0.135;0.8 -0.16 -0.5];
Bemcf=[0.205;-0.44;-0.31];
Bemcw=[0.011;-0.32;-0.25];
C=[0.05 -0.025 0];
D=1.5;
L=1.0e-12 *[ -0.0021 0.1894,0.1809]';
N= -7.0061e-14;%ϵͳ�����˲��������������


f1=zeros(201,1); % �����ź�
afa=zeros(201,1); % ����������� 0.95
d1=zeros(201,1); % �ⲿ�Ŷ�
r1=zeros(201,1); % �в�ֵ
r2=zeros(201,1); % �в�ֵû�й���case
J1=zeros(201,1);  % ��������
J2=zeros(201,1);  % ��������û�й���case
J1q=zeros(201,1);  % ������������
J2q=zeros(201,1);  % ��������û�й���case����
e1=zeros(3,201); % ״̬��
e2=zeros(3,201); % ״̬��û�й���case
x1=zeros(3,201);  % ״̬����
x2=zeros(3,201);  % ״̬����û�й���case
Tem3=0;
Tem4=0;
Cs1=0;
Cs2=0;


for i=1:1:201
    if(i>=31 & i<=61)
        f1(i)=4*sin(i-1);
    else
        f1(i)=0;
    end
    w1(i)=2*exp(-0.01*(i-1))*randsrc(1,1,[-0.01 0.01;1/2 1/2]);
   afa(i)=binornd(1,0.95);
end              % �����źźͶ����������������

x1(:,[1])=[0 0 0]';
x2(:,[1])=[0 0 0]';
e1(:,[1])=[0 0 0]';
e2(:,[1])=[0 0 0]';  %��ʼ״̬

%������״̬����
for k=1:1:200
x1(:,[k+1])=Aemc*x1(:,[k])+Bemcw*w1(k)+Bemcf*f1(k);
x2(:,[k+1])=Aemc*x2(:,[k])+Bemcw*w1(k);
e1(:,[k+1])=(Aemc-afaba*L*C)*e1(:,[k])-(afa(k)-afaba)*L*C*x1(:,[k])+(Bemcw-L*D)*w1(k)+Bemcf*f1(k);
e2(:,[k+1])=(Aemc-afaba*L*C)*e2(:,[k])-(afa(k)-afaba)*L*C*x2(:,[k])+(Bemcw-L*D)*w1(k) ;

r1(k)=afaba*N*C*e1(:,[k])+(afa(k)-afaba)*N*C*x1(:,[k])+N*D*w1(k);
r2(k)=afaba*N*C*e2(:,[k])+(afa(k)-afaba)*N*C*x2(:,[k])+N*D*w1(k);
Tem3=r1(k)*r1(k);
       Tem4=r2(k)*r2(k);
       Cs1=Cs1+Tem3;
       Cs2=Cs2+Tem4;
       J1(k)=Cs1;
       J2(k)=Cs2;
end
r1(201)=afaba*N*C*e1(:,[201])+(afa(201)-afaba)*N*C*x1(:,[201])+N*D*w1(201);
r2(201)=afaba*N*C*e2(:,[201])+(afa(201)-afaba)*N*C*x2(:,[201])+N*D*w1(201);
J1(201)=J1(200)+r1(201)*r1(201);
J2(201)=J2(200)+r2(201)*r2(201);
 for k=1:1:201   
      J1q(k)=J1(k)^(0.5);
      J2q(k)=J2(k)^(0.5);
 end
%   J1q(37)
 % ��ʼ��ͼ
figure(1)
j=1:1:201;
 plot(j-1,J1q);
  hold on;
j=1:1:201;
 plot(j-1,J2q); %��ֵ&��������
 xlabel('Time step (k)');
ylabel('evaluation function');   %�в�


 figure(2)
j=1:1:201;  %r(k)
plot(j-1,r1);
hold on;
j=1:1:201;
plot(j-1,r2,'--');
xlabel('Time step (k)');
ylabel('r(k)');   %�в�


