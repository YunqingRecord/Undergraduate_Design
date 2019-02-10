%FD for discrete system with missing measurements and infinite
%distributed delays
Gama=1.0001;
A=[0.0001 -0.8027;-1.0135 0.1];
Ad=[0.1993 -0.1000;0.0999 -0.1021];
B=[0.1200 -0.595]';
C1=[0.01 -0.01];
C2=0.7827;
G=[2 1.0007]';
D=[0.1 0.1]';
K=[0.1 0.2];
u=0.5*3^(-3);
rufa=0.8;
bata=0.9;
rou1=[rufa*(1-rufa)]^(0.5);
rou2=[bata*(1-bata)]^(0.5);
%matrix variables
setlmis([]);
%definite matrix
P11=lmivar(1,[2,1]);
P22=lmivar(1,[2,1]);
Q11=lmivar(1,[2,1]);
Q12=lmivar(2,[2,2]);
Q22=lmivar(1,[2,1]);
Lw=lmivar(2,[1,2]);
S=lmivar(2,[1,1]); 
% began with here
% 第一个不等式
linearmi1=newlmi;
lmiterm([linearmi1 1 1 Q11],u,1);
lmiterm([linearmi1 1 1 P11],1,-1);
lmiterm([linearmi1 1 2 Q12],u,1);
lmiterm([linearmi1 1 7 P11],bata*K'*D'+A',1);
lmiterm([linearmi1 1 10 Lw],rou1*C1',1);
lmiterm([linearmi1 1 11 P11],-1*rou2*K'*D',1);
lmiterm([linearmi1 1 12 P22],-1*rou2*K'*D',1);
lmiterm([linearmi1 1 14 S],-1*C1',rou1);
%第一行完
lmiterm([linearmi1 2 2 Q22],u,1);
lmiterm([linearmi1 2 2 P22],-1,1);
lmiterm([linearmi1 2 7 P11],-1*bata*K'*D',1);
lmiterm([linearmi1 2 8 P22],A',1);
lmiterm([linearmi1 2 8 Lw],rufa*C1',-1);
lmiterm([linearmi1 2 11 P11],rou2*K'*D',1);
lmiterm([linearmi1 2 12 P22],rou2*K'*D',1);
lmiterm([linearmi1 2 13 S],rufa*C1',1);
%第2行完
lmiterm([linearmi1 3 3 Q11],-1,u^(-1));
lmiterm([linearmi1 3 4 Q12],-1,u^(-1));
lmiterm([linearmi1 3 7 P11],Ad',1);
% 第三行完
lmiterm([linearmi1 4 4 Q22],-1,u^(-1));
lmiterm([linearmi1 4 8 P22],Ad',1);
%第5行完
lmiterm([linearmi1 5 5 0],-1*Gama);
lmiterm([linearmi1 5 7 P11],B',1);
lmiterm([linearmi1 5 8 P22],B',1);
lmiterm([linearmi1 5 8 Lw],C2',-1);
lmiterm([linearmi1 5 13 S],C2',1);
%第6行完
lmiterm([linearmi1 6 6 0],-1*Gama);
lmiterm([linearmi1 6 7 P11],G',1);
lmiterm([linearmi1 6 8 P22],G',1);
lmiterm([linearmi1 6 13 0],-1);
%第7行完
lmiterm([linearmi1 7 7 P11],-1,1);
%第8行完
lmiterm([linearmi1 8 8 P22],-1,1);
%第9行完
lmiterm([linearmi1 9 9 P11],-1,1);
%第10行完
lmiterm([linearmi1 10 10 P22],-1,1);
%第11行完
lmiterm([linearmi1 11 11 P11],-1,1);
lmiterm([linearmi1 12 12 P22],-1,1);
lmiterm([linearmi1 13 13 0],-1);
lmiterm([linearmi1 14 14 0],-1);
% 第一个不等式完;
linearmi2=newlmi;
lmiterm([-linearmi2 1 1 P11],1,1);
lmiterm([linearmi2 1 1 0],0);
% above is P11>0
linearmi3=newlmi;
lmiterm([-linearmi3 1 1 P22],1,1);
lmiterm([linearmi3 1 1 0],0);
% above is P22>0
linearmi4=newlmi;
lmiterm([-linearmi4 1 1 Q11],1,1);
lmiterm([-linearmi4 1 2 Q12],1,1);
lmiterm([-linearmi4 2 2 Q22],1,1);
lmiterm([linearmi4 1 1 0],0);
% above is Q>0
lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys)
Lw=dec2mat(lmisys,xfeas,Lw)
P11=dec2mat(lmisys,xfeas,P11)
P22=dec2mat(lmisys,xfeas,P22)
Q11=dec2mat(lmisys,xfeas,Q11)
Q12=dec2mat(lmisys,xfeas,Q12)
Q22=dec2mat(lmisys,xfeas,Q22)
S=dec2mat(lmisys,xfeas,S)
L=P22^(-1)*Lw'
tmin


