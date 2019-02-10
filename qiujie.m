emc=0.05;
afaba=0.95;
Aemc=[0.8 0.07 0.015;0 -0.5 -0.135;0.8 -0.16 -0.5];
A0=[1 0 0;0 -0.5 -0.135;0.8 -0.16 -0.5];
Bemcf=[0.205;-0.44;-0.31];
Bf0=[0;-0.44;-0.31];
Bemcw=[0.011;-0.32;-0.25];
Bw0=[0;-0.32;-0.25];
C=[0.05 -0.025 0];
D=1.5;%系统矩阵给定完毕
setlmis([]);
P11=lmivar(1,[3,1]);
P12=lmivar(2,[3,3]);
P22=lmivar(1,[3,1]);
S11=lmivar(2,[3,3]);
S12=lmivar(2,[3,3]);
S22=lmivar(2,[3,3]);
W=lmivar(2,[1,3]);
N=lmivar(1,[1,1]);
Gama=lmivar(1,[1,1]);
linearmi1=newlmi;
lmiterm([linearmi1 1 1 P11],-1,1);
lmiterm([linearmi1 1 2 P12],-1,1);
lmiterm([linearmi1 1 5 S11],A0',1);
lmiterm([linearmi1 1 6 S12],A0',1);
lmiterm([linearmi1 1 8 W],C',1);
lmiterm([linearmi1 1 9 W],C',1);
lmiterm([linearmi1 1 10 -N],sqrt((1-afaba)*afaba)*C',1);
 
 
lmiterm([linearmi1 2 2 P22],-1,1);
lmiterm([linearmi1 2 5 S22],A0',1);
lmiterm([linearmi1 2 5 W],-afaba*C',1);
lmiterm([linearmi1 2 6 S22],A0',1);
lmiterm([linearmi1 2 6 W],-afaba*C',1);
lmiterm([linearmi1 2 7 -N],afaba*C',1);
 
 
lmiterm([linearmi1 3 3 Gama],-1,1);
lmiterm([linearmi1 3 5 S11],Bw0',1);
lmiterm([linearmi1 3 5 S22],Bw0',1);
lmiterm([linearmi1 3 5 W],-D',1);
lmiterm([linearmi1 3 6 S12],Bw0',1);
lmiterm([linearmi1 3 6 S22],Bw0',1);
lmiterm([linearmi1 3 6 W],-D',1);
lmiterm([linearmi1 3 7 -N],D',1);
 
 
lmiterm([linearmi1 4 4 Gama],-1,1);
lmiterm([linearmi1 4 5 S11],Bf0',1);
lmiterm([linearmi1 4 5 S22],Bf0',1);
lmiterm([linearmi1 4 6 S12],Bf0',1);
lmiterm([linearmi1 4 6 S22],Bf0',1);
lmiterm([linearmi1 4 7 0],-1);
 
 
lmiterm([linearmi1 5 5 P11],1,1);
lmiterm([linearmi1 5 5 S11],-1,1,'s');
lmiterm([linearmi1 5 6 P12],1,1);
lmiterm([linearmi1 5 6 -S22],-1,1);
lmiterm([linearmi1 5 6 S12],-1,1);
 
lmiterm([linearmi1 6 6 P22],1,1);
lmiterm([linearmi1 6 6 S22],-1,1,'s');
 
lmiterm([linearmi1 7 7 0],-1);
 
lmiterm([linearmi1 8 8 P11],1,1);
lmiterm([linearmi1 8 8 S11],-1,1,'s');
lmiterm([linearmi1 8 9 P12],1,1);
lmiterm([linearmi1 8 9 -S22],-1,1);
lmiterm([linearmi1 8 9 S12],-1,1);
 
lmiterm([linearmi1 9 9 P22],1,1);
lmiterm([linearmi1 9 9 S22],-1,1,'s');
 
lmiterm([linearmi1 10 10 0],-1);
 
 
linearmi1=newlmi;                    %LMI2
lmiterm([linearmi1 1 1 P11],-1,1);
lmiterm([linearmi1 1 2 P12],-1,1);
lmiterm([linearmi1 1 5 S11],A1',1);
lmiterm([linearmi1 1 6 S12],A1',1);
lmiterm([linearmi1 1 8 W],C',1);
lmiterm([linearmi1 1 9 W],C',1);
lmiterm([linearmi1 1 10 -N],sqrt((1-afaba)*afaba)*C',1);
 
 
lmiterm([linearmi1 2 2 P22],-1,1);
lmiterm([linearmi1 2 5 S22],A1',1);
lmiterm([linearmi1 2 5 W],-afaba*C',1);
lmiterm([linearmi1 2 6 S22],A1',1);
lmiterm([linearmi1 2 6 W],-afaba*C',1);
lmiterm([linearmi1 2 7 -N],afaba*C',1);
 
 
lmiterm([linearmi1 3 3 Gama],-1,1);
lmiterm([linearmi1 3 5 S11],Bw1',1);
lmiterm([linearmi1 3 5 S22],Bw1',1);
lmiterm([linearmi1 3 5 W],-D',1);
lmiterm([linearmi1 3 6 S12],Bw1',1);
lmiterm([linearmi1 3 6 S22],Bw1',1);
lmiterm([linearmi1 3 6 W],-D',1);
lmiterm([linearmi1 3 7 -N],D',1);
 
 
lmiterm([linearmi1 4 4 Gama],-1,1);
lmiterm([linearmi1 4 5 S11],Bf1',1);
lmiterm([linearmi1 4 5 S22],Bf1',1);
lmiterm([linearmi1 4 6 S12],Bf1',1);
lmiterm([linearmi1 4 6 S22],Bf1',1);
lmiterm([linearmi1 4 7 0],-1);
 
 
lmiterm([linearmi1 5 5 P11],1,1);
lmiterm([linearmi1 5 5 S11],-1,1,'s');
lmiterm([linearmi1 5 6 P12],1,1);
lmiterm([linearmi1 5 6 -S22],-1,1);
lmiterm([linearmi1 5 6 S12],-1,1);
 
lmiterm([linearmi1 6 6 P22],1,1);
lmiterm([linearmi1 6 6 S22],-1,1,'s');
 
lmiterm([linearmi1 7 7 0],-1);
 
lmiterm([linearmi1 8 8 P11],1,1);
lmiterm([linearmi1 8 8 S11],-1,1,'s');
lmiterm([linearmi1 8 9 P12],1,1);
lmiterm([linearmi1 8 9 -S22],-1,1);
lmiterm([linearmi1 8 9 S12],-1,1);
 
lmiterm([linearmi1 9 9 P22],1,1);
lmiterm([linearmi1 9 9 S22],-1,1,'s');
 
lmiterm([linearmi1 10 10 0],-1);
%第一个LMI完

linearmi2=newlmi;
lmiterm([linearmi2 1 1 P11],-1,1);
lmiterm([linearmi2 1 2 P12],-1,1);
lmiterm([linearmi2 1 5 S11],Aemc',1);
lmiterm([linearmi2 1 6 S12],Aemc',1);
lmiterm([linearmi2 1 8 W],C',1);
lmiterm([linearmi2 1 9 W],C',1);
lmiterm([linearmi2 1 10 -N],sqrt((1-afaba)*afaba)*C',1);


lmiterm([linearmi2 2 2 P22],-1,1);
lmiterm([linearmi2 2 5 S22],Aemc',1);
lmiterm([linearmi2 2 5 W],-afaba*C',1);
lmiterm([linearmi2 2 6 S22],Aemc',1);
lmiterm([linearmi2 2 6 W],-afaba*C',1);
lmiterm([linearmi2 2 7 -N],afaba*C',1);


lmiterm([linearmi2 3 3 Gama],-1,1);
lmiterm([linearmi2 3 5 S11],Bemcw',1);
lmiterm([linearmi2 3 5 S22],Bemcw',1);
lmiterm([linearmi2 3 5 W],-D',1);
lmiterm([linearmi2 3 6 S12],Bemcw',1);
lmiterm([linearmi2 3 6 S22],Bemcw',1);
lmiterm([linearmi2 3 6 W],-D',1);
lmiterm([linearmi2 3 7 -N],D',1);


lmiterm([linearmi2 4 4 Gama],-1,1);
lmiterm([linearmi2 4 5 S11],Bemcf',1);
lmiterm([linearmi2 4 5 S22],Bemcf',1);
lmiterm([linearmi2 4 6 S12],Bemcf',1);
lmiterm([linearmi2 4 6 S22],Bemcf',1);
lmiterm([linearmi2 4 7 0],-1);


lmiterm([linearmi2 5 5 P11],1,1);
lmiterm([linearmi2 5 5 S11],-1,1,'s');
lmiterm([linearmi2 5 6 P12],1,1);
lmiterm([linearmi2 5 6 -S22],-1,1);
lmiterm([linearmi2 5 6 S12],-1,1);

lmiterm([linearmi2 6 6 P22],1,1);
lmiterm([linearmi2 6 6 S22],-1,1,'s');

lmiterm([linearmi2 7 7 0],-1);

lmiterm([linearmi2 8 8 P11],1,1);
lmiterm([linearmi2 8 8 S11],-1,1,'s');
lmiterm([linearmi2 8 9 P12],1,1);
lmiterm([linearmi2 8 9 -S22],-1,1);
lmiterm([linearmi2 8 9 S12],-1,1);

lmiterm([linearmi2 9 9 P22],1,1);
lmiterm([linearmi2 9 9 S22],-1,1,'s');

lmiterm([linearmi2 10 10 0],-1);


linearmi2=newlmi;                    %LMI2
lmiterm([linearmi2 1 1 P11],-1,1);
lmiterm([linearmi2 1 2 P12],-1,1);
lmiterm([linearmi2 1 5 S11],A1',1);
lmiterm([linearmi2 1 6 S12],A1',1);
lmiterm([linearmi2 1 8 W],C',1);
lmiterm([linearmi2 1 9 W],C',1);
lmiterm([linearmi2 1 10 -N],sqrt((1-afaba)*afaba)*C',1);


lmiterm([linearmi2 2 2 P22],-1,1);
lmiterm([linearmi2 2 5 S22],A1',1);
lmiterm([linearmi2 2 5 W],-afaba*C',1);
lmiterm([linearmi2 2 6 S22],A1',1);
lmiterm([linearmi2 2 6 W],-afaba*C',1);
lmiterm([linearmi2 2 7 -N],afaba*C',1);


lmiterm([linearmi2 3 3 Gama],-1,1);
lmiterm([linearmi2 3 5 S11],Bw1',1);
lmiterm([linearmi2 3 5 S22],Bw1',1);
lmiterm([linearmi2 3 5 W],-D',1);
lmiterm([linearmi2 3 6 S12],Bw1',1);
lmiterm([linearmi2 3 6 S22],Bw1',1);
lmiterm([linearmi2 3 6 W],-D',1);
lmiterm([linearmi2 3 7 -N],D',1);


lmiterm([linearmi2 4 4 Gama],-1,1);
lmiterm([linearmi2 4 5 S11],Bf1',1);
lmiterm([linearmi2 4 5 S22],Bf1',1);
lmiterm([linearmi2 4 6 S12],Bf1',1);
lmiterm([linearmi2 4 6 S22],Bf1',1);
lmiterm([linearmi2 4 7 0],-1);


lmiterm([linearmi2 5 5 P11],1,1);
lmiterm([linearmi2 5 5 S11],-1,1,'s');
lmiterm([linearmi2 5 6 P12],1,1);
lmiterm([linearmi2 5 6 -S22],-1,1);
lmiterm([linearmi2 5 6 S12],-1,1);

lmiterm([linearmi2 6 6 P22],1,1);
lmiterm([linearmi2 6 6 S22],-1,1,'s');

lmiterm([linearmi2 7 7 0],-1);

lmiterm([linearmi2 8 8 P11],1,1);
lmiterm([linearmi2 8 8 S11],-1,1,'s');
lmiterm([linearmi2 8 9 P12],1,1);
lmiterm([linearmi2 8 9 -S22],-1,1);
lmiterm([linearmi2 8 9 S12],-1,1);

lmiterm([linearmi2 9 9 P22],1,1);
lmiterm([linearmi2 9 9 S22],-1,1,'s');

lmiterm([linearmi2 10 10 0],-1);

lmisys=getlmis;

[tmin,xfeas]=feasp(lmisys);
tmin

n=decnbr(lmisys);
c=zeros(n,1);
     for j=1:n  
[gaj]=defcx(lmisys,j,Gama);
c(j)=trace(gaj);
     end
% 上述是为了获取c(普适方法)
options=[1e-5,100,0,20,0];
[copt,xopt]=mincx(lmisys,c,options);
P11opt=dec2mat(lmisys,xopt,P11);
P12opt=dec2mat(lmisys,xopt,P12);
P22opt=dec2mat(lmisys,xopt,P22);
Nopt=dec2mat(lmisys,xopt,N);
Wopt=dec2mat(lmisys,xopt,W);
S11opt=dec2mat(lmisys,xopt,S11);
S12opt=dec2mat(lmisys,xopt,S12);
S22opt=dec2mat(lmisys,xopt,S22);
L=[S22opt^(-1)]'*Wopt'
N=Nopt
% the object value is copt
copt