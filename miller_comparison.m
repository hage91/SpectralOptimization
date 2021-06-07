%% different scenarios, prescribe a different peaking and bandwidth 
%% and choose miller suitable
%% ggf müssen wir nochmal die parameter nachoptimieren



c=zeros(1,59); c(19)=1;  %%%entspricht Knoten 26 im Bild
b=zeros(59,1); b(59)=1;

%% miller position
d=zeros(59,1);
d(9)=1; d(10)=-1;

bodeplot(E1,A1,b,c,[1,0,0]);

hold on;

%bodeplot(Eglo,A1,b,c,[0,1,0]);

bodeplot(E1-2.5*10^(-12)*(d*d'),A1,b,c,[0,0,1]);


%err: 0.16223 c_ij: 1.1827e-12 (i,j): 2 25 Re/Im: 6.0039
%err: 0.045315 c_ij: 1.1655e-12 (i,j): 1 25 Re/Im: 3.6529
%err: 0.049598 c_ij: 1.2672e-12 (i,j): 2 16 Re/Im: 2.5749

r1=zeros(59,1); r1(2)=1; r1(25)=-1;
r2=zeros(59,1); r2(1)=1; r2(25)=-1;
r3=zeros(59,1); r3(2)=1; r3(16)=-1;

bodeplot(E1-1.2*10^(-12)*(r1*r1')-1.2*10^(-12)*(r2*r2')-1.2*10^(-12)*(r3*r3'),A1,b,c,[0,0,1]);



%err: 0.089302 c_ij: 1.4853e-12 (i,j): 6 25 Re/Im: 5.7682
%err: 4.2251e-06 c_ij: 1.9036e-12 (i,j): 2 17 Re/Im: 2.3558

r1=zeros(59,1); r1(6)=1; r1(25)=-1;
r2=zeros(59,1); r2(2)=1; r2(17)=-1;

bodeplot(E1-1.4*10^(-12)*(r1*r1')-1.9*10^(-12)*(r2*r2'),A1,b,c,[1,0,1]);

