%function [para_opt] = post_opt(I,Ival,E1,A1,b,c)

%I=[1,3];
%init_val=[1*10^(-13)];
node_list=[1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,26,27,3,4,5,6,7,8,9];
%E1=E_1;
%A1=A_1;
cap=linspace(0*10^(-12),0*10^(-12),1);
cap2=linspace(0*10^(-12),3*10^(-12),40);
cap3=linspace(0*10^(-12),0*10^(-12),1);
s1=length(cap);
s2=length(cap2);
s3=length(cap3);
bd_glo=1;
c=zeros(1,59); c(19)=1;  %% in and output vectors
b=zeros(59,1); b(59)=1;


pk_max=2;   %% contraint on max. peaking

%u1=zeros(59,1); u1(6)=1; u1(26)=-1;    
%u2=zeros(59,1); u2(10)=1; u2(16)=-1; 
%u3=zeros(59,1); u3(6)=1; u3(15)=-1;
%Capacitor between 14 and 8 with cij=3.267617e-12.
%Capacitor between 18 and 23 with cij=9.839690e-13.


%% This includes Miller capacitor
%Capacitor between 1 and 3 with cij=2.318499e-12.
%Capacitor between 14 and 8 with cij=2.885981e-12.
%Capacitor between 17 and 18 with cij=4.962241e-13.
%u1=zeros(59,1); u1(1)=1; u1(21)=-1;    
%u2=zeros(59,1); u2(6)=1; u2(26)=-1; 
%u3=zeros(59,1); u3(9)=1; u3(10)=-1;


%Capacitor between 14 and 8 with cij=1.576760e-12.
%Capacitor between 18 and 23 with cij=9.170497e-13.
%Capacitor between 3 and 8 with cij=4.618399e-13.
u1=zeros(59,1); u1(6)=1; u1(26)=-1;    
u2=zeros(59,1); u2(10)=1; u2(16)=-1; 
u3=zeros(59,1); u3(21)=1; u3(26)=-1;


for j3=1:s3
for j2=1:s2
for j1=1:s1
   
    [bd,p]=get_bandwidth(E1-cap(j1)*(u1*u1')-cap2(j2)*(u2*u2')-cap3(j3)*(u3*u3'),A1,b,c);
    
    if bd>bd_glo && p <pk_max 
        
        bd_glo=bd;
        c_glo=[cap(j1),cap2(j2),cap3(j3)];
        pk_glo=p;
        
    end
end

end
end

