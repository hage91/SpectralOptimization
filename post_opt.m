%function [para_opt] = post_opt(I,Ival,E1,A1,b,c)

%I=[1,3];
%init_val=[1*10^(-13)];
E1=E_RC;
A1=A_RC;
cap=linspace(0,1.5*10^(-12),20);
cap2=linspace(-9*10^(-13),5*10^(-13),20);
s1=length(cap);
s2=length(cap2);
bd_glo=1;
c=zeros(1,12); c(6)=1;  %% in and output vectors
b=zeros(12,1); b(12)=1;
%[s1,s2]=size(I);
%b
%c
pk_max=1;
u=zeros(12,1); u(1)=1; u(9)=-1;
u2=zeros(12,1); u2(6)=1; u2(9)=-1;


for j2=1:s2
for j1=1:s1
   
    [bd,p]=get_bandwidth(E1-cap(j1)*u*u'-cap2(j2)*u2*u2',A1,b,c);
    
    if bd>bd_glo && p <pk_max 
        
        bd_glo=bd;
        c_glo=[cap(j1),cap2(j2)];
        pk_glo=p;
        
    end
end

end

