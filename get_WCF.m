function [S,T,iplus,iminus,invipl,invimin,V1] = get_WCF(Et,At,ity,num)
E=inv(At)*Et;
sz=length(E);     

A=eye(sz);

[U1,V1]=eig(A,E);

[~,I1]=sort(diag(V1));
 
U1perm=zeros(sz,sz);
   for i=1:sz
U1perm(:,i)=U1(:,I1(i));
   end

[U2,V2]=eig(A',E');

[~,I2]=sort(diag(V2));
 
 U2perm=zeros(sz,sz);
   for i=1:sz
U2perm(:,i)=U2(:,I2(i));
   end

 


U1=U1perm;
U2=U2perm';


iplus=zeros(num,1);   %%%%Indices of dominant EigVal
iminus=zeros(num,1);

Icplxpl=zeros(sz,1);
Icplxmin=zeros(sz,1);
cntp=1;
cntm=1;
%% für allgemeine komplexe Zahlen natürlich quatsch...
for i=1:sz
    if imag(V1(i,i))>10^3      
    Icplxpl(cntp)=i;
    cntp=cntp+1;
    end
    if imag(V1(i,i))<-10^3 
    Icplxmin(cntm)=i;   
    cntm=cntm+1;
    end
end




Invcplxmin=zeros(sz,1);
Invcplxpl=zeros(sz,1);

for i=1:sz
    for cnt=1:sz
        if I1(i)==Icplxpl(cnt)
        Invcplxpl(cnt)=i;
        end
    end
   for cnt=1:sz
    if I1(i)==Icplxmin(cnt)
    Invcplxmin(cnt)=i;
    end
   end 
end

%% diagonalize the part at infinity

Ed=U2*E*U1;  
Ad=U2*A*U1;

N1=Ad([(sz-ity+1):sz],:);
N=N1(:,[(sz-ity+1):sz]);

[UN,SN,VN]=svd(N); %%% dann ist UN'*N*VN=SN

H=VN*sqrt(inv(SN));
PV=blkdiag(eye(sz-ity),H);
PU=blkdiag(eye(sz-ity),sqrt(inv(SN))*UN');

U2=PU*U2;
U1=U1*PV;

vh=diag(Ed);
Pfin=blkdiag(inv(diag(sqrt(abs(vh(1:(sz-ity)) )))),eye(ity));

U2=Pfin*U2;
U1=U1*Pfin;

S=U2;
T=U1;

%% swap the columns at conjugate complex eigenvalues
 for i=1:sz
 
if Invcplxmin(i)~=0
 x=S(Invcplxmin(i),:);
S(Invcplxmin(i),:)=S(Invcplxpl(i),:);
S(Invcplxpl(i),:)=x;        
     end
 end

d=diag(S*E*T);

for i=1:sz-ity
    S(i,:)=S(i,:)*d(i)^(-1);
end

%% select a critical eigenvalue 
%% work with function angle() here... select the largest angle...
max=pi*ones(num,1);%zeros(num,1); 
for i=1:sz
    if Icplxpl(i)~=0
        pos1=1;
       for j1=1:num
       if (angle(V1(Icplxpl(i),Icplxpl(i)))<max(j1)) %&& (angle(V1(Icplxpl(i),Icplxpl(i)))>0)
       %if abs(imag(V1(Icplxpl(i),Icplxpl(i)))/real(V1(Icplxpl(i),Icplxpl(i))))>max(j1)
            pos1=j1;  
        end
       end
              if (angle(V1(Icplxpl(i),Icplxpl(i)))<max(pos1)) %&& (angle(V1(Icplxpl(i),Icplxpl(i)))>0)
     %  if abs(imag(V1(Icplxpl(i),Icplxpl(i)))/real(V1(Icplxpl(i),Icplxpl(i))))>max(pos1)
       for j2=1:(pos1-1) 
                 max(j2)=max(j2+1);
                 iplus(j2)=iplus(j2+1);
                 iminus(j2)=iminus(j2+1);
       end
            iplus(pos1)=Icplxpl(i);
            iminus(pos1)=Icplxmin(i);
            max(pos1)=angle(V1(Icplxpl(i),Icplxpl(i)));%abs(imag(V1(Icplxpl(i),Icplxpl(i)))/real(V1(Icplxpl(i),Icplxpl(i))));
       end
    end
end


invipl=zeros(59,1);
invimin=zeros(59,1);

for j=1:num
    for i=1:sz
        if Icplxpl(i)==iplus(j)
        invipl(iplus(j))=Invcplxpl(i);
        invimin(iminus(j))=Invcplxmin(i);
        end
    end
end




end

