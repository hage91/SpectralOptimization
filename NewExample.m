function [Eo,iglo,jglo,tglo,gloquo]=NewExample(Et,At,ratio,domsens,err,cmax)

%%%%%%%
%TO DO bei der definition von fac müssen die sensitiven wieder
%hinzugenommen werden; 
%TO DO 2: Suche über den ganzen Bereich
%TO DO 3: Schleife außen rum sodass einfach alle möglichen Kritischen mit
%jeweiligen tmin gesprichert werden und am ende alle möglichkeiten getestet werden
%dadurch das man tmin speichert muss man keine Parametersimulation machen
%%%%%%%%%%


%%%%%%%%%
% TO DO 2: Projektionsansatz testen, wo der ideale rang eins runterprojeziert wird 
%  auf den unterraum 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung der Transformationsmatrizen S und T
%%%%%%%%%%%%%%%

%%%%%%%%PRECONDITIONING STEP
E=inv(At)*Et;
[sz,s2]=size(E);     %%%get matrix size

%sz=s1;
A=eye(sz);

%E=E1;
%A=A1;

Iout=zeros(100,3);
Ioutcnt=1;

%ratio=0.8; %%%Vorgabe der gerade mit Im/Re

%num=2; %%% Zahl der Sensitiven Eigenwerte

%domsens=5*10^10;  %%%Schwellen-Sensitivität der dominanten Eigenwerte



ity=6; %%%Zahl der Eigenwerte bei oo

num=1; %%% Zahl der Sensitiven Eigenwerte

gloquo=10^(14); %%quotient aus Fehler und Schrittweite 

gminglo=10^(14);
tglo=0;

%Irel=zeros(2,200);  %%%ausgabe: Relevante indizes %größe per hand heuristisch
glcnt=1;

step=5*10^(-14); % size of the constant maximal time step using in approximations

[U1,V1]=eig(A,E);
%plot(v,'o');
%hold on;
iglo=0;
jglo=0;
%for i=1:1000
%v=2*rand(27,1);
%u=rand(27,1);
%w=eig(A,E+10*v*v');
%plot(w,'x')
%hold on;
%end

%%Allgemeines umsortieren sonst keine iteration möglich
  % vector I1 enthält die reihenfolge
   [a1,I1]=sort(diag(V1));
 
U1perm=zeros(sz,sz);
   for i=1:sz
U1perm(:,i)=U1(:,I1(i));
   end

[U2,V2]=eig(A',E');

[a2,I2]=sort(diag(V2));
 
 U2perm=zeros(sz,sz);
   for i=1:sz
U2perm(:,i)=U2(:,I2(i));
   end



U1=U1perm;
U2=U2perm';


iplus=zeros(num,1);   %%%%Indices of sensitive EigVal
iminus=zeros(num,1);
invimin=zeros(num,1);
invipl=zeros(num,1);

Icplxpl=zeros(sz,1);
Icplxmin=zeros(sz,1);
cnt=1;
%%bestimmen komplexe indices
%geht noch nicht für mehrere komplexe indices
for i=1:sz
    if imag(V1(i,i))>10^3      
   % iplus=i; 
    Icplxpl(cnt)=i;
    
    for j=1:sz
        if imag(V1(j,j))==-imag(V1(i,i)) 
           Icplxmin(cnt)=j;    
        end
    end
    
    cnt=cnt+1;
    end
%    if imag(V1(i,i))<-10^3
 %       iminus=i;
  %          Icplxpl(cnt)=i+1;
  %  Icplxmin(cnt)=i;
   %     cnt=cnt+1;
    %end
end

Invcplxmin=zeros(sz,1);

Invcplxpl=zeros(sz,1);
%cnt=1;


for i=1:sz
    for cnt=1:sz
        if I1(i)==Icplxpl(cnt)%iplus
        %invipl=i;
        Invcplxpl(cnt)=i;
        end
    end
   for cnt=1:sz
    if I1(i)==Icplxmin(cnt)%iminus
    % invimin=i;
    Invcplxmin(cnt)=i;
    end
   end 
end



%iplus=Icplxpl(1);
%iminus=Icplxmin(1);
%invipl=Invcplxpl(1);
%invimin=Invcplxmin(1);

Ed=U2*E*U1;  %*10^(13) als Vorfaktor rausgenommen
Ad=U2*A*U1;
%%% Problem könnte noch sein, dass die Diagonalmatrizen eigentlich laut
%%% theorie die eigenwerte enthalten müssten, was sie aber nicht tun wegen
%%% numerischer fehler...

%%% Reskalierung der Transformationsmatrizen, sodass Ad die eigenvektoren
%%% enthält 
%%%% aufgrund numerischer Fehler skalieren wir zunächst nicht auf Basis von
%%%% E

D1=diag(diag(Ad));
D2=diag(diag(Ed));


%%%%%%%%
% TEIL bei unendlich diagonalisieren! Um später die Weierstrassform 
%%%%%%%

N1=Ad([(sz-ity+1):sz],:);
N=N1(:,[(sz-ity+1):sz]);

[UN,SN,VN]=svd(N); %%% dann ist UN'*N*VN=SN

H=VN*sqrt(inv(SN));
PV=blkdiag(eye(sz-ity),H);
PU=blkdiag(eye(sz-ity),sqrt(inv(SN))*UN');

U2=PU*U2;
U1=U1*PV;

vh=diag(Ed);
Pfin=blkdiag(inv(diag(sqrt( abs(vh(1:(sz-ity)) )))),eye(ity));

U2=Pfin*U2;
U1=U1*Pfin;

S=U2;
T=U1;


%%% Die Matrix S simmt noch nicht, die Zeilen zu den konjugiert komplexen
%%% EW müssen getauscht werden::

 for i=1:sz
 
     if Invcplxmin(i)~=0
 x=S(Invcplxmin(i),:);
S(Invcplxmin(i),:)=S(Invcplxpl(i),:);
S(Invcplxpl(i),:)=x;        
%x=S(invimin,:);
%S(invimin,:)=S(invipl,:);
%S(invipl,:)=x;
     end
 end

d=diag(S*E*T);

%%vorher iplus=17 und iminus=18

%S(iplus,:)=S(iplus,:)*d(iplus)^(-1);
%S(iminus,:)=S(iminus,:)*d(iminus)^(-1);

for i=1:sz-ity
    S(i,:)=S(i,:)*d(i)^(-1);
end

%%%% Hier die schleifen geschickt zusammenfassen
%for i=19:20
%    S(i,:)=S(i,:)*d(i)^(-1);
%millerend

Iinv=zeros(sz,1);
for i=1:sz
  for j=1:sz
    if I1(i)==j
    Iinv(j)=i;
    end
  end

end


%%%% select a critical eigenvalue %%% later: optimize over all EVal

max=zeros(num,1); %%die k-größten Quotienten im/re
%%Definieren einer Schwelle um dominante Eigenwerte für die Optimierung auszuwählen 
for i=1:sz
    if Icplxpl(i)~=0
        pos1=1;
       for j1=1:num
        if abs(imag(V1(Icplxpl(i),Icplxpl(i)))/real(V1(Icplxpl(i),Icplxpl(i))))>max(j1)
            pos1=j1;  
        end
       end
       if abs(imag(V1(Icplxpl(i),Icplxpl(i)))/real(V1(Icplxpl(i),Icplxpl(i))))>max(pos1)
       for j2=1:(pos1-1)%%richtiges Einsortieren in Liste
                % if j2<(pos1-1)
                 max(j2)=max(j2+1);
                 iplus(j2)=iplus(j2+1);
                 iminus(j2)=iminus(j2+1);
              end
%            iplus2=iplus;
 %           iminus2=iminus;
            iplus(pos1)=Icplxpl(i);
            iminus(pos1)=Icplxmin(i);
            max(pos1)=abs(imag(V1(Icplxpl(i),Icplxpl(i)))/real(V1(Icplxpl(i),Icplxpl(i))));
       end
    end
end
%         end
 %   end
%end


invipl=zeros(12,1);
invimin=zeros(12,1);

for j=1:num
    for i=1:sz
        if Icplxpl(i)==iplus(j)
        invipl(iplus(j))=Invcplxpl(i);
        invimin(iminus(j))=Invcplxmin(i);
        end
    end
end



%%%%%%%%%%%ACHTUNG HIER WURDE sz-ity ersetzt um den LAUFBEREICH zu
%%%%%%%%%%%vergrößern

for i=1:sz
   for j=i:sz
       v=zeros(sz,1);
        v(j)=-1; v(i)=1; 
        uinf=S*(At\v);
        vinf=transpose(T)*v;
         g1=3*10^8; gmin=g1; tmin=0; mmin=0; %mmin2=0;
       
         Ute=vinf*uinf';
%%% Hier einen sensitivitätscheck einbauen 
%%% erstmal schauen ob sich das dominante paar stark verändert
if abs(Ute(invipl(iplus(1:num)),invipl(iplus(1:num))))>=domsens   %>=0.01
        
   % abs(Ute(invipl(iplus(1:num)),invipl(iplus(1:num))))
         
    sens=zeros(sz,1);        
for b=1:sz
    if abs(Ute(b,b))>domsens
        sens(b)=b;
    end
    
end


%sens
%%% Definieren der Laufbereiche für die dominanten Eigenwerte
Csz=40; %vorher 40
num2=6;
C=zeros(num2,Csz);
Cplx=zeros(num2,3*Csz);
%%%Lazy variant: redefine iplus

iplus2=[4,5,6,7,8,9];


for k=1:num2    
    
   ang=angle(V1(iplus2(k),iplus2(k)));
   av=abs(V1(iplus2(k),iplus2(k)));
%    box1=real(V1(iplus2(k),iplus2(k)));
 %   box2=-imag(V1(iplus2(k),iplus2(k)));
  %  boxm=-min(-box1,-box2);
%C(k,:)=linspace(-boxm,-boxm+10^8,Csz); 
C(k,:)=linspace(av-3*10^7,av+3*10^7,Csz); 
 
if ang==pi
   rat2=linspace(ang,3/4*ang,3);
   Cplx(k,:)=[exp(1i*(rat2(1)+pi/2))*C(k,:),exp(1i*(rat2(2)+pi/2))*C(k,:),exp(1i*(rat2(3)+pi/2))*C(k,:)];
end

if ang<pi
   rat2=linspace(ang,pi,3);
    Cplx(k,:)=[exp(1i*(rat2(1)+pi/2))*C(k,:),exp(1i*(rat2(2)+pi/2))*C(k,:),exp(1i*(rat2(3)+pi/2))*C(k,:)];
end

if iplus2(k)==iplus(1)
rat2=linspace(ratio,2*ratio,3);
    Cplx(k,:)=[exp(1i*(rat2(1)+pi/2))*C(k,:),exp(1i*(rat2(2)+pi/2))*C(k,:),exp(1i*(rat2(3)+pi/2))*C(k,:)];
    
end

if iplus2(k)==iminus(1)
rat2=linspace(ratio,2*ratio,3);
    Cplx(k,:)=[exp(1i*(3/2*pi-rat2(1)))*C(k,:),exp(1i*(3/2*pi-rat2(2)))*C(k,:),exp(1i*(3/2*pi-rat2(3)))*C(k,:)];
    
end
end

%rat2=linspace(ratio,2*ratio,5);
%rat2=linspace(ratio,pi/2,3);

%for ii=1:n(1)
    
    %%%%%%%% DEFINITION OF FACTORS
%    factor=ones(6,1);
 %   mu(:)=Huge(ii,:);
  %  lamb(:)=V1(iplus2(:),iplus2(:));
 %   for i3=1 : 6
 %       factor(i3)=(lamb(i3)-mu(i3))/mu(i3);
 %       for i4=1:5
            
 %       factor(i3)=factor(i3)*
        
        
    
 %   end
%end

rat2=linspace(ratio,2*ratio,3);

 ang=angle(V1(iplus(1),iplus(1)));
   av=abs(V1(iplus(1),iplus(1)));
%    box1=real(V1(iplus2(k),iplus2(k)));
 %   box2=-imag(V1(iplus2(k),iplus2(k)));
  %  boxm=-min(-box1,-box2);
%C(k,:)=linspace(-boxm,-boxm+10^8,Csz); 
C(1,:)=linspace(av-3*10^7,av+3*10^7,Csz); 
 
   % C(1,:)=C(iplus(1),:);


for r=1:Csz %10^8:5*10^7:1*10^9
%for s=1:Csz%1*10^9:5*10^7:6*10^9    
for l1=1:3
%    for l2=1:5
m=zeros(num,1);



%m1=C(1,r)-rat(l1)*C(1,r)*1i; 
%m2=C(1,r)+rat(l1)*C(1,r)*1i;  

m1=exp(1i*(rat2(l1)+pi/2))*(C(1,r)); 
m2=exp(1i*3*pi/2-1i*rat2(l1))*(C(1,r));


%%%%%%mögliche Erweiterung%%%%%%%%%%%%%%%%%%%%%
%%% ganzen uv vector benutzen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Verbesserter Faktor

fac=ones(num,1);

m=zeros(num,1);
m(1)=m1;
%m(2)=m3;


for i1=1:num

 for i2=1:num
 if i2==i1
     fac(i1)=fac(i1)/m(i1)*V1(iminus(i1),iminus(i1))*(V1(iplus(i1),iplus(i1))-m(i1))*(V1(iplus(i1),iplus(i1))-conj(m(i1)))/((V1(iplus(i1),iplus(i1))-V1(iminus(i1),iminus(i1)))*conj(m(i1)));
     
 
 end
 if i2~=i1
     fac(i1)=fac(i1)/m(i2)*V1(iminus(i2),iminus(i2))*V1(iplus(i2),iplus(i2))*(V1(iplus(i1),iplus(i1))-m(i2))*(V1(iplus(i1),iplus(i1))-conj(m(i2)))/((V1(iplus(i1),iplus(i1))-V1(iplus(i2),iplus(i2)))*(V1(iplus(i1),iplus(i1))-V1(iminus(i2),iminus(i2)))*conj(m(i2)));   
 end

 end

end

%%%Sensitiven Faktoren mit dazu nehmen
%for i1=1:num
 %   for l=1:sz
    
  %      if l~=iplus(1:num)
            
   %         if l~=iminus(1:num)
          
    %        fac(i1)=fac(i1)*V1(l,l)/(V1(l,l)+V1(l,l)*step*uinf(Iinv(l))*vinf(Iinv(l)))*(V1(iplus(i1),iplus(i1))-V1(l,l)-step*V1(l,l)*uinf(Iinv(l))*vinf(Iinv(l)))/(V1(iplus(i1),iplus(i1))-V1(l,l));
            
     %       end
            
      %  end
        
    %end
%end


%fac=1;
%for i1=1:num   %%%%V1 überall ersetzen durch EW, weil wir die Indexmenge vereinheitlichen müssen
%for l=1:sz
%if abs(V1(l,l))<10^(16) %%%EW bei oo ausschließen
 %   if sens(Iinv(l))>0  %%%für alle sensitiven Eigenwerte nehme ich die Faktoren in das Produkt auf
 %   if l==iplus(i1)
 %   fac(i1)=fac(i1)/m(i1)*V1(iminus(i1),iminus(i1))*(V1(l,l)-m(i1))*(V1(l,l)-conj(m(i1)))/((V1(l,l)-V1(iminus(i1),iminus(i1)))*conj(m(i1)));
    %%%Hier wird der zweite sensitive Faktor auch berücksichtigt (oder allgemeiner alle sensitiven EW)
 %   end
  %      if l~=iplus(1:num)
            
   %         if l~=iminus(1:num)
          
    %        fac(i1)=fac(i1)*V1(l,l)/(V1(l,l)+V1(l,l)*step*uinf(Iinv(l))*vinf(Iinv(l)))*(V1(iplus(i1),iplus(i1))-V1(l,l)-step*V1(l,l)*uinf(Iinv(l))*vinf(Iinv(l)))/(V1(iplus(i1),iplus(i1))-V1(l,l));
            
     %       end
            
      %  end
      
   % end
    
  %  for i2=1:num    %%%setze die nicht sensitiven, die ich aber vorgebe
  %  fac(i1)=fac(i1)/m(i2)*V1(iminus(i2),iminus(i2))*V1(iplus(i2),iplus(i2))*(V1(iplus(i1),iplus(i1))-m(i2))*(V1(iplus(i1),iplus(i1))-conj(m(i2)))/((V1(iplus(i1),iplus(i1))-V1(iplus(i2),iplus(i2)))*(V1(iplus(i1),iplus(i1))-V1(iminus(i2),iminus(i2)))*conj(m(i2)));    
    
  %  end     
    
%end
%end
%end

%z=zeros(num,1);
%for u=1 :2
z1=fac;%V1(iminus,iminus)*(V1(iplus,iplus)-m1)*(V1(iplus,iplus)-m2)/(m1*m2*(V1(iplus,iplus)-V1(iminus,iminus)));
z2=conj(z1);
%end

 q=[z1;z2];
 %%%TASK: ADD MORE EQUATION, BECAUSE THE POLES WILL ENTER THE RIGHT
 %%%HALFPLANE! AT LEST ZERO EQUATIONS...
            uv=[diag(uinf(invipl(iplus))*transpose(vinf(invipl(iplus))));diag(uinf(invimin(iminus))*transpose(vinf(invimin(iminus))))];    
        
        %t=0.02;
         t=norm(q,2)^2/(uv'*q); %(uv'*q)/(norm(uv,2)^2);
          g=norm(t*uv-q,2);
          if abs(g)<gmin
              if t<0  %vorher >0
              if abs(t)<cmax  %vorher  <0.0001 %%muss hier gewisse größe zulassen, weil sonst nur kleine Veränderungen zugelassen werden!
          if abs(t)>10^(-16) %%%%%%%%%%%mal testen%%%%%%%%%%%%%%%%%%%%%%%%%
          %gmin=g;
          tmin=t;
          gmin=abs(g);
          i;
          j;
          mmin=m1;
         % mmin2=m3;
%          tmin=t;
          end
                 end
              end
          end
          
%end
          
          
        end
              
end 
%end
%end
       % end  
                  
       
        if gmin > 0%0.0001%0.0001    %%% gmin darf nicht zu klein sein
          if gmin < err%0.1%10%0.3 %10^7 %0.8 damit komplexe paare sichtbar sonst kleiner 
                %    ub=sqrt(0.01)*Sinv*uinf;
                    %eig(A,E+ub*ub')
                    
                Iout(Ioutcnt,1)=i;
                   Iout(Ioutcnt,2)=j;
                   Iout(Ioutcnt,3)=tmin;
       
                   Ioutcnt=Ioutcnt+1;
                   
                %%%%%%%%%NEU: CHECKEN WELCHER WIRKLICH VERBESSERT
                
                 b=zeros(sz,1);
        b(i)=1; b(j)=-1;
        Enew=Et+tmin*b*b';
       
       Eval=eig(At,Enew);
       
       max=0;
                for i4=1 : sz
                  % if imag(Eval(i4))~=0
                      
                       if abs(imag(Eval(i4))/real(Eval(i4)))>max
                      
                           max=abs(imag(Eval(i4))/real(Eval(i4)));
                       end
                  % end 
                end
                
                  
                   
                       
  %%%%%AUSGABE                  
      %    disp(['g: ',  num2str(gmin), ' t: ', num2str(tmin),' m: ', num2str(mmin,'%e'), ' (i,j): ', num2str(i), ' ',num2str(j),' max: ', num2str(max)])          
                   %%%%%%%%%%%% neues Abbruchkriterium das verhältnis aus
                   %%%%%%%%%%%%  Abweichung des Gleichungssystem zur
                   %%%%%%%%%%%%  benötigten Schrittweite darf nicht zu groß
                   %%%%%%%%%%%%  sein
               if max<gloquo %gloquo
                   
                  % Iout(Ioutcnt,1)=i;
                  % Iout(Ioutcnt,2)=j;
                  % Iout(Ioutcnt,3)=tmin;
                   gloquo=max;%abs(gmin/tmin);
                   tglo=tmin; %0                
                   iglo=i;
                   jglo=j;
                   gminglo=gmin;
                  % Ioutcnt=Ioutcnt+1;
                    %%realteil darf nicht zu groß sein sonst kleines gmin
                    %%obwohl EW weit abweichen
                    %if real(mmin)>-5*10^(10)%(-1)*imag(mmin)  
                     
                     %   Irel(1,glcnt)=i; Irel(2,glcnt)=j;
                      %  glcnt=glcnt+1;
                        % disp(['g: ',  num2str(gmin), ' t: ', num2str(tmin),' m: ', num2str(mmin,'%e'), ' m2: ', num2str(mmin2,'%e'), ' (i,j): ', num2str(i), ' ',num2str(j)]) 
%                    gmin
 %                   tmin  %18.10e\n
  %                  mmin
   %                 i
    %                j
                    %end
               end
          end
        end    
   end
              %    end
%
                  
                  
  %       end
   end
   end


if iglo+jglo>0
b=zeros(sz,1); b(iglo)=1; b(jglo)=-1;
Eo=Et+tglo*(b*b');
end
if iglo+jglo==0
Eo=Et;
tglo=0; %new
end

disp(['err: ',  num2str(gminglo), ' c_ij: ', num2str(tglo), ' (i,j): ', num2str(iglo), ' ',num2str(jglo),' Im/Re: ', num2str(gloquo)])          
      
end %of function
