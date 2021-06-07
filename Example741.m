function [Eo,iglo,jglo,tglo,gloquo,gminglo]=Example741(Et,At,ang_ratio,domsens,eps,cmax,num,ity)

%%%%ratio ist Winkel von imaginärer Achse zum Sektor
%%%%domsens die kritische sensitivität

%num=2; %%% Zahl der Sensitiven Eigenwerte


[~,sz]=size(Et);
%ity=10; %%%Zahl der Eigenwerte bei oo

%num=2; %%% Zahl der Sensitiven Eigenwerte

gloquo=pi/2; %%quotient aus Fehler und Schrittweite 
gminglo=1;
tglo=0;

%isens=0;
icond=0;

iglo=0;
jglo=0;
tmin=-cmax;

[S,T,iplus,iminus,invipl,invimin,V1] = get_WCF(Et,At,ity,num);
%% hier nochmal vereinfachen, wir brauchen nicht alle ausgaben...
%sensitivität auch direkt?


%% Run over 1..11 nodes
for i=1:11
   for j=i:11
       v=zeros(sz,1);
        v(j)=-1; v(i)=1; 
        uinf=S*(At\v);
        vinf=transpose(T)*v;
        gmin=10; %tmin=0; 
         Ute=vinf*uinf';
%Sensitivity check
if abs(Ute(invipl(iplus(1:num)),invipl(iplus(1:num))))>=domsens   


%% Domain for parameter search for mu_i
%% hier müsste ratio mit eingehen...
Csz=50;  
%num2=num; %%Number of prescribed EVal
C=zeros(num,Csz);
Cplx=zeros(num,Csz,5); %%ggf erweitern
%angnew=zeros(3);


for k=1:num    
   ang=angle(V1(iplus(k),iplus(k)));
   angnew2=linspace(ang_ratio,3/2*pi,4); %%nicht nach oben suchen
   %%klassische suchtrategie
   angnew3=linspace((ang_ratio+ang)/2,ang_ratio,5);
   av=abs(V1(iplus(k),iplus(k)));  
   Vneu=av*exp(1i*(ang_ratio+ang)/2);  %sonst ang
   %chose angle as the mean of ratio and actual value
%C(k,:)=linspace(av-0.1*av,av+2*av,Csz); 
% C(k,:)=linspace(0,10*av,Csz); 
C(k,:)=linspace(av-0.5*av,10*av,Csz); 
%%an sensitivität anpassen l(0)+l'(0)*c_max               
  %  Cplx(k,:,1)=Vneu+exp(1i*angnew2(1))*C(k,:);
   % Cplx(k,:,2)=Vneu+exp(1i*angnew2(2))*C(k,:);
    %Cplx(k,:,3)=Vneu+exp(1i*angnew2(3))*C(k,:);
    %Cplx(k,:,4)=Vneu+exp(1i*angnew2(4))*C(k,:);
    Cplx(k,:,1)=exp(1i*angnew3(2))*C(k,:);
    Cplx(k,:,2)=exp(1i*angnew3(3))*C(k,:);
    Cplx(k,:,3)=exp(1i*angnew3(4))*C(k,:);
    Cplx(k,:,4)=exp(1i*angnew3(5))*C(k,:);
  
    
   % Cplx(k,:,5)=V1(iplus(k),iplus(k))+exp(1i*angnew2(5))*C(k,:);
%If imaginary part negative, then set zero...
end

%% compute the RHS

for r=1:Csz 
for l1=1:5

fac=ones(num,1);


m=zeros(num,1);
for l3=1 :num
m(l3)=Cplx(l3,r,l1);
end

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



z1=fac;
z2=conj(z1);

q=[z1;z2];
uv=[diag(uinf(invipl(iplus))*transpose(vinf(invipl(iplus))));diag(uinf(invimin(iminus))*transpose(vinf(invimin(iminus))))];    
        
t=norm(q,2)^2/(uv'*q);%(uv'*q)/(norm(uv,2)^2)
          g=norm(t*uv-q,2);
          %% nochmal überdenken, ob man hier nicht nur eien schwelle definiert
  %%schwelle für g_min und stattdessen kleinste Kapazität  
          if abs(g)<eps && real(t)>tmin && real(t)<0 && imag(t)<10^(-13)
             %the perturbations are negative, because the matrices have opposite sign    
            if real(t)< -cmax && real(t)> -1.1*cmax
            tmin=-cmax; 
            gmin=abs(g/real(t)*cmax);
            end   
            if real(t)>= -cmax
            tmin=real(t); 
            gmin=abs(g);      
            end
          end
          
          
          
        end
              
end 
                  
 %% Verification step: We check if the eigenvalues are inside 
 % the prescribed sector: result: if in sector then keep otherwise 
 % remove; set tmin to the old tmin
 
 
        if gmin > 0 && gmin< eps % I think this condition is obsolete
                
        b=zeros(sz,1);
        b(i)=1; b(j)=-1;
       
       Eval=eig(At,Et+tmin*(b*b'));
       
       max=pi;  %%komisch; funktioniert aber nicht immer
                for i4=1 : sz
                    
                   
                      if angle(Eval(i4))<max && imag(Eval(i4))>0
                      
                           max=angle(Eval(i4));%abs(imag(Eval(i4))/real(Eval(i4)));
                       
                      end 
                end
                %% problem: hier wird schon das beste im Sektor gewählt
                  icond=icond+1;
               if max>gloquo && tmin<0 %% müsste hier dann auch bandbreite überprüfen
                   
                   if tglo == 0
                  gloquo=max; %abs(gmin/tmin);
                   tglo=tmin;                 
                   iglo=i;
                   jglo=j;
                   gminglo=gmin;
                   end
                   
                   if abs(tglo)>abs(tmin)
                           gloquo=max; %abs(gmin/tmin);
                   tglo=tmin;                 
                   iglo=i;
                   jglo=j;
                   gminglo=gmin;
                       
                   end
               end
        end    

end
end
end

 

%% output
if iglo+jglo>0
b=zeros(sz,1); b(iglo)=1; b(jglo)=-1;
Eo=Et+tglo*(b*b');
end
if iglo+jglo==0
Eo=Et; iglo=1; jglo=1; gloquo=ang_ratio;
%gloquo=(ang_ratio+angle(V1(iplus(1),iplus(1))))/2;
end

 disp(['err: ',  num2str(gminglo), ' c_ij: ', num2str(-tglo), ' (i,j): ', num2str(iglo), ' ',num2str(jglo),' Re/Im: ', num2str(gloquo)])          

end 
