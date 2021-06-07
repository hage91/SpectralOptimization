%% typical value domsens=1*10^10*ones(1,40);
%% eps= 2,  c_max=2.5*10^(-12)
function [ang_ratio,tges,I,Eglo] = onepara(E1,A1,phi_max,phi_min,domsens,eps,c_max,glo_final,steps,num,ity,var)

%phi=zeros(1,length(phi_up));
%phi=[pi/8,pi/6,pi/4*ones(1,1),pi/3*ones(1,2), pi/3,pi/3 pi/2*ones(1,6)];

I=sparse(zeros(27,27)); %% stores the i,j indizes (size does depend on the circuit...)
%Enew=E0;
%gloquo=phi_max;%tan(pi/2-phi_max);
ang_ratio=phi_max;
%ang_ratio_n=2*phi_max;
tges=0;
tglo=1;
n=length(E1);
E=zeros(n,n,40);
E(:,:,1)=E1;


%gmax=4;
for i=1:steps
   if ang_ratio<glo_final && abs(tglo)>0
       %% evetl noch über die Varianz loopen und dann einen schritt mehrfach machen
     % if abs(tglo)>0
      ang_ratio_h=min(ang_ratio+0.001+var*abs(rand(1)),phi_min);
     % end
      % end
       if i==1
           ang_ratio_h=phi_max;
       end
        %  ang_ratio=ang_ratio_n;
%% hier dann einen varianz loop rum....
   [E(:,:,i+1),iglo,jglo,tglo,ang_ratio_h]=Example741(E(:,:,i),A1,ang_ratio_h,domsens(i),eps,c_max,num,ity);
   I(iglo,jglo)=I(iglo,jglo)+tglo;
   tges=tges-tglo;
   Eglo=E(:,:,i+1);
   if abs(tglo)>0
       ang_ratio=ang_ratio_h;
   end
   
   end
end


if tges==0
    Eglo=E1;
end
%ang_ratio=ang_ratio_n
end

