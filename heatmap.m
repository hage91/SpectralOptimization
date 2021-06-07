
C=linspace(-3*10^(-12),3*10^(-12),50);
R=linspace(-10^4,0.1,50);

%%R=C=0 ist das Sommer beispiel;
%%Für R=-10^4 wird der Widerstand entfernt
%% Für C=-10^(-12) die Kapazität entfernt

cm=zeros(1,12); cm(6)=1;  %% in and output vectors
bm=zeros(12,1); bm(12)=1;



uC=zeros(12,1); uC(6)=1;uC(9)=-1; %%zwischen 6 und 9
UR=zeros(12,1); uR(3)=1; uR(9)=-1;   %% zwischen 3 und 9

bd_max=zeros(length(C),length(R));
C_max=0;
R_max=0;

for i=1:length(C)
    for j=1:length(R)
    
        En=E1-C(i)*(uC*uC');
        An=A1+R(j)^(-1)*(uR*uR');
        
        f=0;
        
        z=eig(An,En);
        for k=1:length(z)
        if real(z(k))>0  && real(z(k))<10^16
        
            f=1;
           
            
        end
        end
        if f==0
            [b,p]=get_bandwidth(En,An,bm,cm);
           
            
            if p<10 && b>10^7%bd_max

                
                R_max=R(j);
                C_max=C(i);
                bd_max(i,j)=b;
%          plot(C(i),R(j), '*');
           % hold on;
            end
            
        end
    end
          
       
            
          
end
%    figure
 %        mesh(bd_max);
  contour(bd_max,16);
colormap default;    % change color map  