
%% all input parameters
domsens=5*10^10*ones(1,40);
node_list=[1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,26,27,3,4,5,6,7,8,9];
eps=0.5;  
num=1;
ity=10; % eigenavlues at infty (just dimension of ker E)
c_max=2*10^(-12);
E1=E741;
A1=A741;
%c_min=0.7*10^(-12);
cs=linspace(c_max,c_max,1);
peak_max=5;
c=zeros(1,59); c(19)=1;  %%%entspricht Knoten 26 im Bild
b=zeros(59,1); b(59)=1;  %%%entspricht Knoten 5 im Bild; das ist irgendwie komisch und passt nicht zu meiner nummerierung
steps=15; 
glq=10;
ang_final=1.9;
C_glo=1;
bd_glo=1;
%phi_max=pi/15;
phi_min=5/6*pi;%pi*2.5/6; %% should be slightly below glo_final...
phi_ini=linspace(1.75,1.9,3);
var=linspace(1/10,1/3,4);

%% Optimization 
for iv=1:length(var)
for ic=1:length(cs)

    for i=1:length(phi_ini)

    [gloquo,tges,I,Eg] = onepara(E1,A1,phi_ini(i),phi_min,domsens,eps,cs(ic),ang_final,steps,num,ity,var(iv));
   % gloquo
      [bd,pk]=get_bandwidth(Eg,A1,b,c);
    
      if tges<C_glo 
      %  if bd>bd_glo
            if pk<peak_max && gloquo>ang_final && abs(tges)>0
     % if gloquo<glo_final
    glq=gloquo;
          C_glo=tges;
          bd_glo=bd;
      Eglo=Eg;
      Iglo=I;
      pk_glo=pk;
            end
            end
      %  end

    end
end
end

%% Output
%Bodeplot
bodeplot(E1,A1,b,c,[0,51/255,102/255],'-');
hold on;

bodeplot(Eglo,A1,b,c,[1,102/255,0], '-');

%List of Capacities
Ifull=full(Iglo);
for i=1:27
for j=i:27
    if Ifull(i,j)~=0
    fprintf('Capacitor between %d and %d with cij=%d.\n',node_list(i),node_list(j),-Ifull(i,j));
    end
end
end

