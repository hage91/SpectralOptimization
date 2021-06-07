



%%%Winkelwertefür die Iterationen kann man automatisch machen und an gloqou anpassen...dann den mit ausgeben

%phi=[pi/6*ones(1,10),pi/4*ones(1,5),pi/3*ones(1,10), pi/2*ones(1,15)];
%phi=[pi/6*ones(1,2),pi/4*ones(1,1),pi/3*ones(1,2), pi/3,pi/3 pi/2*ones(1,6)];
phi=[pi/8,pi/6,pi/4*ones(1,1),pi/3*ones(1,2), pi/3,pi/3 pi/2*ones(1,6)];


%phi=[pi/3.5*ones(1,5),pi/3*ones(1,6),pi/(2.2)*ones(1,20)];
%phi=[pi/2*ones(1,20),pi/(2)*ones(1,20)];
%domsens=[10^11*ones(1,12),5*10^10*ones(1,18),5*10^10*ones(1,20)];
domsens=1*10^10*ones(1,40);


%%Strategie: automatisierte sensitivität: Falls man in einem Schritt gloquo
%%nicht wesentlich verändert oder keine Kapazität gefunden wird, dann muss
%%sensitivität verringert werden
I=zeros(27,27);
%Enew=E0;
gloquo=2;
gloquoprev=2;
tges=0;
E=zeros(59,59,40);
E(:,:,1)=E1;
%phi=pi/10;
h1=plot(eig(A1,E1), 'd');
%axis([-2*10^9 0 -10^9 10^9]) 
axis([-2*10^(9) 0 -2*10^9 2*10^9]) 
hold on;
steps=6;

for i=2:steps
   if gloquo> 1     %%%%%%%%%%%%%                          phi(i)
       gloquoprev=gloquo;
   [E(:,:,i),iglo,jglo,tglo,gloquo,gminglo]=Example741(E(:,:,i-1),A1,phi(i-1),domsens(i),2,2.5*10^(-12));
   I(iglo,jglo)=I(iglo,jglo)+tglo;
%disp(['Step', num2str(i), ' err: ',  num2str(gminglo), ' c_ij: ', num2str(tglo), ' (i,j): ', num2str(iglo), ' ',num2str(jglo),' Im/Re: ', num2str(gloquo)])          
   tges=tges-tglo;
   d=zeros(59,1); d(iglo)=1; d(jglo)=-1; %c1=tglo/4; c2=tglo/2; c3=3*tglo/4;
 %  plot(eig(A1,E(:,:,i)-c1*(b*b')), 'x', 'Color', [i/40 0 (40-i)/40]);
  % plot(eig(A1,E(:,:,i)-c2*(b*b')), 'x', 'Color', [i/40 0 (40-i)/40]);
   %plot(eig(A1,E(:,:,i)-c3*(b*b')), 'x', 'Color', [i/40 0 (40-i)/40]);
   plot(eig(A1,E(:,:,i)-tges*(d*d')), 'x', 'Color', [(i-1)/steps 0 (steps-(i-1))/steps]);
   end
   %phi=pi/2;
   %if pi+atan(-gloquo)-0.1*abs((gloquoprev-gloquo))/gloquoprev-0.1< pi/2
   %phi=pi+atan(-gloquo)-0.1*abs((gloquoprev-gloquo))/gloquoprev-0.1; %%etwas kleinerer Winkel als vorher, daher subtrahieren wir willkürlich 0.1
   %end
end

%Ausgabe der Kapazitäten mit Knotenübersetzung (da Variablen permutiert)
p=[1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,26,27,3,4,5,6,7,8,9];
for i=1:27
for j=i:27
    if I(i,j)~=0
    fprintf('Kapazität zwischen %d und %d mit cij=%d.\n',p(i),p(j),-I(i,j));
    end
end
end

%%%Generate a Figure
%h2=plot(eig(A1,EC12),'p');
h3=plot(eig(A1,E(:,:,5)), 'o');
%h4=plot(eig(A1,E741old), 'h');
legend([h1,h3], {'Unperturbed', 'New Method'})
