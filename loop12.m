
%%%Winkelwertefür die Iterationen kann man automatisch machen und an gloqou anpassen...dann den mit ausgeben

%phi=[pi/6*ones(1,10),pi/4*ones(1,5),pi/3*ones(1,10), pi/2*ones(1,6)];
phi=[pi/7,pi/6,pi/5,pi/4,pi/3,pi/2,pi/2];
%phi=[pi/3.5*ones(1,5),pi/3*ones(1,6),pi/(2.2)*ones(1,20)];
%domsens=[10^11*ones(1,12),5*10^10*ones(1,18),5*10^10*ones(1,20)];
domsens=0.1*10^(10)*ones(1,31);


%%Strategie: automatisierte sensitivität: Falls man in einem Schritt gloquo
%%nicht wesentlich verändert oder keine Kapazität gefunden wird, dann muss
%%sensitivität verringert werden

I=zeros(12,12);
%Enew=E0;
gloquo=2;
%gloquoprev=2;
tges=0;
E=zeros(12,12,31);
E(:,:,1)=Em;
cnt=1;
%phi=pi/10;
h1=plot(eig(Am,Em), 'd');
axis([-5*10^8 0 -.1*10^9 .1*10^9]) 
hold on;

for i=2:31
   if gloquo>0.1      %%%%%%%%%%%%%                          phi(i)
       gloquoprev=gloquo;
   [E(:,:,i),iglo,jglo,tglo,gloquo,gminglo]=NewExample(E(:,:,i-1),Am,phi(i-1),domsens(i-1),1,10*10^(-13));
   I(iglo,jglo)=I(iglo,jglo)+tglo;
%disp(['Step', num2str(i), ' err: ',  num2str(gminglo), ' c_ij: ', num2str(tglo), ' (i,j): ', num2str(iglo), ' ',num2str(jglo),' Im/Re: ', num2str(gloquo)])          
   tges=tges-tglo;
   b=zeros(12,1); b(iglo)=1; b(jglo)=-1; c1=tglo/4; c2=tglo/2; c3=3*tglo/4;
  % plot(real(eig(Am,E(:,:,i)-c1*(b*b'))),imag(eig(Am,E(:,:,i)-c1*(b*b'))), 'x', 'Color', [i/10 0 (10-i)/10]);
   plot(real(eig(Am,E(:,:,i)-c2*(b*b'))),imag(eig(Am,E(:,:,i)-c2*(b*b'))), 'x', 'Color', [ i/10 0 (10-i)/10]);
  % plot(real(eig(Am,E(:,:,i)-c3*(b*b'))),imag(eig(Am,E(:,:,i)-c3*(b*b'))), 'x', 'Color', [i/10 0 (10-i)/10]);
   plot(real(eig(Am,E(:,:,i)-0*(b*b'))),imag(eig(Am,E(:,:,i)-tglo*(b*b'))), 'x', 'Color', [i/10 0 (10-i)/10]);
   cnt=cnt+1;
   end
   %phi=pi/2;
   %if pi+atan(-gloquo)-0.1*abs((gloquoprev-gloquo))/gloquoprev-0.1< pi/2
   %phi=pi+atan(-gloquo)-0.1*abs((gloquoprev-gloquo))/gloquoprev-0.1; %%etwas kleinerer Winkel als vorher, daher subtrahieren wir willkürlich 0.1
   %end
end

%Ausgabe der Kapazitäten mit Knotenübersetzung (da Variablen permutiert)
%p=[1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,26,27,3,4,5,6,7,8,9];
for i=1:12
for j=i:12
    if I(i,j)~=0
    fprintf('Kapazität zwischen %d und %d mit cij=%d.\n',i,j,-I(i,j));
    end
end
end

%%%Generate a Figure
%h2=plot(eig(Am,EC1),'o');
h3=plot(real(eig(Am,E(:,:,cnt))),imag(eig(Am,E(:,:,cnt))), 'd', 'Color', [0 1 0]);
%h4=plot(eig(Am,E741old), 'h');
legend([h1,h3], {'Unperturbed', 'New Method'})

