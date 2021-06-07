%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

c=zeros(1,59); c(19)=1;  %%%entspricht Knoten 26 im Bild
b=zeros(59,1); b(59)=1; 

bodeplot(E741,A741, b,c,[0,51/255,102/255], '-');
hold on;

%% Algorithm with prescribed angle 2
%Capacitor between 1 and 3 with cij=1.421107e-12.

u1=zeros(59,1); u1(6)=1; u1(26)=-1;    
%u2=zeros(59,1); u2(16)=1; u2(27)=-1; 
u2=zeros(59,1); u2(10)=1; u2(16)=-1;
u3=zeros(59,1); u3(6)=1; u3(15)=-1;

%bodeplot(E1-1.17*10^(-12)*(u1*u1')-1.27*10^(-12)*(u2*u2')-1.18*10^(-12)*(u3*u3'),A1, b,c,[0,102/255,102/255]);
%bodeplot(E_noRC-1.8*10^(-12)*u*u',A_noRC, bim,cim,[1,102/255,0]);
%Capacitor between 14 and 8 with cij=1.815543e-12.
%Capacitor between 18 and 23 with cij=8.582694e-13. 
%Capacitor between 23 and 9 with cij=1.416473e-12.


%% Algorithm with prescribed angle 3/4 pi

bodeplot(E741-10.5*10^(-12)*(u1*u1')-0*10^(-12)*(u2*u2')-0*10^(-12)*(u3*u3'),A741, b,c,[1,102/255,0], '-');
%bodeplot(E1-30*10^(-12)*(u1*u1')-0*10^(-12)*(u2*u2')-0*10^(-12)*(u3*u3'),A1,b,c,[1,102/255,0],'-.');

%eig(A1,E1-20*10^(-12)*(u1*u1'))

u1=zeros(59,1); u1(9)=1; u1(10)=-1; % Miller   
%u1=zeros(59,1); u1(10)=1; u1(16)=-1; %%Miller Var with 2.38
u4=zeros(59,1); u4(1)=1; u4(21)=-1;
bodeplot(E741-2.6*10^(-12)*(u1*u1')-0*10^(-12)*(u4*u4'),A741, b,c,[0,102/255,102/255], '-');
%3.13          2.62                 0.3

ylim([-3 17]);
xlim([2*10^6 3*10^8])

title('Bodeplot of OpAmp $\mu A 741$', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Frequency $\omega$', 'FontSize',18,'Interpreter', 'latex');
ylabel('$20\log_{10}|H(i\omega)|$','FontSize',18,'Interpreter', 'latex');
legend({'OpAmp as in Figure 6.10','$c_{8,14}=10.5$ pF','Miller compensation $C_M=2.62$ pF, $c_{1,3}=0.3$ pF'},'FontSize',18, 'Interpreter', 'latex');


a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)


