%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

cim=zeros(1,11); cim(6)=1;  %% in and output vectors
bim=zeros(11,1); bim(11)=1;

bodeplot(E_noRC,A_noRC, bim,cim,[0,51/255,102/255]);
hold on;

%% Algorithm with prescribed angle 2
%Capacitor between 1 and 3 with cij=1.421107e-12.

u=zeros(11,1); u(1)=1; u(3)=-1;
bodeplot(E_noRC_SO_2,A_noRC, bim,cim,[1,102/255,0]);
%bodeplot(E_noRC-1.8*10^(-12)*u*u',A_noRC, bim,cim,[1,102/255,0]);


%% Algorithm with prescribed angle 3/4 pi

bodeplot(E_noRC_SO_pi4,A_noRC, bim,cim,[0,102/255,102/255]);

ylim([-3 12]);
xlim([10^6 2*10^8])

title('Bodeplot of two-stage OpAmp without precompensation', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Frequency $\omega$', 'FontSize',18,'Interpreter', 'latex');
ylabel('$20\log_{10}|H(i\omega)|$','FontSize',18,'Interpreter', 'latex');
legend({'OpAmp as in Figure 6.5','Algorithm 1 with $\varphi_*=2-\frac{\pi}{2}$','Algorithm 1 with $\varphi_*=\frac{\pi}{4}$'},'FontSize',18, 'Interpreter', 'latex');


a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)
