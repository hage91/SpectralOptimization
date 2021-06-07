%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

cim=zeros(1,12); cim(6)=1;  %% in and output vectors
bim=zeros(12,1); bim(12)=1;

bodeplot(E_RC,A_RC, bim,cim,[0,51/255,102/255]);
hold on;

%% Algorithm with prescribed angle 2
%%%%old: Capacitor between 1 and 3 with cij=1.421107e-12.
%Capacitor between 3 and 4 with cij=2.477798e-13.
%Capacitor between 4 and 9 with cij=2.301542e-13.

%u=zeros(11,1); u(1)=1; u(3)=-1;
bodeplot(E_RC_SO_2,A_RC, bim,cim,[1,102/255,0]);


%% Algorithm with prescribed angle 3/4 pi
%um=zeros(12,1); um(1)=1; um(9)=-1;
%Capacitor between 1 and 9 with cij=9.986652e-13.
%Capacitor between 4 and 9 with cij=2.736976e-13.

bodeplot(E_RC_SO_pi4,A_RC, bim,cim,[0,102/255,102/255]);

%% Comparison with Miller, when C_M is increased
%um=zeros(12,1); um(6)=1; um(9)=-1;

%bodeplot(E_RC-1.2*10^(-12)*(um*um'),A_RC, bim,cim,[0.7,0,0.7]);


ylim([-3 5]);
xlim([10^6 2*10^8])

title('Bodeplot of precompensated two-stage OpAmp', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Frequency $\omega$', 'FontSize',18,'Interpreter', 'latex');
ylabel('$20\log_{10}|H(i\omega)|$','FontSize',18,'Interpreter', 'latex');
legend({'OpAmp as in Figure 6.4','Algorithm 1 with $\varphi_*=2-\frac{\pi}{2}$','Algorithm 1 with $\varphi_*=\frac{\pi}{4}$', 'Miller compensation with $C_M=2.2${\rm pF}'},'FontSize',18, 'Interpreter', 'latex');


a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)


