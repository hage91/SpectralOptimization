%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

cim=zeros(1,12); cim(6)=1;  %% in and output vectors
bim=zeros(12,1); bim(12)=1;

u2=zeros(12,1); u2(6)=1; u2(9)=-1;
%% with Miller
bodeplot(E_RC-2.27*10^(-12)*(u2*u2'),A_RC, bim,cim,[0,51/255,102/255]);
hold on;




%% post opt
u3=zeros(12,1); u3(1)=1; u3(9)=-1;
bodeplot(E_RC+0.9*10^(-12)*(u2*u2')-1.26*10^(-12)*(u3*u3'),A_RC, bim,cim,[1,102/255,0]);

%% without Miller
u1=zeros(11,1); u1(1)=1; u1(3)=-1;
cim=zeros(1,11); cim(6)=1;  %% in and output vectors
bim=zeros(11,1); bim(11)=1;

bodeplot(E_noRC-2.07*10^(-12)*(u1*u1'),A_noRC, bim,cim,[0,102/255,102/255]);

%% Comparison with Miller, when C_M is increased
%um=zeros(12,1); um(6)=1; um(9)=-1;

%bodeplot(E_RC-1.2*10^(-12)*(um*um'),A_RC, bim,cim,[0.7,0,0.7]);


ylim([-3 3]);
xlim([10^6 2*10^8])

title('Bodeplot of the compensated two-stage OpAmp', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Frequency $\omega$', 'FontSize',18,'Interpreter', 'latex');
ylabel('$20\log_{10}|H(i\omega)|$','FontSize',18,'Interpreter', 'latex');
%legend({'without compensation','Algorithm 1 with $\varphi_*=2-\frac{\pi}{2}$','Algorithm 1 with $\varphi_*=\frac{\pi}{4}$', 'Miller compensation with $C_M=2.2${\rm pF}'},'FontSize',14, 'Interpreter', 'latex');
legend({'$C_M=2.27$ pF','$C_M=0.1$ pF, $c_{1,9}=$1.26 pF','$c_{1,3}=2.07$ pF'},'FontSize',18, 'Interpreter', 'latex');



