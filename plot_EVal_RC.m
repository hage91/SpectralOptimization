%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

%plot(eig(A_RC,E_RC), '*');
%hold on ;

plot(eig(A_RC,E_RC), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[0,51/255,102/255],'MarkerFaceColor',[0,51/255,102/255]);
hold on;

%% Algorithm with prescribed angle 2
%Capacitor between 1 and 3 with cij=1.421107e-12.

%u=zeros(11,1); u(1)=1; u(3)=-1;

plot(eig(A_RC,E_RC_SO_2), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[1,102/255,0],'MarkerFaceColor',[1,102/255,0]);


%% Algorithm with prescribed angle 3/4 pi

plot(eig(A_RC,E_RC_SO_pi4), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[0,102/255,102/255],'MarkerFaceColor',[0,102/255,102/255]);

%um=zeros(12,1); um(6)=1; um(9)=-1;

%plot(eig(A_RC,E_RC-0.7*10^(-12)*(um*um')), 's','LineWidth',2,...
 %   'MarkerSize',10,'MarkerEdgeColor',[0.7,0,0.7],'MarkerFaceColor',[0.7,0,0.7]);


title('~~~~~Eigenvalues of precompensated two-stage OpAmp', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Re $z$', 'FontSize',18,'Interpreter', 'latex');
ylabel('Im $z$','FontSize',18,'Interpreter', 'latex');
legend({'OpAmp as in Figure 6.4','Algorithm 1 with $\varphi_*=2-\frac{\pi}{2}$','Algorithm 1 with $\varphi_*=\frac{\pi}{4}$'},'FontSize',18, 'Interpreter', 'latex');

xlim([ -5*10^8 0]);

ylim([ -2.5*10^8 2.5*10^8]);

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

grid on;
