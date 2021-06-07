%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

%plot(eig(A_RC,E_RC), '*');
%hold on ;

plot(eig(A_noRC,E_noRC), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[0,51/255,102/255],'MarkerFaceColor',[0,51/255,102/255]);
hold on;

%% Algorithm with prescribed angle 2
%Capacitor between 1 and 3 with cij=1.421107e-12.

%u=zeros(11,1); u(1)=1; u(3)=-1;

plot(eig(A_noRC,E_noRC_SO_2), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[1,102/255,0],'MarkerFaceColor',[1,102/255,0]);


%% Algorithm with prescribed angle 3/4 pi

plot(eig(A_noRC,E_noRC_SO_pi4), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[0,102/255,102/255],'MarkerFaceColor',[0,102/255,102/255]);


title('~~~~Eigenvalues of two stage CMOS without precomp.', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Re $z$', 'FontSize',18,'Interpreter', 'latex');
ylabel('Im $z$','FontSize',18,'Interpreter', 'latex');
legend({'OpAmp as in Figure 6.5','Algorithm 1 with $\varphi_*=2-\frac{\pi}{2}$','Algorithm 1 with $\varphi_*=\frac{\pi}{4}$'},'FontSize',18, 'Interpreter', 'latex');

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

grid on;
