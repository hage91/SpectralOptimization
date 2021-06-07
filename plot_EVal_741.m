%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

%plot(eig(A_RC,E_RC), '*');
%hold on ;

plot(eig(A741,E741), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[0,51/255,102/255],'MarkerFaceColor',[0,51/255,102/255]);
hold on;

%% root plot
%sys = dss(A1,b,c,0,E1);
%plot(zero(sys), 'o','LineWidth',2,...
 %   'MarkerSize',10,'MarkerEdgeColor',[0,51/255,102/255]);



%% high bandwidth position

u1=zeros(59,1); u1(6)=1; u1(26)=-1;  
Evar1=E741-10.5*10^(-12)*(u1*u1');
plot(eig(A741,Evar1), 's','LineWidth',2,...
    'MarkerSize',6,'MarkerEdgeColor',[1,102/255,0],'MarkerFaceColor',[1,102/255,0]);
%% root plot
%sys = dss(A1,b,c,0,Evar1);
%plot(zero(sys), 'o','LineWidth',2,...
%    'MarkerSize',10,'MarkerEdgeColor',[0,102/255,102/255]);



%% Miller position
u1=zeros(59,1); u1(9)=1; u1(10)=-1; % Miller   
u4=zeros(59,1); u4(1)=1; u4(21)=-1; 

Evar2=E741-2.62*10^(-12)*(u1*u1')-0.33*10^(-12)*(u4*u4');
plot(eig(A741,Evar2), 's','LineWidth',2,...
    'MarkerSize',6,'MarkerEdgeColor',[0,102/255,102/255],'MarkerFaceColor',[0,102/255,102/255]);

%sys = dss(A1,b,c,0,Evar2);
%plot(, 'o', 'Color',[1,102/255,0]);
%plot(zero(sys), 'o','LineWidth',2,...
 %   'MarkerSize',10,'MarkerEdgeColor',[1,102/255,0]);


%plot(eig(A_noRC,E_noRC_SO_pi4), 's','LineWidth',2,...
 %   'MarkerSize',10,'MarkerEdgeColor',[0,102/255,102/255],'MarkerFaceColor',[0,102/255,102/255]);


title('Eigenvalues of OpAmp $\mu A741$', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Re $z$', 'FontSize',18,'Interpreter', 'latex');
ylabel('Im $z$','FontSize',18,'Interpreter', 'latex');
legend({'OpAmp as in Figure 6.10','$c_{8,14}=10.5$ pF','Miller compensation $C_M=2.62$ pF, $c_{1,3}=0.3$ pF'},'FontSize',18, 'Interpreter', 'latex');
xlim([-2*10^(8) 0]);

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

grid on;

