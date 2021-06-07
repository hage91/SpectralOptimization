%% TUI colors
%[0,51/255,102/255] blue
%[0,102/255,102/255]  turkey
%[1,102/255,0]     orange

%plot(eig(A_RC,E_RC), '*');
%hold on ;

u1=zeros(11,1); u1(1)=1; u1(3)=-1;


plot(eig(A_noRC,E_noRC-2.07*10^(-12)*(u1*u1')), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[0,102/255,102/255],'MarkerFaceColor',[0,102/255,102/255]);
hold on;

%% Algorithm with prescribed angle 2
%Capacitor between 1 and 3 with cij=1.421107e-12.

%u=zeros(11,1); u(1)=1; u(3)=-1;



%%u=zeros(12,1); u(1)=1; u(9)=-1;
%%u2=zeros(12,1); u2(6)=1; u2(9)=-1;
%E_RC-1.26*10^(-12)*u*u'+0.9*10^(-12)*u2*u2';
u2=zeros(12,1); u2(6)=1; u2(9)=-1;


plot(eig(A_RC,E_RC-2.27*10^(-12)*(u2*u2')), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[0,51/255,102/255],'MarkerFaceColor',[0,51/255,102/255]);



u3=zeros(12,1); u3(1)=1; u3(9)=-1;


plot(eig(A_RC,E_RC+0.9*10^(-12)*(u2*u2')-1.26*10^(-12)*(u3*u3')), 's','LineWidth',2,...
    'MarkerSize',10,'MarkerEdgeColor',[1,102/255,0],'MarkerFaceColor',[1,102/255,0]);

%um=zeros(12,1); um(6)=1; um(9)=-1;

%plot(eig(A_RC,E_RC-0.7*10^(-12)*(um*um')), 's','LineWidth',2,...
 %   'MarkerSize',10,'MarkerEdgeColor',[0.7,0,0.7],'MarkerFaceColor',[0.7,0,0.7]);


title('Eigenvalues of compensated two-stage OpAmp', 'FontSize',18, 'Interpreter', 'latex');
xlabel('Re $z$', 'FontSize',18,'Interpreter', 'latex');
ylabel('Im $z$','FontSize',18,'Interpreter', 'latex');
legend({'$C_M=0.1$ pF, $c_{1,9}=$1.26 pF','$C_M=2.27$ pF','$c_{1,3}=2.07$ pF'},'FontSize',18, 'Interpreter', 'latex');

xlim([ -8*10^8 0]);

ylim([ -2.5*10^8 2.5*10^8]);

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

grid on;

%% with decreasing C_M u=zeros(12,1); u(1)=1; u(9)=-1; peaking <= 0.1
%%u=zeros(12,1); u(1)=1; u(9)=-1;
%%u2=zeros(12,1); u2(6)=1; u2(9)=-1;
%E_RC-1.26*10^(-12)*u*u'+0.9*10^(-12)*u2*u2';

%% with decreasing C_M u=zeros(12,1); u(1)=1; u(9)=-1; peaking <= 0.5
%%u=zeros(12,1); u(1)=1; u(9)=-1;
%%u2=zeros(12,1); u2(6)=1; u2(9)=-1;
%E_RC-1.184*10^(-12)*u*u'+0.9*10^(-12)*u2*u2';



