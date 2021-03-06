function createFigure(YMatrix1, YMatrix2)
%CREATEFIGURE(YMATRIX1, YMATRIX2)
%  YMATRIX1:  matrix of y data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 16-Oct-2017 13:33:34

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','YMinorTick','on',...
    'YScale','log',...
    'FontSize',12,...
    'Position',[0.0878378378378378 0.119966777408638 0.414864864864865 0.850132890365448]);
box(axes1,'on');
hold(axes1,'on');

% Create ylabel
ylabel('$\frac{F(y^k)-F^{\star}}{\vert F^{\star}\vert}$-in logscale',...
    'HorizontalAlignment','center',...
    'FontSize',15,...
    'Interpreter','latex');

% Create xlabel
xlabel('Iterations','HorizontalAlignment','center','FontSize',16,...
    'Interpreter','latex');

% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(YMatrix1,'Parent',axes1,'MarkerSize',3,'LineWidth',3);
set(semilogy1(1),'DisplayName','scvx-PAPA-1.0','LineStyle','-.','Color',[0 0 1]);
set(semilogy1(2),'DisplayName','scvx-PAPA-0.1','LineStyle','--','Color',[0 0 0]);
set(semilogy1(3),'DisplayName','scvx-PAPA-0.01','LineStyle',':','Color',[1 0 1]);
set(semilogy1(4),'DisplayName','scvx-PAPA-rs-1.0','LineStyle','-.','Color',[1 0 0]);
set(semilogy1(5),'DisplayName','scvx-PAPA-rs-0.1','LineStyle','--',...
    'Color',[0.0784313753247261 0.168627455830574 0.549019634723663]);
set(semilogy1(6),'DisplayName','scvx-PAPA-rs-0.01','LineStyle','-',...
    'Color',[0 1 0]);
set(semilogy1(7),'DisplayName','Chambolle-Pock-1.0','LineStyle','--',...
    'Color',[0.494117647409439 0.184313729405403 0.556862771511078]);
set(semilogy1(8),'DisplayName','Chambolle-Pock-0.1', 'LineStyle','-.',...
    'Color',[0.600000023841858 0.200000002980232 0]);
set(semilogy1(9),'DisplayName','Chambolle-Pock-0.01', 'LineStyle','--',...
    'Color',[0.23137255012989 0.443137258291245 0.337254911661148]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.303865159139948 0.26578645896969 0.175275169557903 0.389801702148178],...
    'Interpreter','latex');

% Create axes
axes2 = axes('Parent',figure1,'YGrid','on','XGrid','on','YMinorTick','on',...
    'YScale','log',...
    'FontSize',12,...
    'Position',[0.585135135135135 0.119966777408638 0.395945945945946 0.846810631229236]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes2,[1e-16 0.1]);
box(axes2,'on');
hold(axes2,'on');

% Create ylabel
ylabel('$\frac{F(y^k)-F^{\star}}{\vert F^{\star}\vert}$-in logscale',...
    'HorizontalAlignment','center',...
    'FontSize',15,...
    'Interpreter','latex');

% Create xlabel
xlabel('Iterations','HorizontalAlignment','center','FontSize',16,...
    'Interpreter','latex');

% Create multiple lines using matrix input to semilogy
semilogy2 = semilogy(YMatrix2,'Parent',axes2,'MarkerSize',3,'LineWidth',3);
set(semilogy2(1),'DisplayName','scvx-PAPA-1.0','LineStyle','-.','Color',[0 0 1]);
set(semilogy2(2),'DisplayName','scvx-PAPA-0.1','LineStyle','--','Color',[0 0 0]);
set(semilogy2(3),'DisplayName','scvx-PAPA-0.01','LineStyle',':','Color',[1 0 1]);
set(semilogy2(4),'DisplayName','scvx-PAPA-rs-1.0','LineStyle','-.','Color',[1 0 0]);
set(semilogy2(5),'DisplayName','scvx-PAPA-rs-0.1','LineStyle','--',...
    'Color',[0.0784313753247261 0.168627455830574 0.549019634723663]);
set(semilogy2(6),'DisplayName','scvx-PAPA-rs-0.01','LineStyle','-',...
    'Color',[0 1 0]);
set(semilogy2(7),'DisplayName','Chambolle-Pock-1.0','LineStyle','--',...
    'Color',[0.494117647409439 0.184313729405403 0.556862771511078]);
set(semilogy2(8),'DisplayName','Chambolle-Pock-0.1', 'LineStyle','-.',...
    'Color',[0.600000023841858 0.200000002980232 0]);
set(semilogy2(9),'DisplayName','Chambolle-Pock-0.01', 'LineStyle','--',...
    'Color',[0.23137255012989 0.443137258291245 0.337254911661148]);
