%% The first plot
figure (1)
hold on
plot(Size(:,1),'LineWidth', 3)
plot(Size(:,7),'LineWidth', 3)
plot(Size(:,2)+Size(:,3)+Size(:,4)+Size(:,5)+Size(:,6),'LineWidth', 3)
hold on
legend('Accumulative','Last 5-Frame','LCS filter','Interpreter','latex')
xlabel('Frame','FontSize',32,'LineWidth' , 3,'Interpreter','latex')
ylabel('Memory Dimensionality','FontSize',32,'LineWidth' , 3,'Interpreter','latex')
set(gca,'FontSize',32); 
%% The second plot
% figure (2)
% hold on
% %plot(Size(:,2),'LineWidth', 2)
% plot(Size(:,3),'LineWidth', 3)
% plot(Size(:,4),'LineWidth', 3)
% plot(Size(:,5),'LineWidth', 3)
% plot(Size(:,6),'LineWidth', 3)
% hold on
% legend('$\bf{E}_r$','$\bf{C}_n$','$\bf{C}_r$','\bf{S}','Interpreter','latex')
% xlabel('Frame','FontSize',32,'LineWidth' , 3,'Interpreter','latex')
% ylabel('Memory Dimensionality','FontSize',32,'LineWidth' , 3,'Interpreter','latex')
% set(gca,'FontSize',32); 
%% The second plot
% figure (3)
% hold on
% 
% plot(Size(:,2)+Size(:,3)+Size(:,4)+Size(:,5)+Size(:,6),'LineWidth', 3)
% plot(Size(:,2),'LineWidth', 3)
% %plot(Size(:,3),'LineWidth', 2)
% % plot(Size(:,4),'LineWidth', 2)
% % plot(Size(:,5),'LineWidth', 2)
% % plot(Size(:,6),'LineWidth', 2)
% hold on
% legend('LCS Filter','$\bf{E}_n$','Interpreter','latex')
% xlabel('Frame','FontSize',32,'LineWidth' , 3,'Interpreter','latex')
% ylabel('Memory Dimensionality','FontSize',32,'LineWidth' , 3,'Interpreter','latex')
% set(gca,'FontSize',32); 