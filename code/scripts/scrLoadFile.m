clear all
file=dir('C:\projectsENS\PolymerChainDynamicsResults\20_8_2015\*.mat');
%file=dir('D:\Zha\Project\PolymerChainDynamicsResults\20_8_2015\*.mat');
m=zeros(64,5000,numel(file));
% for i=1:length(file)
% load(file(i).name);
% m(:,:,i) = results.msd;
% end
% figure,
% for sIdx=1:5000
% MSD(sIdx,:)=mean(m(:,sIdx,:),3);
% end
% plot(MSD);

for i=1:length(file)
load(file(i).name);
sTime(i)       =result.stickyTime(end);
changeTime(i)  = result.stickyTime(end)/numel(result.stickyTime);
meanCluster(i) = sum(result.numcluster)/numel(result.numcluster);
end
% figure,
% plot(sTime,'LineWidth',4);hold on
% plot(repmat(sum(sTime)/numel(file),[1 numel(file)]);,'LineWidth',4','color','r');
% legend('Time for each simulation','Mean time')
% title('Time to assemble a only one cluster of 8 chains with ends beads move on the boundary','FontSize',24)
% xlabel('number of simulations','FontSize',30)
% ylabel('Time','FontSize',30,'Rotation',0)
% set(gca,'FontSize',30)
% figure,
% plot(changeTime,'LineWidth',4);hold on
% plot(repmat(sum(changeTime)/numel(changeTime),[1 numel(changeTime)]),'LineWidth',4,'color','r');
% legend('mean time for each simulation','mean time')
% title('mean time to form/break a cluster of 8 chains with ends beads move on the boundary','FontSize',24)
% xlabel('number of simulations','FontSize',30)
% ylabel('mean Time to form/break a cluster','FontSize',30,'Rotation',90)
% set(gca,'FontSize',30)

figure,
plot(meanCluster,'LineWidth',4);hold on
plot(repmat(sum(meanCluster)/numel(meanCluster),[1 numel(meanCluster)]),'LineWidth',4,'color','r');
legend('mean time for each simulation','mean time','FontSize',24)
title('mean number of clusters of 8 chains with ends beads move on the boundary','FontSize',24)
xlabel('number of simulations','FontSize',30)
ylabel('Mean cluster','FontSize',30,'Rotation',90)
set(gca,'FontSize',30)



