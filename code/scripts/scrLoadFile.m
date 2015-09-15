clear all
%file=dir('C:\projectsENS\PolymerChainDynamicsResults\21_8_2015\*.mat');
file=dir('D:\Zha\Project\PolymerChainDynamicsResults\2_9_2015\*.mat');
% m=zeros(64,5000,numel(file));
% for i=1:length(file)
% load(file(i).name);
% m(:,:,i) = results.msd;
% end
% figure,
% for sIdx=1:5000
% MSD(sIdx,:)=mean(m(:,sIdx,:),3);
% end
% plot(MSD);
% 

for i=1:length(file)
load(file(i).name);
sTime(i)       =result.stickyTime(end);
%changeTime(i)  = result.stickyTime(end)/numel(result.stickyTime);
for j=1:numel(result.clustercomponents)
A(j) = sum(result.clustercomponents{j})/result.numcluster(j); 
end
meanCluster(i) = sum(A)/numel(A);
stickyTime{i}(1) = result.stickyTime(1);
stickyTime{i}(2:numel(result.stickyTime))=diff(result.stickyTime);
end
 figure,
plot(sTime,'LineWidth',4);hold on
plot(repmat(sum(sTime)/numel(file),[1 numel(file)]),'LineWidth',4','color','r');
legend('Time for each simulation','Mean time')
title('Time to assemble a only one cluster of 8 chains with ends beads move on the boundary','FontSize',24)
xlabel('number of simulations','FontSize',30)
ylabel('Time','FontSize',30,'Rotation',0)
set(gca,'FontSize',30)

figure,
[h,bins] = hist(sTime,50);
bar(bins,h)
legend('mean period time','FontSize',24)
xlabel('time','FontSize',30)
ylabel('Frequency','FontSize',30,'Rotation',90)
set(gca,'FontSize',30)
title('Distribution of Time to assemble a only one cluster of 8 chains with ends beads move on the boundary','FontSize',24)
legend('Time for each simulation')
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

figure,
[h,bins] = hist(meanCluster,50);
bar(bins,h)
legend('mean period time','FontSize',24)
xlabel('time','FontSize',30)
ylabel('Frequency','FontSize',30,'Rotation',90)
set(gca,'FontSize',30)
title('Distribution of Mean cluster of 8 chains with ends beads move on the boundary','FontSize',24)
legend('Mean cluster')


for stIdx=1:numel(file)
    meanPeriodTime(stIdx) = mean(stickyTime{stIdx});
end

figure,
[h,bins] = hist(meanPeriodTime,50);
bar(bins,h)
legend('mean period time','FontSize',24)
xlabel('time','FontSize',30)
ylabel('quantity','FontSize',30,'Rotation',90)
set(gca,'FontSize',30)

