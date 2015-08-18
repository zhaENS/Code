clear all
file=dir('C:\projectsENS\PolymerChainDynamicsResults\MSD64beadsMoveOntheBoundary\*.mat');
m=zeros(64,5000,numel(file));
for i=1:length(file)
load(file(i).name);
m(:,:,i) = results.msd;
end
figure,
for sIdx=1:5000
MSD(sIdx,:)=mean(m(:,sIdx,:),3);
end

plot(MSD);