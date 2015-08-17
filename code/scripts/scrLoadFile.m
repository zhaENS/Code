clear all
file=dir('D:\Zha\Project\PolymerChainDynamicsResults\17_8_2015\*.mat');
msd=[];
for i=1:length(file)
load(file(i).name);
msd(:,:,i) = results.msd;
end
figure,
for sIdx=1:1000
    for bIdx=1:64
MSD(sIdx,bIdx)=mean(sum(msd(bIdx,sIdx,:),3));
    end
end

plot(MSD);