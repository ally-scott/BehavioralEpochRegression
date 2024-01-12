%% load stuff
clear all
close all
clc

cd('C:\Users\Ally\Desktop\FP-analysis-variableReward\FP_analysis\FP-analysis\Parker encoding model\encoding_results\satge5criteria\dff\');
load('allSubjResultskernel_Shifted_all.mat');

%% Create graph of all subjects ( SEM area)
%reorganise z-score data for all animals

DSkernel_Shifted_all=[];
PEkernel_Shifted_all=[];
Lickkernel_Shifted_all=[]; 
modelsum_Shifted_all=[];
DSzblueAllTrials_Shifted_all=[];
DS_rat={};
DS_timeLock=[];
DS_trial=[];
DS_event={};

PE_rat={};
PE_timeLock=[];
PE_trial=[];
PE_event={};

lox_rat={};
lox_timeLock=[];
lox_trial=[];
lox_event={};


fns = fieldnames(kernel_Shifted_all.kernels_DSTrials);

for subject= 1:length(fns);
% create matrix with all DS onset DSkernels for all animals (1st sheet in
% 3D matrix)


DSkernel_Shifted_all= cat(2,DSkernel_Shifted_all,kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1));
DS_rat=cat(2,DS_rat,repmat({fns{subject}},size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1))));
DS_timeLock=cat(2,DS_timeLock,repmat([timeLock'],[1 size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1),2)]));
DS_trial=cat(2,DS_trial,repmat((1:size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1),2)),[size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1),1) 1]));
DS_event=cat(2,DS_event,repmat({'DS'}, size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1))));
% create matrix with all PE DSkernels for all animals(2nd sheet in
% 3D matrix)
PEkernel_Shifted_all= cat(2,PEkernel_Shifted_all,kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,2));
PE_rat=cat(2,PE_rat,repmat({fns{subject}},size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,2))));
PE_timeLock=cat(2,PE_timeLock,repmat([timeLock'],[1 size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,2),2)]));
PE_trial=cat(2,PE_trial,repmat((1:size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,2),2)),[size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1),1) 1]));
PE_event=cat(2,PE_event,repmat({'PE'}, size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,2))));
% create matrix with all lick DSkernels for all animals(3rd sheet in
% 3D matrix)

Lickkernel_Shifted_all= cat(2,Lickkernel_Shifted_all,kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,3));
lox_rat=cat(2,lox_rat,repmat({fns{subject}},size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,3))));
lox_timeLock=cat(2,lox_timeLock,repmat([timeLock'],[1 size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,3),2)]));
lox_trial=cat(2,lox_trial,repmat((1:size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,3),2)),[size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,1),1) 1]));
lox_event=cat(2,lox_event,repmat({'lox'}, size(kernel_Shifted_all.kernels_DSTrials.(fns{subject})(:,:,3))));



% modelsum matrix
modelsum_Shifted_all= cat(2,modelsum_Shifted_all,kernel_Shifted_all.gcamp_model_sum.(fns{subject})(:,:));

%DSzalltrials matrix
DSzblueAllTrials_Shifted_all= cat(2,DSzblueAllTrials_Shifted_all,kernel_Shifted_all.DSzblueAllTrials.(fns{subject})(:,:)');
end
%% gramm plots

%DS kernel
figure()
g=gramm('x',timeLock,'y',DSkernel_Shifted_all');
g.stat_summary('geom','area','type','sem')
g.set_names('x','Time from DS onset(sec)','y','Regression Coefficient','group','DSonset');
g.axe_property('YLim',[-1 1.5]);
g.set_title('DS Onset Kernel');
g.draw()

         gcf;
        [filepath,name,ext] = fileparts(file_name);
        figsave_name=strcat('Avg_DSKernel');
        cd(strcat(figsave_folder,'\Avg kernels across animals\'));;
        savefig(figsave_name);

%PE kernel
figure()
g=gramm('x',timeLock,'y',PEkernel_Shifted_all')
g.stat_summary('geom','area','type','sem')
g.set_names('x','Time from PE(sec)','y','Regression Coefficient');
g.set_color_options('map','d3_10');
g.axe_property('YLim',[-1 1.5]);
g.set_title('Port Entry Kernel');
g.draw()

         gcf;
        [filepath,name,ext] = fileparts(file_name);
        figsave_name=strcat('Avg_PEKernel');
        cd(strcat(figsave_folder,'\Avg kernels across animals\'));;
        savefig(figsave_name);


%Lick kernel
figure()
g=gramm('x',timeLock,'y',Lickkernel_Shifted_all')
g.stat_summary('geom','area','type','sem')
g.set_names('x','Time from Initial Lick(sec)','y','Regression Coefficient');
g.set_color_options('map','brewer2');
g.axe_property('YLim',[-1 1.5]);
g.set_title('Initial Lick Kernel');
g.draw()

         gcf;
        [filepath,name,ext] = fileparts(file_name);
        figsave_name=strcat('Avg_LickKernel');
        cd(strcat(figsave_folder,'\Avg kernels across animals\'));;
        savefig(figsave_name);




%% create table
% b one column, rat one colum, timeLock one column, trial and eventype

% squeeze
auckernel=table();
b=cat(2,DSkernel_Shifted_all,PEkernel_Shifted_all,Lickkernel_Shifted_all)';
b = b(:)';
auckernel.b =b';
rat=cat(2,DS_rat,PE_rat,lox_rat)';
rat = rat(:)';
auckernel.rat = rat';
time=cat(2,DS_timeLock,PE_timeLock,lox_timeLock)';
time = time(:)';
auckernel.time = time';
trial=cat(2,DS_trial,PE_trial,lox_trial)';
trial = trial(:)';
auckernel.trial = trial';
event=cat(2,DS_event,PE_event,lox_event)';
event = event(:)';
auckernel.event = event';

auckernel = sortrows(auckernel,'trial','ascend');
auckernel = sortrows(auckernel,'rat','ascend');
auckernel = sortrows(auckernel,'event','ascend');

%% calculate means for animals 

kernelaucavgt=table();
kernelaucavgt= groupsummary(auckernel,{'event','rat','time'},'mean',{'b'});
kernelaucavgt=kernelaucavgt(kernelaucavgt.time<=5 & kernelaucavgt.time>=0,:);

%% table with time bins in 0.1 sec
timebinaucavgt=kernelaucavgt(kernelaucavgt.time<=5 & kernelaucavgt.time>=0,:);
x=1;
ma=table();
for t=1:length(timebinaucavgt.mean_b)/201
    if x==1
        b=0;
    else
        b=1;
       
    end
 ratevent= x:t*(length(timebinaucavgt.mean_b)/(length(subjects)*length(cons)));
%regression coefficients
 ma.b_bins(((t-1)*50)+1:(t*50))= nanmean(reshape(timebinaucavgt.mean_b(x+b:ratevent(end)-1),[],50,1));
% % metadata

% time
 ma.bin(((t-1)*50)+1:(t*50))= nanmean(reshape(timebinaucavgt.time(x+b:ratevent(end)-1),[],50,1));
 
 % rat
 binrat=(reshape(timebinaucavgt.rat(x+b:ratevent(end)-1),[],50,1));
 ma.rat(((t-1)*50)+1:(t*50))= binrat(1,:);
% eventtype
 binevent=(reshape(timebinaucavgt.event(x+b:ratevent(end)-1),[],50,1));
 ma.event(((t-1)*50)+1:(t*50))= binevent(1,:);
% sex and led addded in R




x=t*((length(timebinaucavgt.mean_b)/(length(subjects)*length(cons))));

end
%% calculate AUC
for ratdaymean=1:(length(kernelaucavgt.time)/201)
    if ratdaymean==1
kernelaucavgt.aucvalues(1:ratdaymean*201)= trapz(kernelaucavgt.time((1:ratdaymean*201)), kernelaucavgt.mean_b(~isnan(1:ratdaymean*201)));
   else
kernelaucavgt.aucvalues((((ratdaymean-1)*201)+1):(((ratdaymean-1)*201)+201))= trapz(kernelaucavgt.time((((ratdaymean-1)*201)+1):(((ratdaymean-1)*201)+201)), kernelaucavgt.mean_b((((ratdaymean-1)*201)+1):(((ratdaymean-1)*201)+201)));       
    end
end

kernelaucavgtab= groupsummary(kernelaucavgt,{'event','rat'},'mean',{'aucvalues'});
%% save table for R stats\
cd('G:\Shared drives\Richard Lab\Papers and Abstracts\Papers in Progress\Alexandra Scott\GADVPFP\Figs\fig4\');
writetable(kernelaucavgtab,'criteriadaykernel465_dff.csv')
writetable(ma,'criteriaday465kernelbins_dff.csv')
