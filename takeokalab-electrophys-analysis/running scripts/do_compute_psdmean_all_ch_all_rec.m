clear all
%% Parameters

% data and channels to be excluded in particular cases
%MASTER
% data{1} = 'D:\Experiments\Data\SC\20200907';
% chlist{1} = [];
% data{2} = 'D:\Experiments\Data\SC\20200921\Master_3';
% chlist{2} = [];
% data{3} = 'D:\Experiments\Data\SC\20201002\Master 7';
% chlist{3} = [];
% % 
% %BAD MASTER

%  data{1} = 'D:\Experiments\Data\SC\20200916\Master 2';
% chlist{1} = [7 8 9 10];
% bad = 0;

bad = 1;
data{1} = 'D:\Experiments\Data\SC\20200921\Master_4';
chlist{1} = [];
data{2} = 'D:\Experiments\Data\SC\20200923\Master_5';
chlist{2} = [];

%YOKED
% data{1} = 'D:\Experiments\Data\SC\20200914';
% chlist{1} = [8,11];
% data{2} = 'D:\Experiments\Data\SC\20200916\Yoked 3';
% chlist{2} = [7,8,9,10];
% data{3} = 'D:\Experiments\Data\SC\20200923\Yoked_4';
% chlist{3} = [9,11,32];
% data{4} = 'D:\Experiments\Data\SC\20200928\Yoked 5';
% chlist{4} = [];
% data{1} = 'D:\Experiments\Data\SC\20201005';
% chlist{1} = [];
bad = 0;
%NEW
% good master
% data{1} = 'D:\Experiments\Data\SC\20201012\Master9';
% chlist{1} = [];
% 
% %good master, but weird result
% data{2} = 'D:\Experiments\Data\SC\20201014\Master10';
% chlist{2} = [];

nchs = 32;

for i  = 1:length(data)
    disp(sprintf('Processing data #%d',i))
    %% UPLOAD MEMORY POINTERS TO DATA
    s = dir(data{i});
    
   %find folder with spontaneus rec before horridge
   list = regexp({s.name},'spont_before','match');
   n = find(~cellfun(@isempty,list));
   s2 = dir(fullfile(data{i},s(n).name));
   dirFlags = [s2.isdir];
   subFold = s(dirFlags);
   %switch difference between recordings in .dat and in .openephys formats:
   %the latter need to be pre-preocessed --> use the function
   %"preprocess_openephys_data()"
   if length(subFold)>2
       D1 = load_open_ephys_binary(fullfile(data{i},s(n).name,'experiment1','recording1','\structure.oebin'), 'continuous',1,'mmap');
   else
       datname = dir(fullfile(data{i},s(n).name,'*.dat'));
       contFile = fullfile(data{i},s(n).name,datname.name);
       file = dir(contFile);
       samples = file.bytes/2/nchs;
       D1.Data = memmapfile(contFile,'Format',{'int16' [32 samples] 'mapped'});
   end
   
   %repeat procedure for spontaneous rec after horridge
   list = regexp({s.name},'spont_after','match');
   n = find(~cellfun(@isempty,list));
   s2 = dir(fullfile(data{i},s(n).name));
   dirFlags = [s2.isdir];
   subFold = s(dirFlags);
   %switch difference between recordings in .dat and in .openephys formats:
   %the latter need to be pre-preocessed --> use the function
   %"preprocess_openephys_data()"
   if length(subFold)>2
       D2 = load_open_ephys_binary(fullfile(data{i},s(n).name,'experiment1','recording1','\structure.oebin'), 'continuous',1,'mmap');
   else
       datname = dir(fullfile(data{i},s(n).name,'*.dat'));
       contFile = fullfile(data{i},s(n).name,datname.name);
       file = dir(contFile);
       samples = file.bytes/2/nchs;
       D2.Data = memmapfile(contFile,'Format',{'int16' [32 samples] 'mapped'});
   end
   
   %repeat procedure for recovery rec for bad masters
   if bad
       list = regexp({s.name},'spont_recovery','match');
       n = find(~cellfun(@isempty,list));
       s3 = dir(fullfile(data{i},s(n).name));
       dirFlags = [s3.isdir];
       subFold = s(dirFlags);
       %switch difference between recordings in .dat and in .openephys formats:
       %the latter need to be pre-preocessed --> use the function
       %"preprocess_openephys_data()"
       if length(subFold)>2
           D3 = load_open_ephys_binary(fullfile(data{i},s(n).name,'experiment1','recording1','\structure.oebin'), 'continuous',1,'mmap');
       else
           datname = dir('*.dat');
           contFile = fullfile(data{i},s(n).name,datname.name);
           file = dir(contFile);
           samples = file.bytes/2/nchs;
           D3.Data = memmapfile(contFile,'Format',{'int16' [32 samples] 'mapped'});
       end
   end
   
   %% upload one channel at the time and compute psdmean 
   psdmean = [];
   f = [];

   for j = 1:nchs
       if ~ismember(j,chlist{i})
           datain{1} = double(D1.Data.Data(1).mapped(j,:));
           datain{2} = double(D2.Data.Data(1).mapped(j,:));
           if bad
               datain{3} = double(D3.Data.Data(1).mapped(j,:));
           end
           [f, psdmeantemp] = lfp_comparison_recordings(datain,'select_channels',1,'filter_signal',[0.5 100],'plot_results',0);
            psdmean{j} = psdmeantemp';
       else
           psdmean{j} = NaN(1,size(psdmean,2));
       end
   end
   if ~bad
       save(fullfile(data{i},'psdmean_before_after_allch.mat'), 'psdmean', 'f', '-v7.3')
   else
       save(fullfile(data{i},'psdmean_before_after_recovery_allch.mat'), 'psdmean', 'f', '-v7.3')
   end 
end
       
       
       
       
       
       
       
       
       