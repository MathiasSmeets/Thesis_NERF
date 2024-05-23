function channels_of_interest = get_channels_of_interest(neuron, offset, learner)
% get channels of a mouse according to its identified cluster
% offset resembles the amount of channels around the channel that detects a neuron
% learner resembles if it is a learner or control mouse

recording_database_updated_20200805_JCfreshStart_1;

if learner
    load("X:\Mathias\switch_data\depth\channel_m.mat"); % channel_m
    load("X:\Mathias\switch_data\correlations\template_cluster_m.mat"); % template_cluster
    
else
    load("X:\Mathias\switch_data\depth\channel_y.mat"); % channel_m
    load("X:\Mathias\switch_data\correlations\template_cluster_y.mat") % template_cluster

end

neurons = template_cluster{neuron};
channels_of_interest = channel_m{neuron}(neurons);

channels_to_add = [];
for i = 1:offset
    for j = 1:numel(channels_of_interest)
        channels_to_add = [channels_to_add, max(channels_of_interest(j)-i,1)];
        channels_to_add = [channels_to_add, min(channels_of_interest(j)+i,384)];
    end
end

channels_of_interest = unique([channels_of_interest, channels_to_add]);






