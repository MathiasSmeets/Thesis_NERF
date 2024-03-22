function [EventSizeM, EventSizeY, max_m, max_y, waitingeventM, waitingeventY, switcheventM, switcheventY] = get_start_end_events_np2()
max_m = 0;
max_y = 0;
for i = 1:3
    if i==1 %05/07/2022 //Ptf1a //waiting bin: 3842, 7382];
        lim_start=2;
        lim_event = length(events(evtid).onsets);
        waiting_start = events(evtid).onsets(lim_start)+600;
        waiting_stop = events(evtid).onsets(lim_start)+2400;
        switch_start=1;
    elseif i==2 %26/07/22 // Ptf1a// waiting bin: 4438, 7443];
        lim_start=2;
        lim_event = length(events(evtid).onsets);
        waiting_start = events(evtid).onsets(lim_start)+600;
        waiting_stop = events(evtid).onsets(lim_start)+2400;
        switch_start=1;    
        
        
    elseif i == 3
        lim_start=2;
        lim_event = length(events(evtid).onsets);
        waiting_start = events(evtid).onsets(lim_start)+600;
        waiting_stop = events(evtid).onsets(lim_start)+2400;
        switch_start=1;  
    else% Initialization parameter in case of problem
        lim_start=2;
        lim_event=length(events(evtid).onsets);
    end

    EventSizeY(i,1)=lim_start;
    EventSizeY(i,2)=lim_event;
    waitingeventY(i,1) = waiting_start;
    waitingeventY(i,2) = waiting_stop;
    switcheventY(i,1) = switch_start;

    if max_y < lim_event-lim_start
        max_y = lim_event-lim_start;
    end
end
EventSizeM = 0;
waitingeventM = 0;
switcheventM = 0;
end
