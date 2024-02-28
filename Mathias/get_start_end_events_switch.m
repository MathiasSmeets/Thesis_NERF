function [EventSizeM, EventSizeY, max_m, max_y, waitingeventM, waitingeventY, switcheventM, switcheventY] = get_start_end_events_switch()
max_m = 0;
max_y = 0;
for i = 1:12
    if i==1 %25/07/22 // Ptf1a
        lim_start=1;
        lim_event=324;
        waiting_start = 3665;
        waiting_stop = 8700;
        switch_start = 563; % switch end should be 10 min after this

    elseif i==2 %20220819 // Ptf1a
        lim_start=1;
        lim_event=336;
        waiting_start = 3659;
        waiting_stop = 7270;
        switch_start =  578;

    elseif i==3 %19/08/22 // Eng (waiting: 3676-7264)
        lim_start=1;
        lim_event=749;
        waiting_start = 3655;
        waiting_stop = 7270;
        switch_start =  578;

    elseif i==4 %21/09/22 // Ptf1a (waiting: 3709-7382)
        lim_start=2;
        lim_event=1070;
        waiting_start = 3709;
        waiting_stop = 7382;
        switch_start =  1329;
    elseif i==5 %23/09/22 // Ptf1a
        lim_start=1;
        lim_event=128;
        waiting_start = 3921;
        waiting_stop = 7924;
        switch_start =  1479;

    elseif i==6 %11/10/22 // Ptf1a
        lim_start=1;
        lim_event=103;
        waiting_start = 3940;
        waiting_stop = 7563;
        switch_start = 863;
    elseif i==7 %16/11/22 // Ptf1a
        lim_start=3;
        lim_event=419;
        waiting_start = 3040;
        waiting_stop = 6902;
        switch_start = 545;
    elseif i==8 %24/11/22 // ENG
        lim_start=1;
        lim_event=164;
        waiting_start = 3075;
        waiting_stop = 7212;
        switch_start = 291;
    elseif i==9 %28/11/22 // ENG
        lim_start=1;
        lim_event=225;
        waiting_start = 3207;
        waiting_stop = 7213;
        switch_start = 383;
    elseif i==10 %25/05/22 // horridge only
        lim_start=2;
        lim_event=3372;
        waiting_start = 1;
        waiting_stop = 10;
        switch_start = 1;
    elseif i==11 %22/06/22 // horridge only
        lim_start=1;
        lim_event=610;
        waiting_start = 1;
        waiting_stop = 10;
        switch_start = 1;
    elseif i==12 %15/11/22 // horridge only
        lim_start=1;
        lim_event=183;
        waiting_start = 1;
        waiting_stop = 10;
        switch_start = 1;

        % Initialization parameter in case of problem
    else
        lim_start=2;
        lim_event=length(events(evtid).onsets);
    end

    EventSizeM(i,1)=lim_start;
    EventSizeM(i,2)=lim_event;
    waitingeventM(i,1) = waiting_start;
    waitingeventM(i,2) = waiting_stop;
    switcheventM(i,1) = switch_start;

    if max_m < lim_event-lim_start
        max_m = lim_event-lim_start;
    end
end

for i = 1:12
    if i==1 %05/07/2022 //Ptf1a //waiting bin: 3842, 7382];
        lim_start=1;
        lim_event=712;
        waiting_start = 3842;
        waiting_stop = 7382;
        switch_start=956;
    elseif i==2 %26/07/22 // Ptf1a// waiting bin: 4438, 7443];
        lim_start=1;
        lim_event=568;
        waiting_start = 4438;
        waiting_stop = 7443;
        switch_start=817;
    elseif i==3 %05/09/22 // Ptf1a// waiting bin: 3842, 7265];
        lim_start=1;
        lim_event=279;
        waiting_start = 3847;
        waiting_stop = 7262;
        switch_start=756;
    elseif i==4 %22/09/22 // Ptf1a// waiting bin: 4148, 7806];
        lim_start=1;
        lim_event=1035;
        waiting_start = 4148;
        waiting_stop = 7803;
        switch_start=2319;
    elseif i ==5 %27/09/2022 // Ptf1a // waiting bin: 6963, 10563
        lim_start=1;
        lim_event=127;
        waiting_start = 6963;
        waiting_stop = 10563;
        switch_start=1460;
    elseif i ==6 %06/10/2022 // Ptf1a
        lim_start=1;
        lim_event=271;
        waiting_start = 3915;
        waiting_stop = 7623;
        switch_start=1007;
    elseif i ==7 %18/10/2022 // Ptf1a
        lim_start=1;
        lim_event = 102;
        waiting_start = 3010;
        waiting_stop = 6600;
        switch_start=227;
    elseif i ==8 %06/09/2022 // ENG
        lim_start=107;
        lim_event = 862;
        waiting_start = 3687;
        waiting_stop = 7561;
        switch_start=1116;
    elseif i ==9 %13/04/2022 // ENG
        lim_start=1;
        lim_event = 1787;
        waiting_start = 1;
        waiting_stop = 10;
        switch_start=1;
    elseif i ==10 %27/04/2022 // ENG
        lim_start=1;
        lim_event = 787;
        waiting_start = 1;
        waiting_stop = 10;
        switch_start=1;
    elseif i ==11 %29/04/2022 // ENG
        lim_start=1;
        lim_event = 239;
        waiting_start = 1;
        waiting_stop = 10;
        switch_start=1;
    elseif i ==12 %07/06/2022 // ENG
        lim_start=1;
        lim_event = 2664;
        waiting_start = 1;
        waiting_stop = 10;
        switch_start=1;

                % Initialization parameter in case of problem
    else
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
end
