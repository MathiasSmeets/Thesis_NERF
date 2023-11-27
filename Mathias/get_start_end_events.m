function [EventSizeM, EventSizeY] = get_start_end_events()


for i = 7:34
    lim_start =0;
    lim_event =0;
    if i==7 % WT
        %End of master horridge number 8 at 32min
        lim_start=1; %length(find(events.onsets<(22*60)));
        lim_event=1604; %length(find(events.onsets<(32*60) & events.onsets>(22*60)));    
    elseif i==8 % WT
        %End of master horridge number 8 at 32min
        lim_start=1; %length(find(events.onsets<(22*60)));
        lim_event=156; %length(find(events.onsets<(32*60) & events.onsets>(22*60)));
    elseif i==9 % Tlx3
        %Start horridge 22 End of master horridge number 9 at 32min
        lim_start=259;
        lim_event=1036;
    elseif i==10 %09/11/21 // Tlx3
        lim_start=1;
        lim_event=1394;
    elseif i==11 %17/11/21 // Ptf1a
        lim_start=1; 
        lim_event=2035;
    elseif i==12 % Ptf1a
        lim_start=21; %OK checked CB;
        lim_event=578;
    elseif i==13 %06/12/21 // Ptf1a
        lim_start=1;
        lim_event=981;
    elseif i==14 %04/01/22 // Ptf1a
        lim_start=3;
        lim_event=202;
    elseif i==15 %11/01/22 // Ptf1a
        lim_start=2;
        lim_event=1291;
    elseif i==16 %18/01/22 // Ptf1a
        lim_start=3;
        lim_event=1291;
    elseif i==17 %25/01/22 // Ptf1a
        lim_start=133; % OK checked CB;;
        lim_event=1345; %length(find(events.onsets<(60*25) & events.onsets>(18*60)));
    elseif i==18 %07/02/22 // Ptf1a
        lim_start=3;%find(events.onsets>(17*60),1)+1;
        lim_event=1228;%length(find(events.onsets<(60*25) & events.onsets>(17*60)));
    elseif i==19 %14/02/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        lim_event=286; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==20 %21/02/22 // Ptf1a
        lim_start=4; %find(events.onsets>(15*60),1)+1;
        lim_event=885; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==21 %10/03/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        lim_event=569; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==22 %14/03/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        lim_event=293; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==23 %24/03/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        lim_event=1787; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==24 %19/04/22 // Ptf1a
        lim_start=2; %find(events.onsets>(15*60),1)+1;
        lim_event=788; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==25 %28/04/22 // Ptf1a
        lim_start=1; %find(events.onsets>(15*60),1)+1;
        lim_event=239; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==26 %25/05/22 // Ptf1a
        lim_start=2; %find(events.onsets>(15*60),1)+1;
        lim_event=3372; %length(find(events.onsets<(60*25) & events.onsets>(15*60)));
    elseif i==27 %22/06/22 // Ptf1a
        lim_start=1; 
        lim_event=610;
    elseif i==28 %25/07/22 // Ptf1a
        lim_start=1; 
        lim_event=324;
    elseif i==29 %19/08/22 // Ptf1a
        lim_start=1; 
        lim_event=336;
    elseif i==30 %21/09/22 // Ptf1a
        lim_start=2; 
        lim_event=1070;
    elseif i==31 %23/09/22 // Ptf1a
        lim_start=1; 
        lim_event=128;
    elseif i==32 %11/10/22 // Ptf1a 
        lim_start=1; 
        lim_event=103;         
    elseif i==33 %15/11/22 // Ptf1a 
        lim_start=1; 
        lim_event=183; %corrected
    elseif i==34 %16/11/22 // Ptf1a   (Ran as a Yoked the Copy pasted)
        lim_start=3;
        lim_event=419;
        
    end
    EventSizeM(i,1)=lim_start;
    EventSizeM(i,2)=lim_event;
end



for i =7:34
    lim_start =0;
    lim_event =0;
    
    if i==7 %31/08/21       
        lim_start=1; %1; horridge fail first time %length(find(events.onsets<(17*60)));
        lim_event=1603; %length(find(events.onsets<(60*25) & events.onsets>(60*17)));
    elseif i==8
        %End of yoked horridge number 8 at 27min / start at 17
        lim_start=97; %1; horridge fail first time %length(find(events.onsets<(17*60)));
        lim_event=252; %length(find(events.onsets<(60*25) & events.onsets>(60*17)));
    elseif i==9 %16/11/21
        lim_start=1;
        lim_event=1382;
    elseif i==10 %22/11/21
        lim_start=1; %length(find(events.onsets<(15*60)));
        lim_event=1826; %length(find(events.onsets<(60*25) & events.onsets>(60*15)));
    elseif i==11 %02/12/21
        lim_start=1; %11;%length(find(events.onsets<(15.1*60)));
        lim_event=558;%length(find(events.onsets<(60*25) & events.onsets>(60*15.1)));
    elseif i==12 %07/12/21
        lim_start=1;
        lim_event=998;
    elseif i==13 %05/01/22
        lim_start=4;
        lim_event=204; %197; %length(find(events.onsets<(60*25) & events.onsets>(60*17)));
    elseif i==14 %14/01/22
        lim_start=2;
        lim_event=1072;
    elseif i==15 %21/01/22
        lim_start=2;
        lim_event=1271;
    elseif i==16 %01/02/22
        %test=events.onsets;
        lim_start=1;%find(test>(14*60) & test<(16*60),1)+1;
        lim_event=1213;%length(find(test<(60*23) & test>(60*15)));
    elseif i==17 %08/02/22
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=1224;%length(find(test<(60*25) & test>(60*17)));
    elseif i==18 %17/02/22
        lim_start=3; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=287;%length(find(test<(60*25) & test>(60*17)));
    elseif i==19 %01/03/22
        lim_start=76; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=957; %length(find(test<(60*25) & test>(60*17)));
    elseif i==20 %11/03/22
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=569; %length(find(test<(60*25) & test>(60*17)));
    elseif i==21 %23/03/22
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=293; %length(find(test<(60*25) & test>(60*17)));
    elseif i==22 %13/04/22
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=1786; %length(find(test<(60*25) & test>(60*17)));
    elseif i==23 %27/04/22
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=787; %length(find(test<(60*25) & test>(60*17)));
    elseif i==24 %29/04/22
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=239; %length(find(test<(60*25) & test>(60*17)));     
    elseif i==25 %07/06/22
        lim_start=1; %find(test>(15*60) & test<(18*60),1)+1;
        lim_event=2648; %length(find(test<(60*25) & test>(60*17)));     
    elseif i==26 %05/07/22 // Ptf1a
        lim_start=1;
        lim_event=610;
    elseif i==27 %26/07/22 // Ptf1a
        lim_start=1;
        lim_event=568;
    elseif i==28 %05/09/22 // Ptf1a
        lim_start=1;
        lim_event=279;
    elseif i==29 %22/09/22 // Ptf1a
        lim_start=1;
        lim_event=1035;
    elseif i==30 %27/09/22 // Ptf1a   
        lim_start=1;
        lim_event=127;
    elseif i==31 %06/10/22 // Ptf1a   
        lim_start=2;
        lim_event=221;
    elseif i==32 %18/10/22 // Ptf1a   
        lim_start=1;
        lim_event=183; %corrected
    elseif i==33 %07/12/22 // Ptf1a
        lim_start=107;
        lim_event=287;
    elseif i==34 %19/01/23 // Ptf1a
        lim_start=1;
        lim_event=416;
    end
    EventSizeY(i,1)=lim_start;
    EventSizeY(i,2)=lim_event;

end

end