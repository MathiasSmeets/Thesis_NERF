function colornew = piuchiaropiuscuro(colorold,B)
% B in [0,1];
% B = 0.7
%color1 = lines(length(colortemp1))

if B < 0
    B = -B;
end
if B > 1
    B = 1/B;
end

colornew = ((colorold-B)*(1-B))+B;

% B = 0.7;
% x = [1:10]
% y = ones(1,length(x))
% XY =  x'*y;
% 
% figure
% hold on
% %XY =  x'*y;
% for kgt = 1:length(colortemp1)
% plot(XY(kgt,:),XY(:,kgt),'Color',color1(kgt,:))
% plot(XY(kgt,:)+.1,XY(:,kgt),'Color',color2(kgt,:))
% end
% 
% 
