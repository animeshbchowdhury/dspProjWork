function [w] = calculateFrequency(d11,z)

%     diff(1,:) = abs(d(z-6)-d(z-5));
%     diff(2,:) = abs(d(z-5)-d(z-4));
%     diff(3,:) = abs(d(z-4)-d(z-3));
%     diff(4,:) = abs(d(z-3)-d(z-2));
    diff(1,:) = abs(round(d11(z-3)-d11(z-2)))
    diff(2,:) = abs(round(d11(z-2)-d11(z-1)))
    diff(3,:) = abs(round(d11(z-1)-d11(z)))
    
%     w = int16(diff(2,:));
w = mode(diff);
