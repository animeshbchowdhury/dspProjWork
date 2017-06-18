function [VelocityMarker] = velocityPendCenter(pend_centers)

time = 1/30;

[r,c] = size(pend_centers);

for i = 2:r
    VelocityMarker(i,c-1) = (pend_centers(i,c-1) - pend_centers(i-1,c-1))/time;
    VelocityMarker(i,c) = (pend_centers(i,c) - pend_centers(i-1,c))/time;
end

