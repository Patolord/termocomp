function [a,b]=calc_ab_vw(tc,pc)
r=8.31451;
a=(27/64)*((r*tc)^2)/pc;
b=r*tc/(pc*8);
end
