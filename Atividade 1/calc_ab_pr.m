function [a,b]=calc_ab_pr(t,tc,pc,w)
r=8.31451;
ac=0.45724*((r*tc)^2)/pc;
b=0.0778*r*tc/pc;
kw=0.37464+1.5422*w-0.26992*(w^2);
tr=t/tc;
alfa=(1+kw*(1-tr^0.5))^2;
a=ac*alfa;
end