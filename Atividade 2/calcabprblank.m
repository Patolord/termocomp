function [a,b]=calcabpr(t,tc,pc,w)
r=8.31451;
ac=0.45724*((r*tc)^2)/pc;
b=0.0778*r*tc/pc;
nc=length(t);
a=zeros(nc,1);
kw=0.37464+1.5422*w-0.26992*(w^2);
for k=1:nc
tr=t(k)/tc;
alfa=(1+kw*(1-tr^0.5))^2;
a(k)=ac*alfa;
end
end

