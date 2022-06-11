function [psata] = calc_psat_antoine(pa,t)
a=pa(1,1);
b=pa(1,2);
c=pa(1,3);
psat=10^(a-b/(t+c));
alfa=(1E5);
%psat em Pascal
psata=psat.*alfa;
end


