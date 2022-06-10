function [plv]=calculaelvfipr(t,a,b,chute0)
%
[zl,zv,A,B]=resolveeospr(a,b,t,chute0);
[fiv,fil]=calcfipr(A,B,zl,zv);
rfi=fiv/fil;
dfi=1-rfi;
%

tol=1e-5;

while abs(dfi)>=tol
chute1=chute0/(rfi);

%
[zl,zv,A,B]=resolveeospr(a,b,t,chute1);
[fiv,fil]=calcfipr(A,B,zl,zv);
rfi=fiv/fil;
dfi=1-rfi;
%
chute0=chute1;
end
plv=chute0;
end

function [fiv,fil]=calcfipr(A,B,zl,zv)
fil=exp(zl-1-log(zl-B)+(A/(B*(2^1.5)))*log((zl+B*(1-2^0.5))/(zl+B*(1+2^0.5))));
fiv=exp(zv-1-log(zv-B)+(A/(B*(2^1.5)))*log((zv+B*(1-2^0.5))/(zv+B*(1+2^0.5))));
end