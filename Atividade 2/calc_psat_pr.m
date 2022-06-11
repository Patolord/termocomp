function [psat]=calc_psat_pr(t,tc,pc,w,chutep)
[a,b]=calcabpr(t,tc,pc,w);
[psat]=calculaelvfipr(t,a,b,chutep);
end

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

function [zl,zv,A,B]=resolveeospr(a,b,t,chute)
r=8.31451;
A=a*chute/(r*t)^2;
B=b*chute/(r*t);
alfa=-1+B;
beta=A-3*(B^2)-2*B;
gama=-A*B+B^2+B^3;
syms z;
solz=vpasolve(z^3+alfa*(z^2)+beta*z+gama==0);
[zl,zv]=escolhez(solz);
end

function [zl,zv]=escolhez(solz)
nz=length(solz);
nzf=0;
for k=1:nz
    l=isreal(solz(k));
    if l==1
        nzf=nzf+1;
    else
    end
end
solzf=zeros(nzf,1);
k=0;
for j=1:nz
    l=isreal(solz(j));
    if l==1
        k=k+1;
        solzf(k)=solz(j);
    else
    end
end
if nzf==1
zl=solzf(1);
zv=zl;
else
zl=min(solzf);
zv=max(solzf);
end
end

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

