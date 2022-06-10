function [psat]=calcpsatsrk(t,tc,pc,w,chutep)
[a,b]=calcabsrk(t,tc,pc,w);
[psat]=calculaelvfisrk(t,a,b,chutep);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b]=calcabsrk(t,tc,pc,w)
r=8.31451;
ac=0.42747*((r*tc)^2)/pc;
b=0.08664*r*tc/pc;
kw=0.480+1.574*w-0.176*(w^2);
tr=t/tc;
alfat=(1+kw*(1-tr^0.5))^2;
a=ac*alfat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plv]=calculaelvfisrk(t,a,b,chute0)
[zl,zv,A,B]=resolveeossrk(a,b,t,chute0);
[fiv,fil]=calcfisrk(A,B,zl,zv);
rfi=fiv/fil;
dfi=1-rfi;
tol=1E-5;
it=0;
while abs(dfi)>=tol
chute1=chute0*fil/fiv;
it=it+1;
[zl,zv,A,B]=resolveeossrk(a,b,t,chute1);
[fiv,fil]=calcfisrk(A,B,zl,zv);
rfi=fiv/fil;
dfi=1-rfi;
chute0=chute1;
end
plv=chute0;
fv=plv*fiv;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fiv,fil]=calcfisrk(A,B,zl,zv)
fil=exp(zl-1-log(zl-B)-(A/B)*log((zl+B)/zl));
fiv=exp(zv-1-log(zv-B)-(A/B)*log((zv+B)/zv));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zl,zv,A,B]=resolveeossrk(a,b,t,chute)
r=8.314;
A=a*chute/(r*t)^2;
B=b*chute/(r*t);
alfa=-1;
beta=A-B-B^2;
gama=-A*B;
[solz]=resolvepolz(alfa,beta,gama);
[zl,zv]=escolhez(solz);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [solz]=resolvepolz(alfa,beta,gama)
syms x;
eq=x^3+alfa*(x^2)+beta*x+gama==0;
solz=vpasolve(eq,x);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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