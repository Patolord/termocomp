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
