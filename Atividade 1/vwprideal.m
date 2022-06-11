clear
clc
%load dos dados de tc(K), pc(Pa) e w test
%H2O
%load tcpcwh2o.txt
%tc=tcpcwh2o(1);
%pc=tcpcwh2o(2);
%w=tcpcwh2o(3);

%C3H8
 load tcpcwc3h8.txt
 tc=tcpcwc3h8(1);
 pc=tcpcwc3h8(2);
 w=tcpcwc3h8(3);

%Constante universal dos gases (J/mol.K)
r=8.3145;

%temperatura da isoterma e intervalos e volumes molares
%H20
np=500;
%t=600;
%t=700;
%%t=647.3;
%vvw=linspace(4e-5,8e-4,np);
%vpr=linspace(2.5e-5,8e-4,np);

%C3H8
t=tc*0.7;
%t=tc;
%t=tc*1.2;
%vvw=linspace(0.98e-4,8e-4,np);
%vpr=linspace(0.66e-4,8e-4,np);

vvw=linspace(0.95e-4,8e-4,np);
vpr=linspace(0.6e-4,8e-4,np);

%cálculo de a e b para VW e PR
[avw,bvw]=calc_ab_vw(tc,pc);
[apr,bpr]=calc_ab_pr(t,tc,pc,w);

%definição do número de pontos do diagrama e tamanho dos vetores para
%armazenar as pressões calculadas

pcalcvw=zeros(np,1);
pcalcpr=zeros(np,1);
pcalcid=zeros(np,1);

%cálculo da pressão como função do volume molar
for i=1:np
pcalcvw(i)=r*t/(vvw(i)-bvw)-avw/(vvw(i)^2);
pcalcid(i)=r*t/vpr(i);
pcalcpr(i)=r*t/(vpr(i)-bpr)-apr/(vpr(i)^2+2*bpr*vpr(i)-bpr^2);
end

%plot das isotermas
figure(1)
plot(vvw,pcalcvw,'b',vpr,pcalcid,'r',vpr,pcalcpr,'g')
ylabel('p(Pa)')
xlabel('v(m3/mol)')
legend('VW','Gás ideal','PR')


