%Leitura dos parâmetros da EdE e dados experimentais
load Dados/tcpcwch4.txt
load Dados/parantoinech4.txt
load Dados/psatch4.txt

% pexp=psatch4(:,2).*1000;
% texp=psatch4(:,1);
% nexp=length(texp);
% tc=tcpcwch4(1,1);
% pc=tcpcwch4(1,2);
% w=tcpcwch4(1,3);
% pa=parantoinech4;
% 
% tmin=psatch4(1,1);
% tmax=psatch4(nexp,1);
%  
% load tcpcwacetona.txt
% load parantoineacetona.txt
% load psatacetona.txt
% 
% pexp=psatacetona(:,2).*1000;
% texp=psatacetona(:,1);
% nexp=length(texp);
% 
% tc=tcpcwacetona(1,1);
% pc=tcpcwacetona(1,2);
% w=tcpcwacetona(1,3);
% pa=parantoineacetona;
% 
% tmin=psatacetona(1,1);
% tmax=psatacetona(nexp,1);

load Dados/tcpcwh2o.txt
load Dados/parantoineh2o.txt
load Dados/psath2o.txt

pexp=psath2o(:,2).*1000;
texp=psath2o(:,1);
nexp=length(texp);

tc=tcpcwh2o(1,1);
pc=tcpcwh2o(1,2);
w=tcpcwh2o(1,3);
pa=parantoineh2o;

tmin=psath2o(1,1);
tmax=psath2o(nexp,1);

%definição do intervalo de temperatura e determinação dos valores
% de a(T) e b
if tmax<tc
nt=100;
[tcalc]=linspace(tmin,tmax,nt);
[apr,bpr]=calcabpr(tcalc,tc,pc,w);
psata=zeros(nt,1);
psatpr=zeros(nt,1);
psatsrk=zeros(nt,1);

%Cálculo da pressão de saturação via Antoine
for k=1:nt
psata(k)=calcpsatantoine(pa,tcalc(k));
end

%chute para a pressão de saturação em Pa
chutepr=psata(1);

for k=1:nt
disp(tcalc(k))
[psatpr(k)]=calculaelvfipr(tcalc(k),apr(k),bpr,chutepr);
[psatsrk(k)]=calcpsatsrk(tcalc(k),tc,pc,w,chutepr);

chutepr=psatpr(k);
%chutepr=1e5;
%chutepr=psata(k);
end
else
disp('ajustar t')
end

figure(1)
plot(tcalc,psata,'b',tcalc,psatpr,'g',texp,pexp,'r*',tcalc,psatsrk,'k')
xlabel('t(K)')
ylabel('psat(Pa)')
legend('Antoine','PR','EXP','SRK')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psata] = calcpsatantoine(pa,t)
a=pa(1,1);
b=pa(1,2);
c=pa(1,3);
psat=10^(a-b/(t+c));
alfa=(1E5);
%psat em Pascal
psata=psat.*alfa;
end


