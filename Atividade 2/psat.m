%Leitura dos parâmetros da EdE e dados experimentais
load Dados/tcpcwch4.txt
load Dados/parantoinech4.txt
load Dados/psatch4.txt

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
psata=zeros(nt,1);
psatpr=zeros(nt,1);
psatsrk=zeros(nt,1);

%Cálculo da pressão de saturação via Antoine
for k=1:nt
psata(k)=calc_psat_antoine(pa,tcalc(k));
end

%chute para a pressão de saturação em Pa
chutepr=psata(1);

for k=1:nt
disp(tcalc(k))
[psatpr(k)]=calc_psat_pr(tcalc(k),tc,pc,w,chutepr);
[psatsrk(k)]=calc_psat_srk(tcalc(k),tc,pc,w,chutepr);

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


