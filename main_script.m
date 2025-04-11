clearvars
close all
clc

%% Caricamento e plot dei dati originali

load glucose_data_fm.mat

figure
plot(ts, ys, 'LineWidth', 1.5, 'Color', 'Blue')
title('segnale glicemico CGM')
xlabel('tempo [min]')
ylabel('segnale CGM [mg/dl]')
xlim([ts(1) ts(end)])


%% Deconvoluzione regolarizzata

ns = length(ys);

% Burn in
ts_og = ts;
Tburn = ts(end);
ts = ts+Tburn; % prendo come burn in la durata dell'intera finestra dei dati

% Griglia virtuale (Tv=1)
Tv = 1;
tv = (Tv:Tv:ts(end))';
nv = length(tv);

% Mappatura dei virtually-missing samples
vt = zeros(ts(end), 1);
for i=1:ns
    vt(tv==ts(i)) = 1;
end

% Creazione di Gv e di G
NUM = [1];
DEN = [1 -1];
id = zeros(nv, 1);
id(1) = 1;
gv = filter(NUM, DEN, id);
r = zeros(1, nv);
r(1) = gv(1);
Gv = toeplitz(gv, r);
G = Gv(vt==1, :);

m = 2;

% Creazione di F
if m==0
    F = eye(nv);
else
    c = [1; -1; [zeros(nv-2, 1)]]';
    r = [1; [zeros(nv-1, 1)]];

    delta = toeplitz(c, r);

    F = eye(nv);
    for k=1:m
        F = F*delta;
    end
end

% Matrice di covarianza dell'errore di misura: sigma2*B
sigma2 = 9;
B = eye(ns);

% Creazione di vettori e matrici relativi alla varianza del rumore
IB2 = ones(ns, 1)./diag(B).^(1/2); % vettore
IB = diag(diag(B).^-1); % B^-1

% Creazione di vettori e matrici che rimangono invariate
b = G'*IB*ys;
GT_IB_G = G'*IB*G;
FT_F = (F')*F;
GT_IB = G'*IB;

% Criterio di consistenza 1
convergenza = 0;
gamma_min = 1e-10; 
gamma_max = 1e+10;
while convergenza==0
    gamma = 10^((log10(gamma_min)+log10(gamma_max))/2);
    
    xk = zeros(nv, 1);
    rk = b;
    pk = rk;
    for k=1:nv
        
        temp1 = filter(F(:, 1), 1, pk);
        Q2 = gamma*flipud(filter(F(:, 1), 1, flipud(temp1)));
        
        temp2 = filter(NUM, DEN, pk);
        temp3 = temp2(vt==1);
        temp4 = diag(IB).*temp3;
        temp4_riempito = zeros(nv, 1);
        temp4_riempito(vt==1) = temp4;
        Q1 = flipud(filter(NUM, DEN, flipud(temp4_riempito)));
        
        Qpk = Q1 + Q2;

        alphak = (pk'*rk)/(pk'*Qpk);
        xk = xk + alphak*pk;

        rk = rk - alphak*Qpk;
        betak = -(rk'*Qpk/(pk'*Qpk));
        pk = rk + betak*pk;
    end
    
    u_hat = xk;
    
    WESS = sum((filter(F(:, 1), 1, u_hat)).^2);
    D = (GT_IB_G+gamma*FT_F)\GT_IB;
    q = trace(G*D);
    lambda2 = sigma2/gamma;
    if WESS> lambda2*q
        gamma_max = gamma; % devo abbassare gamm 
    else
        gamma_min = gamma; % devo alzare gamm 
    end
    if abs((WESS-lambda2*q)/WESS) < 0.01
        convergenza=1;
    end
end

% Predizione
yp = Gv*u_hat;

% Residui pesati
resp = IB2.*(ys-yp(vt==1));

% Residui normalizzati
resn = resp/sqrt(sigma2);

% Rimozione burn in
ind_in = (Tburn/Tv)+1;
ind_fin = length(tv);
tv = tv(ind_in:ind_fin)-Tburn;
u_hat = u_hat(ind_in:ind_fin);
yp = yp(ind_in:ind_fin);
ts = ts-Tburn;

%% Plot dei risultati

figure
plot(tv, yp, 'LineWidth', 1.5, 'Color', 'Red')
title('stima del segnale glicemico CGM senza rumore')
xlabel('tempo [min]')
ylabel('segnale CGM [mg/dl]')
xlim([tv(1) tv(end)])

figure
plot(tv, u_hat, 'LineWidth', 1.5, 'Color', 'Green')
title('stima derivata prima del segnale glicemico CGM')
xlabel('tempo [min]')
ylabel('derivata prima CGM [mg/(dl·min)]')
xlim([tv(1) tv(end)])

figure
subplot(211)
hold on
plot(ts, ys, 'LineWidth', 1.5, 'Color', 'Blue')
plot(tv, yp, 'LineWidth', 1.5, 'Color', 'Red', 'LineStyle', '--')
hold off
title('segnale glicemico misurato vs riconvoluzione')
ylabel('segnale CGM [mg/dl]')
xlabel('tempo [min]')
xlim([tv(1) tv(end)])
legend('segnale originale', 'riconvoluzione')
subplot(212)
hold on
plot([tv(1) tv(end)], [-1 -1], 'color', 'black', 'linestyle', '--', 'LineWidth', 1)
plot([tv(1) tv(end)] , [1 1], 'color', 'black', 'linestyle', '--', 'LineWidth', 1)
plot([tv(1) tv(end)] , [0 0], 'color', 'black', 'linestyle', '-.', 'LineWidth', 1)
plot(ts, resn, 'ro', 'LineWidth', 1)
plot(ts, resn, 'r-', 'LineWidth', 1)
hold off
title('residui normalizzati')
xlabel('tempo [min]')
ylabel('(ys-yp)/sd [adimensionale]')
ylim([-3 3])
