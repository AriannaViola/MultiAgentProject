%% 
%  PROGETTO 3 - Consenso event-triggered per dinamiche lineari generali
%  Driver principale: esegue i tre esperimenti, produce figure e metriche.
%
%  Riferimenti:
%   [G] Garcia, Cao, Casbeer, Automatica 2014   (protocollo, Th.2, Th.3)
%   [I] Isidori, Lectures in Feedback Design, 5.4-5.5  (design del guadagno)
%   [D] Dimarogonas, Frazzoli, Johansson, IEEE TAC 2012  (single integrator)
%
%  File del progetto:
%   make_config.m  design_gain.m  sim_continuous.m  sim_event_triggered.m
%   disagreement.m  conv_time.m
% 
clear; clc; close all;

%%  Esperimento 1: oscillatore (LTI generale, marginalmente stabile)
cfg = make_config('oscillator');
g   = design_gain(cfg);                        % Riccati -> P, F, c
c1  = g.c + 0.5*g.c;                            % guadagno usato in continuo e ET
fprintf('[oscillatore] lambda2=%.3f  c=%.3f  c1=%.3f\n', g.lam2, g.c, c1);

[tc, Xc] = sim_continuous(cfg, g, c1);
et = sim_event_triggered(cfg, g, struct('rule','garcia','model','model','rho',0.5));
N = cfg.N; n = g.n;

dc = disagreement(Xc, N, n);
de = disagreement(et.X, N, n);
tol = 1e-2*dc(1);
fprintf('  trigger/agente: %s  (tot %d)\n', mat2str(et.counts(:).'), sum(et.counts));
fprintf('  t_conv  continuo=%.2fs  ET=%.2fs\n', ...
        conv_time(dc,tc,tol), conv_time(de,et.t,tol));

figure('Name','Convergenza');
semilogy(tc, dc, 'LineWidth',1.6); hold on;
semilogy(et.t, de, '--', 'LineWidth',1.6); grid on;
xlabel('t [s]'); ylabel('||disagreement||_F');
legend('continuous','event-triggered'); title('Convergence at consensus (oscillator)');

figure('Name','Stati & uscite');
subplot(1,2,1); plot(et.t, squeeze(et.X(1,:,:)).'); grid on;
xlabel('t [s]'); title('States x_{i,1}(t) (ET)');
subplot(1,2,2);
Yet = zeros(N, numel(et.t));
for i=1:N, Yet(i,:) = cfg.C*squeeze(et.X(:,i,:)); end
plot(et.t, Yet.'); grid on; xlabel('t [s]'); title('Outputs y_i = C x_i (ET)');

%% figura eventi: raster a barrette + zoom + istogramma inter-event
figure('Name','Eventi','Position',[100 100 1200 380]);

% pannello 1: raster completo (ogni barretta verticale = un evento)
subplot(1,3,1); hold on;
for i = 1:N
    te = et.trig{i}(:).';                       % istanti di evento dell'agente i
    if isempty(te), continue; end
    plot([te; te], repmat([i-0.35; i+0.35], 1, numel(te)), 'k-', 'LineWidth', 0.6);
end
grid on; ylim([0.5 N+0.5]); yticks(1:N);
xlim([0 cfg.T]); xlabel('t [s]'); ylabel('agent');
title('Trigger instants');

% pannello 2: zoom sui primi secondi (eventi ben separati)
tz = 3;                                          % finestra di zoom [s]
subplot(1,3,2); hold on;
for i = 1:N
    te = et.trig{i}(:).';
    te = te(te <= tz);
    if isempty(te), continue; end
    plot([te; te], repmat([i-0.35; i+0.35], 1, numel(te)), 'k-', 'LineWidth', 1.0);
end
grid on; ylim([0.5 N+0.5]); yticks(1:N);
xlim([0 tz]); xlabel('t [s]'); ylabel('agent');
title(sprintf('Zoom [0, %g] s: discrete events', tz));

% pannello 3: istogramma inter-event (massa lontana da 0 = no Zeno)
iet = [];
for i = 1:N
    iet = [iet, diff(et.trig{i}(:)).'];          %#ok<AGROW>  forza a RIGA
end
iet = iet(iet > 0);                              % scarta eventuali dt nulli/negativi

subplot(1,3,3); cla;
if ~isempty(iet)
    hH = histogram(iet, 25, 'FaceColor',[0.2 0.2 0.2], 'EdgeColor','w');
    hold on;                                     % <-- senza questo la riga cancella l'istogramma
    hL = xline(min(iet), 'r--', 'LineWidth', 1.5);   % tau_min > 0
    hold off;
    legend([hH hL], {'inter-event','\tau_{min} > 0'}, 'Location','northeast');
end
grid on;
xlabel('inter-event time [s]'); ylabel('count');
title('Inter-event distribution');

%% sweep su rho (il "sigma": curva di trade-off)
rhos = 0.1:0.1:0.9; tot = zeros(size(rhos)); tcv = zeros(size(rhos));
for q = 1:numel(rhos)
    e = sim_event_triggered(cfg, g, struct('rule','garcia','model','model','rho',rhos(q)));
    d = disagreement(e.X, N, n);
    tot(q) = sum(e.counts); tcv(q) = conv_time(d, e.t, tol);
end
figure('Name','Sweep rho');
yyaxis left;  plot(rhos, tot, '-o','LineWidth',1.5); ylabel('total triggers');
yyaxis right; plot(rhos, tcv, '-s','LineWidth',1.5); ylabel('convergence time [s]');
grid on; xlabel('\rho  (parametro \sigma)'); title('Trade-off communication / convergence');

%%  Esperimento 2: agenti instabili - model-based vs ZOH 
cfg2 = make_config('garcia');
P_paper = [ 4.8436  5.4783 -1.1082;       % P fornita nel paper [G]
            5.4783  7.0514 -1.4299;
           -1.1082 -1.4299  0.3778];
g2 = design_gain(cfg2, [], [], P_paper);
% stessa regola di trigger (Th.3 anti-Zeno), cambia SOLO il modello
em = sim_event_triggered(cfg2, g2, struct('rule','garcia_zeno','model','model','rho',0.5,'phi',0.3));
ez = sim_event_triggered(cfg2, g2, struct('rule','garcia_zeno','model','zoh','rho',0.5,'phi',0.3));
fprintf('\n[instabili] trasmissioni model-based=%s (tot %d)\n', mat2str(em.counts(:).'), sum(em.counts));
fprintf('            trasmissioni ZOH       =%s (tot %d)\n', mat2str(ez.counts(:).'), sum(ez.counts));

dm = disagreement(em.X, cfg2.N, g2.n);
dz = disagreement(ez.X, cfg2.N, g2.n);
figure('Name','Model vs ZOH');
subplot(1,2,1);
semilogy(em.t, dm, 'LineWidth',1.6); hold on; semilogy(ez.t, dz, '--', 'LineWidth',1.6);
grid on; xlabel('t [s]'); ylabel('||disagreement||_F'); legend('model-based','ZOH');
title('Disagreement (instable agents)');
subplot(1,2,2);
bar(1:cfg2.N, [em.counts(:) ez.counts(:)]); grid on;
xlabel('agent'); ylabel('n. transmissions'); legend('model-based','ZOH');
title('Communication for each agent');

%% Validazione: single integrator (Dimarogonas) 
cfg3 = make_config('single_integrator');
g3 = design_gain(cfg3);
es = sim_event_triggered(cfg3, g3, struct('rule','garcia','model','model','rho',0.5));
xbar0 = mean(cfg3.x0);
xfin  = squeeze(es.X(1,:,end));
fprintf('\n[single integrator] media iniziale=%.5f\n', xbar0);
fprintf('                    stati finali   =%s\n', mat2str(xfin,5));
fprintf('                    -> convergenza alla media iniziale: %s\n', ...
        ternary(max(abs(xfin-xbar0))<1e-3,'OK','NO'));

%% Esperimento 4: consenso d'uscita observer-based (sola misura y) 
cfg4 = make_config('oscillator');     % C = [1 0]: si misura solo la posizione
g4 = design_gain(cfg4);
Lobs = design_observer(cfg4, [-3 -4]);% osservatore piu' veloce del consenso
eo = sim_et_output(cfg4, g4, Lobs, struct('rho',0.5));
N4 = cfg4.N;
fprintf('\n[uscita] trigger/agente: %s  (tot %d)\n', mat2str(eo.counts(:).'), sum(eo.counts));
fprintf('         errore osservatore finale=%.2e\n', max(eo.XE(:,end)));

Yout = zeros(N4, numel(eo.t));
for i=1:N4, Yout(i,:) = cfg4.C*squeeze(eo.X(:,i,:)); end
do = disagreement(eo.X, N4, g4.n);
figure('Name','Consenso uscita');
subplot(1,3,1); plot(eo.t, Yout.'); grid on;
xlabel('t [s]'); title('$\mathbf{Outputs}\ \mathbf{y_i = C x_i}$','Interpreter','latex');

subplot(1,3,2); semilogy(eo.t, eo.XE.'); grid on;
xlabel('t [s]'); title('$\mathbf{Observer\ error}\ \|\mathbf{x}_i - \hat{\mathbf{x}}_i\|$','Interpreter','latex');

subplot(1,3,3); semilogy(eo.t, do); grid on;
xlabel('t [s]'); title('$\mathbf{States\ disagreement}\ \|\mathbf{x} - \bar{\mathbf{x}}\|_F$','Interpreter','latex');

%%  Esperimento 5: confronto COMPUTAZIONALE (comunicazione vs computazione) 
% (a) single-integrator: l'event-triggered riduce ENTRAMBI gli assi
cfgS = make_config('single_integrator'); gS = design_gain(cfgS);
stepsS = round(cfgS.T/1e-3); NS = cfgS.N;
es5 = sim_event_triggered(cfgS, gS, struct('rule','garcia','model','model','rho',0.5));
commS = [NS*stepsS, sum(es5.counts)];
compS = [NS*stepsS, sum(es5.ctrl_updates)];
fprintf('\n[computazione | single-integrator]\n');
fprintf('  continuo:        comunicazione=%6d  computazione=%6d\n', commS(1), compS(1));
fprintf('  event-triggered: comunicazione=%6d  computazione=%6d\n', commS(2), compS(2));

figure('Name','Computazione');
bar([commS; compS].'); set(gca,'YScale','log'); grid on;
set(gca,'XTickLabel',{'continuous','event-triggered'});
ylabel('log'); legend('communication','computation');
title('Communication and computation (single-integrator)');

% (b) oscillatore (LTI generale): confronto a tre - tabella a console
cfgO = make_config('oscillator'); gO = design_gain(cfgO);
stepsO = round(cfgO.T/1e-3); NO = cfgO.N;
e_mb5  = sim_event_triggered(cfgO, gO, struct('rule','garcia','model','model','rho',0.5));
e_zoh5 = sim_event_triggered(cfgO, gO, struct('rule','garcia','model','zoh','rho',0.5));
fprintf('[computazione | oscillatore]\n');
fprintf('  %-22s %10s %10s\n','controllore','comunic.','comput.');
fprintf('  %-22s %10d %10d\n','continuo', NO*stepsO, NO*stepsO);
fprintf('  %-22s %10d %10d\n','ET model-based', sum(e_mb5.counts), sum(e_mb5.ctrl_updates));
fprintf('  %-22s %10d %10d\n','ET controllo tenuto', sum(e_zoh5.counts), sum(e_zoh5.ctrl_updates));

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
