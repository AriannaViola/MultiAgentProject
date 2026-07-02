function res = sim_event_triggered(cfg, g, opt)
%SIM_EVENT_TRIGGERED  Consenso event-triggered per dinamiche lineari generali.
%
%   res = sim_event_triggered(cfg, g, opt)
%
%   opt (struct, tutti i campi opzionali):
%     .rule   'garcia'      -> Th.2 [G] eq.(10)-(12), consenso asintotico
%             'garcia_zeno' -> Th.3 [G] eq.(30), consenso limitato, anti-Zeno
%             'sigma'       -> regola dell'enunciato/[D]:  ||e_i|| > sigma*||y_i||
%     .model  'model' -> stimatore dinamico  dot(y_i)=A y_i   (Garcia)
%             'zoh'   -> zero-order-hold      dot(y_i)=0
%     .c2,.b,.rho,.phi,.sigma,.dt  parametri (default sotto)
%
%   Protocollo: u_i = c1 F z_i, F=-B'P, z_i = sum_{j in N_i}(y_i - y_j),
%   c1 = c + c2. All'evento dell'agente i: broadcast x_i, reset y_i <- x_i.
%   e_i = y_i - x_i  e' l'errore di modello.
if nargin < 3, opt = struct(); end
opt = set_defaults(opt, g);
A = cfg.A; B = cfg.B; Adj = cfg.Adj; N = cfg.N; n = g.n; m = g.m;
F = g.F; M = g.M; c = g.c;  c2 = opt.c2; c1 = c + c2;  b = opt.b;
dt = opt.dt; steps = round(cfg.T/dt);

nbr = arrayfun(@(i) find(Adj(i,:)>0), 1:N, 'uni', 0);
deg = cellfun(@numel, nbr);

x = cfg.x0.'; y = cfg.x0.';            % (n x N), e_i(0)=0
X = zeros(n, N, steps+1); X(:,:,1) = x;
trig = cell(N,1);
u_prev = nan(m,N); ctrl_updates = zeros(N,1);   % n. ricalcoli controllo (computazione)

for s = 1:steps
    tnow = (s-1)*dt;
    %  verifica trigger per agente 
    for i = 1:N
        di = deg(i); if di==0, continue; end
        zi = di*y(:,i) - sum(y(:,nbr{i}),2);
        ei = y(:,i) - x(:,i);
        switch opt.rule
            case {'garcia','garcia_zeno'}
                Ki  = 2*c*di^2*(1+b) + ((c2-c)/b)*di + c*di*(N-1)*(b+3/b);
                psi = 2*(c2-c)*di*(zi.'*M*ei) + Ki*(ei.'*M*ei);     % (11)
                si  = 2*c2 - b*di*(c2-c);                          % (12)
                thr = opt.rho*si*(zi.'*M*zi);                      % (10) rhs
                if strcmp(opt.rule,'garcia_zeno'), thr = thr + opt.phi; end  % (30)
                fire = psi > thr;
            case 'sigma'
                fire = norm(ei) > opt.sigma*max(norm(y(:,i)),1e-9);
            otherwise
                error('rule sconosciuta');
        end
        if fire
            y(:,i) = x(:,i);           % broadcast + reset modello
            trig{i}(end+1) = tnow;    
        end
    end
    %  conteggio aggiornamenti del controllo (computazione) 
    u_nom = zeros(m,N);
    for i = 1:N
        di = deg(i);
        if di==0, zc = zeros(n,1); else, zc = di*y(:,i)-sum(y(:,nbr{i}),2); end
        u_nom(:,i) = c1*F*zc;
    end
    if s==1, ctrl_updates = ctrl_updates + 1;
    else,    ctrl_updates = ctrl_updates + (vecnorm(u_nom-u_prev).' > 1e-9); end
    u_prev = u_nom;
    % passo RK4 su (x,y) 
    [k1x,k1y] = deriv(x,y);
    [k2x,k2y] = deriv(x+dt/2*k1x, y+dt/2*k1y);
    [k3x,k3y] = deriv(x+dt/2*k2x, y+dt/2*k2y);
    [k4x,k4y] = deriv(x+dt*k3x,   y+dt*k3y);
    x = x + dt/6*(k1x+2*k2x+2*k3x+k4x);
    y = y + dt/6*(k1y+2*k2y+2*k3y+k4y);
    X(:,:,s+1) = x;
end

res.t = (0:steps)*dt;
res.X = X;                              % (n x N x T)
res.trig = trig;
res.counts = cellfun(@numel, trig);
res.ctrl_updates = ctrl_updates;
res.c1 = c1; res.c2 = c2; res.opt = opt;

    % dinamica (controllo dai modelli) 
    function [dx,dy] = deriv(xx,yy)
        dx = zeros(n,N); dy = zeros(n,N);
        for ii = 1:N
            dii = deg(ii); if dii==0, zi2 = zeros(n,1);
            else, zi2 = dii*yy(:,ii) - sum(yy(:,nbr{ii}),2); end
            dx(:,ii) = A*xx(:,ii) + B*(c1*F*zi2);
            if strcmp(opt.model,'model'), dy(:,ii) = A*yy(:,ii); else, dy(:,ii) = 0; end
        end
    end
end

function opt = set_defaults(opt, g)
def = struct('rule','garcia','model','model','c2',0.5*g.c,'b',1.0, ...
             'rho',0.5,'phi',0.3,'sigma',0.3,'dt',1e-3);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opt,fn{k}) || isempty(opt.(fn{k})), opt.(fn{k}) = def.(fn{k}); end
end
end
