function res = sim_et_output(cfg, g, Lobs, opt)
%SIM_ET_OUTPUT  Consenso d'uscita event-triggered con sola misura y_i = C x_i.
%
%   Struttura: ogni agente misura solo y_i e ricostruisce x_i con un osservatore
%   locale; il consenso event-triggered model-based (Garcia) viene
%   eseguito sulle STIME xhat_i. Si trasmette/modella xhat_i, e il trigger usa
%   l'errore stima-modello e_i = chi_i - xhat_i. Per C = I si riduce al caso a
%   stato pieno di sim_event_triggered.
%
%   Stati per agente:  x_i (vero), xhat_i (stima), chi_i (modello condiviso).
%     dot x_i    = A x_i + B u_i
%     dot xhat_i = A xhat_i + B u_i + Lobs (y_i - C xhat_i)
%     dot chi_i  = A chi_i                       (reset chi_i <- xhat_i all'evento)
%     u_i = c1 F z_i,  z_i = sum_j (chi_i - chi_j)
if nargin < 4, opt = struct(); end
def = struct('rho',0.5,'c2',0.5*g.c,'b',1.0,'dt',1e-3,'xhat_offset',2.0);
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(opt,fn{k})||isempty(opt.(fn{k})), opt.(fn{k})=def.(fn{k}); end
end
A=cfg.A; B=cfg.B; C=cfg.C; Adj=cfg.Adj; N=cfg.N; n=g.n; m=g.m;
F=g.F; M=g.M; c=g.c; c2=opt.c2; c1=c+c2; b=opt.b; dt=opt.dt;
steps = round(cfg.T/dt);
nbr = arrayfun(@(i) find(Adj(i,:)>0), 1:N, 'uni', 0);
deg = cellfun(@numel, nbr);

x    = cfg.x0.';
xhat = cfg.x0.' + opt.xhat_offset*randn(n,N);   % stima iniziale errata
chi  = xhat;                                     % e_i(0)=0
X  = zeros(n,N,steps+1); X(:,:,1)=x;
XE = zeros(N,steps+1);   XE(:,1)=vecnorm(x-xhat);
trig = cell(N,1);

for s=1:steps
    tnow=(s-1)*dt;
    % --- trigger su e_i = chi_i - xhat_i (regola Garcia) ---
    for i=1:N
        di=deg(i); if di==0, continue; end
        zi = di*chi(:,i) - sum(chi(:,nbr{i}),2);
        ei = chi(:,i) - xhat(:,i);
        Ki = 2*c*di^2*(1+b) + ((c2-c)/b)*di + c*di*(N-1)*(b+3/b);
        psi = 2*(c2-c)*di*(zi.'*M*ei) + Ki*(ei.'*M*ei);
        si  = 2*c2 - b*di*(c2-c);
        if psi > opt.rho*si*(zi.'*M*zi)
            chi(:,i) = xhat(:,i);          % broadcast stima + reset modello
            trig{i}(end+1) = tnow;         
        end
    end
    % RK4 su (x, xhat, chi) 
    u = ctrl(chi);
    [a1,b1,c1k] = deriv(x,         xhat,         chi);
    [a2,b2,c2k] = deriv(x+dt/2*a1, xhat+dt/2*b1, chi+dt/2*c1k);
    [a3,b3,c3k] = deriv(x+dt/2*a2, xhat+dt/2*b2, chi+dt/2*c2k);
    [a4,b4,c4k] = deriv(x+dt*a3,   xhat+dt*b3,   chi+dt*c3k);
    x    = x    + dt/6*(a1+2*a2+2*a3+a4);
    xhat = xhat + dt/6*(b1+2*b2+2*b3+b4);
    chi  = chi  + dt/6*(c1k+2*c2k+2*c3k+c4k);
    X(:,:,s+1)=x; XE(:,s+1)=vecnorm(x-xhat);
end
res.t=(0:steps)*dt; res.X=X; res.XE=XE; res.trig=trig;
res.counts=cellfun(@numel,trig); res.Lobs=Lobs;

    function u = ctrl(ch)
        u=zeros(m,N);
        for ii=1:N
            dii=deg(ii);
            if dii==0, zi2=zeros(n,1); else, zi2=dii*ch(:,ii)-sum(ch(:,nbr{ii}),2); end
            u(:,ii)=c1*F*zi2;
        end
    end
    function [dx,dxh,dch]=deriv(xx,xh,ch)
        uu = ctrl(ch);
        dx=zeros(n,N); dxh=zeros(n,N); dch=zeros(n,N);
        for ii=1:N
            yi=C*xx(:,ii);
            dx(:,ii)  = A*xx(:,ii) + B*uu(:,ii);
            dxh(:,ii) = A*xh(:,ii) + B*uu(:,ii) + Lobs*(yi - C*xh(:,ii));
            dch(:,ii) = A*ch(:,ii);
        end
    end
end
