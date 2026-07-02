function [t, X] = sim_continuous(cfg, g, c1, dt)
%SIM_CONTINUOUS  Baseline continua: stessa legge di controllo dell'event-
%   triggered ma con comunicazione continua (errore di modello e_i == 0).
%   Dinamica impilata:  dot x = [ I (x) A  -  c1 (L (x) B B' P) ] x .
if nargin < 4, dt = 1e-3; end
A = cfg.A; B = cfg.B; L = cfg.L; N = cfg.N; n = g.n;
Acl = kron(eye(N),A) - c1*kron(L, B*B.'*g.P);
steps = round(cfg.T/dt);
x = cfg.x0.'; x = x(:);                % col(x_1,...,x_N)
X = zeros(N*n, steps+1); X(:,1) = x;
for s = 1:steps
    k1 = Acl*x; k2 = Acl*(x+dt/2*k1); k3 = Acl*(x+dt/2*k2); k4 = Acl*(x+dt*k3);
    x = x + dt/6*(k1+2*k2+2*k3+k4); X(:,s+1) = x;
end
t = (0:steps)*dt;
end
