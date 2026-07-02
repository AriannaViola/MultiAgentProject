function g = design_gain(cfg, theta, q, P_override)
%DESIGN_GAIN  Progetto del guadagno via Riccati ([G] eq.6-8 / [I] 5.37).
%   Risolve  (A+theta I)'P + P(A+theta I) - 2 P B B' P + q I = 0
%   pone  F = -B'P,  c = 1/lambda2(L).  Ritorna struct g con P,F,M,c,lam2,n,m.
%
%   P_override: usa una P data invece di risolvere la Riccati
%   (serve a riprodurre esattamente l'esempio del paper [G]).
if nargin < 2 || isempty(theta), theta = 0.5; end
if nargin < 3 || isempty(q),     q = 1.0;     end
A = cfg.A; B = cfg.B; [n,m] = size(B);

if nargin >= 4 && ~isempty(P_override)
    P = P_override;
else
    % care risolve A_'P + P A_ - P B R^{-1} B' P + Q = 0 ; con R=0.5 I -> 2 P B B' P
    P = care(A + theta*eye(n), B, q*eye(n), 0.5*eye(m));
end
F = -B.'*P;
M = P*B*B.'*P;                         % ricorre nelle soglie
ev = sort(eig(cfg.L));
lam2 = ev(2);                          % autovalore di Fiedler
g = struct('P',P,'F',F,'M',M,'lam2',lam2,'c',1/lam2,'n',n,'m',m);
end
