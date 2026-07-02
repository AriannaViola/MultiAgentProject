function cfg = make_config(name)
%MAKE_CONFIG  Preset di problema (grafo, dinamica LTI, condizioni iniziali).
%   cfg = make_config('oscillator')        oscillatore armonico, 5 agenti
%   cfg = make_config('garcia')            esempio instabile del paper [G], 3 agenti
%   cfg = make_config('single_integrator') validazione Dimarogonas [D], 5 agenti
%
%   Campi di cfg: A,B,C, N, Adj, L, x0, T
rng(0);
switch name
    case 'oscillator'
        A = [0 1; -1 0];  B = [0; 1];  C = [1 0];      % uscita = posizione
        N = 5;
        edges = [1 2; 2 3; 3 4; 4 5; 5 1; 1 3];        % ciclo-5 + corda
        x0 = 4*randn(N,2);  T = 20;
    case 'garcia'
        A = [0.48 0.29 -0.30; 0.13 0.23 0; 0 -1.20 -1.00];
        B = [2 0; -1.5 1; 0 1];  C = eye(3);
        N = 3;  edges = [1 2; 2 3];                     % path 1-2-3
        x0 = 4*randn(N,3);  T = 12;
    case 'single_integrator'
        A = 0; B = 1; C = 1;
        N = 5;  edges = [1 2; 2 3; 3 4; 4 5; 5 1; 1 3];
        x0 = 5*randn(N,1);  T = 15;
    otherwise
        error('preset sconosciuto: %s', name);
end
Adj = zeros(N);
for k = 1:size(edges,1)
    Adj(edges(k,1),edges(k,2)) = 1; Adj(edges(k,2),edges(k,1)) = 1;
end
L = diag(sum(Adj,2)) - Adj;
cfg = struct('name',name,'A',A,'B',B,'C',C,'N',N,'Adj',Adj,'L',L,'x0',x0,'T',T);
end
