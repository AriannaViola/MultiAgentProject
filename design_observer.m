function Lobs = design_observer(cfg, poles)
%DESIGN_OBSERVER  Guadagno d'osservatore L tale che A - L C abbia gli autovalori
%   'poles' (richiede (A,C) osservabile). Struttura observer-based di Isidori 5.5.
A = cfg.A; C = cfg.C;
K = place(A.', C.', poles);     % colloca eig(A' - C' K) in 'poles'
Lobs = K.';                     % => eig(A - L C) = poles
end
