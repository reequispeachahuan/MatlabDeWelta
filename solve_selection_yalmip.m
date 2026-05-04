function [delta_value, objective_value] = solve_selection_yalmip(M_prime_full, Nr, alpha, Nt, lambda, sigma_x, H_eff_r_q_full)
%SOLVE_SELECTION_YALMIP Solve comparator-selection relaxation with YALMIP.
% Keeps the same model as the original CVX block:
%   maximize 1/2 * log_det(I + c * H' * diag(delta) * H)
%   s.t. delta(1:2Nr)=1, 0<=delta<=1, sum(delta)=2Nr+alpha.

delta = sdpvar(M_prime_full, 1);
objective_matrix = eye(2 * Nt) + 1 / lambda * (sigma_x^2 / 2) * (H_eff_r_q_full' * diag(delta) * H_eff_r_q_full);
constraints = [delta(1:2*Nr) == 1, 0 <= delta <= 1, sum(delta) == 2 * Nr + alpha];

% Avoid LMILAB (unsupported/poor SDP behavior in YALMIP). Try better solvers.
candidate_solvers = {'mosek', 'sdpt3', 'sedumi', 'sdpnal'};
diagnostics = [];
for k = 1:numel(candidate_solvers)
    ops = sdpsettings('solver', candidate_solvers{k}, 'verbose', 0);
    diagnostics = optimize(constraints, -0.5 * logdet(objective_matrix), ops);
    if diagnostics.problem == 0
        break;
    end
end

if isempty(diagnostics) || diagnostics.problem ~= 0
    error(sprintf(['YALMIP optimization failed (all preferred SDP solvers failed). ' ...
           'Install/use MOSEK, SDPT3, or SeDuMi. Last code=%d, info=%s'], diagnostics.problem, diagnostics.info));
end

delta_value = value(delta);
objective_value = 0.5 * log(det(value(objective_matrix)));
end
