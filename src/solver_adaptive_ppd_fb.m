function [xsol, L1_v, L1_vp, L2_v, L2_vp, delta_v, sol_v, no_sub_itr_v, v1, v2, sol_reweight_v] = solver_adaptive_ppd_fb(y, epsilont, epsilonts, epsilon, epsilons, A, At, T, pU, W, Psi, Psit, Psiw, Psitw, param)
%
% [xsol, L1_v, L1_vp, L2_v, L2_vp] = solver_adaptive_ppd_fb(y, epsilon, A, At, T, W, Psi, Psit, param) solves:
%
%   min ||Psit x||_1   s.t.  ||y-A x||_2 <= epsilon and x>=0
%
%
% y contains the measurements. A is the forward measurement operator and
% At the associated adjoint operator. Psit is a sparfying transform and Psi
% its adjoint. PARAM a Matlab structure containing the following fields:
%
%   General parameters:
%
%   - verbose: 0 no log, 1 print main steps, 2 print all steps.
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - rel_obj: minimum relative change of the objective value (default:
%   1e-4)
%       The algorithm stops if
%           | ||x(t)||_1 - ||x(t-1)||_1 | / ||x(t)||_1 < rel_obj,
%       where x(t) is the estimate of the solution at iteration t.
%
%   - param.weights: weights (default = 1) for a weighted L1-norm defined
%       as sum_i{weights_i.*abs(x_i)}
%
% The impelmentation simulates a distributed setup for parallel processing
%
% Authors: Alex Onose, Rafael Carrillo, Arwa Dabbech
% E-mail: a.onose@hw.ac.uk, rafael.carrillo@epfl.ch, a.dabbech@hw.ac.uk

% number of nodes
R = length(y);
P = length(Psit);

if W{1} ~= ':'
    % oversampling vectorized data length
    No = size(W{1}, 1);
else
    No = size(T{1}' * y{1} , 1);
end

% number of pixels
[Ny, Nx] = size(At(zeros(No, 1)));


%% optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'nu1')
    param.nu1 = ones(P, 1);
else
    if numel(param.nu1) == 1
        param.nu1 = ones(P, 1) * param.nu1;
    end
end
if ~isfield(param, 'nu2')
    param.nu2 = zeros(R, 1);
    % maximum eigenvalue of operato A^T A
    for q = 1:R
        Tw = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
        Tw(:, W{q}) = T{q};
        fprintf('\nComputing operator norm: block %i \n', q)
        param.nu2(q) = op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
        clear Tw;
    end
    fprintf('\n');
else
    if numel(param.nu2) == 1
        param.nu2 = ones(R, 1) * param.nu2;
    end
end
if ~isfield(param, 'sigma1'), param.sigma1 = 1./param.nu1; end
if ~isfield(param, 'sigma2'), param.sigma2 = 1./param.nu2; end
if ~isfield(param, 'gamma'), param.gamma = 1e-3; end
if ~isfield(param, 'tau'), param.tau = 0.49; end
if ~isfield(param, 'weights')
    param.weights = cell(P, 1);
    for k = 1:P
        param.weights{k} = ones(size(Psit{k}(At(zeros(No, 1))), 1), 1);
    end
else
    if ~iscell(param.weights)
        weights = param.weights;
        param.weights = cell(P, 1);
        for k = 1:P
            param.weights{k} = weights;
        end
    end
end
if isfield(param, 'initsol')
    xsol = param.initsol;
else
    % start from zero solution
    xsol = zeros(Ny, Nx);
end
if isfield(param, 'initv1')
    norm1 = cell(P, 1);
    % r1 = cell(P, 1);
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    if iscell(param.initv1)
        % initial dual variables
        v1 = param.initv1;
    else
        for k = 1:P
            % initial L1 part variable
            v1{k} = param.initv1;
        end
    end
    for k = 1:P
        % initial L1 descent step
        u1{k} = [];%zeros(size(Psi{k}(v1{k})));
    end
else
    norm1 = cell(P, 1);
    r1 = cell(P, 1);
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    for k = 1:P
        % start from zero solution
        v1{k} = zeros(size(Psit{k}(xsol)));
        
        % initial L1 descent step
        u1{k} = [];%zeros(size(Psi{k}(v1{k})));
    end
end
if isfield(param, 'initv2')
    block_l2norm = cell(R, 1);
    r2 = cell(R, 1);
    u2 = cell(R, 1);
    v2 = cell(R, 1);
    % initial L2 ball variables
    % initial dual variables
    if iscell(param.initv2)
        % initial dual variables
        v2 = param.initv2;
    else
        for q = 1:R
            % initial L1 part variable
            v2{q} = param.initv2;
        end
    end
    for q = 1:R
        % initial L1 part descent step
        u2{q} = [];%zeros(size(T{q}, 1), 1);
    end
    %vy2 = v2;
else
    block_l2norm = cell(R, 1);
    r2 = cell(R, 1);
    u2 = cell(R, 1);
    v2 = cell(R, 1);
    for q = 1:R
        % initial L1 part variable
        v2{q} = zeros(length(y{q}), 1);
        % initial L1 part descent step
        u2{q} = [];% zeros(size(T{q}, 2), 1);
    end
end
if ~isfield(param, 'lambda0'), param.lambda0 = 1; end
if ~isfield(param, 'lambda1'), param.lambda1 = 1; end
if ~isfield(param, 'lambda2'), param.lambda2 = 1; end
if ~isfield(param, 'sol_steps')
    param.sol_steps = inf;
else
    if param.sol_steps(end) ~= inf
        param.sol_steps = [param.sol_steps inf];
    end
end
if ~isfield(param, 'reweight_steps')
    param.reweight_steps = inf;
else
    if param.reweight_steps(end) ~= inf
        param.reweight_steps = [param.reweight_steps inf];
    end
end
if ~isfield(param, 'best_bound_steps')
    param.best_bound_steps = inf;
else
    if param.best_bound_steps(end) ~= inf
        param.best_bound_steps = [param.best_bound_steps inf];
    end
end
if ~isfield(param, 'omega1'), param.omega1 = 1; end
if ~isfield(param, 'omega2'), param.omega2 = 1; end

if ~isfield(param, 'global_stop_bound'), param.global_stop_bound = 1; end

if ~isfield(param, 'use_proj_elipse_fb')
    param.use_proj_elipse_fb = 1;
end
if ~isfield(param, 'elipse_proj_max_iter')
    param.elipse_proj_max_iter = 2000;
end
if ~isfield(param, 'elipse_proj_min_iter')
    param.elipse_proj_min_iter = 1;
end
if ~isfield(param, 'elipse_proj_eps')
    param.elipse_proj_eps = 1e-8;
end
if ~isfield(param, 'use_reweight_steps')
    param.use_reweight_steps = 0;
end
if ~isfield(param, 'use_reweight_eps')
    param.use_reweight_eps = 0;
end
if ~isfield(param, 'reweightBound')
    param.reweightBound = 0;
end

%-----------------------%
if ~isfield(param, 'use_adapt_bound_eps') %
    param.use_adapt_bound_eps = 0;
end
if ~isfield(param, 'hard_thres')  %
    param.hard_thres = 0;
end
%-----------------------%
%% set up log variables
L1_v = zeros(param.max_iter, 1);
L1_vp = zeros(param.max_iter, P);
L2_v = zeros(param.max_iter, 1);
L2_vp = zeros(param.max_iter, R);
no_sub_itr_v = cell(param.max_iter, 1);

delta_v = zeros(param.max_iter, 1);

sol_steps = param.sol_steps;
sol_step_count = 1;

reweight_steps = param.reweight_steps;
reweight_step_count = 1;
reweight_last_step_iter = 1;
reweight_step_count_since_best_bound_search = 0;


sol_v =[];% zeros(length(sol_steps)-1, Ny, Nx);

sol_reweight_v = [];% zeros(0, Ny, Nx);


%% useful functions for the projection
% thresholding negative values
%hardt = @(z) max(real(z), min(-param.im0, 0));
hardt = @(z) max(real(z), 0);

%soft thresholding operator
soft = @(z, T) sign(z) .* max(abs(z)-T, 0);

phi_= 0.5*(sqrt(5)-1);
%% initialization

% initial primal gradient like step
%g1 = zeros(size(xsol));
%g2 = zeros(size(xsol));

% solution flag: 0 - max iteration reached; 1 - solution found
flag = 0;

%% store useful variables
% step size for the dual variables
sigma1 = param.sigma1;
sigma2 = param.sigma2;

% step size primal
tau = param.tau;

% relaxation parameters
lambda0 = param.lambda0;
lambda1 = param.lambda1;
lambda2 = param.lambda2;


% weights
weights = param.weights;
param.weights = [];

% omega sizes
omega1 = param.omega1;
omega2 = param.omega2;

gamma = param.gamma;

reweight_alpha = param.reweight_alpha;
reweight_alpha_ff = param.reweight_alpha_ff;
reweightBound = param.reweightBound;

[proj] = solver_find_elipse_point(y, pU, A, T, xsol, v2, W, epsilont, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);

A = afclean(A);
At = afclean(At);

for k = 1:P
    Psi{k} = afclean(Psi{k});
    Psit{k} = afclean(Psit{k});
end
Psiw = [];
Psiwt = [];
util_create_pool(P);

%% main loop: sequential + simulated parallel
% xsol     - current solution estimate
% prev_sol - previous solution estimate
% r1       - L1 part of criterion: r1 = Psi' * sol
% v1       - L1 dual variable
% r2       - L2 part of criterion: r2 = T * A * sol
% v2       - L2 dual variable


% epsilon update counter
t_update=cell(R,1);
count_eps_update = 0; %
for q = 1:R %
    t_update{q}=1;
end


% initial primal gradient like step
%g2 = zeros(size(xsol));

%%%%%%%%%%%%%%%%%%%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear u2 vy2 r1 r2;


g1 = zeros(size(xsol));
use_proj_elipse_fb =  param.use_proj_elipse_fb;
elipse_proj_max_iter = param.elipse_proj_max_iter;
elipse_proj_min_iter =param.elipse_proj_min_iter;
elipse_proj_eps = param.elipse_proj_eps;


for t = 1:param.max_iter
    tm = tic;
    %% primal update
    
    
    prev_xsol = xsol;
    xsol = (hardt(xsol - tau * g1));
  
    norm_prevsol = norm(prev_xsol(:));
  
    % solution relative change
    if (norm_prevsol == 0)
        rel_sol_norm_change = 1;
    else
        rel_sol_norm_change = norm(xsol(:) - prev_xsol(:))/norm_prevsol;
    end
    
    prev_xsol = 2*xsol - prev_xsol;
    
    
    %% L1 prox update: dual variables update
    % parallel for all bases
    tic;
    for k = 1:P
        
        f(k) = parfeval(@run_par_waverec, 3, v1{k}, Psit{k}, Psi{k}, prev_xsol, gamma, weights{k}, sigma1(k), lambda1);
        
    end
    toc;
    
    
    %% L2 ball projection update: dual variables update
    
    % non gridded measurements of current solution
    ns = A(prev_xsol);
    clear prev_xsol;
    
    
    no_sub_itr = cell(R, 1);
    
    ns = ns(:);
    uu = zeros(No, 1);
    
    % parallel for all R blocks
    
    for q = 1:R
        
        %tic
        residual_data_block =T{q}*ns;
        block_l2norm{q} = norm(residual_data_block-y{q});
        
        [proj{q}, no_sub_itr{q}] = solver_proj_elipse_fb(1 ./ pU{q} .* v2{q}, residual_data_block, y{q}, pU{q}, epsilont{q}, proj{q}, elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
        v2{q} = v2{q} + pU{q} .* (residual_data_block -  proj{q});
        
        
        uu = uu +T{q}'*v2{q} ; %  uu(W{q}) = uu(W{q}) + u2{q};
        %toc
        
    end
    
    
    % ADAPTIVE bound update on each block
    if( param.use_adapt_bound_eps ==1)
        for q = 1:R
            if (block_l2norm{q} < epsilonts{q}) && (block_l2norm{q}> (1-param.adapt_bound_tol)*epsilont{q})
                t_update{q} = param.max_iter;
                
            elseif   ((t>t_update{q}+param.adapt_bound_steps) || (t_update{q} == param.max_iter )) && ...
                    (block_l2norm{q}< (1-param.adapt_bound_tol)*epsilont{q}) &&  ...
                    (rel_sol_norm_change <  param.adapt_bound_rel_obj)
                
                t_update{q} = t;
                
                epsilont{q} = block_l2norm{q} + (-block_l2norm{q} + epsilont{q})*(1-phi_);
                epsilonts{q}= (1+param.adapt_bound_tol)*epsilont{q};
                epsilon = norm(cell2mat(epsilont));
                epsilons= (1+param.adapt_bound_tol)*epsilon;
                count_eps_update = count_eps_update +1;
                
                if param.verbose>=1
                    fprintf('========================================================\n')
                    fprintf('Updated bound of block %0.0f,  down to = %0.1f\n',  epsilont{q});
                    fprintf('========================================================\n')
                    
                end
                
                
            end
            
            if (block_l2norm{q} > epsilonts{q})  && ...
                    ((( (t>t_update{q}+param.adapt_bound_steps) && (rel_sol_norm_change < param.adapt_bound_rel_obj)) || ...
                    ((t_update{q} == param.max_iter ) && (rel_sol_norm_change <  param.reweight_rel_obj))    ) ||...
                    (t ==param.adapt_bound_start))
                
                t_update{q} = t;
                epsilont{q} =epsilont{q} + (block_l2norm{q} - epsilont{q})*phi_;% (norm2{q} + epsilont{q})/2;
                epsilonts{q}= (1+param.adapt_bound_tol)*epsilont{q};
                epsilon = norm(cell2mat(epsilont));
                epsilons= (1+param.adapt_bound_tol)*epsilon;
                count_eps_update = count_eps_update +1;
                
                
                if param.verbose>=1
                    fprintf('=========================================================\n')
                    fprintf('Updated bound of block %0.0f, up to = %0.1f\n', q, epsilont{q});
                    fprintf('=========================================================\n')
                    
                end
                
                
            end
        end
        
    end
    % ADAPTIVE bound update on each block
    
    %% update the primal gradient
    
    g1 = (sigma2(1) ./sigma1(1))*(At(uu));
    clear uu;
    for k = 1:P
        [idx, v1_, u1_, norm1_] = fetchNext(f);
        v1{idx} = v1_;
        norm1{idx} = norm1_;
        g1 = g1 + u1_;
    end
    
    clear u1 u1_ v1_ idx;
    
    %spmd
    g1 = sigma1(1) * g1;
    %end
    tm = toc(tm);
    %% stopping criterion and logs
    
    % log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t);
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(block_l2norm)));
        fprintf(' Global residual bound         = %e\n', epsilon);
        fprintf(' Distributed residual L2 ball  = %e\n', norm(cell2mat(epsilont)));
        fprintf(' Distributed residual L2 bound = %e\n', norm(cell2mat(epsilonts)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        if (param.verbose >= 2)
            for q = 1:R
                fprintf('   Residual %i                     = %e\n', q, block_l2norm{q});
                fprintf('   Residual L2 ball %i             = %e\n', q, epsilont{q});
                fprintf('   Residual L2 bound %i            = %e\n\n', q, epsilonts{q});
            end
        end
        fprintf('Time for iteration %i: %3.3f\n\n\n',t, tm);
        
    end
    if (param.verbose <= 0.5)
        fprintf('.\n');fprintf('\b');
        if mod(t, 50) == 0
            fprintf('\n');
        end
    end
    if (param.verbose >= 0.5)
        L1_v(t) = sum(cell2mat(norm1));
        L2_v(t) = norm(cell2mat(block_l2norm));
        no_sub_itr_v{t} = no_sub_itr;
        for q = 1:R
            L2_vp(t, q) = block_l2norm{q};
        end
        for k = 1:P
            L1_vp(t, k) = norm1{k};
        end
        delta_v(t) = rel_sol_norm_change;
        
    end
    
    if (param.use_reweight_steps || param.use_reweight_eps) && ( t < param.reweight_max_reweight_itr)
        
        if (param.use_reweight_steps && t == reweight_steps(reweight_step_count)) || ...
                (param.use_reweight_eps && ...
                norm(cell2mat(block_l2norm)) <= norm(cell2mat(epsilonts)) && ...
                param.reweight_min_steps_rel_obj < t - reweight_last_step_iter && ...
                rel_sol_norm_change < param.reweight_rel_obj)
            
            % parallel for all bases
            for k = 1:P
                d_val = abs(Psit{k}(xsol));
                weights{k} = reweight_alpha(k) ./ (reweight_alpha(k) + (d_val));
            end
            
            reweight_alpha = max(reweight_alpha_ff .* reweight_alpha,reweightBound);
            sigma1(:) = 1/op_norm_wave_par(weights, Psit, Psi, [Ny, Nx], 1e-8, 200, 0, P);
            
            fprintf('\n\n\n\n\n\n\n Performed reweight no %d \n\n\n\n\n', reweight_step_count);
            
            if reweight_step_count > param.total_reweights
                param.reweight_max_reweight_itr = t+1;
                fprintf('\n\n\n\n\n\n\n No more reweights \n\n\n\n\n');
                
            end
            sol_v{reweight_step_count} =xsol;
            reweight_step_count = reweight_step_count + 1;
            reweight_last_step_iter = t;
            reweight_step_count_since_best_bound_search = reweight_step_count_since_best_bound_search + 1;
            
        end
    end
    
    
    
    
    
    % global stopping criteria
    if rel_sol_norm_change < param.rel_obj && ...
            ((param.global_stop_bound && norm(cell2mat(block_l2norm)) <= norm(cell2mat(epsilonts))) || ...
            (~param.global_stop_bound && prod(cell2mat(block_l2norm) <= cell2mat(epsilonts))))
        flag = 1;
        break;
    end
end


% final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('\nSolution found\n');
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(block_l2norm)));
        fprintf(' Global residual bound         = %e\n', epsilon);
        fprintf(' Distributed residual L2 ball  = %e\n', norm(cell2mat(epsilont)));
        fprintf(' Distributed residual L2 bound = %e\n', norm(cell2mat(epsilonts)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        for q = 1:R
            fprintf('   Residual %i                     = %e\n', q, block_l2norm{q});
            fprintf('   Residual L2 ball %i             = %e\n', q, epsilont{q});
            fprintf('   Residual L2 bound %i            = %e\n\n', q, epsilonts{q});
        end
    else
        fprintf('\nMaximum number of iterations reached\n');
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(block_l2norm)));
        fprintf(' Global residual bound         = %e\n', epsilon);
        fprintf(' Distributed residual L2 ball  = %e\n', norm(cell2mat(epsilont)));
        fprintf(' Distributed residual L2 bound = %e\n', norm(cell2mat(epsilonts)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        for q = 1:R
            fprintf('   Residual %i                     = %e\n', q, block_l2norm{q});
            fprintf('   Residual L2 ball %i             = %e\n', q, epsilont{q});
            fprintf('   Residual L2 bound %i            = %e\n\n', q, epsilonts{q});
        end
    end
end
fprintf('\n');

xsol = hardt(real(xsol));


% trim the log vectors to the size of the actual iterations performed
if (param.verbose >= 0.5)
    L1_v = L1_v(1:t);
    L1_vp = L1_vp(1:t, :);
    L2_v = L2_v(1:t);
    L2_vp = L2_vp(1:t, :);
    delta_v = delta_v(1:t);
    no_sub_itr_v = {no_sub_itr_v{1:t}}';
end

end


function [v1_, u1_, norm1_] = run_par_waverec(v1_, Psit, Psi, prev_xsol, gamma, weights_, sigma1_, lambda1)

r1_ = Psit(prev_xsol);
z = v1_ + r1_;
T = (gamma  / sigma1_) * weights_;
v1_ = z - sign(z) .* max(abs(z)-T, 0);
u1_ =  Psi(v1_);%
%u1_ =Psi(weights_ .* v1_);

% local L1 norm of current solution
norm1_ = sum(abs(r1_));
end
