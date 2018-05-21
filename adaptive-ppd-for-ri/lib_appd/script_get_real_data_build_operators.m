
% script that loads or generates the input data for the test

        fprintf('Using real data ... \n\n');
        
        Nx = param_real_data.image_size_Nx;
        Ny = param_real_data.image_size_Ny;
        N  = Nx * Ny;


        %% load noisy real data
        
        [yfull, uw, vw, sigma_noise, nWw,ant1,ant2,time] = util_load_real_data(visibility_file_name, param_real_data);
        
        
        %% compute weights for Projection into the ellipsoids
        param_precond.N = N; % number of pixels in the image
        param_precond.Nox = ox*Nx; % number of pixels in the image
        param_precond.Noy = oy*Ny; % number of pixels in the image
        aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
        
        %% set the blocks structure
        if param_block_structure.use_manual_partitioning == 1
            if flag_single_data_set
               param_block.pos = length(uw);
            end
            out_block =util_antenna_based_block_sp_ar(uw,ant1,param_block);
%            out_block =util_time_based_block_sp_ar(time,param_block);
            partition = out_block.partition;
        else
            partition = length(uw); % one block
        end
        
        param_block_structure.partition =partition; 
        [u, v, ~, uvidx, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
        
        %% measurement operator and sparsity operator initialization 
        fprintf('Initializing the NUFFT operator\n\n');
        tstart = tic;
        [A, At, T, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
        tend = toc(tstart);

        fprintf('Initialization runtime: %ds\n\n', ceil(tend));
        R = length(v);     
        y = cell(R, 1);
        for q = 1:R
             y{q} = yfull(uvidx{q});
        end
        
        %clear unnecessary vars at this stage
        clear uw vw u v  uvidx ;
        clear ant1 ant2 time ;         
   
        % sparsity operator definition
        [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        
        % compute the spectral norm of the measurement operator
        fprintf('Computing operator norms ...\n');
        fprintf('Natural W ...\n');
        evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
        
        fprintf('Preconditioning ...\n');
        evl_precond = op_norm(@(x) sqrt(cell2mat(aW)) .* (Gw * A(x)), @(x) At(Gw' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
       
        evl_blocks = zeros(R, 1);
        %clear unnecessary vars at this stage
            
   

