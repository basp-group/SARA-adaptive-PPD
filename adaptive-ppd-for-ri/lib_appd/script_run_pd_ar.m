    % wavelet mode is a global variable which does not get transfered
    % to the workers; we need to set it manually for each worker
    dwtmode('per');    
    tstart_a = tic;
    fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond_wave_par\n');
    [result_st.sol, result_st.L1_v, result_st.L1_vp, ...
        result_st.L2_v,result_st.L2_vp, result_st.delta_v,...
        result_st.sol_v,result_st.snr_v, result_st.no_sub_itr_v, ~, ~, ...
        result_st.sol_reweight_v] ...
        = solver_adaptive_ppd_fb(yT,...
        epsilonT, epsilonTs, epsilon, epsilons, ...
        A, At, T, aW, W, Psi, Psit, Psiw, Psitw, ...
        param_pdfb_precond);

    tend = toc(tstart_a);
    fprintf(' pdfb_bpcon_par_sim_rescaled_precond_wave_par runtime: %ds\n\n', ceil(tend));

    result_st.time = tend;
    if ~use_real_visibilities
        error = im - result_st.sol;
        result_st.snr = 20 * log10(norm(im(:))/norm(error(:)));
    end
    result_st.no_itr = length(result_st.L1_v);

    wcoef = [];
    for q = 1:length(Psit)
        wcoef = [wcoef; Psit{q}(result_st.sol)];
    end
    result_st.sparsity = sum(abs(wcoef) > 1e-3)/length(wcoef);
    



