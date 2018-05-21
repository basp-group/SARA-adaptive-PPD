
% wavelet mode is a global variable which does not get transfered
    % to the workes; we need to set it manually for each worker
    dwtmode('per');    
    tstart_a = tic;
    fprintf(' Running Adaptive PPD \n');
    [result_st.sol, result_st.L1_v, result_st.L1_vp, ...
        result_st.L2_v,result_st.L2_vp, result_st.delta_v,...
        result_st.sol_v, result_st.no_sub_itr_v, ~, ~, ...
        result_st.sol_reweight_v] ...
        = solver_adaptive_ppd_fb(y,...
        epsilonT, epsilonTs, epsilon, epsilons, ...
        A, At, T, aW, W, Psi, Psit, Psiw, Psitw, ...
        param_pdfb_precond);

    tend = toc(tstart_a);
    fprintf(' Adaptive PPD runtime: %ds\n\n', ceil(tend));

    result_st.time = tend;
   
    result_st.no_itr = length(result_st.L1_v);

    wcoef = [];
    for q = 1:length(Psit)
        wcoef = [wcoef; Psit{q}(result_st.sol)];
    end
    result_st.sparsity = sum(abs(wcoef) > 1e-8)/length(wcoef);
    result_st.residualImage= At(Gw'*(yfull-Gw*(A(result_st.sol))));       
    result_st.ox = ox;
    result_st.oy = oy;
    result_st.Kx = Kx;
