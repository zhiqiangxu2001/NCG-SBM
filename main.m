function [likelihood_bound,modularity,conductance,time_cost] = main(dataset,cluster_num,initdata)
% dataset: adjacency matrix named Mat
% cluster_num: cluster number
% initdata: initial value for natural parameter tOmega

    load(dataset);  % introduce Mat variable.
    Option = set_opt(Mat,cluster_num);
    if isempty(initdata)
        tOmega = init(Option.N,Option.K-1);
    else
        load(initdata); % introduce tOmega variable
    end
    Rslt = ncg_sbm(tOmega,Mat,Option);    
    likelihood_bound = Rslt.f;
    modularity = zeros(size(likelihood_bound));
    conductance = zeros(size(likelihood_bound));
    time_cost = zeros(size(likelihood_bound));
    for i=1:length(modularity)
        [modularity(i),conductance(i)] = evaluate_quality(Mat,Rslt.zcurve(:,i));
        time_cost(i)  = Rslt.tcurve(i)-Rslt.tcurve(1)+1;
    end

end

function Option = set_opt(X,K)

    Option.N = size(X,1);
    Option.K = K;
    
    Option.alpha = ones(1,K);                                                  % hyperparameter of Multinomial(Z_i|alpha)
    Option.beta = [1,1];                                                       % hyperparameter of Beta(phi_kk|beta)

    Option.gtol = 1e-6;
    Option.ftol = 1e-6;    
    Option.maxiter = 200;    
    
    Option.en = sum(X);
    Option.triu = triu(true(K),1);
    
    p = 1e-10;
    Option.coef = [log(p/(1-p)),log(1-p)];

end