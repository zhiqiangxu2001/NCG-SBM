function [modularity,conductance] = evaluate_quality(adj_mat,partition)
%
% Author - Zhiqiang Xu, 03/2014
%
% Email  - zxu1@e.ntu.edu.sg
% 
% Description - evaluates the structure quality of the resulting clustering;
%               modularity and conductant are used to measure the structure quality.
%
% Input  - adj_mat    : adjacency matrix, size NxN
%        - partition  : the final cluster labels, size Nx1
%
% Output - modularity : scalar
%        - conductance : scalar
% -------------------------------------------------------------------------
      
    cltId = unique(partition);
    N = length(partition);
    K = length(cltId);   

	% ----modularity for both weighted and unweighted graphs------
	
    M = sum(adj_mat(:));    
    
    rho = zeros(K);    
    clt = true(N,K);
    for k = 1:K
        clt(:,k) = logical(partition==cltId(k));
        rho(k,k)   = sum(sum(adj_mat(clt(:,k),clt(:,k))))/M;
    end    
    for i = 1:(K-1)        
        for j = (i+1):K            
            rho(i,j) = sum(sum(adj_mat(clt(:,i),clt(:,j))))/M;
            rho(j,i) = rho(i,j);
        end
    end
    
    alpha = sum(rho,2);
    modularity = sum(diag(rho)-alpha.^2);
    conductance = mean((sum(rho)-diag(rho)')./sum(rho));

end