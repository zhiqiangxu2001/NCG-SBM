function Rslt = ncg_sbm(tOmega,Mat,Option)
% approximate Riemmannian nonlinear conjugate gradient in canonical
% parameters (i.e., natural parameters)

    tic;

    Elbo  = zeros(Option.maxiter,1);
    zcurve = zeros(Option.N,Option.maxiter);    
    tcurve = zeros(Option.maxiter,1);     
    [tOmegaK,tTheta,expsum] = transform(tOmega);
   
    i  = 1;
    it = 1;
    steplength = 1;
    while true
        s0 = sum(tTheta);
        s1 = tTheta'*Mat*tTheta;
        s2 = s0'*s0-tTheta'*tTheta;
        
        tAlpha = Option.alpha + s0;
        tBeta  = [Option.beta(1)+0.5*diag(s1),Option.beta(2)+0.5*(diag(s2)-diag(s1))];
        
        Elbo(i) = lower_bound(tAlpha,tBeta,tTheta,s1,s2,Option);
        
        if i>1 
            if Elbo(i) < Elbo(i-1) 
                % backtracking.
                if it == Option.maxiter
                    disp('stop: maxiter exceeded');
                    [~,tTheta] = transform(tOmega_old);
                    i = i-1;
                    break;
                end
                steplength = 0.5*steplength;
                tOmega = tOmega_old + abs((Elbo(i)-Elbo(i-1))/Elbo(i))*steplength*d;
                [tOmegaK,tTheta,expsum] = transform(tOmega);
                it = it+1;
                continue;
            end
            if abs((Elbo(i)-Elbo(i-1))/Elbo(i)) <= Option.ftol
                disp('stop: NCG converged (ftol)');
                break;
            end
        end        
        if it == Option.maxiter
            disp('stop: maxiter exceeded');
            break;
        end
        [~,zcurve(:,i)] = max(tTheta,[],2);
        tcurve(i) = toc;        
        
        c1 = psi(tAlpha);
        c2 = psi([tBeta sum(tBeta,2)])';
        v1 = Mat*tTheta;
        v2 = bsxfun(@minus,s0,tTheta);
        g = bsxfun(@minus,Option.coef(2)*(Option.N-1)-1,tOmegaK)+bsxfun(@plus,log(expsum),c1);
        g = g + Option.coef(1)*bsxfun(@minus,Option.en',v1)-Option.coef(2)*v2;
        g = g + bsxfun(@times,c2(1,:),v1);
        g = g + bsxfun(@times,c2(2,:),v2-v1);
        g = g - bsxfun(@times,c2(3,:),v2);
        
        tG = bsxfun(@minus,g(:,1:Option.K-1),g(:,end));                       % natural gradient in canonical (natural) parameter is equal to the gradient in expectationa parameter.
        tGnorm2 = sum(tTheta(:).*(g(:).^2))-sum((tTheta(:).*g(:)).^2);
        if tGnorm2 <= Option.gtol
            disp('stop: NCG converged (gtol)');
            break;
        end
        if it>1
            % Fletcher-Reeves
            beta = tGnorm2/tGnorm2_old;
            if ~islegal(beta)
                beta = 0;
            end
        end
        tGnorm2_old = tGnorm2;        
              
        if it>1 && beta>0
            d = tG + beta*d;    % CG step
        else            
            d = tG;             % steepest ascend step
        end
        tOmega_old = tOmega;
        tOmega = tOmega + steplength*d;    % the step length is fixed to 1 before decreasing.     
        [tOmegaK,tTheta,expsum] = transform(tOmega);
               
        i = i+1;
        it = it+1;
    end  
    [~,zcurve(:,i)] = max(tTheta,[],2);
    tcurve(i) = toc;    

    Rslt.Option = Option;
    Rslt.f = Elbo(1:i);
    Rslt.zcurve = zcurve(:,1:i);
    Rslt.tcurve = tcurve(1:i);
    Rslt.t = toc;
end