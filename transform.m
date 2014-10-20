function [tOmegaK,tTheta,expsum] = transform(tOmega)

    tOmegaK = [tOmega,zeros(size(tOmega,1),1)];
    tOmegaK = bsxfun(@minus,tOmegaK,max(tOmegaK,[],2));
    tTheta = exp(tOmegaK);
    if nargout > 2
        expsum = sum(tTheta,2);
        tTheta = bsxfun(@rdivide,tTheta,expsum);
    else
        tTheta = bsxfun(@rdivide,tTheta,sum(tTheta,2));
    end

end