function Elbo = lower_bound(tAlpha,tBeta,tTheta,s1,s2,Option)

    Elbo = sum(s1(Option.triu))*Option.coef(1)+sum(s2(Option.triu))*Option.coef(2);
    Elbo = Elbo + sum(gammaln(tAlpha))-gammaln(sum(tAlpha)) - (sum(gammaln(Option.alpha))-gammaln(sum(Option.alpha)));
    Elbo = Elbo + sum(betaln(tBeta(:,1),tBeta(:,2))) - Option.K*betaln(Option.beta(1),Option.beta(2));
    Elbo = Elbo - sum(tTheta(logical(tTheta>0)).*log(tTheta(logical(tTheta>0))));

end