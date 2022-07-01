function [ fit ] = acmThreeStepRegression( data, factors, param )
% Calculates fitted model parameters via threee step linear regression.

%% Calculate VAR.
fit.phi     = ((factors.m.X(:, 1:end-1)*factors.m.X(:, 1:end-1)')\factors.m.X(:, 1:end-1)*factors.m.X(:, 2:end)')';
fit.mu      = zeros(param.K, 1);
fit.vHat    = factors.m.X(:, 2:end) - fit.phi*factors.m.X(:, 1:end-1) - repmat(fit.mu, 1, param.m.T-1);
temp.s0     = cov(fit.vHat');
fit.s0      = temp.s0(:);
fit.s1      = zeros(param.K^2, param.K);
clear temp

%% Calculate excess returns.
temp.lBondPrices  	= -(data.m.ylds/100).*(repmat(param.allMatur, param.m.T, 1));
data.m.rx           = temp.lBondPrices(2:end, 1:end-1) - temp.lBondPrices(1:end-1, 2:end) + ...
                            repmat(temp.lBondPrices(1:end-1, 1), 1, param.N - 1);
data.m.rx           = data.m.rx(:, param.treasRtns-1)';
clear temp

%% Lambda estimation
% OLS
temp.Y              = data.m.rx';
temp.X              = [ones(param.m.T-1, 1), factors.m.X(:, 1:end-1)' fit.vHat'];
temp.b              = (temp.X'*temp.X)\temp.X'*temp.Y; 
temp.e              = temp.Y - temp.X*temp.b; 

fit.omega           = eye(param.rxN).*var(temp.e(:));
temp.b              = temp.b';
fit.beta           	= temp.b(:, (end-param.K+1):end);
fit.betaStar        = zeros(param.rxN, param.K^2);
for i = 1:param.rxN
    fit.betaStar(i, :)  = kron(fit.beta(i, :), fit.beta(i, :))';
end

% Convexity-adjusted lambda
temp.factors 	= [ones(param.m.T - 1, 1), factors.m.X(:, 1:end-1)'];
temp.factors 	= temp.factors - fit.vHat'*((fit.vHat*fit.vHat')\fit.vHat*temp.factors);
temp.Y        	= ((temp.factors'*temp.factors)\temp.factors'*(data.m.rx + (1/2)*...
                 	(repmat(fit.betaStar*fit.s0 + diag(fit.omega), 1, param.m.T-1)))')';
temp.X        	= fit.beta;
temp.Lambda  	= (temp.X'*temp.X)\temp.X'*temp.Y; 
fit.lambda0   	= temp.Lambda(:,1);
fit.lambda1   	= temp.Lambda(:,2:end);
fit.muStar      = fit.mu - fit.lambda0;
fit.phiStar     = fit.phi - fit.lambda1;
clear temp

%% Calculate short rate parameters.
temp.Y        	= data.m.ylds(:,1)/1200;
temp.X        	= [ones(param.m.T,1), factors.m.X'];
temp.Delta  	= (temp.X'*temp.X)\temp.X'*temp.Y;
fit.delta0      = temp.Delta(1)';
fit.delta1      = temp.Delta(2:(1+param.K)); 
clear temp

end

