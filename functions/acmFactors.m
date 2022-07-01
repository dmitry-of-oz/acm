function [ factors ] = acmFactors( data, param )
% Calculate ACM state variables.
% Calculate principal components at monthly and daily frequencies
temp.pcYlds                     = data.m.ylds(:, param.pcMatur);
temp.meanYlds                   = mean(temp.pcYlds);
temp.pcYldsDemeaned             = temp.pcYlds - repmat(temp.meanYlds, param.m.T, 1);
[temp.svd, ~, ~]                = svd(temp.pcYldsDemeaned'*temp.pcYldsDemeaned);
factors.pcLoadings              = temp.svd(:, 1:param.K);
factors.m.X                     = temp.pcYldsDemeaned*factors.pcLoadings;

% Scale factors.
[factors.m.X, ~, temp.sigmaX]   = zscore(factors.m.X);
factors.pcLoadings              = factors.pcLoadings./repmat(temp.sigmaX, length(param.pcMatur), 1);

% Enforce average positive loadings.
temp.neg                        = mean(factors.pcLoadings) < 0;
factors.m.X(:, temp.neg)        = -factors.m.X(:,temp.neg);
factors.pcLoadings(:, temp.neg) = -factors.pcLoadings(:,temp.neg);
factors.m.X                     = factors.m.X';

%     Calculate daily factors.
factors.d.X                     = (data.d.ylds(:, param.pcMatur) - repmat(temp.meanYlds, param.d.T, 1))*factors.pcLoadings;
factors.d.X                     = factors.d.X';
clear temp

end

