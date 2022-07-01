function fit = acmRecursion( fit, param )
% Calculate ACM recursive pricing factors.
fit.A     	= zeros(param.N,1);
fit.B     	= zeros(param.N, param.K);
fit.A(1,1)	= -fit.delta0;
fit.B(1,:)	= -fit.delta1';    
fit.Arf   	= fit.A;
fit.Brf     = fit.B;

temp.s0term = zeros(param.N,1);
temp.s1term = zeros(param.N,param.K);

% Main recursion
for n = 2:param.N
	temp.Bpb            = kron(fit.B(n-1,:), fit.B(n-1,:));
	temp.s0term(n)      = (1/2)*(temp.Bpb*fit.s0 + fit.omega(1,1));
    temp.s1term(n,:)  	= (1/2)*(temp.Bpb*fit.s1);
	fit.A(n)            = fit.A(n-1) + fit.B(n-1,:)*(fit.mu - fit.lambda0) + temp.s0term(n) + fit.A(1);    
	fit.B(n,:)          = fit.B(n-1,:)*(fit.phi - fit.lambda1) + temp.s1term(n,:) + fit.B(1,:);

	% Risk neutral
	temp.Brfpb          = kron(fit.Brf(n-1,:), fit.Brf(n-1,:));
	fit.Arf(n)          = fit.Arf(n-1) + fit.Brf(n-1,:)*fit.mu + (1/2)*(temp.Brfpb*fit.s0 + fit.omega(1,1)) + fit.Arf(1);    
	fit.Brf(n,:)        = fit.Brf(n-1,:)*fit.phi + (1/2)*temp.Brfpb*fit.s1 + fit.Brf(1,:);
end
clear temp


end

