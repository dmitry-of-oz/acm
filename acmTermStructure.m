%%  ACM Term Structure Model
clear
addpath('functions')

%% Initialize paramters.
param.pcMatur             	= 3:120;
param.N                     = param.pcMatur(end);
param.treasRtns           	= [6, 12:12:120];
param.rxN                   = length(param.treasRtns);
param.allMatur            	= (1:120)/12';
param.N                     = length(param.allMatur);
param.K                     = 5;

%% Import GSW parameters and extract zero-coupon yields from GSW parameters.
[~, ~, temp.raw, temp.dateNums]	= xlsread('data\gswParams.xls','gswParams','','',@convertSpreadsheetDates);
temp.raw                        = temp.raw(2:end,:);
temp.dateNums                   = temp.dateNums(2:end,:);
temp.R                          = ~cellfun(@isequalwithequalnans,temp.dateNums,temp.raw) & cellfun('isclass',temp.raw,'char'); % Find Excel dates
temp.raw(temp.R)              	= temp.dateNums(temp.R);
temp.R                          = cellfun(@(x) ~isnumeric(x) || isnan(x),temp.raw); % Find non-numeric cells
temp.raw(temp.R)             	= {0.0}; % Replace non-numeric cells
temp.data                       = cell2mat(temp.raw);

temp.data(isnan(temp.data(:,2)), :)   = [];
temp.data(isnan(temp.data))           = 0;

temp.DATE                       = temp.data(:,1);
temp.BETA0                      = temp.data(:,2);
temp.BETA1                      = temp.data(:,3);
temp.BETA2                      = temp.data(:,4);
temp.BETA3                      = temp.data(:,5);
temp.TAU1                       = temp.data(:,6);
temp.TAU2                       = temp.data(:,7);

clearvars temp.data temp.raw temp.dateNums temp.R temp.columnIndices;

param.d.T       = size(temp.DATE,1);
data.d.ylds     = nan(param.d.T, param.N);
for t=1:param.d.T
        temp.nTau1          = ((1:param.N)./12)./temp.TAU1(t);
        temp.nTau2          = ((1:param.N)./12)./temp.TAU2(t);
        data.d.ylds(t,:)	= temp.BETA0(t) + ...
                                temp.BETA1(t).*((1-exp(-temp.nTau1))./temp.nTau1) + ...
                            	temp.BETA2(t).*((1-exp(-temp.nTau1))./temp.nTau1 - exp(-temp.nTau1)) + ...
                            	temp.BETA3(t).*((1-exp(-temp.nTau2))./temp.nTau2 - exp(-temp.nTau2));
end
data.d.dates                = datenum(temp.DATE);

temp.monthIdx            	= find(month2(data.d.dates(1:(end-1))) ~= month2(data.d.dates(2:end)));
data.m.dates                = data.d.dates(temp.monthIdx);
data.m.ylds                 = data.d.ylds(temp.monthIdx, :);
param.m.T                   = length(data.m.dates);
clear temp

%% Smooth the 1m GSW (risk free rate) for dates before Jan 1982
[~, ~, temp.raw, temp.dateNums]	= xlsread('data\ffr.xls','ffr','','',@convertSpreadsheetDates);
temp.raw                        = temp.raw(2:end,:);
temp.dateNums                   = temp.dateNums(2:end,:);
temp.R                          = ~cellfun(@isequalwithequalnans,temp.dateNums,temp.raw) & cellfun('isclass',temp.raw,'char');
temp.raw(temp.R)               	= temp.dateNums(temp.R);
temp.data                       = cell2mat(temp.raw);
temp.ffrDate                    = temp.data(:,1);
temp.effectiveFFR               = temp.data(:,2);
clearvars raw dateNums R columnIndices;

temp.idx                  	= find(datenum('29-Jan-1982') == data.m.dates);
temp.Y                   	= data.m.ylds(temp.idx:end-length(data.m.ylds(temp.idx+length(temp.effectiveFFR(temp.idx:end, 1)):end, 1)), 1);
temp.X                      = [ones(length(temp.Y),1), temp.effectiveFFR(temp.idx:end, 1)];
temp.beta                   = temp.X\temp.Y;
temp.predict              	= [ones(length(temp.effectiveFFR(:,1)),1), temp.effectiveFFR(:,1)]*temp.beta;
data.m.ylds(1:temp.idx-1,1)	= temp.predict(1:temp.idx-1,1);
data.m.effectiveFFR       	= temp.effectiveFFR(:,1);
clear temp

%% Calculate state variables.
factors         = acmFactors(data, param);

%% Calculate fitted model parameters via threee step linear regression.
fit             = acmThreeStepRegression(data, factors, param);

%% Calculate recursive pricing factors.
fit             = acmRecursion(fit, param);

%% Calculate fitted yields, risk neutral yields, and term premia
% Fitted yields
temp.multiplier      	= repmat(param.allMatur', 1, param.m.T)/100;
fit.m.ylds              = -((repmat(fit.A, 1, param.m.T) + fit.B*factors.m.X)./temp.multiplier)';
fit.m.yldsRn            = -((repmat(fit.Arf, 1, param.m.T) + fit.Brf*factors.m.X)./temp.multiplier)';
fit.m.termPremium       = fit.m.ylds - fit.m.yldsRn;

% Daily estimates
temp.multiplier         = repmat(param.allMatur'/100, 1, param.d.T);
fit.d.ylds              = -((repmat(fit.A, 1, param.d.T) + fit.B*factors.d.X)./temp.multiplier)';
fit.d.yldsRn            = -((repmat(fit.Arf, 1 ,param.d.T) + fit.Brf*factors.d.X)./temp.multiplier)';
fit.d.termPremium       = fit.d.ylds - fit.d.yldsRn;  

clear temp

%% Plot fit and decomposition
figure('Position', [0, 100, 1000, 400])
subplot(1,2,1)
plot(data.m.dates, data.m.ylds(:,120), 'LineWidth', 1); hold on;
plot(data.m.dates, fit.m.ylds(:,120), 'r--', 'LineWidth', 1); hold off;
datetick('x', 'yyyy')
xlim([data.m.dates(1) data.m.dates(end)]);
legend('Actual Yield', 'Fitted Yield')
title('10-Year Model Fit')

subplot(1,2,2)
plot(data.m.dates, fit.m.ylds(:,120), 'b', 'LineWidth', 1); hold on;
plot(data.m.dates, fit.m.yldsRn(:,120), 'g', 'LineWidth', 1); hold on;
plot(data.m.dates, fit.m.termPremium(:,120), 'r', 'LineWidth', 1); hold off;
datetick('x', 'yyyy')
xlim([data.m.dates(1) data.m.dates(end)]);
legend('Yield', 'Risk Neutral Yield', 'Term Premium')
title('10-Year Yield Decomposition')

%% Output decomposition
temp.outputMaturs   = 12:12:120;
temp.excelDates     = data.d.dates - datenum('30-Dec-1899');
temp.header         = cell(1, length(temp.outputMaturs)+1);
temp.header{1}      = 'date';

for i=1:length(temp.outputMaturs)
    temp.header{i+1}     = [num2str(temp.outputMaturs(i)/12) 'y yield'];
end
temp.output         = [temp.excelDates, fit.d.ylds(:,temp.outputMaturs)];
xlswrite('output\acmDecomposition.xlsx', [temp.header; num2cell(temp.output)], 'ACM Fitted Yields')

for i=1:length(temp.outputMaturs)
    temp.header{i+1}     = [num2str(temp.outputMaturs(i)/12) 'y term premium'];
end
temp.output         = [temp.excelDates, fit.d.termPremium(:,temp.outputMaturs)];
xlswrite('output\acmDecomposition.xlsx', [temp.header; num2cell(temp.output)], 'ACM Term Premium')

for i=1:length(temp.outputMaturs)
    temp.header{i+1}     = [num2str(temp.outputMaturs(i)/12) 'y risk neutral yield'];
end
temp.output         = [temp.excelDates, fit.d.yldsRn(:,temp.outputMaturs)];
xlswrite('output\acmDecomposition.xlsx', [temp.header; num2cell(temp.output)], 'ACM Risk Neutral Yields')
