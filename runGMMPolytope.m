function [colorcodes, binaryMatrix, confidenceScores, trustedCells] = runGMMPolytope(input, sensitivity, confidence)
% runGMMthresholding performs GMM-based clustering on log-transformed nx9 
% Polytope epitope tag MFI data:
% INPUTS:
%   input            - n x 9 matrix of log10-transformed MFI data
%   sensitivity      - scalar (e.g., 0.66) to binarize expression calls
%   confidence       - scalar (e.g., 0.25) to decide if a cell is trustworthy
%
% OUTPUT:
%   colorcodes       - n x 1 vector encoding marker combinations
%   binaryMatrix     - n x 9 array with call decisions based on sensitivity
%   confidenceScores - n x 1 vector with average confidence per cell
%   trustedCells     - n x 1 logical vector, true for cells with high confidence
%
% EXAMPLE USAGE:
%   % Generate example bimodal data, then log10-transform it
%   data = arrayfun(@(x) [randn(500,1)*100 + 1000; randn(500,1)*100 + 10000], 1:9, 'UniformOutput', false);
%   data = log10(cell2mat(data));  % Concatenate and log-transform
%   sensitivity = 0.5;              % Set sensitivity threshold
%   GMM_codes = runGMMthresholding(data, sensitivity); % Run GMM thresholding

if nargin < 2
        sensitivity = 0.66; % default threshold
end

if nargin < 3
        confidence = 0.25; % default threshold
end

%%%%%%% Validate Input %%%%%%%
if any(isnan(input(:))) || any(isinf(input(:)))
    error('Input data contains NaN or Inf values. Please check your data.');
end

isLog = median(input(:)) < 10 && max(input(:)) < 20;
if ~isLog
    warning('Input data appears to be raw. Applying log10 transformation.');
    input = log10(input);
end

%%%%%%% Thresholding on bimodal MFI distribution %%%%%%%
thresholds = zeros(1,9);
for i = 1:9
    [y, x] = ksdensity(input(:,i));
    p = find(islocalmax(y));
    if numel(p) < 2
        error('Less than two peaks detected—data may not be bimodal.');
    end
    [~, s] = sort(y(p), 'descend');
    r = sort(p(s(1:2)));
    [~, m] = min(y(r(1):r(2)));
    thresholds(i) = x(r(1) + m - 1);
end
input_cc = input >= thresholds;

%%%%%%% GMM Preparation %%%%%%%
combinations = nchoosek(1:9, 3);
GMM_container = cell(1, 9);
covtype = 'full';
sharedcov = true;
clusternum = 8;

%%%%%%% Run GMM on all combinations %%%%%%%
for i = 1:84
    idx = combinations(i,:);
    initial = input_cc(:,idx) .* [1 3 5];
    initial = sum(initial, 2);
    initial(initial > 1) = initial(initial > 1) - 1;
    initial(initial < 6) = initial(initial < 6) + 1;

    %gmModel = fitgmdist(input(:,idx), clusternum, 'CovarianceType', covtype, 'SharedCovariance', sharedcov, 'Start', initial, 'MaxIter', 500);
    % Create a statset object with convergence options
    options = statset('MaxIter', 800, 'TolFun', 1e-8);
    % Now pass this 'Options' object into fitgmdist
    gmModel = fitgmdist(input(:,idx), clusternum, ...
    'CovarianceType', covtype, ...
    'SharedCovariance', sharedcov, ...
    'Start', initial, ...
    'Options', options);
   
    idx_gm = cluster(gmModel, input(:,idx));

    GMM_container{idx(1)} = [GMM_container{idx(1)}, ismember(idx_gm, [2 4 6 8])];
    GMM_container{idx(2)} = [GMM_container{idx(2)}, ismember(idx_gm, [3 7 4 8])];
    GMM_container{idx(3)} = [GMM_container{idx(3)}, ismember(idx_gm, [5 7 6 8])];
end

%%%%%%% Get likelihood matrix %%%%%%%
likelihood = [mean(GMM_container{1,1},2), ...
    mean(GMM_container{1,2},2), mean(GMM_container{1,3},2), ...
    mean(GMM_container{1,4},2), mean(GMM_container{1,5},2), ...
    mean(GMM_container{1,6},2), mean(GMM_container{1,7},2), ...
    mean(GMM_container{1,8},2), mean(GMM_container{1,9},2)];

%%%%%% Binarize based on sensitivity %%%%%%
binaryMatrix = likelihood > sensitivity;

%%%%%% Compute confidence scores %%%%%%
confidenceScores = abs(likelihood - 0.5) * 2;  % 0 = ambiguous, 1 = clear;
trustedCells = all(confidenceScores >= confidence, 2);

%%%%%%% Generate color code list %%%%%%%
str_GMM_colorcodes = string(binaryMatrix .* (1:9));
str_GMM_colorcodes = join(str_GMM_colorcodes, '', 2);
str_GMM_colorcodes = erase(str_GMM_colorcodes, "0");
colorcodes = double(str_GMM_colorcodes);
colorcodes(isnan(colorcodes)) = 0;

fprintf(strcat("GMM thresholding completed successfully. ", num2str(length(unique(colorcodes))), " unique color codes detected"));
end
