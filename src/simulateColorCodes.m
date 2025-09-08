function [colorcodes, summary] = simulateColorCodes(n, recombinationRates)
% SIMULATECOLORCODES Simulates barcode recombination outcomes.
% 
%   [colorcodes, summary] = simulateColorCodes(n, recombinationRates)
%   Inputs:
%     n - number of cells (e.g., 1000)
%     recombinationRates - 1x10 vector with probability for each recomb.
%     step count (e.g., [0.05 0.15 0.225 0.225 0.17 0.10 0.05 0.02 0.01 0])
%   Outputs:
%     colorcodes - a column vector of resulting color codes (integers)
%     summary - a matrix where column 1 is unique color codes, column 2 is counts

    if length(recombinationRates) ~= 10
        error('recombinationRates must be a 1x10 vector.');
    end

    % Normalize the recombination rates to sum to 1
    recombinationRates = recombinationRates / sum(recombinationRates);

    % Generate initial barcode matrix
    BC = ones(n,1) * 123456789;  % Each cell starts with the same barcode

    % Simulate 9 rounds of potential recombination
    for rounds = 1:9
        for i = 1:n
            bc = num2str(BC(i,rounds));
            if bc == "0" || length(bc) == 1
                BC(i, rounds+1) = 0;
            else
                sz = length(bc);
                loxp = sort(randi(sz, 1, 2));
                bc = eraseBetween(bc, loxp(1), loxp(2));
                if isempty(bc)
                    BC(i, rounds+1) = 0;
                else
                    BC(i, rounds+1) = str2double(bc);
                end
            end
        end
    end

    % Sample recombination step counts for each cell
    stepCounts = randsample(0:9, n, true, recombinationRates);

    % Extract the corresponding barcodes after each cell's final recomb step
    colorcodes = zeros(n, 1);
    for i = 1:n
        colorcodes(i) = BC(i, stepCounts(i) + 1);  % +1 because round 0 is col 1
    end

    % Compute summary of unique codes and their counts
    uniqueCodes = unique(colorcodes);
    counts = histc(colorcodes, uniqueCodes);
    summary = [uniqueCodes, counts];

    % Sort by counts in descending order
    summary = sortrows(summary, -2);
end
