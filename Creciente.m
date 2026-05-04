% Parameters
Nr = 14; % Number of receiving antennas
Nt = 2;  % Number of transmitting antennas

% Calculate alpha based on the number of receiving antennas
alpha = Nr * (2 * Nr - 1); % Total number of rows in B0
disp(['Value of alpha: ', num2str(alpha)]);

% Total positions
total_positions = 2 * Nr; 
% Initialize B0
B0 = zeros(alpha, total_positions); % Preallocate B0 with zeros

% Set to track unique combinations
unique_combinations = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Generate combinations with +1 and -1
row_count = 0; % Initialize row counter

while row_count < alpha
    % Randomly select positions for +1 and -1
    pos_plus = randi(total_positions); % Random position for +1
    pos_minus = randi(total_positions); % Random position for -1
    
    %
    if pos_plus ~= pos_minus 
        row = zeros(1, total_positions); % Create a row of zeros
        row(pos_plus) = 1;  % Assign +1
        row(pos_minus) = -1; % Assign -1
        
        % 
        if pos_plus > pos_minus
            row(pos_plus) = -1;
            row(pos_minus) = 1; % Swap positions if needed
        end
        
        % Create a unique key for the current row
        row_key = num2str(row);
        
        % Check if the combination is unique
        if ~isKey(unique_combinations, row_key)
            row_count = row_count + 1; % Increment row counter
            
            % Store the generated row in B0
            B0(row_count, :) = row; % Add the row to B0
            unique_combinations(row_key) = 1; % Mark this combination as used
        end
    end
end

% Display the resulting B0 matrix
disp('Matrix B0:');
disp(B0);

% Check dimensions of B0
disp(['Dimensions of B0: ', num2str(size(B0))]);

% Calculate the number of comparators based on the number of rows in B0
numComparators = size(B0, 1); % Total number of comparators is equal to the number of rows in B0
disp(['Total number of comparators based on B0: ', num2str(numComparators)]);

% Calculate the number of comparators based on alpha
numComparators_alpha = alpha; % Total number of comparators is equal to alpha
disp(['Total number of comparators based on alpha: ', num2str(numComparators_alpha)]);

% Plot the total number of comparators as a function of the number of receiving antennas
Nr_values = 2:20; % Number of receiving antennas (from 2 to 20)
numComparatorsArray = zeros(size(Nr_values)); % Initialize the vector to store the number of comparators

% Calculate the number of comparators for each number of receiving antennas
for idx = 1:length(Nr_values)
    n = Nr_values(idx);
    numComparatorsArray(idx) = nchoosek(2*n, 2); % αf = C(2, 2Nr)
end

% Plot the total number of comparators
figure;
plot(Nr_values, numComparatorsArray, '-o', 'MarkerFaceColor', 'b');
grid on;
xlabel('Number of Receiving Antennas (Nr)');
ylabel('Total Number of Comparators');
title('Number of Comparators in a Fully Connected Network');
xlim([2 20]);
ylim([0 max(numComparatorsArray) + 10]);
