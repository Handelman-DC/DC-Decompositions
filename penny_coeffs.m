% Specify the folder where your text files are located
folder = 'penny';

% Initialize an empty matrix to store your data
dataMatrix_penny = zeros(127*127, 35);

% Loop through the files in the folder
for i = 0:126
    for j = 0:126
%         disp(i)
        % Generate the filename for the current file
        filename = fullfile(folder, sprintf('coeff_%d_%d.txt', i, j));

        % Check if the file exists
        if exist(filename, 'file')
            % Read the data from the current file and store it in the matrix
            data = dlmread(filename);
            
            % Reshape the data into a row vector
            dataVector = reshape(data', 1, []);

            % Store the data in the matrix
            dataMatrix_penny((i*127 + j + 1), :) = dataVector;
        else
            fprintf('File %s does not exist.\n', filename);
        end
    end
end