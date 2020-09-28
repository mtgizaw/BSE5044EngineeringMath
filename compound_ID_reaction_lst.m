% compound_ID_reaction_lst.m

% Written by Martha Gizaw
% 9/28/2020

% Input:
% reaction.lst

% Output:
% reactionCompoundIDs.txt

% Description: This MATLAB script counts the number of times a compound ID 
% between C00001 and C00100 appears in the reaction.lst database. That 
% result will copy to the file reactionCompoundIDs.txt.

clear
clc

% Locate and open the file 'reaction.lst.'
[file1, path1] = uigetfile('*.lst', 'Select the reaction.lst database.');
fid1 = fopen(fullfile(path1, file1), 'r');

% Create the output file and name it 'reactionCompoundIDs.txt.'
file2 = 'reactionCompoundIDs.txt';
fid2 = fopen(fullfile(path1,file2),'w');

% Initialize the array compoundID.
compoundID = zeros(0, 100);

% Format the entire reaction.lst database as a string.
lst2Str = fscanf(fid1, '%s');

for n = 1:100
    % Determine the last three digits of a compound ID based on the value 
    % of n.
    if n < 10
        lastThree = strcat('00', num2str(n));
    elseif n >= 10 && n < 100
        lastThree = strcat('0', num2str(n));
    elseif n == 100
        lastThree = num2str(n);
    end

    % Calculate the frequency of each compound ID from C00001 to C00100
    % throughout the entire reaction.lst file.
    compoundID(n) = size(strfind(lst2Str, strcat('C00', lastThree)), 2);
end

% Calculate the sum of all values in the compoundID array. This sum
% represents the frequency of the compound ID between C00001 and C00100
% appearing in reaction.lst.
totalFreq = sum(compoundID);

% Print the total number of compound ID apprearances to the command window.
fprintf(['The compound IDs from C00001 to C00100 appear in the reaction.lst database ' num2str(totalFreq) ' times. \n']);

% Print the total number of compound ID apprearances to the output text
% file.
fprintf(fid2, ['The compound IDs from C00001 to C00100 appear in the reaction.lst database ' num2str(totalFreq) ' times. \n']);

% Close both the input and output files.
fid1 = fclose(fid1);
fid2 = fclose(fid2);

% Print "Done!" to the command window.
fprintf('Done!');