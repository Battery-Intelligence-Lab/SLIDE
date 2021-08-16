% This script is a unit test written to test results.
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.

% THIS FILE IS FOR UNIT TESTS PURPOSE.

clear all; close all; clc;

pathVar;

fileList = dir(fullfile(pathvar.results_folder, '*','*.csv'));

n = length(fileList);

file_diff = [];
fprintf('Results checking is started.\n');
for i= 1:n
    
    [~, y] = fileparts(fileList(i).folder);
    
    origPath = fullfile(pathvar.orig_results_folder, y, fileList(i).name);
    file_diff_now = 0;
    if(exist(origPath,'file')==2)
        filePath = fullfile(fileList(i).folder, fileList(i).name);
        try
        testFile = csvread(filePath);
        origFile = csvread(origPath);
        
        file_diff_now = norm(origFile-testFile,1);
        file_diff = [file_diff, file_diff_now];
        catch
            continue;
        end
        

        
    else
        fprintf('%s\n', [origPath, ' is skipped.']);
    end
    
    
    if (mod(i,25)==0 || i==n)
        fprintf('Test of %d/%d files has finished.\n',i,n);
        
    end
    
    if (file_diff_now>eps)
        fprintf('File difference: %4.4f/%4.4f.\n',file_diff_now, norm(origFile,1));
        disp(origPath);
        disp(i);
        [row,col,~] = ind2sub(size(origFile), find(abs(origFile - testFile)>eps));
        row'
        col'
        %%pause();
    end
    
end


if(all(file_diff<eps))
    fprintf('Test has passed!\n');
else
    fprintf('Test has failed!\n');
end





