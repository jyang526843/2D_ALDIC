% ==============================================
% function ReadImage
% ==============================================
function [gridxy] = funReadImageRefUpdate(filename)

% fprintf('Choose method to load images:  \n')
% fprintf('     0: Select images folder;  \n')
% fprintf('     1: Use prefix of image names;  \n')
% fprintf('     2: Manually select images.  \n')
% prompt = 'Input here: ';
% LoadImgMethod = input(prompt);
 
% ==============================================
% Choose ZOI
fprintf('***** Update incremental reference image ***** \n');
disp('--- Choose ZOI two boundary points from the left-top to the right-bottom ---')
imshow( (imread(filename)));

gridx = zeros(1,2); gridy = zeros(1,2);
[gridx(1), gridy(1)] = ginput(1);
fprintf('The left-top coordinates are (%4.3f,%4.3f)\n',gridx(1), gridy(1))

[gridx(2), gridy(2)] = ginput(1);
fprintf('The right-bottom coordinates are (%4.3f,%4.3f)\n',gridx(2), gridy(2))

gridxy.gridx = round(gridx); gridxy.gridy = round(gridy);

% % Choose subset size
% prompt = '--- What is the subset size? --- Input here: ';
% winsize = input(prompt);
% 
% % Choose subset size
% prompt = '--- What is the subset step? --- Input here: ';
% winstepsize = input(prompt);
 


end