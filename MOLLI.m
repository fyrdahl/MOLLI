% Toy example of MOLLI-reconstruction
% -------------------------------------------
% This is just for educational purposes. It probably contains errors and
% it's most certainly not very time efficient.
%
% Alexander Fyrdahl, Karolinska Instituet, 2018

clear variables; clc;

% Load data; this is MOCO data from the scanner
load('data.mat');

% Crop image speed up calculation
im = im(25:105,90:210,:);
im_size = [size(im,1), size(im,2)];

% Standard 3p-MOLLI
molli3p = @(p,x) p(1)-p(2)*exp(-x./p(3));
init = [1000 1000 1000];

opts = optimset('Display','off');

% For masking negative part of relaxation curve
flip_mask = @(n) ones(numel(inv_time),1)-2*tril(ones(1,numel(inv_time)),n-1)';

fit = zeros(3,prod(im_size));
err = zeros(1,prod(im_size));

im = reshape(im, prod(im_size), size(im,3));

% Loop over all voxels and perform pixel-wise fit
for ii = 1:prod(im_size)

    fprintf('processing voxel %d of %d\n', ii, prod(im_size));

    S = im(ii,:)';
    
    thisFit = zeros(3,4);
    fval = zeros(1,4);

    % The provided data are magntiude only, so I perform multiple fits with the
    % n = [0:4] first values flipped. The lowest residual wins.

    [thisFit(:,1),fval(1)] = fminsearch(@(p) sum((S.*flip_mask(0) - molli3p(p,inv_time)).^2), init, opts);
    [thisFit(:,2),fval(2)] = fminsearch(@(p) sum((S.*flip_mask(1) - molli3p(p,inv_time)).^2), init, opts);
    [thisFit(:,3),fval(3)] = fminsearch(@(p) sum((S.*flip_mask(2) - molli3p(p,inv_time)).^2), init, opts);
    [thisFit(:,4),fval(4)] = fminsearch(@(p) sum((S.*flip_mask(3) - molli3p(p,inv_time)).^2), init, opts);

    [err(ii),idx] = min(fval);
    fit(:,ii) = thisFit(:,idx);
    
end
%% Display the map

A = fit(1,:);
B = fit(2,:);
C = fit(3,:); % T1*

t1 = reshape(C.*(B./A-1),im_size); %% Look-Locker correction
err = reshape(err,im_size);

c = 1300;
w = 1300;
low = c - w / 2;
high = c + w / 2;

figure(99);
imshow(fliplr(permute(t1,[1 2])),[low high]);
colormap(jet); % I'm guessing the MyoMaps colormap is copyrighted, so let's go with something else...
axis image;