function [outputArg1,outputArg2] = plot_binary(A)
% plot a binary matrix, 1 as black, 0 as white.
h = figure();
tiledlayout(1,1, 'TileSpacing', 'none', 'Padding', 'none');
A_show = ones(size(A));
A_show(:, :, 1) = 1 - A;
A_show(:, :, 2) = 1 - A;
A_show(:, :, 3) = 1 - A;
imagesc(A_show)
% set(gca,'visible','off')
truesize(h,[size(A, 1) * 20,  size(A, 2) * 20]);
title('binary data')

end