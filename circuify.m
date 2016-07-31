function [x, r42v, r42h, r43v, r43h, r44v, r44h, r45v, r45h] = circuify(img)
  mindetail = 0.001;
  surfacewidth = 0.07;
  N = surfacewidth / mindetail;

  x = imgcondition(imgread(img), N);

  k42 = foursidedkernel(2);
  r42v = crosscorr(x, k42);
  r42h = crosscorr(x, k42.');
  k43 = foursidedkernel(3);
  r43v = crosscorr(x, k43);
  r43h = crosscorr(x, k43.');
  k44 = foursidedkernel(4);
  r44v = crosscorr(x, k44);
  r44h = crosscorr(x, k44.');
  k45 = foursidedkernel(5);
  r45v = crosscorr(x, k45);
  r45h = crosscorr(x, k45.');
end

function [x] = imgread(img)
  [x, map] = imread(img);
  if ~isempty(map)
    x = ind2rgb(x, map);
  end
end

function [y] = imgcondition(x, N)
  y = rgb2gray(x);
  y = imcomplement(y);
  y = mat2gray(y);
  y = crop2square(y);
  y = imresize(y, [N, N]);
end

function [y] = crop2square(x)
  [S1, S2] = size(x);
  S = min(S1, S2);
  y = x(((S1 - S) / 2 + 1):((S1 + S) / 2), ((S2 - S) / 2 + 1):((S2 + S) / 2));
end

function [k] = foursidedkernel(M)
  k = zeros(2 * M + 3, 2 * M + 3);
  k(3:2:(end - 2), [1, end]) = 1;
  k([1, end], 3:2:(end - 2)) = 1;
  k(3:(end - 2), [3, (end - 2)]) = 1;
  k([3, (end - 2)], 3:(end - 2)) = 1;
  k = normalizekernel(k);
end

function [k] = twosidedkernel(M)
  k = zeros(2 * M - 1, 7);
  k(1:2:end, [1, 7]) = 1;
  k(:, [3, 5]) = 1;
  k([1, end], 4) = 1;
  k = normalizekernel(k);
end

function [k] = resistencekernel()
  k = normalizekernel([1, 0, 1]);
end

function [y] = normalizekernel(x)
  m = mean(x(:));
  y = x - m;
  y = y / sum(sum(x .* y));
end

function [r] = crosscorr(x, k)
  [N, M] = size(k);
  r = xcorr2(x, k);
  r = r(N:end, M:end);
  r((end - N + 2):end, :) = 0;
  r(:, (end - M + 2):end) = 0;
end
