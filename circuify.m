function [x, y] = circuify(img)
  mindetail = 0.001;
  surfacewidth = 0.07;
  N = round(surfacewidth / mindetail);
  x = imgread(img, N);
  y = zeros(N, N);
  y = addfoursided(x, y);
end

function [x] = imgread(img, N)
  [x, map] = imread(img);
  if ~isempty(map)
    x = ind2rgb(x, map);
  end
  x = imgcondition(x, N);
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

function [y] = addfoursided(x, y)
  k = foursidedkernels(length(x));
  r = crosscorrs(x, k);
  [rmax, m, n1, n2] = maxcorr(r);
  rmin = rmax / 2;
  while rmax > rmin
    y = addkernel(y, k{m}, n1, n2);
    r = updatecorrs(r, k, m, n1, n2);
    [rmax, m, n1, n2] = maxcorr(r);
  end
end

function [k] = foursidedkernels(N)
  k = {};
  m = 2;
  ktmp = foursidedkernel(m);
  while length(ktmp) < N / 5
    k{m - 1} = ktmp;
    m = m + 1;
    ktmp = foursidedkernel(m);
  end
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

function [r] = crosscorrs(x, k)
  M = length(k);
  r = cell(M, 1);
  for n = 1:M
    r{n} = crosscorr(x, k{n});
  end
end

function [r] = crosscorr(x, k)
  [N, M] = size(k);
  r = xcorr2(x, k);
  r = r(N:end, M:end);
  r((end - N + 2):end, :) = 0;
  r(:, (end - M + 2):end) = 0;
end

function [rmax, m, n1, n2] = maxcorr(r)
  M = length(r);
  n1 = zeros(M, 1);
  n2 = zeros(M, 1);
  rmax = zeros(M, 1);
  for m = 1:M
    [rtmp, n1tmp] = max(r{m});
    [rmax(m), n2(m)] = max(rtmp);
    n1(m) = n1tmp(n2(m));
  end
  [rmax, m] = max(rmax);
  n1 = n1(m);
  n2 = n2(m);
end

function [y] = addkernel(y, k, n1, n2)
  [N1, N2] = size(k);
  y(n1:(n1 + N1 - 1), n2:(n2 + N2 - 1)) = k > 0;
end

function [r] = updatecorrs(r, k, mk, n1, n2)
  for mr = 1:length(r)
    [N1r, N2r] = size(k{mr});
    [N1k, N2k] = size(k{mk});
    r{mr}(max(0, n1 - N1r):(n1 + N1k), max(0, n2 - N2r):(n2 + N2k)) = 0;
  end
end
