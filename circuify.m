function [x, rcv, rch, rrv, rrh] = circuify(img)
  mindetail = 0.001;
  surfacewidth = 0.07;
  N = surfacewidth / mindetail;

  x = imgcondition(imgread(img), N);

  kc = normalizekernel([1, 0, 1, 1, 1, 0, 1; ...
                        0, 0, 1, 1, 1, 0, 0; ...
                        1, 0, 1, 1, 1, 0, 1]);
  rcv = crosscorr(x, kc);
  rch = crosscorr(x, kc.');


  kr = normalizekernel([1, 0, 1]);
  rrv = crosscorr(x, kr);
  rrh = crosscorr(x, kr.');
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
