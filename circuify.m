function [rcv, rch, rrv, rrh] = circuify(img)
  mindetail = 0.0005;
  surfacewidth = 0.05;
  N = surfacewidth / mindetail;
  x = imgcondition(imgread(img), N);
  kc = normalizekernel([1, 0, 1, 1, 1, 0, 1; ...
                        0, 0, 1, 1, 1, 0, 0; ...
                        1, 0, 1, 1, 1, 0, 1]);
  kr = normalizekernel([1, 0, 1]);
  rcv = xcorr2(x, kc);
  rch = xcorr2(x, kc.');
  rrv = xcorr2(x, kr);
  rrh = xcorr2(x, kr.');
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
  y = y / max(y(:));
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
