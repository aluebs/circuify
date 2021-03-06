function [cnxn, cntx, orly] = circuify(img)
  mindetail = 0.001;
  surfacewidth = 0.07;
  N = round(surfacewidth / mindetail);
  x = imgread(img, N);
  [cntx, orly, mask] = addallkernels(x);
  cnxn = addconnections(x, cntx, mask);
  cntx = cntx > 0;
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

function [cntx, orly, mask] = addallkernels(x)
  [N1, N2] = size(x);
  cntx = zeros(N1, N2);
  orly = zeros(N2, N1);
  mask = ones(N1, N2);
  Nc = 0;
  [k4, c4, o4] = foursidedkernels(min(N1, N2));
  [cntx, orly, mask, Nc] = addkernels(x, cntx, orly, mask, k4, c4, o4, Nc);
  [cntx, orly, mask, Nc] = addtwosided(x, cntx, orly, mask, Nc);
  [kr, cr, or] = resistencekernels();
  [cntx, orly, mask, Nc] = addkernels(x, cntx, orly, mask, kr, cr, or, Nc);
end

function [cntx, orly, mask, Nc] = addkernels(x, cntx, orly, mask, k, c, o, Nc)
  r = crosscorrs(x, k);
  r = updatecorrs(r, k, mask);
  [rmax, m, n1, n2] = maxcorr(r);
  rmin = rmax / 2;
  while rmax > rmin
    [cntx, orly, mask, Nc] = ...
        addkernel(cntx, orly, c{m}, o{m}, n1, n2, mask, Nc);
    r = updatecorrs(r, k, mask);
    [rmax, m, n1, n2] = maxcorr(r);
  end
end

function [k, c, o] = foursidedkernels(N)
  k = {};
  c = {};
  o = {};
  m = 2;
  [ktmp, ctmp, otmp] = foursidedkernel(m);
  while length(ktmp) < N / 5
    k{m - 1} = ktmp;
    c{m - 1} = ctmp;
    o{m - 1} = otmp;
    m = m + 1;
    [ktmp, ctmp, otmp] = foursidedkernel(m);
  end
end

function [k, c, o] = foursidedkernel(M)
  c = zeros(2 * M + 3, 2 * M + 3);
  c(3:2:(end - 2), [1, end]) = 1;
  c([1, end], 3:2:(end - 2)) = 1;
  o = zeros(2 * M + 3, 2 * M + 3);
  o(3:(end - 2), [3, (end - 2)]) = 1;
  o([3, (end - 2)], 3:(end - 2)) = 1;
  k = normalizekernel(c | o);
end

function [cntx, orly, mask, Nc] = addtwosided(x, cntx, orly, mask, Nc)
  [k, c, o] = twosidedkernels();
  r = crosscorrs(x, k);
  r = updatecorrs(r, k, mask);
  [rmax, m, n1, n2] = maxcorr(r);
  rmin = rmax / 2;
  while rmax > rmin
    [cntx, orly, mask, Nc] = ...
        addtwosidedkernel(cntx, orly, r{m}, m, n1, n2, mask, Nc);
    r = updatecorrs(r, k, mask);
    [rmax, m, n1, n2] = maxcorr(r);
  end
end

function [k, c, o] = twosidedkernels()
  k = cell(2, 1);
  c = cell(2, 1);
  o = cell(2, 1);
  [k{1}, c{1}, o{1}] = twosidedkernel(2);
  k{2} = k{1}.';
  c{2} = c{1}.';
  o{2} = o{1}.';
end

function [k, c, o] = twosidedkernel(M)
  c = zeros(2 * M - 1, 7);
  o = zeros(2 * M - 1, 7);
  c(1:2:end, [1, 7]) = 1;
  o(:, [3, 5]) = 1;
  o([1, end], 4) = 1;
  k = normalizekernel(c | o);
end

function [cntx, orly, mask, Nc] = addtwosidedkernel(cntx, orly, r, m, n1, n2, mask, Nc)
  [N1, N2] = size(r);
  switch m
    case 1
      s1 = 2;
      s2 = 0;
    case 2
      s1 = 0;
      s2 = 2;
  end
  if ((n1 + s1) > N1) || ...
     ((n2 + s2) > N2) || ...
     ((n1 > s1) && (n2 > s2) && (r(n1 - s1, n2 - s2) > r(n1 + s1, n2 + s2)))
    s1 = -s1;
    s2 = -s2;
  end
  rmin = r(n1, n2) / 2;
  ns = 0;
  n1p = n1 + (ns + 1) * s1;
  n2p = n2 + (ns + 1) * s2;
  while r(n1p, n2p) > rmin && ...
        all(mask((n1p - 1):(n1p + 4), (n2p - 1):(n2p + 4)));
    ns = ns + 1;
    n1p = n1 + (ns + 1) * s1;
    n2p = n2 + (ns + 1) * s2;
  end
  [~, c, o] = twosidedkernel(ns + 2);
  if m == 2
      c = c.';
      o = o.';
  end
  if s1 < 0
    n1 = n1 + ns * s1;
  end
  if s2 < 0
    n2 = n2 + ns * s2;
  end
  [cntx, orly, mask, Nc] = addkernel(cntx, orly, c, o, n1, n2, mask, Nc);
end

function [k, c, o] = resistencekernels()
  k = cell(2, 1);
  c = cell(2, 1);
  o = cell(2, 1);
  [k{1}, c{1}, o{1}] = resistencekernel();
  k{2} = k{1}.';
  c{2} = c{1}.';
  o{2} = o{1}.';
end

function [k, c, o] = resistencekernel()
  c = [1, 0, 1];
  o = [0, 0, 0];
  k = normalizekernel(c | o);
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

function [cntx, orly, mask, Nc] = addkernel(cntx, orly, c, o, n1, n2, mask, Nc)
  [N1, N2] = size(c);
  Nc = Nc + 1;
  cntx(n1:(n1 + N1 - 1), n2:(n2 + N2 - 1)) = c * Nc;
  orly(n1:(n1 + N1 - 1), n2:(n2 + N2 - 1)) = o;
  mask(n1:(n1 + N1 - 1), n2:(n2 + N2 - 1)) = 0;
end

function [r] = updatecorrs(r, k, mask)
  [row, col] = find(mask == 0);
  for m = 1:length(r)
    [N1r, N2r] = size(r{m});
    [N1k, N2k] = size(k{m});
    for l = 1:length(row)
      r{m}(max(0, row(l) - N1k):min(N1r, row(l) + 1), ...
           max(0, col(l) - N2k):min(N2r, col(l) + 1)) = 0;
    end
  end
end

function [cnxn] = addconnections(x, cntx, mask)
  [N1, N2] = size(x);
  v = directionalderivate(x, 1, 0, mask);
  h = directionalderivate(x, 0, 1, mask);
  du = directionalderivate(x, 1, 1, mask);
  dd = directionalderivate(x, 1, -1, mask);
  mask = dilatemask(cntx, mask);
  cnxn = zeros(N1, N2);
  while sum(sum(x .* (cntx > 0))) > 0
    [n1, n2] = maxcntx(x, cntx);
    cmp = cntx(n1, n2);
    cntx(n1, n2) = 0;
    cmpmask = removesamecmp(cntx, cmp, mask);
    [d1, d2] = getbestd(n1, n2, v, h, du, dd, cmpmask);
  end
end

function [dd] = directionalderivate(x, dx, dy, mask)
  [N1, N2] = size(x);
  dd = zeros(N1, N2);
  for n1 = 1:N1
    for n2 = 1:N2
      count = 0;
      if iswithinlimits(n1 - dx, N1) && iswithinlimits(n2 - dy, N2)
        dd(n1, n2) = dd(n1, n2) + abs(x(n1, n2) - x(n1 - dx, n2 - dy));
        count = count + 1;
      end
      if iswithinlimits(n1 + dx, N1) && iswithinlimits(n2 + dy, N2)
        dd(n1, n2) = dd(n1, n2) + abs(x(n1 + dx, n2 + dy) - x(n1, n2));
        count = count + 1;
      end
      if count > 0
        dd(n1, n2) = dd(n1, n2) / count;
      end
    end
  end
  dd = dd .* x .* mask / norm([dx, dy]);
end

function [mask] = dilatemask(cntx, mask)
  [N1, N2] = size(mask);
  for n1 = 1:N1
    for n2 = 1:N2
      numcntx = 0;
      for d1 = -1:1
        for d2 = -1:1
          if iswithinlimits(n1 + d1, N1) && iswithinlimits(n2 + d2, N2)
            numcntx = numcntx + (cntx(n1 + d1, n2 + d2) > 0);
          end
        end
      end
    end
  end
end

function [n1, n2] = maxcntx(x, cntx)
  cntxweights = x .* (cntx > 0);
  [~, n] = max(cntxweights(:));
  [n1, n2] = ind2sub(size(cntxweights), n);
end

function [cmpmask] = removesamecmp(cntx, cmp, mask)
  [N1, N2] = size(cntx);
  cmpmask = mask;
  [n1, n2] = find(cntx == cmp);
  for n = 1:length(n1)
    for d1 = -1:1
      for d2 = -1:1
        if iswithinlimits(n1(n) + d1, N1) && iswithinlimits(n2(n) + d2, N2)
          cmpmask(n1(n), n2(n)) = 0;
        end
      end
    end
  end
end

function [d1, d2] = getbestd(n1, n2, v, h, du, dd, mask)
  dl1 = [-1, -1, -1, 0, 0, 1, 1, 1];
  dl2 = [-1, 0, 1, -1, 1, -1, 0, 1];
  [d1, d2] = getbestdfromlist(n1, n2, dl1, dl2, v, h, du, dd, mask);
end

function [d1, d2] = getbestdfromlist(n1, n2, dl1, dl2, v, h, du, dd, mask)
  d1 = 0;
  d2 = 0;
  dmax = 0;
  for n = 1:length(dl1)
    d = getd(n1, n2, dl1(n), dl2(n), v, h, du, dd, mask);
    if d > dmax
      dmax = d;
      d1 = dl1(n);
      d2 = dl2(n);
    end
  end
end

function [d] = getd(n1, n2, d1, d2, v, h, du, dd, mask)
  [N1, N2] = size(mask);
  d = inf;
  if iswithinlimits(n1 + d1, N1) && iswithinlimits(n2 + d2, N2)
    if d1 == 0
      der = v .* mask;
    elseif d2 == 0
      der = h .* mask;
    elseif d1 * d2 > 0
      der = dd .* mask;
    else
      der = du .* mask;
    end
    d = der(n1 + d1, n2 + d2);
  end
end

function [b] = iswithinlimits(x, N)
  b = (x > 0) && (x <= N);
end
