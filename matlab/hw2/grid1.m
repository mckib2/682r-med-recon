function [ m ] = grid1(d,k,w,n,o,kw)

% function m = grid1(d,k,n)
%     d -- k-space data
%     k -- k-trajectory, scaled -0.5 to 0.5
%     w -- preweighting function
%     n -- image size
%     o -- oversampling factor
%    kw -- kaiser bessel kernel parameters, [ W,beta ]

% convert to single column
d = d(:);
k = k(:);

if ~isempty(w)
    w = w(:);
else
    w = ones(size(d));
end

if isempty(o)
    o = 1;
end

% Oversampling adjustment
n = n*o;

% cconvert k-space samples to matrix indices
nx = (n/2+1) + n*real(k);
ny = (n/2+1) + n*imag(k);

% zero out output array
m = zeros(n,n);

% loop over samples in kernel
for lx = -o:o
  for ly = -o:o

    % find nearest samples
    nxt = round(nx+lx);
    nyt = round(ny+ly);

    % compute weighting for triangular kernel
    kwx = max(1-abs(nx-nxt),0);
    kwy = max(1-abs(ny-nyt),0);

    % map samples outside the matrix to the edges
    nxt = max(nxt,1); nxt = min(nxt,n);
    nyt = max(nyt,1); nyt = min(nyt,n);

    % use sparse matrix to turn k-space trajectory into 2D matrix
    m = m + sparse(nxt,nyt,d.*kwx.*kwy.*w,n,n);
  end
end

% zero out edge samples, since these may be due to samples outside
% the matrix
m(:,1) = 0; m(:,n) = 0;
m(1,:) = 0; m(n,:) = 0;


