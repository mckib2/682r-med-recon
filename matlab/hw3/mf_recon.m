function [ im_mf ] = mf_recon(d,k,w,n,te,tad,fmin,fmax,fstep,opt)

    if nargin < 10
        opt = 'simple';
    end

    osf = 2; wg = 5;
    idy = round((.5*n*(osf-1)+1):(.5*n*(osf+1)));
    idx = fliplr(idy);

    [ ns,ni ] = size(d);
    t = te + (0:ns-1)*tad/ns;
    t = repmat(t.',1,ni);

    fs = fmin:fstep:fmax;
    im_mf = zeros(n,n,numel(fs));

    % Discrete Frequency Algorithm
    if strcmp(opt,'simple')
        for ii = 1:numel(fs)
            tmp = gridkb(d.*exp(-1j*2*pi*fs(ii)*t),k,w,n,osf,wg,'image');
            im_mf(:,:,ii) = tmp(idx,idy);
        end
    elseif strcmp(opt,'phasefun')
        krt = abs(k);
        [ m_hat,~,da ] = gridkb(d,k,w,n,osf,wg,'k-space');

        for ii = 1:numel(fs)
            tmp = m_hat.*P(fs(ii),krt,size(m_hat,1));
            tmp = fftshift(fft2(fftshift(tmp)))./da;
            im_mf(:,:,ii) = tmp(idx,idy);
        end
    end
end

function out = P(f,krt,n)
    kx = linspace(min(min(krt)),max(max(krt)),n);
    [ KX,KY ] = meshgrid(kx,kx);
    KR = sqrt(KX.^2 + KY.^2);
    
    tkr = interp1(krt(:,1),1:size(krt,1),KR,'linear');
    tkr(isnan(tkr)) = 1;
        
    out = exp(1j*2*pi*f*tkr);
end
