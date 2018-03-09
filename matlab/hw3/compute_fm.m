function [ fm ] = compute_fm(im1,te1,im2,te2)
    dphi = angle(im2.*conj(im1));
    dt = abs(te1 - te2);
    w = dphi/dt;
    fm = 1/(2*pi)*w;
end