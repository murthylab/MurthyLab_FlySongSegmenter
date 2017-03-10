function f = scal2frq(a,wname,delta)
%SCAL2FRQ Scale to frequency.
if isempty(a) , f = a; return; end
if nargin == 2, delta = 1; end
err = (min(size(a))>1) | (min(a)<eps);
if err
    error('Wavelet:FunctionArgVal:Invalid_ArgVal', ...
        'Invalid Value for Scales a !')
end
if delta <= 0
    error('Wavelet:FunctionArgVal:Invalid_ArgVal', ...
        'Invalid Value for Delta !')
end

% Compute pseudo-frequencies
f = centfrq(wname)./(a.*delta);     
