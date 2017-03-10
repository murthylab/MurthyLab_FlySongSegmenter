function y = wconv1(x,f,shape)
%WCONV1 1-D Convolution.
%   Y = WCONV1(X,F) performs the 1-D convolution of the 
%   vectors X and F.
%   Y = WCONV1(...,SHAPE) returns a subsection of the
%   convolution with size specified by SHAPE (See CONV2).


if nargin<3
    shape = 'full';
end
y = conv2(x(:)',f(:)',shape); 
if size(x,1)>1
    y = y';
end
