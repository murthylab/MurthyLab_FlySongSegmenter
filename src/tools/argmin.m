function index = argmin(varargin)
%ARGMIN  Index of maximum element.
% ARGMIN(X) returns an index I such that X(I) == MIN(X(:)).
%
% See also MIN, ARGMAX.

[~,index] = min(varargin{:});
