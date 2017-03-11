function index = argsort(varargin)
%ARGSORT   Index of sorted elements.
% ARGSORT(X) returns an index I such that X(I) is sorted.
%
% See also SORT, ARGMIN, ARGMAX.

[~,index] = sort(varargin{:});
