function varargout = mapFun(fun,x, args)
% [x1, x2,...] = mapFun(fun,x,args)
% applies function FUN to each plane of a stack / column(?) of a matrix
% PARAM
%  FUN - handle to function
%   x  - 2d/3d array
% args - array for passing single arguments or cell array for passing multiple arguments
% RETURNS
%  x1..xN - variable number of output arguments
% see also BSXFUN, ARRAYFUN, CELLFUN

% TODO: make parallel - however, current implementation doesn't work well with PARFOR

%% some rudimentary error checking
% check if first argument is a function handle
%assert(isa(fun, 'function_handle'), 'First argument needs to be a valid function handle.') 
% check if the handle points to an actual function
%assert(exist(func2str(fun))>2, 'First argument needs to point to a valid function')
% check if the function returns the requested number of output arguments
%assert(nargout<=nargout(fun), ['Too many output requested: The function ' func2str(fun) ' returns ' int2str(nargout(fun)) '. You requested ' int2str(nargout) '.'])
   
%%
if nargin==2
   args = {};
end

if ~iscell(args)
   args = {args};
end

if ndims(x)==3
   for i = size(x,3):-1:1
      [tmp{1:nargout}] = fun(x(:,:,i), args{:});
      y(:,:,i) = fun(x(:,:,i), args{:});
      for arg = 1:nargout
         varargout{arg}(:,:,i) = tmp{arg};
      end
   end
elseif ndims(x)==2
   for i = size(x,2):-1:1
      %y(:,i) = fun(x(:,i), args{:});
      [tmp{1:nargout}] = fun(x(:,i), args{:});
      for arg = 1:nargout
         varargout{arg}(:,i) = tmp{arg};
      end
   end
else
   error('first input must be either 2D or 3D matrix!!')
end
