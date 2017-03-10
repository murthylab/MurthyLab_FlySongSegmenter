function all_ha = gcas(varargin)
% get handles to the axes of all subplots in a figure
all_ha = findobj(gcf, 'type', 'axes', 'tag', '' );
if nargin==0
   return
elseif nargin
   all_ha = all_ha(varargin{1});
end