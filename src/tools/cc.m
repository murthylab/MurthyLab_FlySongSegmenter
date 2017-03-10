try % to make it work in case there's no graphics system
   if ~isempty(findall(0,'Type','Figure'))
      clf
%       close all
   end
end
clear
clc
clear global