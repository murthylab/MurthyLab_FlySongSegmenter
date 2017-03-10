function [fZ, scale] = scaleZ2M(Z,M,t)
%USAGE
%fZ = scaleZ2M(Z,M)
[n_samples,total_length] = size(Z);
%rescale data to mean
%equivalent to following, on whole array
%a = mean(M.*Z(n,:))/mean(Z(n,:).^2);
%Z(n,:) = a*Z(n,:);

Ma = repmat(M,n_samples,1);
num = mean(Ma'.*Z');
den = mean(Z'.^2);
a = num./den;
ar = repmat(a',1,total_length);
Z = ar.* Z;

MZ = mean(Z,1);
% (MZ/M)
% if mxType ~= 0
%     scale = max([(MZ/M) 0.45]);
% elseif mxType == 2
%     scale = max([(MZ/M) 0.15]);
% else
scale = max([(MZ/M) 0.30]);
% end
% scale = max([0.5 (MZ/M)]);

fZ = Z/scale;