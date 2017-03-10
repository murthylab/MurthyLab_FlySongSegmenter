function s = deblankl(x)
%DEBLANKL Convert string to lowercase without blanks.
%   S = DEBLANKL(X) is the string X converted to lowercase 

if ~isempty(x)
    s = lower(x);
    s = s(s~=' ');
else
    s = [];
end
