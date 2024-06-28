function [X]=str2num_mjk(str)
[X,OK]=str2num(str);
X(OK==0)=NaN;
end