function [ysmooth] = mvg(y)
% 5-point moving average

ysmooth = zeros(length(y),1);

ysmooth(1)      = y(1);
ysmooth(2)      = (y(1) + y(2) + y(3))/3;
for i = 3:length(y)-2
    ysmooth(i) = (y(i-2) + y(i-1) + y(i) + y(i+1) + y(i+2))/5;
end
ysmooth(end-1)  = (y(end-2) + y(end-1) + y(end))/3;
ysmooth(end)    = y(end);