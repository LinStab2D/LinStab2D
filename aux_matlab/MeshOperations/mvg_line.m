function []=mvg_line(val)

lo = gco;
lo.YData = smooth(lo.YData,val);