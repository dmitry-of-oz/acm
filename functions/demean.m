function Yd = demean(Y)

Yd      = bsxfun(@minus, Y, mean(Y));