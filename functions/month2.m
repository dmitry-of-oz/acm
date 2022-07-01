function month1 = month2(arg1)

b       = datestr(arg1);
c       = b(:,4:6)';

c       = c(:)';
month1      = nan(size(arg1,1),1);
allmonths   = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
                'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

for i=1:12
    month1((strfind(c,allmonths{i})-1)/3+1) = i;
end

