function beh = behavior(beh, filtstr)
% Filter the spikes at behavior struct with a filter
% string

filtstr = string(filtstr);
filtstr = replace(filtstr, '$', 'beh.');

filt = true(height(beh), 1);
for f = filtstr(:)'
   filt = filt & eval(f);
end

beh = beh(filt,:);
