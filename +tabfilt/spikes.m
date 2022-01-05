function s = spikes(s, filtstr, varargin)
% Filter the spikes at behavior struct with a filter
% string

ip = inputParser;
ip.addParameter('lean', true);
ip.parse(varargin{:})
Opt = ip.Results;

filtstr = string(filtstr);
filtstr = replace(filtstr, '$', 's.beh.');

filt = true(height(s.beh), 1);
for f = filtstr(:)'
   filt = filt & eval(f);
end

s.beh = s.beh(filt,:);

if Opt.lean
    fieldset = setdiff(string(fieldnames(s)), "beh");
    s = rmfield(s, fieldset);
end
