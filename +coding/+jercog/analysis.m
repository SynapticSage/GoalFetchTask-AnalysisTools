function spikes = analysis(animal, day, spikes, Opt)

% Get sparse behavior at spike times for the next few analyses
beh = gbehavior.lookup(animal, [], day);
[spikes, beh, ~] = units.atBehavior(beh, spikes,...
    'useGPU', false,....
    'merge', true,...
    'query', Opt.behFilter);
Opt.checkpointVar = 'jercogtmp';
spikes = jercog(spikes, beh, Opt);
