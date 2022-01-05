config.clean_string = '$area == "PFC" | $area == "CA1"';
cellTable = util.table.query(spikes.cellTable, config.clean_string);
config.sig = 0.05/height(spikes.cellTable);


                              
% ,---.,   .,---.|    ,---.,---.   
% |---||\  ||  _.|    |--- `---.
% |   || \ ||   ||    |        |
% `   '`  `'`---'`---'`---'`---'
                              
ca1 = spikes.cellTable.area == "CA1";
pfc = spikes.cellTable.area == "PFC";
neurons.set =  1:height(spikes.cellTable);
neurons.N =  height(spikes.cellTable);
neurons.Npfc =  sum(pfc);
neurons.Nca1 =  sum(ca1);


% ------ OVERALL -----------------------------
% Q1 : What percent of neurons are angular?
sigNeurons.all = any(sarel.stops.rayleigh.currentAngle.pval < config.sig,2);
disp(sum(sigNeurons.all))
% Q2 : How many goals per neuron?
sigNeurons = sarel.stops.rayleigh.currentAngle.pval < config.sig;
disp(sum(sigNeurons))
% Q3 : When a neuron fires to a goal, same angle across goals or differnt>?
sigNeurons = sarel.stops.rayleigh.currentAngle.pval < config.sig;
disp(sum(sigNeurons))


% ------ REGIONAL ----------------------------
% Q1 : What percent of neurons are angular?
sigNeurons.ca1 = any(sarel.stops.rayleigh.currentAngle.pval(ca1,:) < config.sig,2);
sigNeurons.pfc = any(sarel.stops.rayleigh.currentAngle.pval(pfc,:) < config.sig,2);
disp(sum(sigNeurons.ca1))
disp(sum(sigNeurons.ca1)/neurons.Nca1)
disp(sum(sigNeurons.pfc))
disp(sum(sigNeurons.pfc)/neurons.Npfc)
% Q2 : How many goals per neuron?
sigNeurons.goal.all = sum(sarel.stops.rayleigh.currentAngle.pval < config.sig,2);
sigNeurons.goal.ca1 = sum(sarel.stops.rayleigh.currentAngle.pval(ca1,:) < config.sig,2);
sigNeurons.goal.pfc = sum(sarel.stops.rayleigh.currentAngle.pval(pfc,:) < config.sig,2);
disp(sum(sigNeurons))
% Q3 : When a neuron fires to a goal, same angle across goals or differnt>?
sigNeurons = sarel.stops.rayleigh.currentAngle.pval < config.sig;
disp(sum(sigNeurons))
% Q4 : Cells with high rayleigh for multiple goals share angles across goals? Or different?

                                    
%,--. |,---.--.--,---.,   .,---.,---.
%|   ||`---.  |  |---||\  ||    |--- 
%|   ||    |  |  |   || \ ||    |    
%`--' ``---'  `  `   '`  `'`---'`---'
%Distance


