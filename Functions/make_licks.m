


Fp.animals = {'JZ3'};
Fp.add_params = {'wtrackdays', 'lickbouts', 'valid_ntrodes', 'nonMU_cells'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
    'excludetime', Fp.timefilter);


for ani = 1:length(F)
    animdef = animaldef(F(ani).animal{3});
    DIO = loaddatastruct(animdef{2}, F(ani).animal{3}, 'DIO');
    task = loaddatastruct(animdef{2}, F(ani).animal{3}, 'task');
    get_licks(F(ani).animal{3}, F(ani).epochs{1}, DIO, task)
end
