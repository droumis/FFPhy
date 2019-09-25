files = dir('*DIO*.mat');

if ~exist('pin')
  pin = 48;
  % pin = 15;
end

pot_threshold = 20000;

for i = 1:length(files)
   ff = files(i).name;
   load(ff);
   for j = 1:length(DIO)-1
      if ~isempty(DIO{j})
         error('Multiple non-empty DIO cells!');
      end
   end
   for j = 1:length(DIO{end})
      if ~isfield(DIO{end}{j}{pin},'timesincelast')
         if length(DIO{end}{j}{pin}.pulsetimes) > 0
            DIO{end}{j}{pin}.timesincelast = diff(DIO{end}{j}{pin}.pulsetimes(:,1));
            DIO{end}{j}{pin}.timesincelast = [Inf; DIO{end}{j}{pin}.timesincelast];
         else
            DIO{end}{j}{pin}.timesincelast = [];
         end
      end
      frequency = 10000./DIO{end}{j}{pin}.timesincelast;
      frequency2 = frequency;
      if length(frequency) > 1
         frequency = [frequency(2:end); nan];
      end
      DIO{end}{j}{pin}.frequency = frequency;
      index_in_sequence = -ones(length(DIO{end}{j}{pin}.timesincelast),1);
      first_pulse = -1;
      for k = 1:length(DIO{end}{j}{pin}.timesincelast)
         if (DIO{end}{j}{pin}.timesincelast(k) < pot_threshold)
            if (first_pulse < 0)
               first_pulse = k;
            else
               frequency2(k) = frequency(first_pulse-1);
            end
            index_in_sequence(k) = index_in_sequence(k-1) + 1;
         else
            first_pulse = -1;
         end
      end
      pot_inds = find(index_in_sequence >= 0);
      index_in_sequence(pot_inds) = index_in_sequence(pot_inds) + 1;
      pot_starts = find(index_in_sequence == 1);
      index_in_sequence(pot_starts - 1) = 0;
      DIO{end}{j}{pin}.index_in_sequence = index_in_sequence;
      DIO{end}{j}{pin}.frequency2 = frequency2;
   end
   save(ff,'DIO','rawdio','diopulses');
end

