function out_snippet = interp_snippet(in_snippet, up_rate)
% function out_snippet = interp_snippet(in_snippet, up_rate)


out_snippet = zeros(40*up_rate, 4);

for i = 1:4
  fi = fft(double(in_snippet(:,i)));
  fo = zeros(40*up_rate, 1);
  fo(1:20) = fi(1:20);
  fo(end-18:end) = fi(22:end);
  out_snippet(:,i) = up_rate * ifft(fo);
end


