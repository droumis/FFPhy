% GABRIELR: Investigate edge lengths of random graphs

function gabrielr(N,distrib_size)

  distrib = zeros(distrib_size,1);
  nb = 1;
  ne = 2;

  while (ne < distrib_size)
    crds = rand(N,2);
    [connect,dist] = gabriel(crds,1);
    t = trilow(dist);
    d = t(t>0);
    lend = length(d);

    ne = nb+lend-1;
    if (ne > distrib_size)
      ex = ne-distrib_size;
      lend = lend - ex;
      ne = distrib_size;
    end;
    distrib(nb:ne) = d(1:lend);
    nb = ne+1;
  end;

  figure;
  histgram(distrib);

  return;
