% ROTXY: Variant of PCA: rotates the x/y plane so as to maximize variation in z.
%         Spatial structure.
%

function fac = rotxy(X)
  [n,p] = size(X);

  x = X(:,1);
  y = X(:,2);
  z = X(:,3);
  
if (1)
  X = [x y ones(n,1)];               % Multiple regression of z on (x,y)
  Sxx = X'*X./(n-1);
  sxy = X'*z./(n-1);
  b = inv(Sxx)*sxy;                  
  R2 = var(X*b)/var(z);
R2
end;

fac = pcacorr(corrcoef([x y z]));       % Get initial loadings
fac

%for theta = linspace(0,2*pi,20)
for i = 1:10
i
  cos2theta = fac(3,1)^2 - fac(3,2)^2
  sin2theta = 2*fac(3,1)*fac(3,2)
acos2theta = acos(cos2theta)
asin2theta = asin(sin2theta)
  lae = sqrt(cos2theta^2 + sin2theta^2)
  costheta =  sqrt(1-cos2theta/lae)/2
  sintheta = -sqrt(1+cos2theta/lae)/2
  fac(:,1) = fac(:,1)*costheta - fac(:,2)*sintheta;
  fac(:,2) = fac(:,2)*sintheta + fac(:,2)*costheta;
fac
end;


if (0)
  fac = pcacorr(corrcoef(X));       % Get initial loadings

  sin1 = 1;
  while (abs(sin1)>0.0001)
    cos2sum = sum(fac(:,1).^2 - fac(:,2).^2);
    sin2sum = sum(2*fac(:,1).*fac(:,2));

    lae = sqrt(cos2sum^2 + sin2sum^2);
    cos1 =  0.5*sqrt(1-cos2sum/lae);
    sin1 = -0.5*sqrt(1+cos2sum/lae);

    fx = fac(:,1)*cos1 - fac(:,2)*sin1;
    fy = fac(:,1)*sin1 + fac(:,2)*cos1;
    fac(:,1:2) = [fx fy]
  end;
end;
  

  return;