function[L,d,B]=mresp(Z,nx)
% MRESP.M
%function[L,d,B]=mresp(Z,nx)
% L normalised eigen vector matrix of vn'  L= ny.ny
% d eigenvalues in col vector
% B coeff matrix nx.ny
% Z =[X ,Y]  X= indep var, Y= dep var
% nx=no.  of cols containing indep variables
% calc of B in Y=X*B + E , corr(B), 95%CL on B
% Also calcs eigen vectors /values for Y=X*B + E
% for max predictability i.e. L.Y' most predictable L.B' best predictor
% Y=n.ny, X=n.nx, B=nx.ny, L=ny.ny ,# of parameters for each y is nx
      [n,p]=size(Z);
      X=Z(:,1:nx);
      Y=Z(:,nx+1:p);
      B=X\Y ;       %solve for B in Y=X*B + E
      disp(' Coefficient Matrix B in Y= X*B + E ,  B=nx.ny')
      disp(B)
      pause
      xm=mean(X);
      XMEAN=X - ones(n,1)*xm;
      PSI= XMEAN'*XMEAN;
      BPB=B'*PSI*B;
      E=Y-X*B;
      ETE=E'*E;
% extra commands to calc corr(B),(%)CL
      xtxinv=inv(X'*X);
      corrb=zeros(nx,nx);
      for k=1:nx
      for j=k:nx
      den=sqrt(xtxinv(k,k)*xtxinv(j,j));
      corrb(k,j)=xtxinv(k,j)/den;
      corrb(j,k)=corrb(k,j);
      end
      end
      disp(' correlation matrix for coeficient matrix B')
      disp(corrb)
      pause
      % calc diag of std dev for each y
      dete=sqrt(diag(ETE)); % ny.1
      dxtxinv=sqrt(diag(xtxinv)); % nx.1
      p=dxtxinv*dete'*2/sqrt(n-nx); % nx.1 * 1.ny = nx.ny
      disp(' 95%CL for B  B= B +/- elements in this matrix')
      disp(p)
      pause
% Eigenvalue Calcs
      disp(' eigen values ' )
      [v,di]=eig(BPB,ETE);
      d=diag(di);
      disp(d)
      pause
      maxvec=max(abs(v));
      disp(' L= most predictable linear combinations of Y')
      disp(' normalised eigenvector Matrix L  (rows)')
      L=(v./(ones(n-nx,1)*maxvec))';
      disp(L)
      pause
      disp('LB(trans) = linear combinations of X that are best Predictors')
      disp(L*B')
      return

