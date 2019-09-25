function [t, x, y, z]= auxEvalModel(m, traj, nt, nx, ny)
%function z= auxEvalModel(m, traj, t, x, y)
%
% t     real time [sec], starting at 0 
%
%output
%   z   z{traj}{t}(x,y);


ntraj= length(traj);
tiny= 1e-7;
%nx= length(x);
%ny= length(y);

tp= 0.5; % tension parameter
M= [-tp, 2-tp, tp-2, tp; 2*tp, tp-3, 3-2*tp, -tp; -tp, 0, tp, 0; 0, 1, 0, 0];

switch m.name
case 'CSplines'
%    keyboard
    if(m.outputCompress) 
        theta= double(m.parameters.param)/m.parameters.convertFac + ...
        m.parameters.convertMin;
    else
        theta= double(m.parameters.param);
    end
    tstep= m.outputInterval;
    tmax= tstep*size(theta,2);
    t= linspace(tstep,tmax,nt);
    
    y= {};
    for itraj=1:ntraj   % loop over requested trajectories
        cpx= m.cpx{itraj}; ncpx= length(cpx);
        
        x= linspace(cpx(2)+tiny, cpx(ncpx-1)-tiny, nx);

        for it=1:nt     % loop over requested snapshot times
            tindex= ceil(t(it)/tstep);

            % evaluate bicubic cardinal spline in z
            y{itraj,it}= zeros(nx,1);
            for ix=1:nx
                icpx= max(find(cpx < x(ix)));
                ux= (x(ix) - cpx(icpx))/ (cpx(icpx+1) - cpx(icpx));
                C(:,1)= theta(icpx-1:icpx+2, 1);
                y{itraj,it}(ix)= [ux^3 ux^2 ux 1]*M*C;
            end
        end    % for it
    end % for itraj
case 'CardSplines2d'
    if isfield(m, 'periodic') & strcmp(m.periodic, 'y')
        periodicY= 1;
%        disp('periodic');
    else
        periodicY= 0;
    end

    tstep= m.outputInterval;
    tmax= tstep*size(m.theta,2);
    t= linspace(tstep,tmax,nt);
    
    zMin= m.convertMin;
    zFac= m.convertFac;
    for itraj=1:ntraj   % loop over requested trajectories
        cpx= m.cpx{itraj}; ncpx= length(cpx);
        cpy= m.cpy{itraj}; ncpy= length(cpy);
        ncpxy= ncpx*ncpy;
        
        % check whether x and y are within bounds
%        if sum(x < cpx(2)) |  sum(x > cpx(ncpx-2)) | ...
%            sum(y < cpy(2)) |  sum(y > cpy(ncpy-2))  ...
%                error('control points out of bounds');
%        end
        x= linspace(cpx(2)+tiny, cpx(ncpx-1)-tiny, nx);

        if periodicY
            y= linspace(cpy(1)+tiny, cpy(ncpy)-tiny, ny);
        else
            y= linspace(cpy(2)+tiny, cpy(ncpy-1)-tiny, ny);
        end

        for it=1:nt     % loop over requested snapshot times
            tindex= ceil(t(it)/tstep);
%            keyboard
            lintheta= double(m.theta(ncpxy*traj(itraj)+1:ncpxy*(traj(itraj)+1),tindex))./zFac + zMin;
%            theta= reshape(lintheta, ncpx, ncpy);
            % data was saved in row major format
            theta= reshape(lintheta, ncpy, ncpx)';

%            subplot(4,1,itraj);
%            imagesc(cpx,cpy,theta);
%            axis xy
%            colorbar
%            break

            % evaluate bicubic cardinal spline in z
            z{itraj,it}= zeros(nx,ny);
            for ix=1:nx
                icpx= max(find(cpx < x(ix)));
                ux= (x(ix) - cpx(icpx))/ (cpx(icpx+1) - cpx(icpx));
                for iy=1:ny
                    icpy= max(find(cpy < y(iy)));
                    uy= (y(iy) - cpy(icpy))/ (cpy(icpy+1) - cpy(icpy));

                    % find control point values
                    if periodicY
                        ylo= icpy-1; yhi= icpy+2;
                        if ylo==0
                            C(:,2:4)= theta(icpx-1:icpx+2, 1:yhi);
                            C(:,1)= theta(icpx-1:icpx+2, ncpy-1);
                        elseif yhi >= ncpy
                            yover= yhi-ncpy+1;
                            C(:,1:4-yover)= theta(icpx-1:icpx+2, ylo:ncpy-1);
                            C(:,4-yover+1:4)= theta(icpx-1:icpx+2, 1:yover);
                        else
                            C= theta(icpx-1:icpx+2, icpy-1:icpy+2);
                        end
                    else
                        C= theta(icpx-1:icpx+2, icpy-1:icpy+2);
                    end

                    z{itraj,it}(ix,iy)= ...
                                [uy^3 uy^2 uy 1]*M*([ux^3 ux^2 ux 1]*M*C)';
                end
            end
        end    % for it
    end % for itraj
otherwise
    error(['unknown model ' m.name]);
end

