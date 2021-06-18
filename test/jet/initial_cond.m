function [rho] = initial_cond
nlat   = 2000;
nz     = 100;
lat_c  =     30;
H      =  -4000;
width  =   2000;
L_jet  =   1600;
radius =  10825;

rho_0  =  27.75;
S_b    =   9.8e-6;
z_surf =   -300;

drho_surf = [0 1.5];
dz        = [300 700];
z_int     = [-400 -1000];

lat_width = rad2deg (width/radius);

lat = lat_c - lat_width/2 : lat_width/nlat : lat_c + lat_width/2;
z   = H : abs(H)/nz : 0;
rho = zeros(nz+1, nlat+1);
sm = zeros(nlat+1);
for i = 1:nlat+1
    sm(i) = smoothing (lat(i));
    for k = 1:nz+1
        rho(k,i) = density_ic (lat(i), z(k));
    end
end
c1 = [28 27.5 27 26.5 26 25.5 25 24.5 24];
caxis([min(c1) max(c1)]);
contourf(lat,z,rho,c1);colormap(autumn);colorbar;
figure;plot(rho(:,1000),z,'linewidth',2);grid on;ylabel('z');xlabel('rho');

    function [rho] = density_ic(phi, depth)
        drho_N = rho_NS (1, depth);
        drho_S = rho_NS (2, depth);
        
        %rho = rho_0 - S_b * (depth - H) + smoothing(phi) * drho_S + drho_N;
        rho = rho_0 - S_b * (depth - H) + smoothing(phi) * drho_S + (1-smoothing(phi))*drho_N;
        %rho = rho_0 - S_b * (depth - H) + drho_N;
    end

    function rho = rho_NS (hemi, depth)
        rho = - 0.5 * drho_int(hemi) * (1 + tanh((d_NS(depth,hemi) - z_int(hemi))/dz(hemi))) ...
            - drho_surf(hemi) / (2*tanh(1)) * (1 + tanh ((z_surf - depth) / z_surf));
    end

    function d = d_NS (depth, hemi)
        d = z_int (hemi) + (depth - z_int (hemi)) ...
            * sqrt (1 + 0.5 * ((depth - z_int (hemi) + abs (depth - z_int (hemi))) ...
            / (1.3*dz(hemi)))^2);
    end

    function  drho = drho_int(hemi)
        if hemi == 1
            %drho = 1.41;
            drho = 1.40 * (1 + tanh((d_NS(0,2)-z_int(2))/dz(2))) / ...
                (1 + tanh((d_NS(0,1)-z_int(1))/dz(1)));
        else
            drho = 1.40;
        end
    end

    function sm = smoothing (phi)
        y = pi * (width/L_jet * (phi - (lat_c - lat_width/2)) / lat_width + 0.5 * (1 - width/L_jet));
        
        if y < 0
            sm = 1;
        elseif y > pi
            sm = 0;
        else
            sm = 1 - (y - sin (y) * cos (y)) / pi;
        end
    end
end

