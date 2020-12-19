function make_in_file(FileNamePrefix, fpn, f, zs, zrr, r, bty, z, c)

% make file or erase current contents of file:
fid = fopen(fpn, 'w');
fclose(fid);

% open file to append contents
fid = fopen(fpn, 'at');

% line 1: Prefix
fprintf(fid, FileNamePrefix);
fprintf(fid, '\n');

% line 2: f, zs, zrr (don't know purpose of zs and zrr) 
% fprintf(fid, '300\t!Freq (Hz)\n');
fprintf(fid, '%.6f\t%.6f\t%.6f\t\tf, zs, zrr\n', f, zs, zrr)

% line 3: rmax, dr, ndr (not clear on purpose of ndr)
fprintf(fid, '%.6f\t%.6f\t%d\t\trmax, dr, ndr\n', max(r), 1, 1)

% line 4: zmax dz, ndz, zmplot (don't know purpose of ndz or zmplot)
%fprintf(fid, '%.6f\t%.6f\t%d\t%.6f\t\tzmax, dz, ndz, zmplot\n', max([max(bty), max(20)]), 1, 1, 805)
fprintf(fid, '%.6f\t%.6f\t%d\t%.6f\t\tzmax, dz, ndz, zmplot\n', 930, 1, 1, 805)

% line 5: c0, np, irot, theta  (don't know purpose of any of these)
fprintf(fid, '%.6f\t%d\t%d\t%.6f\t\tc0, np, irot, theta\n', 1500, 8, 0, 0)

%% Bathymetry
for nb = 1:length(bty)
    fprintf(fid, '%.6f\t%.6f\n', r(nb), bty(nb));
end

% End of bathymetry
fprintf(fid, '-1\t-1\n')

%% Sound speed profile
for nz = 1:length(c)
    fprintf(fid, '%.6f\t%.6f\n', z(nz), c(nz));
end

% End of SSP
fprintf(fid, '-1\t-1\n')

%% compressive SSP
zcomp = [0, 90, 200, 350];
ccomp = [1500, 1700, 1850, 1850];
for nz = 1:length(zcomp)
    fprintf(fid, '%.6f\t%.6f\n', zcomp(nz), ccomp(nz))
end

% End of compressive SSP
fprintf(fid, '-1\t-1\n')

%% shear SSP
zsh = [0, 90, 200, 350];
csh = [10, 24, 72, 72];
for nz = 1:length(zcomp)
    fprintf(fid, '%.6f\t%.6f\n', zsh(nz), csh(nz))
end

% End of shear SSP
fprintf(fid, '-1\t-1\n')

%% Density profile in substrate
zd = [0, 90, 200, 350];
rhob = [1.35, 1.5, 1.7, 1.7];
for nz = 1:length(zcomp)
    fprintf(fid, '%.6f\t%.6f\n', zd(nz), rhob(nz))
end

% End of density
fprintf(fid, '-1\t-1\n')

%% compressive attenuation profile in substrate
zac = [0, 90, 200, 350];
attnp = [.2, 1, .5, 10];

for nz = 1:length(zcomp)
    fprintf(fid, '%.6f\t%.6f\n', zac(nz), attnp(nz))
end

% End of compr. atten.
fprintf(fid, '-1\t-1\n')

%% shear attenuation profile in substrate
zas = [0, 90, 200, 350];
attns = [1, 1.5, 1.5, 10];

for nz = 1:length(zcomp)
    fprintf(fid, '%.6f\t%.6f\n', zas(nz), attns(nz))
end

% End of shear atten.
fprintf(fid, '-1\t-1\n')

%% END OF FILE
fclose(fid);
end