% Nonparaxial Fourier 3D Beam-Propagation Simulation in Inhomogeneous Medium
% Copyright (C) 2024 by Dr. Dashan Dong (dongdashan@icloud.com)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% Last Modified: 2024/04/15

% This script is used to simulate the propagation of a Gauss beam in a
% inhomogeneous medium using the Beam Propagation Method (BPM). BPM is
% a class of algorithms designed for calculating the optical field
% distribution in space or in time given initial conditions.
%
% Here the complex envelope of the optical field is calculated, which
% includes both the amplitude and phase of the field.
%
% Please note that this script is only a simple demonstration of the BPM
% algorithm and may not be suitable for all applications. For more
% information on BPM, please refer to the following literatures.
%
% 1. M.D. Feit and J.A. Fleck, "Beam nonparaxiality, filament formation,
% and beam breakup in the self-focusing of optical beams," J. Opt. Soc. Am.
% B 5, 633-640 (1988).
% 2. A. Taflove and S.C. Hagness, "Computational Electrodynamics: The
% Finite-Difference Time-Domain Method," 3rd ed. (Artech House, 2005).
% 3. J.W. Goodman, "Introduction to Fourier Optics," 3rd ed. (Roberts and
% Company, 2005).
% 4. A. Yariv and P. Yeh, "Photonics: Optical Electronics in Modern
% Communications," 6th ed. (Oxford University Press, 2007).

%% Initialization
clear;
clc;

%% Propagation Field Parameters
N_x = 512;
N_y = 512;

d_x = 50e-6 / 512;
L_x = d_x * N_x;		% X_Length of the field
Axis_x = linspace( -L_x/2, L_x/2, N_x );

d_y = 50e-6 / 512;
L_y = d_y * N_y;		% Y_Length of the field
Axis_y = linspace( -L_y/2, L_y/2, N_y );

[Grid_x, Grid_y] = meshgrid( Axis_x, Axis_y );
[Grid_th, Grid_r] = cart2pol( Grid_x, Grid_y );

d_z = 100e-9;	% Pixel Size in Z direction
z_o = -5e-6;			% Distance between Field_origin_z & Field_input_z
z_e = 15e-6;			% Distance between Field_input_z & Field_end_z
L_z = z_e - z_o;		% Z_Length of the field
N_z = floor( L_z / d_z ) + 1;
Axis_z = linspace( z_o, z_o + d_z * (N_z - 1), N_z );

[X, Y, Z] = meshgrid( Axis_x, Axis_y, Axis_z );
[TH, PHI, R] = cart2sph( X, Y, Z );

%% Propagation Environment Parameters
n_medium = 1.33;	% Refractive Index of immerse medium
n_sphere = 1.40;	% Refractive Index of sphere ball
n = n_medium * ones( N_y, N_x, N_z );
n( R<2.5e-6 ) = n_sphere;
n_delta = n - n_medium;
n_0z = squeeze( n( floor(N_x/2), floor(N_y/2), : ) );

%% Frequency Domain
dk_x = 2*pi / L_x;
dk_y = 2*pi / L_y;

k_x = dk_x .* ((1:1:N_x) - ceil((N_x+1)/2));
k_y = dk_y .* ((1:1:N_y) - ceil((N_y+1)/2));

[K_x, K_y] = meshgrid( k_x, k_y );
[K_th, K_r] = cart2pol( K_x, K_y );
K_r = ifftshift( K_r );

%% Beam Parameters
lambda = 561.0e-9;		% Wavelength
k0 = 2 * pi / lambda;
k = 2 * pi * n_medium / lambda;	% Wavevector
theta = 0.0;			% Propagation angle (in degree)
phi = 0.0;				% Rotation angle (in degree)

waistRadius = 15e-6;    % Beam radius for Gauss beam

%% Input & Output field Initialization

N_start = floor(-z_o / d_z)+1;

Section_0 = exp( -(Grid_r / waistRadius).^2 + 1i * sind( theta ) * ( cosd( phi ) * k * Grid_x + sind( phi ) * k * Grid_y ) );		% Input complex field

Section_yx = complex(Section_0);

Section_zx = nan( N_z, N_x );
Section_zx( N_start, : ) = Section_0( floor(N_y/2), : );

Section_zy = nan( N_z, N_y );
Section_zy( N_start, : ) = Section_0( :, floor(N_x/2) );

scrsz = get(groot,'ScreenSize');
fig = figure('Position',[50 50 scrsz(3)-100 scrsz(4)-200]);

pos_XZ = 0;
pos_YZ = 0;

[Grid_zx_x, Grid_zx_z] = meshgrid(Axis_x, Axis_z);
[Grid_zy_y, Grid_zy_z] = meshgrid(Axis_y, Axis_z);

tile = tiledlayout(fig, 1, 2);
axAmp = nexttile;
view(axAmp,3);
grid on;
box on;
hold on;
ampXZ = surf(Grid_zx_x, pos_XZ * ones(N_z,N_x), Grid_zx_z, log10(sqrt(abs( Section_zx ))), 'EdgeColor','none', 'FaceAlpha', 0.2);
ampYZ = surf(pos_YZ * ones(N_z,N_y), Grid_zy_y, Grid_zy_z, log10(sqrt(abs( Section_zy ))), 'EdgeColor','none', 'FaceAlpha', 0.2);
ampXY = surf(Grid_x, Grid_y, Axis_z(N_start)*ones(N_y, N_x), log10(sqrt(abs( Section_yx ))), "EdgeColor", 'none', 'FaceAlpha', 0.5);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
clim( [-1 1] );
colormap(axAmp, hot);
cbAmp = colorbar(axAmp);
cbAmp.Label.String = 'dB';

axPhs = nexttile;
view(axPhs,3);
grid on;
box on;
hold on;
phsXZ = surf(Grid_zx_x, pos_XZ * ones(N_z,N_x), Grid_zx_z, angle( Section_zx ), 'EdgeColor','none', 'FaceAlpha', 0.2);
phsYZ = surf(pos_YZ * ones(N_z,N_y), Grid_zy_y, Grid_zy_z, angle( Section_zy ), 'EdgeColor','none', 'FaceAlpha', 0.2);
phsXY = surf(Grid_x, Grid_y, Axis_z(N_start)*ones(N_y, N_x), angle( Section_yx ), "EdgeColor", 'none', 'FaceAlpha', 0.5);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
clim( [-pi pi] );
colormap(axPhs, jet);
cbPhs = colorbar(axPhs);
cbPhs.Label.String = 'rad';

drawnow;

frame = getframe(fig);
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,'test.gif','gif','LoopCount',Inf,'DelayTime',0.05);

%% FFT_BPM Caculation
temp = Section_0;

covKernal = exp( -1i * ( K_x.^2 + K_y.^2 ) ./ (k + sqrt((k)^2 - K_x.^2 - K_y.^2)) * d_z );
covKernal_back = ifftshift( conj(covKernal) );  %% conjuate for backpropagation
covKernal_forw = ifftshift( covKernal );

%% Free space backpropagation to the entrance plane
for m = N_start:-1:1
    disp( [ num2str(m) ' / ' num2str(N_z) ] );

    %% Propagation
    temp = fft2(temp) .* covKernal_back;
    temp( K_r > k ) = 0;

    %% Phase modulation
    temp = ifft2(temp) .* exp( 1i * k * 0 * d_z);

    Section_yx = temp;
    Section_zx( m, : ) = Section_yx( floor(N_y/2), : );
    Section_zy( m, : ) = Section_yx( :, floor(N_x/2) );

    ampXZ.CData = log10( sqrt(abs( Section_zx )) );
    ampYZ.CData = log10( sqrt(abs( Section_zy )) );
    ampXY.CData = log10( sqrt(abs( Section_yx )) );
    ampXY.ZData = Axis_z(m)*ones(N_y, N_x);

    phsXZ.CData = angle( Section_zx );
    phsYZ.CData = angle( Section_zy );
    phsXY.CData = angle( Section_yx );
    phsXY.ZData = Axis_z(m)*ones(N_y, N_x);

    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,'test.gif','gif','WriteMode','append','DelayTime',0.05);
end

%% Forward propagation in inhomogeneous medium
for m = 2:N_z
    disp( [ num2str(m) ' / ' num2str(N_z) ] );

    %% Propagation
    temp = fft2(temp) .* covKernal_forw;
    temp( K_r > k ) = 0;

    %% Phase modulation
    temp = ifft2(temp) .* exp( 1i * k * n_delta( :, :, m ) * d_z);

    Section_yx = temp;
    Section_zx( m, : ) = Section_yx( floor(N_y/2), : );
    Section_zy( m, : ) = Section_yx( :, floor(N_x/2) );

    ampXZ.CData = log10( sqrt(abs( Section_zx )) );
    ampYZ.CData = log10( sqrt(abs( Section_zy )) );
    ampXY.CData = log10( sqrt(abs( Section_yx )) );
    ampXY.ZData = Axis_z(m)*ones(N_y, N_x);

    phsXZ.CData = angle( Section_zx );
    phsYZ.CData = angle( Section_zy );
    phsXY.CData = angle( Section_yx );
    phsXY.ZData = Axis_z(m)*ones(N_y, N_x);

    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,'test.gif','gif','WriteMode','append','DelayTime',0.05);
end

%% Free space backpropagation to the focal plane
for m = N_z:-1:N_start
    disp( [ num2str(m) ' / ' num2str(N_z) ] );

    %% Propagation
    temp = fft2(temp) .* covKernal_back;
    temp( K_r > k ) = 0;

    %% Phase modulation
    temp = ifft2(temp) .* exp( 1i * k * 0 * d_z);

    Section_yx = temp;
    Section_zx( m, : ) = Section_yx( floor(N_y/2), : );
    Section_zy( m, : ) = Section_yx( :, floor(N_x/2) );

    ampXZ.CData = log10( sqrt(abs( Section_zx )) );
    ampYZ.CData = log10( sqrt(abs( Section_zy )) );
    ampXY.CData = log10( sqrt(abs( Section_yx )) );
    ampXY.ZData = Axis_z(m)*ones(N_y, N_x);

    phsXZ.CData = angle( Section_zx );
    phsYZ.CData = angle( Section_zy );
    phsXY.CData = angle( Section_yx );
    phsXY.ZData = Axis_z(m)*ones(N_y, N_x);

    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,'test.gif','gif','WriteMode','append','DelayTime',0.05);
end