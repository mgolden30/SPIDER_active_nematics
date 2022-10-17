addpath("../bin")

fprintf('Starting n cross product n regression...\n');

num_lib         = 8;
num_windows     = num_lib*10;
nw = num_windows; %some shorthand
nl = num_lib;

labels = cell(nl, 1);

res_aves = [];
res_stds = [];
c_aves = [];
c_stds = [];

corners = zeros(3, num_windows);
size_vec = round( 2./dx_vec ); %integrate over length&time scales
size_vec(3) = 50; %manually change for rice

Ly = size(U,1);
Lx = size(U,2);
Lt = size(U,3);
for i=1:num_windows
    boundary = 10;
    corners(:,i) = [ randi(Ly-size_vec(1)-2*boundary) randi(Lx-size_vec(2)-2*boundary) randi(Lt-size_vec(3)-2*boundary) ] + [1 1 1]*boundary;
end

%The polynomial weights uses these powers: x^m(1-x)^n
m = 4*[1 1 1];
n = m;

G = zeros(num_windows, num_lib);

un = U.*n1 + V.*n2;
B1 = n1.*n1_x + n2.*n1_y;
B2 = n1.*n2_x + n2.*n2_y;
S  = n1_x + n2_y;

[B1_x, B1_y, ~] = gradient3d(B1, dx_vec);
[B2_x, B2_y, ~] = gradient3d(B2, dx_vec);



vorticity = V_x - U_y;


G_full = [];
lengths = size_vec.*dx_vec;
%weights = { 0*x+1, cos( pi*x/lengths(2) ), cos( pi*y/lengths(1)), cos( pi*t/lengths(3)) };
mask0 = mask;

for ii = 1:1%numel(weights)
a = 1;

clear functions            %We will recompute derivatives of the mask
mask = mask0;%.*weights{ii}; %Multiply your weight by anything you want


G(:, a) = gpu_integrate_data_cross_n( n1_t, n2_t, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\partial_t n_i";
a = a+1;

G(:,a) = gpu_integrate_data_cross_n( U.*(n1_x + n2_y), V.*(n1_x + n2_y), n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i S";
a=a+1;

G(:,a) = gpu_integrate_data_cross_n( n1_xx + n2_xy, n1_xy + n2_yy, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\nabla_i S";
a=a+1;

%Using Ryskin's definition of \Omega_{ij}
G(:,a) = gpu_integrate_data_cross_n( -vorticity.*n2, vorticity.*n1, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "\Omega_{ij} n_j";
a = a+1;

G(:, a) = gpu_integrate_data_cross_n( A_11.*n1 + A_12.*n2, A_12.*n1 - A_11.*n2, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "A'_{ij} n_j";
a = a+1;

G(:, a) = gpu_integrate_data_cross_n( U.*n1_x + V.*n1_y, U.*n2_x + V.*n2_y, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_j \nabla_j n_i";
a = a+1;


%{
G(:, a) = gpu_integrate_data_cross_n( un.*B1, un.*B2, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j u_j B_i";
a = a+1;
%}

G(:, a) = gpu_integrate_data_cross_n( S.*B1, S.*B2, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S B_i";
a = a+1;

G(:, a) = gpu_integrate_data_cross_n( n1.*B1_x + n2.*B1_y, n1.*B2_x + n2.*B2_y, n1, n2, dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j \nabla_j B_i";
a = a+1;

G_full = [G_full; G];
end %end weight iteration

G = G_full;
mask = mask0;
fprintf("done integrating.\n")





function vals = gpu_integrate_data_cross_n( u, v, n1, n2, dx_vec, corners, size_vec, m, n, mask )
  to_be_integrated = n1.*v - n2.*u;
  vals = gpu_integrate_data( to_be_integrated, [], dx_vec, corners, size_vec, m, n, mask );
  
  vals = gather(vals);
end

