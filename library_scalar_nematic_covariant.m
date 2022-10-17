res_aves = [];
res_stds = [];
fprintf('Starting scalar regression for libary A2...\n');

num_lib         = 17;
num_windows     = num_lib*40;
nw = num_windows; %some shorthand
nl = num_lib;

labels = cell(nl,1);
corners = zeros(3, num_windows);
size_vec= round( 2./dx_vec );

L = 10;
for i=1:num_windows
    corners(:,i) = [ randi(Ly-size_vec(1)-2*L) randi(Lx-size_vec(2)-2*L) randi(Lt-size_vec(3)-2*L) ] + [1 1 1]*L;
end

%The polynomial weights uses these powers: x^m(1-x)^n
m = 4*[1 1 1];
n = m;

G = zeros(2*num_windows, num_lib);

%u dot n shows up a lot.
un = U.*n1 + V.*n2;
S  = n1_x + n2_y;
B1 = -(n1.*n1_x + n2.*n1_y);
B2 = -(n1.*n2_x + n2.*n2_y);

[U_yy, U_xx, U_tt, U_xy, U_xt, U_yt] = secondgradient3d( U, dx_vec );
[V_yy, V_xx, V_tt, V_xy, V_xt, V_yt] = secondgradient3d( V, dx_vec );
[B1_x, B1_y,~] = gradient3d(B1, dx_vec);
[B2_x, B2_y,~] = gradient3d(B2, dx_vec);



G_full = [];
lengths = size_vec.*dx_vec;
weights = { 0*x+1, cos( x/lengths(2) ), cos(y/lengths(1)), cos(t/lengths(3)) };
mask0 = mask;

for ii = 1:numel(weights)
a = 1;

clear functions            %We will recompute derivatives of the mask
mask = mask0.*weights{ii}; %Multiply your weight by anything you want




G(:,a) = gpu_integrate_data_A2( n1_x + n2_y, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S";
a = a+1;

G(:,a) = gpu_integrate_data_A2( un, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(u_i n_i)";
a = a+1;

G(:,a) = gpu_integrate_data_A2( U.*n1_t + V.*n2_t, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i \partial_t n_i";
a = a+1;

G(:,a) = gpu_integrate_data_A2( U_t.*n1 + V_t.*n2, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i \partial_t u_i";
a = a+1;

G(:,a) = gpu_integrate_data_A2( un.*S.^2, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i u_i S^2";
a = a+1;

G(:,a) = gpu_integrate_data_A2( S.*(U.*B1 + V.*B2), n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S u_i B_i";
a = a+1;

G(:,a) = gpu_integrate_data_A2( S.*(U_x + V_y), n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S \nabla_i u_i";
a = a+1;

G(:,a) = gpu_integrate_data_A2( S.*Ann, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "S n_i n_j A'_{ij}";
a = a+1;

G(:,a) = gpu_integrate_data_A2( un.*( n1.*(n1_xx + n2_xy) + n2.*(n1_xy + n2_yy) ), n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(n_i u_i) n_j \nabla_j S";
a = a+1;

G(:,a) = gpu_integrate_data_A2( U.*(n1_xx + n2_xy) + V.*(n1_xy + n2_yy), n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i \nabla_i S";
a = a+1;

scalar = n1.*n1.*n1.*U_xx + n1.*n2.*n2.*U_yy + 2*n1.*n1.*n2.*U_xy + ...
         n2.*n1.*n1.*V_xx + n2.*n2.*n2.*V_yy + 2*n2.*n1.*n2.*V_xy;
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i n_j n_k \nabla_i \nabla_j u_k";
a = a+1;

scalar = n1.*(U_xx + U_yy) + n2.*(V_xx + V_yy);
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i \nabla^2 u_i";
a = a+1;

scalar = n1.*(U_xx + V_xy) + n2.*(U_xy + V_yy);
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i \nabla_i \nabla_j u_j";
a = a+1;

scalar = un.*(B1.^2 + B2.^2);
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i n_i B^2";
a = a+1;

scalar = n1.*B1.*U_x + n1.*B2.*U_y + n2.*B1.*V_x + n2.*B2.*V_y;
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_i B_j \nabla_j u_i";
a = a+1;

scalar = n1.*B1.*U_x + n2.*B1.*U_y + n1.*B2.*V_x + n2.*B2.*V_y;
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "n_j B_i \nabla_j u_i";
a = a+1;

%{
scalar = un.*(n1.*(n1_xx + n2_xy) + n2.*(n2_yy + n1_xy) + S.^2);
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "(u_i n_i) \nabla_j B_j";
a = a+1;
%}

scalar = U.*(n1.*B1_x + n2.*B1_y) + V.*(n1.*B2_x + n2.*B2_y);
G(:,a) = gpu_integrate_data_A2( scalar, n1, n2, [], dx_vec, corners, size_vec, m, n, mask );
labels{a} = "u_i n_j \nabla_j B_i";
a = a+1;


G_full = [G_full; G];
end %end weight iteration

G = G_full;
mask = mask0;
fprintf("done integrating.\n")




function g = gpu_integrate_data_A2( data, n1, n2, derivs, dx_vec, corners, size_vec, m, n, mask )
  l = round(numel(corners)/3);
  g = zeros(2*l,1);
  
  g(1:2:end) = gpu_integrate_data( data.*n1, derivs, dx_vec, corners, size_vec, m, n, mask );
  g(2:2:end) = gpu_integrate_data( data.*n2, derivs, dx_vec, corners, size_vec, m, n, mask );
  
  g = gather(g);
end