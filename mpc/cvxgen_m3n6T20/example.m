clear

if 0
  m = 3;
  n = 6;
  T = 20;
  dt = .1;
  Q = randn(n,n); Q = Q'*Q;
  R = randn(m,m); R = R'*R;
  A = eye + dt*randn(n,n);
  B = dt*randn(n,m);
  umax = 20;
  x0 = randn(n,1)/dt;
else
  load mpc_data
end

params.A = A;
params.B = B;
params.Q = Q;
params.R = R;
params.u_max = umax;
params.x_0 = x0;

tic
[vars, status] = csolve(params);
toc
