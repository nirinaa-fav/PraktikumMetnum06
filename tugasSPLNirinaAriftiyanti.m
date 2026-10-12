1;  % tugasSPL Nirina Ariftiyanti

% ---------- Eliminasi Gauss ----------
function [x, waktu] = Gauss(A, b)
  tic;
  n = length(b);
  Ab = [A b(:)];
  for k = 1:n-1
    for i = k+1:n
      m = Ab(i,k) / Ab(k,k);
      Ab(i,:) = Ab(i,:) - m * Ab(k,:);
    end
  end
  x = zeros(n,1);
  for i = n:-1:1
    x(i) = (Ab(i,n+1) - Ab(i,i+1:n) * x(i+1:n)) / Ab(i,i);
  end
  waktu = toc;
end

% ---------- Eliminasi Gauss-Jordan ----------
function [x, waktu] = GaussJordan(A, b)
  tic;
  n = length(b);
  Ab = [A b(:)];
  for k = 1:n
    Ab(k,:) = Ab(k,:) / Ab(k,k);
    for i = 1:n
      if i ~= k
        Ab(i,:) = Ab(i,:) - Ab(i,k) * Ab(k,:);
      end
    end
  end
  x = Ab(:,n+1);
  waktu = toc;
end

% ---------- Dekomposisi LU (Doolittle) ----------
function [x, waktu] = DekomposisiLU(A, b)
  tic;
  n = length(b);
  L = eye(n);
  U = A;
  for k = 1:n-1
    for i = k+1:n
      L(i,k) = U(i,k) / U(k,k);
      U(i,:) = U(i,:) - L(i,k) * U(k,:);
    end
  end
  y = zeros(n,1);                 % Ly = b
  for i = 1:n
    y(i) = b(i) - L(i,1:i-1) * y(1:i-1);
  end
  x = zeros(n,1);                 % Ux = y
  for i = n:-1:1
    x(i) = (y(i) - U(i,i+1:n) * x(i+1:n)) / U(i,i);
  end
  waktu = toc;
end

% ---------- Iterasi Jacobi ----------
function [x, iter, waktu] = Jacobi(A, b, tol, max_iter)
  tic;
  [n, m] = size(A);
  for i = 1:n
    xlama(i) = b(i) / A(i,i);
  end
  xlama = xlama';
  C = -A;
  for i = 1:n
    C(i,i) = 0.0;
    C(i,:) = C(i,:) / A(i,i);
    d(i,1) = b(i) / A(i,i);
  end
  i = 1;
  while (i <= max_iter)
    xbaru = C * xlama + d;
    if (abs(xbaru - xlama) <= tol)
      x = xbaru;
      iter = i;
      waktu = toc;
      disp('Jacobi method konverge');
      return;
    else
      xlama = xbaru;
    end
    i = i + 1;
  end
  x = xlama;
  iter = max_iter;
  waktu = toc;
  disp('Jacobi method not konverge');
end

% ---------- Iterasi Gauss-Seidel ----------
function [x, iter, waktu] = IterasiGaussSeidel(A, b, tol, max_iter)
  tic;
  [n, m] = size(A);
  for i = 1:n
    if (i == 1)
      x(i) = b(i) / A(i,i);
    else
      x(i) = 0;
    end
  end
  x = x';
  C = -A;
  for i = 1:n
    C(i,i) = 0;
    C(i,:) = C(i,:) / A(i,i);
    d(i,1) = b(i) / A(i,i);
  end
  i = 1;
  while (i <= max_iter)
    xlama = x;
    for j = 1:n
      x(j) = C(j,:) * x + d(j);
    end
    if (abs(xlama - x) <= tol)
      iter = i;
      waktu = toc;
      disp('Gauss Seidel method konverge');
      return;
    end
    i = i + 1;
  end
  iter = max_iter;
  waktu = toc;
  disp('Gauss Seidel method not konverge');
end

% Urutan variabel: p q r s t
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];
b = [-1; 2; 6; 2; -1];

max_iter = 20;      % maksimum iterasi
tol      = 0.001;   % toleransi error

% ---- Cek diagonal dominan ----
disp('=== Cek diagonal dominan (|aii| >= sum |aij|, j~=i) ===');
for i = 1:5
  dg  = abs(A(i,i));
  sm  = sum(abs(A(i,:))) - dg;
  printf('Baris %d: %g vs %g\n', i, dg, sm);
end

% ---- SOAL 1 ----
disp(' ');
disp('=== SOAL 1: Himpunan Penyelesaian ===');
[x_gauss,  t_gauss]         = Gauss(A, b);
[x_gj,     t_gj]            = GaussJordan(A, b);
[x_lu,     t_lu]            = DekomposisiLU(A, b);
[x_jacobi, it_j,  t_jacobi] = Jacobi(A, b, tol, max_iter);
[x_gs,     it_gs, t_gs]     = IterasiGaussSeidel(A, b, tol, max_iter);

nama = {'p','q','r','s','t'};
printf('\n%-4s %12s %12s %12s %12s %12s\n', ...
       'Var','Gauss','Gauss-Jordan','LU','Jacobi','Gauss-Seidel');
for i = 1:5
  printf('%-4s %12.6f %12.6f %12.6f %12.6f %12.6f\n', nama{i}, ...
         x_gauss(i), x_gj(i), x_lu(i), x_jacobi(i), x_gs(i));
end
printf('\nIterasi Jacobi       : %d\n', it_j);
printf('Iterasi Gauss-Seidel : %d\n', it_gs);

% ---- SOAL 2 ----
% Nilai eksak diambil dari metode analitik/langsung (Gauss-Jordan)
x_eksak = x_gj;
galat_jacobi = abs(x_eksak - x_jacobi(:));
galat_gs     = abs(x_eksak - x_gs(:));

disp(' ');
disp('=== SOAL 2: Galat = abs(nilai eksak - nilai hampiran) ===');
printf('%-4s %12s %16s %16s\n', 'Var','Eksak','Galat Jacobi','Galat G-Seidel');
for i = 1:5
  printf('%-4s %12.6f %16.8f %16.8f\n', nama{i}, x_eksak(i), ...
         galat_jacobi(i), galat_gs(i));
end
printf('Galat maks Jacobi       : %.8f\n', max(galat_jacobi));
printf('Galat maks Gauss-Seidel : %.8f\n', max(galat_gs));

% ---- SOAL 3 ----
disp(' ');
disp('=== SOAL 3: Running Time (detik) ===');
printf('Gauss        : %.8f\n', t_gauss);
printf('Gauss-Jordan : %.8f\n', t_gj);
printf('Dekomp. LU   : %.8f\n', t_lu);
printf('Jacobi       : %.8f\n', t_jacobi);
printf('Gauss-Seidel : %.8f\n', t_gs);


