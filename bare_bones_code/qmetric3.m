function distance = qmetric3(tiles1, tiles2)
% QMETRIC 
%
% QMETRIC 
%
%   usage: distance = qmetric(tiles1, tiles2);
%

t1 = tiles1(:, 1);
f1 = tiles1(:, 2);
q1 = tiles1(:, 3);
z1 = tiles1(:, 4);

t2 = tiles2(:, 1);
f2 = tiles2(:, 2);
q2 = tiles2(:, 3);
z2 = tiles2(:, 4);

dq = q1 - q2;
dz = z1 - z2;

df1 = 2 * sqrt(pi) * f1 ./ q1;
df2 = 2 * sqrt(pi) * f2 ./ q2;
dt1 = 1 ./ df1;
dt2 = 1 ./ df2;

df1 = df1 * 0.5;
df2 = df2 * 0.5;
dt1 = dt1 * 0.5;
dt2 = dt2 * 0.5;

dt = max(0, abs(t2 - t1) - (dt1 + dt2) / 2);
df = max(0, abs(f2 - f1) - (df1 + df2) / 2);

% distance = sqrt(dt.^2 + df.^2); % 0.5, 1

% distance = sqrt(100*dt + 100*df + sqrt(dq) + 8*dz); % 7.6, 2a

% distance = sqrt(100*dt + 100*df + 2*sqrt(dq) + 8*dz); % 8.05 2b

distance = sqrt(100*dt + 100*df + 2* sqrt(dq) + 10*dz); % 8.7 2c
