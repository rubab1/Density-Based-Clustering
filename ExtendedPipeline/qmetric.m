function distance = qmetric(tiles1, tiles2)
% QMETRIC 
%
% QMETRIC 
%
%   usage: distance = qmetric(tiles1, tiles2);
%

t1 = tiles1(:, 1);
f1 = tiles1(:, 2);
q1 = tiles1(:, 3);

t2 = tiles2(:, 1);
f2 = tiles2(:, 2);
q2 = tiles2(:, 3);

t = (t1 + t2) / 2;
f = sqrt(f1 .* f2);
q = sqrt(q1 .* q2);

dt = t2 - t1;
df = f2 - f1;
dq = q2 - q1;

gtt = 4 * pi^2 * f.^2 ./ q.^2;
gff = (2 + q.^2) ./ (4 .* f.^2);
gqq = 1 ./ (2 * q.^2);
gfq = - 1 ./ (2 * f .* q);
           
% distance = sqrt(dt.^2 + df.^2);

distance = sqrt(gtt .* dt.^2 + ...
                gff .* df.^2 + ...
                gqq .* dq.^2 + ...
                2 * gfq .* df .* dq);

% df1 = 2 * sqrt(pi) * f1 ./ q1;
% df2 = 2 * sqrt(pi) * f2 ./ q2;
% dt1 = 1 ./ df1;
% dt2 = 1 ./ df2;

% df1 = df1 * 0.5;
% df2 = df2 * 0.5;
% dt1 = dt1 * 0.5;
% dt2 = dt2 * 0.5;

% dt = max(0, abs(t2 - t1) - (dt1 + dt2) / 2);
% df = max(0, abs(f2 - f1) - (df1 + df2) / 2);

% distance = sqrt(dt.^2 + df.^2);
