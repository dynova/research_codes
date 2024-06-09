format bank

x = readNPY('inferences/Dallas4-n5-parLog.npy');
y = readNPY('inferences/Dallas4-n5-logLLog.npy');
[aa,bb] = max(y);
mle = x(bb,:);
z = [mle(1) mle(21) mle(1)+mle(2) mle(12) mle(7) mle(1)+mle(2)+mle(3) mle(13) mle(8) mle(1)+mle(2)+mle(3)+mle(4) mle(14) mle(9) sum(mle(1:5)) mle(15) mle(10) sum(mle(1:6)) mle(16) mle(11) mle(17) mle(19) mle(17)+mle(18) mle(20) mle(22) mle(23)];
disp(z');

x = readNPY('inferences/Houston2-n4-parLog.npy');
y = readNPY('inferences/Houston2-n4-logLLog.npy');
[aa,bb] = max(y);
mle = x(bb,:);
z = [mle(1) mle(18) mle(1)+mle(2) mle(10) mle(6) mle(1)+mle(2)+mle(3) mle(11) mle(7) mle(1)+mle(2)+mle(3)+mle(4) mle(12) mle(8) sum(mle(1:5)) mle(13) mle(9) mle(14) mle(16) mle(14)+mle(15) mle(17) mle(19) mle(20)];
disp(z');

x = readNPY('inferences/NYC0-n4-parLog.npy');
y = readNPY('inferences/NYC0-n4-logLLog.npy');
[aa,bb] = max(y);
mle = x(bb,:);
z = [mle(1) mle(18) mle(1)+mle(2) mle(10) mle(6) mle(1)+mle(2)+mle(3) mle(11) mle(7) mle(1)+mle(2)+mle(3)+mle(4) mle(12) mle(8) sum(mle(1:5)) mle(13) mle(9) mle(14) mle(16) mle(14)+mle(15) mle(17) mle(19) mle(20)];
disp(z');

x = readNPY('inferences/Phoenix29-n5-parLog.npy');
y = readNPY('inferences/Phoenix29-n5-logLLog.npy');
[aa,bb] = max(y);
mle = x(bb,:);
z = [mle(1) mle(21) mle(1)+mle(2) mle(12) mle(7) mle(1)+mle(2)+mle(3) mle(13) mle(8) mle(1)+mle(2)+mle(3)+mle(4) mle(14) mle(9) sum(mle(1:5)) mle(15) mle(10) sum(mle(1:6)) mle(16) mle(11) mle(17) mle(19) mle(17)+mle(18) mle(20) mle(22) mle(23)];
disp(z');