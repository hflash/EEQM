OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg m2[1];
creg m0[1];
creg m1[1];
u2(0,pi) q[0];
u2(0,pi) q[1];
u2(0,pi) q[2];
cx q[0],q[2];
u1(5.65442695349013) q[2];
cx q[0],q[2];
cx q[0],q[1];
cx q[1],q[2];
u1(-11.3088853229068) q[2];
cx q[1],q[2];
cx q[0],q[1];
u3(1.7132487,-pi/2,pi/2) q[0];
u3(1.7132487,-pi/2,-2.8261453) q[1];
u3(1.7132487,-pi/2,pi/2) q[2];
