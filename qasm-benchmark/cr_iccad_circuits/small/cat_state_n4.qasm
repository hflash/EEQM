OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[4];
u2(0,pi) q[0];
cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[3];

