OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[4];
u2(0,0) q[0];
cx q[0],q[1];
u2(-pi,-pi) q[0];
u2(0,-pi) q[1];
cx q[0],q[1];
u2(0,pi) q[0];
u2(0,0) q[2];
cx q[2],q[3];
u2(-pi,-pi) q[2];
u2(0,-pi) q[3];
cx q[2],q[3];
u2(0,pi) q[2];

