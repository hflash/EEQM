OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[3];
u3(1.91063,0,0) q[0];
u2(-pi/2,-pi) q[1];
cx q[0],q[1];
u2(pi/4,-pi) q[1];
cx q[0],q[1];
u1(pi/2) q[0];
u2(-pi/2,pi/4) q[1];
u2(0,pi) q[2];
cx q[1],q[2];
u1(-pi/4) q[2];
cx q[0],q[2];
u1(pi/4) q[2];
cx q[1],q[2];
u1(pi/4) q[1];
u1(-pi/4) q[2];
cx q[0],q[2];
cx q[0],q[1];
u1(pi/4) q[0];
u1(-pi/4) q[1];
cx q[0],q[1];
u3(pi,0,-pi) q[0];
u3(pi,0,-pi) q[1];
cx q[0],q[1];
u2(0,-3*pi/4) q[2];

