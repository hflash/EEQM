OPENQASM 2.0;
include "qelib1.inc";
qreg q[80];
creg c[80];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
h q[47];
h q[48];
h q[49];
h q[50];
h q[51];
h q[52];
h q[53];
h q[54];
h q[55];
h q[56];
h q[57];
h q[58];
h q[59];
h q[60];
h q[61];
h q[62];
h q[63];
h q[64];
h q[65];
h q[66];
h q[67];
h q[68];
h q[69];
h q[70];
h q[71];
h q[72];
h q[73];
h q[74];
h q[75];
h q[76];
h q[77];
h q[78];
h q[79];
cx q[10], q[41];
cx q[28], q[3];
cx q[31], q[1];
cx q[22], q[75];
cx q[17], q[47];
cx q[68], q[36];
cx q[25], q[55];
cx q[14], q[12];
cx q[69], q[16];
cx q[76], q[65];
cx q[79], q[0];
cx q[4], q[49];
cx q[35], q[63];
cx q[58], q[30];
cx q[54], q[8];
cx q[33], q[42];
cx q[78], q[52];
cx q[18], q[44];
cx q[64], q[62];
cx q[20], q[45];
rx(0.2) q[2];
rx(0.2) q[13];
rx(0.2) q[15];
rx(0.2) q[43];
rx(0.2) q[57];
rx(0.2) q[59];
rx(0.2) q[67];
rz(0.4) q[41];
rz(0.4) q[3];
rz(0.4) q[1];
rz(0.4) q[75];
rz(0.4) q[47];
rz(0.4) q[36];
rz(0.4) q[55];
rz(0.4) q[12];
rz(0.4) q[16];
rz(0.4) q[65];
rz(0.4) q[0];
rz(0.4) q[49];
rz(0.4) q[63];
rz(0.4) q[30];
rz(0.4) q[8];
rz(0.4) q[42];
rz(0.4) q[52];
rz(0.4) q[44];
rz(0.4) q[62];
rz(0.4) q[45];
rx(1) q[2];
rx(1) q[13];
rx(1) q[15];
rx(1) q[43];
rx(1) q[57];
rx(1) q[59];
rx(1) q[67];
cx q[10], q[41];
cx q[28], q[3];
cx q[31], q[1];
cx q[22], q[75];
cx q[17], q[47];
cx q[68], q[36];
cx q[25], q[55];
cx q[14], q[12];
cx q[69], q[16];
cx q[76], q[65];
cx q[79], q[0];
cx q[4], q[49];
cx q[35], q[63];
cx q[58], q[30];
cx q[54], q[8];
cx q[33], q[42];
cx q[78], q[52];
cx q[18], q[44];
cx q[64], q[62];
cx q[20], q[45];
rx(0.2) q[2];
rx(0.2) q[13];
rx(0.2) q[15];
rx(0.2) q[43];
rx(0.2) q[57];
rx(0.2) q[59];
rx(0.2) q[67];
cx q[28], q[37];
cx q[56], q[22];
cx q[29], q[55];
cx q[38], q[69];
cx q[11], q[47];
cx q[14], q[46];
cx q[6], q[10];
cx q[50], q[65];
cx q[27], q[25];
cx q[58], q[4];
cx q[35], q[19];
cx q[16], q[52];
cx q[26], q[42];
cx q[31], q[18];
cx q[74], q[36];
cx q[21], q[45];
cx q[73], q[44];
cx q[51], q[64];
rx(0.2) q[0];
rx(0.2) q[1];
rx(0.2) q[8];
rx(0.2) q[17];
rx(0.2) q[20];
rx(0.2) q[33];
rx(0.2) q[49];
rx(0.2) q[54];
rx(0.2) q[62];
rx(0.2) q[78];
rx(0.2) q[79];
rz(0.4) q[37];
rz(0.4) q[22];
rz(0.4) q[55];
rz(0.4) q[69];
rz(0.4) q[47];
rz(0.4) q[46];
rz(0.4) q[10];
rz(0.4) q[65];
rz(0.4) q[25];
rz(0.4) q[4];
rz(0.4) q[19];
rz(0.4) q[52];
rz(0.4) q[42];
rz(0.4) q[18];
rz(0.4) q[36];
rz(0.4) q[45];
rz(0.4) q[44];
rz(0.4) q[64];
cx q[79], q[0];
cx q[54], q[8];
cx q[28], q[37];
cx q[56], q[22];
cx q[29], q[55];
cx q[38], q[69];
cx q[11], q[47];
cx q[14], q[46];
cx q[6], q[10];
cx q[50], q[65];
cx q[27], q[25];
cx q[58], q[4];
cx q[35], q[19];
cx q[16], q[52];
cx q[26], q[42];
cx q[31], q[18];
cx q[74], q[36];
cx q[21], q[45];
cx q[73], q[44];
cx q[51], q[64];
rz(1.2) q[0];
rz(1.2) q[8];
cx q[48], q[56];
cx q[6], q[24];
cx q[11], q[40];
cx q[39], q[69];
cx q[5], q[10];
cx q[29], q[60];
cx q[47], q[32];
cx q[50], q[28];
cx q[4], q[61];
cx q[70], q[22];
cx q[23], q[35];
cx q[46], q[68];
cx q[52], q[76];
cx q[58], q[21];
rx(0.2) q[14];
rx(0.2) q[26];
rx(0.2) q[31];
rx(0.2) q[36];
rx(0.2) q[44];
rx(0.2) q[51];
rx(0.2) q[55];
rx(0.2) q[64];
rx(0.2) q[65];
rx(0.2) q[73];
cx q[79], q[0];
cx q[54], q[8];
rz(0.4) q[56];
rz(0.4) q[24];
rz(0.4) q[40];
rz(0.4) q[69];
rz(0.4) q[10];
rz(0.4) q[60];
rz(0.4) q[32];
rz(0.4) q[28];
rz(0.4) q[61];
rz(0.4) q[22];
rz(0.4) q[35];
rz(0.4) q[68];
rz(0.4) q[76];
rz(0.4) q[21];
cx q[31], q[1];
cx q[64], q[62];
rx(1) q[0];
rx(1) q[8];
rx(1) q[54];
rx(1) q[79];
cx q[48], q[56];
cx q[6], q[24];
cx q[11], q[40];
cx q[39], q[69];
cx q[5], q[10];
cx q[29], q[60];
cx q[47], q[32];
cx q[50], q[28];
cx q[4], q[61];
cx q[70], q[22];
cx q[23], q[35];
cx q[46], q[68];
cx q[52], q[76];
cx q[58], q[21];
rz(1.2) q[1];
rz(1.2) q[62];
cx q[79], q[0];
cx q[54], q[8];
cx q[72], q[56];
cx q[48], q[3];
cx q[47], q[27];
cx q[71], q[23];
cx q[35], q[12];
cx q[70], q[30];
cx q[6], q[74];
cx q[45], q[68];
cx q[60], q[69];
cx q[39], q[34];
cx q[22], q[76];
rx(0.2) q[5];
rx(0.2) q[24];
rx(0.2) q[28];
rx(0.2) q[29];
rx(0.2) q[32];
rx(0.2) q[40];
rx(0.2) q[50];
rx(0.2) q[52];
rx(0.2) q[58];
cx q[31], q[1];
cx q[64], q[62];
rz(0.4) q[0];
rz(0.4) q[8];
rz(0.4) q[56];
rz(0.4) q[3];
rz(0.4) q[27];
rz(0.4) q[23];
rz(0.4) q[12];
rz(0.4) q[30];
rz(0.4) q[74];
rz(0.4) q[68];
rz(0.4) q[69];
rz(0.4) q[34];
rz(0.4) q[76];
cx q[78], q[52];
cx q[51], q[64];
rx(1) q[1];
rx(1) q[62];
cx q[79], q[0];
cx q[54], q[8];
cx q[72], q[56];
cx q[48], q[3];
cx q[47], q[27];
cx q[71], q[23];
cx q[35], q[12];
cx q[70], q[30];
cx q[6], q[74];
cx q[45], q[68];
cx q[60], q[69];
cx q[39], q[34];
cx q[22], q[76];
rz(1.2) q[52];
rz(1.2) q[64];
rx(0.2) q[0];
rx(0.2) q[8];
rx(0.2) q[54];
rx(0.2) q[79];
cx q[37], q[56];
cx q[41], q[72];
cx q[27], q[63];
cx q[47], q[9];
cx q[71], q[61];
cx q[35], q[21];
cx q[3], q[46];
cx q[19], q[76];
rx(0.2) q[6];
rx(0.2) q[12];
rx(0.2) q[22];
rx(0.2) q[23];
rx(0.2) q[30];
rx(0.2) q[34];
rx(0.2) q[39];
rx(0.2) q[45];
rx(0.2) q[60];
rx(0.2) q[68];
rx(0.2) q[69];
rx(0.2) q[74];
cx q[78], q[52];
cx q[51], q[64];
rz(0.4) q[56];
rz(0.4) q[72];
rz(0.4) q[63];
rz(0.4) q[9];
rz(0.4) q[61];
rz(0.4) q[21];
rz(0.4) q[46];
rz(0.4) q[76];
cx q[68], q[36];
cx q[14], q[12];
cx q[58], q[30];
cx q[20], q[45];
rx(1) q[51];
rx(1) q[64];
rx(1) q[78];
cx q[37], q[56];
cx q[41], q[72];
cx q[27], q[63];
cx q[47], q[9];
cx q[71], q[61];
cx q[35], q[21];
cx q[3], q[46];
cx q[19], q[76];
rz(1.2) q[36];
rz(1.2) q[12];
rz(1.2) q[30];
rz(1.2) q[45];
cx q[64], q[62];
cx q[7], q[72];
cx q[77], q[37];
cx q[56], q[75];
cx q[41], q[11];
cx q[61], q[25];
rx(0.2) q[3];
rx(0.2) q[9];
rx(0.2) q[19];
rx(0.2) q[21];
rx(0.2) q[27];
rx(0.2) q[35];
rx(0.2) q[46];
rx(0.2) q[47];
rx(0.2) q[63];
rx(0.2) q[71];
rx(0.2) q[76];
cx q[68], q[36];
cx q[14], q[12];
cx q[58], q[30];
cx q[20], q[45];
rz(0.4) q[62];
rz(0.4) q[72];
rz(0.4) q[37];
rz(0.4) q[75];
rz(0.4) q[11];
rz(0.4) q[25];
cx q[28], q[3];
cx q[17], q[47];
cx q[76], q[65];
cx q[14], q[46];
cx q[35], q[63];
cx q[74], q[36];
cx q[21], q[45];
rx(1) q[20];
cx q[64], q[62];
cx q[7], q[72];
cx q[77], q[37];
cx q[56], q[75];
cx q[41], q[11];
cx q[61], q[25];
rz(1.2) q[3];
rz(1.2) q[47];
rz(1.2) q[65];
rz(1.2) q[46];
rz(1.2) q[63];
rz(1.2) q[36];
rz(1.2) q[45];
cx q[51], q[64];
rx(0.2) q[62];
cx q[37], q[53];
cx q[41], q[4];
cx q[18], q[72];
cx q[77], q[10];
cx q[61], q[42];
rx(0.2) q[7];
rx(0.2) q[11];
rx(0.2) q[25];
rx(0.2) q[56];
rx(0.2) q[75];
cx q[28], q[3];
cx q[17], q[47];
cx q[76], q[65];
cx q[14], q[46];
cx q[35], q[63];
cx q[74], q[36];
cx q[21], q[45];
rz(0.4) q[64];
rz(0.4) q[53];
rz(0.4) q[4];
rz(0.4) q[72];
rz(0.4) q[10];
rz(0.4) q[42];
cx q[22], q[75];
cx q[25], q[55];
cx q[11], q[47];
cx q[50], q[65];
cx q[35], q[19];
cx q[46], q[68];
rx(1) q[14];
rx(1) q[17];
rx(1) q[36];
cx q[51], q[64];
cx q[37], q[53];
cx q[41], q[4];
cx q[18], q[72];
cx q[77], q[10];
cx q[61], q[42];
rz(1.2) q[75];
rz(1.2) q[55];
rz(1.2) q[47];
rz(1.2) q[65];
rz(1.2) q[19];
rz(1.2) q[68];
rx(0.2) q[51];
rx(0.2) q[64];
cx q[38], q[53];
cx q[66], q[41];
cx q[70], q[37];
rx(0.2) q[4];
rx(0.2) q[10];
rx(0.2) q[18];
rx(0.2) q[42];
rx(0.2) q[61];
rx(0.2) q[72];
rx(0.2) q[77];
cx q[22], q[75];
cx q[25], q[55];
cx q[11], q[47];
cx q[50], q[65];
cx q[35], q[19];
cx q[46], q[68];
rz(0.4) q[53];
rz(0.4) q[41];
rz(0.4) q[37];
cx q[56], q[22];
cx q[29], q[55];
cx q[4], q[49];
cx q[11], q[40];
cx q[27], q[25];
cx q[33], q[42];
cx q[47], q[32];
cx q[18], q[44];
cx q[23], q[35];
cx q[45], q[68];
rx(1) q[65];
cx q[38], q[53];
cx q[66], q[41];
cx q[70], q[37];
rz(1.2) q[22];
rz(1.2) q[55];
rz(1.2) q[49];
rz(1.2) q[40];
rz(1.2) q[25];
rz(1.2) q[42];
rz(1.2) q[32];
rz(1.2) q[44];
rz(1.2) q[35];
rz(1.2) q[68];
cx q[48], q[70];
cx q[16], q[41];
rx(0.2) q[37];
rx(0.2) q[38];
rx(0.2) q[53];
rx(0.2) q[66];
cx q[56], q[22];
cx q[29], q[55];
cx q[4], q[49];
cx q[11], q[40];
cx q[27], q[25];
cx q[33], q[42];
cx q[47], q[32];
cx q[18], q[44];
cx q[23], q[35];
cx q[45], q[68];
rz(0.4) q[70];
rz(0.4) q[41];
cx q[28], q[37];
cx q[29], q[60];
cx q[58], q[4];
cx q[47], q[27];
cx q[71], q[23];
cx q[35], q[12];
cx q[26], q[42];
cx q[31], q[18];
cx q[73], q[44];
rx(1) q[32];
rx(1) q[33];
rx(1) q[40];
rx(1) q[45];
rx(1) q[49];
rx(1) q[55];
rx(1) q[68];
cx q[48], q[70];
cx q[16], q[41];
rz(1.2) q[37];
rz(1.2) q[60];
rz(1.2) q[4];
rz(1.2) q[27];
rz(1.2) q[23];
rz(1.2) q[12];
rz(1.2) q[42];
rz(1.2) q[18];
rz(1.2) q[44];
cx q[68], q[36];
cx q[20], q[45];
rx(0.2) q[16];
rx(0.2) q[41];
rx(0.2) q[48];
rx(0.2) q[70];
cx q[28], q[37];
cx q[29], q[60];
cx q[58], q[4];
cx q[47], q[27];
cx q[71], q[23];
cx q[35], q[12];
cx q[26], q[42];
cx q[31], q[18];
cx q[73], q[44];
rz(0.4) q[36];
rz(0.4) q[45];
cx q[10], q[41];
cx q[48], q[56];
cx q[69], q[16];
cx q[50], q[28];
cx q[4], q[61];
cx q[70], q[22];
cx q[27], q[63];
cx q[47], q[9];
cx q[58], q[21];
rx(1) q[12];
rx(1) q[23];
rx(1) q[26];
rx(1) q[29];
rx(1) q[31];
rx(1) q[44];
rx(1) q[73];
cx q[68], q[36];
cx q[20], q[45];
rz(1.2) q[41];
rz(1.2) q[56];
rz(1.2) q[16];
rz(1.2) q[28];
rz(1.2) q[61];
rz(1.2) q[22];
rz(1.2) q[63];
rz(1.2) q[9];
rz(1.2) q[21];
cx q[31], q[1];
cx q[14], q[12];
rx(0.2) q[20];
cx q[10], q[41];
cx q[48], q[56];
cx q[69], q[16];
cx q[50], q[28];
cx q[4], q[61];
cx q[70], q[22];
cx q[27], q[63];
cx q[47], q[9];
cx q[58], q[21];
rz(0.4) q[1];
rz(0.4) q[12];
cx q[72], q[56];
cx q[38], q[69];
cx q[48], q[3];
cx q[6], q[10];
cx q[16], q[52];
cx q[70], q[30];
cx q[71], q[61];
cx q[35], q[21];
rx(1) q[9];
rx(1) q[27];
rx(1) q[28];
rx(1) q[47];
rx(1) q[50];
rx(1) q[58];
rx(1) q[63];
cx q[31], q[1];
cx q[14], q[12];
rz(1.2) q[56];
rz(1.2) q[69];
rz(1.2) q[3];
rz(1.2) q[10];
rz(1.2) q[52];
rz(1.2) q[30];
rz(1.2) q[61];
rz(1.2) q[21];
cx q[17], q[47];
rx(0.2) q[1];
cx q[72], q[56];
cx q[38], q[69];
cx q[48], q[3];
cx q[6], q[10];
cx q[16], q[52];
cx q[70], q[30];
cx q[71], q[61];
cx q[35], q[21];
rz(0.4) q[47];
cx q[37], q[56];
cx q[6], q[24];
cx q[39], q[69];
cx q[5], q[10];
cx q[41], q[72];
cx q[52], q[76];
cx q[61], q[25];
cx q[3], q[46];
rx(1) q[21];
rx(1) q[30];
rx(1) q[35];
rx(1) q[71];
cx q[17], q[47];
rz(1.2) q[56];
rz(1.2) q[24];
rz(1.2) q[69];
rz(1.2) q[10];
rz(1.2) q[72];
rz(1.2) q[76];
rz(1.2) q[25];
rz(1.2) q[46];
cx q[35], q[63];
cx q[58], q[30];
cx q[21], q[45];
rx(0.2) q[17];
cx q[37], q[56];
cx q[6], q[24];
cx q[39], q[69];
cx q[5], q[10];
cx q[41], q[72];
cx q[52], q[76];
cx q[61], q[25];
cx q[3], q[46];
rz(0.4) q[63];
rz(0.4) q[30];
rz(0.4) q[45];
cx q[7], q[72];
cx q[77], q[37];
cx q[56], q[75];
cx q[41], q[11];
cx q[6], q[74];
cx q[60], q[69];
cx q[39], q[34];
cx q[22], q[76];
cx q[61], q[42];
rx(1) q[3];
rx(1) q[5];
rx(1) q[24];
rx(1) q[25];
rx(1) q[46];
rx(1) q[52];
cx q[35], q[63];
cx q[58], q[30];
cx q[21], q[45];
rz(1.2) q[72];
rz(1.2) q[37];
rz(1.2) q[75];
rz(1.2) q[11];
rz(1.2) q[74];
rz(1.2) q[69];
rz(1.2) q[34];
rz(1.2) q[76];
rz(1.2) q[42];
cx q[28], q[3];
cx q[25], q[55];
cx q[14], q[46];
cx q[78], q[52];
cx q[7], q[72];
cx q[77], q[37];
cx q[56], q[75];
cx q[41], q[11];
cx q[6], q[74];
cx q[60], q[69];
cx q[39], q[34];
cx q[22], q[76];
cx q[61], q[42];
rz(0.4) q[3];
rz(0.4) q[55];
rz(0.4) q[46];
rz(0.4) q[52];
cx q[37], q[53];
cx q[41], q[4];
cx q[18], q[72];
cx q[77], q[10];
cx q[19], q[76];
rx(1) q[6];
rx(1) q[7];
rx(1) q[11];
rx(1) q[22];
rx(1) q[34];
rx(1) q[39];
rx(1) q[42];
rx(1) q[56];
rx(1) q[60];
rx(1) q[61];
rx(1) q[69];
rx(1) q[74];
rx(1) q[75];
cx q[28], q[3];
cx q[25], q[55];
cx q[14], q[46];
cx q[78], q[52];
rz(1.2) q[53];
rz(1.2) q[4];
rz(1.2) q[72];
rz(1.2) q[10];
rz(1.2) q[76];
cx q[22], q[75];
cx q[29], q[55];
cx q[11], q[47];
cx q[27], q[25];
cx q[33], q[42];
cx q[46], q[68];
cx q[74], q[36];
rx(0.2) q[14];
rx(0.2) q[78];
cx q[37], q[53];
cx q[41], q[4];
cx q[18], q[72];
cx q[77], q[10];
cx q[19], q[76];
rz(0.4) q[75];
rz(0.4) q[55];
rz(0.4) q[47];
rz(0.4) q[25];
rz(0.4) q[42];
rz(0.4) q[68];
rz(0.4) q[36];
cx q[38], q[53];
cx q[66], q[41];
cx q[70], q[37];
rx(1) q[4];
rx(1) q[10];
rx(1) q[18];
rx(1) q[19];
rx(1) q[72];
rx(1) q[76];
rx(1) q[77];
cx q[22], q[75];
cx q[29], q[55];
cx q[11], q[47];
cx q[27], q[25];
cx q[33], q[42];
cx q[46], q[68];
cx q[74], q[36];
rz(1.2) q[53];
rz(1.2) q[41];
rz(1.2) q[37];
cx q[56], q[22];
cx q[76], q[65];
cx q[4], q[49];
cx q[11], q[40];
cx q[29], q[60];
cx q[47], q[32];
cx q[35], q[19];
cx q[18], q[44];
cx q[26], q[42];
cx q[45], q[68];
rx(0.2) q[33];
rx(0.2) q[36];
rx(0.2) q[55];
cx q[38], q[53];
cx q[66], q[41];
cx q[70], q[37];
rz(0.4) q[22];
rz(0.4) q[65];
rz(0.4) q[49];
rz(0.4) q[40];
rz(0.4) q[60];
rz(0.4) q[32];
rz(0.4) q[19];
rz(0.4) q[44];
rz(0.4) q[42];
rz(0.4) q[68];
cx q[48], q[70];
cx q[16], q[41];
rx(1) q[37];
rx(1) q[38];
rx(1) q[53];
rx(1) q[66];
cx q[56], q[22];
cx q[76], q[65];
cx q[4], q[49];
cx q[11], q[40];
cx q[29], q[60];
cx q[47], q[32];
cx q[35], q[19];
cx q[18], q[44];
cx q[26], q[42];
cx q[45], q[68];
rz(1.2) q[70];
rz(1.2) q[41];
cx q[28], q[37];
cx q[50], q[65];
cx q[58], q[4];
cx q[47], q[27];
cx q[23], q[35];
cx q[31], q[18];
cx q[73], q[44];
rx(0.2) q[26];
rx(0.2) q[29];
rx(0.2) q[32];
rx(0.2) q[40];
rx(0.2) q[45];
rx(0.2) q[49];
rx(0.2) q[68];
cx q[48], q[70];
cx q[16], q[41];
rz(0.4) q[37];
rz(0.4) q[65];
rz(0.4) q[4];
rz(0.4) q[27];
rz(0.4) q[35];
rz(0.4) q[18];
rz(0.4) q[44];
rx(1) q[16];
rx(1) q[41];
rx(1) q[48];
rx(1) q[70];
cx q[28], q[37];
cx q[50], q[65];
cx q[58], q[4];
cx q[47], q[27];
cx q[23], q[35];
cx q[31], q[18];
cx q[73], q[44];
cx q[10], q[41];
cx q[48], q[56];
cx q[69], q[16];
cx q[50], q[28];
cx q[4], q[61];
cx q[70], q[22];
cx q[71], q[23];
cx q[35], q[12];
cx q[27], q[63];
cx q[47], q[9];
cx q[58], q[21];
rx(0.2) q[31];
rx(0.2) q[44];
rx(0.2) q[65];
rx(0.2) q[73];
rz(0.4) q[41];
rz(0.4) q[56];
rz(0.4) q[16];
rz(0.4) q[28];
rz(0.4) q[61];
rz(0.4) q[22];
rz(0.4) q[23];
rz(0.4) q[12];
rz(0.4) q[63];
rz(0.4) q[9];
rz(0.4) q[21];
cx q[10], q[41];
cx q[48], q[56];
cx q[69], q[16];
cx q[50], q[28];
cx q[4], q[61];
cx q[70], q[22];
cx q[71], q[23];
cx q[35], q[12];
cx q[27], q[63];
cx q[47], q[9];
cx q[58], q[21];
cx q[72], q[56];
cx q[38], q[69];
cx q[48], q[3];
cx q[6], q[10];
cx q[16], q[52];
cx q[70], q[30];
cx q[71], q[61];
cx q[35], q[21];
rx(0.2) q[9];
rx(0.2) q[12];
rx(0.2) q[23];
rx(0.2) q[27];
rx(0.2) q[28];
rx(0.2) q[47];
rx(0.2) q[50];
rx(0.2) q[58];
rx(0.2) q[63];
rz(0.4) q[56];
rz(0.4) q[69];
rz(0.4) q[3];
rz(0.4) q[10];
rz(0.4) q[52];
rz(0.4) q[30];
rz(0.4) q[61];
rz(0.4) q[21];
cx q[72], q[56];
cx q[38], q[69];
cx q[48], q[3];
cx q[6], q[10];
cx q[16], q[52];
cx q[70], q[30];
cx q[71], q[61];
cx q[35], q[21];
cx q[37], q[56];
cx q[6], q[24];
cx q[39], q[69];
cx q[5], q[10];
cx q[41], q[72];
cx q[52], q[76];
cx q[61], q[25];
cx q[3], q[46];
rx(0.2) q[21];
rx(0.2) q[30];
rx(0.2) q[35];
rx(0.2) q[71];
rz(0.4) q[56];
rz(0.4) q[24];
rz(0.4) q[69];
rz(0.4) q[10];
rz(0.4) q[72];
rz(0.4) q[76];
rz(0.4) q[25];
rz(0.4) q[46];
cx q[37], q[56];
cx q[6], q[24];
cx q[39], q[69];
cx q[5], q[10];
cx q[41], q[72];
cx q[52], q[76];
cx q[61], q[25];
cx q[3], q[46];
cx q[7], q[72];
cx q[77], q[37];
cx q[56], q[75];
cx q[41], q[11];
cx q[6], q[74];
cx q[60], q[69];
cx q[39], q[34];
cx q[22], q[76];
cx q[61], q[42];
rx(0.2) q[3];
rx(0.2) q[5];
rx(0.2) q[24];
rx(0.2) q[25];
rx(0.2) q[46];
rx(0.2) q[52];
rz(0.4) q[72];
rz(0.4) q[37];
rz(0.4) q[75];
rz(0.4) q[11];
rz(0.4) q[74];
rz(0.4) q[69];
rz(0.4) q[34];
rz(0.4) q[76];
rz(0.4) q[42];
cx q[7], q[72];
cx q[77], q[37];
cx q[56], q[75];
cx q[41], q[11];
cx q[6], q[74];
cx q[60], q[69];
cx q[39], q[34];
cx q[22], q[76];
cx q[61], q[42];
cx q[37], q[53];
cx q[41], q[4];
cx q[18], q[72];
cx q[77], q[10];
cx q[19], q[76];
rx(0.2) q[6];
rx(0.2) q[7];
rx(0.2) q[11];
rx(0.2) q[22];
rx(0.2) q[34];
rx(0.2) q[39];
rx(0.2) q[42];
rx(0.2) q[56];
rx(0.2) q[60];
rx(0.2) q[61];
rx(0.2) q[69];
rx(0.2) q[74];
rx(0.2) q[75];
rz(0.4) q[53];
rz(0.4) q[4];
rz(0.4) q[72];
rz(0.4) q[10];
rz(0.4) q[76];
cx q[37], q[53];
cx q[41], q[4];
cx q[18], q[72];
cx q[77], q[10];
cx q[19], q[76];
cx q[38], q[53];
cx q[66], q[41];
cx q[70], q[37];
rx(0.2) q[4];
rx(0.2) q[10];
rx(0.2) q[18];
rx(0.2) q[19];
rx(0.2) q[72];
rx(0.2) q[76];
rx(0.2) q[77];
rz(0.4) q[53];
rz(0.4) q[41];
rz(0.4) q[37];
cx q[38], q[53];
cx q[66], q[41];
cx q[70], q[37];
cx q[48], q[70];
cx q[16], q[41];
rx(0.2) q[37];
rx(0.2) q[38];
rx(0.2) q[53];
rx(0.2) q[66];
rz(0.4) q[70];
rz(0.4) q[41];
cx q[48], q[70];
cx q[16], q[41];
rx(0.2) q[16];
rx(0.2) q[41];
rx(0.2) q[48];
rx(0.2) q[70];
