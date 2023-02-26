#pragma warning (disable : 4996)

#include "input_output.h"

// Strips blanks from string src[il:i2] and returns the result in res
void strstrip(char src[], char res[], int i1, int i2){
	int j = 0;
	for (int i = i1; i <= i2 and src[i] != 0; ++i)
		if (src[i] != ' ' and src[i] != '\n') res[j++] = src[i];
	res[j] = 0;	
}

// Retrieves the number of atoms from a pdb file
int GetDimPDB(const char* file) {
	FILE* pdb;
	char line[100], string[10];
	int n=0;
	pdb = fopen(file, "r");
	if (!pdb) {
		std::cout << "GetDimPDB: pdb file missing!\n"; _getch(); exit(1);
	}
	strcpy_s(string, "");
	while (!strstr(string, "END")) { // while not finding "END"
		fgets(line, sizeof(line), pdb); // get line
		sscanf_s(line, "%s", string); // first string in line
		if (!strcmp(string, "ATOM") || !strcmp(string, "HETATM")) n += 1;
	}
	fclose(pdb);
	return n;
}

// Reads data from a pdb file in a list of atoms.
void ReadPDB0(const char* file, Atom atoms[], int n) {
	FILE* pdb;
	int i=0;
	char line[100], string[10];
	pdb = fopen(file, "r");
	if (!pdb) {
		std::cout<<("ReadPDBO: pdb file missing!\n"); _getch(); exit(1); }
		strcpy_s(string, "");
		while (!strstr(string, "END")) { // while not finding "END"
			fgets(line, sizeof(line), pdb); // get line
			sscanf_s(line, "%s", string); // first string in line
			if (!strcmp(string, "ATOM") || !strcmp(string, "HETATM")) {
				i += 1;
				if (i > n) {
					printf("ReadPDBO: inconsistent natm ! %i\n", n); exit(1);
				}
				strstrip(line, atoms[i].name, 12, 15);
				strstrip(line, atoms[i].resi, 17, 20);
				strstrip(line, atoms[i].segm, 21, 21);
				strstrip(line, string, 22, 25); atoms[i].ires = strtol(string, NULL, 10);
				strstrip(line, string, 30, 37); atoms[i].r.x = strtof(string, NULL);
				strstrip(line, string, 38, 45); atoms[i].r.y = strtof(string, NULL);
				strstrip(line, string, 46, 53); atoms[i].r.z = strtof(string, NULL);
				strstrip(line, string, 54, 59); atoms[i].occp = strtof(string, NULL);
				strstrip(line, string, 60, 65); atoms[i].beta = strtod(string, NULL);
				strstrip(line, atoms[i].symb, 76, 77);
			}
		}
	fclose(pdb);
}

// Writes data from a list of natm atoms to a pdb file (mode = "w" or "a").
void WritePDB0(const char* file, Atom atoms[], int n, const char* mode) {
	FILE* pdb;
	int i;
	const char* strfrm ="ATOM \t %5d%-4s%-4s%1s%4d \t %8.3f%8.3f%8.3f%6.2f%6.2f \t %-2s\n";
	pdb = fopen(file, mode);
	for (i = 1; i <= n; i++)
		fprintf_s(pdb, strfrm, i, atoms[i].name, atoms[i].resi, atoms[i].segm, atoms[i].ires, atoms[i].r.x, atoms[i].r.y, atoms[i].r.z, atoms[i].occp, atoms[i].beta, atoms[i].symb);
	fprintf_s(pdb, "END\n");
	fclose(pdb);
}

// Retrieves the main dimensions from a psf file
int GetDimPSF0(const char* file)
{
	FILE* psf;
	int n;
	char line[100];
	psf = fopen(file, "r");
	if (!psf) { printf("GetDimPSF: input file missing!\n"); _getch(); exit(1); }
	strcpy_s(line, ""); // search for header "ÎNATOM"
	while (!strstr(line, "!NATOM")) fgets(line, sizeof(line), psf);
		sscanf_s(line, "%d", &n);
	fclose(psf);
	return n;
}


// Reads in data from a psf file
void ReadPSF0(const char* file, Atom atoms[], int n){
	FILE* psf;
	int n1;
	char line[100];
	psf = fopen(file, "r");
	if (!psf) {
		printf("ReadPSF: input file missing!\n"); _getch(); exit(1); }
		strcpy_s(line, ""); // ATOMS - scan for header "ÎNATOM"
		while (!strstr(line, "!NATOM")) fgets(line, sizeof(line), psf); // get lines
			sscanf_s(line, "%d", &n1);
		if (n1 != n) {
			printf("ReadPSF: inconsistent natm ! %i %i\n", n, n1); exit(1);
		}
		for (int i = 1; i <= n; i++) { // read a line for each atom
			fscanf_s(psf, "%*ld %s %ld %s %s %s %lf %lf %*d", atoms[i].segm, &atoms[i].ires, atoms[i].resi, atoms[i].name, atoms[i].type, &atoms[i].chrg, &atoms[i].mass);
		}
		fclose(psf);
	}

// Reads non-bonded force field parameters from a CHARMM-type file
void ReadPAR0(const char* file, Atom atoms[], int natm) {
	FILE* par;
	int found;
	char line[100], str0[20];
	par = fopen(file, "r");
	if (!par) {
		printf("ReadPARO: input file missing! \n"); _getch(); exit(1);
	}
	// NONBONDED
	for (int i = 1; i <= natm; i++) {
		rewind(par); //set pointer to beginning of file
		strcpy_s(str0, "");
		while (strcmp(str0, "NONBONDED")) { //skip lines up to "NONBONDED"
			fgets(line, sizeof(line), par);
			sscanf_s(line, "%s", str0);
		}
		found = 0;
		while (fgets(line, sizeof(line), par)) {
			sscanf_s(line, "%s", str0);
			if (!strcmp(atoms[i].type, str0)) {//check match of searched type
				sscanf_s(line, "%*s %*f %lf %lf", &atoms[i].epsLJ, &atoms[i].RminLJ);
				atoms[i].epsLJ = sqrt(fabs(atoms[i].epsLJ));
				found = 1;
				break;
			}
		}
		if (!found) {
			printf("ReadPAR0: Nonbonded parameters missing for atom type %s\n", atoms[i].type);
			_getch(); exit(1);
		}
	}
	fclose(par);
}

// Reads input
void Input0(int &iopRun, real &Rcut, real &Temp, real &dt, int &nstep, int &nout) {
	FILE* in;
	in = fopen("mdsim.dat", "r");
	if (in == NULL) {
		printf("MDsim.dat missing !\n"); _getch(); exit(1); }
	fscanf_s(in, "%*15c%d", &iopRun);
	fscanf_s(in, "%*15c%lf", &Rcut);
	fscanf_s(in, "%*15c%lf", &Temp);
	fscanf_s(in, "%*15c%lf", &dt);
	fscanf_s(in, "%*15c%d", &nstep);
	fscanf_s(in, "%*15c%d", &nout);
	fclose(in);
}

// Prints main output
void Output0(int istep, int iopRun, int natm, real Rcut, real Temp, real dt, int nstep, int nouț, real ELJ, real Ekin) {
	FILE* out;
	real Etot;
	if (istep == 0) {
		out = fopen("mdsim.out", "w");
		fprintf_s(out," iopRun   =%6d\n", iopRun);
		fprintf_s(out," Rcut (A) =%6.21f\n", Rcut);
		fprintf_s(out," Temp(K)  =%6.21f\n", Temp);
		fprintf_s(out," dt (ps)  =%6g\n", dt);
		fprintf_s(out," nstep    =%6d\n", nstep);
		fprintf_s(out," nout     =%6d\n", nouț);
		fprintf_s(out, "\n");
	}
	else {
		out = fopen("mdsim.out", "a");
	}
	Etot = ELJ + Ekin;
	fprintf_s(out, "%7d%10.4f%15.4f%15.4f%15.4f\n",istep, istep * dt, Etot, ELJ, Ekin);
	fclose(out);
}
