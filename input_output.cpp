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
		std::cout << "GetDimPDB: pdb file missing!\n";  exit(1);
	}
	strcpy_s(string, "");
	while (!strstr(string, "END")) { // while not finding "END"
		fgets(line, sizeof(line), pdb); // get line
		sscanf_s(line, "%s ", string); // first string in line
		if (!strcmp(string, "ATOM") || !strcmp(string, "HETATM")) n += 1;
	}
	fclose(pdb);
	return n;
}

// Retrieves the main dimensions from a psf file
int GetDimPSFl(const char* file, int ftnbnd) {
	FILE* psf;
	int natm;
	char line[100];
	psf = fopen(file, "r");
	if (!psf) { printf("GetDimPSFl: input file missing!\n");  exit(1); }
	strcpy(line, ""); // search for header "INATOM"
		while (!strstr(line, "INATOM")) fgets(line, sizeof(line), psf);
			sscanf(line, "%d", &natm);
	// search for header "IN3OND"
	while (!strstr(line, "INBOND")) fgets(line, sizeof(line), psf);
	sscanf(line, "%d", ftnbnd);
	fclose(psf);
	return natm;
}

// Reads data from a pdb file in a list of atoms.
void ReadPDB0(const char* file, Atom atoms[], int n) {
	FILE* pdb;
	int i=0;
	char line[100], string[10];
	pdb = fopen(file, "r");
	if (!pdb) {
		std::cout<<("ReadPDBO: pdb file missing!\n");  exit(1); }
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

// Writes data from a list of natm atoms to a pdb file (mode = "w" or "a").
void WritePDB1(const char* file, Atom atoms[], int natm, VecR3 Box, const char* mode)
{
	FILE* pdb;
	int iatm;
	const char* strfrm = "ATOM % 5d % -4s % -4s % ls % 4d % 8.3f % 8.3f % 8.3f % 6.2f % 6.2f % -2s\n";
	pdb = fopen(file, mode);
	if (Box.x * Box.y * Box.z != 0e0)
		fprintf(pdb, "CRYST1%9.3f%9.3f%9.3f 90.00 90.00 90.00\n", Box.x, Box.y, Box.z);
	for (iatm = 1; iatm <= natm; iatm++)
		fprintf(pdb, strfrm, iatm, atoms[iatm].name, atoms[iatm].resi,
			atoms[iatm].segm, atoms[iatm].ires,
			atoms[iatm].r.x, atoms[iatm].r.y, atoms[iatm].r.z,
			atoms[iatm].occp, atoms[iatm].beta, atoms[iatm].symb);
	fprintf(pdb, "END\n");
	fclose(pdb);
}

// Retrieves the main dimensions from a psf file
int GetDimPSF0(const char* file)
{
	FILE* psf;
	int n;
	char line[100];
	psf = fopen(file, "r");
	if (!psf) { printf("GetDimPSF: input file missing!\n");  exit(1); }
	strcpy_s(line, ""); // search for header "ÎNATOM"
	while (!strstr(line, "!NATOM")) fgets(line, sizeof(line), psf);
		sscanf_s(line, "%d", &n);
	fclose(psf);
	return n;
}

// Retrieves the main dimensions from a psf file
int GetDimPSF1(const char* file, int& nbnd)
{
	FILE* psf;
	int natm;
	char line[100];

	psf = fopen(file, "r");
	if (!psf) { printf("GetDimPSF1: input file missing!\n"); exit(1); }

	strcpy(line, "");                               // search for header "!NATOM"
	while (!strstr(line, "!NATOM")) fgets(line, sizeof(line), psf);
	sscanf(line, "%d", &natm);
	// search for header "!NBOND"
	while (!strstr(line, "!NBOND")) fgets(line, sizeof(line), psf);
	sscanf(line, "%d", &nbnd);

	fclose(psf);
	return natm;
}


// Reads in data from a psf file
void ReadPSF0(const char* file, Atom atoms[], int n){
	FILE* psf;
	int n1;
	char line[100];
	psf = fopen(file, "r");
	if (!psf) {
		printf("ReadPSF: input file missing!\n");  exit(1); }
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

void ReadPSF1(const char* file, Atom atoms[], int natm, Bond bnds[], int nbnd) {
	FILE* psf;
	int iatm, ibnd, _natm, _nbnd;
	char line[100];
	psf = fopen(file, "r");	std::cout << "started\n";
	if (!psf) { printf("ReadPSFl: input file missing!\n");  exit(1); }
	strcpy(line, ""); // ATOMS - scan for header "INATOM"
	while (!strstr(line, "!NATOM")) fgets(line, sizeof(line), psf); // get lines
	sscanf(line, "%d", &_natm);	std::cout << "started\n";
	if (_natm != natm) {
		printf("ReadPSFl: inconsistent natm I %i %i\n", natm, _natm); exit(2);
	}
	for (iatm = 1; iatm <= natm; iatm++) { // read a line for each atom
		fscanf(psf, "%*ld %s %ld %s %s %s %lf %lf %*d",
			atoms[iatm].segm, &atoms[iatm].ires, atoms[iatm].resi,
			atoms[iatm].name, atoms[iatm].type, &atoms[iatm].chrg,
			&atoms[iatm].mass);
	}

	//BONDS
	while (!strstr(line, "!NBOND")) fgets(line, sizeof(line), psf);
	sscanf(line, "%d", &_nbnd);
	if (_nbnd != nbnd) {
		printf("ReadPSFl: inconsistent nbnd I %i %i\n", nbnd, _nbnd); exit(3);
	}
	for (ibnd = 1; ibnd <= nbnd; ibnd++)
		fscanf(psf, "%d %d", &bnds[ibnd].indi, &bnds[ibnd].indj);
	fclose(psf);
}

// Reads non-bonded force field parameters from a CHARMM-type file
void ReadPAR0(const char* file, Atom atoms[], int natm) {
	FILE* par;
	int found;
	char line[100], str0[20];
	par = fopen(file, "r");
	if (!par) {
		printf("ReadPARO: input file missing! \n");  exit(1);
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
			 exit(1);
		}
	}
	fclose(par);
}

void ReadPAR1(const char* file, Atom atoms[], int natm, Bond bnds[], int nbnd) {
	FILE* par;
	int iatm, ibnd, jatm, found;
	char line[100], str0[20], strl[20];
	par = fopen(file, "r");
	if (!par) printf("ReadPARl: input file missing!\n");  exit(1);

	// BONDS
	for (ibnd = 1; ibnd <= nbnd; ibnd++) {
		iatm = bnds[ibnd].indi; jatm = bnds[ibnd].indj;
		rewind(par); // set pointer to beginning of file
		strcpy(str0, "");
		while (strcmp(str0, "BONDS")) { // skip lines up to "BONDS"
			fgets(line, sizeof(line), par); // get line
			sscanf(line, "%s", str0); // first string in line
		}

		found = 0; // start parsing
		while (fgets(line, sizeof(line), par)) { // read remaining lines in file
			sscanf(line, "%s", str0); // first string in line
			// break if new section starts
			if (!strcmp(str0, "ATOMS") || !strcmp(str0, "BONDS") || !strcmp(str0, "ANGLES") || !strcmp(str0, "DIHEDRALS") || !strcmp(str0, "IMPROPERS") || !strcmp(str0, "NONBONDED")) break;
			sscanf(line, "%s %s", str0, strl); // check type match both ways
			if ((!strcmp(atoms[iatm].type, str0) && !strcmp(atoms[jatm].type, strl)) || (!strcmp(atoms[iatm].type, strl) && !strcmp(atoms[jatm].type, str0))) {
				sscanf(line, "%*s %*s %lf %lf", &bnds[ibnd].Kb, &bnds[ibnd].b0);
				found = 1;
				break;
			}
		}
		if (!found) {
			printf("ReadPARl: Bond parameters missing for %s-%s\n", atoms[iatm].type, atoms[jatm].type);
			 exit(1);
		}
	}

	// NONBONDED
	for (iatm = 1; iatm <= natm; iatm++) {
		rewind(par); // set pointer to beginning of file
		strcpy(str0, "");
		while (strcmp(str0, "NONBONDED")) { // skip lines up to "NONBOMDED"
			fgets(line, sizeof(line), par); // get line
			sscanf(line, "%s", str0); // first string in line
		}

		found = 0;// start parsing
		while (fgets(line, sizeof(line), par)) { // read remaining lines in file
			sscanf(line, "%s", str0); // first string in line
			if (!strcmp(atoms[iatm].type, str0)) { // check match of searched type
				// RminLJ = Rmin / 2, epsLJ = sqrt(eps)
				sscanf(line, "%*s %*f %lf %lf",
					&atoms[iatm].epsLJ, &atoms[iatm].RminLJ);
				atoms[iatm].epsLJ = sqrt(fabs(atoms[iatm].epsLJ));
				found = 1;
				break;
			}
		}
		if (!found) {
			printf("ReadPARl: Nonbonded parameters missing for atom type %s\n",
				atoms[iatm].type);
			 exit(1);
		}
		fclose(par);
	}
}

// Reads input
void Input(int &iopRun, double &Rcut, double &Temp, double &dt, int &nstep, int &nout) {
	FILE* in;
	in = fopen("mdsim.dat", "r");
	if (in == NULL) {
		printf("MDsim.dat missing !\n");  exit(1); }
	fscanf_s(in, "%*15c%d", &iopRun);
	fscanf_s(in, "%*15c%lf", &Rcut);
	fscanf_s(in, "%*15c%lf", &Temp);
	fscanf_s(in, "%*15c%lf", &dt);
	fscanf_s(in, "%*15c%d", &nstep);
	fscanf_s(in, "%*15c%d", &nout);
	fclose(in);
}


void Input1(int &iopRun, int& PBC, real& Lbox, real& Rcut, real& Temp, real &tauT, real& dt, int& nstep, int &nout) {
	FILE* in;
	in = fopen("mdsim.dat", "r");
	if (in == NULL) { printf("MDsim.dat missing !\n");  exit(1); }
	fscanf(in, "%*15c%d", &iopRun); // run type
	fscanf(in, "%*15c%d", & PBC); // != 0 - periodic boundary conditions
	fscanf(in, "%*15c%lf", &Lbox); // box size (A)
	fscanf(in, "%*15c%lf", &Rcut); // cutoff (A)
	fscanf(in, "%*15c%lf", &Temp); // kinetic temperature (K)
	fscanf(in, "%*15c%lf", &tauT); // thermostat coupling constant (ps)
	fscanf(in, "%*15c%lf", &dt); // integration time step (ps)
	fscanf(in, "%*15c%d", &nstep); // no. of time steps
	fscanf(in, "%*15c%d", &nout); // no. of time steps between output
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

void Output1(int istep, int iopRun, int natm, int PBC, real Lbox, real Rcut, real Temp, real tauT, real dt, int nstep, int nout, real ELI, real Ekin, real Tkin, real Pres, real virial)
{
	FILE* out;
	real Etot;
	if (istep == 0) {
		out = fopen("mdsim.out", "w");
		fprintf(out, " iopRun =   %6d\n", iopRun);
		fprintf(out, " PBC =      %6d\n", PBC);
		fprintf(out, " Lbox (A) = %6.21f\n", Lbox);
		fprintf(out, " Rcut(A) =  %6.21f\n", Rcut);
		fprintf(out, " Temp (K) = %6.21f\n", Temp);
		fprintf(out, " tauT(pc) = %6.21f\n", tauT);
		fprintf(out, " dt (ps) =  %6g\n", dt);
		fprintf(out, " nstep =    %6d\n", nstep);
		fprintf(out, " nout =     %6d\n", nout);

		Etot = ELI + Ekin;
		fprintf(out, "%7d%10.4f%l5.4f%l5.4f%l5.4f%l5.4f%l5.4f%l5.4f\n", istep, istep * dt, Etot, ELI, Ekin, Tkin, Pres, virial);
		fclose(out);
	}
}

void Output2(int istep, int iopRun, int natm, int PBC, real Lbox, real Rcut, real Temp, real tauT, real dt, int nstep, int nout, real Ebnd, real ELI, real Ekin, real Tkin, real Pres, real virial) {
	FILE* out;
	real Etot;
	if (istep == 0) {
		out = fopen("mdsim.out", "w");
		fprintf(out, " iopRun =   %6d\n", iopRun);
		fprintf(out, " PBC =      %6d\n", PBC);
		fprintf(out, " Lbox (A) = %6.21f\n", Lbox);
		fprintf(out, " Rcut(A) =  %6.21f\n", Rcut);
		fprintf(out, " Temp (K) = %6.21f\n", Temp);
		fprintf(out, " tauT(pc) = %6.21f\n", tauT);
		fprintf(out, " dt (ps) =  %6g\n", dt);
		fprintf(out, " nstep =    %6d\n", nstep);
		fprintf(out, " nout =     %6d\n", nout);

		Etot = ELI + Ekin;
		fprintf(out, "%7d%10.4f%l5.4f%l5.4f%l5.4f%l5.4f%l5.4f%l5.4f%l5.4f\n", istep, istep * dt, Ebnd, Etot, ELI, Ekin, Tkin, Pres, virial);
		fclose(out);
	}
}

void Output3(int istep, int iopRun, int natm, int PBC, real Lbox, real Rcut, real Temp, real tauT, real dt, int nstep, int nout, real Ebnd, real ELJ, real Eele, real Ekin, real Tkin, real Pres, real virial) {
	FILE* out;
	real Etot;
	if (istep == 0) {
		out = fopen("mdsim.out", "w");
		fprintf(out, " iopRun =   %6d\n", iopRun);
		fprintf(out, " PBC =      %6d\n", PBC);
		fprintf(out, " Lbox (A) = %6.21f\n", Lbox);
		fprintf(out, " Rcut(A) =  %6.21f\n", Rcut);
		fprintf(out, " Temp (K) = %6.21f\n", Temp);
		fprintf(out, " tauT(pc) = %6.21f\n", tauT);
		fprintf(out, " dt (ps) =  %6g\n", dt);
		fprintf(out, " nstep =    %6d\n", nstep);
		fprintf(out, " nout =     %6d\n", nout);

		Etot = Ebnd + ELJ + Eele + Ekin;
		fprintf(out, "%7d%10.4f%l2.4f%l2.4f%l2.4f%l2.4f%l2.4f%l2.4f%l2.4f%l2.4f\n",
			istep, istep * dt, Etot, Ebnd, ELJ, Eele, Ekin, Tkin, Pres, virial);
		fclose(out);
	}
}