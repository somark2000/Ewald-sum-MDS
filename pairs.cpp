#include "pairs.h"


// Sets up atom pair list for short-range interactions eliminating bonded pairs.
int GetPairNum0(Atom atoms[], int natm, Bond bnds[], int nbnd) {
	int iatm, ibnd, jatm, npair;
	npair = 0;
	for (iatm = 1; iatm <= natm - 1; iatm++) { // loop over all atom pairs
		for (jatm = iatm + 1; jatm <= natm; jatm++) {
			for (ibnd = 1; ibnd <= nbnd; ibnd++) // check if pair forms bond
				if (ComPair(iatm, jatm, bnds[ibnd].indi, bnds[ibnd].indj)) break;
			if (ibnd <= nbnd) continue; // found matching bond? - skip pair
			npair++;
		}
	}
	return npair;
}

// Sets up atom pair list for short-range interactions eliminating bonded pairs.
void PairList0(Atom atoms[], int natm, Bond bnds[], int nbnd, Pair pairs[], int& npair)
{
	int iatm, ibnd, jatm;
	npair = 0;
	for (iatm = 1; iatm <= natm - 1; iatm++) { // loop over all atom oairs
		for (jatm = iatm + 1; jatm <= natm; jatm++) {
			for (ibnd = 1; ibnd <= nbnd; ibnd++) // check if pair forms bond
				if (ComPair(iatm, jatm, bnds[ibnd].indi, bnds[ibnd].indj)) break;
			if (ibnd <= nbnd) continue; // found matching bond? - skip pair
			pairs[npair] = Pair();// define new atom pair
			pairs[npair].indi = iatm;
			pairs[npair].indj = jatm;
		}
	}
}