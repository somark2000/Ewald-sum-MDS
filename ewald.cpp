#include "ewald.h"

// k-space Ewald contributions to potential energy and forces.
void EwaldK0(Atom atoms[], int natm, VecR3 Box, real alpha, real& Eele)
{
    const int nkx = 7, nky = 7, nkz = 7;              // number of k-numbers >= 0
    const real gkmin = 1e-10;                     // threshold for Green function

    static complex** eikx, ** eiky, ** eikz, * eikr;      // Fourier basis functions
    static real alph4, Eself, factE, factF, factg, kx1, ky1, kz1;
    static int init = 1;
    complex eikxy, e_ikr, qe_ikr, potk, potki, rhok;
    real Ek, fk, fgkz, gk, kx, ky, kz, k2, Vbox;
    int iatm, ikx, iky, ikz;
    const real pi2 = pi / 2e0;
    // initialization
    if (init) {
        eikx = Mat<complex>(1, natm, -nkx, nkx);
        eiky = Mat<complex>(1, natm, -nky, nky);
        eikz = Mat<complex>(1, natm, 0, nkz);
        eikr = Vec<complex>(1, natm);

        kx1 = pi2 / Box.x;                                // fundamental k numbers
        ky1 = pi2 / Box.y;
        kz1 = pi2 / Box.z;
        Vbox = Box.x * Box.y * Box.z;
        factg = 4e0 * pi * fcoul / Vbox;                      // 1 / (eps0 * Vbox)
        alph4 = 1e0 / (4e0 * alpha * alpha);

        for (iatm = 1; iatm <= natm; iatm++) {// Fourier basis functions for k = 0
            eikx[iatm][0] = complex(1e0, 0e0);
            eiky[iatm][0] = complex(1e0, 0e0);
            eikz[iatm][0] = complex(1e0, 0e0);
        }

        Eself = 0e0;                                                // self energy
        for (iatm = 1; iatm <= natm; iatm++)
            Eself += atoms[iatm].chrg * atoms[iatm].chrg;
        Eself *= fcoul * alpha / sqrt(pi);

        init = 0;
    }
    // generate Fourier basis functions recursively
    for (iatm = 1; iatm <= natm; iatm++) {    // exp(ikx, y, z) for fundamental k
        kx = kx1 * atoms[iatm].r.x; eikx[iatm][1] = complex(cos(kx), sin(kx));
        ky = ky1 * atoms[iatm].r.y; eiky[iatm][1] = complex(cos(ky), sin(ky));
        kz = kz1 * atoms[iatm].r.z; eikz[iatm][1] = complex(cos(kz), sin(kz));
        eikx[iatm][-1] = conj(eikx[iatm][1]);       // exp(-ikx) for fundamental k
        eiky[iatm][-1] = conj(eiky[iatm][1]);       // exp(-iky) for fundamental k

        for (ikx = 2; ikx <= nkx; ikx++) {              // recurrence for exp(ikx)
            eikx[iatm][ikx] = eikx[iatm][ikx - 1] * eikx[iatm][1];
            eikx[iatm][-ikx] = conj(eikx[iatm][ikx]);                  // exp(-ikx)
        }
        for (iky = 2; iky <= nky; iky++) {              // recurrence for exp(iky)
            eiky[iatm][iky] = eiky[iatm][iky - 1] * eiky[iatm][1];
            eiky[iatm][-iky] = conj(eiky[iatm][iky]);                  // exp(-iky)
        }
        for (ikz = 2; ikz <= nkz; ikz++)                // recurrence for exp(ikz)
            eikz[iatm][ikz] = eikz[iatm][ikz - 1] * eikz[iatm][1];
    }
    // k-space Ewald sum contributions
    Ek = 0e0;
    for (ikz = 0; ikz <= nkz; ikz++) {   // invariance k -> -k: take only kz >= 0
        kz = ikz * kz1;
        fgkz = (ikz ? 2e0 * factg : factg);           // factor for Green function
        for (iky = -nky; iky <= nky; iky++) {
            ky = iky * ky1;
            for (ikx = -nkx; ikx <= nkx; ikx++) {
                kx = ikx * kx1;
                k2 = kx * kx + ky * ky + kz * kz;
                if (k2 > 0e0) {
                    gk = fgkz * exp(-alph4 * k2) / k2;       // Green function / Vbox

                    rhok = complex(0e0, 0e0);        // k-space charge density / Vbox
                    for (iatm = 1; iatm <= natm; iatm++) {
                        eikr[iatm] = eikx[iatm][ikx] * eiky[iatm][iky] * eikz[iatm][ikz];
                        rhok += atoms[iatm].chrg * conj(eikr[iatm]);
                    }

                    if (fabs(gk) >= gkmin) {
                        potk = gk * rhok;                  // k-space potential / Vbox
                        // k-space contributions to energy
                        Ek += rhok.real() * potk.real() + rhok.imag() * potk.imag();
                        // k-space contributions to forces
                        for (iatm = 1; iatm <= natm; iatm++) {
                            fk = atoms[iatm].chrg * (potk * eikr[iatm]).imag();
                            atoms[iatm].f.x += fk * kx;
                            atoms[iatm].f.y += fk * ky;
                            atoms[iatm].f.z += fk * kz;
                        }
                    }
                }
            }
        }
    }
    Ek *= 0.5e0;

    Eele += Ek - Eself;
}