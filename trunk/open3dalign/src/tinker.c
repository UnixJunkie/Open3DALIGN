/*

tinker.c

is part of

Open3DALIGN
-----------

An open-source software aimed at unsupervised molecular alignment

Copyright (C) 2010-2014 Paolo Tosco, Thomas Balle

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For further information, please contact:

Paolo Tosco, PhD
Dipartimento di Scienza e Tecnologia del Farmaco
Universita' degli Studi di Torino
Via Pietro Giuria, 9
10125 Torino (Italy)
Phone:  +39 011 670 7680
Mobile: +39 348 553 7206
Fax:    +39 011 670 7687
E-mail: paolo.tosco@unito.it

*/

#include <include/o3header.h>
#include <include/prog_exe_info.h>
#include <include/proc_env.h>


int compare_bond_list(const void *a, const void *b)
{
  int result;
  const BondList **da = (const BondList **)a;
  const BondList **db = (const BondList **)b;


  result = ((*da)->a[0] > (*db)->a[0]) - ((*da)->a[0] < (*db)->a[0]);
  return (result ? result : ((*da)->a[1] > (*db)->a[1]) - ((*da)->a[1] < (*db)->a[1]));
}


int fill_tinker_types(AtomInfo **atom)
{
  int i;
  int j;
  int k;
  int n;
  int c;
  int o;
  int h;
  int fully_assigned;
  int n_hydro;
  int n_carb;
  int n_nitro;
  int n_oxy;
  int n_sulf;
  int n_alpha;
  int n_os_alpha;
  int n_beta;
  int n_os_beta;
  int n_atoms;
  int missing;
  AtomInfo *temp_atom = NULL;
  AtomInfo *alpha[MAX_BONDS];
  AtomInfo *beta[MAX_BONDS];


  n_atoms = 0;
  while (atom[n_atoms]->atom_type != -1) {
    ++n_atoms;
  }
  fully_assigned = 0;
  while (fully_assigned < 3) {
    for (i = 0; i < n_atoms; ++i) {
      if (atom[i]->tinker_type) {
        continue;
      }
      for (j = 0, n_hydro = 0, n_nitro = 0, n_oxy = 0, n_sulf = 0,
        n_carb = 0; j < atom[i]->n_bonded; ++j) {
        if (!strcmp(atom[atom[i]->bonded[j].num]->element, "H")) {
          ++n_hydro;
        }
        else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "C")) {
          ++n_carb;
        }
        else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "N")) {
          ++n_nitro;
        }
        else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "O")) {
          ++n_oxy;
        }
        else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "S")) {
          ++n_sulf;
        }
      }
      switch (atom[i]->atom_type) {
        case 1:
        /*
        ALKYL CARBON, SP3
        */
        atom[i]->tinker_type = 1;
        break;

        case 2:
        for (j = 0; j < atom[i]->n_bonded; ++j) {
          if ((atom[atom[i]->bonded[j].num]->atom_type == 2)
            || (atom[atom[i]->bonded[j].num]->atom_type == 4)
            || (atom[atom[i]->bonded[j].num]->atom_type == 30)) {
            /*
            VINYLIC CARBON, SP2
            */
            atom[i]->tinker_type = 2;
            break;
          }
        }
        if (!(atom[i]->tinker_type)) {
          /*
          GENERIC SP2 CARBON
          */
          atom[i]->tinker_type = 3;
        }
        break;

        case 3:
        for (j = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          if ((!strcmp(atom[atom[i]->bonded[j].num]->element, "N"))
            && (atom[i]->bonded[j].order == 2)) {
            if (n_nitro == 3) {
              for (k = 0; (k < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++k) {
                temp_atom = atom[atom[i]->bonded[k].num];
                if ((temp_atom == atom[atom[i]->bonded[j].num])
                  || (strcmp(temp_atom->element, "N"))) {
                  continue;
                }
                for (n = 0; n < temp_atom->n_bonded; ++n) {
                  if (temp_atom->bonded[n].order == 2) {
                    n = -1;
                    break;
                  }
                }
                if (n == -1) {
                  /*
                  SP2 CARBON IN C=N
                  */
                  atom[i]->tinker_type = 5;
                }
              }
              if (!(atom[i]->tinker_type)) {
                /*
                GUANIDINE C=N
                */
                atom[i]->tinker_type = 6;
              }
            }
            else {
              /*
              SP2 CARBON IN C=N
              */
              atom[i]->tinker_type = 5;
             }
          }
          else if ((!strcmp(atom[atom[i]->bonded[j].num]->element, "O"))
            && (atom[i]->bonded[j].order == 2)) {
            if ((n_carb == 2) || ((n_carb == 1) && (n_hydro == 1))
              || (n_hydro == 2)) {
              /*
              KETONE/ALD CARBONYL C
              */
              atom[i]->tinker_type = 7;
            }
            else if (((n_carb == 1) && (n_nitro == 1))
              || ((n_hydro == 1) && (n_nitro == 1))
              || ((n_sulf == 1) && (n_nitro == 1))) {
              /*
              AMIDE CARBONYL C
              */
              atom[i]->tinker_type = 8;
            }
            else if (n_nitro == 2) {
              /*
              UREA CARBONYL C
              */
              atom[i]->tinker_type = 10;
            }
            else if (((n_carb == 1) && (n_oxy == 2))
              || ((n_hydro == 1) && (n_oxy == 2))) {
              for (k = 0; k < atom[i]->n_bonded; ++k) {
                if (atom[atom[i]->bonded[k].num]->atom_type == 51) {
                  k = -1;
                  break;
                }
              }
              if (k > 0) {
                /*
                CARBOX AC/ESTER CARB C
                */
                atom[i]->tinker_type = 11;
              }
            }
            else if ((n_nitro == 1) && (n_oxy == 2)) {
              for (n = 0, c = 0; n < atom[i]->n_bonded; ++n) {
                temp_atom = atom[atom[i]->bonded[n].num];
                c = (temp_atom->atom_type == 10);
                if (c) {
                  break;
                }
              }
              if (c) {
                for (n = 0; n < temp_atom->n_bonded; ++n) {
                  if ((!strcmp(atom[temp_atom->bonded[n].num]->element, "O"))
                    || (!strcmp(atom[temp_atom->bonded[n].num]->element, "S"))
                    || (atom[temp_atom->bonded[n].num]->atom_type == 3)) {
                    c = 0;
                    break;
                  }
                }
              }
              for (n = 0; (n < atom[i]->n_bonded) && c; ++n) {
                if (!strcmp(atom[atom[i]->bonded[n].num]->element, "O")
                  && (atom[i]->bonded[n].order == 1)) {
                  temp_atom = atom[atom[i]->bonded[n].num];
                  for (k = 0; k < temp_atom->n_bonded; ++k) {
                    if (!strcmp(atom[temp_atom->bonded[k].num]->element, "N")) {
                      c = 0;
                      break;
                    }
                  }
                }
              }
              if (c) {  
                /*
                CARBAMATE CARBONYL C
                */
                atom[i]->tinker_type = 12;
              }
              else {
                /*
                AMIDE CARBONYL C
                */
                atom[i]->tinker_type = 8;
              }
            }
            else if (n_oxy == 3) {
              /*
              CARBONIC AC/EST CARB C
              */
              atom[i]->tinker_type = 13;
            }
            else if (((n_carb == 1) && (n_sulf == 1))
              || ((n_hydro == 1) && (n_sulf == 1))
              || (n_sulf == 2)) {
              /*
              THIOESTER CARB C=O
              note: also applies to dithiocarbamates
              */
              atom[i]->tinker_type = 14;
            }
            if (!(atom[i]->tinker_type)) {
              /*
              GENERAL CARBONYL CARBON
              */
              atom[i]->tinker_type = 4;
            }
          }
          else if ((!strcmp(atom[atom[i]->bonded[j].num]->element, "S"))
            && (atom[i]->bonded[j].order == 2)) {
            if (atom[atom[i]->bonded[j].num]->atom_type == 16) {
              if (((n_hydro == 1) && (n_oxy == 1))
                || ((n_carb == 1) && (n_oxy == 1))
                || (n_carb == 2) || (n_hydro ==2)
                || ((n_carb == 1) && (n_hydro == 1))
                || (n_oxy == 2)) {
                /*
                THIOESTER C=S
                note: applies also to thioketones/thioaldehydes,
                thiocarbonates and thiocarbamates ionized on nitrogen
                */
                atom[i]->tinker_type = 15;
              }
              else if (((n_hydro == 1) && (n_nitro == 1))
                || ((n_carb == 1) && (n_nitro == 1))
                || ((n_nitro == 1) && (n_oxy == 1))
                || ((n_nitro == 1) && (n_sulf == 2))
                || (n_nitro == 2)) {
                for (n = 0; (n < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++n) {
                  if (!strcmp(atom[atom[i]->bonded[n].num]->element, "N")) {
                    if ((atom[atom[i]->bonded[n].num]->atom_type == 10)
                      || (atom[atom[i]->bonded[n].num]->atom_type == 39)) {
                      /*
                      THIOAMIDE, C=S
                      note: applies also to thioureas and
                      thiocarbamates non-ionized on nitrogen
                      */
                      atom[i]->tinker_type = 16;
                    }
                    else if (atom[atom[i]->bonded[n].num]->atom_type == 62) {
                      /*
                      THIOESTER C=S
                      note: applies also to thioketones/thioaldehydes,
                      thiocarbonates and thiocarbamates ionized on nitrogen
                      */
                      atom[i]->tinker_type = 15;
                    }
                  }
                }
                if (!(atom[i]->tinker_type)) {
                  /*
                  THIOCARB.AC/EST.CARB.C
                  note: applies also to dithiocarbonates
                  */
                  atom[i]->tinker_type = 19;
                }
              }
              else if (((n_hydro == 1) && (n_sulf == 2))
                || ((n_carb == 1) && (n_sulf == 2))
                || (n_sulf == 3)) {
                /*
                THIOCARB.AC/EST.CARB.C
                note: applies also to dithiocarbonates
                */
                atom[i]->tinker_type = 19;
              }
            }
            else if ((atom[atom[i]->bonded[j].num]->atom_type == 18)
              || (atom[atom[i]->bonded[j].num]->atom_type == 73)) {
              /*
              CARBON IN >C=SO2
              */
              atom[i]->tinker_type = 17;
            }
            else if (atom[atom[i]->bonded[j].num]->atom_type == 74) {
              /*
              C IN >C=S=O
              */
              atom[i]->tinker_type = 18;
            }
          }
          else if (atom[atom[i]->bonded[j].num]->atom_type == 75) {
            /*
            C=P
            */
            atom[i]->tinker_type = 20;
          }
        }
        break;

        case 4:
        if ((atom[i]->n_bonded == 2) && (atom[i]->bonded[0].order == 2)
          && (atom[i]->bonded[1].order == 2)) {
          /*
          ALLENIC CARBON
          */
          atom[i]->tinker_type = 22;
        }
        else {
          /*
          ACETYLENIC CARBON
          */
          atom[i]->tinker_type = 21;
        }
        break;

        case 5:
        if ((atom[i]->n_bonded == 1) && (!strcmp
          (atom[atom[i]->bonded[0].num]->element, "C"))) {
          /*
          H ATTACHED TO C
          */
          atom[i]->tinker_type = 23;
        }
        else {
          /*
          H ATTACHED TO SI
          */
          atom[i]->tinker_type = 24;
        }
        break;

        case 6:
        if (atom[i]->n_bonded == 2) {
          for (j = 0, c = 0, n = 0, o = 0, h = 0; j < atom[i]->n_bonded; ++j) {
            if ((atom[atom[i]->bonded[j].num]->atom_type == 1)
              || (atom[atom[i]->bonded[j].num]->atom_type == 20)
              || (atom[atom[i]->bonded[j].num]->atom_type == 22)) {
              ++c;
            }
            else if ((atom[atom[i]->bonded[j].num]->atom_type == 8)
              || (atom[atom[i]->bonded[j].num]->atom_type == 6)
              || (atom[atom[i]->bonded[j].num]->atom_type == 9)
              || (atom[atom[i]->bonded[j].num]->atom_type == 10)
              || (atom[atom[i]->bonded[j].num]->atom_type == 19)
              || (atom[atom[i]->bonded[j].num]->atom_type == 40)
              || (atom[atom[i]->bonded[j].num]->atom_type == 43)) {
              ++n;
            }
            else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "H")) {
              ++h;
            }
            else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "O")) {
              ++o;
            }
          }
          if (((c == 1) && (h == 1))
            || ((n == 1) && (c == 1))
            || ((o == 1) && (c == 1))
            || (c == 2) || (o == 2)) {
            /*
            O ALCOHOL/ETHER
            note: also applies to N-(hydroxy,alkoxy)amines,
            N-(hydroxy,alkoxy)amides, silylethers,
            alkoxyimines and peroxides
            */
            atom[i]->tinker_type = 25;
            break;
          }
          for (j = 0, c = 0, n = 0, o = 0, h = 0; j < atom[i]->n_bonded; ++j) {
            if ((atom[atom[i]->bonded[j].num]->atom_type == 2)
              || (atom[atom[i]->bonded[j].num]->atom_type == 4)
              || (atom[atom[i]->bonded[j].num]->atom_type == 30)) {
              ++c;
            }
            else if (atom[atom[i]->bonded[j].num]->atom_type == 37) {
              temp_atom = atom[atom[i]->bonded[j].num];
              for (k = 0; k < temp_atom->n_bonded; ++k) {
                if (atom[temp_atom->bonded[k].num] == atom[i]) {
                  continue;
                }
                if (atom[temp_atom->bonded[k].num]->atom_type == 38) {
                  k = -1;
                  break;
                }
              }
              if (k > 0) {
                ++c;
              }
              else {
                ++n;
              }
            }
            else if ((atom[atom[i]->bonded[j].num]->atom_type == 63)
              || (atom[atom[i]->bonded[j].num]->atom_type == 64)) {
              temp_atom = atom[atom[i]->bonded[j].num];
              for (k = 0; k < temp_atom->n_bonded; ++k) {
                if (atom[temp_atom->bonded[k].num] == atom[i]) {
                  continue;
                }
                if (!strcmp(atom[temp_atom->bonded[k].num]->element, "N")
                  && (temp_atom->bonded[k].order == 2)) {
                  k = -1;
                  break;
                }
              }
              if (k > 0) {
                ++c;
              }
            }
            else if ((atom[atom[i]->bonded[j].num]->atom_type == 3)
              || (atom[atom[i]->bonded[j].num]->atom_type == 57)) {
              temp_atom = atom[atom[i]->bonded[j].num];
              for (k = 0; k < temp_atom->n_bonded; ++k) {
                if ((!strcmp(atom[temp_atom->bonded[k].num]->element, "N"))
                  && (temp_atom->bonded[k].order == 2)) {
                  ++n;
                  break;
                }
              }
            }
            else if ((atom[atom[i]->bonded[j].num]->atom_type == 1)
              || (atom[atom[i]->bonded[j].num]->atom_type == 19)
              || (atom[atom[i]->bonded[j].num]->atom_type == 20)
              || (atom[atom[i]->bonded[j].num]->atom_type == 22)) {
              ++n;
            }
            else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "H")) {
              ++h;
            }
            else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "O")) {
              ++o;
            }
          }
          if (((c == 1) && (h == 1))
            || ((c == 1) && (n == 1))
            || ((c == 1) && (o == 1))
            || (c == 2) || (o == 2)) {
            /*
            ENOLIC OR PHENOLIC O
            */
            atom[i]->tinker_type = 27;
            break;
          }
          for (j = 0, missing = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
            if (atom[atom[i]->bonded[j].num]->tinker_type) {
              switch (atom[atom[i]->bonded[j].num]->tinker_type) {
                case 85:
                /*
                DIVALENT O-PO3
                */
                atom[i]->tinker_type = 37;
                break;

                case 86:
                /*
                DIVALENT O-PO2
                */
                atom[i]->tinker_type = 38;
                break;

                case 87:
                case 88:
                case 90:
                temp_atom = atom[atom[i]->bonded[j].num];
                for (k = 0, n = 0; k < temp_atom->n_bonded; ++k) {
                  if (!strcmp(atom[temp_atom->bonded[k].num]->element, "O")) {
                    ++n;
                  }
                }
                if (n == 3) {
                  /*
                  DIVALENT O-PO2
                  */
                  atom[i]->tinker_type = 38;
                }
                else if (n == 2) {
                  /*
                  DIVALENT O-P
                  */
                  atom[i]->tinker_type = 39;
                }
                else if (n == 1) {
                  /*
                  DIVALENT O-P
                  */
                  atom[i]->tinker_type = 40;
                }
                break;

                case 4:
                case 8:
                case 11:
                case 12:
                case 13:
                /*
                ESTER/CARBOX ACID -O-
                note: applies also to carbamates,
                carbonates and anhydrides; if it is
                a mixed carbonic/phosphoric anhydride,
                then it is considered as a phosphoric
                ester oxygen
                */
                atom[i]->tinker_type = 26;
                break;
              }
            }
            else if (!fully_assigned) {
              missing = 1;
            }
          }
          for (j = 0; (!missing) && (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
            if (atom[atom[i]->bonded[j].num]->tinker_type) {
              switch (atom[atom[i]->bonded[j].num]->tinker_type) {
                case 5:
                case 133:
                case 170:
                case 171:
                /*
                DIVALENT OXYGEN
                */
                atom[i]->tinker_type = 28;
                break;

                case 15:
                case 16:
                case 141:
                /*
                THIOESTER/THIOACID -O-
                note: applies also to thiocarbamates
                */
                atom[i]->tinker_type = 29;
                break;

                case 150:
                /*
                DIVALENT NITRATE O
                */
                atom[i]->tinker_type = 30;
                break;

                case 151:
                /*
                DIVALENT NITRITE O
                */
                atom[i]->tinker_type = 31;
                break;

                case 68:
                /*
                DIVALENT O-SO3
                */
                atom[i]->tinker_type = 32;
                break;

                case 66:
                /*
                DIVALENT O-SO2
                */
                atom[i]->tinker_type = 33;
                break;

                case 60:
                if ((atom[atom[i]->bonded[j].num]->n_bonded == 2)
                  && (!strcmp(atom[atom[atom[i]->bonded[j].num]->bonded[0].num]->element, "O"))
                  && (!strcmp(atom[atom[atom[i]->bonded[j].num]->bonded[1].num]->element, "O"))) {
                  /*
                  DIVALENT O-S
                  */
                  atom[i]->tinker_type = 34;
                }
                else {
                  /*
                  GENERAL DIVALENT OX-S
                  */
                  atom[i]->tinker_type = 36;
                }
                break;

                case 62:
                /*
                DIVALENT O-SULFOXIDE S
                */
                atom[i]->tinker_type = 35;
                break;
              }
            }
            else {
              missing = 1;
            }
          }
          if ((!missing) && (!(atom[i]->tinker_type))) {
            /*
            GENERAL DIVALENT O
            */
            atom[i]->tinker_type = 41;
          }
        }
        break;

        case 7:
        for (j = 0, missing = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          if (atom[atom[i]->bonded[j].num]->tinker_type) {
            switch  (atom[atom[i]->bonded[j].num]->tinker_type) {
              case 8:
              case 10:
              case 12:
              /*
              CARBONYL O, AMIDES
              note: applies also to ureas
              and carbamates
              */
              atom[i]->tinker_type = 43;
              break;

              case 7:
              /*
              CARBONYL O,ALD/KETONES
              */
              atom[i]->tinker_type = 44;
              break;

              case 11:
              case 13:
              /*
              CARBONYL O,CARB.AC/EST
              note: applies also to carbonates
              */
              atom[i]->tinker_type = 45;
              break;

              case 70:
              /*
              TERMIN OonTETRAC S
              */
              atom[i]->tinker_type = 112;
              break;

              case 151:
              /*
              NITROSO OXYGEN
              */
              atom[i]->tinker_type = 46;
              break;

              case 62:
              /*
              O=S IN SULFOXIDES
              */
              atom[i]->tinker_type = 47;
              break;

              case 69:
              /*
              TERM O-S SULFONE
              */
              atom[i]->tinker_type = 113;
              break;

              case 186:
              /*
              3COR SinTHIOSULFINATE
              */
              atom[i]->tinker_type = 116;
              break;

              case 187:
              /*
              O=S ON S= TO,E.G.,C
              */
              atom[i]->tinker_type = 48;
              break;
            }
          }
          else {
            missing = 1;
          }
          if ((!missing) && (!(atom[i]->tinker_type))) {
            /*
            GENERAL C=O
            */
            atom[i]->tinker_type = 42;
          }
        }
        break;

        case 8:
        atom[i]->tinker_type = 49;
        break;

        case 9:
        for (j = 0; j < atom[i]->n_bonded; ++j) {
          if ((!strcmp(atom[atom[i]->bonded[j].num]->element, "C"))
            && (atom[i]->bonded[j].order == 2)) {
            /*
            N=C IMINES
            */
            atom[i]->tinker_type = 50;
            break;
          }
          else if ((!strcmp(atom[atom[i]->bonded[j].num]->element, "N"))
            && (atom[i]->bonded[j].order == 2)) {
            /*
            N=N AZO COMPOUNDS
            */
            atom[i]->tinker_type = 51;
            break;
          }
        }
        break;

        case 10:
        for (j = 0, missing = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          if ((!(atom[atom[i]->bonded[j].num]->tinker_type)) && (!fully_assigned)) {
            missing = 1;
          }
          if ((atom[atom[i]->bonded[j].num]->tinker_type == 8)
            || (atom[atom[i]->bonded[j].num]->tinker_type == 10)
            || (atom[atom[i]->bonded[j].num]->atom_type == 41)
            || (atom[atom[i]->bonded[j].num]->tinker_type == 12)) {
            /*
            N AMIDE N-C=O
            */
            atom[i]->tinker_type = 52;
          }
        }
        for (j = 0; (!missing) && (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          if ((atom[atom[i]->bonded[j].num]->tinker_type == 16)
            || (atom[atom[i]->bonded[j].num]->tinker_type == 17)
            || (atom[atom[i]->bonded[j].num]->tinker_type == 141)) {
            /*
            N THIOAMIDE N-C=S
            */
            atom[i]->tinker_type = 53;
          }
        }
        for (j = 0; (!missing) && (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          switch  (atom[atom[i]->bonded[j].num]->tinker_type) {
            case 50:
            /*
            NITROGEN IN N-N=C
            */
            atom[i]->tinker_type = 54;
            break;

            case 51:
            /*
            NITROGEN IN N-N=N
            */
            atom[i]->tinker_type = 55;
            break;
          }
        }
        break;

        case 11:
        /*
        FLUORINE
        */
        atom[i]->tinker_type = 56;
        break;

        case 12:
        /*
        CHLORINE
        */
        atom[i]->tinker_type = 57;
        break;

        case 13:
        /*
        BROMINE
        */
        atom[i]->tinker_type = 58;
        break;

        case 14:
        /*
        IODINE
        */
        atom[i]->tinker_type = 59;
        break;

        case 15:
        /*
        S THIOETHER/MERCAPTAN
        */
        atom[i]->tinker_type = 60;
        break;

        case 16:
        /*
        TERMINAL S=C
        */
        atom[i]->tinker_type = 61;
        break;

        case 17:
        for (j = 0; j < atom[i]->n_bonded; ++j) {
          if ((atom[atom[i]->bonded[j].num]->atom_type == 48)
            || (atom[atom[i]->bonded[j].num]->atom_type == 43)) {
            /*
            S,TRICOORD,S=N
            */
            atom[i]->tinker_type = 63;
          }
        }
        if (!(atom[i]->tinker_type)) {
          /*
          SULFUR IN SULFOXIDES
          */
          atom[i]->tinker_type = 62;
          break;
        }
        break;

        case 18:
        if (n_oxy == 4) {
          /*
          SULFATE SULFUR
          */
          atom[i]->tinker_type = 68;
        }
        else if (n_oxy == 3) {
          /*
          SULFONATE SULFUR
          */
          atom[i]->tinker_type = 66;
        }
        else if (n_oxy == 2) {
          if (n_nitro >= 1) {
            /*
            S IN SULFONAMIDES
            */
            atom[i]->tinker_type = 65;
          }
          else if (atom[i]->n_bonded == 3) {
            /*
            SULFONE S=C
            */
            atom[i]->tinker_type = 69;
          }
          else {
            /*
            S SULFONES -SO2-
            */
            atom[i]->tinker_type = 64;
          }
        }
        else if ((n_oxy ==1) && (n_nitro == 1)) {
          /*
          S SULFONE ANALOGS
          */
          atom[i]->tinker_type = 70;
        }
        break;

        case 19:
        /*
        SILICON
        */
        atom[i]->tinker_type = 71;
        break;

        case 20:
        /*
        C IN 4-MEMBERED RINGS
        */
        atom[i]->tinker_type = 72;
        break;

        case 21:
        if (atom[atom[i]->bonded[0].num]->tinker_type || fully_assigned) {
          if (atom[atom[i]->bonded[0].num]->tinker_type == 25) {
            /*
            H-O ALCOHOLS
            */
            atom[i]->tinker_type = 73;
          }
          else if (atom[atom[i]->bonded[0].num]->tinker_type == 27) {
            /*
            H-O ENOLS/PHENOLS
            */
            atom[i]->tinker_type = 103;
          }
          else if (atom[atom[i]->bonded[0].num]->tinker_type == 124) {
            /*
            H HYDROXIDE ANION
            */
            atom[i]->tinker_type = 75;
          }
          else {
            /*
            H-O GENERAL OXYGEN
            */
            atom[i]->tinker_type = 74;
          }
        }
        break;

        case 22:
        /*
        C 3-MEMBERED RING
        */
        atom[i]->tinker_type = 76;
        break;

        case 23:
        if (atom[atom[i]->bonded[0].num]->atom_type == 8) {
          /*
          H-N(SP3)
          */
          atom[i]->tinker_type = 77;
        }
        else if (atom[atom[i]->bonded[0].num]->atom_type == 39) {
          /*
          H-N PYRROLE
          */
          atom[i]->tinker_type = 79;
        }
        else if ((atom[atom[i]->bonded[0].num]->atom_type == 67)
          || (atom[atom[i]->bonded[0].num]->atom_type == 68)) {
          /*
          H-N N-OXIDE
          */
          atom[i]->tinker_type = 80;
        }
        else if (atom[atom[i]->bonded[0].num]->atom_type == 62) {
          /*
          H ON DICOORD, N(-)
          */
          atom[i]->tinker_type = 81;
        }
        else {
          /*
          H-N GENERAL NITROGEN
          */
          atom[i]->tinker_type = 82;
        }
        break;

        case 24:
        if (atom[atom[i]->bonded[0].num]->tinker_type == 26) {
          /*
          H-O CARBOXYLIC ACIDS
          */
          atom[i]->tinker_type = 83;
        }
        else if (atom[atom[i]->bonded[0].num]->tinker_type) {
          /*
          H-O ATTACHED TO P
          */
          atom[i]->tinker_type = 84;
        }
        break;

        case 25:
        if (n_oxy == 4) {
          /*
          P PO4/PHOSPHODIESTER
          */
          atom[i]->tinker_type = 85;
        }
        else if (n_oxy == 3) {
          /*
          TETRACOORDINATE PO3
          */
          atom[i]->tinker_type = 86;
        }
        else if (n_oxy == 2) {
          /*
          TETRACOORDINATE PO2
          */
          atom[i]->tinker_type = 87;
        }
        else if (n_oxy == 1) {
          for (j = 0, n = 0; j < atom[i]->n_bonded; ++j) {
            if ((atom[atom[i]->bonded[j].num]->atom_type == 32)
              || (atom[atom[i]->bonded[j].num]->atom_type == 72)) {
              ++n;
            }
          }
          if (n == 2) {
            /*
            TETRACOORDINATE PO2
            */
            atom[i]->tinker_type = 87;
          }
          else {
            /*
            TETRACOORDINATE PO
            */
            atom[i]->tinker_type = 88;
          }
        }
        if (!(atom[i]->tinker_type)) {
          /*
          GENERALTETRACOORD P
          */
          atom[i]->tinker_type = 89;
        }
        break;

        case 26:
        atom[i]->tinker_type = 90;
        break;

        case 27:
        if (atom[atom[i]->bonded[0].num]->tinker_type == 51) {
          /*
          AZO HYDROGEN
          */
          atom[i]->tinker_type = 91;
        }
        else if (atom[atom[i]->bonded[0].num]->tinker_type) {
          /*
          IMINE HYDROGEN
          */
          atom[i]->tinker_type = 92;
        }
        break;

        case 28:
        if (atom[atom[i]->bonded[0].num]->tinker_type) {
          switch (atom[atom[i]->bonded[0].num]->tinker_type) {
            case 52:
            /*
            AMIDE HYDROGEN
            */
            atom[i]->tinker_type = 93;
            break;

            case 53:
            /*
            THIOAMIDE HYDROGEN
            */
            atom[i]->tinker_type = 94;
            break;

            case 136:
            /*
            H-N ENAMINES
            */
            atom[i]->tinker_type = 95;
            break;

            case 137:
            /*
            H-N H-N-C=N
            */
            atom[i]->tinker_type = 96;
            break;

            case 54:
            /*
            H-N H-N-N=C
            */
            atom[i]->tinker_type = 97;
            break;

            case 55:
            /*
            H-N H-N-N=N
            */
            atom[i]->tinker_type = 98;
            break;

            case 143:
            case 144:
            case 153:
            /*
            H-N SULFONAMIDE
            */
            atom[i]->tinker_type = 99;
            break;

            case 145:
            case 146:
            /*
            H-N PHOSPHONAMIDE
            */
            atom[i]->tinker_type = 100;
            break;

            case 139:
            case 147:
            /*
            H-N TRIPLY BONDED C
            */
            atom[i]->tinker_type = 101;
            break;

            default:
            /*
            GENERAL H ON SP2 N
            */
            atom[i]->tinker_type = 102;
            break;
          }
        }
        break;

        case 29:
        if (atom[atom[i]->bonded[0].num]->tinker_type) {
          switch (atom[atom[i]->bonded[0].num]->tinker_type) {
            case 27:
            /*
            H-O ENOLS/PHENOLS
            */
            atom[i]->tinker_type = 103;
            break;

            case 28:
            /*
            H-O HO-C=N
            */
            atom[i]->tinker_type = 104;
            break;
          }
        }
        break;

        case 30:
        /*
        C=C 4-RING OLEFIN
        */
        atom[i]->tinker_type = 105;
        break;

        case 31:
        /*
        H-OH WATER
        */
        atom[i]->tinker_type = 106;
        break;

        case 32:
        if (atom[i]->n_bonded == 1) {
          temp_atom = atom[atom[i]->bonded[0].num];
          if (temp_atom->tinker_type) {
            switch (temp_atom->tinker_type) {
              case 50:
              /*
              OXIDE O ON CSP2(-)
              */
              atom[i]->tinker_type = 125;
              break;

              case 52:
              /*
              ALKOXIDE O(-)
              */
              atom[i]->tinker_type = 124;
              break;

              case 140:
              /*
              O CARBOX ANION
              */
              atom[i]->tinker_type = 107;
              break;

              case 174:
              case 175:
              case 176:
              case 198:
              case 199:
              case 200:
              /*
              N-OXIDE OXYGEN
              */
              atom[i]->tinker_type = 108;
              break;

              case 149:
              /*
              NITRO OXYGEN
              */
              atom[i]->tinker_type = 109;
              break;

              case 150:
              temp_atom = atom[atom[i]->bonded[0].num];
              for (j = 0, n = 0; j < temp_atom->n_bonded; ++j) {
                if (!strcmp(atom[temp_atom->bonded[j].num]->element, "O")
                  && (atom[temp_atom->bonded[j].num]->n_bonded == 1)) {
                  ++n;
                }
              }
              if (n == 3) {
                /*
                NITRATE ANION O
                */
                atom[i]->tinker_type = 111;
              }
              else if (n == 2) {
                /*
                NITRO O IN NITRATE
                */
                atom[i]->tinker_type = 110;
              }
              break;

              case 70:
              /*
              TERMIN OonTETRAC S
              */
              atom[i]->tinker_type = 112;
              break;

              case 64:
              case 65:
              case 66:
              case 68:
              case 69:
              case 185:
              for (n = 0, c = 0; n < temp_atom->n_bonded; ++n) {
                if ((!strcmp(atom[temp_atom->bonded[n].num]->element, "O"))
                  && (atom[temp_atom->bonded[n].num]->n_bonded == 1)) {
                  ++c;
                }
              }
              if (c == 4) {
                /*
                TERM O IN SO4(-3)
                */
                atom[i]->tinker_type = 115;
              }
              else if (c == 3) {
                /*
                TERM O IN SULFONATES
                */
                atom[i]->tinker_type = 114;
              }
              else {
                /*
                TERM O-S SULFONE
                */
                atom[i]->tinker_type = 113;
              }
              break;

              case 186:
              /*
              TERM O(-)THIOSULFINATE
              */
              atom[i]->tinker_type = 116;
              break;

              case 85:
              case 86:
              case 87:
              case 88:
              temp_atom = atom[atom[i]->bonded[0].num];
              for (j = 0, n = 0; j < temp_atom->n_bonded; ++j) {
                if (strcmp(atom[temp_atom->bonded[j].num]->element, "O")
                  && strcmp(atom[temp_atom->bonded[j].num]->element, "S")) {
                  continue;
                }
                if (atom[temp_atom->bonded[j].num]->n_bonded == 1) {
                  ++n;
                }
              }
              switch (n) {
                case 1:
                /*
                TERM O IN PHOSPHOXIDES
                */
                atom[i]->tinker_type = 117;
                break;

                case 2:
                /*
                TERM O IN PHOSPHINATES
                */
                atom[i]->tinker_type = 118;
                break;

                case 3:
                /*
                TERM O IN PHOSPHONATES
                */
                atom[i]->tinker_type = 119;
                break;

                case 4:
                /*
                TERM O IN PO4/DIESTER
                */
                atom[i]->tinker_type = 120;
                break;
              }
              break;

              case 190:
              /*
              O IN CLO4(-)
              */
              atom[i]->tinker_type = 121;
              break;
            }
          }
        }
        break;

        case 33:
        /*
        H ON O ATTACHED TO S
        */
        atom[i]->tinker_type = 122;
        break;

        case 34:
        /*
        QUATERNARY NSP3(+)
        */
        atom[i]->tinker_type = 123;
        break;

        case 35:
        if (atom[i]->n_bonded == 1) {
          if ((atom[atom[i]->bonded[0].num]->atom_type == 2)
            || (atom[atom[i]->bonded[0].num]->atom_type == 3)
            || (atom[atom[i]->bonded[0].num]->atom_type == 5)
            || (atom[atom[i]->bonded[0].num]->atom_type == 9)
            || (atom[atom[i]->bonded[0].num]->atom_type == 30)
            || (atom[atom[i]->bonded[0].num]->atom_type == 37)
            || (atom[atom[i]->bonded[0].num]->atom_type == 63)
            || (atom[atom[i]->bonded[0].num]->atom_type == 64)
            || (atom[atom[i]->bonded[0].num]->atom_type == 78)) {
            /*
            OXIDE O ON CSP2(-)
            */
            atom[i]->tinker_type = 125;
          }
          else {
            /*
            ALKOXIDE O(-)
            */
            atom[i]->tinker_type = 124;
          }
        }
        break;

        case 36:
        if (atom[i]->n_bonded == 1) {
          switch (atom[atom[i]->bonded[0].num]->atom_type) {
            case 34:
            /*
            H ON QUATERNARY N
            */
            atom[i]->tinker_type = 126;
            break;

            case 80:
            /*
            H ON IMIDAZOLIUM N
            */
            atom[i]->tinker_type = 127;
            break;

            case 58:
            /*
            H ON PYRIDINE N(+)
            */
            atom[i]->tinker_type = 128;
            break;

            case 55:
            /*
            H ON AMIDINIUM N
            */
            atom[i]->tinker_type = 129;
            break;

            case 54:
            /*
            H ON IMINE N(+)
            */
            atom[i]->tinker_type = 130;
            break;

            case 56:
            /*
            H-NH-R GUANIDINIUM
            */
            atom[i]->tinker_type = 131;
            break;

            case 81:
            /*
            H ON N5+,N5A+ORN5B+
            */
            atom[i]->tinker_type = 132;
            break;
          }
        }
        break;

        case 37:
        /*
        C BENZENE/AROMATIC
        */
        atom[i]->tinker_type = 133;
        break;

        case 38:
        /*
        N PYRIDINE
        */
        atom[i]->tinker_type = 134;
        break;

        case 39:
        /*
        N PYRROLE
        */
        atom[i]->tinker_type = 135;
        break;

        case 40:
        for (j = 0, missing = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          temp_atom = atom[atom[i]->bonded[j].num];
          if (!(temp_atom->tinker_type)) {
            missing = 1;
            break;
          }
          if ((temp_atom->tinker_type == 2)
            || (temp_atom->tinker_type == 3)
            || (temp_atom->tinker_type == 5)
            || (temp_atom->tinker_type == 6)
            || (temp_atom->tinker_type == 20)
            || (temp_atom->tinker_type == 105)
            || (temp_atom->tinker_type == 133)
            || (temp_atom->tinker_type == 170)
            || (temp_atom->tinker_type == 171)) {
            for (n = 0; (n < temp_atom->n_bonded) && (!(atom[i]->tinker_type)); ++n) {
              if (atom[temp_atom->bonded[n].num] == atom[i]) {
                continue;
              }
              if (((!strcmp(atom[temp_atom->bonded[n].num]->element, "N"))
                && (temp_atom->bonded[n].order == 2))
                || (atom[temp_atom->bonded[n].num]->atom_type == 38)
                || (atom[temp_atom->bonded[n].num]->atom_type == 54)
                || (atom[temp_atom->bonded[n].num]->atom_type == 58)) {
                /*
                NITROGEN ON N-C=N
                */
                atom[i]->tinker_type = 137;
              }
              else if (!strcmp(atom[temp_atom->bonded[n].num]->element, "P")) {
                /*
                NITROGEN IN N-C=P
                */
                atom[i]->tinker_type = 138;
                break;
              }
            }
            for (n = 0; (n < temp_atom->n_bonded) && (!(atom[i]->tinker_type)); ++n) {
              if (!strcmp(atom[temp_atom->bonded[n].num]->element, "C")) {
                if (temp_atom->bonded[n].order == 3) {
                  /*
                  N-C-C TRIPLE BOND
                  */
                  atom[i]->tinker_type = 139;
                  break;
                }
              }
            }
          }
        }
        if (((!missing) || fully_assigned) && (!(atom[i]->tinker_type))) {
          /*
          NITROGEN ON N-C=C
          */
          atom[i]->tinker_type = 136;
          break;
        }
        break;

        case 41:
        if (n_oxy >= 2) {
          /*
          CARBOXYLATE ANION C
          */
          atom[i]->tinker_type = 140;
        }
        else if (n_sulf >= 2) {
          /*
          C THIOCARBOXYLATE
          */
          atom[i]->tinker_type = 141;
        }
        break;

        case 42:
        /*
        N  TRIPLE BONDED
        */
        atom[i]->tinker_type = 142;
        break;

        case 43:
        for (j = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          switch (atom[atom[i]->bonded[j].num]->tinker_type) {
            case 65:
            /*
            N IN SULFONAMIDES
            */
            atom[i]->tinker_type = 143;
            break;

            case 66:
            /*
            N SULFONAMIDES SO3
            */
            atom[i]->tinker_type = 144;
            break;

            case 87:
            /*
            N IN PHOSPHONAMIDES
            */
            atom[i]->tinker_type = 145;
            break;

            case 86:
            /*
            N PHOSPHONAMIDES,PO3
            */
            atom[i]->tinker_type = 146;
            break;

            case 21:
            /*
            N-(CYANO GROUP)
            */
            atom[i]->tinker_type = 147;
            break;
          }
        }
        break;

        case 44:
        /*
        S THIOPHENE
        */
        atom[i]->tinker_type = 148;
        break;

        case 45:
        if (n_oxy == 2) {
          /*
          NITRO GROUP N
          */
          atom[i]->tinker_type = 149;
        }
        else if (n_oxy == 3) {
          /*
          NITRATE GROUP N
          */
          atom[i]->tinker_type = 150;
        }
        break;

        case 46:
        /*
        NITROSO NITROGEN
        */
        atom[i]->tinker_type = 151;
        break;

        case 47:
        /*
        TERM N IN AZIDO/DIAZO
        */
        atom[i]->tinker_type = 152;
        break;

        case 48:
        /*
        DIVALNforMONOVALOinSO2
        */
        atom[i]->tinker_type = 153;
        break;

        case 49:
        /*
        OXONIUM(TRICOORD)O
        */
        atom[i]->tinker_type = 154;
        break;

        case 50:
        /*
        H ON O+ O
        */
        atom[i]->tinker_type = 155;
        break;

        case 51:
        /*
        OXENIUM (DICOORD) O
        */
        atom[i]->tinker_type = 156;
        break;

        case 52:
        /*
        H ON OXENIUM OXYGEN
        */
        atom[i]->tinker_type = 157;
        break;

        case 53:
        /*
        N IN C=N=N OR -N=N=N
        */
        atom[i]->tinker_type = 158;
        break;

        case 54:
        for (j = 0, missing = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          if (!(atom[atom[i]->bonded[j].num]->tinker_type)) {
            missing = 1;
          }
          if ((atom[atom[i]->bonded[j].num]->tinker_type == 164)
            || (atom[atom[i]->bonded[j].num]->tinker_type == 5)) {
            /*
            IMINIUM NITROGEN
            */
            atom[i]->tinker_type = 159;
          }
        }
        if ((!missing) && (!(atom[i]->tinker_type))) {
          /*
          N(+)=N
          */
          atom[i]->tinker_type = 160;
        }
        break;

        case 55:
        /*
        N IN +N=C-N
        */
        atom[i]->tinker_type = 161;
        break;

        case 56:
        /*
        N GUANIDINIUM
        */
        atom[i]->tinker_type = 162;
        break;

        case 57:
        for (j = 0, missing = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
          if (!(atom[atom[i]->bonded[j].num]->tinker_type)) {
            missing = 1;
          }
          if (atom[atom[i]->bonded[j].num]->tinker_type == 162) {
            /*
            C GUANIDINIUM
            */
            atom[i]->tinker_type = 163;
          }
        }
        if ((!missing) && (!(atom[i]->tinker_type))) {
          /*
          C IN +N=C-N
          */
          atom[i]->tinker_type = 164;
        }
        break;

        case 58:
        /*
        PYRIDINIUM-TYPE N
        */
        atom[i]->tinker_type = 165;
        break;

        case 59:
        /*
        O FURAN
        */
        atom[i]->tinker_type = 166;
        break;

        case 60:
        /*
        ISONITRILE CARBON
        */
        atom[i]->tinker_type = 167;
        break;

        case 61:
        /*
        ISONITRILE N/DIAZO N
        */
        atom[i]->tinker_type = 168;
        break;

        case 62:
        /*
        DEPROT.SULFONAMIDE N-
        */
        atom[i]->tinker_type = 169;
        break;

        case 63:
        /*
        A-Cin5MEMB HETEROAROM
        */
        atom[i]->tinker_type = 170;
        break;

        case 64:
        /*
        B-Cin5MEMB HETEROAROM
        */
        atom[i]->tinker_type = 171;
        break;

        case 65:
        /*
        A-AROMHETEROCYC5RING N
        */
        atom[i]->tinker_type = 172;
        break;

        case 66:
        /*
        B-AROMHETEROCYC5RING N
        */
        atom[i]->tinker_type = 173;
        break;

        case 67:
        /*
        NSP2-OXIDE N
        */
        atom[i]->tinker_type = 174;
        break;

        case 68:
        /*
        NSP3-OXIDE N
        */
        atom[i]->tinker_type = 175;
        break;

        case 69:
        /*
        PYRIDINE N-OXIDE N
        */
        atom[i]->tinker_type = 176;
        break;

        case 70:
        /*
        OXYGEN ON WATER
        */
        atom[i]->tinker_type = 177;
        break;

        case 71:
        if (!strcmp(atom[atom[i]->bonded[0].num]->element, "P")) {
          /*
          H ATACHto3/4CORD P
          */
          atom[i]->tinker_type = 180;
        }
        else if (atom[atom[i]->bonded[0].num]->tinker_type == 63) {
          /*
          H ATACHto4VAL,3COR S=N
          */
          atom[i]->tinker_type = 179;
        }
        else if (!strcmp(atom[atom[i]->bonded[0].num]->element, "S")) {
          /*
          H ATACHtoDIVAL,DICOR S
          */
          atom[i]->tinker_type = 178;
        }
        break;

        case 72:
        temp_atom = atom[atom[i]->bonded[0].num];
        if (!strcmp(temp_atom->element, "P")) {
          for (j = 0, n = 0; (j < temp_atom->n_bonded) && (!n); ++j) {
            n = ((!strcmp(atom[temp_atom->bonded[j].num]->element, "O"))
              && (atom[temp_atom->bonded[j].num]->n_bonded == 1));
          }
          if (!n) {
            /*
            TERM S BONDED TO P
            */
            atom[i]->tinker_type = 181;
          }
        }
        else if (atom[atom[i]->bonded[0].num]->tinker_type == 141) {
          /*
          TER SinTHIOCARBOXYLATE
          */
          atom[i]->tinker_type = 182;
        }
        else if (atom[atom[i]->bonded[0].num]->tinker_type == 186) {
          /*
          TER S IN THIOSULFINATE
          */
          atom[i]->tinker_type = 184;
        }
        if ((atom[atom[i]->bonded[0].num]->tinker_type) && (!(atom[i]->tinker_type))) {
          /*
          TERMINAL SULFUR
          */
          atom[i]->tinker_type = 183;
        }
        break;

        case 73:
        if (n_oxy == 2) {
          for (j = 0; (j < atom[i]->n_bonded) && (!(atom[i]->tinker_type)); ++j) {
            if (atom[atom[i]->bonded[j].num]->atom_type == 3) {
              /*
              SULFONE S=C
              */
              atom[i]->tinker_type = 69;
            }
          }
          if (!(atom[i]->tinker_type)) {
            /*
            S IN SO2(-)
            */
            atom[i]->tinker_type = 185;
          }
        }
        else if ((n_oxy >= 1) && (n_sulf >= 1)) {
          /*
          3COR SinTHIOSULFINATE
          */
          atom[i]->tinker_type = 186;
        }
        break;

        case 74:
        /*
        SULFINYL S,EG.inC=S=O
        */
        atom[i]->tinker_type = 187;
        break;

        case 75:
        /*
        P=C
        */
        atom[i]->tinker_type = 188;
        break;

        case 76:
        /*
        N(-)in,E.G,3-/4AZOLE
        */
        atom[i]->tinker_type = 189;
        break;

        case 77:
        /*
        CL IN CLO4(-)
        */
        atom[i]->tinker_type = 190;
        break;

        case 78:
        /*
        GENERAL C 5MEMBHETAR
        */
        atom[i]->tinker_type = 191;
        break;

        case 79:
        /*
        GENERAL N 5MEMHETERCYC
        */
        atom[i]->tinker_type = 192;
        break;

        case 80:
        /*
        C IMIDAZOLIUM N-C-N
        */
        atom[i]->tinker_type = 193;
        break;

        case 81:
        case 82:
        for (j = 0; (j < atom[i]->n_bonded) && (atom[i]->atom_type == 81); ++j) {
          if (atom[atom[i]->bonded[j].num]->atom_type == 80) {
            /*
            IMIDAZOLIUM N
            */
            atom[i]->tinker_type = 194;
            break;
          }
        }
        if (!(atom[i]->tinker_type)) {
          memset(alpha, 0, sizeof(AtomInfo *) * MAX_BONDS);
          memset(beta, 0, sizeof(AtomInfo *) * MAX_BONDS);
          for (j = 0, n = 0, n_alpha = 0, n_os_alpha = 0; j < atom[i]->n_bonded; ++j) {
            if ((atom[atom[i]->bonded[j].num]->atom_type == 81)
              || (atom[atom[i]->bonded[j].num]->atom_type == 82)
              || ((atom[atom[i]->bonded[j].num]->atom_type >= 63)
              && (atom[atom[i]->bonded[j].num]->atom_type <= 66))
              || (atom[atom[i]->bonded[j].num]->atom_type == 44)
              || (atom[atom[i]->bonded[j].num]->atom_type == 39)
              || (atom[atom[i]->bonded[j].num]->atom_type == 59)
              || (atom[atom[i]->bonded[j].num]->atom_type == 78)
              || (atom[atom[i]->bonded[j].num]->atom_type == 79)) {
              alpha[n] = atom[atom[i]->bonded[j].num];
              ++n;
            }
            if ((atom[atom[i]->bonded[j].num]->atom_type == 81)
              || (atom[atom[i]->bonded[j].num]->atom_type == 82)
              || (atom[atom[i]->bonded[j].num]->atom_type == 65)
              || (atom[atom[i]->bonded[j].num]->atom_type == 66)
              || (atom[atom[i]->bonded[j].num]->atom_type == 44)
              || (atom[atom[i]->bonded[j].num]->atom_type == 39)
              || (atom[atom[i]->bonded[j].num]->atom_type == 59)
              || (atom[atom[i]->bonded[j].num]->atom_type == 79)) {
              ++n_alpha;
            }
            if ((atom[atom[i]->bonded[j].num]->atom_type == 44)
              || (atom[atom[i]->bonded[j].num]->atom_type == 39)
              || (atom[atom[i]->bonded[j].num]->atom_type == 59)) {
              ++n_os_alpha;
            }  
          }
          for (j = 0, n_beta = 0, n_os_beta = 0; j < n; ++j) {
            for (k = 0; k < alpha[j]->n_bonded; ++k) {
              if (atom[alpha[j]->bonded[k].num] == atom[i]) {
                continue;
              }
              if ((atom[alpha[j]->bonded[k].num]->atom_type == 81)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 82)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 65)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 66)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 44)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 39)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 59)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 79)) {
                ++n_beta;
              }
              if ((atom[alpha[j]->bonded[k].num]->atom_type == 44)
                || (atom[atom[i]->bonded[j].num]->atom_type == 39)
                || (atom[alpha[j]->bonded[k].num]->atom_type == 59)) {
                ++n_os_beta;
              }
            }
          }
          if ((n_beta && (!n_os_alpha)) || (n_beta && n_os_beta)) {
            if (atom[i]->atom_type == 81) {
              /*
              POSITIVE N5B N
              */
              atom[i]->tinker_type = 196;
            }
            else {
              /*
              N-OXIDEN 5-RING B-POS
              */
              atom[i]->tinker_type = 199;
            }
          }
          if (n_alpha && (!(atom[i]->tinker_type))) {
            if (atom[i]->atom_type == 81) {
              /*
              POSITIVE N5A N
              */
              atom[i]->tinker_type = 195;
            }
            else {
              /*
              N-OXIDEN 5-RING A-POS
              */
              atom[i]->tinker_type = 198;
            }
          }
          if (!(atom[i]->tinker_type)) {
            if (atom[i]->atom_type == 81) {
              /*
              POSITIVE N5 N
              */
              atom[i]->tinker_type = 197;
            }
            else {
              /*
              N-OXIDN GENER5RINGPOS
              */
              atom[i]->tinker_type = 200;
            }
          }
        }
        break;

        case 87:
        /*
        IRON +2 CATION
        */
        atom[i]->tinker_type = 201;
        break;
        
        case 88:
        /*
        IRON +3 CATION
        */
        atom[i]->tinker_type = 202;
        break;

        case 89:
        /*
        FLUORIDE ANION
        */
        atom[i]->tinker_type = 203;
        break;

        case 90:
        /*
        CHLORIDE ANION
        */
        atom[i]->tinker_type = 204;
        break;

        case 91:
        /*
        BROMIDE ANION
        */
        atom[i]->tinker_type = 205;
        break;

        case 92:
        /*
        LITHIUM CATION
        */
        atom[i]->tinker_type = 206;
        break;

        case 93:
        /*
        SODIUM CATION
        */
        atom[i]->tinker_type = 207;
        break;

        case 94:
        /*
        POTASSIUM CATION
        */
        atom[i]->tinker_type = 208;
        break;

        case 95:
        /*
        DIPOSITIVE ZINC
        */
        atom[i]->tinker_type = 210;
        break;

        case 96:
        /*
        CALCIUM(+2) CATION
        */
        atom[i]->tinker_type = 211;
        break;

        case 97:
        /*
        COPPER(+1) CATION
        */
        atom[i]->tinker_type = 212;
        break;

        case 98:
        /*
        COPPER(+2) CATION
        */
        atom[i]->tinker_type = 213;
        break;

        case 99:
        /*
        MAGNESIUM(+2) CATION
        */
        atom[i]->tinker_type = 214;
        break;
      }
    }
    for (i = 0; i < n_atoms; ++i) {
      if (!(atom[i]->tinker_type)) {
        break;
      }
    }
    if (i == n_atoms) {
      fully_assigned = 10;
    }
    ++fully_assigned;
  }

  return ((fully_assigned < 10) ? FL_UNKNOWN_ATOM_TYPE : 0);
}


int v_intersection(int *v1, int *v2)
{
  /*
  returns how many elements
  are in common between v1 and v2
  */
  int len1 = 0;
  int len2 = 0;
  int i1;
  int i2;
  int is = 0;
  
  
  while (v1[len1] != -1) {
    ++len1;
  }
  while (v2[len2] != -1) {
    ++len2;
  }
  for (i1 = 0; i1 < len1; ++i1) {
    for (i2 = 0; i2 < len2; ++i2) {
      if (v1[i1] == v2[i2]) {
        ++is;
        break;
      }
    }
  }
  
  return is;
}


int v_union(int *v_union, int *v1, int *v2)
{
  /*
  returns the union between
  v1 and v2 (no duplicates)
  */
  int len1 = 0;
  int len2 = 0;
  int flag = 0;
  int i = 0;


  while (v1[len1] != -1) {
    v_union[len1] = v1[len1];
    ++len1;
  }
  while (v2[len2] != -1) {
    i = 0;
    flag = 0;
    while ((v_union[i] != -1) && (!flag)) {
      flag = (v_union[i] == v2[len2]);
      ++i;
    }
    if (!flag) {
      v_union[i] = v2[len2];
    }
    ++len2;
  }
  
  return (i + 1);
}


int bond_in_aromatic_ring(RingInfo **ring, int *a)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int found[2];
  
  
  while (ring[i]) {
    for (j = 0, found[0] = 0, found[1] = 0;
      ((!found[0]) || (!found[1])) && (j < ring[i]->size); ++j) {
      for (k = 0; k < 2; ++k) {
        if (!found[k]) {
          found[k] = (ring[i]->atom_id[j] == a[k]);
        }
      }
    }
    if (found[0] && found[1] && ring[i]->arom && (!((ring[i]->ele - 2) % 4))) {
      return 1;
    }
    ++i;
  }
  
  return 0;
}


int check_bond_type(AtomInfo **atom, int *tinker_types, int *a, int i, BondInfo *bond_info, int value)
{
  int j = 0;
  int k = 0;
  
  
  bond_info[i].num = 0;
  while (tinker_types[j] && (!(bond_info[i].num))) {
    if (atom[a[i]]->tinker_type == tinker_types[j]) {
      bond_info[i].num = value;
      bond_info[i].order = 0;
      for (k = 0; k < atom[a[i]]->n_bonded; ++k) {
        if (atom[a[i]]->bonded[k].num == (a[1 - i])) {
          bond_info[i].order = atom[a[i]]->bonded[k].order;
          break;
        }
      }
      break;
    }
    ++j;
  }
  
  return bond_info[i].num;
}


int check_tinker_type(int type, int *type_array)
{
  int i = 0;
  
  
  while (type_array[i] && (type_array[i] != type)) {
    ++i;
  }
  
  return (type_array[i] ? i : -1);
}

  
int is_in_ring(RingInfo *ring, int atom_id)
{
  int i;
  
  
  for (i = 0; i < ring->size; ++i) {
    if (ring->atom_id[i] == atom_id) {
      return 1;
    }
  }
  
  return 0;
}


int fill_tinker_bond_info(O3Data *od, FileDescriptor *inp_fd, AtomInfo **atom, BondList **bond_list, int object_num)
{
  int i;
  int j;
  int k;
  int n = 0;
  int n_pairs = -1;
  int n_atoms = 0;
  int n_arom_atoms = 0;
  int n_arom_rings = 0;
  int n_bonds = 0;
  int print_bond = 0;
  int fully_aromatic = 0;
  int is = 0;
  int root;
  int front_node;
  int **path = NULL;
  static int sp2_any_tinker_types[] = {
    2, 3, 4, 5, 6, 7, 8, 10, 11, 12,
    13, 14, 15, 16, 19, 21, 50, 51,
    105, 160, 164, 174, 133, 134, 135,
    148, 165, 166, 170, 171, 172, 173,
    176, 191, 192, 193, 194, 195, 196,
    197, 198, 199, 200, 0 };
  static int sp2_non_arom_tinker_types[] = {
    2, 3, 4, 5, 6, 7, 8, 10, 11, 12,
    13, 14, 15, 16, 19, 21, 50, 51,
    105, 160, 164, 174, 0 };
  int sp2_arom_tinker_types[] = {
    133, 134, 135, 148, 165, 166, 170,
    171, 172, 173, 176, 189, 191, 192,
    193, 194, 195, 196, 197, 198, 199,
    200, 0 };
  int sp2_arom_tinker_types_ele[] = {
    1, 1, 2, 2, 1, 2, 1,
    1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1,
    1, 0 };
  BondInfo bond_info[2];
  NodeInfo *fnode = NULL;
  NodeInfo *bnode = NULL;
  NodeInfo *tnode = NULL;
  RingInfo **ring = NULL;


  n_atoms = od->al.mol_info[object_num]->n_atoms;
  n_bonds = od->al.mol_info[object_num]->n_bonds;
  if (!(path = (int **)malloc(n_atoms * sizeof(int *)))) {
    free_node(fnode, path, ring, n_atoms);
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    return FL_OUT_OF_MEMORY;
  }
  memset(path, 0, n_atoms * sizeof(int *));
  for (i = 0; i < n_atoms; ++i) {
    if (atom[i]->ring & AROMATIC) {
      atom[i]->ring |= AROMATIC_TEMP;
      ++n_arom_atoms;
    }
  }
  n_arom_rings = 0;
  while (n_arom_atoms) {
    if (!(ring = (RingInfo **)realloc(ring,
      (n_arom_rings + 2) * sizeof(RingInfo *)))) {
      free_node(fnode, path, ring, n_atoms);
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      return FL_OUT_OF_MEMORY;
    }
    ring[n_arom_rings] = NULL;
    ring[n_arom_rings + 1] = NULL;
    for (i = 0; i < n_atoms; ++i) {
      /*
      count how many aromatic atoms are bound to each aromatic atom
      and put their IDs in the arom_bonded vector
      */
      if (!(atom[i]->ring & AROMATIC_TEMP)) {
        continue;
      }
      atom[i]->n_arom_bonded = 0;
      for (j = 0; j < atom[i]->n_bonded; ++j) {
        if (atom[atom[i]->bonded[j].num]->ring & AROMATIC_TEMP) {
          atom[i]->arom_bonded[atom[i]->n_arom_bonded] =
            atom[i]->bonded[j].num;
          ++(atom[i]->n_arom_bonded);
        }
      }
    }
    for (i = 0; i < n_atoms; ++i) {
      /*
      for each aromatic atom allocate a path vector
      and initialize its positions to -1
      */
      if (!(atom[i]->ring & AROMATIC_TEMP)) {
        continue;
      }
      if (!(path[i] = (int *)realloc(path[i],
        (n_arom_atoms + 1) * sizeof(int)))) {
        free_node(fnode, path, ring, n_atoms);
        return FL_OUT_OF_MEMORY;
      }
      for (j = 0; j < (n_arom_atoms + 1); ++j) {
        path[i][j] = -1;
      }
    }
    if (!(fnode = (NodeInfo *)malloc(sizeof(NodeInfo)))) {
      free_node(fnode, path, ring, n_atoms);
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      return FL_OUT_OF_MEMORY;
    }
    fnode->next = NULL;
    bnode = fnode;
    /*
    allocate a new RingInfo structure and allocate a vector
    for the atom IDs it contains, initializing its positions to -1
    */
    if (!(ring[n_arom_rings] = (RingInfo *)malloc
      (sizeof(RingInfo)))) {
      free_node(fnode, path, ring, n_atoms);
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      return FL_OUT_OF_MEMORY;
    }
    ring[n_arom_rings]->size = 0;
    ring[n_arom_rings]->arom = 0;
    if (!(ring[n_arom_rings]->atom_id = (int *)malloc
      ((n_arom_atoms + 1) * sizeof(int)))) {
      free_node(fnode, path, ring, n_atoms);
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      return FL_OUT_OF_MEMORY;
    }
    for (j = 0; j < (n_arom_atoms + 1); ++j) {
      ring[n_arom_rings]->atom_id[j] = -1;
    }
    /*
    loop over aromatic atoms and stop as soon as one with
    only two connected aromatic atoms is found
    */
    for (i = 0, root = -1; i < n_atoms; ++i) {
      if ((atom[i]->ring & AROMATIC_TEMP) && (atom[i]->n_arom_bonded == 2)) {
        root = i;
        break;
      }
    }
    /*
    if no aromatic atom is found with only two connected aromatic atoms
    check if this is a completely aromatic system (e.g., fullerene)
    */
    if (root == -1) {
      for (i = 0, fully_aromatic = 1; (i < n_atoms) && fully_aromatic; ++i) {
        fully_aromatic = (check_tinker_type(atom[i]->tinker_type, sp2_arom_tinker_types) != -1);
      }
      break;
    }
    /*
    start from the root atom and loop over aromatic atoms
    bound to the root atom
    */
    for (j = 0; j < atom[root]->n_arom_bonded; ++j) {
      /*
      for each aromatic atom bound to the root atom
      allocate a new NodeInfo structure
      and link it to the previous one
      */
      if (j) {
        if (!(bnode->next = (NodeInfo *)malloc(sizeof(NodeInfo)))) {
          free_node(fnode, path, ring, n_atoms);
          O3_ERROR_LOCATE(od->al.task_list[object_num]);
          return FL_OUT_OF_MEMORY;
        }
        bnode = bnode->next;
        bnode->next = NULL;
      }
      bnode->source = root;
      bnode->atom_id = atom[root]->arom_bonded[j];
      path[bnode->atom_id][0] = root;
      path[bnode->atom_id][1] = bnode->atom_id;
    }
    is = 0;
    /*
    fnode is the first node of the chain
    */
    while (fnode) {
      front_node = fnode->atom_id;
      for (j = 0; j < atom[front_node]->n_arom_bonded; ++j) {
        n = atom[front_node]->arom_bonded[j];
        if (n == fnode->source) {
          continue;
        }
        if  (path[n][0] != -1) {
          is = v_intersection(path[front_node], path[n]);
          if (is == 1) {
            ring[n_arom_rings]->size =
              v_union(ring[n_arom_rings]->atom_id,
              path[front_node], path[n]);
            tnode = fnode;
            while (tnode) {
              tnode = fnode->next;
              free(fnode);
              fnode = tnode;
            }
            break;
          }
        }
        else {
          i = 0;
          while (path[n][i] != -1) {
            ++i;
          }
          path[n][i] = n;
          k = 0;
          while (path[front_node][k] != -1) {
            path[n][i + k + 1] = path[front_node][k];
            ++k;
          }
          if (!(bnode->next = (NodeInfo *)malloc(sizeof(NodeInfo)))) {
            free_node(fnode, path, ring, n_atoms);
            O3_ERROR_LOCATE(od->al.task_list[object_num]);
            return FL_OUT_OF_MEMORY;
          }
          bnode = bnode->next;
          bnode->next = NULL;
          bnode->source = front_node;
          bnode->atom_id = n;
        }
      }
      if (ring[n_arom_rings]->size) {
        break;
      }
      if (fnode) {
        tnode = fnode->next;
        free(fnode);
        fnode = tnode;
      }
    }
    /*
    if no ring was found starting from that root atom, then we are
    in the case where two aromatic atoms have remained between two
    aromatic rings, after the other atoms of the ring have been
    removed in the previous step (e.g., COX2 inhibitors after
    removing three of the 5 atoms of the central pentatomic
    aromatic ring)
    */
    if (!(ring[n_arom_rings]->size)) {
      atom[root]->ring &= (~AROMATIC_TEMP);
    }
    for (i = 0, n_arom_atoms = 0; i < n_atoms; ++i) {
      if ((atom[i]->ring & AROMATIC_TEMP) && (atom[i]->n_arom_bonded <= 2)) {
        if ((atom[i]->n_arom_bonded < 2) || is_in_ring(ring[n_arom_rings], i)) {
          atom[i]->ring &= (~AROMATIC_TEMP);
        }
      }
      if (atom[i]->ring & AROMATIC_TEMP) {  
        ++n_arom_atoms;
      }
    }
    free_node(fnode, NULL, NULL, 0);
    ++n_arom_rings;
  }
  for (i = 0; i < n_arom_rings; ++i) {
    ring[i]->arom = 1;
    ring[i]->ele = 0;
    for (j = 0, n = 0; (j < ring[i]->size) && (ring[i]->arom); ++j) {
      ring[i]->arom = (check_tinker_type
        (atom[ring[i]->atom_id[j]]->tinker_type,
        sp2_non_arom_tinker_types) == -1);
      if (ring[i]->arom) {
        k = check_tinker_type
          (atom[ring[i]->atom_id[j]]->tinker_type,
          sp2_arom_tinker_types);
        if (k != -1) {
          /*
          special case for protonated imidazole
          or deprotonated azoles
          */
          if (((sp2_arom_tinker_types[k] == 194)
            || (sp2_arom_tinker_types[k] == 189)) && (n < 2)) {
            ++n;
          }
          ring[i]->ele += sp2_arom_tinker_types_ele[k];
        }
      }
    }
    /*
    special case for protonated imidazole
    or deprotonated azoles
    */
    if (ring[i]->arom && n && (!(n % 2))) {
      ++(ring[i]->ele);
    }
  }
  qsort(bond_list, n_bonds, sizeof(BondList *), compare_bond_list);
  for (n = 0; n < n_bonds; ++n) {
    print_bond = 0;
    if (bond_list[n]->order == 1) {
      for (i = 0; i < 2; ++i) {
        check_bond_type(atom, sp2_any_tinker_types, bond_list[n]->a, i, bond_info, SP2_ANY);
      }
      print_bond = (bond_info[0].num && bond_info[1].num);
    }
    else if (bond_list[n]->order == AROMATIC) {
      for (i = 0; i < 2; ++i) {
        if (!check_bond_type(atom, sp2_non_arom_tinker_types,
          bond_list[n]->a, i, bond_info, SP2_NON_AROMATIC)) {
          check_bond_type(atom, sp2_arom_tinker_types,
            bond_list[n]->a, i, bond_info, AROMATIC);
        }
      }
      print_bond = (!fully_aromatic)
        && ((((bond_info[0].num == SP2_NON_AROMATIC) && (bond_info[1].num == SP2_NON_AROMATIC))
        || ((bond_info[0].num == SP2_NON_AROMATIC) && (bond_info[1].num == AROMATIC))
        || ((bond_info[0].num == AROMATIC) && (bond_info[1].num == SP2_NON_AROMATIC))
        || ((bond_info[0].num == AROMATIC) && (bond_info[1].num == AROMATIC)
        && (!(bond_in_aromatic_ring(ring, bond_list[n]->a)))))
        && (bond_info[0].order == 1) && (bond_info[1].order == 1));
    }
    if (print_bond) {
      if (n_pairs < 0) {
        n_pairs = 0;
      }
      if (!n_pairs) {
        fprintf(inp_fd->handle, "mmff-pibond");
      }
      fprintf(inp_fd->handle, "%6d%6d",
        bond_list[n]->a[0] + 1, bond_list[n]->a[1] + 1);
      ++n_pairs;
      if (n_pairs >= 4) {
        n_pairs = 0;
        fprintf(inp_fd->handle, "\n");
      }
    }
  }
  if (n_pairs > 0) {
    fprintf(inp_fd->handle, "\n");
  }
  free_node(fnode, path, ring, n_atoms);
  
  return 0;
}


int write_tinker_energy(FileDescriptor *fd, double energy)
{
  char buffer[BUF_LEN];
  int i = 0;
  
  
  memset(buffer, 0, BUF_LEN);
  if (!(fd->handle = fopen(fd->name, "rb+"))) {
    return FL_CANNOT_WRITE_INP_FILE;
  }
  else {
    if (!fgets(buffer, BUF_LEN, fd->handle)) {
      return FL_CANNOT_WRITE_INP_FILE;
    }
    else {
      rewind(fd->handle);
      buffer[BUF_LEN - 1] = '\0';
      sprintf(&buffer[8], "%.4lf", energy);
      i = 8;
      while (buffer[i]) {
        ++i;
      }
      while (i < 53) {
        buffer[i] = '_';
        ++i;
      }
      fprintf(fd->handle, "%s", buffer);
    }
    fclose(fd->handle);
  }
  
  return 0;
}


int tinker_analyze(O3Data *od, char *work_dir, char *xyz, int object_num, int conf_num)
{
  char buffer[BUF_LEN];
  int pid;
  int error = 0;
  double energy = 0.0;
  FileDescriptor inp_fd;
  FileDescriptor out_fd;
  FileDescriptor log_fd;
  ProgExeInfo prog_exe_info;


  memset(buffer, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&out_fd, 0, sizeof(FileDescriptor));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.proc_env = fill_env(od, minimal_env, od->qmd.tinker_exe_path, 0);
  if (!(prog_exe_info.proc_env)) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    return FL_OUT_OF_MEMORY;
  }
  sprintf(inp_fd.name, "%s%c%04d_ana.key", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (conf_num != -1) {
    sprintf(buffer, "_%06d", conf_num);
  }
  sprintf(out_fd.name, "%s%c%04d%s.out", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id, buffer);
  sprintf(log_fd.name, "%s%c%04d%s.log", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id, buffer);
  od->al.task_list[object_num]->data[TEMPLATE_OBJECT_NUM] = object_num;
  od->al.task_list[object_num]->data[TEMPLATE_CONF_NUM] = conf_num;
  if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  fprintf(inp_fd.handle,
    "parameters %s%c%s\n"
    "noversion\n"
    "dielectric %lf\n"
    "%s", od->qmd.tinker_prm_path,
    SEPARATOR, TINKER_MMFF94S_PRM_FILE,
    od->qmd.diel_const,
    ((od->qmd.options & QMD_GBSA) ?
    "solvate gbsa\n"
    "solvateterm\n" : ""));
  fclose(inp_fd.handle);
  sprintf(buffer, "%s%c%04d.bnd", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (!fcopy(buffer, inp_fd.name, "ab")) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  prog_exe_info.stdout_fd = &out_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.exedir = work_dir;
  sprintf(prog_exe_info.command_line,
    "%s%c%s -k %04d_ana.key %s e", od->qmd.tinker_exe_path, SEPARATOR,
    TINKER_ANALYZE_EXE, od->al.mol_info[object_num]->object_id, xyz);
  pid = ext_program_exe(&prog_exe_info, &error);
  ext_program_wait(&prog_exe_info, pid);
  free_proc_env(prog_exe_info.proc_env);
  if (error) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
  }
  if (!error) {
    if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      error = FL_CANNOT_READ_OUT_FILE;
    }
    else if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], log_fd.name);
      error = FL_CANNOT_READ_OUT_FILE;
    }
  }
  if (!error) {
    if (!fgrep(out_fd.handle, buffer, "Total Potential Energy")) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      error = FL_ABNORMAL_TERMINATION;
    }
    else if (fgets(buffer, BUF_LEN, log_fd.handle)) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], log_fd.name);
      error = FL_ABNORMAL_TERMINATION;
    }
  }
  if (!error) {
    sscanf(&buffer[25], "%lf", &energy);
    sprintf(inp_fd.name, "%s%c%s", work_dir, SEPARATOR, xyz);
    if ((error = write_tinker_energy(&inp_fd, energy))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    }
  }
  if (out_fd.handle) {
    fclose(out_fd.handle);
    out_fd.handle = NULL;
    if (!error) {
      remove(out_fd.name);
    }
  }
  if (log_fd.handle) {
    fclose(log_fd.handle);
    log_fd.handle = NULL;
    if (!error) {
      remove(log_fd.name);
    }
  }
  
  return error;
}


int tinker_minimize(O3Data *od, char *work_dir, char *xyz, char *xyz_min, int object_num, int conf_num)
{
  char buffer[BUF_LEN];
  int pid;
  int error = 0;
  double energy = 0.0;
  FileDescriptor inp_fd;
  FileDescriptor out_fd;
  FileDescriptor log_fd;
  ProgExeInfo prog_exe_info;


  memset(buffer, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&out_fd, 0, sizeof(FileDescriptor));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.proc_env = fill_env(od, minimal_env, od->qmd.tinker_exe_path, 0);
  if (!(prog_exe_info.proc_env)) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    return FL_OUT_OF_MEMORY;
  }
  sprintf(inp_fd.name, "%s%c%04d_min.key", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (conf_num != -1) {
    sprintf(buffer, "_%06d", conf_num);
  }
  sprintf(out_fd.name, "%s%c%04d%s.out", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id, buffer);
  sprintf(log_fd.name, "%s%c%04d%s.log", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id, buffer);
  od->al.task_list[object_num]->data[TEMPLATE_OBJECT_NUM] = object_num;
  od->al.task_list[object_num]->data[TEMPLATE_CONF_NUM] = conf_num;
  if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  fprintf(inp_fd.handle,
    "parameters %s%c%s\n"
    "noversion\n"
    "dielectric %lf\n"
    "maxiter %d\n"
    "%s", od->qmd.tinker_prm_path,
    SEPARATOR, TINKER_MMFF94S_PRM_FILE,
    od->qmd.diel_const, od->qmd.min_maxiter,
    ((od->qmd.options & QMD_GBSA) ?
    "solvate gbsa\n"
    "solvateterm\n" : ""));
  fclose(inp_fd.handle);
  sprintf(buffer, "%s%c%04d.bnd", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (!fcopy(buffer, inp_fd.name, "ab")) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  prog_exe_info.stdout_fd = &out_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.exedir = work_dir;
  sprintf(prog_exe_info.command_line,
    "%s%c%s -k %04d_min.key %s %lf %s", od->qmd.tinker_exe_path,
    SEPARATOR, od->qmd.minimizer, od->al.mol_info[object_num]->object_id,
    xyz, od->qmd.min_grad, xyz_min);
  pid = ext_program_exe(&prog_exe_info, &error);
  ext_program_wait(&prog_exe_info, pid);
  free_proc_env(prog_exe_info.proc_env);
  if (error) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
  }
  if (!error) {
    if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      error = FL_CANNOT_READ_OUT_FILE;
    }
    else if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], log_fd.name);
      error = FL_CANNOT_READ_OUT_FILE;
    }
  }
  if (!error) {
    if ((!fgrep(out_fd.handle, buffer, TINKER_MINIMIZE_TERMINATION))
      && (!fgrep(out_fd.handle, buffer, "Incomplete Convergence"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      error = FL_ABNORMAL_TERMINATION;
    }
    else if (fgets(buffer, BUF_LEN, log_fd.handle)) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], log_fd.name);
      error = FL_ABNORMAL_TERMINATION;
    }
  }
  if (!error) {
    if (!fgrep(out_fd.handle, buffer, "Final Function Value")) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      error = FL_ABNORMAL_TERMINATION;
    }
  }
  if (!error) {
    sscanf(&buffer[24], "%lf", &energy);
    sprintf(inp_fd.name, "%s%c%s", work_dir, SEPARATOR, xyz_min);
    if ((error = write_tinker_energy(&inp_fd, energy))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    }
  }
  if (out_fd.handle) {
    fclose(out_fd.handle);
    out_fd.handle = NULL;
    if (!error) {
      remove(out_fd.name);
    }
  }
  if (log_fd.handle) {
    fclose(log_fd.handle);
    log_fd.handle = NULL;
    if (!error) {
      remove(log_fd.name);
    }
  }
  
  return error;
}


int tinker_dynamic(O3Data *od, char *work_dir, char *xyz, int object_num, int conf_num, unsigned long seed)
{
  char buffer[BUF_LEN];
  int pid;
  int error = 0;
  FileDescriptor inp_fd;
  FileDescriptor out_fd;
  FileDescriptor log_fd;
  ProgExeInfo prog_exe_info;


  memset(buffer, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&out_fd, 0, sizeof(FileDescriptor));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.proc_env = fill_env(od, minimal_env, od->qmd.tinker_exe_path, 0);
  if (!(prog_exe_info.proc_env)) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    return FL_OUT_OF_MEMORY;
  }
  sprintf(inp_fd.name, "%s%c%04d_dyn.key", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  sprintf(out_fd.name, "%s%c%04d_%06d.out", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id, conf_num);
  sprintf(log_fd.name, "%s%c%04d_%06d.log", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id, conf_num);
  od->al.task_list[object_num]->data[TEMPLATE_OBJECT_NUM] = object_num;
  od->al.task_list[object_num]->data[TEMPLATE_CONF_NUM] = conf_num;
  if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  fprintf(inp_fd.handle,
    "parameters %s%c%s\n"
    "dielectric %lf\n"
    "maxiter %d\n"
    "%s"
    "randomseed %ld\n"
    "integrate verlet\n"
    "thermostat berendsen\n",
    od->qmd.tinker_prm_path,
    SEPARATOR, TINKER_MMFF94_PRM_FILE,
    od->qmd.diel_const, od->qmd.min_maxiter,
    ((od->qmd.options & QMD_GBSA) ?
    "solvate gbsa\n"
    "solvateterm\n" : ""),
    (long)seed);
  fclose(inp_fd.handle);
  inp_fd.handle = NULL;
  sprintf(buffer, "%s%c%04d.bnd", work_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (!fcopy(buffer, inp_fd.name, "ab")) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  prog_exe_info.stdout_fd = &out_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.exedir = work_dir;
  sprintf(prog_exe_info.command_line,
    "%s%c%s -k %04d_dyn.key %s %d %lf %lf 2 %lf",
    od->qmd.tinker_exe_path, SEPARATOR,
    TINKER_DYNAMIC_EXE, od->al.mol_info[object_num]->object_id,
    xyz, (int)safe_rint(od->qmd.window * 1.0e03 / od->qmd.time_step),
    od->qmd.time_step, od->qmd.window, od->qmd.temperature);
  pid = ext_program_exe(&prog_exe_info, &error);
  ext_program_wait(&prog_exe_info, pid);
  free_proc_env(prog_exe_info.proc_env);
  if (error) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
  }
  if (!error) {
    if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      error = FL_CANNOT_READ_OUT_FILE;
    }
    else if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], log_fd.name);
      error = FL_CANNOT_READ_OUT_FILE;
    }
  }
  if (!error) {
    if ((!fgrep(out_fd.handle, buffer, TINKER_DYNAMIC_TERMINATION))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      error = FL_ABNORMAL_TERMINATION;
    }
    else if (fgets(buffer, BUF_LEN, log_fd.handle)) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], log_fd.name);
      error = FL_ABNORMAL_TERMINATION;
    }
  }
  if (out_fd.handle) {
    fclose(out_fd.handle);
    out_fd.handle = NULL;
    if (!error) {
      remove(out_fd.name);
    }
  }
  if (log_fd.handle) {
    fclose(log_fd.handle);
    log_fd.handle = NULL;
    if (!error) {
      remove(log_fd.name);
    }
  }
  
  return error;
}


int write_tinker_xyz_bnd(O3Data *od, AtomInfo **atom, BondList **bond_list, int n_atoms, int object_num, char *xyz_name, char *bnd_name)
{
  int i;
  int j;
  int result = 0;
  FileDescriptor out_fd;
  
  
  memset(&out_fd, 0, sizeof(FileDescriptor));
  if (xyz_name) {
    strcpy(out_fd.name, xyz_name);
    if (!(out_fd.handle = fopen(out_fd.name, "wb+"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      return FL_CANNOT_WRITE_TEMP_FILE;
    }
    /*
    write number of atoms and leave room which will then
    be overwritten by the energy value (nice hack)
    then write coordinates, TINKER types and connectivity table
    */
    fprintf(out_fd.handle,
      "%6d  _____________________________________________\n", n_atoms);
    for (i = 0; i < n_atoms; ++i) {
      fprintf(out_fd.handle,
        "%6d  %-3s%12.6lf%12.6lf%12.6lf%6d", i + 1,
        atom[i]->element, atom[i]->coord[0],
        atom[i]->coord[1], atom[i]->coord[2],
        atom[i]->tinker_type);
      for (j = 0; j < atom[i]->n_bonded; ++j) {
        fprintf(out_fd.handle, "%6d", atom[i]->bonded[j].num + 1);
      }
      fprintf(out_fd.handle, "\n");
    }
    fclose(out_fd.handle);
  }
  if (bnd_name) {
    /*
    in a separate .bnd file write MMFF charges and bond info
    */
    strcpy(out_fd.name, bnd_name);
    if (!(out_fd.handle = fopen(out_fd.name, "wb+"))) {
      O3_ERROR_LOCATE(od->al.task_list[object_num]);
      O3_ERROR_STRING(od->al.task_list[object_num], out_fd.name);
      return FL_CANNOT_WRITE_TEMP_FILE;
    }
    for (i = 0; i < n_atoms; ++i) {
      fprintf(out_fd.handle,
        "charge -%d %lf\n", i + 1,
        ((atom[i]->tinker_type == 201) || (atom[i]->tinker_type == 202)
        || (atom[i]->tinker_type == 212) || (atom[i]->tinker_type == 213)
        ? (double)(atom[i]->sdf_charge) : atom[i]->charge));
    }
    result = fill_tinker_bond_info(od, &out_fd, atom, bond_list, object_num);
    fclose(out_fd.handle);
  }

  return result;
}


void read_tinker_xyz_n_atoms_energy(char *line, int *n_atoms, double *energy)
{
  char buffer[BUF_LEN];
  int i = 0;
  
  
  memset(buffer, 0, BUF_LEN);
  memcpy(buffer, line, BUF_LEN - 1);
  while (buffer[i]) {
    if (buffer[i] == '_') {
      buffer[i] = ' ';
    }
    ++i;
  }
  sscanf(buffer, "%d %lf", n_atoms, energy);
}
