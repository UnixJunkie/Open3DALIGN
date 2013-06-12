/*

conf.c

is part of

Open3DALIGN
-----------

An open-source software aimed at unsupervised molecular alignment

Copyright (C) 2010-2013 Paolo Tosco, Thomas Balle

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


ConfInfo *alloc_conf(int n_atoms)
{
  ConfInfo *conf;
  
  
  conf = (ConfInfo *)malloc(sizeof(ConfInfo));
  if (!conf) {
    return NULL;
  }
  memset(conf, 0, sizeof(ConfInfo));
  if (!(conf->coord = (double *)malloc(n_atoms * 3 * sizeof(double)))) {
    free(conf);
    return NULL;
  }
  memset(conf->coord, 0, n_atoms * 3 * sizeof(double));
  
  return conf;
}


int alloc_lap_info(LAPInfo *li, int max_n_atoms)
{
  int i;
  int alloc_fail = 0;


  for (i = 0; i < O3_MAX_SLOT; ++i) {
    if ((li->array[i] = (int *)malloc(max_n_atoms * sizeof(int)))) {
      memset(li->array[i], 0, max_n_atoms * sizeof(int));
    }
    else {
      alloc_fail = 1;
    }
  }
  if (!(li->cost = (int **)alloc_array(max_n_atoms, max_n_atoms * sizeof(int)))) {
    alloc_fail = 1;
  }
  if (!(li->diff = (double **)alloc_array(max_n_atoms, max_n_atoms * sizeof(double)))) {
    alloc_fail = 1;
  }
  
  return alloc_fail;
}


void free_conf(ConfInfo *conf)
{
  if (conf) {
    if (conf->coord) {
      free(conf->coord);
    }
    free(conf);
  }
}

    
void free_lap_info(LAPInfo *li)
{
  int i;
  
  
  if (li) {
    if (li->array) {
      for (i = 0; i < O3_MAX_SLOT; ++i) {
        if (li->array[i]) {
          free(li->array[i]);
        }
      }
    }
    if (li->cost) {
      free(li->cost);
    }
    if (li->diff) {
      free(li->diff);
    }
  }
}


void free_node(NodeInfo *fnode, int **path, RingInfo **ring, int n_atoms)
{
  int i;
  NodeInfo *tnode = NULL;
  
  
  if (path) {
    for (i = 0; i < n_atoms; ++i) {
      if (path[i]) {
        free(path[i]);
      }
    }
    free(path);
  }
  if (ring) {
    i = 0;
    while (ring[i]) {
      if (ring[i]->atom_id) {
        free(ring[i]->atom_id);
      }
      free(ring[i]);
      ++i;
    }
    free(ring);
  }
  while (fnode) {
    tnode = fnode->next;
    free(fnode);
    fnode = tnode;
  }
}


int check_conf_db(O3Data *od, char *conf_dir, int type, int *wrong_object_num, int *wrong_conf_num)
{
  char buffer[BUF_LEN];
  int i;
  int object_num;
  int line;
  int conf_num;
  int found = 0;
  int result = 0;
  FileDescriptor multi_conf_sdf_fd;
  MolInfo temp_mol_info;
  
  
  memset(&multi_conf_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&temp_mol_info, 0, sizeof(MolInfo));
  memset(buffer, 0, BUF_LEN);
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    if (type == TEMPLATE_DB) {
      for (i = 0, found = 0; (!found) && (i < od->pel.numberlist[OBJECT_LIST]->size); ++i) {
        found = (object_num == (od->pel.numberlist[OBJECT_LIST]->pe[i] - 1));
      }
      if (!found) {
        continue;
      }
    }
    sprintf(multi_conf_sdf_fd.name, "%s%c%04d.sdf", conf_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (!(multi_conf_sdf_fd.handle = fopen(multi_conf_sdf_fd.name, "rb"))) {
      if (type != ANY_DB) {
        *wrong_object_num = object_num;
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), multi_conf_sdf_fd.name);
        return CANNOT_READ_ORIGINAL_SDF;
      }
      else {
        continue;
      }
    }
    line = 0;
    conf_num = 0;
    while ((!result) && fgets(buffer, BUF_LEN, multi_conf_sdf_fd.handle)) {
      ++line;
      buffer[BUF_LEN - 1] = '\0';
      if (line == 4) {
        ++conf_num;
        if (get_n_atoms_bonds(&temp_mol_info,
          multi_conf_sdf_fd.handle, buffer)) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), multi_conf_sdf_fd.name);
          result = PREMATURE_EOF;
        }
        else if ((temp_mol_info.n_atoms != od->al.mol_info[object_num]->n_atoms)
          || (temp_mol_info.n_bonds != od->al.mol_info[object_num]->n_bonds)) {
          *wrong_conf_num = conf_num;
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), multi_conf_sdf_fd.name);
          result = N_ATOM_BOND_MISMATCH;
        }
      }
      else if (line > 4) {
        if (!strncmp(buffer, SDF_DELIMITER, 4)) {
          line = 0;
        }
      }
    }
    fclose(multi_conf_sdf_fd.handle);
    multi_conf_sdf_fd.handle = NULL;
    if (result) {
      *wrong_object_num = object_num;
      return result;
    }
    od->pel.conf_population[type]->pe[object_num] = conf_num;
  }
  
  return 0;
}
