/*

compare.c

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
#ifdef WIN32
#include <windows.h>
#endif


int compare(O3Data *od, O3Data *od_comp, int type, int verbose)
{
  char buffer[BUF_LEN];
  int i;
  int object_num;
  int n_threads;
  int result;
  double heavy_msd = 0.0;
  double ave_heavy_msd = 0.0;
  double t_vec1[RT_VEC_SIZE];
  double t_vec2[RT_VEC_SIZE];
  double rt_mat[RT_MAT_SIZE];
  ConfInfo *template_conf = NULL;
  ConfInfo *fitted_conf = NULL;
  FileDescriptor mol_comp_fd;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  ThreadInfo **ti;


  ti = od->mel.thread_info;
  memset(&mol_comp_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  memset(t_vec1, 0, RT_VEC_SIZE * sizeof(double));
  t_vec1[3] = 1.0;
  if (alloc_threads(od)) {
    return OUT_OF_MEMORY;
  }
  if (!(od->al.task_list = (TaskInfo **)alloc_array(od->grid.object_num, sizeof(TaskInfo)))) {
    return OUT_OF_MEMORY;
  }
  if (!(od->vel.heavy_msd_list = double_vec_alloc(od->grid.object_num))) {
    return OUT_OF_MEMORY;
  }
  for (i = 0, od->field.max_n_heavy_atoms = 0; i < od->grid.object_num; ++i) {
    if (!(od->al.mol_info[i]->atom = (AtomInfo **)
      alloc_array(od->al.mol_info[i]->n_atoms + 1, sizeof(AtomInfo)))) {
      return OUT_OF_MEMORY;
    }
    if (!(od_comp->al.mol_info[i]->atom = (AtomInfo **)
      alloc_array(od_comp->al.mol_info[i]->n_atoms + 1, sizeof(AtomInfo)))) {
      return OUT_OF_MEMORY;
    }
    result = fill_atom_info(od, &(od->task), od->al.mol_info[i]->atom, NULL, i, O3_MMFF94);
    if (result) {
      return result;
    }
    result = fill_atom_info(od_comp, &(od->task), od_comp->al.mol_info[i]->atom, NULL, i, O3_MMFF94);
    if (result) {
      return result;
    }
    if ((!i) || (od->al.mol_info[i]->n_heavy_atoms > od->field.max_n_heavy_atoms)) {
      od->field.max_n_heavy_atoms = od->al.mol_info[i]->n_heavy_atoms;
      od_comp->field.max_n_heavy_atoms = od_comp->al.mol_info[i]->n_heavy_atoms;
    }
  }
  if ((type & BLOCK_COMPARE)) {
    if (!(od->al.rt_list = (RotoTransList **)alloc_array
      (od->grid.object_num, sizeof(RotoTransList)))) {
      return OUT_OF_MEMORY;
    }
  }
  od_comp->al.task_list = od->al.task_list;
  #ifndef WIN32
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
  #endif
  n_threads = fill_thread_info(od, od->grid.object_num);
  for (i = 0; i < n_threads; ++i) {
    memcpy(&(ti[i]->od), od, sizeof(O3Data));
    memcpy(&(ti[i]->od_comp), od_comp, sizeof(O3Data));
    ti[i]->model_type = type;
    ti[i]->thread_num = i;
    /*
    create the i-th thread
    */
    #ifndef WIN32
    od->error_code = pthread_create(&(od->thread_id[i]),
      &thread_attr, (void *(*)(void *))compare_thread, ti[i]);
    if (od->error_code) {
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)compare_thread,
      ti[i], 0, &(od->dwThreadIdArray[i]));
    if (!(od->hThreadArray[i])) {
      return CANNOT_CREATE_THREAD;
    }
    #endif
  }
  #ifndef WIN32
  /*
  free the pthread attribute memory
  */
  pthread_attr_destroy(&thread_attr);
  /*
  wait for all threads to have finished
  */
  for (i = 0; i < n_threads; ++i) {
    od->error_code = pthread_join(od->thread_id[i],
      &(od->thread_result[i]));
    if (od->error_code) {
      return CANNOT_JOIN_THREAD;
    }
  }
  #else
  WaitForMultipleObjects(n_threads, od->hThreadArray, TRUE, INFINITE);
  for (i = 0; i < n_threads; ++i) {
    CloseHandle(od->hThreadArray[i]);
  }
  #endif
  for (i = 0; (i < od->grid.object_num) && (!(od->al.task_list[i]->code)); ++i);
  /*
  if errors did not occur
  */
  if (type & BLOCK_COMPARE) {
    if (!(template_conf = alloc_conf(od->field.max_n_atoms))) {
      return OUT_OF_MEMORY;
    }
    if (!(fitted_conf = alloc_conf(od->field.max_n_atoms))) {
      free_conf(template_conf);
      return OUT_OF_MEMORY;
    }
    rms_algorithm_multi(od, od_comp, rt_mat, &ave_heavy_msd);
    for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
      template_conf->atom = od->al.mol_info[object_num]->atom;
      template_conf->n_atoms = od->al.mol_info[object_num]->n_atoms;
      template_conf->n_heavy_atoms = od_comp->al.mol_info[object_num]->n_heavy_atoms;
      fitted_conf->atom = od->al.mol_info[object_num]->atom;
      fitted_conf->n_atoms = od_comp->al.mol_info[object_num]->n_atoms;
      fitted_conf->n_heavy_atoms = od_comp->al.mol_info[object_num]->n_heavy_atoms;
      for (i = 0, heavy_msd = 0.0; i < od_comp->al.mol_info[object_num]->n_atoms; ++i) {
        /*
        apply the rotation/translation matrix to comp coordinates
        */
        cblas_dcopy(3, od_comp->al.mol_info[object_num]->atom[i]->coord, 1, t_vec1, 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans,
          RT_VEC_SIZE, RT_VEC_SIZE, 1.0, rt_mat, RT_VEC_SIZE,
          t_vec1, 1, 0.0, t_vec2, 1);
        cblas_dcopy(3, t_vec2, 1, &(fitted_conf->coord[i * 3]), 1);
        cblas_dcopy(3, od->al.mol_info[object_num]->atom[i]->coord,
          1, &(template_conf->coord[i * 3]), 1);
      }
      overall_msd(od->al.rt_list[object_num]->sdm, od->al.rt_list[object_num]->pairs,
        fitted_conf, template_conf, &(od->vel.heavy_msd_list->ve[object_num]));
      if ((result = write_aligned_mol(od, od_comp, &(od->task), fitted_conf, object_num))) {
        free_conf(template_conf);
        free_conf(fitted_conf);
        return result;
      }
    }
    free_conf(template_conf);
    free_conf(fitted_conf);
  }
  if (od->file[ASCII_IN]->name[0]) {
    if (!(od->file[MOLFILE_IN]->handle = fopen(od->file[MOLFILE_IN]->name, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), od->file[MOLFILE_IN]->name);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    if (!(od->file[ASCII_IN]->handle = fopen(od->file[ASCII_IN]->name, "wb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), od->file[ASCII_IN]->name);
      return CANNOT_WRITE_ALIGNED_SDF;
    }
    for (i = 0; i < od->grid.object_num; ++i) {
      sprintf(mol_comp_fd.name, "%s%c%04d_aligned.mol",
        od_comp->field.mol_dir, SEPARATOR, od->al.mol_info[i]->object_id);
      if (!(mol_comp_fd.handle = fopen(mol_comp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), mol_comp_fd.name);
        return CANNOT_READ_TEMP_FILE;
      }
      while (fgets(buffer, BUF_LEN, mol_comp_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        remove_newline(buffer);
        fprintf(od->file[ASCII_IN]->handle, "%s\n", buffer);
      }
      result = 0;
      while ((!result) && fgets(buffer, BUF_LEN, od->file[MOLFILE_IN]->handle)) {
        buffer[BUF_LEN - 1] = '\0';
        result = (strstr(buffer, MOL_DELIMITER) ? 1 : 0);
      }
      if (!result) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), od->file[MOLFILE_IN]->name);
        return CANNOT_READ_ORIGINAL_SDF;
      }
      result = 0;
      while ((!result) && fgets(buffer, BUF_LEN, od->file[MOLFILE_IN]->handle)) {
        buffer[BUF_LEN - 1] = '\0';
        remove_newline(buffer);
        fprintf(od->file[ASCII_IN]->handle, "%s\n", buffer);
        result = ((!strncmp(buffer, SDF_DELIMITER, 4)) ? 1 : 0);
      }
      if (!result) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), od->file[MOLFILE_IN]->name);
        return CANNOT_READ_ORIGINAL_SDF;
      }
      fclose(mol_comp_fd.handle);
      mol_comp_fd.handle = NULL;
    }
    fclose(od->file[MOLFILE_IN]->handle);
    od->file[MOLFILE_IN]->handle = NULL;
    fclose(od->file[ASCII_IN]->handle);
    od->file[ASCII_IN]->handle = NULL;
  }
  if (verbose) {
    tee_printf(od,
      "\n"
      "Number of objects:   %d\n\n",
      od->grid.object_num);
    tee_printf(od,
      "---------------------------------------------------------------------------------------------\n"
      "%5s    %-36s%-36s%12s\n"
      "---------------------------------------------------------------------------------------------\n",
      "N", "Reference name", "Name", "RMSD");
    ave_heavy_msd = 0.0;
    for (i = 0; i < od->grid.object_num; ++i) {
      heavy_msd = sqrt(od->vel.heavy_msd_list->ve[i]);
      tee_printf(od, "%5d    %-36s%-36s%12.4lf\n",
        i + 1, od->al.mol_info[i]->object_name,
        od_comp->al.mol_info[i]->object_name, heavy_msd);
      ave_heavy_msd += heavy_msd;
    }
    if (od->grid.object_num > 1) {
      ave_heavy_msd /= (double)(od->grid.object_num);
      tee_printf(od,
        "---------------------------------------------------------------------------------------------\n"
        "%-5s    %-36s%-36s%12.4lf\n",
        "", "", "Average", ave_heavy_msd);
    }
    tee_printf(od,
      "---------------------------------------------------------------------------------------------\n");
  }
  for (i = 0; i < od->grid.object_num; ++i) {
    free_array(od->al.mol_info[i]->atom);
    od->al.mol_info[i]->atom = NULL;
    free_array(od_comp->al.mol_info[i]->atom);
    od_comp->al.mol_info[i]->atom = NULL;
  }
  double_vec_free(od->vel.heavy_msd_list);
  od->vel.heavy_msd_list = NULL;
  free_array(od->al.rt_list);
  od->al.rt_list = NULL;
  
  return 0;
}


#define O3_TEMPLATE      0
#define O3_COMP        1
#define O3_FITTED_LAP      2
#define O3_FITTED_SYST      3
#define O3_PROGRESS      4
#define O3_CAND        5
#define O3_MAX_CONF      6
#define O3_TEMP_SDM1      0
#define O3_TEMP_SDM2      1
#define O3_BEST_SDM_LAP      2
#define O3_BEST_SDM_SYST    3
#define O3_MAX_SDM      4

#ifndef WIN32
void *compare_thread(void *pointer)
#else
DWORD compare_thread(void *pointer)
#endif
{
  char *used[2] = { NULL, NULL };
  int object_num;
  int error = 0;
  int n_atoms;
  int i;
  int j;
  int pairs;
  int pairs_lap;
  int pairs_syst;
  int alloc_fail = 0;
  int **h[2] = { NULL, NULL };
  double msd_lap = 0.0;
  double msd_syst = 0.0;
  LAPInfo li;
  AtomPair *sdm[O3_MAX_SDM] = { NULL, NULL, NULL, NULL };
  AtomPair *best_sdm = NULL;
  ConfInfo *conf[O3_MAX_CONF] = { NULL, NULL, NULL, NULL, NULL, NULL };
  ConfInfo *fitted_conf = NULL;
  ThreadInfo *ti;
  

  ti = (ThreadInfo *)pointer;
  for (i = 0; i < O3_MAX_CONF; ++i) {
    if (!(conf[i] = alloc_conf(ti->od.field.max_n_atoms))) {
      alloc_fail = 1;
    }
  }
  if (alloc_lap_info(&li, ti->od.field.max_n_atoms)) {
    alloc_fail = 1;
  }
  for (i = 0; i < O3_MAX_SDM; ++i) {
    if ((sdm[i] = (AtomPair *)malloc(square(ti->od.field.max_n_atoms) * sizeof(AtomPair)))) {
      memset(sdm[i], 0, square(ti->od.field.max_n_atoms) * sizeof(AtomPair));
    }
    else {
      alloc_fail = 1;
    }
  }
  for (i = 0; i < 2; ++i) {
    if (!(used[i] = malloc(ti->od.field.max_n_atoms))) {
      alloc_fail = 1;
    }
    if (!(h[i] = (int **)alloc_array(ti->od.field.max_n_atoms, MAX_H_BINS * sizeof(int)))) {
      alloc_fail = 1;
    }
  }
  for (object_num = ti->start; (!error) && (object_num <= ti->end); ++object_num) {
    if (alloc_fail) {
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      error = 1;
      continue;
    }
    conf[O3_TEMPLATE]->atom = ti->od.al.mol_info[object_num]->atom;
    n_atoms = ti->od.al.mol_info[object_num]->n_atoms;
    for (i = 0; i < O3_MAX_CONF; ++i) {
      if (i) {
        conf[i]->atom = ti->od_comp.al.mol_info[object_num]->atom;
      }
      conf[i]->n_atoms = n_atoms;
      conf[i]->n_heavy_atoms = ti->od.al.mol_info[object_num]->n_heavy_atoms;
    }
    for (j = 0; j < n_atoms; ++j) {
      cblas_dcopy(3, conf[O3_TEMPLATE]->atom[j]->coord, 1, &(conf[O3_TEMPLATE]->coord[j * 3]), 1);
      cblas_dcopy(3, conf[O3_COMP]->atom[j]->coord, 1, &(conf[O3_COMP]->coord[j * 3]), 1);
    }
    conf[O3_TEMPLATE]->h = h[O3_TEMPLATE];
    compute_conf_h(conf[O3_TEMPLATE]);
    conf[O3_COMP]->h = h[O3_COMP];
    compute_conf_h(conf[O3_COMP]);
    superpose_conf_syst(conf[O3_COMP], conf[O3_TEMPLATE], conf[O3_FITTED_SYST],
      conf[O3_PROGRESS], conf[O3_CAND], sdm[O3_TEMP_SDM1], sdm[O3_TEMP_SDM2],
      sdm[O3_BEST_SDM_SYST], used, NULL, ANGLE_STEP,
      &msd_syst, NULL, &pairs_syst);
    overall_msd(sdm[O3_BEST_SDM_SYST], pairs_syst,
      conf[O3_FITTED_SYST], conf[O3_TEMPLATE], &msd_syst);
    superpose_conf_lap(&li, conf[O3_COMP], conf[O3_TEMPLATE], conf[O3_FITTED_LAP],
      conf[O3_PROGRESS], sdm[O3_TEMP_SDM1], sdm[O3_BEST_SDM_LAP], used, NULL,
      &msd_lap, NULL, &pairs_lap);
    overall_msd(sdm[O3_BEST_SDM_LAP], pairs_lap,
      conf[O3_FITTED_LAP], conf[O3_TEMPLATE], &msd_lap);
    if ((msd_syst - msd_lap) > MSD_THRESHOLD) {
      ti->od.vel.heavy_msd_list->ve[object_num] = msd_lap;
      best_sdm = sdm[O3_BEST_SDM_LAP];
      pairs = pairs_lap;
      fitted_conf = conf[O3_FITTED_LAP];
    }
    else {
      ti->od.vel.heavy_msd_list->ve[object_num] = msd_syst;
      best_sdm = sdm[O3_BEST_SDM_SYST];
      pairs = pairs_syst;
      fitted_conf = conf[O3_FITTED_SYST];
    }
    if (ti->model_type & BLOCK_COMPARE) {
       if (!(ti->od.al.rt_list[object_num]->sdm =
         (AtomPair *)malloc(pairs * sizeof(AtomPair)))) {
         ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        error = 1;
        continue;
      }
      memcpy(ti->od.al.rt_list[object_num]->sdm, best_sdm, pairs * sizeof(AtomPair));
      ti->od.al.rt_list[object_num]->pairs = pairs;
    }
    if (ti->od.file[ASCII_IN]->name[0]) {
      if (!(ti->model_type & BLOCK_COMPARE)) {
        if ((ti->od.al.task_list[object_num]->code = write_aligned_mol
          (&(ti->od), &(ti->od_comp), ti->od.al.task_list[object_num],
          fitted_conf, object_num))) {
          error = 1;
        }
      }
    }
  }
  for (i = 0; i < O3_MAX_SDM; ++i) {
    if (sdm[i]) {
      free(sdm[i]);
    }
  }
  for (i = 0; i < 2; ++i) {
    if (used[i]) {
      free(used[i]);
    }
    free_array(h[i]);
  }
  for (i = 0; i < O3_MAX_CONF; ++i) {
    free_conf(conf[i]);
  }
  free_lap_info(&li);

  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


int write_aligned_mol(O3Data *od, O3Data *od_comp, TaskInfo *task, ConfInfo *fitted_conf, int object_num)
{
  char buffer[BUF_LEN];
  int atom_num;
  int result = 0;
  FileDescriptor inp_fd;
  FileDescriptor mol_comp_fd;


  memset(buffer, 0, BUF_LEN);
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&mol_comp_fd, 0, sizeof(FileDescriptor));
  sprintf(mol_comp_fd.name, "%s%c%04d.mol", od_comp->field.mol_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (!(mol_comp_fd.handle = fopen(mol_comp_fd.name, "rb"))) {
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], mol_comp_fd.name);
    return FL_CANNOT_READ_MOL_FILE;
  }
  sprintf(inp_fd.name, "%s%c%04d_aligned.mol", od_comp->field.mol_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
    fclose(mol_comp_fd.handle);
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], inp_fd.name);
    return FL_CANNOT_WRITE_SDF_FILE;
  }
  if (find_conformation_in_sdf(mol_comp_fd.handle, inp_fd.handle, 0)) {
    fclose(mol_comp_fd.handle);
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], mol_comp_fd.name);
    return FL_CANNOT_READ_MOL_FILE;
  }
  atom_num = 0;
  while (fgets(buffer, BUF_LEN, mol_comp_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    remove_newline(buffer);
    if (atom_num < fitted_conf->n_atoms) {
      if (replace_coord(od_comp->al.mol_info[object_num]->sdf_version,
        buffer, &(fitted_conf->coord[atom_num * 3]))) {
        break;
      }
      ++atom_num;
    }
    fprintf(inp_fd.handle, "%s\n", buffer);
  }
  if (atom_num < fitted_conf->n_atoms) {
    result = FL_CANNOT_READ_MOL_FILE;
    O3_ERROR_LOCATE(od->al.task_list[object_num]);
    O3_ERROR_STRING(od->al.task_list[object_num], mol_comp_fd.name);
  }
  fclose(mol_comp_fd.handle);
  fclose(inp_fd.handle);
  
  return result;
}
