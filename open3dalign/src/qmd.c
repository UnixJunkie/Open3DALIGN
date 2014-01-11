/*

qmd.c

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


int qmd(O3Data *od)
{
  int i;
  int n_threads;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  ThreadInfo **ti;


  if (alloc_threads(od)) {
    return OUT_OF_MEMORY;
  }
  if (!(od->al.task_list = (TaskInfo **)
    alloc_array(od->grid.object_num, sizeof(TaskInfo)))) {
    return OUT_OF_MEMORY;
  }
  od->mel.random_seed_array = (unsigned long *)
    malloc(od->qmd.runs * sizeof(unsigned long));
  if (!(od->mel.random_seed_array)) {
    return OUT_OF_MEMORY;
  }
  for (i = 0; i < od->grid.object_num ; ++i) {
    od->al.mol_info[i]->done = 0;
  }
  set_random_seed(od, od->random_seed);
  for (i = 0; i < od->qmd.runs; ++i) {
    od->mel.random_seed_array[i] = (unsigned long)
      (genrand_real(od) * (double)0x4FFFFFFFUL);
  }
  #ifndef WIN32
  pthread_mutex_init(od->mel.mutex, NULL);
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
  #else
  if (!(*(od->mel.mutex) = CreateMutex(NULL, FALSE, NULL))) {
    return CANNOT_CREATE_THREAD;
  }
  #endif  
  n_threads = fill_thread_info(od, od->grid.object_num);
  ti = od->mel.thread_info;
  for (i = 0; i < n_threads; ++i) {
    /*
    create the i-th thread
    */
    #ifndef WIN32
    od->error_code = pthread_create(&(od->thread_id[i]),
      &thread_attr, (void *(*)(void *))qmd_thread, ti[i]);
    if (od->error_code) {
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)qmd_thread,
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
  pthread_mutex_destroy(od->mel.mutex);
  #else
  WaitForMultipleObjects(n_threads, od->hThreadArray, TRUE, INFINITE);
  for (i = 0; i < n_threads; ++i) {
    CloseHandle(od->hThreadArray[i]);
  }
  CloseHandle(*(od->mel.mutex));
  #endif
  free(od->mel.random_seed_array);
  od->mel.random_seed_array = NULL;
  
  return 0;
}


#define O3_CURR        0
#define O3_FITTED      1
#define O3_FITTED_LAP      1
#define O3_FITTED_SYST      3
#define O3_PROGRESS      2
#define O3_CAND        3
#define O3_MAX_CONF      4
#define O3_MAX_SDM      3

#ifndef WIN32
void *qmd_thread(void *pointer)
#else
DWORD qmd_thread(void *pointer)
#endif
{
  char buffer[BUF_LEN];
  char buffer2[BUF_LEN];
  char work_dir[BUF_LEN];
  char *used[2] = { NULL, NULL };
  int i;
  int j;
  int n;
  int run = 0;
  int atom_num = 0;
  int object_num;
  int n_atoms;
  int n_conf = 0;
  int min_pos = 0;
  int result;
  int pairs = 0;
  int alloc_fail = 0;
  int restart = 0;
  int maybe_restart = 0;
  int assigned = 1;
  double heavy_msd_lap = 0.0;
  double heavy_msd_syst = 0.0;
  double min_heavy_msd = 0.0;
  double user_msd;
  double glob_min = 0.0;
  double energy = 0.0;
  double rt_mat[RT_MAT_SIZE];
  LAPInfo li;
  AtomPair *sdm[O3_MAX_SDM] = { NULL, NULL, NULL };
  AtomInfo **atom = NULL;
  BondList **bond_list = NULL;
  ConfInfo *conf[O3_MAX_CONF] = { NULL, NULL, NULL, NULL };
  ConfInfo **conf_array = NULL;
  ConfInfo *fitted_conf = NULL;
  ThreadInfo *ti;
  FileDescriptor mol_fd;
  FileDescriptor inp_fd;


  ti = (ThreadInfo *)pointer;
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  memset(buffer2, 0, BUF_LEN);
  memset(work_dir, 0, BUF_LEN);
  memset(&li, 0, sizeof(LAPInfo));
  /*
  allocate memory for AtomInfo structure array
  */
  if (!(atom = (AtomInfo **)alloc_array(ti->od.field.max_n_atoms + 1, sizeof(AtomInfo)))) {
    alloc_fail = 1;
  }
  if (!(bond_list = (BondList **)alloc_array(ti->od.field.max_n_bonds + 1, sizeof(BondList)))) {
    alloc_fail = 1;
  }
  if ((conf_array = (ConfInfo **)malloc((ti->od.qmd.runs + 2) * sizeof(ConfInfo *)))) {
    memset(conf_array, 0, (ti->od.qmd.runs + 2) * sizeof(ConfInfo *));
  }
  else {
    alloc_fail = 1;
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
  }
  for (i = 0; i < O3_MAX_CONF; ++i) {
    if (!(conf[i] = alloc_conf(ti->od.field.max_n_atoms))) {
      alloc_fail = 1;
    }
  }
  while (assigned) {
    object_num = 0;
    assigned = 0;
    while ((!assigned) && (object_num < ti->od.grid.object_num)) {
      if (!(ti->od.al.mol_info[object_num]->done)) {
        #ifndef WIN32
        pthread_mutex_lock(ti->od.mel.mutex);
        #else
        WaitForSingleObject(ti->od.mel.mutex, INFINITE);
        #endif
        if (!(ti->od.al.mol_info[object_num]->done)) {
          ti->od.al.mol_info[object_num]->done = 1;
          assigned = 1;
        }
        #ifndef WIN32
        pthread_mutex_unlock(ti->od.mel.mutex);
        #else
        ReleaseMutex(ti->od.mel.mutex);
        #endif
      }
      if (!assigned) {
        ++object_num;
      }
    }
    if (!assigned)  {
      break;
    }
    if (alloc_fail) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
      continue;
    }
    if ((ti->od.al.task_list[object_num]->code =
      fill_atom_info(&(ti->od), ti->od.al.task_list[object_num],
        atom, bond_list, object_num, O3_MMFF94))) {
      continue;
    }
    if ((ti->od.al.task_list[object_num]->code = fill_tinker_types(atom))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      continue;
    }
    n_atoms = ti->od.al.mol_info[object_num]->n_atoms;
    for (i = 0; i < O3_MAX_CONF; ++i) {
      conf[i]->atom = atom;
      conf[i]->n_atoms = n_atoms;
      conf[i]->n_heavy_atoms = ti->od.al.mol_info[object_num]->n_heavy_atoms;
    }
    sprintf(work_dir, "%s%c%04d", ti->od.qmd.qmd_dir, SEPARATOR,
      ti->od.al.mol_info[object_num]->object_id);
    sprintf(buffer, "%s.sdf", work_dir);
    /*
    if we are restarting a run, the %04d.sdf file might already exist
    if so, directly skip to the next object
    */
    if (fexist(buffer)) {
      continue;
    }
    /*
    if we are restarting a run, the folder might already exist
    if so, then it is necessary to check where the previous run
    stopped
    */
    maybe_restart = dexist(work_dir);
    if (!maybe_restart) {
      #ifndef WIN32
      result = mkdir(work_dir, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(work_dir);
      #endif
      if (result) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], work_dir);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_CREATE_SCRDIR;
        continue;
      }
    }
    /*
    prepare the initial XYZ geometry
    */
    sprintf(buffer, "%s%c%04d%c%04d.xyz", ti->od.qmd.qmd_dir, SEPARATOR,
      ti->od.al.mol_info[object_num]->object_id, SEPARATOR,
      ti->od.al.mol_info[object_num]->object_id);
    sprintf(buffer2, "%s%c%04d%c%04d.bnd", ti->od.qmd.qmd_dir, SEPARATOR,
      ti->od.al.mol_info[object_num]->object_id, SEPARATOR,
      ti->od.al.mol_info[object_num]->object_id);
    ti->od.al.task_list[object_num]->code = write_tinker_xyz_bnd
      (&(ti->od), atom, bond_list, n_atoms, object_num, buffer, buffer2);
    if (ti->od.al.task_list[object_num]->code) {
      continue;
    }
    run = 0;
    restart = 0;
    if (maybe_restart) {
      for (restart = 1; restart && (run <= ti->od.qmd.runs); ++run) {
        /*
        check all existing XYZ files; as soon as a malformed one is found,
        exit the loop
        directly from scratch
        */
        sprintf(inp_fd.name, "%s%c%04d%c%04d_%06d.xyz", ti->od.qmd.qmd_dir, SEPARATOR,
           ti->od.al.mol_info[object_num]->object_id, SEPARATOR,
           ti->od.al.mol_info[object_num]->object_id, run);
        /*
        if the file already exists, check that it is not malformed
        if it is ok, skip ahead
        */
        if ((restart = fexist(inp_fd.name))) {
          inp_fd.handle = fopen(inp_fd.name, "rb");
          restart = (inp_fd.handle ? 1 : 0);
          if (restart) {
            i = 0;
            while ((i < (n_atoms + 1)) && fgets(buffer2, BUF_LEN, inp_fd.handle)) {
              ++i;
            }
            fclose(inp_fd.handle);
            restart = ((i == (n_atoms + 1)) ? 1 : 0);
          }
        }
      }
    }
    if (!restart) {
      /*
      if a malformed file was found, start 3 runs before that point,
      if possible, otherwise just start from scratch
      */
      restart = run - 3;
      if (restart < 2) {
        restart = -1;
      }
      for (run = restart + 1; maybe_restart && (run < ti->od.qmd.runs); ++run) {
        /*
        removed all files eventually created during previous runs
        which are malformed or not needed anymore
        */
        sprintf(buffer,
          "%s%c%04d%c%04d_%06d.xyz", ti->od.qmd.qmd_dir,
          SEPARATOR, ti->od.al.mol_info[object_num]->object_id,
          SEPARATOR, ti->od.al.mol_info[object_num]->object_id, run);
        remove(buffer);
        sprintf(buffer, "%s%c%04d%c%04d_%06d.001", ti->od.qmd.qmd_dir, SEPARATOR,
          ti->od.al.mol_info[object_num]->object_id,
          SEPARATOR, ti->od.al.mol_info[object_num]->object_id, run - 1);
        remove(buffer);
        sprintf(buffer, "%s%c%04d%c%04d_%06d.dyn", ti->od.qmd.qmd_dir, SEPARATOR,
          ti->od.al.mol_info[object_num]->object_id,
          SEPARATOR, ti->od.al.mol_info[object_num]->object_id, run - 1);
        remove(buffer);
      }
    }
    else {
      /*
      all files are already in place and well-formed,
      so no QMD runs have to be performed (the check
      "if (run > restart)" in the following will be
      always false, so no runs will be carried out)
      */
      restart = run + 1;
    }
    for (run = 0, n_conf = 0; (!(ti->od.al.task_list[object_num]->code))
      && (run <= ti->od.qmd.runs); ++run) {
      sprintf(buffer2, "%04d_%06d.xyz", ti->od.al.mol_info[object_num]->object_id, run);
      if (run > restart) {
        if (!run) {
          /*
          if this is the first run, the starting geometry has to be optimized
          */
          sprintf(buffer, "%04d.xyz", ti->od.al.mol_info[object_num]->object_id);
        }
        else {
          /*
          otherwise MD is carried out on the previous optimized geometry,
          then the last geometry of the MD trajectory is optimized
          */
          sprintf(buffer, "%04d_%06d.xyz", ti->od.al.mol_info[object_num]->object_id, run - 1);
          if ((ti->od.al.task_list[object_num]->code = tinker_dynamic
            (&(ti->od), work_dir, buffer, object_num, run,
            ti->od.mel.random_seed_array[run - 1]))) {
            continue;
          }
          sprintf(buffer, "%04d_%06d.001", ti->od.al.mol_info[object_num]->object_id, run - 1);
        }
        if ((ti->od.al.task_list[object_num]->code = tinker_minimize
          (&(ti->od), work_dir, buffer, buffer2, object_num, run))) {
          continue;
        }
        if (!run) {
          continue;
        }
        /*
        the raw MD geometry can be removed
        */
        sprintf(buffer, "%s%c%04d%c%04d_%06d.001", ti->od.qmd.qmd_dir, SEPARATOR,
          ti->od.al.mol_info[object_num]->object_id,
          SEPARATOR, ti->od.al.mol_info[object_num]->object_id, run - 1);
        remove(buffer);
        /*
        as well as the .dyn file
        */
        sprintf(buffer, "%s%c%04d%c%04d_%06d.dyn", ti->od.qmd.qmd_dir, SEPARATOR,
          ti->od.al.mol_info[object_num]->object_id,
          SEPARATOR, ti->od.al.mol_info[object_num]->object_id, run - 1);
        remove(buffer);
      }
      /*
      open the optimized geometry
      */
      sprintf(inp_fd.name,
        "%s%c%04d%c%04d_%06d.xyz", ti->od.qmd.qmd_dir,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id, run);
      if (!(inp_fd.handle = fopen(inp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_TEMP_FILE;
        continue;
      }
      j = 0;
      n = 0;
      while ((n < n_atoms) && fgets(buffer, BUF_LEN, inp_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        /*
        if this is the first line, read number of atoms and energy
        */
        if (!j) {
          read_tinker_xyz_n_atoms_energy(buffer, &n, &energy);
          j = 1;
          if (n != n_atoms) {
            n = 0;
            ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
            break;
          }
          n = 0;
        }
        else {
          /*
          otherwise read and store coordinates
          */
          sscanf(buffer, "%*s %*s %lf %lf %lf", &(conf[O3_CURR]->coord[n * 3]),
            &(conf[O3_CURR]->coord[n * 3 + 1]), &(conf[O3_CURR]->coord[n * 3 + 2]));
          ++n;
        }
      }
      fclose(inp_fd.handle);
      inp_fd.handle = NULL;
      if (n != n_atoms) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
        continue;
      }
      if (!(conf[O3_CURR]->h = (int **)alloc_array(conf[O3_CURR]->n_heavy_atoms, MAX_H_BINS * sizeof(int)))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      compute_conf_h(conf[O3_CURR]);
      min_heavy_msd = MAX_CUTOFF;
      min_pos = 0;
      user_msd = square(ti->od.qmd.rmsd);
      i = 0;
      /*
      both superposition techniques are attempted;
      if by either method the current conformation turns out to be
      similar enough as per the user-defined RMSD criterion to an
      already existing one, energies are compared, and the one
      having the lowest energy is kept; otherwise, the current
      conformation is added to conf_array
      */
      while (((min_heavy_msd - user_msd) > MSD_THRESHOLD) && (i < n_conf)) {
        superpose_conf_syst(conf_array[i], conf[O3_CURR], conf[O3_FITTED],
          conf[O3_PROGRESS], conf[O3_CAND], sdm[0], sdm[1], sdm[2], 
          used, rt_mat, ANGLE_STEP, &heavy_msd_syst, NULL, &pairs);
        if ((min_heavy_msd - heavy_msd_syst) > MSD_THRESHOLD) {
          min_heavy_msd = heavy_msd_syst;
          min_pos = i;
        }
        ++i;
      }
      i = 0;
      while (((min_heavy_msd - user_msd) > MSD_THRESHOLD) && (i < n_conf)) {
        superpose_conf_lap(&li, conf_array[i], conf[O3_CURR], conf[O3_FITTED],
          conf[O3_PROGRESS], sdm[0], sdm[1],
          used, rt_mat, &heavy_msd_lap, NULL, &pairs);
        if ((min_heavy_msd - heavy_msd_lap) > MSD_THRESHOLD) {
          min_heavy_msd = heavy_msd_lap;
          min_pos = i;
        }
        ++i;
      }
      n = -1;
      if ((min_heavy_msd - user_msd) > MSD_THRESHOLD) {
        /*
        this conformation is new: let's store it
        */
        if (!(conf_array[n_conf] = alloc_conf(n_atoms))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
          break;
        }
        n = n_conf;
        ++n_conf;
      }
      else if (n_conf && ((conf_array[min_pos]->energy - energy) > ALMOST_ZERO)) {
        /*
        this conformation is very similar to an existing one,
        but its energy is lower, so the old geometry is replaced
        */
        n = min_pos;
      }
      if (n != -1) {
        conf_array[n]->n_atoms = n_atoms;
        conf_array[n]->n_heavy_atoms = conf[O3_CURR]->n_heavy_atoms;
        conf_array[n]->atom = atom;
        if (conf_array[n]->h) {
          free_array(conf_array[n]->h);
        }
        conf_array[n]->h = conf[O3_CURR]->h;
        conf_array[n]->energy = energy;
        conf_array[n]->n_conf = n + 1;
        cblas_dcopy(n_atoms * 3, conf[O3_CURR]->coord, 1, conf_array[n]->coord, 1);
        /*
        conformations are sorted by increasing energy, then
        the array is truncated as soon as the user-defined
        threshold above the global minimum is exceeded
        */
        qsort(conf_array, n_conf, sizeof(ConfInfo *), compare_conf_energy);
        if (conf_array[0]) {
          glob_min = conf_array[0]->energy;
        }
        i = 1;
        while (i < n_conf) {
          energy = conf_array[i]->energy - glob_min;
          if (energy > ti->od.qmd.range) {
            break;
          }
          ++i;
        }
        if (i < n_conf) {
          for (j = i; j < n_conf; ++j) {
            free_array(conf_array[j]->h);
            conf_array[j]->h = NULL;
            free_conf(conf_array[j]);
            conf_array[j] = NULL;
          }
          n_conf = i;
        }
      }
      else if (conf[O3_CURR]->h) {
        free_array(conf[O3_CURR]->h);
        conf[O3_CURR]->h = NULL;
      }
    }
    if (ti->od.qmd.options & QMD_KEEP_INITIAL) {
      /*
      after all runs have been carried out, eventually
      the starting conformation is added to the conformational
      array, its energy is computed, and the array is sorted
      once again by increasing energy
      */
      if (!(conf_array[n_conf] = alloc_conf(n_atoms))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      if (!(conf_array[n_conf]->h = (int **)alloc_array
        (conf[O3_CURR]->n_heavy_atoms, MAX_H_BINS * sizeof(int)))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      compute_conf_h(conf_array[n_conf]);
      conf_array[n_conf]->n_atoms = n_atoms;
      conf_array[n_conf]->n_heavy_atoms = conf[O3_CURR]->n_heavy_atoms;
      conf_array[n_conf]->atom = atom;
      sprintf(buffer, "%04d.xyz", ti->od.al.mol_info[object_num]->object_id);
      ti->od.al.task_list[object_num]->code = tinker_analyze
        (&(ti->od), work_dir, buffer, object_num, -1);
      if (ti->od.al.task_list[object_num]->code) {
        continue;
      }
      sprintf(inp_fd.name, "%s%c%04d%c%04d.xyz", ti->od.qmd.qmd_dir,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
      if (!(inp_fd.handle = fopen(inp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_TEMP_FILE;
        continue;
      }
      n = 0;
      if (fgets(buffer, BUF_LEN, inp_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        read_tinker_xyz_n_atoms_energy(buffer, &n, &energy);
      }
      if (n != n_atoms) {
        energy = 0.0;
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
      }
      fclose(inp_fd.handle);
      inp_fd.handle = NULL;
      if (ti->od.al.task_list[object_num]->code) {
        continue;
      }
      conf_array[n_conf]->energy = energy;
      conf_array[n_conf]->n_conf = n_conf + 1;
      for (i = 0; i < n_atoms; ++i) {
        cblas_dcopy(3, atom[i]->coord, 1, &(conf_array[n_conf]->coord[i * 3]), 1);
      }
      ++n_conf;
      qsort(conf_array, n_conf, sizeof(ConfInfo *), compare_conf_energy);
    }
    if (!(ti->od.al.task_list[object_num]->code)) {
      /*
      finally, each conformation is aligned onto the global minimum
      using the method which yields the lowest RMSD
      */
      for (i = 1; i < n_conf; ++i) {
        superpose_conf_syst(conf_array[i], conf_array[0], conf[O3_FITTED_SYST],
          conf[O3_PROGRESS], conf[O3_CAND], sdm[0], sdm[1],
          sdm[2], used, rt_mat, ANGLE_STEP, &heavy_msd_syst, NULL, &pairs);
        superpose_conf_lap(&li, conf_array[i], conf_array[0], conf[O3_FITTED_LAP],
          conf[O3_PROGRESS], sdm[0], sdm[1],
          used, rt_mat, &heavy_msd_lap, NULL, &pairs);
        fitted_conf = (((heavy_msd_syst - heavy_msd_lap) > MSD_THRESHOLD)
          ? conf[O3_FITTED_LAP] : conf[O3_FITTED_SYST]);
        cblas_dcopy(n_atoms * 3, fitted_conf->coord, 1, conf_array[i]->coord, 1);
      }
      sprintf(mol_fd.name, "%s%c%04d.mol", ti->od.field.mol_dir,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
      if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_MOL_FILE;
        continue;
      }
      sprintf(inp_fd.name, "%s%c%04d.sdf", ti->od.qmd.qmd_dir, SEPARATOR,
        ti->od.al.mol_info[object_num]->object_id);
      if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        fclose(mol_fd.handle);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
        continue;
      }
      if (conf_array[0]) {
        glob_min = conf_array[0]->energy;
      }
      for (i = 0; (!(ti->od.al.task_list[object_num]->code)) && (i < n_conf); ++i) {
        /*
        SDF header is written
        */
        fprintf(inp_fd.handle,
          "%04d_%06d\n"
          "  %-18s3D\n"
          "%04d_%06d\n",
          ti->od.al.mol_info[object_num]->object_id, i + 1,
          PACKAGE_NAME, ti->od.al.mol_info[object_num]->object_id, i + 1);
        if (find_conformation_in_sdf(mol_fd.handle, inp_fd.handle, SKIP_MOL_NAME)) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
          ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_MOL_FILE;
          continue;
        }
        atom_num = 0;
        while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          if (atom_num < conf[O3_FITTED]->n_atoms) {
            /*
            coordinates from the original MOL file
            are replaced by those of the current conformation
            */
            if (replace_coord(ti->od.al.mol_info[object_num]->sdf_version,
              buffer, &(conf_array[i]->coord[atom_num * 3]))) {
              break;
            }
            ++atom_num;
          }
          fprintf(inp_fd.handle, "%s\n", buffer);
        }
        if (atom_num < n_atoms) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
          ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_MOL_FILE;
          continue;
        }
        /*
        absolute and relative energies are written
        */
        fprintf(inp_fd.handle,
          ">  <MMFF94S %sENERGY>\n"
          "%.4lf\n\n"
          ">  <DELTA>\n"
          "%.4lf\n\n"
          SDF_DELIMITER "\n", ((ti->od.qmd.options & QMD_GBSA) ? "GBSA " : ""),
          conf_array[i]->energy, conf_array[i]->energy - glob_min);
        rewind(mol_fd.handle);
      }
      fclose(inp_fd.handle);
      fclose(mol_fd.handle);
    }
    i = 0;
    while (conf_array[i]) {
      free_array(conf_array[i]->h);
      conf_array[i]->h = NULL;
      free_conf(conf_array[i]);
      conf_array[i] = NULL;
      ++i;
    }
    if ((!(ti->od.al.task_list[object_num]->code)) && (ti->od.qmd.options & QMD_REMOVE_FOLDER)) {
      /*
      remove folder with intermediate XYZ files if the user wants so
      */
      sprintf(buffer, "%s%c%04d", ti->od.qmd.qmd_dir, SEPARATOR,
        ti->od.al.mol_info[object_num]->object_id);
      remove_recursive(buffer);
    }
  }
  free_array(atom);
  free_array(bond_list);
  free_lap_info(&li);
  for (i = 0; i < O3_MAX_SDM; ++i) {
    if (sdm[i]) {
      free(sdm[i]);
    }
  }
  for (i = 0; i < 2; ++i) {
    if (used[i]) {
      free(used[i]);
    }
  }
  if (conf_array) {
    free(conf_array);
  }
  for (i = 0; i < O3_MAX_CONF; ++i) {
    free_conf(conf[i]);
  }

  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}  


int energy(O3Data *od)
{
  int i;
  int result = 0;
  int n_threads;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  ThreadInfo **ti;


  if (alloc_threads(od)) {
    return OUT_OF_MEMORY;
  }
  if (!(od->al.task_list = (TaskInfo **)
    alloc_array(od->grid.object_num, sizeof(TaskInfo)))) {
    return OUT_OF_MEMORY;
  }
  for (i = 0; i < od->grid.object_num ; ++i) {
    od->al.mol_info[i]->done = 0;
  }
  #ifndef WIN32
  pthread_mutex_init(od->mel.mutex, NULL);
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
  #else
  if (!(*(od->mel.mutex) = CreateMutex(NULL, FALSE, NULL))) {
    return CANNOT_CREATE_THREAD;
  }
  #endif  
  n_threads = fill_thread_info(od, od->grid.object_num);
  ti = od->mel.thread_info;
  for (i = 0; i < n_threads; ++i) {
    /*
    create the i-th thread
    */
    #ifndef WIN32
    od->error_code = pthread_create(&(od->thread_id[i]),
      &thread_attr, (void *(*)(void *))energy_thread, ti[i]);
    if (od->error_code) {
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)energy_thread,
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
  pthread_mutex_destroy(od->mel.mutex);
  #else
  WaitForMultipleObjects(n_threads, od->hThreadArray, TRUE, INFINITE);
  for (i = 0; i < n_threads; ++i) {
    CloseHandle(od->hThreadArray[i]);
  }
  CloseHandle(*(od->mel.mutex));
  #endif
  if (!(od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT)) {
    for (i = 0, result = 1; result && (i < od->grid.object_num); ++i) {
      sprintf(od->qmd.src, "%s%c%04d.sdf", od->align.align_scratch,
        SEPARATOR, od->al.mol_info[i]->object_id);
      result = fcopy(od->qmd.src, od->qmd.dest, (i ? "ab" : "wb"));
    }
    result = (result ? 0 : ERROR_MERGING_FILES);
  }
  
  return result;
}


#ifndef WIN32
void *energy_thread(void *pointer)
#else
DWORD energy_thread(void *pointer)
#endif
{
  char buffer[BUF_LEN];
  char buffer2[BUF_LEN];
  char *used[2] = { NULL, NULL };
  int i;
  int j;
  int n;
  int object_num;
  int atom_num = 0;
  int conf_num = 0;
  int conf_max = 1;
  int n_conf;
  int n_atoms = 0;
  int alloc_fail = 0;
  int alloc_new_conf = 0;
  int assigned = 1;
  int minimize = 0;
  int min_pos = 0;
  int pairs = 0;
  double heavy_msd_lap = 0.0;
  double heavy_msd_syst = 0.0;
  double original_heavy_msd_lap = 0.0;
  double original_heavy_msd_syst = 0.0;
  double min_heavy_msd = 0.0;
  double user_msd;
  double glob_min = 0.0;
  double energy = 0.0;
  double rt_mat[RT_MAT_SIZE];
  LAPInfo li;
  AtomPair *sdm[O3_MAX_SDM] = { NULL, NULL, NULL };
  AtomInfo **atom = NULL;
  BondList **bond_list = NULL;
  ConfInfo *conf[O3_MAX_CONF] = { NULL, NULL, NULL, NULL };
  ConfInfo **conf_array = NULL;
  ConfInfo *fitted_conf = NULL;
  ThreadInfo *ti;
  FileDescriptor sdf_fd;
  FileDescriptor mol_fd;


  ti = (ThreadInfo *)pointer;
  memset(&sdf_fd, 0, sizeof(FileDescriptor));
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  memset(buffer2, 0, BUF_LEN);
  minimize = (strcmp(ti->od.qmd.minimizer, TINKER_ANALYZE_EXE) ? 1 : 0);
  /*
  allocate memory for AtomInfo structure array
  */
  if (!(atom = (AtomInfo **)alloc_array(ti->od.field.max_n_atoms + 1, sizeof(AtomInfo)))) {
    alloc_fail = 1;
  }
  if (!(bond_list = (BondList **)alloc_array(ti->od.field.max_n_bonds + 1, sizeof(BondList)))) {
    alloc_fail = 1;
  }
  for (i = 0; i < O3_MAX_CONF; ++i) {
    if (!(conf[i] = alloc_conf(ti->od.field.max_n_atoms))) {
      alloc_fail = 1;
    }
  }
  if ((ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT)
    && (ti->od.qmd.options & (QMD_ALIGN | QMD_REMOVE_DUPLICATES))) {
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
    }
  }
  while (assigned) {
    object_num = 0;
    assigned = 0;
    while ((!assigned) && (object_num < ti->od.grid.object_num)) {
      if (!(ti->od.al.mol_info[object_num]->done)) {
        #ifndef WIN32
        pthread_mutex_lock(ti->od.mel.mutex);
        #else
        WaitForSingleObject(ti->od.mel.mutex, INFINITE);
        #endif
        if (!(ti->od.al.mol_info[object_num]->done)) {
          ti->od.al.mol_info[object_num]->done = 1;
          assigned = 1;
        }
        #ifndef WIN32
        pthread_mutex_unlock(ti->od.mel.mutex);
        #else
        ReleaseMutex(ti->od.mel.mutex);
        #endif
      }
      if (!assigned) {
        ++object_num;
      }
    }
    if (!assigned)  {
      break;
    }
    if (alloc_fail) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
      continue;
    }
    if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
      sprintf(sdf_fd.name, "%s%c%04d.sdf", ti->od.qmd.src,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
      /*
      if this object is missing, no harm
      */
      if (!(sdf_fd.handle = fopen(sdf_fd.name, "rb"))) {
        continue;
      }
    }
    if ((ti->od.al.task_list[object_num]->code =
      fill_atom_info(&(ti->od), ti->od.al.task_list[object_num],
        atom, bond_list, object_num, O3_MMFF94))) {
      continue;
    }
    if ((ti->od.al.task_list[object_num]->code = fill_tinker_types(atom))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      continue;
    }
    n_atoms = ti->od.al.mol_info[object_num]->n_atoms;
    /*
    prepare the .bnd file
    */
    sprintf(buffer, "%s%c%04d.bnd", ti->od.align.align_scratch, SEPARATOR,
      ti->od.al.mol_info[object_num]->object_id);
    ti->od.al.task_list[object_num]->code = write_tinker_xyz_bnd
      (&(ti->od), atom, bond_list, n_atoms, object_num, NULL, buffer);
    if (ti->od.al.task_list[object_num]->code) {
      continue;
    }
    for (i = 0; i < O3_MAX_CONF; ++i) {
      conf[i]->atom = atom;
      conf[i]->n_atoms = n_atoms;
      conf[i]->n_heavy_atoms = ti->od.al.mol_info[object_num]->n_heavy_atoms;
    }
    if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
      conf_max = ti->od.pel.conf_population[ANY_DB]->pe[object_num];
      if ((conf_array = (ConfInfo **)realloc(conf_array,
        ((conf_max + 1) * sizeof(ConfInfo *))))) {
        memset(conf_array, 0, (conf_max + 1) * sizeof(ConfInfo *));
      }
    }
    for (conf_num = 0, n_conf = 0; conf_num < conf_max; ++conf_num) {
      if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
        if (find_conformation_in_sdf(sdf_fd.handle, NULL, 0)) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[object_num], sdf_fd.name);
          ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_SDF_FILE;
          continue;
        }
        atom_num = 0;
        while ((atom_num < n_atoms)
          && fgets(buffer, BUF_LEN, sdf_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          /*
          read coordinates of the current conformation
          */
          sscanf(buffer, "%lf %lf %lf", &(atom[atom_num]->coord[0]),
            &(atom[atom_num]->coord[1]), &(atom[atom_num]->coord[2]));
          ++atom_num;
        }
        if (atom_num < n_atoms) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[object_num], sdf_fd.name);
          ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_SDF_FILE;
          continue;
        }
        while (fgets(buffer, BUF_LEN, sdf_fd.handle)
          && strncmp(buffer, SDF_DELIMITER, 4));
      }
      /*
      prepare the XYZ geometry
      */
      sprintf(buffer, "%s%c%04d_%06d.xyz", ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[object_num]->object_id, conf_num + 1);
      ti->od.al.task_list[object_num]->code = write_tinker_xyz_bnd
        (&(ti->od), atom, bond_list, n_atoms, object_num, buffer, NULL);
      if (ti->od.al.task_list[object_num]->code) {
        continue;
      }
      /*
      call TINKER tool
      */
      if (minimize) {
        sprintf(buffer, "%04d_%06d.xyz",
          ti->od.al.mol_info[object_num]->object_id, conf_num + 1);
        sprintf(buffer2, "%04d_%06d_min.xyz",
          ti->od.al.mol_info[object_num]->object_id, conf_num + 1);
      }
      else {
        sprintf(buffer2, "%04d_%06d.xyz",
          ti->od.al.mol_info[object_num]->object_id, conf_num + 1);
      }
      ti->od.al.task_list[object_num]->code = (minimize
        ? tinker_minimize(&(ti->od), ti->od.align.align_scratch,
        buffer, mol_fd.name, object_num, conf_num + 1)
        : tinker_analyze(&(ti->od), ti->od.align.align_scratch,
        buffer2, object_num, conf_num + 1));
      if (ti->od.al.task_list[object_num]->code) {
        continue;
      }
      /*
      open the optimized (or analyzed) geometry
      */
      sprintf(mol_fd.name, "%s%c%s", ti->od.align.align_scratch, SEPARATOR, buffer2);
      if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_TEMP_FILE;
        continue;
      }
      j = 0;
      n = 0;
      while ((n < n_atoms) && fgets(buffer, BUF_LEN, mol_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        /*
        if this is the first line, read number of atoms and energy
        */
        if (!j) {
          read_tinker_xyz_n_atoms_energy(buffer, &n, &energy);
          j = 1;
          if (n != n_atoms) {
            n = 0;
            ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
            break;
          }
          n = 0;
        }
        else {
          /*
          otherwise read and store coordinates
          */
          sscanf(buffer, "%*s %*s %lf %lf %lf", &(conf[O3_CURR]->coord[n * 3]),
            &(conf[O3_CURR]->coord[n * 3 + 1]), &(conf[O3_CURR]->coord[n * 3 + 2]));
          ++n;
        }
      }
      fclose(mol_fd.handle);
      remove(mol_fd.name);
      if (minimize) {
        sprintf(mol_fd.name, "%s%c%s", ti->od.align.align_scratch, SEPARATOR, buffer2);
        remove(mol_fd.name);
      }
      if (n != n_atoms) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
        continue;
      }
      if (!(ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT)) {
        continue;
      }
      if (ti->od.qmd.options & (QMD_ALIGN | QMD_REMOVE_DUPLICATES)) {
        if (!(conf[O3_CURR]->h = (int **)alloc_array
          (conf[O3_CURR]->n_heavy_atoms, MAX_H_BINS * sizeof(int)))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
          continue;
        }
        compute_conf_h(conf[O3_CURR]);
      }
      n = n_conf;
      alloc_new_conf = 1;
      if (ti->od.qmd.options & QMD_REMOVE_DUPLICATES) {
        min_heavy_msd = MAX_CUTOFF;
        min_pos = 0;
        user_msd = square(ti->od.qmd.rmsd);
        i = 0;
        /*
        both superposition techniques are attempted;
        if by either method the current conformation turns out to be
        similar enough as per the user-defined RMSD criterion to an
        already existing one, energies are compared, and the one
        having the lowest energy is kept; otherwise, the current
        conformation is added to conf_array
        */
        while (((min_heavy_msd - user_msd) > MSD_THRESHOLD) && (i < n_conf)) {
          superpose_conf_syst(conf_array[i], conf[O3_CURR], conf[O3_FITTED],
            conf[O3_PROGRESS], conf[O3_CAND], sdm[0], sdm[1], sdm[2], 
            used, rt_mat, ANGLE_STEP, &heavy_msd_syst,
            &original_heavy_msd_syst, &pairs);
          if (ti->od.qmd.options & QMD_DONT_SUPERPOSE) {
            heavy_msd_syst = original_heavy_msd_syst;
          }
          if ((min_heavy_msd - heavy_msd_syst) > MSD_THRESHOLD) {
            min_heavy_msd = heavy_msd_syst;
            min_pos = i;
          }
          ++i;
        }
        i = 0;
        while (((min_heavy_msd - user_msd) > MSD_THRESHOLD) && (i < n_conf)) {
          superpose_conf_lap(&li, conf_array[i], conf[O3_CURR], conf[O3_FITTED],
            conf[O3_PROGRESS], sdm[0], sdm[1],
            used, rt_mat, &heavy_msd_lap,
            &original_heavy_msd_lap, &pairs);
          if (ti->od.qmd.options & QMD_DONT_SUPERPOSE) {
            heavy_msd_lap = original_heavy_msd_lap;
          }
          if ((min_heavy_msd - heavy_msd_lap) > MSD_THRESHOLD) {
            min_heavy_msd = heavy_msd_lap;
            min_pos = i;
          }
          ++i;
        }
        n = -1;
        alloc_new_conf = 0;
        if ((min_heavy_msd - user_msd) > MSD_THRESHOLD) {
          /*
          this conformation is new: let's store it
          */
          alloc_new_conf = 1;
          n = n_conf;
        }
        else if (n_conf && ((conf_array[min_pos]->energy - energy) > ALMOST_ZERO)) {
          /*
          this conformation is very similar to an existing one,
          but its energy is lower, so the old geometry is replaced
          */
          n = min_pos;
        }
      }
      if (alloc_new_conf) {
        if (!(conf_array[n_conf] = alloc_conf(n_atoms))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
          break;
        }
        ++n_conf;
      }
      if (n != -1) {
        conf_array[n]->n_atoms = n_atoms;
        conf_array[n]->n_heavy_atoms = conf[O3_CURR]->n_heavy_atoms;
        conf_array[n]->atom = atom;
        if (conf_array[n]->h) {
          free_array(conf_array[n]->h);
        }
        conf_array[n]->h = conf[O3_CURR]->h;
        conf_array[n]->energy = energy;
        conf_array[n]->n_conf = n + 1;
        cblas_dcopy(n_atoms * 3, conf[O3_CURR]->coord, 1, conf_array[n]->coord, 1);
        /*
        conformations are sorted by increasing energy, then
        the array is truncated as soon as the user-defined
        threshold above the global minimum is exceeded
        */
        qsort(conf_array, n_conf, sizeof(ConfInfo *), compare_conf_energy);
        if (conf_array[0]) {
          glob_min = conf_array[0]->energy;
        }
        i = 1;
        while (i < n_conf) {
          energy = conf_array[i]->energy - glob_min;
          if (energy > ti->od.qmd.range) {
            break;
          }
          ++i;
        }
        if (i < n_conf) {
          for (j = i; j < n_conf; ++j) {
            free_array(conf_array[j]->h);
            conf_array[j]->h = NULL;
            free_conf(conf_array[j]);
            conf_array[j] = NULL;
          }
          n_conf = i;
        }
      }
      else if (conf[O3_CURR]->h) {
        free_array(conf[O3_CURR]->h);
        conf[O3_CURR]->h = NULL;
      }
    }
    if (sdf_fd.handle) {
      fclose(sdf_fd.handle);
      sdf_fd.handle = NULL;
      conf_max = n_conf;
    }
    if (!(ti->od.al.task_list[object_num]->code)) {
      if ((ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT)
        && (ti->od.qmd.options & QMD_ALIGN)) {
        /*
        finally, if the user wants so,
        each conformation is aligned onto the global minimum
        using the method which yields the lowest RMSD
        */
        for (i = 1; i < n_conf; ++i) {
          superpose_conf_syst(conf_array[i], conf_array[0], conf[O3_FITTED_SYST],
            conf[O3_PROGRESS], conf[O3_CAND], sdm[0], sdm[1],
            sdm[2], used, rt_mat, ANGLE_STEP, &heavy_msd_syst, NULL, &pairs);
          superpose_conf_lap(&li, conf_array[i], conf_array[0], conf[O3_FITTED_LAP],
            conf[O3_PROGRESS], sdm[0], sdm[1],
            used, rt_mat, &heavy_msd_lap, NULL, &pairs);
          fitted_conf = (((heavy_msd_syst - heavy_msd_lap) > MSD_THRESHOLD)
            ? conf[O3_FITTED_LAP] : conf[O3_FITTED_SYST]);
          cblas_dcopy(n_atoms * 3, fitted_conf->coord, 1, conf_array[i]->coord, 1);
        }
      }
      sprintf(mol_fd.name, "%s%c%04d.mol", ti->od.field.mol_dir,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
      if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_MOL_FILE;
        continue;
      }
      sprintf(sdf_fd.name, "%s%c%04d.sdf",
        (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT)
        ? ti->od.qmd.dest : ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[object_num]->object_id);
      if (!(sdf_fd.handle = fopen(sdf_fd.name, "wb+"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], sdf_fd.name);
        fclose(mol_fd.handle);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
        continue;
      }
      if ((ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) && conf_array[0]) {
        glob_min = conf_array[0]->energy;
      }
      for (i = 0; (!(ti->od.al.task_list[object_num]->code)) && (i < conf_max); ++i) {
        /*
        SDF header is written
        */
        if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
          fprintf(sdf_fd.handle,
            "%04d_%06d\n"
            "  %-18s3D\n"
            "%04d_%06d\n",
            ti->od.al.mol_info[object_num]->object_id, i + 1,
            PACKAGE_NAME, ti->od.al.mol_info[object_num]->object_id, i + 1);
        }
        if (find_conformation_in_sdf(mol_fd.handle, sdf_fd.handle,
          (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) ? SKIP_MOL_NAME : 0)) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
          ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_MOL_FILE;
          continue;
        }
        atom_num = 0;
        while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          if (atom_num < conf[O3_FITTED]->n_atoms) {
            /*
            coordinates from the original MOL file
            are replaced by those of the current conformation
            */
            if (replace_coord(ti->od.al.mol_info[object_num]->sdf_version,
              buffer, &(((ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT)
              ? conf_array[i] : conf[O3_CURR])->coord[atom_num * 3]))) {
              break;
            }
            ++atom_num;
          }
          fprintf(sdf_fd.handle, "%s\n", buffer);
        }
        if (atom_num < n_atoms) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[object_num], mol_fd.name);
          ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_MOL_FILE;
          continue;
        }
        /*
        absolute and relative energies are written
        */
        fprintf(sdf_fd.handle,
          ">  <MMFF94S %sENERGY>\n"
          "%.4lf\n\n"
          ">  <DELTA>\n"
          "%.4lf\n\n"
          SDF_DELIMITER "\n", ((ti->od.qmd.options & QMD_GBSA) ? "GBSA " : ""),
          conf_array[i]->energy, conf_array[i]->energy - glob_min);
        rewind(mol_fd.handle);
      }
      fclose(sdf_fd.handle);
      fclose(mol_fd.handle);
    }
    i = 0;
    while ((ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT)
      && conf_array[i]) {
      if (ti->od.qmd.options & (QMD_ALIGN | QMD_REMOVE_DUPLICATES)) {
        free_array(conf_array[i]->h);
        conf_array[i]->h = NULL;
      }
      free_conf(conf_array[i]);
      conf_array[i] = NULL;
      ++i;
    }
  }
  free_array(atom);
  free_array(bond_list);
  if (conf_array) {
    free(conf_array);
  }
  for (i = 0; i < O3_MAX_CONF; ++i) {
    free_conf(conf[i]);
  }
  if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
    if (ti->od.qmd.options & (QMD_ALIGN | QMD_REMOVE_DUPLICATES)) {
      free_lap_info(&li);
      for (i = 0; i < O3_MAX_SDM; ++i) {
        if (sdm[i]) {
          free(sdm[i]);
        }
      }
      for (i = 0; i < 2; ++i) {
        if (used[i]) {
          free(used[i]);
        }
      }
    }
  }

  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}  
