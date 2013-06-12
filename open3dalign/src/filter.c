/*

filter.c

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
#include <include/prog_exe_info.h>
#include <include/proc_env.h>
#ifdef WIN32
#include <windows.h>
#endif


int filter(O3Data *od)
{
  char buffer[BUF_LEN];
  int i;
  int n_conf;
  int n_conf_overall;
  int found;
  int last_checked_found;
  int n_retained_conf;
  int deleted_object_num;
  int deleted_conf_num;
  int last_object_num;
  int last_conf_num;
  int template_num;
  int template_object_num;
  int old_template_object_num;
  int template_conf_num;
  int n_threads;
  FileDescriptor temp_fd;
  FileDescriptor sdf_fd;
  FileDescriptor filter_log_fd;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  ThreadInfo **ti;


  ti = od->mel.thread_info;
  memset(buffer, 0, BUF_LEN);
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&sdf_fd, 0, sizeof(FileDescriptor));
  memset(&filter_log_fd, 0, sizeof(FileDescriptor));
  for (template_num = 0, i = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
    template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    sprintf(buffer, "%s%c%04d.sdf", od->align.filter_conf_dir,
      SEPARATOR, od->al.mol_info[template_object_num]->object_id);
    od->al.mol_info[template_object_num]->done = 0;
    if (fexist(buffer)) {
      if (od->align.filter_type & FILTER_INTRA_CONF_DB_BIT) {
        od->al.mol_info[template_object_num]->done = OBJECT_ALREADY_DONE;
      }
      ++i;
    }
  }
  if (i == od->pel.numberlist[OBJECT_LIST]->size) {
    return NOTHING_TO_DO_FILTER;
  }
  if (alloc_threads(od)) {
    return OUT_OF_MEMORY;
  }
  od->align.n_tasks = od->pel.numberlist[OBJECT_LIST]->size;
  if (!(od->al.task_list = (TaskInfo **)alloc_array(od->align.n_tasks, sizeof(TaskInfo)))) {
    return OUT_OF_MEMORY;
  }
  /*
  count how many conformations we have overall
  */
  for (template_num = 0, n_conf_overall = 0;
    template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
    template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    n_conf_overall += od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num];
  }
  if (!(od->al.phar_conf_list = (PharConfInfo **)alloc_array(n_conf_overall, sizeof(PharConfInfo)))) {
    return OUT_OF_MEMORY;
  }
  for (template_num = 0, n_conf_overall = 0;
    template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
    template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    for (template_conf_num = 0; template_conf_num < od->pel.conf_population
      [TEMPLATE_DB]->pe[template_object_num]; ++template_conf_num) {
      od->al.phar_conf_list[n_conf_overall]->object_num = template_object_num;
      od->al.phar_conf_list[n_conf_overall]->conf_num = template_conf_num;
      ++n_conf_overall;
    }
  }
  n_threads = fill_thread_info(od, od->align.n_tasks);
  for (i = 0; i < od->pel.numberlist[OBJECT_LIST]->size; ++i) {
    od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM] =
      od->pel.numberlist[OBJECT_LIST]->pe[i] - 1;
    od->al.task_list[i]->data[TEMPLATE_CONF_NUM] = -1;
    od->al.task_list[i]->data[MOVED_OBJECT_NUM] = -1;
    od->al.task_list[i]->data[MOVED_CONF_NUM] = -1;
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
  for (i = 0; i < n_threads; ++i) {
    /*
    create the i-th thread
    */
    #ifndef WIN32
    od->error_code = pthread_create(&(od->thread_id[i]),
      &thread_attr, (void *(*)(void *))filter_extract_split_phar_thread, ti[i]);
    if (od->error_code) {
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)filter_extract_split_phar_thread,
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
  for (i = 0; (i < od->align.n_tasks)
    && (!(od->al.task_list[i]->code)); ++i);
  if (i != od->align.n_tasks) {
    return ERROR_IN_FILTER_EXTRACT_SPLIT;
  }
  /*
  sort pharmacophores in decreasing n_phar_points order
  */
  qsort(od->al.phar_conf_list, n_conf_overall, sizeof(PharConfInfo *), compare_n_phar_points);
  for (template_num = 0, i = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
    template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    sprintf(buffer, "%s%c%04d.sdf", od->align.filter_conf_dir,
      SEPARATOR, od->al.mol_info[template_object_num]->object_id);
    od->al.mol_info[template_object_num]->done = 0;
    if (fexist(buffer)) {
      if (od->align.filter_type & FILTER_INTRA_CONF_DB_BIT) {
        od->al.mol_info[template_object_num]->done = OBJECT_ALREADY_DONE;
      }
      ++i;
    }
  }
  if (od->align.filter_type & FILTER_INTRA_CONF_DB_BIT) {
    #ifndef WIN32
    pthread_mutex_init(od->mel.mutex, NULL);
    pthread_attr_init(&thread_attr);
    pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
    #else
    if (!(*(od->mel.mutex) = CreateMutex(NULL, FALSE, NULL))) {
      return CANNOT_CREATE_THREAD;
    }
    #endif
    for (i = 0; i < n_threads; ++i) {
      /*
      create the i-th thread
      */
      #ifndef WIN32
      od->error_code = pthread_create(&(od->thread_id[i]),
        &thread_attr, (void *(*)(void *))filter_intra_thread, ti[i]);
      if (od->error_code) {
        return CANNOT_CREATE_THREAD;
      }
      #else
      od->hThreadArray[i] = CreateThread(NULL, 0,
        (LPTHREAD_START_ROUTINE)filter_intra_thread,
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
    for (i = 0; (i < od->align.n_tasks)
      && (!(od->al.task_list[i]->code)); ++i);
    if (i != od->align.n_tasks) {
      return ERROR_IN_FILTER_INTRA;
    }
  }
  else if (od->align.filter_type & FILTER_INTER_CONF_DB_BIT) {
    n_conf = 0;
    last_checked_found = 0;
    i = 0;
    sprintf(filter_log_fd.name, "%s%cfilter_inter.log", od->align.filter_conf_dir, SEPARATOR);
    if ((filter_log_fd.handle = fopen(filter_log_fd.name, "rb"))) {
      if (fgrep(filter_log_fd.handle, buffer, "LAST_CHECKED")) {
        rewind(filter_log_fd.handle);
        while ((!last_checked_found) && (i < n_conf_overall) && (n_conf < n_conf_overall)
          && fgets(buffer, BUF_LEN, filter_log_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          if (strncmp(buffer, "LAST_CHECKED", 12)) {
            sscanf(buffer, "%d %d", &deleted_object_num, &deleted_conf_num);
            i = 0;
            found = 0;
            while ((i < n_conf_overall)
              && (!(found = (od->al.phar_conf_list[i]->object_num == deleted_object_num)
              && (od->al.phar_conf_list[i]->conf_num == deleted_conf_num)))) {
              ++i;
            }
            if (found) {
              found = 0;
              od->al.phar_conf_list[i]->delete = 1;
            }
          }
          else {
            sscanf(&buffer[13], "%d %d", &last_object_num, &last_conf_num);
            n_conf = 0;
            last_checked_found = 0;
            while ((n_conf < n_conf_overall)
              && (!(last_checked_found = ((od->al.phar_conf_list[n_conf]->object_num == last_object_num)
              && (od->al.phar_conf_list[n_conf]->conf_num == last_conf_num))))) {
              ++n_conf;
            }
            ++n_conf;
          }
        }
        if (!last_checked_found) {
          /*
          the log file cannot be trusted; start from scratch
          */
          n_conf = 0;
          for (i = 0; i < n_conf_overall; ++i) {
            od->al.phar_conf_list[i]->delete = 0;
          }
        }
      }
      fclose(filter_log_fd.handle);
    }
    if (od->al.task_list) {
      free_array(od->al.task_list);
    }
    if (!(od->al.task_list = (TaskInfo **)alloc_array(od->n_proc, sizeof(TaskInfo)))) {
      return OUT_OF_MEMORY;
    }
    while (n_conf < n_conf_overall) {
      if (od->al.phar_conf_list[n_conf]->delete) {
        ++n_conf;
        continue;
      }
      for (i = 0, n_retained_conf = 0; i < n_conf_overall; ++i) {
        if (i == n_conf) {
          continue;
        }
        if (!(od->al.phar_conf_list[i]->delete)) {
          ++n_retained_conf;
        }
      }
      n_threads = fill_thread_info(od, n_retained_conf);
      od->align.n_tasks = n_threads;
      #ifndef WIN32
      pthread_attr_init(&thread_attr);
      pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
      #endif
      for (i = 0; i < n_threads; ++i) {
        ti[i]->data[DATA_N_CONF] = n_conf;
        ti[i]->data[DATA_N_CONF_OVERALL] = n_conf_overall;
        /*
        create the i-th thread
        */
        #ifndef WIN32
        od->error_code = pthread_create(&(od->thread_id[i]),
          &thread_attr, (void *(*)(void *))filter_inter_thread, ti[i]);
        if (od->error_code) {
          return CANNOT_CREATE_THREAD;
        }
        #else
        od->hThreadArray[i] = CreateThread(NULL, 0,
          (LPTHREAD_START_ROUTINE)filter_inter_thread,
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
      for (i = 0; (i < od->align.n_tasks)
        && (!(od->al.task_list[i]->code)); ++i);
      if (i != od->align.n_tasks) {
        return ERROR_IN_FILTER_INTER;
      }
      sprintf(filter_log_fd.name, "%s%cfilter_inter_new.log", od->align.filter_conf_dir, SEPARATOR);
      if (!(filter_log_fd.handle = fopen(filter_log_fd.name, "wb"))) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), filter_log_fd.name);
        return CANNOT_WRITE_TEMP_FILE;
      }
      for (i = 0; i < n_conf_overall; ++i) {
        if (od->al.phar_conf_list[i]->delete) {
          fprintf(filter_log_fd.handle, "%d\t%d\n",
            od->al.phar_conf_list[i]->object_num,
            od->al.phar_conf_list[i]->conf_num);
        }
      }
      fprintf(filter_log_fd.handle, "LAST_CHECKED\t%d\t%d\n",
        od->al.phar_conf_list[n_conf]->object_num,
        od->al.phar_conf_list[n_conf]->conf_num);
      fclose(filter_log_fd.handle);
      sprintf(buffer, "%s%cfilter_inter.log", od->align.filter_conf_dir, SEPARATOR);
      rename(filter_log_fd.name, buffer);
      ++n_conf;
    }
    free_array(od->al.task_list);
    od->al.task_list = NULL;
    old_template_object_num = -1;
    temp_fd.handle = NULL;
    sdf_fd.handle = NULL;
    for (i = 0; i < n_conf_overall; ++i) {
      template_object_num = od->al.phar_conf_list[i]->object_num;
      if (!(od->al.phar_conf_list[i]->delete)) {
        if (template_object_num != old_template_object_num) {
          if (temp_fd.handle) {
            fclose(temp_fd.handle);
            temp_fd.handle = NULL;
          }
          if (sdf_fd.handle) {
            fclose(sdf_fd.handle);
            sdf_fd.handle = NULL;
          }
          old_template_object_num = template_object_num;
          sprintf(temp_fd.name, "%s%c%04d.sdf", od->align.template_conf_dir,
            SEPARATOR, od->al.mol_info[template_object_num]->object_id);
          if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
            O3_ERROR_LOCATE(&(od->task));
            O3_ERROR_STRING(&(od->task), temp_fd.name);
            return FL_CANNOT_READ_SDF_FILE;
          }
          sprintf(sdf_fd.name, "%s%c%04d.sdf", od->align.filter_conf_dir,
            SEPARATOR, od->al.mol_info[template_object_num]->object_id);
          if (!(sdf_fd.handle = fopen(sdf_fd.name, "wb"))) {
            O3_ERROR_LOCATE(&(od->task));
            O3_ERROR_STRING(&(od->task), sdf_fd.name);
            return FL_CANNOT_WRITE_SDF_FILE;
          }
        }
        template_conf_num = od->al.phar_conf_list[i]->conf_num;
        if (find_conformation_in_sdf(temp_fd.handle, sdf_fd.handle, template_conf_num)) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), temp_fd.name);
          return FL_CANNOT_READ_SDF_FILE;
        }
        while (fgets(buffer, BUF_LEN, temp_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          fprintf(sdf_fd.handle, "%s\n", buffer);
          if (!strncmp(buffer, SDF_DELIMITER, 4)) {
            break;
          }
        }
        rewind(temp_fd.handle);
      }
    }
    if (temp_fd.handle) {
      fclose(temp_fd.handle);
      temp_fd.handle = NULL;
    }
    if (sdf_fd.handle) {
      fclose(sdf_fd.handle);
      sdf_fd.handle = NULL;
    }
    sprintf(filter_log_fd.name, "%s%cfilter_inter.log", od->align.filter_conf_dir, SEPARATOR);
    remove(filter_log_fd.name);
    for (template_num = 0;
      template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
      template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
      sprintf(buffer, "%s%c%04d.phar",
        od->align.filter_conf_dir, SEPARATOR,
        od->al.mol_info[template_object_num]->object_id);
      remove(buffer);
      sprintf(buffer, "%s%c%04d_phar_conf",
        od->align.filter_conf_dir, SEPARATOR,
        od->al.mol_info[template_object_num]->object_id);
      remove_recursive(buffer);
    }
  }
  free_array(od->al.phar_conf_list);
  od->al.phar_conf_list = NULL;
  
  return 0;
}


int filter_extract_split_phar_thread(void *pointer)
{
  char buffer[BUF_LEN];
  int error;
  int template_object_num;
  int template_num;
  int alloc_fail = 0;
  int assigned = 0;
  int n_phar_conf;
  int next_conf;
  int n_conf;
  int found;
  int phar_exist_ok;
  int pid;
  int result;
  FileDescriptor temp_fd;
  FileDescriptor phar_conf_fd;
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  
  
  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&phar_conf_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.exedir = ti->od.align.pharao_exe_path;
  if (!(prog_exe_info.proc_env = fill_env
    (&(ti->od), babel_env, ti->od.align.pharao_exe_path, 0))) {
    alloc_fail = 1;
  }
  prog_exe_info.stdout_fd = &temp_fd;
  prog_exe_info.stderr_fd = &temp_fd;
  prog_exe_info.sep_proc_grp = 1;
  error = 0;
  assigned = 1;
  while ((!error) && assigned) {
    template_num = 0;
    assigned = 0;
    while ((!assigned) && (template_num < ti->od.pel.numberlist[OBJECT_LIST]->size)) {
      template_object_num = ti->od.pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
      if (!(ti->od.al.mol_info[template_object_num]->done)) {
        #ifndef WIN32
        pthread_mutex_lock(ti->od.mel.mutex);
        #else
        WaitForSingleObject(ti->od.mel.mutex, INFINITE);
        #endif
        if (!(ti->od.al.mol_info[template_object_num]->done)) {
          ti->od.al.mol_info[template_object_num]->done = OBJECT_ASSIGNED;
          assigned = 1;
        }
        #ifndef WIN32
        pthread_mutex_unlock(ti->od.mel.mutex);
        #else
        ReleaseMutex(ti->od.mel.mutex);
        #endif
      }
      if (!assigned) {
        ++template_num;
      }
    }
    if (!assigned)  {
      break;
    }
    template_object_num = ti->od.pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    ti->od.al.task_list[template_num]->code = 0;
    ti->od.al.task_list[template_num]->data[TEMPLATE_OBJECT_NUM] = template_object_num;
    ti->od.al.task_list[template_num]->data[TEMPLATE_CONF_NUM] = -1;
    if (alloc_fail) {
      O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
      ti->od.al.task_list[template_num]->code = FL_OUT_OF_MEMORY;
      error = 1;
      continue;
    }
    sprintf(phar_conf_fd.name, "%s%c%04d.phar",
      ti->od.align.filter_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    phar_exist_ok = 0;
    if (fexist(phar_conf_fd.name)) {
      /*
      if XXXX.phar already exist, check that it is not malformed
      by checking that the number of pharmacophores matches the
      number of conformers in XXXX.sdf
      */
      if (!(phar_conf_fd.handle = fopen(phar_conf_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
        O3_ERROR_STRING(ti->od.al.task_list[template_num], phar_conf_fd.name);
        ti->od.al.task_list[template_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
        error = 1;
        continue;
      }
      n_phar_conf = 0;
      while (fgets(buffer, BUF_LEN, phar_conf_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        if (!strncmp(buffer, SDF_DELIMITER, 4)) {
          ++n_phar_conf;
        }
      }
      fclose(phar_conf_fd.handle);
      phar_conf_fd.handle = NULL;
      phar_exist_ok = (n_phar_conf == ti->od.pel.conf_population
        [TEMPLATE_DB]->pe[template_object_num]);
    }
    if (!phar_exist_ok) {
      sprintf(temp_fd.name, "%s%c%04d_phar_extract.log",
        ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[template_object_num]->object_id);
      sprintf(prog_exe_info.command_line,
        "%s -q %s %s -d %s%c%04d.sdf --dbType MOL -p %s",
        ti->od.align.pharao_exe,
        ((ti->od.align.type & ALIGN_TOGGLE_HYBRID_BIT) ? "" : PHARAO_NO_HYBRID),
        ((ti->od.align.type & ALIGN_TOGGLE_MERGE_BIT) ? PHARAO_MERGE : ""),
        ti->od.align.template_conf_dir, SEPARATOR,
        ti->od.al.mol_info[template_object_num]->object_id,
        phar_conf_fd.name);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[template_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      /*
      check if the Pharao computation was OK
      */
      if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
        O3_ERROR_STRING(ti->od.al.task_list[template_num], temp_fd.name);
        ti->od.al.task_list[template_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
        error = 1;
      }
      else if (fgrep(temp_fd.handle, buffer, "Error")) {
        O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
        ti->od.al.task_list[template_num]->code = FL_PHARAO_ERROR;
        error = 1;
      }
      if (temp_fd.handle) {
        fclose(temp_fd.handle);
        temp_fd.handle = NULL;
      }
      remove(temp_fd.name);
    }
    sprintf(buffer, "%s%c%04d_phar_conf",
      ti->od.align.filter_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    if (!dexist(buffer)) {
      #ifndef WIN32
      result = mkdir(buffer, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(buffer);
      #endif
      if (result) {
        O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
        O3_ERROR_STRING(ti->od.al.task_list[template_num], buffer);
        ti->od.al.task_list[template_num]->code = FL_CANNOT_CREATE_SCRDIR;
        error = 1;
        continue;
      }
    }
    /*
    now split XXXX.phar into n_phar_conf XXXX_YYYYYY.phar files
    and put them in a XXXX_phar_conf folder inside template_conf_dir
    */
    sprintf(temp_fd.name, "%s%c%04d.phar",
      ti->od.align.filter_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
      O3_ERROR_STRING(ti->od.al.task_list[template_num], temp_fd.name);
      ti->od.al.task_list[template_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
      error = 1;
      continue;
    }
    n_phar_conf = -1;
    next_conf = 1;
    n_conf = 0;
    phar_conf_fd.handle = NULL;
    while ((!error) && fgets(buffer, BUF_LEN, temp_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      if (next_conf) {
        if (phar_conf_fd.handle) {
          fclose(phar_conf_fd.handle);
        }
        next_conf = 0;
        ++n_phar_conf;
        n_conf = 0;
        found = 0;
        while (!(found = ((ti->od.al.phar_conf_list[n_conf]->object_num == template_object_num)
          && (ti->od.al.phar_conf_list[n_conf]->conf_num == n_phar_conf)))) {
          ++n_conf;
        }
        if (!found) {
          O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
          O3_ERROR_STRING(ti->od.al.task_list[template_num], temp_fd.name);
          ti->od.al.task_list[template_num]->code = FL_CANNOT_FIND_CONF;
          error = 1;
          continue;
        }
        ti->od.al.phar_conf_list[n_conf]->n_phar_points = -2;
        sprintf(phar_conf_fd.name, "%s%c%04d_phar_conf%c%04d_%06d.phar",
          ti->od.align.filter_conf_dir,
          SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
          SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
          n_phar_conf + 1);
        if (!(phar_conf_fd.handle = fopen(phar_conf_fd.name, "wb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
          O3_ERROR_STRING(ti->od.al.task_list[template_num], phar_conf_fd.name);
          ti->od.al.task_list[template_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
          error = 1;
          continue;
        }
      }
      fprintf(phar_conf_fd.handle, "%s\n", buffer);
      next_conf = (!strncmp(buffer, SDF_DELIMITER, 4));
      ++(ti->od.al.phar_conf_list[n_conf]->n_phar_points);
    }
    if (phar_conf_fd.handle) {
      fclose(phar_conf_fd.handle);
      phar_conf_fd.handle = NULL;
    }
    fclose(temp_fd.handle);
  }
  free_proc_env(prog_exe_info.proc_env);
  
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


int filter_intra_thread(void *pointer)
{
  char buffer[BUF_LEN];
  char *delete_list = NULL;
  int i;
  int error;
  int template_object_num;
  int template_num;
  int last_template_num = 0;
  int last_conf_num = 0;
  int delete_template_num;
  int delete_conf_num;
  int delete_file_ok = 0;
  int alloc_fail = 0;
  int assigned = 0;
  int n_phar_conf;
  int n_deleted;
  int n_appended;
  int line_n;
  int last_checked_line_n;
  int pid;
  double tanimoto;
  FileDescriptor temp_fd;
  FileDescriptor single_phar_conf_fd;
  FileDescriptor multi_phar_conf_fd;
  FileDescriptor filter_log_fd;
  FileDescriptor score_fd;
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  
  
  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&single_phar_conf_fd, 0, sizeof(FileDescriptor));
  memset(&multi_phar_conf_fd, 0, sizeof(FileDescriptor));
  memset(&filter_log_fd, 0, sizeof(FileDescriptor));
  memset(&score_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.exedir = ti->od.align.pharao_exe_path;
  if (!(prog_exe_info.proc_env = fill_env
    (&(ti->od), babel_env, ti->od.align.pharao_exe_path, 0))) {
    alloc_fail = 1;
  }
  prog_exe_info.stdout_fd = &temp_fd;
  prog_exe_info.stderr_fd = &temp_fd;
  prog_exe_info.sep_proc_grp = 1;
  error = 0;
  assigned = 1;
  while ((!error) && assigned) {
    template_num = 0;
    assigned = 0;
    while ((!assigned) && (template_num < ti->od.pel.numberlist[OBJECT_LIST]->size)) {
      template_object_num = ti->od.pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
      if (!(ti->od.al.mol_info[template_object_num]->done)) {
        #ifndef WIN32
        pthread_mutex_lock(ti->od.mel.mutex);
        #else
        WaitForSingleObject(ti->od.mel.mutex, INFINITE);
        #endif
        if (!(ti->od.al.mol_info[template_object_num]->done)) {
          ti->od.al.mol_info[template_object_num]->done = OBJECT_ASSIGNED;
          assigned = 1;
        }
        #ifndef WIN32
        pthread_mutex_unlock(ti->od.mel.mutex);
        #else
        ReleaseMutex(ti->od.mel.mutex);
        #endif
      }
      if (!assigned) {
        ++template_num;
      }
    }
    if (!assigned)  {
      break;
    }
    template_object_num = ti->od.pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    n_phar_conf = ti->od.pel.conf_population[TEMPLATE_DB]->pe[template_object_num];
    ti->od.al.task_list[template_num]->code = 0;
    ti->od.al.task_list[template_num]->data[TEMPLATE_OBJECT_NUM] = template_object_num;
    ti->od.al.task_list[template_num]->data[TEMPLATE_CONF_NUM] = -1;
    if (alloc_fail) {
      O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
      ti->od.al.task_list[template_num]->code = FL_OUT_OF_MEMORY;
      error = 1;
      continue;
    }
    if (!(delete_list = malloc(n_phar_conf))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
      ti->od.al.task_list[template_num]->code = FL_OUT_OF_MEMORY;
      error = 1;
      continue;
    }
    memset(delete_list, 0, n_phar_conf);
    /*
    check if we are restarting a previous filtering operation
    */
    sprintf(filter_log_fd.name, "%s%c%04d_phar_conf%cfilter_intra.log",
      ti->od.align.filter_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id, SEPARATOR);
    delete_file_ok = 0;
    n_deleted = 0;
    if ((filter_log_fd.handle = fopen(filter_log_fd.name, "rb"))) {
      if ((last_checked_line_n = fgrep(filter_log_fd.handle, buffer, "LAST_CHECKED"))) {
        sscanf(&buffer[13], "%d %d", &last_template_num, &last_conf_num);
        --last_template_num;
        --last_conf_num;
        if ((last_template_num == template_num)
          && (last_conf_num >= 0) && (last_conf_num < n_phar_conf)) {
          rewind(filter_log_fd.handle);
          delete_file_ok = 1;
          line_n = 0;
          while ((line_n < last_checked_line_n) && delete_file_ok
            && fgets(buffer, BUF_LEN, filter_log_fd.handle)) {
            ++line_n;
            delete_file_ok = 0;
            buffer[BUF_LEN - 1] = '\0';
            if (!strncmp(buffer, "LAST_CHECKED", 12)) {
              delete_file_ok = 1;
              continue;
            }
            sscanf(buffer, "%d %d", &delete_template_num, &delete_conf_num);
            --delete_template_num;
            --delete_conf_num;
            if ((delete_template_num == template_num)
              && (delete_conf_num >= 0) && (delete_conf_num < n_phar_conf)) {
              delete_file_ok = 1;
              delete_list[delete_conf_num] = 1;
              ++n_deleted;
            }
          }
        }
      }
      fclose(filter_log_fd.handle);
    }
    if (!delete_file_ok) {
      last_conf_num = 0;
      memset(delete_list, 0, n_phar_conf);
      remove(filter_log_fd.name);
    }
    sprintf(temp_fd.name, "%s%c%04d_phar_compare.log",
      ti->od.align.align_scratch, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    if (!(filter_log_fd.handle = fopen(filter_log_fd.name, "ab"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
      O3_ERROR_STRING(ti->od.al.task_list[template_num], filter_log_fd.name);
      ti->od.al.task_list[template_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
      error = 1;
      continue;
    }
    sprintf(multi_phar_conf_fd.name, "%s%c%04d_multi_phar_conf.phar",
      ti->od.align.align_scratch, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    sprintf(score_fd.name, "%s%c%04d_phar_compare.scores",
      ti->od.align.align_scratch, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    while ((!error) && (last_conf_num < n_phar_conf)) {
      if (!delete_list[last_conf_num]) {
        sprintf(single_phar_conf_fd.name, "%s%c%04d_phar_conf%c%04d_%06d.phar",
          ti->od.align.filter_conf_dir,
          SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
          SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
          last_conf_num + 1);
        remove(multi_phar_conf_fd.name);
        for (i = 0, n_appended = 0; ((!error) && (i < n_phar_conf)); ++i) {
          if ((i == last_conf_num) || delete_list[i]) {
            continue;
          }
          sprintf(buffer, "%s%c%04d_phar_conf%c%04d_%06d.phar",
            ti->od.align.filter_conf_dir,
            SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
            SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id, i + 1);
          if (!fcopy(buffer, multi_phar_conf_fd.name, "ab")) {
            O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
            O3_ERROR_STRING(ti->od.al.task_list[template_num], multi_phar_conf_fd.name);
            ti->od.al.task_list[template_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
            error = 1;
            continue;
          }
          ++n_appended;
        }
        if (error) {
          continue;
        }
        if (n_appended) {
          sprintf(prog_exe_info.command_line,
            "%s -q -r %s --refType PHAR -d %s --dbType PHAR -s %s",
            ti->od.align.pharao_exe, single_phar_conf_fd.name,
            multi_phar_conf_fd.name, score_fd.name);
          pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[template_num]->code));
          ext_program_wait(&prog_exe_info, pid);
          /*
          check if the Pharao computation was OK
          */
          if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
            O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
            O3_ERROR_STRING(ti->od.al.task_list[template_num], temp_fd.name);
            ti->od.al.task_list[template_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
            error = 1;
          }
          else if (fgrep(temp_fd.handle, buffer, "Error")) {
            O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
            ti->od.al.task_list[template_num]->code = FL_PHARAO_ERROR;
            error = 1;
          }
          if (temp_fd.handle) {
            fclose(temp_fd.handle);
            temp_fd.handle = NULL;
          }
          remove(temp_fd.name);
          if (error) {
            continue;
          }
          /*
          check the Pharao scores
          */
          if (!(score_fd.handle = fopen(score_fd.name, "rb"))) {
            O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
            O3_ERROR_STRING(ti->od.al.task_list[template_num], score_fd.name);
            ti->od.al.task_list[template_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
            error = 1;
            continue;
          }
          i = 0;
          while ((i < n_phar_conf) && fgets(buffer, BUF_LEN, score_fd.handle)) {
            sscanf(buffer, "%*s %*s %*s %*s %*s %*s %*s %*s %lf", &tanimoto);
            while ((i < n_phar_conf) && (delete_list[i] || (i == last_conf_num))) {
              ++i;
            }
            if (i == n_phar_conf) {
              break;
            }
            if (tanimoto > ti->od.align.level) {
              delete_list[i] = 1;
              ++n_deleted;
              fprintf(filter_log_fd.handle, "%d\t%d\n",
                ti->od.al.mol_info[template_object_num]->object_id, i + 1);
            }
            ++i;
          }
          fclose(score_fd.handle);
        }
        fprintf(filter_log_fd.handle, "LAST_CHECKED\t%d\t%d\n",
          ti->od.al.mol_info[template_object_num]->object_id, last_conf_num + 1);
        fflush(filter_log_fd.handle);
      }
      ++last_conf_num;
    }
    fclose(filter_log_fd.handle);
    if (error) {
      continue;
    }
    sprintf(multi_phar_conf_fd.name, "%s%c%04d.sdf",
      ti->od.align.template_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    if (!(multi_phar_conf_fd.handle = fopen(multi_phar_conf_fd.name, "rb"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
      O3_ERROR_STRING(ti->od.al.task_list[template_num], multi_phar_conf_fd.name);
      ti->od.al.task_list[template_num]->code = FL_CANNOT_READ_SDF_FILE;
      error = 1;
      continue;
    }
    sprintf(single_phar_conf_fd.name, "%s%c%04d.sdf",
      ti->od.align.filter_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    if (!(single_phar_conf_fd.handle = fopen(single_phar_conf_fd.name, "wb"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
      O3_ERROR_STRING(ti->od.al.task_list[template_num], single_phar_conf_fd.name);
      ti->od.al.task_list[template_num]->code = FL_CANNOT_WRITE_SDF_FILE;
      error = 1;
      continue;
    }
    for (i = 0; i < n_phar_conf; ++i) {
      if (find_conformation_in_sdf(multi_phar_conf_fd.handle,
        (delete_list[i] ? NULL : single_phar_conf_fd.handle), 0)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[template_num]);
        O3_ERROR_STRING(ti->od.al.task_list[template_num], multi_phar_conf_fd.name);
        ti->od.al.task_list[template_num]->code = FL_CANNOT_READ_SDF_FILE;
        error = 1;
        break;
      }
      while (fgets(buffer, BUF_LEN, multi_phar_conf_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        remove_newline(buffer);
        if (!delete_list[i]) {
          fprintf(single_phar_conf_fd.handle, "%s\n", buffer);
        }
        if (!strncmp(buffer, SDF_DELIMITER, 4)) {
          break;
        }
      }
    }
    fclose(multi_phar_conf_fd.handle);
    fclose(single_phar_conf_fd.handle);
    remove(filter_log_fd.name);
    sprintf(buffer, "%s%c%04d.phar",
      ti->od.align.filter_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    remove(buffer);
    sprintf(buffer, "%s%c%04d_phar_conf",
      ti->od.align.filter_conf_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id);
    remove_recursive(buffer);
  }
  free_proc_env(prog_exe_info.proc_env);
  
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}    


int filter_inter_thread(void *pointer)
{
  char buffer[BUF_LEN];
  int i;
  int j;
  int j_start;
  int error;
  int run;
  int n_runs;
  int n_calc_per_run;
  int n_calc_excess;
  int n_conf = 0;
  int n_conf_overall;
  int n_retained_conf;
  int n_appended;
  int template_object_num;
  int template_conf_num;
  int pid;
  double tanimoto;
  double tversky_ref;
  double tversky_db;
  FileDescriptor temp_fd;
  FileDescriptor single_phar_conf_fd;
  FileDescriptor multi_phar_conf_fd;
  FileDescriptor filter_log_fd;
  FileDescriptor score_fd;
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  
  
  ti = (ThreadInfo *)pointer;
  ti->od.al.task_list[ti->thread_num]->code = 0;
  memset(buffer, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&single_phar_conf_fd, 0, sizeof(FileDescriptor));
  memset(&multi_phar_conf_fd, 0, sizeof(FileDescriptor));
  memset(&filter_log_fd, 0, sizeof(FileDescriptor));
  memset(&score_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.exedir = ti->od.align.pharao_exe_path;
  if (!(prog_exe_info.proc_env = fill_env
    (&(ti->od), babel_env, ti->od.align.pharao_exe_path, 0))) {
    O3_ERROR_LOCATE(ti->od.al.task_list[ti->thread_num]);
    ti->od.al.task_list[ti->thread_num]->code = FL_OUT_OF_MEMORY;
    #ifndef WIN32
    pthread_exit(pointer);
    #else
    return 0;
    #endif
  }
  prog_exe_info.stdout_fd = &temp_fd;
  prog_exe_info.stderr_fd = &temp_fd;
  prog_exe_info.sep_proc_grp = 1;
  error = 0;
  n_runs = ti->n_calc / MAX_CONF_PER_PHARAO_RUN
    + ((ti->n_calc % MAX_CONF_PER_PHARAO_RUN ? 1 : 0));
  n_calc_per_run = ti->n_calc / n_runs;
  n_calc_excess = ti->n_calc % n_runs;
  n_retained_conf = ti->start;
  n_conf = ti->data[DATA_N_CONF];
  n_conf_overall = ti->data[DATA_N_CONF_OVERALL];
  j = 0;
  while ((j < n_conf_overall) && n_retained_conf) {
    while (ti->od.al.phar_conf_list[j]->delete) {
      ++j;
    }
    --n_retained_conf;
    ++j;
  }
  template_object_num = ti->od.al.phar_conf_list[n_conf]->object_num;
  template_conf_num = ti->od.al.phar_conf_list[n_conf]->conf_num;
  ti->od.al.task_list[ti->thread_num]->data[TEMPLATE_OBJECT_NUM] = template_object_num;
  ti->od.al.task_list[ti->thread_num]->data[TEMPLATE_CONF_NUM] = template_conf_num;
  sprintf(multi_phar_conf_fd.name, "%s%cmulti_phar_conf_t%02d.phar",
    ti->od.align.align_scratch, SEPARATOR, ti->thread_num);
  sprintf(temp_fd.name, "%s%cmulti_phar_compare_t%02d.log",
    ti->od.align.align_scratch, SEPARATOR, ti->thread_num);
  sprintf(score_fd.name, "%s%cphar_compare_t%02d.scores",
    ti->od.align.align_scratch, SEPARATOR, ti->thread_num);
  sprintf(single_phar_conf_fd.name, "%s%c%04d_phar_conf%c%04d_%06d.phar",
    ti->od.align.filter_conf_dir,
    SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
    SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
    template_conf_num + 1);
  for (run = 0; (!error) && (run < n_runs); ++run) {
    remove(multi_phar_conf_fd.name);
    j_start = j;
    for (i = 0, n_appended = 0; (!error) && (j < n_conf_overall)
      && (i < (n_calc_per_run + ((run == (n_runs - 1)) ? n_calc_excess : 0))); ++i, ++j) {
      while ((j < n_conf_overall) && (ti->od.al.phar_conf_list[j]->delete)) {
        ++j;
      }
      if ((j == n_conf_overall) || (n_conf == j)) {
        continue;
      }
      template_object_num = ti->od.al.phar_conf_list[j]->object_num;
      template_conf_num = ti->od.al.phar_conf_list[j]->conf_num;
      sprintf(buffer, "%s%c%04d_phar_conf%c%04d_%06d.phar",
        ti->od.align.filter_conf_dir,
        SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
        SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
        template_conf_num + 1);
      if (!fcopy(buffer, multi_phar_conf_fd.name, "ab")) {
        O3_ERROR_LOCATE(ti->od.al.task_list[ti->thread_num]);
        O3_ERROR_STRING(ti->od.al.task_list[ti->thread_num], multi_phar_conf_fd.name);
        ti->od.al.task_list[ti->thread_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
        error = 1;
        continue;
      }
      ++n_appended;
    }
    if (error) {
      continue;
    }
    if (n_appended) {
      sprintf(prog_exe_info.command_line,
        "%s -q -r %s --refType PHAR -d %s --dbType PHAR -s %s",
        ti->od.align.pharao_exe, single_phar_conf_fd.name,
        multi_phar_conf_fd.name, score_fd.name);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[ti->thread_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      /*
      check if the Pharao computation was OK
      */
      if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[ti->thread_num]);
        O3_ERROR_STRING(ti->od.al.task_list[ti->thread_num], temp_fd.name);
        ti->od.al.task_list[ti->thread_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
        error = 1;
      }
      else if (fgrep(temp_fd.handle, buffer, "Error")) {
        O3_ERROR_LOCATE(ti->od.al.task_list[ti->thread_num]);
        ti->od.al.task_list[ti->thread_num]->code = FL_PHARAO_ERROR;
        error = 1;
      }
      if (temp_fd.handle) {
        fclose(temp_fd.handle);
        temp_fd.handle = NULL;
      }
      remove(temp_fd.name);
      if (error) {
        continue;
      }
      /*
      check the Pharao scores
      */
      if (!(score_fd.handle = fopen(score_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[ti->thread_num]);
        O3_ERROR_STRING(ti->od.al.task_list[ti->thread_num], score_fd.name);
        ti->od.al.task_list[ti->thread_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
        error = 1;
        continue;
      }
      i = j_start;
      while ((i < n_conf_overall) && fgets(buffer, BUF_LEN, score_fd.handle)) {
        sscanf(buffer, "%*s %*s %*s %*s %*s %*s %*s %*s %lf %lf %lf",
          &tanimoto, &tversky_ref, &tversky_db);
        while ((i < n_conf_overall)
          && ((ti->od.al.phar_conf_list[i]->delete) || (i == n_conf))) {
          ++i;
        }
        if (i == n_conf_overall) {
          break;
        }
        if ((tanimoto > ti->od.align.level)
          || (tversky_db > ti->od.align.level)) {
          ti->od.al.phar_conf_list[i]->delete = 1;
        }
        ++i;
      }
      fclose(score_fd.handle);
    }
  }
  free_proc_env(prog_exe_info.proc_env);
  
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}    
