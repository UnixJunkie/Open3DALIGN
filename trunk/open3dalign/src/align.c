/*

align.c

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
#ifdef WIN32
#include <windows.h>
#endif


double score_alignment(O3Data *od, ConfInfo *template_conf, ConfInfo *fitted_conf, AtomPair *sdm, int pairs)
{
  int i;
  /*
  int not_fitted;
  int size_diff;
  double  penalty;
  */
  double score;
  double dist;
  double charge_diff;
  double charge_sum;
  
  for (i = 0, score = 0.0; i < pairs; ++i) {
    dist = squared_euclidean_distance
      (&(template_conf->coord[sdm[i].a[0] * 3]),
      &(fitted_conf->coord[sdm[i].a[1] * 3]));
    charge_diff = fabs(template_conf->atom[sdm[i].a[0]]->charge
      - fitted_conf->atom[sdm[i].a[1]]->charge);
    charge_sum = fabs(template_conf->atom[sdm[i].a[0]]->charge
      + fitted_conf->atom[sdm[i].a[1]]->charge);
    score += ((O3_SCORING_FUNCTION_ALPHA + (1.0 + O3_CHARGE_COEFF * charge_sum)
      / (1.0 + charge_diff)) * exp(-O3_SCORING_FUNCTION_BETA * dist));
  }
  /*
  not_fitted = fitted_conf->n_heavy_atoms - pairs;
  size_diff = fitted_conf->n_heavy_atoms - template_conf->n_heavy_atoms;
  if (size_diff < 0) {
    size_diff = 0;
  }
  penalty = 0.5 * (double)(not_fitted - size_diff);
  score += (50.0 - penalty);
  */
  
  return score;
}


int get_alignment_score(O3Data *od, FileDescriptor *fd, int object_num,
  double *score, int *best_template_object_num)
{
  char buffer[BUF_LEN];
  int found;
  int pos;
  
  
  memset(buffer, 0, BUF_LEN);
  if (!(fd->handle = fopen(fd->name, "rb"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), fd->name);
    return CANNOT_READ_TEMP_FILE;
  }
  if (find_conformation_in_sdf(fd->handle, NULL, object_num)) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), fd->name);
    fclose(fd->handle);
    return CANNOT_READ_TEMP_FILE;
  }
  pos = ftell(fd->handle);
  found = 0;
  if (score) {
    while ((!found) && fgets(buffer, BUF_LEN, fd->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      found = (strstr(buffer, ((od->align.type & ALIGN_PHARAO_BIT)
        ? "PHARAO_TANIMOTO" : "O3A_SCORE")) ? 1 : 0);
    }
    if ((!found) || (!fgets(buffer, BUF_LEN, fd->handle))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), fd->name);
      fclose(fd->handle);
      return CANNOT_READ_TEMP_FILE;
    }
    sscanf(buffer, "%lf", score);
    fseek(fd->handle, pos, SEEK_SET);
  }
  if (best_template_object_num) {
    found = 0;
    while ((!found) && fgets(buffer, BUF_LEN, fd->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      found = (strstr(buffer, "BEST_TEMPLATE_ID") ? 1 : 0);
    }
    if ((!found) || (!fgets(buffer, BUF_LEN, fd->handle))) {
      *best_template_object_num = -1;
    }
    else {
      buffer[BUF_LEN - 1] = '\0';
      sscanf(buffer, "%d", best_template_object_num);
      if ((*best_template_object_num > 0)
        && (*best_template_object_num <= od->grid.object_num)) {
        --(*best_template_object_num);
      }
      else {
        *best_template_object_num = -1;
      }
    }
  }
  fclose(fd->handle);

  return 0;
}


int align_random(O3Data *od)
{
  char buffer[BUF_LEN];
  int i;
  int j;
  int x;
  int object_num;
  int result;
  double centroid[3] = { 0.0, 0.0, 0.0 };
  double min_coord[3] = { 0.0, 0.0, 0.0 };
  double max_coord[3] = { 0.0, 0.0, 0.0 };
  double rand_trans[3] = { 0.0, 0.0, 0.0 };
  double rand_rot[3] = { 0.0, 0.0, 0.0 };
  double t_mat1[RT_MAT_SIZE];
  double t_mat2[RT_MAT_SIZE];
  double rt_mat[RT_MAT_SIZE];
  double t_vec1[RT_VEC_SIZE];
  double t_vec2[RT_VEC_SIZE];
  FileDescriptor mol_fd;
  FileDescriptor out_sdf_fd;
  AtomInfo **atom = NULL;
  
  
  memset(buffer, 0, BUF_LEN);
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(&out_sdf_fd, 0, sizeof(FileDescriptor));
  if (!(atom = (AtomInfo **)alloc_array
    (od->field.max_n_atoms + 1, sizeof(AtomInfo)))) {
    O3_ERROR_LOCATE(&(od->task));
    return OUT_OF_MEMORY;
  }
  set_random_seed(od, od->random_seed);
  sprintf(out_sdf_fd.name, "%s%c%04d-%04d_random.sdf",
    od->align.align_dir, SEPARATOR,
    od->al.mol_info[0]->object_id,
    od->al.mol_info[od->grid.object_num - 1]->object_id);
  if (!(out_sdf_fd.handle = fopen(out_sdf_fd.name, "wb"))) {
    free(atom);
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), out_sdf_fd.name);
    return CANNOT_WRITE_ALIGNED_SDF;
  }
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    result = fill_atom_info(od, &(od->task), atom, NULL, object_num, O3_MMFF94);
    if (result) {
      return result;
    }
    sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      fclose(out_sdf_fd.handle);
      free(atom);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    if (find_conformation_in_sdf(mol_fd.handle, out_sdf_fd.handle, 0)) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      fclose(mol_fd.handle);
      fclose(out_sdf_fd.handle);
      free(atom);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    memset(centroid, 0, 3 * sizeof(double));
    /*
    compute centroid, random rotation and random translation
    */
    for (x = 0; x < 3; ++x) {
      for (i = 0, j = 0; i < od->al.mol_info[object_num]->n_atoms; ++i) {
        /*
        if (!strcmp(atom[i]->element, "H")) {
          continue;
        }
        */
        if ((!j) || (atom[i]->coord[x] < min_coord[x])) {
          min_coord[x] = atom[i]->coord[x];
        }
        if ((!j) || (atom[i]->coord[x] > max_coord[x])) {
          max_coord[x] = atom[i]->coord[x];
        }
        centroid[x] += atom[i]->coord[x];
        ++j;
      }
      centroid[x] /= (double)(od->al.mol_info[object_num]->n_atoms);
      rand_rot[x] = angle2rad(360.0 * genrand_real(od));
    }
    /*
    initialize t_mat1 and t_mat2
    */
    memset(t_mat1, 0, RT_MAT_SIZE * sizeof(double));
    for (i = 0; i < RT_MAT_SIZE; i += 5) {
      t_mat1[i] = 1.0;
    }
    memcpy(t_mat2, t_mat1, RT_MAT_SIZE * sizeof(double));
    /*
    translate molecule to origin
    */
    cblas_daxpy(3, -1.0, centroid, 1, &t_mat1[3 * RT_VEC_SIZE], 1);
    /*
    translate molecule randomly
    */
    for (x = 0; x < 3; ++x) {
      rand_trans[x] = (max_coord[x] - min_coord[x])
        * RANDOM_TRANS_COEFF * (1.0 + (genrand_real(od) - 0.5));
      t_mat2[3 * RT_VEC_SIZE + x] = centroid[x] + rand_trans[x];
    }
    memset(t_vec1, 0, RT_VEC_SIZE * sizeof(double));
    t_vec1[3] = 1.0;
    /*
    rototranslate molecule randomly
    */
    prepare_rototrans_matrix(rt_mat, t_mat1, t_mat2, rand_rot);
    for (i = 0; i < od->al.mol_info[object_num]->n_atoms; ++i) {
      cblas_dcopy(3, atom[i]->coord, 1, t_vec1, 1);
      cblas_dgemv(CblasColMajor, CblasNoTrans,
        RT_VEC_SIZE, RT_VEC_SIZE, 1.0, rt_mat, RT_VEC_SIZE,
        t_vec1, 1, 0.0, t_vec2, 1);
      if (!fgets(buffer, BUF_LEN, mol_fd.handle)) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), mol_fd.name);
        fclose(mol_fd.handle);
        fclose(out_sdf_fd.handle);
        return CANNOT_READ_ORIGINAL_SDF;
      }
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      if (replace_coord(od->al.mol_info[object_num]->sdf_version, buffer, t_vec2)) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), mol_fd.name);
        fclose(mol_fd.handle);
        fclose(out_sdf_fd.handle);
        free(atom);
        return CANNOT_READ_ORIGINAL_SDF;
      }
      fprintf(out_sdf_fd.handle, "%s\n", buffer);
    }
    while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      fprintf(out_sdf_fd.handle, "%s\n", buffer);
    }
    fprintf(out_sdf_fd.handle, SDF_DELIMITER"\n");
    fclose(mol_fd.handle);
  }
  fclose(out_sdf_fd.handle);
  free(atom);

  return 0;
}


int is_in_list(IntPerm *list, int elem)
{
  int i;
  int found;
  
  
  for (i = 0, found = 0; (!found) && (i < list->size); ++i) {
    found = (elem == list->pe[i]);
  }
  
  return found;
}


int add_to_list(IntPerm **list, int elem)
{
  if (!(*list = int_perm_resize
    (*list, (*list)->size + 1))) {
    return OUT_OF_MEMORY;
  }
  (*list)->pe[(*list)->size - 1] = elem;
  
  return 0;
}


int remove_from_list(IntPerm **list, int elem)
{
  int i;
  int found;
  
  
  for (i = 0, found = -1; (found == -1) && (i < (*list)->size); ++i) {
    found = (((*list)->pe[i] == elem) ? i : -1);
  }
  if (found != -1) {
    for (i = found + 1; i < (*list)->size; ++i) {
      (*list)->pe[i - 1] = (*list)->pe[i];
    }
    if (!(*list = int_perm_resize
      (*list, (*list)->size - 1))) {
      return OUT_OF_MEMORY;
    }
  }
  
  return 0;
}


int preload_best_templates(O3Data *od, FileDescriptor *fd, int *skip, double *score)
{
  int object_num;
  int template_num = 0;
  int template_object_num;
  double file_score;
  
  
  if (!alignment_exists(od, fd)) {
    return 0;
  }
  for (object_num = 0, *skip = 1;
    *skip && (object_num < od->grid.object_num); ++object_num) {
    /*
    read template_object_num and score from file for this object
    */
    if (get_alignment_score(od, fd, object_num, &file_score, &template_object_num)) {
      *skip = 0;
      break;
    }
    if (template_object_num == -1) {
      *skip = 0;
      break;
    }
    /*
    check if it is already included in the current template list,
    (OBJECT_LIST, which initially is empty);
    if not, then add it at the end of the current template list
    */
    if (!is_in_list(od->pel.numberlist[OBJECT_LIST], template_object_num + 1)) {
      if (add_to_list(&(od->pel.numberlist[OBJECT_LIST]), template_object_num + 1)) {
        return OUT_OF_MEMORY;
      }
    }
    /*
    put the score in mol_info->score, multiplied by the gold coefficient
    */
    od->al.mol_info[object_num]->score = file_score * od->align.gold;
    *score += file_score;
    /*
    populate the per_object_template list with template_object_num
    */
    od->mel.per_object_template[object_num] = template_object_num;
  }
  if (!(*skip)) {
    /*
    if something went wrong, reset mol_info->score to 0.0 and
    clean up per_object_template
    */
    for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
      od->mel.per_object_template[object_num] = -1;
      od->al.mol_info[object_num]->score = 0.0;
    }
  }
  else {
    /*
    sort OBJECT_LIST by object number
    */
    qsort(od->pel.numberlist[OBJECT_LIST]->pe,
      od->pel.numberlist[OBJECT_LIST]->size, sizeof(int), compare_integers);
    for (template_num = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
      /*
      if the same templates read from template_file have not already been supplied by the user,
      then add them to the end of the user-supplied list (ID_LIST)
      */
      template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
      if (!is_in_list(od->pel.numberlist[ID_LIST], template_object_num + 1)) {
        if (add_to_list(&(od->pel.numberlist[ID_LIST]), template_object_num + 1)) {
          return OUT_OF_MEMORY;
        }
      }
    }
    qsort(od->pel.numberlist[ID_LIST]->pe,
      od->pel.numberlist[ID_LIST]->size, sizeof(int), compare_integers);
  }

  return 0;
}


int align_iterative(O3Data *od)
{
  char buffer[BUF_LEN];
  char element[MAX_FF_TYPE_LEN];
  char *dashed_line =
    "---------------------------------------------"
    "---------------------------------------------\n";
  int conv = 0;
  int i = 0;
  int iter = 0;
  int k;
  int n_candidate_templates = 0;
  int orig_n_candidate_templates = 0;
  int pool_size = 0;
  int incr = 0;
  int n_fail = 0;
  int failed_to_improve = 0;
  int n_heavy_atoms;
  int object_num;
  int template_num;
  int added_template_object_num = 0;
  int best_template_object_num = 0;
  int arch_template_object_num = 0;
  int template_object_num = 0;
  int found;
  int skip = 0;
  int n_gold = 0;
  int improve = 0;
  int old_size;
  int result;
  double max_score = 0.0;
  double sum_score = 0.0;
  double sum_score_print = 0.0;
  double prev_score = 0.0;
  double prev_score_print = 0.0;
  double best_template_score = 0.0;
  double best_template_score_print = 0.0;
  double msd;
  double overall_msd;
  double coord[2][3];
  FileDescriptor temp_fd;
  FileDescriptor best_fd;
  FileDescriptor *comp_fd[2];
  
  
  memset(buffer, 0, BUF_LEN);
  memset(element, 0, MAX_FF_TYPE_LEN);
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&best_fd, 0, sizeof(FileDescriptor));
  if (!(od->al.score_matrix = (double **)alloc_array
    (od->grid.object_num, od->grid.object_num * sizeof(double)))) {
    return OUT_OF_MEMORY;
  }
  if (!(od->al.candidate_template_object_list = (TemplateInfo **)alloc_array
    (od->grid.object_num, sizeof(TemplateInfo)))) {
    return OUT_OF_MEMORY;
  }
  if (!(od->mel.per_object_template = (int *)malloc
    (od->grid.object_num * sizeof(int)))) {
    return OUT_OF_MEMORY;
  }
  if (!(od->mel.per_object_template_temp = (int *)malloc
    (od->grid.object_num * sizeof(int)))) {
    return OUT_OF_MEMORY;
  }
  if (!(od->mel.score_temp = (double *)malloc
    (od->grid.object_num * sizeof(double)))) {
    return OUT_OF_MEMORY;
  }
  if (!(od->pel.best_template_object_list = int_perm_resize
    (od->pel.best_template_object_list, od->grid.object_num))) {
    return OUT_OF_MEMORY;
  }
  int_perm_resize(od->pel.best_template_object_list, 0);
  if (!(od->pel.failed_template_object_list = int_perm_resize
    (od->pel.failed_template_object_list, od->grid.object_num))) {
    return OUT_OF_MEMORY;
  }
  int_perm_resize(od->pel.failed_template_object_list, 0);
  for (i = 0; i < od->grid.object_num; ++i) {
    od->mel.per_object_template[i] = -1;
    od->al.mol_info[i]->score = 0.0;
  }
  if (!(od->pel.numberlist[ID_LIST] = int_perm_resize
    (od->pel.numberlist[ID_LIST], od->pel.numberlist[OBJECT_LIST]->size))) {
    return OUT_OF_MEMORY;
  }
  /*
  user specifies objects which may act as superposition templates
  through the OBJECT_LIST; the latter is copied into ID_LIST
  */
  memcpy(od->pel.numberlist[ID_LIST]->pe, od->pel.numberlist[OBJECT_LIST]->pe,
    od->pel.numberlist[OBJECT_LIST]->size * sizeof(int));
  if (od->align.template_file[0]) {
    /*
    if a template file is supplied, then OBJECT_LIST is populated with
    the n_gold templates, while ID_LIST is populated with the objects
    that the user initially put into OBJECT_LIST, to which the n_gold
    templates are added (if not already present)
    */
    int_perm_resize(od->pel.numberlist[OBJECT_LIST], 0);
    strcpy(temp_fd.name, od->align.template_file);
    if (preload_best_templates(od, &temp_fd, &skip, &sum_score)) {
      return OUT_OF_MEMORY;
    }
    if (!skip) {
      /*
      if something went wrong with preload_best_templates(),
      then the original user-supplied OBJECT_LIST is restored
      */
      if (!(od->pel.numberlist[OBJECT_LIST] = int_perm_resize
        (od->pel.numberlist[OBJECT_LIST], od->pel.numberlist[ID_LIST]->size))) {
        return OUT_OF_MEMORY;
      }
      memcpy(od->pel.numberlist[OBJECT_LIST]->pe, od->pel.numberlist[ID_LIST]->pe,
        od->pel.numberlist[ID_LIST]->size * sizeof(int));
    }
    else {
      /*
      instead, if all goes well best templates read from
      template_file are put into OBJECT_LIST; they are also added
      to ID_LIST in case they were not already there
      */
      n_gold = od->pel.numberlist[OBJECT_LIST]->size;
    }
  }
  tee_printf(od, "%6s%16s%16s    Best template object IDs\n%s",
    "Iter", "O3A_SCORE", "Delta", dashed_line);
  tee_flush(od);
  /*
  initially, best_template_object_list->size = 0
  */
  tee_printf(od, "n_gold = %d\n", n_gold);
  tee_flush(od);
  sum_score = 0.0;
  while ((!conv) && (iter < od->align.max_iter)
    && ((od->pel.best_template_object_list->size - n_gold) < od->pel.numberlist[ID_LIST]->size)) {
    if (!skip) {
      /*
      loop over all ID_LIST (user-supplied templates + n_gold templates)
      and copy in candidate_template_object_list
      those which are not already in best_template_object_list
      */
      for (template_num = 0, n_candidate_templates = 0;
        template_num < od->pel.numberlist[ID_LIST]->size; ++template_num) {
        template_object_num = od->pel.numberlist[ID_LIST]->pe[template_num] - 1;
        if (!is_in_list(od->pel.best_template_object_list, template_object_num)) {
          od->al.candidate_template_object_list
            [n_candidate_templates]->num = template_object_num;
          od->al.candidate_template_object_list
            [n_candidate_templates]->score =
            od->al.mol_info[template_object_num]->score;
          ++n_candidate_templates;
        }
      }
      /*
      sort the list by decreasing alignment score, so the first templates
      in the list will be those already better placed
      */
      qsort(od->al.candidate_template_object_list, n_candidate_templates,
        sizeof(TemplateInfo *), compare_template_score);
      tee_printf(od, "n_candidate_templates = %d\n", n_candidate_templates);
      orig_n_candidate_templates = n_candidate_templates;
      /*
      if this is not the first iteration, trim the list to keep
      only best-aligned templates; instead if it is the first
      it means that skip = 0 and so we have no clue about which
      templates might be the most promising and we need to try
      all of them
      */
      best_template_score = od->al.candidate_template_object_list[0]->score;
      pool_size = 1;
      /*
      keep as candidate templates those which have a score
      > (DEFAULT_BEST_PERCENT_ITER_ALIGN * best_score)
      */
      while ((pool_size < n_candidate_templates)
        && (od->al.candidate_template_object_list[pool_size]->score
        > (DEFAULT_BEST_PERCENT_ITER_ALIGN * best_template_score))) {
        ++pool_size;
      }
      /*
      increment the candidate template pool if previous attempts to find
      a better template failed
      */
      pool_size += incr;
      if (od->al.candidate_template_object_list[pool_size]->score
        < (DEFAULT_MIN_PERCENT_ITER_ALIGN * best_template_score)) {
        conv = 1;
        continue;
      }
      /*
      if this is the first iteration, and if there is no template_file,
      then keep all possible candidate templates wince we have no clue about
      those better placed, otherwise trim the list to those
      which have a score > (DEFAULT_BEST_PERCENT_ITER_ALIGN * best_score)
      plus eventual increments to do previous failures
      */
      if (iter && (n_candidate_templates > pool_size)) {
        n_candidate_templates = pool_size;
      }
      tee_printf(od, "first %d new templates: ", n_candidate_templates);
      for (i = 0; i < n_candidate_templates; ++i) {
        tee_printf(od, "%d,%.2lf ", od->al.candidate_template_object_list[i]->num + 1, od->al.candidate_template_object_list[i]->score);
      }
      tee_printf(od, "\ncurrent best templates: ");
      for (i = 0; i < od->pel.best_template_object_list->size; ++i) {
        tee_printf(od, "%d ", od->pel.best_template_object_list->pe[i] + 1);
      }
      /*
      trim OBJECT_LIST to n_candidate_templates and
      copy the most promising templates into OBJECT_LIST
      */
      old_size = od->pel.numberlist[OBJECT_LIST]->size;
      if (!(od->pel.numberlist[OBJECT_LIST] = int_perm_resize
        (od->pel.numberlist[OBJECT_LIST], old_size + n_candidate_templates))) {
        return OUT_OF_MEMORY;
      }
      for (i = 0; i < n_candidate_templates; ++i) {
        od->pel.numberlist[OBJECT_LIST]->pe[old_size + i] =
          od->al.candidate_template_object_list[i]->num + 1;
      }
      /*
      sort numberlist[OBJECT_LIST] (this can be removed later)
      */
      /*
      qsort(od->pel.numberlist[OBJECT_LIST]->pe, od->pel.numberlist[OBJECT_LIST]->size,
        sizeof(int), compare_integers);
      */
    }
    tee_printf(od, "\nOBJECT_LIST: ");
    for (i = 0; i < od->pel.numberlist[OBJECT_LIST]->size; ++i) {
      tee_printf(od, "%d ", od->pel.numberlist[OBJECT_LIST]->pe[i]);
    }
    tee_printf(od, "\n");
    tee_flush(od);
    /*
    do the alignment
    */
    if (!skip) {
      if ((result = align(od))) {
        tee_printf(od, "result = %d\n", result);
        tee_flush(od);
        return result;
      }
      free_array(od->al.task_list);
      od->al.task_list = NULL;
    }
    /*
    read all alignment scores for each template and
    copy them into score_matrix
    */
    for (template_num = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
      template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
      sprintf(temp_fd.name, "%s%c%04d-%04d_on_%04d.sdf",
        od->align.align_dir, SEPARATOR,
        od->al.mol_info[0]->object_id,
        od->al.mol_info[od->grid.object_num - 1]->object_id,
        od->al.mol_info[template_object_num]->object_id);
      for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
        if ((result = get_alignment_score(od, &temp_fd, object_num,
          &(od->al.score_matrix[object_num][template_object_num]), NULL))) {
          return result;
        }
        for (i = 0, found = 0; (!skip) && (!found) && (i < n_gold); ++i) {
          found = (template_object_num == od->pel.best_template_object_list->pe[i]);
        }
        if (found || skip) {
          /*
          if the gold template scores are being read,
          then multiply them by the gold coefficient
          */
          od->al.score_matrix[object_num][template_object_num] *= od->align.gold;
        }
      }
    }
    for (object_num = 0, best_template_score = 0.0, best_template_score_print = 0.0;
      object_num < od->grid.object_num; ++object_num) {
      best_template_object_num = od->mel.per_object_template[object_num];
      if (best_template_object_num != -1) {
        for (i = 0, found = 0; (!found) && (i < n_gold); ++i) {
          found = (best_template_object_num
            == od->pel.best_template_object_list->pe[i]);
        }
        best_template_score += od->al.score_matrix
          [object_num][best_template_object_num];
        best_template_score_print += od->al.score_matrix
          [object_num][best_template_object_num]
          / (found ? od->align.gold : 1.0);
      }
    }
    prev_score = best_template_score;
    prev_score_print = best_template_score_print;
    added_template_object_num = -1;
    if (!skip) {
      for (template_num = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
        /*
        loop over new templates and check if any of them improves the
        previous alignment score, and if so, which one improves it most
        */
        template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
        for (i = 0, found = 0; (!found) && (i < od->pel.failed_template_object_list->size); ++i) {
          /*
          if this template has already failed before, then skip it
          */
          found = (template_object_num == od->pel.failed_template_object_list->pe[i]);
        }
        if (found) {
          continue;
        }
        for (object_num = 0, improve = 0, sum_score = 0.0;
          object_num < od->grid.object_num; ++object_num) {
          /*
          loop over objects, not counting superimposition on self
          */
          found = (object_num == template_object_num);
          /*
          we do not want to change template to the objects which are gold templates
          */
          for (i = 0; (!found) && (i < n_gold); ++i) {
            found = (object_num == od->pel.best_template_object_list->pe[i]);
          }
          arch_template_object_num = od->mel.per_object_template[object_num];
          max_score = ((arch_template_object_num == -1)
            ? 0.0 : od->al.score_matrix[object_num][arch_template_object_num]);
          /*
          does this potential new template perform better than the previous one?
          */
          if ((!found) && ((od->al.score_matrix
            [object_num][template_object_num] - max_score) > ALMOST_ZERO)) {
            max_score = od->al.score_matrix[object_num][template_object_num];
            improve = 1;
            /*
            tee_printf(od, "1) object %d, replacing %d with %d, score %.4lf\n",
              object_num + 1, od->mel.per_object_template[object_num] + 1,
              template_object_num + 1, max_score);
            */
          }
          /*
          else {
            tee_printf(od, "1) object %d, keeping %d, score %.4lf\n",
              object_num + 1, od->mel.per_object_template[object_num] + 1,
              max_score);
          }
          */
          sum_score += max_score;
        }
        /*
        if at least one object improved its alignment score
        with respect to the previous template and the overall score
        is greater than the best improvement found so far, then store
        the new best score and the template to be added
        */
        if (improve && (sum_score > best_template_score)) {
          best_template_score = sum_score;
          added_template_object_num = template_object_num;
        }
      }
      if (added_template_object_num != -1) {
        /*
        if some good new template was found, then add it to the best_template_object_list
        */
        tee_printf(od, "adding %d to list\n", added_template_object_num + 1);
        if (!is_in_list(od->pel.best_template_object_list, added_template_object_num)) {
          add_to_list(&(od->pel.best_template_object_list), added_template_object_num);
        }
      }
      else {
        /*
        none of the templates succeeded in improving the previous alignment
        */
        tee_printf(od, "  no templates improve  ");
        conv = 1;
      }
      sum_score = 0.0;
      sum_score_print = 0.0;
    }
    else {
      /*
      if skip = 1, initial templates were read from template_file
      and now are copied into best_template_object_list
      */
      for (template_num = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
        template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
        add_to_list(&(od->pel.best_template_object_list), template_object_num);
      }
    }
    /*
    per_object_templates are backed up into per_object_template_temp
    */
    memcpy(od->mel.per_object_template_temp, od->mel.per_object_template,
      od->grid.object_num * sizeof(int));
    if (!conv) {
      ++iter;
      sprintf(best_fd.name, "%s%c%04d-%04d_on_best_template_iter_%04d.sdf",
        od->align.align_dir, SEPARATOR,
        od->al.mol_info[0]->object_id,
        od->al.mol_info[od->grid.object_num - 1]->object_id, iter);
      if (!(best_fd.handle = fopen(best_fd.name, "wb"))) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), best_fd.name);
        return CANNOT_WRITE_ALIGNED_SDF;
      }
      /*
      loop over objects
      */
      /*
      tee_printf(od, "1) od->pel.best_template_object_list->size = %d\n", od->pel.best_template_object_list->size);
      */
      for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
        /*
        loop over best_template_object_list, the first n_gold
        templates are those eventually read from template_file
        */
        for (i = 0, found = 0; (!found) && (i < n_gold); ++i) {
          found = (object_num == od->pel.best_template_object_list->pe[i]);
        }
        arch_template_object_num = od->mel.per_object_template[object_num];
        max_score = od->al.score_matrix[object_num][arch_template_object_num];
        /*
        if there is a template which improves overall alignment and
        for this object it performs better than the previous one
        then replace it
        */
        if ((added_template_object_num != -1) && (!found)
          && (object_num != added_template_object_num)
          && ((od->al.score_matrix[object_num][added_template_object_num]
          - max_score) > ALMOST_ZERO)) {
          max_score = od->al.score_matrix[object_num][added_template_object_num];
          best_template_object_num = added_template_object_num;
          /*
            tee_printf(od, "2) object %d, replacing %d with %d, score %.4lf\n",
            object_num + 1, od->mel.per_object_template[object_num] + 1,
            added_template_object_num + 1, max_score);
          */
        }
        else {
          best_template_object_num = od->mel.per_object_template[object_num];
          /*
          tee_printf(od, "2) object %d, keeping %d, score %.4lf\n",
            object_num + 1, od->mel.per_object_template[object_num] + 1,
            max_score);
          */
        }
        /*
        at the end of this loop we know which is the best template
        for this object, and we know the corresponding alignment score
        (max_score)
        */
        sum_score += max_score;
        for (i = 0, found = 0; (!found) && (i < n_gold); ++i) {
          found = (best_template_object_num
            == od->pel.best_template_object_list->pe[i]);
        }
        sum_score_print += max_score / (found ? od->align.gold : 1.0);
        /*
        the previous alignment score for this object
        is saved in score_temp; the new score goes in mol_info->score
        the current alignment template is stored in per_object_template
        */
        od->mel.score_temp[object_num] = od->al.mol_info[object_num]->score;
        od->al.mol_info[object_num]->score = max_score;
        od->mel.per_object_template[object_num] = best_template_object_num;
        /*
        tee_printf(od, "object_num = %d, old_pot,score = %d,%.4lf, new_pot,score = %d,%.4lf\n", object_num + 1,
          od->mel.per_object_template_temp[object_num] + 1, od->mel.score_temp[object_num],
          od->mel.per_object_template[object_num] + 1, od->al.mol_info[object_num]->score);
        */
        /*
        we pick from the file which contains the poses aligned on
        best_template_object_num the poses for objects aligned
        on that template, and we copy the corresponding section
        in the current iteration output file
        */
        sprintf(temp_fd.name, "%s%c%04d-%04d_on_%04d.sdf",
          od->align.align_dir, SEPARATOR,
          od->al.mol_info[0]->object_id,
          od->al.mol_info[od->grid.object_num - 1]->object_id,
          od->al.mol_info[best_template_object_num]->object_id);
        if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), temp_fd.name);
          fclose(best_fd.handle);
          return CANNOT_READ_ORIGINAL_SDF;
        }
        if (find_conformation_in_sdf(temp_fd.handle, best_fd.handle, object_num)) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), temp_fd.name);
          fclose(best_fd.handle);
          fclose(temp_fd.handle);
          return CANNOT_READ_TEMP_FILE;
        }
        found = 0;
        while ((!found) && fgets(buffer, BUF_LEN, temp_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          found = (!strncmp(buffer, SDF_DELIMITER, 4));
          if (!found) {
            fprintf(best_fd.handle, "%s\n", buffer);
          }
        }
        fclose(temp_fd.handle);
        if (!found) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), temp_fd.name);
          fclose(best_fd.handle);
          return CANNOT_READ_TEMP_FILE;
        }
        fprintf(best_fd.handle,
          ">  <BEST_TEMPLATE_ID>\n"
          "%d\n\n"
          SDF_DELIMITER"\n",
          best_template_object_num + 1);
      }
      fclose(best_fd.handle);
    }
    if ((iter > 1) && (!conv)) {
      /*
      conv = (((sum_score - found) / sum_score)
        < MIN_IMPROVEMENT_ITER_ALIGN);
      */
      conv = ((sum_score - prev_score) < ALMOST_ZERO);
      if (conv) {
        tee_printf(od, "sum_score = %.4lf, prev_score = %.4lf  score does not improve  ", sum_score, prev_score);
      }
    }
    if ((!skip) && (!conv) && od->align.template_file[0]) {
      strcpy(temp_fd.name, od->align.template_file);
      comp_fd[0] = &temp_fd;
      comp_fd[1] = &best_fd;
      for (k = 0; k <= 1; ++k) {
        if (!(comp_fd[k]->handle = fopen(comp_fd[k]->name, "rb"))) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), comp_fd[k]->name);
          if (k) {
            fclose(comp_fd[0]->handle);
          }
          return CANNOT_READ_ORIGINAL_SDF;
        }
      }
      for (object_num = 0, overall_msd = 0.0; object_num < od->grid.object_num; ++object_num) {
        for (k = 0; k <= 1; ++k) {
          if (find_conformation_in_sdf(comp_fd[k]->handle, NULL, 0)) {
            O3_ERROR_LOCATE(&(od->task));
            O3_ERROR_STRING(&(od->task), comp_fd[k]->name);
            fclose(comp_fd[0]->handle);
            fclose(comp_fd[1]->handle);
            return CANNOT_READ_TEMP_FILE;
          }
        }
        i = 0;
        n_heavy_atoms = 0;
        msd = 0.0;
        while (i < od->al.mol_info[object_num]->n_atoms) {
          for (k = 0; k <= 1; ++k) {
            if (!fgets(buffer, BUF_LEN, comp_fd[k]->handle)) {
              O3_ERROR_LOCATE(&(od->task));
              O3_ERROR_STRING(&(od->task), comp_fd[k]->name);
              break;
            }
            buffer[BUF_LEN - 1] = '\0';
            parse_sdf_coord_line(od->al.mol_info[object_num]->sdf_version,
              buffer, element, coord[k], NULL);
          }
          if (strcmp(element, "H")) {
            msd += squared_euclidean_distance(coord[0], coord[1]);
            ++n_heavy_atoms;
          }
          ++i;
        }
        for (k = 0; k <= 1; ++k) {
          found = 0;
          while ((!found) && fgets(buffer, BUF_LEN, comp_fd[k]->handle)) {
            buffer[BUF_LEN - 1] = '\0';
            found = (!strncmp(buffer, SDF_DELIMITER, 4));
          }
          if (!found) {
            O3_ERROR_LOCATE(&(od->task));
            O3_ERROR_STRING(&(od->task), comp_fd[k]->name);
            i = 0;
            break;
          }
        }
        if (i != od->al.mol_info[object_num]->n_atoms) {
          fclose(comp_fd[0]->handle);
          fclose(comp_fd[1]->handle);
          return CANNOT_READ_TEMP_FILE;
        }
        msd /= (double)n_heavy_atoms;
        overall_msd += msd;
      }
      fclose(comp_fd[0]->handle);
      fclose(comp_fd[1]->handle);
      conv = (overall_msd < square(DEFAULT_RMSD_ITER_ALIGN));
      if (conv) {
        tee_printf(od, "  same RMSD  ");
      }
    }
    /*
    OBJECT_LIST is reinitialized
    */
    int_perm_resize(od->pel.numberlist[OBJECT_LIST], 0);
    if (!conv) {
      int_perm_resize(od->pel.failed_template_object_list, 0);
      for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
        /*
        loop over all objects
        */
        if (od->mel.per_object_template[object_num] != od->mel.per_object_template_temp[object_num]) {
          /*
          if the template for this object changed from the last time,
          then remove the pose file where this object was acting as the template;
          in fact, if the template for this object changed, also the pose of
          this object will likely change, and therefore the whole dataset
          must be realigned on this object
          this will happen on the next align() call (vide infra)
          */
          sprintf(temp_fd.name, "%s%c%04d-%04d_on_%04d.sdf",
            od->align.align_dir, SEPARATOR,
            od->al.mol_info[0]->object_id,
            od->al.mol_info[od->grid.object_num - 1]->object_id,
            od->al.mol_info[object_num]->object_id);
          remove(temp_fd.name);
          /*
          if this object is also included in the best_template_object_list
          then add it to OBJECT_LIST, since the pose database corresponding
          to this object used as template has been removed since it was not
          valid anymore
          */
          if (is_in_list(od->pel.best_template_object_list, object_num)) {
            add_to_list(&(od->pel.numberlist[OBJECT_LIST]), object_num + 1);
          }
        }
      }
      /*
      at the end of the loop OBJECT_LIST is populated with objects
      whose template object changed since the last align() call
      */
      n_fail = 0;
      incr = 0;
    }
    else {
      /*
      if the last iteration produced no improvement, then we
      remove added_template_object_num from best_template_object_list
      */
      tee_printf(od, "removing %d from list\n", added_template_object_num + 1);
      remove_from_list(&(od->pel.best_template_object_list), added_template_object_num);
      add_to_list(&(od->pel.failed_template_object_list), added_template_object_num);
      memcpy(od->mel.per_object_template, od->mel.per_object_template_temp,
        od->grid.object_num * sizeof(int));
      /*
      the previous score is restored from the backup copy
      */
      for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
        od->al.mol_info[object_num]->score = od->mel.score_temp[object_num];
      }
      failed_to_improve = 1;
      /*
      if max_fail failure count is not reached, then increase
      the template pool size by one as well as the failure count
      */
      if (n_fail < od->align.max_fail) {
        conv = 0;
        ++incr;
        ++n_fail;
      }
    }
    if (!failed_to_improve) {
      /*
      if there was some improvement from the last time,
      then replace template_file with the new alignment
      */
      strcpy(od->align.template_file, best_fd.name);
    }
    else {
      /*
      if there was no improvement, then the score is the same as before
      */
      sum_score = prev_score;
      sum_score_print = prev_score_print;
    }
    tee_printf(od, "%6d%16.2lf%16.2lf",
      iter, sum_score_print, sum_score_print - prev_score_print);
    /*
    qsort(od->mel.best_template_object_list, i, sizeof(int), compare_integers);
    */
    /*
    for (j = 0; j < i; ++j) {
      tee_printf(od, "%d%c", od->mel.best_template_object_list[j] + 1,
        ((j == (i - 1)) ? '\n' : ','));
    }
    */
    if (failed_to_improve) {
      /*
      if there was no improvement and the original candidate_template_list
      was larger than the current template pool 
      */
      if (orig_n_candidate_templates >= (pool_size + /*iter + */incr)) {
        tee_printf(od, "    Increasing template pool size");
      }
      else {
        conv = 1;
      }
    }
    failed_to_improve = 0;
    tee_printf(od, "\n");
    tee_flush(od);
    skip = 0;
  }
  --iter;
  sum_score = prev_score;
  sum_score_print = prev_score_print;
  strcpy(best_fd.name, od->align.template_file);
  sprintf(temp_fd.name, "%s%c%04d-%04d_on_best_template_final.sdf",
    od->align.align_dir, SEPARATOR,
    od->al.mol_info[0]->object_id,
    od->al.mol_info[od->grid.object_num - 1]->object_id);
  if (best_fd.name[0]) {
    if (!fcopy(best_fd.name, temp_fd.name, "wb")) {
      return CANNOT_WRITE_ALIGNED_SDF;
    }
  }
  else {
    if ((result = join_mol_to_sdf(od, &(od->task),
      &temp_fd, od->align.candidate_dir))) {
      return result;
    }
  }
  tee_printf(od, "%s\n"
    "%16s    %s\n"
    "%s"
    "%16.2lf    %s\n\n",
    dashed_line, "Final O3A_SCORE", "Best alignment",
    dashed_line, sum_score_print, temp_fd.name);
  tee_flush(od);
  free_array(od->al.score_matrix);
  od->al.score_matrix = NULL;
  free_array(od->al.candidate_template_object_list);
  od->al.candidate_template_object_list = NULL;
  int_perm_free(od->pel.best_template_object_list);
  od->pel.best_template_object_list = NULL;
  free(od->mel.per_object_template);
  od->mel.per_object_template = NULL;
  free(od->mel.per_object_template_temp);
  od->mel.per_object_template_temp = NULL;
  free(od->mel.score_temp);
  od->mel.score_temp = NULL;
    
  return 0;
}


int align(O3Data *od)
{
  char buffer[BUF_LEN];
  char template_conf_string[MAX_NAME_LEN];
  char *dashed_line =
    "----------------------------------------------------\n";
  char *temp_dir_suffix[2] =
    { "template", "candidate" };
  char *temp_dir_name = NULL;
  char *filename = NULL;
  char first_char;
  int i;
  int j;
  int found;
  int object_num;
  int bit[2] = { ALIGN_MULTICONF_TEMPLATE_BIT,
    ALIGN_MULTICONF_CANDIDATE_BIT };
  int template_num = 0;
  int template_object_num;
  int template_conf_num;
  int best_template_object_num = 0;
  int best_template_conf_num = 0;
  int done_array_pos;
  int result;
  int n_threads;
  double score;
  double best_score;
  double overall_score;
  FileDescriptor mol_fd;
  FileDescriptor temp_fd;
  FileDescriptor best_fd;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  ThreadInfo **ti;
  void *align_func = NULL;


  ti = od->mel.thread_info;
  memset(buffer, 0, BUF_LEN);
  memset(template_conf_string, 0, MAX_NAME_LEN);
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&best_fd, 0, sizeof(FileDescriptor));
  for (i = 0; i < 2; ++i) {
    first_char = (i ? od->align.candidate_file[0] : od->align.template_file[0]);
    temp_dir_name = (i ? od->align.candidate_dir : od->align.template_dir);
    filename = (i ? od->align.candidate_file : od->align.template_file);
    if (first_char && (!(od->align.type & bit[i]))) {
      /*
      if an external SDF file is supplied and a single template/candidate
      is present for each object, then break the SDF file in multiple MOL files
      */
      sprintf(temp_dir_name, "%s%c%s",
        od->align.align_scratch, SEPARATOR, temp_dir_suffix[i]);
      if (dexist(temp_dir_name)) {
        remove_recursive(temp_dir_name);
      }
      #ifndef WIN32
      result = mkdir(temp_dir_name, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(temp_dir_name);
      #endif
      if (result) {
        return CANNOT_CREATE_DIRECTORY;
      }
      strcpy(temp_fd.name, filename);
      if ((result = break_sdf_to_mol(od, &(od->task), &temp_fd, temp_dir_name))) {
        return result;
      }
    }
    else {
      strcpy(temp_dir_name, od->field.mol_dir);
    }
  }
  if (od->align.type & (ALIGN_PHARAO_BIT | ALIGN_MIXED_BIT)) {
    /*
    if PHARAO or MIXED alignment types are chosen
    */
    if (!(od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT)) {
      sprintf(temp_fd.name, "%s%c%04d-%04d.sdf",
        od->align.align_scratch, SEPARATOR,
        od->al.mol_info[0]->object_id,
        od->al.mol_info[od->grid.object_num - 1]->object_id);
      /*
      a single candidate is present for each object,
      then concatenate all MOL files of currently loaded
      objects into a single SDF file as needed by Pharao
      */
      if ((result = join_mol_to_sdf(od, &(od->task),
        &temp_fd, od->align.candidate_dir))) {
        return result;
      }
    }
    /*
    if multi-conformational templates are used, then as many tasks
    as template conformations are needed, otherwise just 1 per template
    */
    for (i = 0, od->align.n_tasks = 0; i < od->pel.numberlist[OBJECT_LIST]->size; ++i) {
      template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[i] - 1;
      od->align.n_tasks += ((od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
        ? od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1);
    }
    if (!(od->al.task_list = (TaskInfo **)alloc_array(od->align.n_tasks, sizeof(TaskInfo)))) {
      O3_ERROR_LOCATE(&(od->task));
      return OUT_OF_MEMORY;
    }
    n_threads = fill_thread_info(od, od->align.n_tasks);
    /*
    evenly distribute Pharao templates among threads
    */
    if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
      for (i = 0, template_num = 0, template_conf_num = 0; i < n_threads; ++i) {
        for (j = ti[i]->start; j <= ti[i]->end; ++j) {
          template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
          od->al.task_list[j]->data[TEMPLATE_OBJECT_NUM] = template_object_num;
          od->al.task_list[j]->data[TEMPLATE_CONF_NUM] = template_conf_num;
          od->al.task_list[j]->data[MOVED_OBJECT_NUM] = -1;
          od->al.task_list[j]->data[MOVED_CONF_NUM] = -1;
          ++template_conf_num;
          if (template_conf_num == od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num]) {
            template_conf_num = 0;
            ++template_num;
          }
        }
      }
    }
    else {
      for (i = 0, j = 0; i < od->pel.numberlist[OBJECT_LIST]->size; ++i) {
        template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[i] - 1;
        od->al.task_list[j]->data[TEMPLATE_OBJECT_NUM] = template_object_num;
        od->al.task_list[j]->data[TEMPLATE_CONF_NUM] = -1;
        od->al.task_list[j]->data[MOVED_OBJECT_NUM] = -1;
        od->al.task_list[j]->data[MOVED_CONF_NUM] = -1;
        ++j;
      }
    }
    if (od->align.type & (ALIGN_MULTICONF_CANDIDATE_BIT | ALIGN_MIXED_BIT)) {
      /*
      if there are multiple candidates or the MIXED alignment
      type was chosen, then the first step will be to
      extract pharmacophores from templates with Pharao
      */
      #ifndef WIN32
      pthread_attr_init(&thread_attr);
      pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
      #endif
      for (i = 0; i < n_threads; ++i) {
        /*
        create the i-th thread
        */
        #ifndef WIN32
        od->error_code = pthread_create(&(od->thread_id[i]),
          &thread_attr, (void *(*)(void *))phar_extract_thread, ti[i]);
        if (od->error_code) {
          O3_ERROR_LOCATE(&(od->task));
          return CANNOT_CREATE_THREAD;
        }
        #else
        od->hThreadArray[i] = CreateThread(NULL, 0,
          (LPTHREAD_START_ROUTINE)phar_extract_thread,
          ti[i], 0, &(od->dwThreadIdArray[i]));
        if (!(od->hThreadArray[i])) {
          O3_ERROR_LOCATE(&(od->task));
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
          O3_ERROR_LOCATE(&(od->task));
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
      /*
      if errors occurred
      */
      if (i != od->align.n_tasks) {
        return ERROR_EXTRACTING_PHAR;
      }
      free_array(od->al.task_list);
      od->al.task_list = NULL;
    }
  }
  if ((!(od->align.type & ALIGN_PHARAO_BIT)) || ((od->align.type & ALIGN_PHARAO_BIT)
    && (od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT))) {
    /*
    if type != PHARAO or type == PHARAO with multiple candidates
    then there will be a task for each template conformation
    */
    od->align.n_tasks = od->grid.object_num;
    if (!(od->al.task_list = (TaskInfo **)alloc_array(od->align.n_tasks, sizeof(TaskInfo)))) {
      O3_ERROR_LOCATE(&(od->task));
      return OUT_OF_MEMORY;
    }
    for (i = 0, template_num = 0; i < od->pel.numberlist[OBJECT_LIST]->size; ++i) {
      template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[i] - 1;
      template_num += ((od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
        ? od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1);
    }
    if (!(od->al.done_objects = alloc_array(template_num, od->grid.object_num))) {
      O3_ERROR_LOCATE(&(od->task));
      return OUT_OF_MEMORY;
    }
    /*
    done_array_pos ranges from 0 to the overall number of
    template conformations (computed over all template objects)
    done_objects is a (done_array_pos, object_num) byte matrix
    */
    for (template_num = 0, done_array_pos = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
      template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
      for (template_conf_num = 0; template_conf_num < ((od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
        ? od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1); ++template_conf_num) {
        if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
          sprintf(template_conf_string, "_%06d", template_conf_num + 1);
        }
        sprintf(temp_fd.name, "%s%c%04d-%04d_on_%04d%s.sdf",
          od->align.align_dir, SEPARATOR,
          od->al.mol_info[0]->object_id,
          od->al.mol_info[od->grid.object_num - 1]->object_id,
          od->al.mol_info[template_object_num]->object_id,
          template_conf_string);
        if (alignment_exists(od, &temp_fd)) {
          for (i = 0; i < od->grid.object_num; ++i) {
            od->al.done_objects[done_array_pos][i] =
              OBJECT_ASSIGNED | OBJECT_FINISHED | OBJECT_COPIED;
          }
        }
        ++done_array_pos;
      }
    }
  }
  if (od->align.type & ALIGN_PHARAO_BIT) {
    align_func = (void *)((od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT)
      ? align_multi_pharao_thread : align_single_pharao_thread);
  }
  else {
    /*
    if atom-based alignment was chosen, the most
    efficient way to parallelise is by object
    */
    for (i = 0; i < od->grid.object_num; ++i) {
      od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM] = -1;
      od->al.task_list[i]->data[TEMPLATE_CONF_NUM] = -1;
      od->al.task_list[i]->data[MOVED_OBJECT_NUM] = i;
      od->al.task_list[i]->data[MOVED_CONF_NUM] = -1;
    }
    /*
    allocate and fill AtomInfo structures for each object
    */
    for (i = 0, od->field.max_n_heavy_atoms = 0; i < od->grid.object_num; ++i) {
      if (!(od->al.mol_info[i]->atom = (AtomInfo **)
        alloc_array(od->al.mol_info[i]->n_atoms + 1, sizeof(AtomInfo)))) {
        O3_ERROR_LOCATE(&(od->task));
        return OUT_OF_MEMORY;
      }
      result = fill_atom_info(od, &(od->task), od->al.mol_info[i]->atom, NULL, i, O3_MMFF94);
      if (result) {
        return result;
      }
      if ((!i) || (od->al.mol_info[i]->n_heavy_atoms > od->field.max_n_heavy_atoms)) {
        od->field.max_n_heavy_atoms = od->al.mol_info[i]->n_heavy_atoms;
      }
    }
    align_func = (void *)align_atombased_thread;
  }
  n_threads = fill_thread_info(od, od->align.n_tasks);
  #ifndef WIN32
  pthread_mutex_init(od->mel.mutex, NULL);
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
  #else
  if (!(*(od->mel.mutex) = CreateMutex(NULL, FALSE, NULL))) {
    O3_ERROR_LOCATE(&(od->task));
    return CANNOT_CREATE_THREAD;
  }
  #endif
  for (i = 0; i < n_threads; ++i) {
    /*
    create the i-th thread
    */
    #ifndef WIN32
    od->error_code = pthread_create(&(od->thread_id[i]),
      &thread_attr, (void *(*)(void *))align_func, ti[i]);
    if (od->error_code) {
      O3_ERROR_LOCATE(&(od->task));
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)align_func,
      ti[i], 0, &(od->dwThreadIdArray[i]));
    if (!(od->hThreadArray[i])) {
      O3_ERROR_LOCATE(&(od->task));
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
      O3_ERROR_LOCATE(&(od->task));
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
  /*
  if errors occurred
  */
  if (i != od->align.n_tasks) {
    return ERROR_IN_ALIGNMENT;
  }
  if ((od->align.type & ALIGN_ATOMBASED_BIT)
    || ((od->align.type & ALIGN_PHARAO_BIT) && (od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT))) {
    for (i = 0; i < od->grid.object_num; ++i) {
      free_array(od->al.mol_info[i]->atom);
      od->al.mol_info[i]->atom = NULL;
    }
    free_array(od->al.done_objects);
    od->al.done_objects = NULL;
  }
  if (!(od->align.type & ALIGN_ITERATIVE_TEMPLATE_BIT)) {
    tee_printf(od, "%8s%8s%16s%20s\n%s",
      "Template", "ID", "Conformer", ((od->align.type & ALIGN_PHARAO_BIT)
      ? "PHARAO_TANIMOTO" : "O3A_SCORE"), dashed_line);
    /*
    loop over user-requested template numbers
    */
    for (template_num = 0; template_num < od->pel.numberlist[OBJECT_LIST]->size; ++template_num) {
      template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
      /*
      loop over available conformations for each template
      */
      for (template_conf_num = 0; template_conf_num < ((od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
        ? od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1); ++template_conf_num) {
        /*
        set the aligned SDF name depending on whether the template
        includes multiple conformations or not
        */
        if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
          sprintf(template_conf_string, "_%06d", template_conf_num + 1);
        }
        sprintf(temp_fd.name, "%s%c%04d-%04d_on_%04d%s.sdf",
          od->align.align_dir, SEPARATOR,
          od->al.mol_info[0]->object_id,
          od->al.mol_info[od->grid.object_num - 1]->object_id,
          od->al.mol_info[template_object_num]->object_id,
          template_conf_string);
        for (object_num = 0, overall_score = 0.0; object_num < od->grid.object_num; ++object_num) {
          if ((result = get_alignment_score(od, &temp_fd, object_num, &score, NULL))) {
            return result;
          }
          overall_score += score;
        }
        tee_printf(od, "%8d%8d%16d%20.2lf\n", template_object_num + 1,
          od->al.mol_info[template_object_num]->object_id,
          template_conf_num + 1, overall_score);
      }
    }
    tee_printf(od, "%s\n", dashed_line);
  }
  if (od->align.type & ALIGN_KEEP_BEST_TEMPLATE_BIT) {
    sprintf(best_fd.name, "%s%c%04d-%04d_on_best_template.sdf",
      od->align.align_dir,
      SEPARATOR, od->al.mol_info[0]->object_id,
      od->al.mol_info[od->grid.object_num - 1]->object_id);
    if (!(best_fd.handle = fopen(best_fd.name, "wb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), best_fd.name);
      return CANNOT_WRITE_ALIGNED_SDF;
    }
    for (object_num = 0, overall_score = 0.0; object_num < od->grid.object_num; ++object_num) {
      best_template_object_num = -1;
      best_template_conf_num = -1;
      best_score = 0.0;
      for (i = 0; i < od->pel.numberlist[OBJECT_LIST]->size; ++i) {
        template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[i] - 1;
        for (template_conf_num = 0; template_conf_num < ((od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
          ? od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1); ++template_conf_num) {
          if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
            sprintf(template_conf_string, "_%06d", template_conf_num + 1);
          }
          sprintf(temp_fd.name, "%s%c%04d-%04d_on_%04d%s.sdf",
            od->align.align_dir, SEPARATOR,
            od->al.mol_info[0]->object_id,
            od->al.mol_info[od->grid.object_num - 1]->object_id,
            od->al.mol_info[template_object_num]->object_id,
            template_conf_string);
          if ((result = get_alignment_score(od, &temp_fd, object_num, &score, NULL))) {
            fclose(best_fd.handle);
            return result;
          }
          if ((best_template_object_num == -1) || ((score - best_score) > ALMOST_ZERO)) {
            best_score = score;
            best_template_object_num = template_object_num;
            if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
              best_template_conf_num = template_conf_num;
            }
          }
        }
      }
      overall_score += best_score;
      if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
        sprintf(template_conf_string, "_%06d", best_template_conf_num + 1);
      }
      sprintf(temp_fd.name, "%s%c%04d-%04d_on_%04d%s.sdf",
        od->align.align_dir, SEPARATOR,
        od->al.mol_info[0]->object_id,
        od->al.mol_info[od->grid.object_num - 1]->object_id,
        od->al.mol_info[best_template_object_num]->object_id,
        template_conf_string);
      if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), temp_fd.name);
        fclose(best_fd.handle);
        return CANNOT_READ_TEMP_FILE;
      }
      if (find_conformation_in_sdf(temp_fd.handle, best_fd.handle, object_num)) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), temp_fd.name);
        fclose(best_fd.handle);
        fclose(temp_fd.handle);
        return CANNOT_READ_TEMP_FILE;
      }
      found = 0;
      while ((!found) && fgets(buffer, BUF_LEN, temp_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        remove_newline(buffer);
        found = (!strncmp(buffer, SDF_DELIMITER, 4));
        if (!found) {
          fprintf(best_fd.handle, "%s\n", buffer);
        }
      }
      fprintf(best_fd.handle, ">  <BEST_TEMPLATE_ID>\n"
        "%d\n\n", best_template_object_num + 1);
      if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
        fprintf(best_fd.handle, ">  <BEST_TEMPLATE_CONF>\n"
          "%d\n\n", best_template_conf_num + 1);
      }
      fprintf(best_fd.handle, SDF_DELIMITER"\n");
      fclose(temp_fd.handle);
      if (!found) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), temp_fd.name);
        fclose(best_fd.handle);
        return CANNOT_READ_TEMP_FILE;
      }
    }
    tee_printf(od, "%-25s%20.2lf\n\n",
      "Best template combination score:", overall_score);
    fclose(best_fd.handle);
  }
  
  return 0;
}


#define O3_TEMPLATE      0
#define O3_MOVED      1
#define O3_BEST        2
#define O3_FITTED      3
#define O3_PROGRESS      4
#define O3_CAND        5
#define O3_TEMP1      7
#define O3_TEMP2      8
#define O3_GLOBAL      9


int alignment_exists(O3Data *od, FileDescriptor *sdf_fd)
{
  char buffer[BUF_LEN];
  int object_num;
  int line;
  int result;
  MolInfo temp_mol_info;


  memset(buffer, 0, BUF_LEN);
  memset(&temp_mol_info, 0, sizeof(MolInfo));
  result = 1;
  if (fexist(sdf_fd->name)) {
    if ((sdf_fd->handle = fopen(sdf_fd->name, "rb"))) {
      object_num = 0;
      line = 0;
      result = 0;
      while ((!result) && fgets(buffer, BUF_LEN, sdf_fd->handle)) {
        ++line;
        buffer[BUF_LEN - 1] = '\0';
        if (line == 4) {
          if (get_n_atoms_bonds(&temp_mol_info, sdf_fd->handle, buffer)) {
            result = PREMATURE_EOF;
          }
          else if ((temp_mol_info.n_atoms != od->al.mol_info[object_num]->n_atoms)
            || (temp_mol_info.n_bonds != od->al.mol_info[object_num]->n_bonds)) {
            result = N_ATOM_BOND_MISMATCH;
          }
        }
        else if (line > 4) {
          if (!strncmp(buffer, SDF_DELIMITER, 4)) {
            line = 0;
            ++object_num;
          }
        }
      }
      result = ((!result) && (object_num != od->grid.object_num));
      fclose(sdf_fd->handle);
      sdf_fd->handle = NULL;
    }
    if (result) {
      remove(sdf_fd->name);
    }
  }
  
  return (result ? 0 : 1);
}


int join_aligned_files(O3Data *od, int done_array_pos, char *error_filename)
{
  char buffer[BUF_LEN];
  char buffer2[BUF_LEN];
  char template_conf_string[MAX_NAME_LEN];
  int i = -1;
  int error = 0;
  int template_num;
  int template_object_num;
  int template_conf_num;
  int temp_done_array_pos = 0;
  int finished = 0;
  int assigned = 0;


  memset(buffer, 0, BUF_LEN);
  memset(buffer2, 0, BUF_LEN);
  memset(template_conf_string, 0, MAX_NAME_LEN);
  for (template_num = 0, temp_done_array_pos = 0; (!error) && (temp_done_array_pos <= done_array_pos)
    && (template_num < od->pel.numberlist[OBJECT_LIST]->size); ++template_num) {
    template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    for (template_conf_num = 0; (!error) && (temp_done_array_pos <= done_array_pos)
      && (template_conf_num < ((od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
      ? od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1)); ++template_conf_num) {
      for (i = 0, finished = 1; finished && (i < od->grid.object_num); ++i) {
        finished = (od->al.done_objects[temp_done_array_pos][i] & OBJECT_FINISHED);
      }
      if (finished) {
        assigned = 0;
        #ifndef WIN32
        pthread_mutex_lock(od->mel.mutex);
        #else
        WaitForSingleObject(*(od->mel.mutex), INFINITE);
        #endif
        if (!(od->al.done_objects[temp_done_array_pos][0] & OBJECT_COPIED)) {
          od->al.done_objects[temp_done_array_pos][0] |= OBJECT_COPIED;
          assigned = 1;
        }
        #ifndef WIN32
        pthread_mutex_unlock(od->mel.mutex);
        #else
        ReleaseMutex(*(od->mel.mutex));
        #endif
        if (assigned) {
          if (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
            sprintf(template_conf_string, "_%06d", template_conf_num + 1);
          }
          sprintf(buffer, "%s%c%04d-%04d_on_%04d%s.sdf",
            od->align.align_dir, SEPARATOR,
            od->al.mol_info[0]->object_id,
            od->al.mol_info[od->grid.object_num - 1]->object_id,
            od->al.mol_info[template_object_num]->object_id,
            template_conf_string);
          for (i = 0; i < od->grid.object_num; ++i) {
            sprintf(buffer2, "%s%c%04d%c%04d_on_%04d%s.sdf",
              od->align.align_scratch, SEPARATOR,
              od->al.mol_info[i]->object_id, SEPARATOR,
              od->al.mol_info[i]->object_id,
              od->al.mol_info[template_object_num]->object_id,
              template_conf_string);
            if (!fcopy(buffer2, buffer, (i ? "ab" : "wb"))) {
              error = 1;
              strcpy(error_filename, buffer2);
              break;
            }
            remove(buffer2);
          }
        }
      }
      ++temp_done_array_pos;
    }
  }
  
  return (error ? i : -1);
}


#ifndef WIN32
void *align_atombased_thread(void *pointer)
#else
DWORD align_atombased_thread(void *pointer)
#endif
{
  char *used[2];
  char buffer[BUF_LEN];
  char template_conf_string[MAX_NAME_LEN];
  int i;
  int k;
  int error;
  int conf_found = 0;
  int conf_file_error = 0;
  int moved_object_num;
  int template_num;
  int template_object_num;
  int template_conf_num;
  int moved_conf_num;
  int best_conf_num;
  int iter;
  int flag;
  int n_equiv = 0;
  int largest_n_heavy_atoms = 0;
  int coeff = 0;
  int options = 0;
  int pid;
  int sdm_threshold_iter;
  int alloc_fail = 0;
  int assigned = 0;
  int done_array_pos = 0;
  int pairs[O3_MAX_SLOT];
  double rt_mat[RT_MAT_SIZE];
  double n_equiv_heavy_msd;
  double pairs_heavy_msd[O3_MAX_SLOT];
  double score[O3_MAX_SLOT];
  double original_heavy_msd;
  double sdm_threshold_dist;
  double centroid[2][3];
  FileDescriptor mol_fd;
  FileDescriptor out_sdf_fd;
  FileDescriptor moved_fd;
  FileDescriptor pharao_sdf_fd;
  FileDescriptor phar_fd;
  FileDescriptor scores_fd;
  FileDescriptor temp_fd;
  LAPInfo li;
  AtomInfo **template_atom = NULL;
  AtomInfo **moved_atom = NULL;
  ConfInfo *conf[O3_MAX_SLOT];
  AtomPair *sdm[O3_MAX_SLOT];
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  
  
  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, BUF_LEN);
  memset(template_conf_string, 0, MAX_NAME_LEN);
  memset(pairs, 0, O3_MAX_SLOT * sizeof(int));
  memset(score, 0, O3_MAX_SLOT * sizeof(double));
  memset(conf, 0, O3_MAX_SLOT * sizeof(ConfInfo *));
  memset(sdm, 0, O3_MAX_SLOT * sizeof(AtomPair *));
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&pharao_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&phar_fd, 0, sizeof(FileDescriptor));
  memset(&scores_fd, 0, sizeof(FileDescriptor));
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(&out_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&moved_fd, 0, sizeof(FileDescriptor));
  memset(&li, 0, sizeof(LAPInfo));
  /*
  allocate two char vectors to store
  already used pairs
  */
  if (alloc_lap_info(&li, ti->od.field.max_n_heavy_atoms)) {
    alloc_fail = 1;
  }
  for (i = 0; i < O3_MAX_SLOT; ++i) {
    /*
    allocate MAX_SLOT conformations and SDM matrices
    */
    if (!(conf[i] = alloc_conf(ti->od.field.max_n_atoms))) {
      alloc_fail = 1;
    }
    if ((sdm[i] = (AtomPair *)malloc(square(ti->od.field.max_n_heavy_atoms) * sizeof(AtomPair)))) {
      memset(sdm[i], 0, square(ti->od.field.max_n_heavy_atoms) * sizeof(AtomPair));
    }
    else {
      alloc_fail = 1;
    }
  }
  for (i = 0; i < 2; ++i) {
    if (conf[i]) {
      /*
      for the first two conformations, allocate a
      (max_n_heavy_atoms, MAX_H_BINS) int matrix
      and also a (max_n_heavy_atoms) byte vector
      */
      if (!(conf[i]->h = (int **)alloc_array
        (ti->od.field.max_n_heavy_atoms, MAX_H_BINS * sizeof(int)))) {
        alloc_fail = 1;
      }
    }
    if (!(used[i] = malloc(ti->od.field.max_n_atoms))) {
      alloc_fail = 1;
    }
  }
  if (ti->od.align.type & ALIGN_MIXED_BIT) {
    prog_exe_info.exedir = ti->od.align.pharao_exe_path;
    if (!(prog_exe_info.proc_env = fill_env
      (&(ti->od), babel_env, ti->od.align.pharao_exe_path, 0))) {
      alloc_fail = 1;
    }
    prog_exe_info.stdout_fd = &temp_fd;
    prog_exe_info.stderr_fd = &temp_fd;
    prog_exe_info.sep_proc_grp = 1;
  }
  /*
  loop over all template objects
  */
  for (template_num = 0, done_array_pos = 0, error = 0; (!error)
    && (template_num < ti->od.pel.numberlist[OBJECT_LIST]->size); ++template_num) {
    template_object_num = ti->od.pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    template_atom = ti->od.al.mol_info[template_object_num]->atom;
    if (conf[O3_TEMPLATE]) {
      conf[O3_TEMPLATE]->atom = template_atom;
      conf[O3_TEMPLATE]->n_atoms = ti->od.al.mol_info[template_object_num]->n_atoms;
      conf[O3_TEMPLATE]->n_heavy_atoms = ti->od.al.mol_info[template_object_num]->n_heavy_atoms;
    }
    /*
    loop over conformations of this template object
    */
    for (template_conf_num = 0; (!error) && (template_conf_num < ((ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
      ? ti->od.pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1)); ++template_conf_num) {
      if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
        sprintf(template_conf_string, "_%06d", template_conf_num + 1);
      }
      conf_found = 0;
      if (conf[O3_TEMPLATE]) {
        if ((ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) || (ti->od.align.template_file[0])) {
          if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
            sprintf(temp_fd.name, "%s%c%04d.sdf", ti->od.align.template_conf_dir, SEPARATOR,
              ti->od.al.mol_info[template_object_num]->object_id);
          }
          else {
            sprintf(temp_fd.name, "%s%c%04d.mol", ti->od.align.template_dir, SEPARATOR,
              ti->od.al.mol_info[template_object_num]->object_id);
          }
          conf_file_error = 0;
          /*
          Read template conformation coordinates from a file,
          either an external file (for single-conformation template aligment)
          or a SDF conformational database (for multi-conformation template aligment)
          */
          if ((temp_fd.handle = fopen(temp_fd.name, "rb"))) {
            conf_found = (!find_conformation_in_sdf(temp_fd.handle, NULL, template_conf_num));
            if (conf_found) {
              i = 0;
              while ((i < ti->od.al.mol_info[template_object_num]->n_atoms)
                && fgets(buffer, BUF_LEN, temp_fd.handle)) {
                buffer[BUF_LEN - 1] = '\0';
                parse_sdf_coord_line(ti->od.al.mol_info[template_object_num]->sdf_version,
                  buffer, NULL, &(conf[O3_TEMPLATE]->coord[i * 3]), NULL);
                ++i;
              }
              conf_found = (i == ti->od.al.mol_info[template_object_num]->n_atoms);
            }
            fclose(temp_fd.handle);
          }
          else {
            conf_file_error = 1;
          }
        }
        else {
          /*
          read coordinates from the currently loaded structure
          */
          for (i = 0; i < ti->od.al.mol_info[template_object_num]->n_atoms; ++i) {
            cblas_dcopy(3, template_atom[i]->coord, 1, &(conf[O3_TEMPLATE]->coord[i * 3]), 1);
          }
        }
        /*
        compute H array for the template
        */
        compute_conf_h(conf[O3_TEMPLATE]);
      }
      assigned = 1;
      /*
      loop over all objects and find one not already assigned to some thread
      */
      while ((!error) && assigned) {
        moved_object_num = 0;
        assigned = 0;
        while ((!assigned) && (moved_object_num < ti->od.grid.object_num)) {
          if (!(ti->od.al.done_objects[done_array_pos][moved_object_num])) {
            sprintf(buffer, "%s%c%04d%c%04d_on_%04d%s.sdf",
              ti->od.align.align_scratch, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id,
              ti->od.al.mol_info[template_object_num]->object_id,
              template_conf_string);
            #ifndef WIN32
            pthread_mutex_lock(ti->od.mel.mutex);
            #else
            WaitForSingleObject(ti->od.mel.mutex, INFINITE);
            #endif
            if (!(ti->od.al.done_objects[done_array_pos][moved_object_num])) {
              ti->od.al.done_objects[done_array_pos][moved_object_num] = OBJECT_ASSIGNED;
              assigned = 1;
            }
            #ifndef WIN32
            pthread_mutex_unlock(ti->od.mel.mutex);
            #else
            ReleaseMutex(ti->od.mel.mutex);
            #endif
          }
          if (!assigned) {
            ++moved_object_num;
          }
        }
        /*
        if no object is left, stop searching
        */
        if (!assigned)  {
          break;
        }
        ti->od.al.task_list[moved_object_num]->code = 0;
        ti->od.al.task_list[moved_object_num]->data[TEMPLATE_OBJECT_NUM] = template_object_num;
        ti->od.al.task_list[moved_object_num]->data[TEMPLATE_CONF_NUM] = -1;
        ti->od.al.task_list[moved_object_num]->data[MOVED_CONF_NUM] = -1;
        if (alloc_fail) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          ti->od.al.task_list[moved_object_num]->code = FL_OUT_OF_MEMORY;
          error = 1;
          continue;
        }
        /*
        if it is a multi-conformational template
        */
        if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
          ti->od.al.task_list[moved_object_num]->data[TEMPLATE_CONF_NUM] = template_conf_num;
          if (conf_file_error) {
            O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
            /*
            sprintf(temp_fd.name, "%s%c%04d.sdf", ti->od.align.template_conf_dir, SEPARATOR,
              ti->od.al.mol_info[template_object_num]->object_id);
            */
            O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], temp_fd.name);
            ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_SDF_FILE;
            error = 1;
            continue;
          }
          if (!conf_found) {
            O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
            /*
            sprintf(temp_fd.name, "%s%c%04d.sdf", ti->od.align.template_conf_dir, SEPARATOR,
              ti->od.al.mol_info[template_object_num]->object_id);
            */
            O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], temp_fd.name);
            ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_FIND_CONF;
            error = 1;
            continue;
          }
        }
        /*
        if this is a mixed alignment
        */
        if (ti->od.align.type & ALIGN_MIXED_BIT) {
          if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
            /*
            align the conformational database of object "moved_object_num"
            on the current template object,conformation pair
            */
            sprintf(temp_fd.name, "%s%c%04d_on_%04d%s.log",
              ti->od.align.align_scratch, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id,
              ti->od.al.mol_info[template_object_num]->object_id,
              template_conf_string);
            sprintf(pharao_sdf_fd.name, "%s%c%04d_on_%04d%s_pharao.sdf",
              ti->od.align.align_scratch, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id,
              ti->od.al.mol_info[template_object_num]->object_id,
              template_conf_string);
            sprintf(phar_fd.name, "%s%c%04d%s.phar",
              ti->od.align.align_scratch, SEPARATOR,
              ti->od.al.mol_info[template_object_num]->object_id,
              template_conf_string);
            sprintf(scores_fd.name, "%s%c%04d_on_%04d%s.scores",
              ti->od.align.align_scratch, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id,
              ti->od.al.mol_info[template_object_num]->object_id,
              template_conf_string);
            sprintf(prog_exe_info.command_line,
              "%s -q %s %s -r %s --refType PHAR "
              "-d %s%c%04d.sdf --dbType MOL -s %s -o %s",
              ti->od.align.pharao_exe,
              ((ti->od.align.type & ALIGN_TOGGLE_HYBRID_BIT) ? "" : PHARAO_NO_HYBRID),
              ((ti->od.align.type & ALIGN_TOGGLE_MERGE_BIT) ? PHARAO_MERGE : ""),
              phar_fd.name, ti->od.align.candidate_conf_dir, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id,
              scores_fd.name, pharao_sdf_fd.name);
            pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[moved_object_num]->code));
            ext_program_wait(&prog_exe_info, pid);
            /*
            check if the Pharao computation was OK
            */
            if (ti->od.al.task_list[moved_object_num]->code) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              error = 1;
              continue;
            }
            if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], temp_fd.name);
              ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
              error = 1;
            }
            else if (fgrep(temp_fd.handle, buffer, "Error")) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], temp_fd.name);
              ti->od.al.task_list[moved_object_num]->code = FL_PHARAO_ERROR;
              error = 1;
            }
            if (temp_fd.handle) {
              fclose(temp_fd.handle);
              temp_fd.handle = NULL;
            }
            if (!error) {
              if (!(pharao_sdf_fd.handle = fopen(pharao_sdf_fd.name, "rb"))) {
                O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
                O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
                ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
                error = 1;
              }  
            }
            if (error) {
              continue;
            }
            remove(temp_fd.name);
            remove(scores_fd.name);
          }
          else {
            sprintf(pharao_sdf_fd.name, "%s%c%04d-%04d_on_%04d%s_pharao%c%04d.sdf",
              ti->od.align.align_scratch, SEPARATOR,
              ti->od.al.mol_info[0]->object_id,
              ti->od.al.mol_info[ti->od.grid.object_num - 1]->object_id,
              ti->od.al.mol_info[template_object_num]->object_id,
              template_conf_string,
              SEPARATOR, ti->od.al.mol_info[moved_object_num]->object_id);
            if (!(pharao_sdf_fd.handle = fopen(pharao_sdf_fd.name, "rb"))) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
              ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
              error = 1;
            }
            if (error) {
              continue;
            }
          }
        }
        /*
        create a scratch folder for the object currently assigned to this thread
        */
        sprintf(buffer, "%s%c%04d", ti->od.align.align_scratch, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id);
        if (!dexist(buffer)) {
          #ifndef WIN32
          error = mkdir(buffer, S_IRWXU | S_IRGRP | S_IROTH);
          #else
          error = mkdir(buffer);
          #endif
        }
        if (error) {
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_CREATE_SCRDIR;
          if (pharao_sdf_fd.handle) {
            fclose(pharao_sdf_fd.handle);
            pharao_sdf_fd.handle = NULL;
          }
          continue;
        }
        /*
        open in the scratch folder a SDF file for the assigned object,
        which will be aligned on the current template
        */
        sprintf(out_sdf_fd.name, "%s%c%04d_on_%04d%s.sdf",
          buffer, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id,
          ti->od.al.mol_info[template_object_num]->object_id,
          template_conf_string);
        if (!(out_sdf_fd.handle = fopen(out_sdf_fd.name, "wb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], out_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_WRITE_SDF_FILE;
          error = 1;
          if (pharao_sdf_fd.handle) {
            fclose(pharao_sdf_fd.handle);
            pharao_sdf_fd.handle = NULL;
          }
          continue;
        }
        /*
        if we are dealing with a multi-conformational candidate
        */
        if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
          /*
          find where the conformational database for this candidate object
          is stored and open it
          */
          sprintf(moved_fd.name, "%s%c%04d.sdf", ti->od.align.candidate_conf_dir,
            SEPARATOR, ti->od.al.mol_info[moved_object_num]->object_id);
          if (!(moved_fd.handle = fopen(moved_fd.name, "rb"))) {
            O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
            O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], moved_fd.name);
            ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_SDF_FILE;
            error = 1;
            if (pharao_sdf_fd.handle) {
              fclose(pharao_sdf_fd.handle);
              pharao_sdf_fd.handle = NULL;
            }
            fclose(out_sdf_fd.handle);
            out_sdf_fd.handle = NULL;
            continue;
          }
        }
        /*
        get AtomInfo for the candidate object
        */
        moved_atom = ti->od.al.mol_info[moved_object_num]->atom;
        for (i = 1; i <= 5; ++i) {
          conf[i]->atom = moved_atom;
          conf[i]->n_atoms = ti->od.al.mol_info[moved_object_num]->n_atoms;
          conf[i]->n_heavy_atoms = ti->od.al.mol_info[moved_object_num]->n_heavy_atoms;
        }
        /*
        loop over conformations of the candidate object
        */
        for (moved_conf_num = 0, best_conf_num = 0, score[O3_GLOBAL] = 0.0, pairs[O3_GLOBAL] = 0;
          (!error) && (moved_conf_num < ((ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT)
          ? ti->od.pel.conf_population[CANDIDATE_DB]->pe[moved_object_num] : 1)); ++moved_conf_num) {
          /*
          if the candidate has multiple conformations, get the relevant one
          from the SDF conformational database
          */
          if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
            ti->od.al.task_list[moved_object_num]->data[MOVED_CONF_NUM] = moved_conf_num;
            ti->od.al.task_list[moved_object_num]->code =
              find_conformation_in_sdf(moved_fd.handle, NULL, 0);
            if (ti->od.al.task_list[moved_object_num]->code) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], moved_fd.name);
              error = 1;
              fclose(moved_fd.handle);
              moved_fd.handle = NULL;
              if (pharao_sdf_fd.handle) {
                fclose(pharao_sdf_fd.handle);
                pharao_sdf_fd.handle = NULL;
              }
              fclose(out_sdf_fd.handle);
              out_sdf_fd.handle = NULL;
              continue;
            }
            k = 0;
            while ((k < ti->od.al.mol_info[moved_object_num]->n_atoms)
              && fgets(buffer, BUF_LEN, moved_fd.handle)) {
              buffer[BUF_LEN - 1] = '\0';
              parse_sdf_coord_line(ti->od.al.mol_info[moved_object_num]->sdf_version,
                buffer, NULL, &(conf[O3_MOVED]->coord[k * 3]), NULL);
              ++k;
            }
            if (k != ti->od.al.mol_info[moved_object_num]->n_atoms) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], moved_fd.name);
              ti->od.al.task_list[moved_object_num]->code =  FL_CANNOT_READ_SDF_FILE;
              error = 1;
              fclose(moved_fd.handle);
              moved_fd.handle = NULL;
              if (pharao_sdf_fd.handle) {
                fclose(pharao_sdf_fd.handle);
                pharao_sdf_fd.handle = NULL;
              }
              fclose(out_sdf_fd.handle);
              out_sdf_fd.handle = NULL;
              continue;
            }
            /*
            skip to the beginning of the next conformation
            */
            i = 0;
            while ((!i) && fgets(buffer, BUF_LEN, moved_fd.handle)) {
              i = (!strncmp(buffer, SDF_DELIMITER, 4));
            }
          }
          else {
            /*
            if the candidate has a single conformations, retrieve it from
            the candidate file (if supplied by the user)
            */
            if (ti->od.align.candidate_file[0]) {
              sprintf(moved_fd.name, "%s%c%04d.mol",
                ti->od.align.candidate_dir, SEPARATOR,
                ti->od.al.mol_info[moved_object_num]->object_id);
              if (!(moved_fd.handle = fopen(moved_fd.name, "rb"))) {
                O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
                O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], moved_fd.name);
                ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_SDF_FILE;
                error = 1;
                if (pharao_sdf_fd.handle) {
                  fclose(pharao_sdf_fd.handle);
                  pharao_sdf_fd.handle = NULL;
                }
                fclose(out_sdf_fd.handle);
                out_sdf_fd.handle = NULL;
                continue;
              }
              ti->od.al.task_list[moved_object_num]->code =
                find_conformation_in_sdf(moved_fd.handle, NULL, 0);
              if (ti->od.al.task_list[moved_object_num]->code) {
                O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
                O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], moved_fd.name);
                error = 1;
                fclose(moved_fd.handle);
                moved_fd.handle = NULL;
                if (pharao_sdf_fd.handle) {
                  fclose(pharao_sdf_fd.handle);
                  pharao_sdf_fd.handle = NULL;
                }
                fclose(out_sdf_fd.handle);
                out_sdf_fd.handle = NULL;
                continue;
              }
              k = 0;
              while ((k < ti->od.al.mol_info[moved_object_num]->n_atoms)
                && fgets(buffer, BUF_LEN, moved_fd.handle)) {
                buffer[BUF_LEN - 1] = '\0';
                parse_sdf_coord_line(ti->od.al.mol_info[moved_object_num]->sdf_version,
                  buffer, NULL, &(conf[O3_MOVED]->coord[k * 3]), NULL);
                ++k;
              }
              fclose(moved_fd.handle);
              moved_fd.handle = NULL;
              if (k != ti->od.al.mol_info[moved_object_num]->n_atoms) {
                O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
                O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], moved_fd.name);
                ti->od.al.task_list[moved_object_num]->code =  FL_CANNOT_READ_SDF_FILE;
                error = 1;
                if (pharao_sdf_fd.handle) {
                  fclose(pharao_sdf_fd.handle);
                  pharao_sdf_fd.handle = NULL;
                }
                fclose(out_sdf_fd.handle);
                out_sdf_fd.handle = NULL;
                continue;
              }
            }
            else {
              /*
              otherwise just get coordinates from the currently loaded objects
              */
              for (i = 0; i < ti->od.al.mol_info[moved_object_num]->n_atoms; ++i) {
                cblas_dcopy(3, moved_atom[i]->coord, 1, &(conf[O3_MOVED]->coord[i * 3]), 1);
              }
            }
          }
          /*
          compute the H array for the candidate conformation
          */
          compute_conf_h(conf[O3_MOVED]);
          for (options = 0, pairs[0] = 0, score[0] = 0.0;
            options <= (ti->od.align.type & ALIGN_TOGGLE_LOOP_BIT ? 0 : 1); ++options) {
            /*
            get object on which we are going to align
            the portion of dataset assigned to this thread
            */
            for (coeff = (ti->od.align.type & ALIGN_TOGGLE_LOOP_BIT ? 5 : 0), pairs[1] = 0, score[1] = 0.0;
              coeff < (ti->od.align.type & ALIGN_TOGGLE_LOOP_BIT ? 6 : 5); ++coeff) {
              /*
              compute the cost matrix for matching candidate to template
              */
              largest_n_heavy_atoms = compute_cost_matrix(&li, conf[O3_MOVED], conf[O3_TEMPLATE],
                MAX_H_BINS, coeff, (options ? MATCH_ATOM_TYPES_BIT : 0));
              /*
              find the lowest cost atom matching
              */
              lap(&li, largest_n_heavy_atoms);
              calc_conf_centroid(conf[O3_TEMPLATE], centroid[O3_TEMPLATE]);
              calc_conf_centroid(conf[O3_MOVED], centroid[O3_MOVED]);
              /*
              copy coordinates from O3_MOVED to O3_CAND
              */
              cblas_dcopy(conf[O3_MOVED]->n_atoms * 3, conf[O3_MOVED]->coord, 1, conf[O3_CAND]->coord, 1);
              for (i = 0; i < conf[O3_MOVED]->n_atoms; ++i) {
                /*
                subtract O3_MOVED centroid from O3_CAND
                (that is, O3_CAND is O3_MOVED with the center translated to origin)
                */
                cblas_daxpy(3, -1.0, centroid[O3_MOVED], 1, &(conf[O3_CAND]->coord[i * 3]), 1);
                /*
                add O3_TEMPLATE centroid from O3_CAND
                (that is, O3_CAND is O3_MOVED with the center translated to origin and
                then translated where the template centroid is)
                */
                cblas_daxpy(3, 1.0, centroid[O3_TEMPLATE], 1, &(conf[O3_CAND]->coord[i * 3]), 1);
              }
              /*
              filter the solution vector keeping only the safest n_equiv matches
              */
              n_equiv = filter_sol_vector(&li, conf[O3_CAND], conf[O3_TEMPLATE], sdm[O3_TEMP2], sdm[O3_TEMP1]);
              /*
              remove highest scores, keeping at least 3 superposition points
              */
              for (i = 3, pairs[2] = 0, score[2] = 0.0; i < n_equiv; ++i) {
                /*
                conf[O3_MOVED] is fitted on conf[O3_TEMPLATE]; fitted coordinates
                are placed in conf[O3_CAND]. If rms_algorithm fails, then
                fitted coordinates are not produced, hence candidate coordinates
                translated on the template centroid are copied in to fitted coordinates
                (better than nothing)
                */
                if (rms_algorithm
                  ((options ? USE_MMFF_WEIGHTS : USE_CHARGE_WEIGHTS), sdm[O3_TEMP1], i,
                  conf[O3_MOVED], conf[O3_TEMPLATE], conf[O3_CAND], rt_mat, &n_equiv_heavy_msd, NULL)) {
                  cblas_dcopy(conf[O3_MOVED]->n_atoms * 3, conf[O3_MOVED]->coord, 1, conf[O3_CAND]->coord, 1);
                  for (k = 0; k < conf[O3_MOVED]->n_atoms; ++k) {
                    cblas_daxpy(3, -1.0, centroid[O3_MOVED], 1, &(conf[O3_CAND]->coord[k * 3]), 1);
                    cblas_daxpy(3, 1.0, centroid[O3_TEMPLATE], 1, &(conf[O3_CAND]->coord[k * 3]), 1);
                  }
                }
                for (sdm_threshold_iter = (ti->od.align.type & ALIGN_TOGGLE_LOOP_BIT ? 2 : 0),
                  pairs[3] = 0, score[3] = 0.0; sdm_threshold_iter < 3; ++sdm_threshold_iter) {
                  pairs[4] = 0;
                  pairs_heavy_msd[4] = MAX_CUTOFF;
                  flag = 1;
                  iter = 0;
                  cblas_dcopy(conf[O3_CAND]->n_atoms * 3, conf[O3_CAND]->coord, 1, conf[O3_PROGRESS]->coord, 1);
                  while (flag && (iter < MAX_SDM_ITERATIONS)) {
                    ++iter;
                    /*
                    call sdm_algorithm
                    */
                    sdm_threshold_dist = SDM_THRESHOLD_START + (double)sdm_threshold_iter * SDM_THRESHOLD_STEP;
                    pairs[5] = sdm_algorithm(sdm[O3_TEMP2], conf[O3_PROGRESS], conf[O3_TEMPLATE], used, 0, sdm_threshold_dist);
                    if (pairs[5] < 3) {
                      break;
                    }
                    /*
                    call rms_algorithm
                    */
                    if (rms_algorithm((options ? USE_MMFF_WEIGHTS : USE_CHARGE_WEIGHTS),
                      sdm[O3_TEMP2], pairs[5], conf[O3_PROGRESS], conf[O3_TEMPLATE], conf[O3_FITTED],
                      rt_mat, &pairs_heavy_msd[5], NULL)) {
                      cblas_dcopy(conf[O3_PROGRESS]->n_atoms * 3, conf[O3_PROGRESS]->coord, 1, conf[O3_FITTED]->coord, 1);
                      break;
                    }
                    /*
                    keep looping until:
                    1) it is possible to increase the number of fitted pairs
                    2) it is not possible to increase the number of fitted pairs
                       anymore, but the msd is improved compared to the previous one
                    */
                    flag = ((pairs[5] > pairs[4]) || ((pairs[5] == pairs[4])
                      && ((pairs_heavy_msd[4] - pairs_heavy_msd[5]) > MSD_THRESHOLD)));
                    if (flag) {
                      pairs[4] = pairs[5];
                      pairs_heavy_msd[4] = pairs_heavy_msd[5];
                      memcpy(sdm[4], sdm[O3_TEMP2], pairs[4] * sizeof(AtomPair));
                      cblas_dcopy(conf[O3_FITTED]->n_atoms * 3, conf[O3_FITTED]->coord, 1, conf[O3_PROGRESS]->coord, 1);
                    }
                  }
                  score[4] = ((pairs[4] >= 3) ? score_alignment(&(ti->od), conf[O3_TEMPLATE],
                    conf[O3_FITTED], sdm[4], pairs[4]) : 0.0);
                  if ((score[4] - score[3]) > ALMOST_ZERO) {
                    pairs[3] = pairs[4];
                    score[3] = score[4];
                    memcpy(sdm[3], sdm[4], pairs[3] * sizeof(AtomPair));
                  }
                }
                if ((score[3] - score[2]) > ALMOST_ZERO) {
                  pairs[2] = pairs[3];
                  score[2] = score[3];
                  memcpy(sdm[2], sdm[3], pairs[2] * sizeof(AtomPair));
                }
              }
              if ((pairs[2] < 3) && (n_equiv >= 3)) {
                pairs[2] = n_equiv;
                memcpy(sdm[2], sdm[O3_TEMP1], pairs[2] * sizeof(AtomPair));
                score[2] = score_alignment(&(ti->od), conf[O3_TEMPLATE],
                  conf[O3_FITTED], sdm[2], pairs[2]);
              }
              if ((score[2] - score[1]) > ALMOST_ZERO) {
                score[1] = score[2];
                pairs[1] = pairs[2];
                memcpy(sdm[1], sdm[2], pairs[1] * sizeof(AtomPair));
              }
            }
            if ((score[1] - score[0]) > ALMOST_ZERO) {
              score[0] = score[1];
              pairs[0] = pairs[1];
              memcpy(sdm[0], sdm[1], pairs[0] * sizeof(AtomPair));
            }
          }
          if (ti->od.align.type & ALIGN_MIXED_BIT) {
            /*
            align moved_object on template_object
            */
            ti->od.al.task_list[moved_object_num]->code =
              find_conformation_in_sdf(pharao_sdf_fd.handle, NULL, 0);
            if (ti->od.al.task_list[moved_object_num]->code) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
              error = 1;
              if (moved_fd.handle) {
                fclose(moved_fd.handle);
                moved_fd.handle = NULL;
              }
              fclose(pharao_sdf_fd.handle);
              pharao_sdf_fd.handle = NULL;
              fclose(out_sdf_fd.handle);
              out_sdf_fd.handle = NULL;
              continue;
            }
            i = 0;
            while ((i < ti->od.al.mol_info[moved_object_num]->n_atoms)
              && fgets(buffer, BUF_LEN, pharao_sdf_fd.handle)) {
              buffer[BUF_LEN - 1] = '\0';
              parse_sdf_coord_line(ti->od.al.mol_info[moved_object_num]->sdf_version,
                buffer, NULL, &(conf[O3_CAND]->coord[i * 3]), NULL);
              ++i;
            }
            if (i != ti->od.al.mol_info[moved_object_num]->n_atoms) {
              O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
              O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
              ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_MOL_FILE;
              error = 1;
              if (moved_fd.handle) {
                fclose(moved_fd.handle);
                moved_fd.handle = NULL;
              }
              fclose(pharao_sdf_fd.handle);
              pharao_sdf_fd.handle = NULL;
              fclose(out_sdf_fd.handle);
              out_sdf_fd.handle = NULL;
              continue;
            }
            for (sdm_threshold_iter = 0, pairs[3] = 0, score[3] = 0.0; sdm_threshold_iter < 3; ++sdm_threshold_iter) {
              pairs[4] = 0;
              pairs_heavy_msd[4] = MAX_CUTOFF;
              flag = 1;
              iter = 0;
              cblas_dcopy(conf[O3_CAND]->n_atoms * 3, conf[O3_CAND]->coord, 1, conf[O3_PROGRESS]->coord, 1);
              while (flag && (iter < MAX_SDM_ITERATIONS)) {
                ++iter;
                /*
                call sdm_algorithm
                */
                sdm_threshold_dist = SDM_THRESHOLD_START + (double)sdm_threshold_iter * SDM_THRESHOLD_STEP;
                pairs[5] = sdm_algorithm(sdm[O3_TEMP2], conf[O3_PROGRESS], conf[O3_TEMPLATE], used, 0, sdm_threshold_dist);
                if (pairs[5] < 3) {
                  break;
                }
                /*
                call rms_algorithm
                */
                if (rms_algorithm((options ? USE_MMFF_WEIGHTS : USE_CHARGE_WEIGHTS),
                  sdm[O3_TEMP2], pairs[5], conf[O3_PROGRESS], conf[O3_TEMPLATE], conf[O3_FITTED],
                  rt_mat, &pairs_heavy_msd[5], NULL)) {
                  cblas_dcopy(conf[O3_PROGRESS]->n_atoms * 3, conf[O3_PROGRESS]->coord, 1, conf[O3_FITTED]->coord, 1);
                  break;
                }
                /*
                keep looping until:
                1) it is possible to increase the number of fitted pairs
                2) it is not possible to increase the number of fitted pairs
                   anymore, but the msd is improved compared to the previous one
                */
                flag = ((pairs[5] > pairs[4]) || ((pairs[5] == pairs[4])
                  && ((pairs_heavy_msd[4] - pairs_heavy_msd[5]) > MSD_THRESHOLD)));
                if (flag) {
                  pairs[4] = pairs[5];
                  pairs_heavy_msd[4] = pairs_heavy_msd[5];
                  memcpy(sdm[4], sdm[O3_TEMP2], pairs[4] * sizeof(AtomPair));
                  cblas_dcopy(conf[O3_CAND]->n_atoms * 3, conf[O3_FITTED]->coord, 1, conf[O3_PROGRESS]->coord, 1);
                }
              }
              score[4] = ((pairs[4] >= 3) ? score_alignment(&(ti->od), conf[O3_TEMPLATE],
                conf[O3_FITTED], sdm[4], pairs[4]) : 0.0);
              if ((score[4] - score[3]) > ALMOST_ZERO) {
                score[3] = score[4];
                pairs[3] = pairs[4];
                memcpy(sdm[3], sdm[4], pairs[3] * sizeof(AtomPair));
              }
            }
            if ((score[3] - score[0]) > 0.1) {
              score[0] = score[3];
              pairs[0] = pairs[3];
              memcpy(sdm[0], sdm[3], pairs[0] * sizeof(AtomPair));
            }
          }
          if (((score[0] - score[O3_GLOBAL]) > ALMOST_ZERO)
            || (!(ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT))) {
            score[O3_GLOBAL] = score[0];
            pairs[O3_GLOBAL] = pairs[0];
            best_conf_num = moved_conf_num;
            cblas_dcopy(conf[O3_MOVED]->n_atoms * 3, conf[O3_MOVED]->coord, 1, conf[O3_BEST]->coord, 1);
            memcpy(sdm[O3_GLOBAL], sdm[0], pairs[O3_GLOBAL] * sizeof(AtomPair));
          }
        }
        if (error) {
          continue;
        }
        if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
          fclose(moved_fd.handle);
          moved_fd.handle = NULL;
          cblas_dcopy(conf[O3_MOVED]->n_atoms * 3, conf[O3_BEST]->coord, 1, conf[O3_MOVED]->coord, 1);
        }
        rms_algorithm(USE_CHARGE_WEIGHTS, sdm[O3_GLOBAL], pairs[O3_GLOBAL], conf[O3_MOVED],
          conf[O3_TEMPLATE], conf[O3_FITTED], rt_mat, &pairs_heavy_msd[O3_GLOBAL],
          &original_heavy_msd);
        sprintf(mol_fd.name, "%s%c%04d.mol", ti->od.align.candidate_dir,
          SEPARATOR, ti->od.al.mol_info[moved_object_num]->object_id);
        if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], mol_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_MOL_FILE;
          fclose(out_sdf_fd.handle);
          out_sdf_fd.handle = NULL;
          error = 1;
          continue;
        }
        ti->od.al.task_list[moved_object_num]->code =
          find_conformation_in_sdf(mol_fd.handle, out_sdf_fd.handle, 0);
        if (ti->od.al.task_list[moved_object_num]->code) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], mol_fd.name);
          fclose(mol_fd.handle);
          mol_fd.handle = NULL;
          fclose(out_sdf_fd.handle);
          out_sdf_fd.handle = NULL;
          error = 1;
          continue;
        }
        k = 0;
        while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          if (k < conf[O3_FITTED]->n_atoms) {
            if (replace_coord(ti->od.al.mol_info[moved_object_num]->sdf_version,
              buffer, &(conf[O3_FITTED]->coord[k * 3]))) {
              break;
            }
            ++k;
          }
          fprintf(out_sdf_fd.handle, "%s\n", buffer);
        }
        if (k < conf[O3_MOVED]->n_atoms) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], mol_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_SDF_FILE;
          fclose(mol_fd.handle);
          mol_fd.handle = NULL;
          fclose(out_sdf_fd.handle);
          out_sdf_fd.handle = NULL;
          error = 1;
          continue;
        }
        while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          fprintf(out_sdf_fd.handle, "%s\n", buffer);
        }
        fclose(mol_fd.handle);
        fprintf(out_sdf_fd.handle, "\n"
          ">  <O3A_SCORE>\n"
          "%.4lf\n\n", score[O3_GLOBAL]);
        fprintf(out_sdf_fd.handle,
          ">  <O3A_ORIGINAL_SCORE>\n"
          "%.4lf\n\n", score_alignment(&(ti->od), conf[O3_TEMPLATE],
          conf[O3_MOVED],sdm[O3_GLOBAL], pairs[O3_GLOBAL]));
        if (ti->od.align.type & ALIGN_PRINT_RMSD_BIT) {
          fprintf(out_sdf_fd.handle,
            ">  <ORIGINAL_RMSD>\n"
            "%.4lf\n\n"
            ">  <ALIGNED_RMSD>\n"
            "%.4lf\n\n",
            sqrt(original_heavy_msd),
            sqrt(pairs_heavy_msd[O3_GLOBAL]));
        }
        if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
          fprintf(out_sdf_fd.handle, ">  <BEST_CANDIDATE_CONF>\n%d\n\n", best_conf_num + 1);
        }
        fprintf(out_sdf_fd.handle, SDF_DELIMITER"\n");
        if (out_sdf_fd.handle) {
          fclose(out_sdf_fd.handle);
          out_sdf_fd.handle = NULL;
        }
        if (pharao_sdf_fd.handle) {
          fclose(pharao_sdf_fd.handle);
          pharao_sdf_fd.handle = NULL;
          if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
            remove(pharao_sdf_fd.name);
          }
        }
        #ifndef WIN32
        pthread_mutex_lock(ti->od.mel.mutex);
        #else
        WaitForSingleObject(ti->od.mel.mutex, INFINITE);
        #endif
        ti->od.al.done_objects[done_array_pos][moved_object_num] |= OBJECT_FINISHED;
        #ifndef WIN32
        pthread_mutex_unlock(ti->od.mel.mutex);
        #else
        ReleaseMutex(ti->od.mel.mutex);
        #endif
      }
      if (error) {
        continue;
      }
      i = join_aligned_files(&(ti->od), done_array_pos, buffer);
      if (i != -1) {
        O3_ERROR_LOCATE(ti->od.al.task_list[i]);
        O3_ERROR_STRING(ti->od.al.task_list[i], buffer);
        ti->od.al.task_list[i]->code = FL_CANNOT_READ_SDF_FILE;
        error = 1;
      }
      if (!error) {
        ++done_array_pos;
      }
    }
  }
  for (i = 0; i < 2; ++i) {
    free_array(conf[i]->h);
    if (used[i]) {
      free(used[i]);
    }
  }
  for (i = 0; i < O3_MAX_SLOT; ++i) {
    free_conf(conf[i]);
    if (sdm[i]) {
      free(sdm[i]);
    }
  }
  free_lap_info(&li);
  free_proc_env(prog_exe_info.proc_env);
  
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


#ifndef WIN32
void *align_single_pharao_thread(void *pointer)
#else
DWORD align_single_pharao_thread(void *pointer)
#endif
{
  char buffer[BUF_LEN];
  char buffer2[BUF_LEN];
  char template_conf_string[MAX_NAME_LEN];
  char db_name[BUF_LEN];
  char pharao_temp_dir[BUF_LEN];
  int i;
  int task_num;
  int template_object_num;
  int template_conf_num;
  int object_num;
  int found = 0;
  int error;
  int result;
  int pid;
  FileDescriptor temp_fd;
  FileDescriptor inp_sdf_fd;
  FileDescriptor out_sdf_fd;
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  
  
  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, BUF_LEN);
  memset(buffer2, 0, BUF_LEN);
  memset(template_conf_string, 0, MAX_NAME_LEN);
  memset(db_name, 0, BUF_LEN);
  memset(pharao_temp_dir, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  memset(&inp_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&out_sdf_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.exedir = ti->od.align.pharao_exe_path;
  prog_exe_info.proc_env = fill_env(&(ti->od), babel_env, ti->od.align.pharao_exe_path, 0);
  prog_exe_info.stdout_fd = &temp_fd;
  prog_exe_info.stderr_fd = &temp_fd;
  prog_exe_info.sep_proc_grp = 1;
  /*
  loop over all tasks
  */
  for (task_num = ti->start, error = 0; (!error) && (task_num <= ti->end); ++task_num) {
    ti->od.al.task_list[task_num]->code = 0;
    if (!(prog_exe_info.proc_env)) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      ti->od.al.task_list[task_num]->code = FL_OUT_OF_MEMORY;
      error = 1;
      continue;
    }
    /*
    get object on which we are going
    to align the dataset
    */
    template_object_num = ti->od.al.task_list[task_num]->data[TEMPLATE_OBJECT_NUM];
    template_conf_num = ti->od.al.task_list[task_num]->data[TEMPLATE_CONF_NUM];
    sprintf(buffer, "%04d-%04d_on_%04d",
      ti->od.al.mol_info[0]->object_id,
      ti->od.al.mol_info[ti->od.grid.object_num - 1]->object_id,
      ti->od.al.mol_info[template_object_num]->object_id);
    sprintf(db_name, "%s%c%04d-%04d.sdf", ti->od.align.align_scratch,
      SEPARATOR, ti->od.al.mol_info[0]->object_id,
      ti->od.al.mol_info[ti->od.grid.object_num - 1]->object_id);
    if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
      sprintf(template_conf_string, "_%06d", template_conf_num + 1);
      sprintf(inp_sdf_fd.name, "%s%c%04d.sdf", ti->od.align.template_conf_dir, SEPARATOR,
        ti->od.al.mol_info[template_object_num]->object_id);
      if (!(inp_sdf_fd.handle = fopen(inp_sdf_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], inp_sdf_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        error = 1;
        continue;
      }
      sprintf(temp_fd.name, "%s%c%04d_%06d.mol", ti->od.align.template_dir,
        SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
        template_conf_num + 1);
      if (!(temp_fd.handle = fopen(temp_fd.name, "wb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], temp_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
        fclose(inp_sdf_fd.handle);
        error = 1;
        continue;
      }
      if (find_conformation_in_sdf(inp_sdf_fd.handle, temp_fd.handle, template_conf_num)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], inp_sdf_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_FIND_CONF;
        fclose(temp_fd.handle);
        fclose(inp_sdf_fd.handle);
        error = 1;
        continue;
      }
      i = 0;
      while ((!i) && fgets(buffer2, BUF_LEN, inp_sdf_fd.handle)) {
        buffer2[BUF_LEN - 1] = '\0';
        remove_newline(buffer2);
        i = (!strncmp(buffer2, MOL_DELIMITER, strlen(MOL_DELIMITER)));
        if (!i) {
          fprintf(temp_fd.handle, "%s\n", buffer2);
        }
      }
      fclose(temp_fd.handle);
      fclose(inp_sdf_fd.handle);
    }
    sprintf(out_sdf_fd.name, "%s%c%s%s.sdf", ti->od.align.align_dir,
      SEPARATOR, buffer, template_conf_string);
    sprintf(temp_fd.name, "%s%c%s%s.log", ti->od.align.align_scratch,
      SEPARATOR, buffer, template_conf_string);
    sprintf(prog_exe_info.command_line,
      "%s -q %s %s -r %s%c%04d%s.mol --refType MOL "
      "-d %s --dbType MOL "
      "-s %s%c%s%s.scores -o %s%c%s%s_pharao.sdf",
      ti->od.align.pharao_exe,
      ((ti->od.align.type & ALIGN_TOGGLE_HYBRID_BIT) ? "" : PHARAO_NO_HYBRID),
      ((ti->od.align.type & ALIGN_TOGGLE_MERGE_BIT) ? PHARAO_MERGE : ""),
      ti->od.align.template_dir, SEPARATOR,
      ti->od.al.mol_info[template_object_num]->object_id,
      template_conf_string, db_name,
      ti->od.align.align_scratch, SEPARATOR, buffer, template_conf_string,
      ti->od.align.align_scratch, SEPARATOR, buffer, template_conf_string);
    /*
    align objects on the current template
    */
    pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[task_num]->code));
    ext_program_wait(&prog_exe_info, pid);
    if (ti->od.al.task_list[task_num]->code) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      error = 1;
      continue;
    }
    /*
    check if the Pharao computation was OK
    */
    if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      O3_ERROR_STRING(ti->od.al.task_list[task_num], temp_fd.name);
      ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
      error = 1;
      continue;
    }
    if (fgrep(temp_fd.handle, buffer2, "Error")) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      O3_ERROR_STRING(ti->od.al.task_list[task_num], temp_fd.name);
      ti->od.al.task_list[task_num]->code = FL_PHARAO_ERROR;
    }
    fclose(temp_fd.handle);
    if (ti->od.al.task_list[task_num]->code) {
      error = 1;
      continue;
    }
    sprintf(pharao_temp_dir, "%s%c%s%s_pharao",
      ti->od.align.align_scratch, SEPARATOR, buffer, template_conf_string);
    if (!dexist(pharao_temp_dir)) {
      #ifndef WIN32
      result = mkdir(pharao_temp_dir, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(pharao_temp_dir);
      #endif
    }
    if (result) {
      ti->od.al.task_list[task_num]->code = FL_CANNOT_CREATE_SCRDIR;
      error = 1;
      continue;
    }
    sprintf(temp_fd.name, "%s.sdf", pharao_temp_dir);
    if ((ti->od.al.task_list[task_num]->code = break_sdf_to_sdf
      (&(ti->od), ti->od.al.task_list[task_num], &temp_fd, pharao_temp_dir))) {
      error = 1;
      continue;
    }
    if (!(out_sdf_fd.handle = fopen(out_sdf_fd.name, "wb"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      O3_ERROR_STRING(ti->od.al.task_list[task_num], out_sdf_fd.name);
      ti->od.al.task_list[task_num]->code = FL_CANNOT_WRITE_SDF_FILE;
      error = 1;
      continue;
    }
    for (object_num = 0; (!(ti->od.al.task_list[task_num]->code))
      && (object_num < ti->od.grid.object_num); ++object_num) {
      sprintf(inp_sdf_fd.name, "%s%c%04d.mol",
        ti->od.align.candidate_dir, SEPARATOR,
        ti->od.al.mol_info[object_num]->object_id);
      if (!(inp_sdf_fd.handle = fopen(inp_sdf_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], inp_sdf_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        continue;
      }
      sprintf(temp_fd.name, "%s%c%04d.sdf", pharao_temp_dir,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
      if (!(temp_fd.handle = fopen(temp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], temp_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        fclose(inp_sdf_fd.handle);
        continue;
      }
      ti->od.al.task_list[task_num]->code = find_conformation_in_sdf
        (inp_sdf_fd.handle, out_sdf_fd.handle, 0);
      if (ti->od.al.task_list[task_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], inp_sdf_fd.name);
        fclose(inp_sdf_fd.handle);
        fclose(temp_fd.handle);
        continue;
      }
      ti->od.al.task_list[task_num]->code = find_conformation_in_sdf
        (temp_fd.handle, NULL, 0);
      if (ti->od.al.task_list[task_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], temp_fd.name);
        fclose(inp_sdf_fd.handle);
        fclose(temp_fd.handle);
        continue;
      }
      i = 0;
      while ((i < ti->od.al.mol_info[object_num]->n_atoms)
        && fgets(buffer, BUF_LEN, inp_sdf_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        if (!fgets(buffer2, BUF_LEN, temp_fd.handle)) {
          break;
        }
        buffer2[BUF_LEN - 1] = '\0';
        remove_newline(buffer2);
        fprintf(out_sdf_fd.handle, "%s\n", buffer2);
        ++i;
      }
      if (i < ti->od.al.mol_info[object_num]->n_atoms) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], inp_sdf_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        fclose(inp_sdf_fd.handle);
        fclose(temp_fd.handle);
        continue;
      }
      found = 0;
      while ((!found) && fgets(buffer, BUF_LEN, inp_sdf_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        remove_newline(buffer);
        fprintf(out_sdf_fd.handle, "%s\n", buffer);
        found = (!strncmp(buffer, MOL_DELIMITER, 4));
      }
      if (!found) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], inp_sdf_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        fclose(inp_sdf_fd.handle);
        fclose(temp_fd.handle);
        continue;
      }
      found = 0;
      while ((!found) && fgets(buffer2, BUF_LEN, temp_fd.handle)) {
        buffer2[BUF_LEN - 1] = '\0';
        remove_newline(buffer2);
        found = (!strncasecmp(buffer2, ">  <PHARAO_TANIMOTO>", 20));
      }
      if (!found) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], temp_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        fclose(inp_sdf_fd.handle);
        fclose(temp_fd.handle);
        continue;
      }
      fprintf(out_sdf_fd.handle, "%s\n", buffer2);
      if (!fgets(buffer2, BUF_LEN, temp_fd.handle)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], temp_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        fclose(inp_sdf_fd.handle);
        fclose(temp_fd.handle);
        continue;
      }
      buffer2[BUF_LEN - 1] = '\0';
      remove_newline(buffer2);
      fprintf(out_sdf_fd.handle, "%s\n%s\n", buffer2, SDF_DELIMITER);
      fclose(temp_fd.handle);
      fclose(inp_sdf_fd.handle);
      remove(temp_fd.name);
    }
    fclose(out_sdf_fd.handle);
  }
  free_proc_env(prog_exe_info.proc_env);
  
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


#ifndef WIN32
void *phar_extract_thread(void *pointer)
#else
DWORD phar_extract_thread(void *pointer)
#endif
{
  char buffer[BUF_LEN];
  char template_conf_string[MAX_NAME_LEN];
  char pharao_temp_dir[BUF_LEN];
  int task_num;
  int template_object_num;
  int template_conf_num;
  int eof;
  int error = 0;
  int result;
  int pid;
  FileDescriptor log_fd;
  FileDescriptor multi_conf_sdf_fd;
  FileDescriptor single_conf_mol_fd;
  FileDescriptor pharao_sdf_fd;
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  
  
  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, BUF_LEN);
  memset(template_conf_string, 0, MAX_NAME_LEN);
  memset(pharao_temp_dir, 0, BUF_LEN);
  memset(&log_fd, 0, sizeof(FileDescriptor));
  memset(&pharao_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&multi_conf_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&single_conf_mol_fd, 0, sizeof(FileDescriptor));
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  prog_exe_info.exedir = ti->od.align.pharao_exe_path;
  prog_exe_info.proc_env = fill_env(&(ti->od), babel_env, ti->od.align.pharao_exe_path, 0);
  prog_exe_info.stdout_fd = &log_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.sep_proc_grp = 1;
  for (task_num = ti->start, error = 0; (!error) && (task_num <= ti->end); ++task_num) {
    ti->od.al.task_list[task_num]->code = 0;
    if (!(prog_exe_info.proc_env)) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      ti->od.al.task_list[task_num]->code = FL_OUT_OF_MEMORY;
      error = 1;
      continue;
    }
    template_object_num = ti->od.al.task_list[task_num]->data[TEMPLATE_OBJECT_NUM];
    template_conf_num = ti->od.al.task_list[task_num]->data[TEMPLATE_CONF_NUM];
    if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
      sprintf(multi_conf_sdf_fd.name, "%s%c%04d.sdf", ti->od.align.template_conf_dir,
        SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id);
      if (!(multi_conf_sdf_fd.handle = fopen(multi_conf_sdf_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], multi_conf_sdf_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        error = 1;
        continue;
      }
      sprintf(single_conf_mol_fd.name, "%s%c%04d_%06d.mol", ti->od.align.align_scratch,
        SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id, template_conf_num + 1);
      if (!(single_conf_mol_fd.handle = fopen(single_conf_mol_fd.name, "wb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], single_conf_mol_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
        fclose(multi_conf_sdf_fd.handle);
        error = 1;
        continue;
      }
      ti->od.al.task_list[task_num]->code = find_conformation_in_sdf
        (multi_conf_sdf_fd.handle, single_conf_mol_fd.handle, template_conf_num);
      if (ti->od.al.task_list[task_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], multi_conf_sdf_fd.name);
        fclose(multi_conf_sdf_fd.handle);
        fclose(single_conf_mol_fd.handle);
        error = 1;
        continue;
      }
      eof = 0;
      while ((!eof) && fgets(buffer, BUF_LEN, multi_conf_sdf_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        remove_newline(buffer);
        eof = (!strncmp(buffer, SDF_DELIMITER, 4));
        if (!eof) {
          fprintf(single_conf_mol_fd.handle, "%s\n", buffer);
        }
      }
      fclose(multi_conf_sdf_fd.handle);
      fclose(single_conf_mol_fd.handle);
      if (!eof) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], multi_conf_sdf_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_SDF_FILE;
        error = 1;
        continue;
      }
      sprintf(log_fd.name, "%s%c%04d_%06d_phar.log", ti->od.align.align_scratch,
        SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id,
        template_conf_num + 1);
      sprintf(prog_exe_info.command_line,
        "%s -q %s %s -d %s --dbType MOL -p %s%c%04d_%06d.phar", ti->od.align.pharao_exe,
        ((ti->od.align.type & ALIGN_TOGGLE_HYBRID_BIT) ? "" : PHARAO_NO_HYBRID),
        ((ti->od.align.type & ALIGN_TOGGLE_MERGE_BIT) ? PHARAO_MERGE : ""),
        single_conf_mol_fd.name, ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[template_object_num]->object_id, template_conf_num + 1);
    }
    else {
      sprintf(log_fd.name, "%s%c%04d_phar.log", ti->od.align.align_scratch,
        SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id);
      sprintf(prog_exe_info.command_line,
        "%s -q %s %s -d %s%c%04d.mol --dbType MOL -p %s%c%04d.phar", ti->od.align.pharao_exe,
        ((ti->od.align.type & ALIGN_TOGGLE_HYBRID_BIT) ? "" : PHARAO_NO_HYBRID),
        ((ti->od.align.type & ALIGN_TOGGLE_MERGE_BIT) ? PHARAO_MERGE : ""),
        ti->od.align.template_dir, SEPARATOR,
        ti->od.al.mol_info[template_object_num]->object_id,
        ti->od.align.align_scratch, SEPARATOR, ti->od.al.mol_info[template_object_num]->object_id);
    }
    pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[task_num]->code));
    ext_program_wait(&prog_exe_info, pid);
    if (ti->od.al.task_list[task_num]->code) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      error = 1;
      continue;
    }
    if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
      remove(single_conf_mol_fd.name);
    }
    if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      O3_ERROR_STRING(ti->od.al.task_list[task_num], log_fd.name);
      ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
      error = 1;
      continue;
    }
    if (fgrep(log_fd.handle, buffer, "Error")) {
      O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
      O3_ERROR_STRING(ti->od.al.task_list[task_num], log_fd.name);
      ti->od.al.task_list[task_num]->code = FL_PHARAO_ERROR;
    }
    fclose(log_fd.handle);
    if (ti->od.al.task_list[task_num]->code) {
      error = 1;
      continue;
    }
    remove(log_fd.name);
    if ((ti->od.align.type & ALIGN_MIXED_BIT)
      && (!(ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT))) {
      sprintf(buffer, "%s%c%04d-%04d.sdf",
        ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[0]->object_id,
        ti->od.al.mol_info[ti->od.grid.object_num - 1]->object_id);
      if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
        sprintf(template_conf_string, "_%06d", template_conf_num + 1);
      }
      sprintf(log_fd.name, "%s%c%04d-%04d_on_%04d%s.log",
        ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[0]->object_id,
        ti->od.al.mol_info[ti->od.grid.object_num - 1]->object_id,
        ti->od.al.mol_info[template_object_num]->object_id,
        template_conf_string);
      sprintf(pharao_temp_dir, "%s%c%04d-%04d_on_%04d%s_pharao",
        ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[0]->object_id,
        ti->od.al.mol_info[ti->od.grid.object_num - 1]->object_id,
        ti->od.al.mol_info[template_object_num]->object_id,
        template_conf_string);
      sprintf(pharao_sdf_fd.name, "%s.sdf", pharao_temp_dir);
      sprintf(prog_exe_info.command_line,
        "%s -q %s %s -r %s%c%04d%s.phar --refType PHAR "
        "-d %s --dbType MOL "
        "-s %s%c%04d-%04d_on_%04d%s.scores -o %s",
        ti->od.align.pharao_exe,
        ((ti->od.align.type & ALIGN_TOGGLE_HYBRID_BIT) ? "" : PHARAO_NO_HYBRID),
        ((ti->od.align.type & ALIGN_TOGGLE_MERGE_BIT) ? PHARAO_MERGE : ""),
        ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[template_object_num]->object_id,
        template_conf_string, buffer,
        ti->od.align.align_scratch, SEPARATOR,
        ti->od.al.mol_info[0]->object_id,
        ti->od.al.mol_info[ti->od.grid.object_num - 1]->object_id,
        ti->od.al.mol_info[template_object_num]->object_id,
        template_conf_string, pharao_sdf_fd.name);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[task_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      if (ti->od.al.task_list[task_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        error = 1;
        continue;
      }
      /*
      check if the Pharao computation was OK
      */
      if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], log_fd.name);
        ti->od.al.task_list[task_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
        error = 1;
        continue;
      }
      if (fgrep(log_fd.handle, buffer, "Error")) {
        O3_ERROR_LOCATE(ti->od.al.task_list[task_num]);
        O3_ERROR_STRING(ti->od.al.task_list[task_num], log_fd.name);
        ti->od.al.task_list[task_num]->code = FL_PHARAO_ERROR;
      }
      fclose(log_fd.handle);
      if (ti->od.al.task_list[task_num]->code) {
        error = 1;
        continue;
      }
      if (!dexist(pharao_temp_dir)) {
        #ifndef WIN32
        result = mkdir(pharao_temp_dir, S_IRWXU | S_IRGRP | S_IROTH);
        #else
        result = mkdir(pharao_temp_dir);
        #endif
      }
      if (result) {
        ti->od.al.task_list[task_num]->code = FL_CANNOT_CREATE_SCRDIR;
        error = 1;
        continue;
      }
      if ((ti->od.al.task_list[task_num]->code = break_sdf_to_sdf
        (&(ti->od), ti->od.al.task_list[task_num], &pharao_sdf_fd, pharao_temp_dir))) {
        error = 1;
        continue;
      }
      remove(log_fd.name);
    }
  }
  free_proc_env(prog_exe_info.proc_env);
  
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


#ifndef WIN32
void *align_multi_pharao_thread(void *pointer)
#else
DWORD align_multi_pharao_thread(void *pointer)
#endif
{
  char buffer[BUF_LEN];
  char buffer2[BUF_LEN];
  char template_conf_string[MAX_NAME_LEN];
  int i;
  int template_num;
  int template_object_num;
  int template_conf_num;
  int moved_object_num;
  int best_conf_num;
  int found = 0;
  int assigned = 0;
  int done_array_pos = 0;
  int error;
  int pid;
  double score;
  double best_score = 0.0;
  FileDescriptor log_fd;
  FileDescriptor conf_sdf_fd;
  FileDescriptor out_sdf_fd;
  FileDescriptor pharao_sdf_fd;
  FileDescriptor phar_fd;
  FileDescriptor scores_fd;
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  
  
  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, BUF_LEN);
  memset(buffer2, 0, BUF_LEN);
  memset(template_conf_string, 0, MAX_NAME_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  memset(&conf_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&out_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&pharao_sdf_fd, 0, sizeof(FileDescriptor));
  memset(&phar_fd, 0, sizeof(FileDescriptor));
  memset(&scores_fd, 0, sizeof(FileDescriptor));
  prog_exe_info.exedir = ti->od.align.pharao_exe_path;
  prog_exe_info.proc_env = fill_env(&(ti->od), babel_env, ti->od.align.pharao_exe_path, 0);
  prog_exe_info.stdout_fd = &log_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.sep_proc_grp = 1;
  for (template_num = 0, done_array_pos = 0, error = 0; (!error)
    && (template_num < ti->od.pel.numberlist[OBJECT_LIST]->size); ++template_num) {
    /*
    loop over all templates
    */
    template_object_num = ti->od.pel.numberlist[OBJECT_LIST]->pe[template_num] - 1;
    for (template_conf_num = 0; (!error) && (template_conf_num < ((ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
      ? ti->od.pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1)); ++template_conf_num) {
      if (ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) {
        sprintf(template_conf_string, "_%06d", template_conf_num + 1);
      }
      /*
      loop over all conformations available
      for the current template
      */
      assigned = 1;
      while ((!error) && assigned) {
        moved_object_num = 0;
        assigned = 0;
        while ((!assigned) && (moved_object_num < ti->od.grid.object_num)) {
          if (!(ti->od.al.done_objects[done_array_pos][moved_object_num])) {
            sprintf(buffer, "%s%c%04d%c%04d_on_%04d%s.sdf",
              ti->od.align.align_scratch, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id, SEPARATOR,
              ti->od.al.mol_info[moved_object_num]->object_id,
              ti->od.al.mol_info[template_object_num]->object_id,
              template_conf_string);
            #ifndef WIN32
            pthread_mutex_lock(ti->od.mel.mutex);
            #else
            WaitForSingleObject(ti->od.mel.mutex, INFINITE);
            #endif
            if (!(ti->od.al.done_objects[done_array_pos][moved_object_num])) {
              ti->od.al.done_objects[done_array_pos][moved_object_num] = OBJECT_ASSIGNED;
              assigned = 1;
            }
            #ifndef WIN32
            pthread_mutex_unlock(ti->od.mel.mutex);
            #else
            ReleaseMutex(ti->od.mel.mutex);
            #endif
          }
          if (!assigned) {
            ++moved_object_num;
          }
        }
        if (!assigned)  {
          break;
        }
        ti->od.al.task_list[moved_object_num]->code = 0;
        if (!(prog_exe_info.proc_env)) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          ti->od.al.task_list[moved_object_num]->code = FL_OUT_OF_MEMORY;
          error = 1;
          continue;
        }
        ti->od.al.task_list[moved_object_num]->data[TEMPLATE_OBJECT_NUM] = template_object_num;
        sprintf(buffer, "%s%c%04d", ti->od.align.align_scratch, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id);
        if (!dexist(buffer)) {
          #ifndef WIN32
          error = mkdir(buffer, S_IRWXU | S_IRGRP | S_IROTH);
          #else
          error = mkdir(buffer);
          #endif
        }
        if (error) {
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_CREATE_SCRDIR;
          continue;
        }
        sprintf(out_sdf_fd.name, "%s%c%04d_on_%04d%s.sdf",
          buffer, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id,
          ti->od.al.mol_info[template_object_num]->object_id,
          template_conf_string);
        if (!(out_sdf_fd.handle = fopen(out_sdf_fd.name, "wb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], out_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_WRITE_SDF_FILE;
          error = 1;
          continue;
        }
        sprintf(conf_sdf_fd.name, "%s%c%04d.sdf",
          ti->od.align.candidate_conf_dir, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id);
        ti->od.al.task_list[moved_object_num]->data[TEMPLATE_CONF_NUM] =
          ((ti->od.align.type & ALIGN_MULTICONF_TEMPLATE_BIT) ? template_conf_num : -1);
        /*
        align the conformational database of object "moved_object_num"
        on the current object,conformation pair
        */
        sprintf(log_fd.name, "%s%c%04d_on_%04d%s.log",
          ti->od.align.align_scratch, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id,
          ti->od.al.mol_info[template_object_num]->object_id,
          template_conf_string);
        sprintf(pharao_sdf_fd.name, "%s%c%04d_on_%04d%s_pharao.sdf",
          ti->od.align.align_scratch, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id,
          ti->od.al.mol_info[template_object_num]->object_id,
          template_conf_string);
        sprintf(phar_fd.name, "%s%c%04d%s.phar",
          ti->od.align.align_scratch, SEPARATOR,
          ti->od.al.mol_info[template_object_num]->object_id,
          template_conf_string);
        sprintf(scores_fd.name, "%s%c%04d_on_%04d%s.scores",
          ti->od.align.align_scratch, SEPARATOR,
          ti->od.al.mol_info[moved_object_num]->object_id,
          ti->od.al.mol_info[template_object_num]->object_id,
          template_conf_string);
        sprintf(prog_exe_info.command_line,
          "%s -q %s %s -r %s --refType PHAR "
          "-d %s --dbType MOL -s %s -o %s",
          ti->od.align.pharao_exe,
          ((ti->od.align.type & ALIGN_TOGGLE_HYBRID_BIT) ? "" : PHARAO_NO_HYBRID),
          ((ti->od.align.type & ALIGN_TOGGLE_MERGE_BIT) ? PHARAO_MERGE : ""),
          phar_fd.name, conf_sdf_fd.name,
          scores_fd.name, pharao_sdf_fd.name);
        pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[moved_object_num]->code));
        ext_program_wait(&prog_exe_info, pid);
        if (ti->od.al.task_list[moved_object_num]->code) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        /*
        check if the Pharao computation was OK
        */
        if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], log_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        if (fgrep(log_fd.handle, buffer2, "Error")) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], log_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_PHARAO_ERROR;
          fclose(out_sdf_fd.handle);
          error = 1;
        }
        fclose(log_fd.handle);
        remove(log_fd.name);
        if (ti->od.al.task_list[moved_object_num]->code) {
          continue;
        }
        /*
        find the best scoring conformation for object "object_num"
        */
        if (!(scores_fd.handle = fopen(scores_fd.name, "rb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], scores_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_PHARAO_OUTPUT;
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        best_conf_num = -1;
        i = 0;
        while (fgets(buffer, BUF_LEN, scores_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          sscanf(buffer, "%*s %*s %*s %*s %*s %*s %*s %*s %lf", &score);
          if ((best_conf_num == -1) || (score > best_score)) {
            best_score = score;
            best_conf_num = i;
          }
          ++i;
        }
        fclose(scores_fd.handle);
        remove(scores_fd.name);
        /*
        read Pharao SDF output
        */
        if (!(pharao_sdf_fd.handle = fopen(pharao_sdf_fd.name, "rb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_SDF_FILE;
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        /*
        look for the conformation having the best Tanimoto score
        */
        ti->od.al.task_list[moved_object_num]->code =
          find_conformation_in_sdf(pharao_sdf_fd.handle, NULL, best_conf_num);
        if (ti->od.al.task_list[moved_object_num]->code) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
          fclose(pharao_sdf_fd.handle);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        /*
        open the conformational SDF database
        */
        if (!(conf_sdf_fd.handle = fopen(conf_sdf_fd.name, "rb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], conf_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_READ_CONF_FILE;
          fclose(pharao_sdf_fd.handle);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        ti->od.al.task_list[moved_object_num]->code = find_conformation_in_sdf
          (conf_sdf_fd.handle, out_sdf_fd.handle, best_conf_num);
        if (ti->od.al.task_list[moved_object_num]->code) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], conf_sdf_fd.name);
          fclose(pharao_sdf_fd.handle);
          fclose(conf_sdf_fd.handle);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        i = 0;
        while ((i < ti->od.al.mol_info[moved_object_num]->n_atoms)
          && fgets(buffer, BUF_LEN, conf_sdf_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          if (!fgets(buffer2, BUF_LEN, pharao_sdf_fd.handle)) {
            break;
          }
          buffer2[BUF_LEN - 1] = '\0';
          remove_newline(buffer2);
          fprintf(out_sdf_fd.handle, "%s\n", buffer2);
          ++i;
        }
        if (i < ti->od.al.mol_info[moved_object_num]->n_atoms) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], conf_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_FIND_CONF;
          fclose(pharao_sdf_fd.handle);
          fclose(conf_sdf_fd.handle);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        found = 0;
        while ((!found) && fgets(buffer, BUF_LEN, conf_sdf_fd.handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          if (!(found = (!strncmp(buffer, SDF_DELIMITER, 4)))) {
            fprintf(out_sdf_fd.handle, "%s\n", buffer);
          }
        }
        if (!found) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], conf_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_FIND_CONF;
          fclose(pharao_sdf_fd.handle);
          fclose(conf_sdf_fd.handle);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        found = 0;
        while ((!found) && fgets(buffer2, BUF_LEN, pharao_sdf_fd.handle)) {
          buffer2[BUF_LEN - 1] = '\0';
          remove_newline(buffer2);
          found = (!strncasecmp(buffer2, ">  <PHARAO_TANIMOTO>", 20));
        }
        if (!found) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_FIND_CONF;
          fclose(pharao_sdf_fd.handle);
          fclose(conf_sdf_fd.handle);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        fprintf(out_sdf_fd.handle, "%s\n", buffer2);
        if (!fgets(buffer2, BUF_LEN, pharao_sdf_fd.handle)) {
          O3_ERROR_LOCATE(ti->od.al.task_list[moved_object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[moved_object_num], pharao_sdf_fd.name);
          ti->od.al.task_list[moved_object_num]->code = FL_CANNOT_FIND_CONF;
          fclose(pharao_sdf_fd.handle);
          fclose(conf_sdf_fd.handle);
          fclose(out_sdf_fd.handle);
          error = 1;
          continue;
        }
        buffer2[BUF_LEN - 1] = '\0';
        fprintf(out_sdf_fd.handle, "%s\n", buffer2);
        if (ti->od.align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
          fprintf(out_sdf_fd.handle, ">  <BEST_CANDIDATE_CONF>\n%d\n\n", best_conf_num + 1);
        }
        fprintf(out_sdf_fd.handle, SDF_DELIMITER"\n");
        fclose(pharao_sdf_fd.handle);
        remove(pharao_sdf_fd.name);
        fclose(conf_sdf_fd.handle);
        fclose(out_sdf_fd.handle);
        #ifndef WIN32
        pthread_mutex_lock(ti->od.mel.mutex);
        #else
        WaitForSingleObject(ti->od.mel.mutex, INFINITE);
        #endif
        ti->od.al.done_objects[done_array_pos][moved_object_num] |= OBJECT_FINISHED;
        #ifndef WIN32
        pthread_mutex_unlock(ti->od.mel.mutex);
        #else
        ReleaseMutex(ti->od.mel.mutex);
        #endif
      }
      i = join_aligned_files(&(ti->od), done_array_pos, buffer);
      if (i != -1) {
        O3_ERROR_LOCATE(ti->od.al.task_list[i]);
        O3_ERROR_STRING(ti->od.al.task_list[i], buffer);
        ti->od.al.task_list[i]->code = FL_CANNOT_READ_SDF_FILE;
        error = 1;
      }
      else {
        ++done_array_pos;
      }
    }
  }
  free_proc_env(prog_exe_info.proc_env);
  
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}
