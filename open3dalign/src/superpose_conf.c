/*

superpose_conf.c

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
#include <include/align.h>


#define WORK_SIZE    BUF_LEN
#define IWORK_SIZE    BUF_LEN
#define EIGENVECTORS    'V'
#define UPPER_DIAG    'U'


void calc_conf_centroid(ConfInfo *conf, double *centroid)
{
  int i;
  int x;
  int n_heavy_atoms;
  
  
  memset(centroid, 0, 3 * sizeof(double));
  for (x = 0; x < 3; ++x) {
    for (i = 0, n_heavy_atoms = 0; i < conf->n_atoms; ++i) {
      if (strcmp(conf->atom[i]->element, "H")) {
        centroid[x] += conf->coord[i * 3 + x];
        ++n_heavy_atoms;
      }
    }
    centroid[x] /= (double)n_heavy_atoms;
  }
}


int superpose_conf_lap(LAPInfo *li, ConfInfo *moved_conf, ConfInfo *template_conf, ConfInfo *fitted_conf,
  ConfInfo *progress_conf, AtomPair *temp_sdm, AtomPair *fitted_sdm, char **used, double *rt_mat,
  double *heavy_msd, double *original_heavy_msd, int *pairs)
{
  int i;
  int flag;
  int iter;
  int pairs_temp;
  int result;
  double fitted_rt_mat[RT_MAT_SIZE];
  double template_centroid[3];
  double moved_centroid[3];
  double pairs_heavy_msd;
  double pairs_heavy_msd_temp;
  
  
  compute_cost_matrix(li, moved_conf, template_conf, MAX_H_BINS, 0, MATCH_CONFORMERS_BIT);
  lap(li, template_conf->n_heavy_atoms);
  calc_conf_centroid(template_conf, template_centroid);
  calc_conf_centroid(moved_conf, moved_centroid);
  cblas_dcopy(moved_conf->n_atoms * 3, moved_conf->coord, 1, fitted_conf->coord, 1);
  for (i = 0; i < moved_conf->n_atoms; ++i) {
    cblas_daxpy(3, -1.0, moved_centroid, 1, &(fitted_conf->coord[i * 3]), 1);
    cblas_daxpy(3, 1.0, template_centroid, 1, &(fitted_conf->coord[i * 3]), 1);
  }
  *pairs = filter_sol_vector(li, fitted_conf, template_conf, temp_sdm, fitted_sdm);
  if (*pairs < 3) {
    for (i = 0, *pairs = 0; i < template_conf->n_atoms; ++i) {
      if (!strcmp(template_conf->atom[i]->element, "H")) {
        continue;
      }
      fitted_sdm[i].a[0] = i;
      fitted_sdm[i].a[1] = i;
      ++(*pairs);
    }
  }
  result = rms_algorithm(DONT_USE_WEIGHTS, fitted_sdm, *pairs, moved_conf,
    template_conf, fitted_conf, fitted_rt_mat, heavy_msd, NULL);
  if (result) {
    *heavy_msd = MAX_CUTOFF;
    if (original_heavy_msd) {
      *original_heavy_msd = MAX_CUTOFF;
    }
    *pairs = 0;
    return result;
  }
  *pairs = 0;
  pairs_heavy_msd = MAX_CUTOFF;
  flag = 1;
  iter = 0;
  cblas_dcopy(template_conf->n_atoms * 3, fitted_conf->coord, 1, progress_conf->coord, 1);
  while (flag && (iter < MAX_SDM_ITERATIONS)) {
    ++iter;
    /*
    call sdm_algorithm
    */
    pairs_temp = sdm_algorithm(temp_sdm, progress_conf, template_conf,
      used, MATCH_ATOM_TYPES_BIT | CENTER_TO_ORIGIN_BIT, SDM_THRESHOLD_START);
    if (pairs_temp < 3) {
      break;
    }
    /*
    call rms_algorithm
    */
    if (rms_algorithm(DONT_USE_WEIGHTS, temp_sdm, pairs_temp, progress_conf,
      template_conf, fitted_conf, fitted_rt_mat, &pairs_heavy_msd_temp, NULL)) {
      break;
    }
    /*
    keep looping until:
    1) it is possible to increase the number of fitted pairs
    2) it is not possible to increase the number of fitted pairs
       anymore, but the msd is improved compared to the previous one
    */
    flag = ((pairs_temp > *pairs) || ((pairs_temp == *pairs)
      && ((pairs_heavy_msd - pairs_heavy_msd_temp) > MSD_THRESHOLD)));
    if (flag) {
      *pairs = pairs_temp;
      pairs_heavy_msd = pairs_heavy_msd_temp;
      memcpy(fitted_sdm, temp_sdm, *pairs * sizeof(AtomPair));
      cblas_dcopy(template_conf->n_atoms * 3, fitted_conf->coord, 1, progress_conf->coord, 1);
    }
  }
  result = rms_algorithm(DONT_USE_WEIGHTS, fitted_sdm, *pairs, moved_conf,
    template_conf, fitted_conf, fitted_rt_mat, &pairs_heavy_msd, NULL);
  overall_msd(fitted_sdm, *pairs, fitted_conf, template_conf, heavy_msd);
  if (original_heavy_msd) {
    overall_msd(fitted_sdm, *pairs, moved_conf, template_conf, original_heavy_msd);
  }
  if (rt_mat) {
    memcpy(rt_mat, fitted_rt_mat, RT_MAT_SIZE * sizeof(double));
  }
  
  return result;
}


int superpose_conf_syst(ConfInfo *moved_conf, ConfInfo *template_conf, ConfInfo *fitted_conf,
  ConfInfo *progress_conf, ConfInfo *cand_conf, AtomPair *sdm, AtomPair *local_best_sdm,
  AtomPair *fitted_sdm, char **used, double *rt_mat, int angle_step,
  double *heavy_msd, double *original_heavy_msd, int *pairs)
{
  int i = 0;
  int pairs1 = 0;
  int pairs2 = 0;
  int iter = 0;
  int flag;
  int result;
  int angle[3];
  double pairs_heavy_msd = 0.0;
  double pairs_heavy_msd1 = 0.0;
  double pairs_heavy_msd2 = 0.0;
  double rt_mat2[RT_MAT_SIZE];
  double rt_mat_sys[RT_MAT_SIZE];
  double fitted_rt_mat[RT_MAT_SIZE];
  double t_mat1[RT_MAT_SIZE];
  double t_mat2[RT_MAT_SIZE];
  double t_vec1[RT_VEC_SIZE];
  double t_vec2[RT_VEC_SIZE];
  double moved_centroid[3];
  double template_centroid[3];
  double rad[3];
  
  
  *pairs = 0;
  /*
  find the centroids of the reference
  and candidate structures
  */
  calc_conf_centroid(moved_conf, moved_centroid);
  calc_conf_centroid(template_conf, template_centroid);
  /*
  now we systematically rotate moved_conf by angle_step
  increments around the three axes; each rotation will constitute
  a new starting point for the superimposition procedure
  */
  memset(t_mat1, 0, RT_MAT_SIZE * sizeof(double));
  for (i = 0; i < RT_MAT_SIZE; i += 5) {
    t_mat1[i] = 1.0;
  }
  memcpy(t_mat2, t_mat1, RT_MAT_SIZE * sizeof(double));
  cblas_daxpy(3, -1.0, moved_centroid, 1, &t_mat1[3 * RT_VEC_SIZE], 1);
  cblas_dcopy(3, template_centroid, 1, &t_mat2[3 * RT_VEC_SIZE], 1);
  memset(t_vec1, 0, RT_VEC_SIZE * sizeof(double));
  t_vec1[3] = 1.0;
  pairs_heavy_msd = MAX_CUTOFF;
  flag = 1;
  for (angle[0] = 0; angle[0] < 360; angle[0] += angle_step) {
    rad[0] = angle2rad(angle[0]);
    for (angle[1] = 0; angle[1] < 360; angle[1] += angle_step) {
      rad[1] = angle2rad(angle[1]);
      for (angle[2] = 0; angle[2] < 180; angle[2] += angle_step) {
        rad[2] = angle2rad(angle[2]);
        prepare_rototrans_matrix(rt_mat_sys, t_mat1, t_mat2, rad);
        for (i = 0; i < moved_conf->n_atoms; ++i) {
          cblas_dcopy(3, &(moved_conf->coord[i * 3]), 1, t_vec1, 1);
          cblas_dgemv(CblasColMajor, CblasNoTrans,
            RT_VEC_SIZE, RT_VEC_SIZE, 1.0, rt_mat_sys, RT_VEC_SIZE,
            t_vec1, 1, 0.0, t_vec2, 1);
          cblas_dcopy(3, t_vec2, 1, &(cand_conf->coord[i * 3]), 1);
        }
        pairs1 = 0;
        pairs_heavy_msd1 = MAX_CUTOFF;
        flag = 1;
        iter = 0;
        while (flag && (iter < MAX_SDM_ITERATIONS)) {
          ++iter;
          /*
          call sdm_algorithm
          */
          pairs2 = sdm_algorithm(sdm, cand_conf, template_conf,
            used, MATCH_ATOM_TYPES_BIT | CENTER_TO_ORIGIN_BIT, 2.0);
          if (pairs2 < 3) {
            break;
          }
          /*
          call rms_algorithm
          */
          result = rms_algorithm(DONT_USE_WEIGHTS, sdm, pairs2, cand_conf,
            template_conf, progress_conf, rt_mat2, &pairs_heavy_msd2, NULL);
          if (result) {
            return result;
          }
          /*
          keep looping until:
          1) it is possible to increase the number of fitted pairs
          2) it is not possible to increase the number of fitted pairs
             anymore, but the msd is improved compared to the previous one
          */
          flag = ((pairs2 > pairs1) || ((pairs2 == pairs1)
            && ((pairs_heavy_msd1 - pairs_heavy_msd2) > MSD_THRESHOLD)));
          if (flag) {
            pairs1 = pairs2;
            pairs_heavy_msd1 = pairs_heavy_msd2;
            memcpy(local_best_sdm, sdm, pairs1 * sizeof(AtomPair));
            cblas_dcopy(cand_conf->n_atoms * 3, progress_conf->coord, 1, cand_conf->coord, 1);
          }
        }
        if ((pairs1 > *pairs) || ((pairs1 == *pairs)
          && ((pairs_heavy_msd - pairs_heavy_msd1) > MSD_THRESHOLD))) {
          *pairs = pairs1;
          pairs_heavy_msd = pairs_heavy_msd1;
          memcpy(fitted_sdm, local_best_sdm, *pairs * sizeof(AtomPair));
        }
      }
    }
  }
  result = rms_algorithm(DONT_USE_WEIGHTS, fitted_sdm, *pairs, moved_conf,
    template_conf, fitted_conf, fitted_rt_mat, &pairs_heavy_msd, NULL);
  if (result) {
    return result;
  }
  overall_msd(fitted_sdm, *pairs, fitted_conf, template_conf, heavy_msd);
  if (original_heavy_msd) {
    overall_msd(fitted_sdm, *pairs, moved_conf, template_conf, original_heavy_msd);
  }
  if (rt_mat) {
    memcpy(rt_mat, fitted_rt_mat, RT_MAT_SIZE * sizeof(double));
  }
  
  return 0;
}


int rms_algorithm(int options, AtomPair *sdm, int pairs, ConfInfo *moved_conf,
  ConfInfo *template_conf, ConfInfo *fitted_conf, double *rt_mat,
  double *heavy_msd, double *original_heavy_msd)
{
  extern char mmff_corr_matrix[99][99];
  char corr;
  char jobz = EIGENVECTORS;
  char uplo = UPPER_DIAG;
  int i;
  int x;
  int n;
  int info = 0;
  #ifndef HAVE_LIBSUNPERF
  int lwork = WORK_SIZE;
  double work[WORK_SIZE];
  #endif
  double moved_pairs_centroid[3];
  double template_pairs_centroid[3];
  double *ev;
  double m[3];
  double p[3];
  double d[4];
  double rt_mat1[RT_MAT_SIZE];
  double rt_mat2[RT_MAT_SIZE];
  double t_mat1[RT_MAT_SIZE];
  double t_mat2[RT_MAT_SIZE];
  double t_vec1[RT_VEC_SIZE];
  double t_vec2[RT_VEC_SIZE];
  double z[16];
  double weight = 1.0;
  
  
  #ifndef HAVE_LIBSUNPERF
  memset(work, 0, WORK_SIZE);
  #endif
  memset(d, 0, 4 * sizeof(double));
  memset(z, 0, 16 * sizeof(double));
  memset(moved_pairs_centroid, 0, 3 * sizeof(double));
  memset(template_pairs_centroid, 0, 3 * sizeof(double));
  /*
  compute centroids taking into account only
  atom pairs used for the superposition
  */
  for (x = 0; x < 3; ++x) {
    for (i = 0; i < pairs; ++i) {
      template_pairs_centroid[x] += template_conf->coord[sdm[i].a[0] * 3 + x];
      moved_pairs_centroid[x] += moved_conf->coord[sdm[i].a[1] * 3 + x];
    }
    template_pairs_centroid[x] /= (double)pairs;
    moved_pairs_centroid[x] /= (double)pairs;
  }
  memset(t_vec1, 0, RT_VEC_SIZE * sizeof(double));
  t_vec1[3] = 1.0;
  memset(t_mat1, 0, RT_MAT_SIZE * sizeof(double));
  for (i = 0; i < RT_MAT_SIZE; i += 5) {
    t_mat1[i] = 1.0;
  }
  memcpy(t_mat2, t_mat1, RT_MAT_SIZE * sizeof(double));
  cblas_daxpy(3, -1.0, moved_pairs_centroid, 1, &t_mat1[3 * RT_VEC_SIZE], 1);
  cblas_dcopy(3, template_pairs_centroid, 1, &t_mat2[3 * RT_VEC_SIZE], 1);
  /*
  find the best rotation matrix
  */
  for (i = 0; i < pairs; ++i) {
    weight = 1.0;
    if (options & USE_CHARGE_WEIGHTS) {
      weight = (1.0 + O3_CHARGE_COEFF * fabs(template_conf->atom[sdm[i].a[0]]->charge + moved_conf->atom[sdm[i].a[1]]->charge))
        / (1.0 + fabs(template_conf->atom[sdm[i].a[0]]->charge - moved_conf->atom[sdm[i].a[1]]->charge));
    }
    else if (options & USE_MMFF_WEIGHTS) {
      corr = mmff_corr_matrix[template_conf->atom[sdm[i].a[0]]->atom_type - 1]
        [moved_conf->atom[sdm[i].a[1]]->atom_type - 1];
      weight = 1.0 / (1.0 + (double)(corr - 1));
    }
    for (x = 0; x < 3; ++x) {
      m[x] = (template_conf->coord[sdm[i].a[0] * 3 + x] - template_pairs_centroid[x])
        - (moved_conf->coord[sdm[i].a[1] * 3 + x] - moved_pairs_centroid[x]);
      p[x] = (template_conf->coord[sdm[i].a[0] * 3 + x] - template_pairs_centroid[x])
        + (moved_conf->coord[sdm[i].a[1] * 3 + x] - moved_pairs_centroid[x]);
    }
    z[0]  += weight * (square(m[0]) + square(m[1]) + square(m[2]));
    z[4]  += weight * (p[1] * m[2] - m[1] * p[2]);
    z[5]  += weight * (square(p[1]) + square(p[2]) + square(m[0]));
    z[8]  += weight * (m[0] * p[2] - p[0] * m[2]);
    z[9]  += weight * (m[0] * m[1] - p[0] * p[1]);
    z[10] += weight * (square(p[0]) + square(p[2]) + square(m[1]));
    z[12] += weight * (p[0] * m[1] - m[0] * p[1]);
    z[13] += weight * (m[0] * m[2] - p[0] * p[2]);
    z[14] += weight * (m[1] * m[2] - p[1] * p[2]);
    z[15] += weight * (square(p[0]) + square(p[1]) + square(m[2]));
  }
  n = 4;
  #ifdef HAVE_LIBMKL
  dsyev(&jobz, &uplo, &n, z, &n, d, work, &lwork, &info);
  #elif HAVE_LIBSUNPERF
  dsyev(jobz, uplo, n, z, n, d, &info);
  #elif HAVE_LIBACCELERATE
  dsyev_(&jobz, &uplo, &n, z, &n, d, work, &lwork, &info);
  #elif HAVE_LIBATLAS
  dsyev_(&jobz, &uplo, &n, z, &n, d, work, &lwork, &info);
  #endif
  if (info) {
    return FL_ABNORMAL_TERMINATION;
  }
  ev = z;
  if (!rt_mat) {
    rt_mat = rt_mat1;
  }
  memset(rt_mat, 0, RT_MAT_SIZE * sizeof(double));
  rt_mat[0]                   = square(ev[0]) + square(ev[1]) - square(ev[2]) - square(ev[3]);
  rt_mat[1]                   = 2.0 * (ev[1] * ev[2] - ev[0] * ev[3]);
  rt_mat[2]                   = 2.0 * (ev[1] * ev[3] + ev[0] * ev[2]);
  rt_mat[RT_VEC_SIZE]         = 2.0 * (ev[1] * ev[2] + ev[0] * ev[3]);
  rt_mat[RT_VEC_SIZE + 1]     = square(ev[0]) + square(ev[2]) - square(ev[1]) - square(ev[3]);
  rt_mat[RT_VEC_SIZE + 2]     = 2.0 * (ev[2] * ev[3] - ev[0] * ev[1]);
  rt_mat[RT_VEC_SIZE * 2]     = 2.0 * (ev[1] * ev[3] - ev[0] * ev[2]);
  rt_mat[RT_VEC_SIZE * 2 + 1] = 2.0 * (ev[2] * ev[3] + ev[0] * ev[1]);
  rt_mat[RT_VEC_SIZE * 2 + 2] = square(ev[0]) + square(ev[3]) - square(ev[1]) - square(ev[2]);
  rt_mat[RT_VEC_SIZE * 3 + 3] = 1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, t_mat2, RT_VEC_SIZE, rt_mat, RT_VEC_SIZE,
    0.0, rt_mat2, RT_VEC_SIZE);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, rt_mat2, RT_VEC_SIZE, t_mat1, RT_VEC_SIZE,
    0.0, rt_mat, RT_VEC_SIZE);
  if (original_heavy_msd) {
    *original_heavy_msd = 0.0;
    for (i = 0; i < pairs; ++i) {
      *original_heavy_msd += squared_euclidean_distance
        (&(template_conf->coord[sdm[i].a[0] * 3]),
        &(moved_conf->coord[sdm[i].a[1] * 3]));
    }
    *original_heavy_msd /= (double)pairs;
  }
  for (i = 0; i < moved_conf->n_atoms; ++i) {
    /*
    apply the rotation/translation matrix to moved_conf coordinates
    */
    cblas_dcopy(3, &(moved_conf->coord[i * 3]), 1, t_vec1, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans,
      RT_VEC_SIZE, RT_VEC_SIZE, 1.0, rt_mat, RT_VEC_SIZE,
      t_vec1, 1, 0.0, t_vec2, 1);
    cblas_dcopy(3, t_vec2, 1, &(fitted_conf->coord[i * 3]), 1);
  }
  *heavy_msd = safe_rint(d[0] / (double)pairs * 1.0e06) / 1.0e06;
  
  return 0;
}


int rms_algorithm_multi(O3Data *od, O3Data *od_comp, double *rt_mat, double *heavy_msd)
{
  char jobz = EIGENVECTORS;
  char uplo = UPPER_DIAG;
  int i;
  int x;
  int n;
  int object_num;
  int pairs;
  int overall_pairs;
  int info = 0;
  #ifndef HAVE_LIBSUNPERF
  int lwork = WORK_SIZE;
  double work[WORK_SIZE];
  #endif
  double *ev;
  double rt_mat2[RT_MAT_SIZE];
  double t_mat1[RT_MAT_SIZE];
  double t_mat2[RT_MAT_SIZE];
  double moved_pairs_centroid[3];
  double template_pairs_centroid[3];
  double m[3];
  double p[3];
  double d[4];
  double z[16];
  AtomPair *sdm = NULL;
  
  
  #ifndef HAVE_LIBSUNPERF
  memset(work, 0, WORK_SIZE);
  #endif
  memset(d, 0, 4 * sizeof(double));
  memset(z, 0, 16 * sizeof(double));
  memset(moved_pairs_centroid, 0, 3 * sizeof(double));
  memset(template_pairs_centroid, 0, 3 * sizeof(double));
  for (object_num = 0, overall_pairs = 0; object_num < od->grid.object_num; ++object_num) {
    pairs = od->al.rt_list[object_num]->pairs;
    sdm = od->al.rt_list[object_num]->sdm;
    /*
    compute centroids taking into account only
    atom pairs used for the superposition
    */
    for (x = 0; x < 3; ++x) {
      for (i = 0; i < pairs; ++i) {
        template_pairs_centroid[x] +=
          od->al.mol_info[object_num]->atom[sdm[i].a[0]]->coord[x];
        moved_pairs_centroid[x] +=
          od_comp->al.mol_info[object_num]->atom[sdm[i].a[1]]->coord[x];
      }
    }
    overall_pairs += pairs;
  }
  for (x = 0; x < 3; ++x) {
    template_pairs_centroid[x] /= (double)overall_pairs;
    moved_pairs_centroid[x] /= (double)overall_pairs;
  }
  memset(t_mat1, 0, RT_MAT_SIZE * sizeof(double));
  for (i = 0; i < RT_MAT_SIZE; i += 5) {
    t_mat1[i] = 1.0;
  }
  memcpy(t_mat2, t_mat1, RT_MAT_SIZE * sizeof(double));
  cblas_daxpy(3, -1.0, moved_pairs_centroid, 1, &t_mat1[3 * RT_VEC_SIZE], 1);
  cblas_dcopy(3, template_pairs_centroid, 1, &t_mat2[3 * RT_VEC_SIZE], 1);
  /*
  find the best rotation matrix
  */
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    pairs = od->al.rt_list[object_num]->pairs;
    sdm = od->al.rt_list[object_num]->sdm;
    for (i = 0; i < pairs; ++i) {
      for (x = 0; x < 3; ++x) {
        m[x] = (od->al.mol_info[object_num]->atom[sdm[i].a[0]]->coord[x] - template_pairs_centroid[x])
          - (od_comp->al.mol_info[object_num]->atom[sdm[i].a[1]]->coord[x] - moved_pairs_centroid[x]);
        p[x] = (od->al.mol_info[object_num]->atom[sdm[i].a[0]]->coord[x] - template_pairs_centroid[x])
          + (od_comp->al.mol_info[object_num]->atom[sdm[i].a[1]]->coord[x] - moved_pairs_centroid[x]);
      }
      z[0]  += (square(m[0]) + square(m[1]) + square(m[2]));
      z[4]  += (p[1] * m[2] - m[1] * p[2]);
      z[5]  += (square(p[1]) + square(p[2]) + square(m[0]));
      z[8]  += (m[0] * p[2] - p[0] * m[2]);
      z[9]  += (m[0] * m[1] - p[0] * p[1]);
      z[10] += (square(p[0]) + square(p[2]) + square(m[1]));
      z[12] += (p[0] * m[1] - m[0] * p[1]);
      z[13] += (m[0] * m[2] - p[0] * p[2]);
      z[14] += (m[1] * m[2] - p[1] * p[2]);
      z[15] += (square(p[0]) + square(p[1]) + square(m[2]));
    }
  }
  n = 4;
  #ifdef HAVE_LIBMKL
  dsyev(&jobz, &uplo, &n, z, &n, d, work, &lwork, &info);
  #elif HAVE_LIBSUNPERF
  dsyev(jobz, uplo, n, z, n, d, &info);
  #elif HAVE_LIBACCELERATE
  dsyev_(&jobz, &uplo, &n, z, &n, d, work, &lwork, &info);
  #elif HAVE_LIBATLAS
  dsyev_(&jobz, &uplo, &n, z, &n, d, work, &lwork, &info);
  #endif
  if (info) {
    return FL_ABNORMAL_TERMINATION;
  }
  ev = z;
  memset(rt_mat, 0, RT_MAT_SIZE * sizeof(double));
  rt_mat[0]                   = square(ev[0]) + square(ev[1]) - square(ev[2]) - square(ev[3]);
  rt_mat[1]                   = 2.0 * (ev[1] * ev[2] - ev[0] * ev[3]);
  rt_mat[2]                   = 2.0 * (ev[1] * ev[3] + ev[0] * ev[2]);
  rt_mat[RT_VEC_SIZE]         = 2.0 * (ev[1] * ev[2] + ev[0] * ev[3]);
  rt_mat[RT_VEC_SIZE + 1]     = square(ev[0]) + square(ev[2]) - square(ev[1]) - square(ev[3]);
  rt_mat[RT_VEC_SIZE + 2]     = 2.0 * (ev[2] * ev[3] - ev[0] * ev[1]);
  rt_mat[RT_VEC_SIZE * 2]     = 2.0 * (ev[1] * ev[3] - ev[0] * ev[2]);
  rt_mat[RT_VEC_SIZE * 2 + 1] = 2.0 * (ev[2] * ev[3] + ev[0] * ev[1]);
  rt_mat[RT_VEC_SIZE * 2 + 2] = square(ev[0]) + square(ev[3]) - square(ev[1]) - square(ev[2]);
  rt_mat[RT_VEC_SIZE * 3 + 3] = 1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, t_mat2, RT_VEC_SIZE, rt_mat, RT_VEC_SIZE,
    0.0, rt_mat2, RT_VEC_SIZE);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, rt_mat2, RT_VEC_SIZE, t_mat1, RT_VEC_SIZE,
    0.0, rt_mat, RT_VEC_SIZE);
  *heavy_msd = safe_rint(d[0] / (double)overall_pairs * 1.0e06) / 1.0e06;
  
  return 0;
}


int sdm_algorithm(AtomPair *sdm, ConfInfo *moved_conf, ConfInfo *template_conf, char **used, int options, double threshold)
{
  int i;
  int j;
  int x;
  int n = 0;
  int pairs = 0;
  int largest_n_atoms;
  double dist = 0.0;
  double template_centroid[3];
  double moved_centroid[3];
  double coord1[3];
  double coord2[3];
  
  
  if (options & CENTER_TO_ORIGIN_BIT) {
    calc_conf_centroid(template_conf, template_centroid);
    calc_conf_centroid(moved_conf, moved_centroid);
  }
  largest_n_atoms = ((moved_conf->n_atoms > template_conf->n_atoms)
    ? moved_conf->n_atoms : template_conf->n_atoms);
  memset(used[0], 0, largest_n_atoms);
  memset(used[1], 0, largest_n_atoms);
  /*
  loop over template_conf atoms
  */
  for (i = 0; i < template_conf->n_atoms; ++i) {
    /*
    skip hydrogens
    */
    if (!strcmp(template_conf->atom[i]->element, "H")) {
      continue;
    }
    /*
    loop over moved_conf atoms
    */
    for (j = 0; j < moved_conf->n_atoms; ++j) {
      /*
      skip hydrogens
      */
      if ((!strcmp(moved_conf->atom[j]->element, "H"))
        || ((options & MATCH_ATOM_TYPES_BIT) && ((template_conf->atom[i]->atom_type !=
        moved_conf->atom[j]->atom_type) || (template_conf->atom[i]->charge !=
        moved_conf->atom[j]->charge)))) {
        continue;
      }
      for (x = 0; x < 3; ++x) {
        coord1[x] = (double)(template_conf->coord[i * 3 + x])
          - ((options & CENTER_TO_ORIGIN_BIT) ? template_centroid[x] : 0.0);
        coord2[x] = (double)(moved_conf->coord[j * 3 + x])
          - ((options & CENTER_TO_ORIGIN_BIT) ? moved_centroid[x] : 0.0);
      }
      dist = squared_euclidean_distance(coord1, coord2);
      /*
      if the distance between these two atoms is lower
      than threshold, then include this
      pair in the SDM matrix
      */
      if (dist < square(threshold)) {
        sdm[n].a[0] = i;
        sdm[n].a[1] = j;
        sdm[n].dist = dist;
        ++n;
      }
    }
  }
  /*
  sort the SDM matrix by increasing distances
  */
  qsort(sdm, n, sizeof(AtomPair), compare_dist);
  /*
  increase the number of pairs which will be used
  by the rms_algorithm until an atom which has been
  included in a previous pair is found
  */
  for (i = 0; i < n; ++i) {
    if (used[0][sdm[i].a[0]] || used[1][sdm[i].a[1]]) {
      break;
    }
    ++pairs;
    used[0][sdm[i].a[0]] = 1;
    used[1][sdm[i].a[1]] = 1;
  }
  
  return pairs;
}


void overall_msd(AtomPair *sdm, int pairs, ConfInfo *moved_conf, ConfInfo *template_conf, double *heavy_msd)
{
  int i;
  int j;
  int x;
  int temp_match = 0;
  int n_heavy_atoms;
  double coord1[3];
  double coord2[3];
  double dist = 0.0;
  double min_dist = 0.0;
  

  for (i = 0; i < template_conf->n_atoms; ++i) {
    template_conf->atom[i]->match = -1;
  }
  for (i = 0; i < pairs; ++i) {
    template_conf->atom[sdm[i].a[0]]->match = sdm[i].a[1];
  }
  /*
  find best matches
  */
  for (i = 0; i < template_conf->n_atoms; ++i) {
    if (template_conf->atom[i]->match != -1) {
      continue;
    }
    min_dist = MAX_CUTOFF;
    cblas_dcopy(3, &(template_conf->coord[i * 3]), 1, coord1, 1);
    for (j = 0; j < moved_conf->n_atoms; ++j) {
      if ((template_conf->atom[i]->tinker_type != moved_conf->atom[j]->tinker_type)
        || strcmp(template_conf->atom[i]->element, template_conf->atom[j]->element)) {
        continue;
      }
      x = 0;
      while ((x != -1) && (x < template_conf->n_atoms)) {
        if (template_conf->atom[x]->match == j) {
          x = -1;
        }
        else {
          ++x;
        }
      }
      if (x == -1) {
        continue;
      }
      cblas_dcopy(3, &(moved_conf->coord[j * 3]), 1, coord2, 1);
      dist = squared_euclidean_distance(coord1, coord2);
      if ((min_dist - dist) > MSD_THRESHOLD) {
        min_dist = dist;
        temp_match = j;
      }
    }
    template_conf->atom[i]->match = temp_match;
  }
  *heavy_msd = 0.0;
  for (i = 0, n_heavy_atoms = 0; i < template_conf->n_atoms; ++i) {
    if (!strcmp(template_conf->atom[i]->element, "H")) {
      continue;
    }
    j = template_conf->atom[i]->match;
    cblas_dcopy(3, &(template_conf->coord[i * 3]), 1, coord1, 1);
    cblas_dcopy(3, &(moved_conf->coord[j * 3]), 1, coord2, 1);
    dist = squared_euclidean_distance(coord1, coord2);
    *heavy_msd += dist;
    ++n_heavy_atoms;
  }
  *heavy_msd = safe_rint(*heavy_msd / (double)n_heavy_atoms * 1.0e06) / 1.0e06;
}


void compute_conf_h(ConfInfo *conf)
{
  int i;
  int y;
  int j;
  int dist;
  
  
  for (i = 0, y = 0; i < conf->n_atoms; ++i) {
    if (!strcmp(conf->atom[i]->element, "H")) {
      continue;
    }
    memset(conf->h[y], 0, MAX_H_BINS * sizeof(int));
    for (j = 0; j < conf->n_atoms; ++j) {
      dist = (int)sqrt(squared_euclidean_distance(&(conf->coord[i * 3]), &(conf->coord[j * 3])));
      if (dist < MAX_H_BINS) {
        ++(conf->h[y][dist]);
      }
    }
    ++y;
  }
}


void lap(LAPInfo *li, int dim)
{
  /*
  input:
  dim        - problem size
  cost       - cost matrix

  output:
  rowsol     - column assigned to row in solution
  colsol     - row assigned to column in solution
  v          - dual variables, column reduction numbers
  */

  int loopcnt;
  int unassignedfound;
  /*
  row vars
  */
  int i;
  int imin;
  int numfree = 0;
  int prvnumfree;
  int f;
  int i0;
  int k;
  int freerow;
  /*
  col vars
  */
  int j;
  int j1;
  int j2 = 0;
  int endofpath;
  int last = 0;
  int low;
  int up;
  /*
  cost vars
  */
  int min = 0;
  int h;
  int umin;
  int usubmin;
  int v2;

  /*
  init how many times a row will be assigned in the column reduction
  */
  memset(li->array[O3_LI_MATCHES], 0, dim * sizeof(int));
  /*
  COLUMN REDUCTION
  reverse order gives better results
  */
  for (j = dim - 1; j >= 0; --j) {
    /*
    find minimum cost over rows
    */
    min = li->cost[0][j];
    imin = 0;
    for (i = 1; i < dim; ++i) {
      if (li->cost[i][j] < min) {
        min = li->cost[i][j];
        imin = i;
      }
    }
    li->array[O3_LI_V][j] = min;

    if (++(li->array[O3_LI_MATCHES][imin]) == 1) {
      /*
      init assignment if minimum row assigned for first time
      */
      li->array[O3_LI_ROWSOL][imin] = j;
      li->array[O3_LI_COLSOL][j] = imin;
    }
    else {
      /*
      row already assigned, column not assigned
      */
      li->array[O3_LI_COLSOL][j] = -1;
    }
  }
  /*
  REDUCTION TRANSFER
  */
  for (i = 0; i < dim; ++i) {
    if (!(li->array[O3_LI_MATCHES][i])) {
      /*
      fill list of unassigned 'free' rows
      */
      li->array[O3_LI_FREE][numfree++] = i;
    }
    else {
      if (li->array[O3_LI_MATCHES][i] == 1) {
        /*
        transfer reduction from rows that are assigned once
        */
        j1 = li->array[O3_LI_ROWSOL][i];
        min = DUMMY_COST;
        for (j = 0; j < dim; ++j) {
          if ((j != j1) && (li->cost[i][j] - li->array[O3_LI_V][j] < min)) {
            min = li->cost[i][j] - li->array[O3_LI_V][j];
          }
        }
        li->array[O3_LI_V][j1] -= min;
      }
    }
  }
  /*
  AUGMENTING ROW REDUCTION
  do-loop to be done twice
  */
  for (loopcnt = 0; loopcnt < 2; ++loopcnt) {
    /*
    scan all free rows
    in some cases, a free row may be replaced with another one to be scanned next
    */
    k = 0;
    prvnumfree = numfree;
    numfree = 0;
    /*
    start list of rows still free after augmenting row reduction
    */
    while (k < prvnumfree) {
      i = li->array[O3_LI_FREE][k];
      k++;
      /*
      find minimum and second minimum reduced cost over columns
      */
      umin = li->cost[i][0] - li->array[O3_LI_V][0];
      j1 = 0;
      usubmin = DUMMY_COST;
      for (j = 1; j < dim; ++j) {
        h = li->cost[i][j] - li->array[O3_LI_V][j];
        if (h < usubmin) {
          if (h >= umin) {
            usubmin = h;
            j2 = j;
          }
          else {
            usubmin = umin;
            umin = h;
            j2 = j1;
            j1 = j;
          }
        }
      }
      i0 = li->array[O3_LI_COLSOL][j1];
      if (umin < usubmin) {
        /*
        change the reduction of the minimum column to increase the minimum
        reduced cost in the row to the subminimum
        */
        li->array[O3_LI_V][j1] -= (usubmin - umin);
      }
      else {
        /*
        minimum and subminimum equal
        */
        if (i0 >= 0) {
          /*
          minimum column j1 is assigned
          swap columns j1 and j2, as j2 may be unassigned
          */
          j1 = j2;
          i0 = li->array[O3_LI_COLSOL][j2];
        }
      }
      /*
      (re-)assign i to j1, possibly de-assigning an i0
      */
      li->array[O3_LI_ROWSOL][i] = j1;
      li->array[O3_LI_COLSOL][j1] = i;

      if (i0 >= 0) {
        /*
        minimum column j1 assigned earlier
        */
        if (umin < usubmin) {
          /*
          put in current k, and go back to that k
          continue augmenting path i - j1 with i0
          */
          li->array[O3_LI_FREE][--k] = i0;
        }
        else {
          /*
          no further augmenting reduction possible
          store i0 in list of free rows for next phase
          */
          li->array[O3_LI_FREE][numfree++] = i0;
        }
      }
    }
  }
  /*
  AUGMENT SOLUTION for each free row
  */
  for (f = 0; f < numfree; ++f) {
    /*
    start row of augmenting path
    */
    freerow = li->array[O3_LI_FREE][f];
    /*
    Dijkstra shortest path algorithm
    runs until unassigned column added to shortest path tree
    */
    for (j = 0; j < dim; ++j) {
      li->array[O3_LI_D][j] = li->cost[freerow][j] - li->array[O3_LI_V][j];
      li->array[O3_LI_PRED][j] = freerow;
      /*
      init column list
      */
      li->array[O3_LI_COLLIST][j] = j;
    }
    /*
    columns in 0..low-1 are ready, now none
    */
    low = 0;
    /*
    columns in low..up-1 are to be scanned for current minimum, now none
    */
    up = 0;
    /*
    columns in up..dim-1 are to be considered later to find new minimum,
    at this stage the list simply contains all columns
    */
    unassignedfound = 0;
    while (!unassignedfound) {
      if (up == low) {
        /*
        no more columns to be scanned for current minimum
        */
        last = low - 1;
        /*
        scan columns for up..dim-1 to find all indices for which new minimum occurs
        store these indices between low..up-1 (increasing up)
        */
        min = li->array[O3_LI_D][li->array[O3_LI_COLLIST][up++]];
        for (k = up; k < dim; ++k) {
          j = li->array[O3_LI_COLLIST][k];
          h = li->array[O3_LI_D][j];
          if (h <= min) {
            if (h < min) {
              /*
              new minimum
              restart list at index low
              */
              up = low;
              min = h;
            }
            /*
            new index with same minimum, put on undex up, and extend list
            */
            li->array[O3_LI_COLLIST][k] = li->array[O3_LI_COLLIST][up];
            li->array[O3_LI_COLLIST][up++] = j;
          }
        }
        /*
        check if any of the minimum columns happens to be unassigned
        if so, we have an augmenting path right away
        */
        for (k = low; k < up; ++k) {
          if (li->array[O3_LI_COLSOL][li->array[O3_LI_COLLIST][k]] < 0) {
            endofpath = li->array[O3_LI_COLLIST][k];
            unassignedfound = 1;
            break;
          }
        }
      }
      if (!unassignedfound) {
        /*
        update 'distances' between freerow and all unscanned columns, via next scanned column
        */
        j1 = li->array[O3_LI_COLLIST][low];
        low++;
        i = li->array[O3_LI_COLSOL][j1];
        h = li->cost[i][j1] - li->array[O3_LI_V][j1] - min;
        for (k = up; k < dim; ++k) {
          j = li->array[O3_LI_COLLIST][k];
          v2 = li->cost[i][j] - li->array[O3_LI_V][j] - h;
          if (v2 < li->array[O3_LI_D][j]) {
            li->array[O3_LI_PRED][j] = i;
            if (v2 == min) {
              /*
              new column found at same minimum value
              */
              if (li->array[O3_LI_COLSOL][j] < 0) {
                /*
                if unassigned, shortest augmenting path is complete
                */
                endofpath = j;
                unassignedfound = 1;
                break;
              }
              else {
                /*
                else add to list to be scanned right away
                */
                li->array[O3_LI_COLLIST][k] = li->array[O3_LI_COLLIST][up];
                li->array[O3_LI_COLLIST][up++] = j;
              }
            }
            li->array[O3_LI_D][j] = v2;
          }
        }
      }
    }
    /*
    update column prices
    */
    for (k = 0; k <= last; ++k) {
      j1 = li->array[O3_LI_COLLIST][k];
      li->array[O3_LI_V][j1] += (li->array[O3_LI_D][j1] - min);
    }
    /*
    reset row and column assignments along the alternating path
    */
    do {
      i = li->array[O3_LI_PRED][endofpath];
      li->array[O3_LI_COLSOL][endofpath] = i;
      j1 = endofpath;
      endofpath = li->array[O3_LI_ROWSOL][i];
      li->array[O3_LI_ROWSOL][i] = j1;
    } while (i != freerow);
  }
}


int compute_cost_matrix(LAPInfo *li, ConfInfo *moved_conf, ConfInfo *template_conf, int n_bins, int coeff, int options)
{
  extern char mmff_corr_matrix[99][99];
  char corr;
  int i;
  int j;
  int k;
  int x;
  int y;
  int largest_n_heavy_atoms = 0;
  double c;
  double h_sum = 0.0;
  
  
  /*
  template_conf is on rows, moved_conf is on columns
  */
  largest_n_heavy_atoms = (template_conf->n_heavy_atoms > moved_conf->n_heavy_atoms)
    ? template_conf->n_heavy_atoms : moved_conf->n_heavy_atoms;
  for (i = 0; i < largest_n_heavy_atoms; ++i) {
    for (j = 0; j < largest_n_heavy_atoms; ++j) {
      li->cost[i][j] = DUMMY_COST;
    }
  }
  c = ((options & MATCH_ATOM_TYPES_BIT) ? 1.0 : 0.0);
  for (i = 0, y = 0; i < template_conf->n_atoms; ++i) {
    if (!strcmp(template_conf->atom[i]->element, "H")) {
      continue;
    }
    for (j = 0, x = 0; j < moved_conf->n_atoms; ++j) {
      if (!strcmp(moved_conf->atom[j]->element, "H")) {
        continue;
      }
      for (k = 0, h_sum = 0.0; (k < n_bins) && (template_conf->h[y][k] + moved_conf->h[x][k]); ++k) {
        h_sum += ((double)square(template_conf->h[y][k] - moved_conf->h[x][k])
          / (double)(template_conf->h[y][k] + moved_conf->h[x][k]));
      }
      if (options & MATCH_CONFORMERS_BIT) {
        if (fabs(template_conf->atom[i]->charge - moved_conf->atom[j]->charge) < ALMOST_ZERO) {
          li->cost[y][x] = (int)safe_rint(h_sum * 1.0e03);
        }
      }
      else {
        corr = mmff_corr_matrix[template_conf->atom[i]->atom_type - 1]
          [moved_conf->atom[j]->atom_type - 1];
        li->cost[y][x] = (int)safe_rint(((double)coeff * CHARGE_WEIGHT * fabs
          (template_conf->atom[i]->charge - moved_conf->atom[j]->charge) +
          c * (double)(5 - coeff) * (double)corr + h_sum) * 1.0e03);
      }
      ++x;
    }
    ++y;
  }
  
  return largest_n_heavy_atoms;
}


int filter_sol_vector(LAPInfo *li, ConfInfo *moved_conf, ConfInfo *template_conf, AtomPair *temp_sdm, AtomPair *sdm)
{
  int i;
  int j;
  int k;
  int n;
  int n_atoms;
  int n_equiv;
  ConfInfo *conf[2];
  
  
  conf[0] = template_conf;
  conf[1] = moved_conf;
  /*
  filter out atom equivalences whose cost is > DUMMY_COST
  */
  for (i = 0, n_equiv = 0; i < template_conf->n_heavy_atoms; ++i) {
    if (li->cost[i][li->array[O3_LI_ROWSOL][i]] >= DUMMY_COST) {
      continue;
    }
    temp_sdm[n_equiv].a[0] = i;
    temp_sdm[n_equiv].a[1] = li->array[O3_LI_ROWSOL][i];
    temp_sdm[n_equiv].cost = li->cost[i][li->array[O3_LI_ROWSOL][i]];
    ++n_equiv;
  }
  for (i = 0, k = 0; i < n_equiv; ++i) {
    if (temp_sdm[i].a[0] == -1) {
      continue;
    }
    memcpy(&sdm[k], &temp_sdm[i], sizeof(AtomPair));
    ++k;
  }
  /*
  find out correspondences between heavy atom indexes
  in the original molecule and heavy atoms in the sdm matrix
  */
  n_equiv = k;
  for (n = 0; n < 2; ++n) {
    for(i = 0; i < n_equiv; ++i) {
      j = 0;
      n_atoms = 0;
      k = sdm[i].a[n];
      while ((n_atoms < conf[n]->n_atoms) && (j <= k)) {
        if (strcmp(conf[n]->atom[n_atoms]->element, "H")) {
          ++j;
          sdm[i].a[n] = n_atoms;
        }
        ++n_atoms;
      }
    }
  }
  for (i = 0; i < n_equiv; ++i) {
    memset(li->diff[i], 0, n_equiv * sizeof(double));
    for (j = 0; j < n_equiv; ++j) {
      /*
      if i == j then the distance will be 0,
      no need to spend time computing it
      */
      if (i == j) {
        continue;
      }
      li->diff[i][j] = fabs(sqrt(squared_euclidean_distance
        (&(template_conf->coord[sdm[i].a[0] * 3]), &(template_conf->coord[sdm[j].a[0] * 3])))
        - sqrt(squared_euclidean_distance
        (&(moved_conf->coord[sdm[i].a[1] * 3]), &(moved_conf->coord[sdm[j].a[1] * 3]))));
    }
  }
  for (i = 0; i < n_equiv; ++i) {
    for (j = 0, sdm[i].score = 0; j < n_equiv; ++j) {
      if (li->diff[i][j] > THRESHOLD_DIFF_DISTANCE) {
        ++(sdm[i].score);
      }
    }
  }
  /*
  sort the SDM matrix by increasing scores
  */
  qsort(sdm, n_equiv, sizeof(AtomPair), compare_score);
  
  return n_equiv;
}
