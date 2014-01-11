/*

parse_o3a_input.c

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
#include <include/error_messages.h>
#include <include/rl_runtime.h>
#ifdef WIN32
#include <include/nice_windows.h>
#endif

#ifdef HAVE_EDITLINE_FUNCTIONALITY
#ifdef HAVE_GNU_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#else
#include <editline/readline.h>
#include <histedit.h>
#endif
#endif


int parse_input(O3Data *od, FILE *input_stream, int run_type)
{
  char *ptr;
  char *context = NULL;
  char prompt[TITLE_LEN];
  char grid_fill[MAX_NAME_LEN];
  char history_file[BUF_LEN];
  char line_orig[BUF_LEN];
  char last_line[BUF_LEN];
  char whole_line[BUF_LEN];
  char file_basename[BUF_LEN];
  char name_list[BUF_LEN];
  char comma_hyphen_list[BUF_LEN];
  char buffer[BUF_LEN];
  char *parameter = NULL;
  char *point;
  char *eof;
  char *conf_dir = NULL;
  unsigned short attr = 0;
  int fail = 0;
  int command = 0;
  int nesting = 0;
  int i;
  int j;
  int space_count;
  int len;
  int found;
  int argn;
  int line_num;
  int overall_line_num = 0;
  int line_complete;
  int spaces;
  int result;
  int import_y_vars = 0;
  int ref_y_var;
  int most_active_percent;
  int criterion;
  int type = 0;
  int state = 0;
  int synonym = 0;
  int list_type = 0;
  int multi_file_type;
  int quote;
  int quote_error;
  int nice_value;
  int options;
  int begin;
  int active_object_num;
  int wrong_object_num;
  int wrong_conf_num;
  int from_file = 0;
  int template_object_num;
  int template_conf_num;
  int template_num;
  int moved_object_num;
  int moved_conf_num;
  long random_seed;
  double outgap = 0.0;
  double trans[3];
  double rot[3];
  CharMat *arg;
  struct timeval start;
  struct timeval end;
  GridInfo temp_grid;
  FileDescriptor log_fd;
  FileDescriptor **source;
  FileDescriptor *current_source;
  O3Data od_comp;
  FILE *touch;
  HIST_ENTRY *last_entry;
  

  memset(comma_hyphen_list, 0, BUF_LEN);
  memset(history_file, 0, BUF_LEN);
  memset(last_line, 0, BUF_LEN);
  memset(&log_fd, 0, sizeof(FileDescriptor));
  memset(&temp_grid, 0, sizeof(GridInfo));
  current_source = od->mel.source;
  while (current_source) {
    ++nesting;
    current_source = current_source->next;
  }
  if (have_editline) {
    if (!(run_type & INTERACTIVE_RUN)) {
      /*
      if this is not an interactive session,
      then we are going to fetch lines by fgets
      so allocate memory for buffer
      */
      od->mel.line = realloc(od->mel.line, BUF_LEN);
      if (!(od->mel.line)) {
        tee_error(od, run_type, overall_line_num,
          E_OUT_OF_MEMORY, PARSE_INPUT_FAILED);
        return PARSE_INPUT_ERROR;
      }
    }
    else {
      /*
      if this is an interactive session,
      then check if the history file exists,
      otherwise create it
      */
      sprintf(history_file, "%s%c%s", od->home_dir, SEPARATOR, HISTORY_FILE);
      if (access(history_file, F_OK) == -1) {
        touch = fopen(history_file, "w+");
        if (touch) {
          fclose(touch);
        }
      }
      /*
      read history from file and get the last entry
      */
      using_history();
      read_history(history_file);
      last_entry = previous_history();
      next_history();
      if (last_entry) {
        strncpy(last_line, last_entry->line, BUF_LEN);
        last_line[BUF_LEN - 1] = '\0';
      }
    }
  }
  else {
    od->mel.line = realloc(od->mel.line, BUF_LEN);
    if (!(od->mel.line)) {
      tee_error(od, run_type, overall_line_num,
        E_OUT_OF_MEMORY, PARSE_INPUT_FAILED);
      return PARSE_INPUT_ERROR;
    }
  }
  /*
  allocate an array of pointers to string buffers
  for command line arguments and values, then store
  the pointers in the MemList substructure
  */
  od->cimal.arg = alloc_char_matrix(od->cimal.arg, MAX_ARG, BUF_LEN);
  if (!(od->cimal.arg)) {
    tee_error(od, run_type, overall_line_num,
      E_OUT_OF_MEMORY, PARSE_INPUT_FAILED);
    return PARSE_INPUT_ERROR;
  }
  arg = od->cimal.arg;
  /*
  main parse loop
  read one line at a time
  */
  while (!fail) {
    /*
    start closing file handles used by previous
    commands (if any); make sure we do not close
    MAIN_INPUT and MAIN_OUTPUT files
    */
    for (i = 0; i < MAX_FILES; ++i) {
      if ((i != MAIN_INPUT) && (i != MAIN_OUTPUT)) {
        if (od->file[i]->handle) {
          fclose(od->file[i]->handle);
          od->file[i]->handle = NULL;
        }
      }
    }
    /*
    if we are running interactively then we are using
    readline which allocates a pointer to the entered
    line; the user is supposed to free it after use
    */
    if (have_editline) {
      if (run_type & INTERACTIVE_RUN) {
        if (od->mel.line) {
          #if (defined HAVE_EDITLINE_FUNCTIONALITY && defined HAVE_GNU_READLINE)
          rl_free(od->mel.line);
          od->mel.line = NULL;
          #endif
          #ifndef HAVE_EDITLINE_FUNCTIONALITY
          if (_dlsym_rl_free) {
            rl_free(od->mel.line);
            od->mel.line = NULL;
          }
          #endif
        }
      }
      else {
        memset(od->mel.line, 0, BUF_LEN);
      }
    }
    else {
      memset(od->mel.line, 0, BUF_LEN);
    }
    memset(line_orig, 0, BUF_LEN);
    line_complete = 0;
    line_num = 0;
    memset(whole_line, 0, BUF_LEN);
    sprintf(prompt, "%s%s", PACKAGE_NAME, SHORT_PROMPT);
    /*
    read line(s) until a line feed which is not
    preceded by a backslash is met
    */
    while (!line_complete) {
      if (run_type & INTERACTIVE_RUN) {
        /*
        if it is an interactive session we use readline
        */
        SET_INK(od, HIGHLIGHT_INK);
        if (have_editline) {
          if (rl_line_buffer && rl_line_buffer[0]) {
            memset(rl_line_buffer, 0, strlen(rl_line_buffer));
          }
          od->mel.line = readline(prompt);
          if (rl_line_buffer) {
            rl_line_buffer[0] = '\0';
          }
          eof = od->mel.line;
        }
        else {
          printf("%s", prompt);
          if ((eof = fgets(od->mel.line, BUF_LEN, input_stream))) {
            od->mel.line[BUF_LEN - 1] = '\0';
            remove_newline(od->mel.line);
          }
        }
        SET_INK(od, NORMAL_INK);
      }
      else {
        /*
        otherwise we use fgets
        */
        if ((eof = fgets(od->mel.line, BUF_LEN, input_stream))) {
          od->mel.line[BUF_LEN - 1] = '\0';
          remove_newline(od->mel.line);
        }
      }
      if (eof) {
        /*
        if line ends with backslash, then wait for
        more input; at the end, put the line together
        */
        len = strlen(od->mel.line);
        /*
        remove spaces eventually present after the last character
        (often present when doing copy/paste from scrolled konsole)
        */
        while (len && isspace(od->mel.line[len - 1])) {
          --len;
        }
        if (len && (od->mel.line[len - 1] == '\\')) {
          ++line_num;
          strcpy(prompt, SHORT_PROMPT);
          od->mel.line[len - 1] = '\0';
          if ((strlen(whole_line) + len) < BUF_LEN) {
            strcat(whole_line, od->mel.line);
          }
          else {
            strncat(whole_line, od->mel.line,
              BUF_LEN - strlen(whole_line));
            line_complete = 1;
          }
        }
        else {
          strncat(whole_line, od->mel.line,
            BUF_LEN - strlen(whole_line));
          line_complete = 1;
        }
      }
      else {
        break;
      }
    }
    if (!eof) {
      /*
      if we are running interactively and the user
      pushes CTRL-D, print a new prompt
      this will prevent editline from quitting
      Open3DALIGN upon window resizing
      */
      if (run_type & INTERACTIVE_RUN) {
        printf("\n");
        continue;
      }
      else {
        break;
      }
    }
    whole_line[BUF_LEN - 1] = '\0';
    i = 0;
    quote = 0;
    /*
    pretreat the line:
    - removing comments
    - converting tabs into spaces
    - converting ' ' into ';' provided it is not
      enclosed between quotes
    - converting '\ ' into ' ' provided it is not
      enclosed between quotes
    */
    len = strlen(whole_line);
    while (i < len) {
      if (isspace(whole_line[i]) && (!quote)) {
        whole_line[i] = ';';
      }
      if (!isprint(whole_line[i])) {
        for (j = i; j < len; ++j) {
          whole_line[j] = whole_line[j + 1];
        }
        --len;
        continue;
      }
      if (whole_line[i] == '#') {
        whole_line[i] = '\0';
        break;
      }
      if (whole_line[i] == '\"') {
        quote = 1 - quote;
      }
      if ((whole_line[i] == '\\') && (!quote)) {
        if ((i + 1) < len) {
          if (whole_line[i + 1] == ' ') {
            for (j = i; j < len; ++j) {
              whole_line[j] = whole_line[j + 1];
            }
            --len;
            ++i;
            continue;
          }
        }
      }
      ++i;
    }
    /*
    if quotes were opened but never closed issue an error
    */
    quote_error = 0;
    if (quote) {
      tee_error(od, run_type, overall_line_num, "Error: missing quote\n\n");
      quote_error = 1;
    }
    argn = 0;
    
    /*
    copy command line arguments in the "arg"
    array until MAX_ARG is reached
    */
    ptr = strtok_r(whole_line, ";", &context);
    while (ptr && (argn < MAX_ARG)) {
      strcpy(arg->me[argn], ptr);
      if (argn) {
        strcat(line_orig, ";");
      }
      strcat(line_orig, ptr);
      ++argn;
      ptr = strtok_r(NULL, ";", &context);
    }
    len = strlen(line_orig);
    i = 0;
    quote = 0;
    /*
    regenerate the original line:
    - converting ';' into ' ' provided it is not
      enclosed between quotes
    - converting ' ' into '\ ' provided it is not
      enclosed between quotes
    so basically we remove from the entered line
    comments and multiple spaces
    */
    while (i < len) {
      if (whole_line[i] == '\"') {
        quote = 1 - quote;
      }
      if ((line_orig[i] == ' ') && (!quote)) {
        for (j = len; j >= i; --j) {
          line_orig[j + 1] = line_orig[j];
        }
        line_orig[i] = '\\';
        ++i;
        ++len;
      }
      else if ((line_orig[i] == ';') && (!quote)) {
        line_orig[i] = ' ';
      }
      ++i;
    }
    if (line_orig[0] && (run_type & INTERACTIVE_RUN)
      && strcmp(last_line, line_orig)) {
      if (have_editline) {
        add_history(line_orig);
        write_history(history_file);
      }
      strcpy(last_line, line_orig);
    }
    /*
    if a quote was missing it is time to bomb out,
    but at least the user finds the line in history
    so he can correct it easily
    */
    if (quote_error) {
      fail = !(run_type & INTERACTIVE_RUN);
      continue;
    }
    /*
    if there are no arguments (because the command line
    is empty or constituted by comments only), we need
    to count it or pretreatment errors will not refer to
    the correct text line
    */
    if (!argn) {
      ++overall_line_num;
      continue;
    }
    od->argn = argn;
    /*
    the first argument is the command itself
    */
    if ((parameter = get_args(od, NULL))) {
      tee_error(od, run_type, overall_line_num, "Error parsing "
        "the following command:\n\n"
        "> %s\n\n"
        "Please check your input near:\n\n", line_orig);
      point = strstr(line_orig, parameter);
      if (point) {
        spaces = 0;
        space_count = 0;
        while (((point - spaces) > line_orig)
          && (spaces <= MAX_SPACES)) {
          if (*(point - spaces) == ' ') {
            ++space_count;
          }
          if (space_count == 2) {
            break;
          }
          ++spaces;
        }
        i = 0;
        space_count = 0;
        while (point[i] && (point[i] != '\n')
          && (point[i] != '\r')
          && (i <= MAX_SPACES)) {
          if (point[i] == ' ') {
            ++space_count;
          }
          if (space_count == 2) {
            break;
          } 
          ++i;
        }
        point[i] = '\0';
        if (spaces > MAX_SPACES) {
          tee_printf(od, "[...] ");
        }
        tee_printf(od, "%s", (char *)(point - spaces));
        if (i > MAX_SPACES) {
          tee_printf(od, "[...] ");
        }
        tee_printf(od, "\n");
        for (i = 0; i < spaces; ++i) {
          tee_printf(od, " ");
        }
        tee_printf(od, "^\n");
      }
      fail = !(run_type & INTERACTIVE_RUN);
      continue;
    }
    overall_line_num += (line_num + 1);

    if (!(od->file[TEMP_OUT]->name[0])) {
      result = open_temp_file(od, od->file[TEMP_OUT], "temp_out");
      if (result) {
        tee_error(od, run_type, overall_line_num,
          E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
          od->file[TEMP_OUT]->name, IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (od->file[TEMP_OUT]->handle) {
        fclose(od->file[TEMP_OUT]->handle);
        od->file[TEMP_OUT]->handle = NULL;
      }
    }
    if (!(od->file[TEMP_LOG]->name[0])) {
      result = open_temp_file(od, od->file[TEMP_LOG], "temp_log");
      if (result) {
        tee_error(od, run_type, overall_line_num,
          E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
          od->file[TEMP_LOG]->name, IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (od->file[TEMP_LOG]->handle) {
        fclose(od->file[TEMP_LOG]->handle);
        od->file[TEMP_LOG]->handle = NULL;
      }
    }
    if (!strcasecmp(arg->me[0], "box")) {
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, BOX_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if ((parameter = get_args(od, "mode"))) {
        if (!strncasecmp(parameter, "get", 3)) {
          print_grid_coordinates(od, &(od->grid));
          continue;
        }
      }
      if (od->field_num) {
        tee_error(od, run_type, overall_line_num,
          E_CANNOT_ALTER_GRID_BOX,
          BOX_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      memset(grid_fill, 0, MAX_NAME_LEN);
      temp_grid.step[0] = 1.0;
      outgap = 5.0;
      from_file = 0;
      if ((parameter = get_args(od, "file"))) {
        grid_fill[0] = 1;
        for (i = 0; i < 6; ++i) {
          grid_fill[i] = 1;
        }
        outgap = -1.0;
        strcpy(od->file[ASCII_IN]->name, parameter);
        if (!(od->file[ASCII_IN]->handle =
          fopen(od->file[ASCII_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            od->file[ASCII_IN]->name, BOX_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        from_file = 1;
      }
      else {
        for (result = 0, i = 0; i < 3; ++i) {
          memset(buffer, 0, BUF_LEN);
          sprintf(buffer, "%c_start", i + 'x');
          if ((parameter = get_args(od, buffer))) {
            grid_fill[i] = 1;
            sscanf(parameter, "%f", &(temp_grid.start_coord[i]));
          }
          synonym = 0;
          sprintf(buffer, "%c_end", i + 'x');
          if ((parameter = get_args(od, buffer))) {
            grid_fill[i + 3] = 1;
            sscanf(parameter, "%f", &(temp_grid.end_coord[i]));
            ++synonym;
          }
          sprintf(buffer, "%c_nodes", i + 'x');
          if ((parameter = get_args(od, buffer))) {
            grid_fill[i + 3] = 1;
            sscanf(parameter, "%d", &(temp_grid.nodes[i]));
            ++synonym;
          }
          if (synonym > 1) {
            tee_error(od, run_type, overall_line_num,
              "Conflicting grid definitions were detected; "
              "please specify either %c_end or %c_nodes.\n",
              i + 'x', i + 'x');
            result = 1;
          }
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            "%s", BOX_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if ((parameter = get_args(od, "step"))) {
          grid_fill[6] = 1;
          sscanf(parameter, "%f", &(temp_grid.step[0]));
        }
        if ((parameter = get_args(od, "outgap"))) {
          grid_fill[7] = 1;
          sscanf(parameter, "%lf", &outgap);
        }
        if (grid_fill[7]) {
          /*
          if the outgap parameter is given, no other parameters
          should be given except step (which defaults to 1.0)
          */
          for (i = 0, j = 0; ((i <= 5) && !j); ++i) {
            j += (int)grid_fill[i];
          }
          if (j) {
            tee_error(od, run_type, overall_line_num,
              "When the outgap parameter is input, "
              "only the step parameter may be "
              "present.\n%s",
              BOX_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          outgap = (double)safe_rint(fabs(outgap) * 1.0e04) / 1.0e04;
        }
        else {
          for (i = 0, j = 0; i < 6; ++i) {
            j += (int)grid_fill[i];
          }
          /*
          if the outgap parameter is not given, full grid
          start and end coordinates should be supplied,
          or none; in the latter case the default outgap (5.0)
          will be used
          */
          if (j && (j < 6)) {
            tee_error(od, run_type, overall_line_num,
              "Full grid box coordinates are needed.\n%s",
              BOX_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          if (j == 6) {
            outgap = -1.0;
          }
        }
        if (temp_grid.step[0] < ALMOST_ZERO) {
          tee_error(od, run_type, overall_line_num,
            "The step value must "
            "be greater than 0.0.\n%s",
            BOX_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "BOX", line_orig);
        tee_flush(od);
        result = create_box(od, &temp_grid, outgap, from_file);
        switch (result) {
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "CoMFA region",
            od->file[ASCII_IN]->name, BOX_FAILED);
          return PARSE_INPUT_ERROR;
        }
        print_grid_coordinates(od, &(od->grid));
        if (od->file[ASCII_IN]->handle) {
          fclose(od->file[ASCII_IN]->handle);
          od->file[ASCII_IN]->handle = NULL;
        }
        update_pymol(od);
        update_jmol(od);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "BOX");
        tee_flush(od);
      }
      else {
        /*
        fake data for DRY_RUN
        */
        od->grid.nodes[0] = 10;
        od->grid.nodes[1] = 10;
        od->grid.nodes[2] = 10;
        od->grid.object_num = 10;
        od->grid.struct_num = 10;
        od->object_num = 10;
      }
    }
    else if (!strcasecmp(arg->me[0], "rototrans")) {
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (od->field_num) {
        tee_error(od, run_type, overall_line_num,
          E_CANNOT_ALTER_OBJECT_COORDINATES,
          ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num,
          E_FILE_FOR_MODIFIED_COORDINATES,
          "transformed", ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      strcpy(file_basename, parameter);
      gettimeofday(&start, NULL);
      result = parse_synonym_lists(od, "ROTOTRANS", ROTOTRANS_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST) | (1 << STRUCT_LIST), &list_type,
        OBJECT_LIST, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      memset(trans, 0, 3 * sizeof(double));
      memset(rot, 0, 3 * sizeof(double));
      for (i = 0; i < 3; ++i) {
        sprintf(buffer, "%c_trans", i + 'x');
        if ((parameter = get_args(od, buffer))) {
          sscanf(parameter, "%lf", &trans[i]);
        }
        sprintf(buffer, "%c_rot", i + 'x');
        if ((parameter = get_args(od, buffer))) {
          sscanf(parameter, "%lf", &rot[i]);
        }
      }
      for (i = 0; i < 3; ++i) {
        rot[i] = angle2rad(rot[i]);
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "ROTOTRANS", line_orig);
        tee_flush(od);
        result = call_obenergy(od, O3_MMFF94);
        switch (result) {
          case FL_CANNOT_CREATE_CHANNELS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PIPE, ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CHDIR:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_PROCESS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PROCESS, "OpenBabel", ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case OPENBABEL_ERROR:
          tee_error(od, run_type, overall_line_num,
            E_PROGRAM_ERROR, "OpenBabel");
          if ((od->file[TEMP_LOG]->handle = fopen
            (od->file[TEMP_LOG]->name, "rb"))) {
            while (fgets(buffer, BUF_LEN,
              od->file[TEMP_LOG]->handle)) {
              buffer[BUF_LEN - 1] = '\0';
              tee_printf(od, "%s", buffer);
            }
            fclose(od->file[TEMP_LOG]->handle);
            od->file[TEMP_LOG]->handle = NULL;
            tee_printf(od, "\n%s", ROTOTRANS_FAILED);
          }
          else {
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", ROTOTRANS_FAILED);
          }
          return PARSE_INPUT_ERROR;
        }
        result = rototrans(od, file_basename, trans, rot);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_ROTOTRANSED_SDF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_SDF_FILE,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_ORIGINAL_SDF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_SDF_FILE,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_READ_MOL_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_MOL_FILE,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_READ_OB_OUTPUT:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_OB_OUTPUT,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_UNKNOWN_ATOM_TYPE:
          tee_error(od, run_type, overall_line_num,
            E_UNKNOWN_ATOM_TYPE, "atom", ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "ROTOTRANS");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "import")) {
      /*
      user wants to IMPORT something; check
      if the file type was specified
      */
      multi_file_type = 0;
      gettimeofday(&start, NULL);
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,
          E_SPECIFY_FILE_TYPE, IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!strcasecmp(parameter, "sdf")) {
        if (!(od->field.babel_exe_path[0])) {
          tee_error(od, run_type, overall_line_num,
            E_OPENBABEL_PATH, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            E_SPECIFY_STRUCTURE_FILE, "import", IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        memset(od->file[MOLFILE_IN]->name, 0, BUF_LEN);
        strcpy(od->file[MOLFILE_IN]->name, parameter);
        absolute_path(od->file[MOLFILE_IN]->name);
        if (!(run_type & DRY_RUN)) {
          if (!fexist(od->file[MOLFILE_IN]->name)) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_READING,
              parameter, IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          if (od->valid & SDF_BIT) {
            od->valid = 0;
            od->object_num = 0;
            od->grid.object_num = 0;
            od->active_object_num = 0;
            od->ext_pred_object_num = 0;
            free_y_var_array(od);
            memset(&(od->grid), 0, sizeof(GridInfo));
            if (od->pel.pymol_old_object_id) {
              int_perm_free(od->pel.pymol_old_object_id);
              od->pel.pymol_old_object_id = NULL;
            }
            if (od->pel.jmol_old_object_id) {
              int_perm_free(od->pel.jmol_old_object_id);
              od->pel.jmol_old_object_id = NULL;
            }
          }
          import_y_vars = 0;
          memset(name_list, 0, BUF_LEN);
          if ((parameter = get_args(od, "y_var_name"))) {
            strcpy(name_list, parameter);
            import_y_vars |= IMPORT_Y_VARS_BIT;
          }
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command,
            "IMPORT SDF", line_orig);
          tee_flush(od);
          if (!(od->file[MOLFILE_IN]->handle =
            fopen(od->file[MOLFILE_IN]->name, "rb"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_READING,
              parameter, IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = open_temp_dir(od, NULL, "mol_dir", od->field.mol_dir);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_DIR_CANNOT_BE_CREATED, od->field.mol_dir,
              IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          memset(&(od->grid), 0, sizeof(GridInfo));
          result = parse_sdf(od, import_y_vars | ALLOC_MOL_INFO_BIT, name_list);
          if (od->file[MOLFILE_IN]->handle) {
            fclose(od->file[MOLFILE_IN]->handle);
            od->file[MOLFILE_IN]->handle = NULL;
          }
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          switch (result) {
            case PREMATURE_EOF:
            tee_error(od, run_type, overall_line_num, 
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "SDF",
              od->file[MOLFILE_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_WRITING_TEMP_FILE,
              od->task.string, IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_FIND_Y_VAR_NAME:
            tee_error(od, run_type, overall_line_num, E_CANNOT_FIND_Y_VAR_NAME,
              od->file[MOLFILE_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case Y_VAR_LOW_SD:
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_LOW_SD, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;
          }
          if (alloc_object_attr(od, 0)) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;
          }
          od->valid |= SDF_BIT;
          if (update_mol(od)) {
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE,
              od->task.string, IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
          }
          update_field_object_attr(od, VERBOSE_BIT);
          update_pymol(od);
          update_jmol(od);
          tee_printf(od, "\n"
            "This set of molecules fits in a grid box "
            "whose bottom left, top right cooordinates "
            "(no outgap) are at least:\n"
            "[(%.4f,%.4f,%.4f), (%.4f,%.4f,%.4f)]\n\n",
            od->min_coord[0], od->min_coord[1],
            od->min_coord[2], od->max_coord[0],
            od->max_coord[1], od->max_coord[2]);
          strcpy(od->default_folder, od->file[MOLFILE_IN]->name);
          get_dirname(od->default_folder);
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "IMPORT SDF");
          tee_flush(od);
        }
        else {
          od->valid |= SDF_BIT;
          od->grid.object_num = 10;
        }
      }
      else if (!strcasecmp(parameter, "dependent")) {
        /*
        user wants to IMPORT dependent variables
        */
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            "Please specify a file from which you wish "
            "to import dependent variables.\n%s",
            IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(od->valid & SDF_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_IMPORT_MOLFILE_FIRST, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        we need a txt file
        */
        memset(od->file[DEP_IN]->name, 0, BUF_LEN);
        strcpy(od->file[DEP_IN]->name, parameter);
        absolute_path(od->file[DEP_IN]->name);
        if (!((od->file[DEP_IN]->handle) =
          fopen(od->file[DEP_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            parameter, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        import the txt file
        */
        if (!(run_type & DRY_RUN)) {
          strcpy(name_list, "all");
          if ((parameter = get_args(od, "y_var_name"))) {
            strcpy(name_list, parameter);
          }
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command,
            "IMPORT DEPENDENT", line_orig);
          tee_flush(od);
          result = import_dependent(od, name_list);
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          /*
          check for errors
          */
          switch (result) {
            case CANNOT_READ_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case NOT_ENOUGH_OBJECTS:
            tee_error(od, run_type, overall_line_num,
              "In the file \"%s\" there must be exactly "
              "%d lines, the first one with variable "
              "names followed by dependent variable values.\n%s",
              od->file[DEP_IN]->name, od->active_object_num + 1,
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case WRONG_NUMBER_OF_Y_VARS:
            tee_error(od, run_type, overall_line_num, E_WRONG_NUMBER_OF_Y_VARS,
              od->file[DEP_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case PREMATURE_DEP_IN_EOF:
            tee_error(od, run_type, overall_line_num,
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, ".txt",
              od->file[DEP_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_FIND_Y_VAR_NAME:
            tee_error(od, run_type, overall_line_num, E_CANNOT_FIND_Y_VAR_NAME,
              od->file[DEP_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case Y_VAR_LOW_SD:
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_LOW_SD, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            default:
            tee_printf(od, M_TOOL_SUCCESS, nesting, command, "IMPORT DEPENDENT");
            tee_flush(od);
          }
        }
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Only \"SDF\", and \"DEPENDENT\" "
          "types are allowed "
          "for the IMPORT keyword.\n%s",
          IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
    }
    else if (!strcasecmp(arg->me[0], "set")) {
      gettimeofday(&start, NULL);
      result = parse_synonym_lists(od, "SET", SET_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST) | (1 << STRUCT_LIST),
        &list_type, 0, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      if ((parameter = get_args(od, "attribute"))) {
        if (!strncasecmp(parameter, "training", 8)) {
          attr = ACTIVE_BIT;
          state = 1;
        }
        else if (!strncasecmp(parameter, "excluded", 8)) {
          attr = ACTIVE_BIT;
          state = 0;
        }
        else if (!strncasecmp(parameter, "test", 4)) {
          attr = PREDICT_BIT;
          state = 1;
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "Only \"TRAINING\", \"EXCLUDED\", "
            "\"TEST\" parameters are allowed.\n%s",
            SET_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Please specify an attribute which should be set.\n%s",
          SET_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SET", line_orig);
        tee_flush(od);
        result = set(od, intlog2(list_type), attr, state, VERBOSE_BIT);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", SET_FAILED);
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, SET_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            "The range you selected is not valid.\n%s", SET_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;

          default:
          update_pymol(od);
          update_jmol(od);
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SET");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "env")) {
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "ENV", line_orig);
        tee_flush(od);
      }
      if ((parameter = get_args(od, "n_cpus"))) {
        determine_best_cpu_number(od, "all");
        determine_best_cpu_number(od, parameter);
        tee_printf(od, M_NUMBER_OF_CPUS, PACKAGE_NAME, od->n_proc);
      }
      else if ((parameter = get_args(od, "random_seed"))) {
        if (!(run_type & DRY_RUN)) {
          sscanf(parameter, "%ld", &random_seed);
          set_random_seed(od, (unsigned long)absval(random_seed));
          tee_printf(od, "The random seed has been set to %ld.\n\n",
            random_seed);
        }
      }
      else if ((parameter = get_args(od, "temp_dir"))) {
        if (dexist(parameter)) {
          strcpy(od->temp_dir, parameter);
          absolute_path(od->temp_dir);
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "The temporary directory has been set to %s.\n\n",
              od->temp_dir);
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            E_DIR_NOT_EXISTING, parameter, ENV_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else if ((parameter = get_args(od, "babel_path"))) {
        memset(od->field.babel_exe_path, 0, BUF_LEN);
        found = 0;
        result = 1;
        if (parameter[0]) {
          strncpy(od->field.babel_exe_path, parameter, BUF_LEN - 1);
          absolute_path(od->field.babel_exe_path);
          if (dexist(od->field.babel_exe_path)) {
            sprintf(buffer, "%s%c%s", od->field.babel_exe_path, SEPARATOR, BABEL_EXE);
            if (fexist(buffer)) {
              sprintf(buffer, "%s%c%s", od->field.babel_exe_path, SEPARATOR, OBENERGY_EXE);
              if (fexist(buffer)) {
                found = 1;
              }
            }
          }
        }
        if (!(run_type & DRY_RUN)) {
          if (found) {
            result = check_babel(od, od->field.babel_exe_path);
            switch (result) {
              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, ENV_FAILED);
              return PARSE_INPUT_ERROR;

              case CANNOT_READ_TEMP_FILE:
              tee_printf(od, E_FILE_CANNOT_BE_OPENED_FOR_READING,
                od->file[TEMP_LOG]->name);
              break;

              case CANNOT_WRITE_TEMP_FILE:
              tee_printf(od, E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
                od->file[TEMP_LOG]->name);
              break;

              case BABEL_NOT_WORKING:
              tee_printf(od, E_OPENBABEL_NOT_WORKING);
              if ((od->file[TEMP_LOG]->handle = fopen(od->file[TEMP_LOG]->name, "rb"))) {
                while (fgets(buffer, BUF_LEN, od->file[TEMP_LOG]->handle)) {
                  buffer[BUF_LEN - 1] = '\0';
                  tee_printf(od, "%s", buffer);
                }
                fclose(od->file[TEMP_LOG]->handle);
                od->file[TEMP_LOG]->handle = NULL;
              }
              break;

              case BABEL_PLUGINS_NOT_FOUND:
              tee_printf(od, E_OPENBABEL_DATA_PLUGINS, "");
              break;
            }
          }
          if ((!found) || (result)) {
            if (found) {
              tee_printf(od, "Working ");
            }
            tee_error(od, run_type, overall_line_num,
              E_OPENBABEL_MISSING_OR_TOO_OLD, ENV_FAILED);
            memset(od->field.babel_exe_path, 0, BUF_LEN);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          tee_printf(od, M_EXECUTABLE_PATH, "OpenBabel",
            od->field.babel_exe_path);
        }
      }
      else if ((parameter = get_args(od, "tinker_path"))) {
        memset(od->qmd.tinker_exe_path, 0, BUF_LEN);
        found = 0;
        if (parameter[0]) {
          strncpy(od->qmd.tinker_exe_path, parameter, BUF_LEN - 1);
          absolute_path(od->qmd.tinker_exe_path);
          if (dexist(od->qmd.tinker_exe_path)) {
            sprintf(buffer, "%s%c%s", od->qmd.tinker_exe_path, SEPARATOR, TINKER_DYNAMIC_EXE);
            if (fexist(buffer)) {
              sprintf(buffer, "%s%c%s", od->qmd.tinker_exe_path, SEPARATOR, TINKER_OPTIMIZE_EXE);
              if (fexist(buffer)) {
                sprintf(buffer, "%s%c%s", od->qmd.tinker_exe_path, SEPARATOR, TINKER_MINIMIZE_EXE);
                if (fexist(buffer)) {
                  found = 1;
                }
              }
            }
          }
        }
        if (found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, M_EXECUTABLE_PATH, "TINKER",
              od->qmd.tinker_exe_path);
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "TINKER binaries could not be found.\n%s\n", ENV_FAILED);
          memset(od->qmd.tinker_exe_path, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else if ((parameter = get_args(od, "pharao"))) {
        memset(od->align.pharao_exe, 0, BUF_LEN);
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->align.pharao_exe);
          }
          else {
            strcpy(od->align.pharao_exe, parameter);
          }
        }
        if (!found) {
          tee_error(od, run_type, overall_line_num,
            "The PHARAO binary %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->align.pharao_exe, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        else {
          absolute_path(od->align.pharao_exe);
          strcpy(od->align.pharao_exe_path, od->align.pharao_exe);
          get_dirname(od->align.pharao_exe_path);
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, M_EXECUTABLE, "PHARAO",
              od->align.pharao_exe);
          }
        }
      }
      else if ((parameter = get_args(od, "pymol"))) {
        memset(od->pymol.pymol_exe, 0, BUF_LEN);
        od->pymol.use_pymol = 0;
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->pymol.pymol_exe);
          }
          else {
            strcpy(od->pymol.pymol_exe, parameter);
          }
        }
        if (parameter[0] && (!found)) {
          tee_error(od, run_type, overall_line_num,
            "The PyMOL binary %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->pymol.pymol_exe, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "PyMOL will not be used.\n\n");
          }
        }
        else {
          od->pymol.use_pymol = 1;
          absolute_path(od->pymol.pymol_exe);
          tee_printf(od, M_EXECUTABLE, "PyMOL",
            od->pymol.pymol_exe);
          strcat(od->pymol.pymol_exe, " " PYMOL_ARGS);
        }
      }
      else if ((parameter = get_args(od, "jmol"))) {
        memset(od->jmol.jmol_exe, 0, BUF_LEN);
        od->jmol.use_jmol = 0;
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->jmol.jmol_exe);
          }
          else {
            strcpy(od->jmol.jmol_exe, parameter);
          }
        }
        if (parameter[0] && (!found)) {
          tee_error(od, run_type, overall_line_num,
            "The Jmol start script %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->jmol.jmol_exe, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "Jmol will not be used.\n\n");
          }
        }
        else {
          od->jmol.use_jmol = 1;
          absolute_path(od->jmol.jmol_exe);
          tee_printf(od, M_EXECUTABLE, "Jmol",
            od->jmol.jmol_exe);
          strcat(od->jmol.jmol_exe, " " JMOL_ARGS);
        }
      }
      else if ((parameter = get_args(od, "nice"))) {
        if (!(run_type & DRY_RUN)) {
          #ifndef WIN32
          sscanf(parameter, "%d", &nice_value);
          #else
          for (nice_value = 0; nice_value < 6; ++nice_value) {
            if (!strcasecmp(parameter, nice_name[nice_value])) {
              break;
            }
          }
          #endif
          set_nice_value(od, nice_value);
        }
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Allowed environmental variables which may be set are: "
          "\"random_seed\", \"temp_dir\", \"n_cpus\", \"nice\", "
          "\"babel_path\", \"tinker_path\", \"pharao\", "
          "\"jmol\" and \"pymol\".\n%s",
          ENV_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "ENV");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "load")) {
      gettimeofday(&start, NULL);
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, LOAD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      /*
      user wants to LOAD a .dat file; check
      if the file name was specified
      */
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify an "PACKAGE_NAME" file "
          "from which data should be loaded.\n%s",
          LOAD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      strcpy(od->file[DAT_IN]->name, parameter);
      absolute_path(od->file[DAT_IN]->name);
      options = VERBOSE_BIT;
      if ((parameter = get_args(od, "mode"))) {
        if (od->grid.object_num
          && (!strncasecmp(parameter, "append", 6))) {
          options |= APPEND_BIT;
        }
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "LOAD", line_orig);
        tee_flush(od);
        if (!(options & APPEND_BIT)) {
          close_files(od, MAX_FILES);
          remove_temp_files(od->package_code);
          if (od->pel.pymol_old_object_id) {
            int_perm_free(od->pel.pymol_old_object_id);
            od->pel.pymol_old_object_id = NULL;
          }
          if (od->pel.jmol_old_object_id) {
            int_perm_free(od->pel.jmol_old_object_id);
            od->pel.jmol_old_object_id = NULL;
          }
        }
        if (!(od->file[DAT_IN]->handle = (FILE *)
          fzopen(od->file[DAT_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            od->file[DAT_IN]->name, LOAD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        load .dat file
        */
        if (!(options & APPEND_BIT)) {
          result = open_temp_dir(od, NULL, "mol_dir", od->field.mol_dir);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_DIR_CANNOT_BE_CREATED, od->field.mol_dir,
              LOAD_FAILED);
              if (od->file[DAT_IN]->handle) {
                fzclose((fzPtr *)(od->file[DAT_IN]->handle));
                od->file[DAT_IN]->handle = NULL;
              }
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
          }
        }
        sprintf(od->file[TEMP_MOLFILE]->name, "%s%cloaded_molfile.sdf",
          od->field.mol_dir, SEPARATOR);
        result = load_dat(od, DAT_IN, VERBOSE_BIT);
        if (od->file[DAT_IN]->handle) {
          fzclose((fzPtr *)(od->file[DAT_IN]->handle));
          od->file[DAT_IN]->handle = NULL;
        }
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, ".dat",
            od->file[DAT_IN]->name, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case WRONG_NUMBER_OF_FIELDS:
          tee_error(od, run_type, overall_line_num,
            "In the .dat file \"%s\" the number of fields "
            "does not match the number of currently loaded "
            "fields.\n%s",
            od->file[DAT_IN]->name, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case WRONG_NUMBER_OF_X_VARS:
          tee_error(od, run_type, overall_line_num,
            "In the .dat file \"%s\" the number of X variables "
            "does not match the number of currently loaded "
            "X variables.\n%s",
            od->file[DAT_IN]->name, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
            
          case WRONG_NUMBER_OF_Y_VARS:
          tee_error(od, run_type, overall_line_num,
            "In the .dat file \"%s\" the number of Y variables "
            "does not match the number of currently loaded "
            "Y variables.\n%s",
            od->file[DAT_IN]->name, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
            
          case GRID_NOT_MATCHING:
          tee_error(od, run_type, overall_line_num,
            E_GRID_NOT_MATCHING, "DAT",
            od->file[DAT_IN]->name, "");
          print_grid_comparison(od);
          tee_error(od, run_type, overall_line_num, "%s\n",
            LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case STATUS_INCONSISTENCY:
          tee_error(od, run_type, overall_line_num,
            "In the .dat file \"%s\" the active/inactive "
            "status of objects is not consistent across the "
            "different fields.\n%s",
            od->file[DAT_IN]->name, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_CHANNELS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PIPE, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CHDIR:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_PROCESS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PROCESS, "OpenBabel", LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case OPENBABEL_ERROR:
          tee_error(od, run_type, overall_line_num,
            E_PROGRAM_ERROR, "OpenBabel");
          if ((od->file[TEMP_LOG]->handle = fopen
            (od->file[TEMP_LOG]->name, "rb"))) {
            while (fgets(buffer, BUF_LEN,
              od->file[TEMP_LOG]->handle)) {
              buffer[BUF_LEN - 1] = '\0';
              tee_printf(od, "%s", buffer);
            }
            tee_printf(od, "\n%s", LOAD_FAILED);
            fclose(od->file[TEMP_LOG]->handle);
            od->file[TEMP_LOG]->handle = NULL;
          }
          else {
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", LOAD_FAILED);
          }
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
        }
        if (update_mol(od)) {
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE,
            od->task.string, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
        }
        update_pymol(od);
        update_jmol(od);
        /*
        the following is for align, so that the default align_dir
        is defined also if user LOADs a DAT file
        instead of IMPORTing a SDF file
        */
        strcpy(od->default_folder, od->file[DAT_IN]->name);
        get_dirname(od->default_folder);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "LOAD");
        tee_flush(od);
      }
      else {
        /*
        fake data for DRY_RUN
        */
        od->grid.nodes[0] = 10;
        od->grid.nodes[1] = 10;
        od->grid.nodes[2] = 10;
        od->grid.object_num = 10;
        od->object_num = 10;
        od->valid |= SDF_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "save")) {
      gettimeofday(&start, NULL);
      /*
      user wants to SAVE a .dat file; check
      if the file name was specified
      */
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num, "Please specify a file "
          "where data should be saved.\n%s",
          SAVE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      memset(od->file[DAT_OUT]->name, 0, BUF_LEN);
      strcpy(od->file[DAT_OUT]->name, parameter);
      absolute_path(od->file[DAT_OUT]->name);
      if (!(run_type & DRY_RUN)) {
        if (!(od->grid.object_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_OBJECTS_PRESENT, SAVE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        save .dat file
        */
        if (!(od->file[DAT_OUT]->handle =
          fopen(od->file[DAT_OUT]->name, "wb+"))) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[DAT_OUT]->name, SAVE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SAVE", line_orig);
        tee_flush(od);
        result = save_dat(od, DAT_OUT);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE,
            od->file[DAT_OUT]->name, SAVE_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SAVE_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SAVE");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "remove_box")) {
      if (!(od->grid.nodes[0])) {
        tee_error(od, run_type, overall_line_num,
          E_NO_GRID_BOX_PRESENT,
          REMOVE_BOX_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "REMOVE_BOX", line_orig);
        tee_flush(od);
        remove_box(od);
        update_pymol(od);
        update_jmol(od);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_BOX");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "remove_object")) {
      gettimeofday(&start, NULL);
      result = parse_synonym_lists(od, "REMOVE_OBJECT", REMOVE_OBJECT_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST) | (1 << STRUCT_LIST), &list_type,
        0, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "REMOVE_OBJECT", line_orig);
        tee_flush(od);
        result = remove_object(od);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, REMOVE_OBJECT_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", REMOVE_OBJECT_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", REMOVE_OBJECT_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          if (!(od->grid.object_num)) {
            tee_printf(od, "Since all objects were removed, all grid "
              "information has been removed as well.\n");
          }
          update_pymol(od);
          update_jmol(od);
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_OBJECT");
          tee_flush(od);
        }
      }
      else {
        od->valid &= SDF_BIT;
      }
    }  
    else if (!strcasecmp(arg->me[0], "remove_y_vars")) {
      gettimeofday(&start, NULL);
      if (od->y_vars) {
        if (!(run_type & DRY_RUN)) {
          strcpy(comma_hyphen_list, "all");
          if ((parameter = get_args(od, "y_var_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          result = parse_comma_hyphen_list_to_array
            (od, comma_hyphen_list, Y_VAR_LIST);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              E_LIST_PARSING, "Y variables", "EXPORT", REMOVE_Y_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = set(od, Y_VAR_LIST, OPERATE_BIT, 1, SILENT);
          if (result) {
            tee_error(od, run_type, overall_line_num, E_Y_VAR_ALLOWED_RANGE,
              od->y_vars, REMOVE_Y_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = remove_y_vars(od);
          switch (result) {
            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_WRITING_TEMP_FILE, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_READ_TEMP_FILE:
            case NOT_ENOUGH_OBJECTS:
            case WRONG_NUMBER_OF_Y_VARS:
            case PREMATURE_DEP_IN_EOF:
            case CANNOT_FIND_Y_VAR_NAME:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;

            case Y_VAR_LOW_SD:
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_LOW_SD, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;
          }
        }
      }
      else {
        tee_printf(od, "No Y variables are present; nothing was done.\n\n");
      }
      tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_Y_VARS");
      tee_flush(od);
    }
    else if (!strcasecmp(arg->me[0], "qmd")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, QMD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, QMD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->qmd.tinker_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_TINKER_PATH, QMD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      sprintf(od->qmd.tinker_prm_path, "%s%c..%cshare%ctinker",
        od->qmd.tinker_exe_path, SEPARATOR, SEPARATOR, SEPARATOR);
      found = 0;
      if (dexist(od->qmd.tinker_prm_path)) {
        sprintf(buffer, "%s%c%s", od->qmd.tinker_prm_path,
          SEPARATOR, TINKER_MMFF94_PRM_FILE);
        if (fexist(buffer)) {
          sprintf(buffer, "%s%c%s", od->qmd.tinker_prm_path,
            SEPARATOR, TINKER_MMFF94S_PRM_FILE);
          found = fexist(buffer);
        }
      }
      if (!found) {
        memset(od->qmd.tinker_prm_path, 0, BUF_LEN);
      }
      if (!(parameter = get_args(od, "prm_dir"))) {
        if (!(od->qmd.tinker_prm_path[0])) {
          tee_error(od, run_type, overall_line_num,
            E_MISSING_DIR, "TINKER MMFF parameter",
            "prm_dir", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else {
        strncpy(od->qmd.tinker_prm_path, parameter, BUF_LEN - 1);
        od->qmd.tinker_prm_path[BUF_LEN - 1] = '\0';
        absolute_path(od->qmd.tinker_prm_path);
      }
      strcpy(od->qmd.minimizer, TINKER_OPTIMIZE_EXE);
      if ((parameter = get_args(od, "minimizer"))) {
        if (!strncasecmp(parameter, "min", 3)) {
          strcpy(od->qmd.minimizer, TINKER_MINIMIZE_EXE);
        }
        else if (strncasecmp(parameter, "opt", 3)) {
          tee_error(od, run_type, overall_line_num,
            "The only allowed values for the \"minimizer\" "
            "parameter are \"MINIMIZE\" and \"OPTIMIZE\".\n%s", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.min_grad = 0.001;
      if ((parameter = get_args(od, "min_grad"))) {
        sscanf(parameter, "%lf", &(od->qmd.min_grad));
        if (od->qmd.min_grad <= 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "per-atom gradient", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.min_maxiter = 1000;
      if ((parameter = get_args(od, "min_maxiter"))) {
        sscanf(parameter, "%d", &(od->qmd.min_maxiter));
        if (od->qmd.min_maxiter <= 0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "minimization maximum iterations", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.rmsd = 0.2;
      if ((parameter = get_args(od, "rmsd"))) {
        sscanf(parameter, "%lf", &(od->qmd.rmsd));
        if (od->qmd.rmsd < 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "rmsd", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.range = 3.0;
      if ((parameter = get_args(od, "range"))) {
        sscanf(parameter, "%lf", &(od->qmd.range));
        if (od->qmd.range < 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "range", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.temperature = 1000.0;
      if ((parameter = get_args(od, "temperature"))) {
        sscanf(parameter, "%lf", &(od->qmd.temperature));
        if (od->qmd.temperature < 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "temperature", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.runs = 200;
      if ((parameter = get_args(od, "runs"))) {
        sscanf(parameter, "%d", &(od->qmd.runs));
        if (od->qmd.runs <= 0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "runs", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (od->qmd.runs > 999999) {
          tee_printf(od, "The number of QMD runs will be "
            "limited to 999999.\n\n");
          od->qmd.runs = 999999;
        }
      }
      od->qmd.window = 10.0;
      if ((parameter = get_args(od, "window"))) {
        sscanf(parameter, "%lf", &(od->qmd.window));
        od->qmd.window = safe_rint
          (od->qmd.window * 1.0e03) / 1.0e03;
        if (od->qmd.window < 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "window", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.time_step = 1.0;
      if ((parameter = get_args(od, "time_step"))) {
        sscanf(parameter, "%lf", &(od->qmd.time_step));
        od->qmd.time_step = safe_rint(od->qmd.time_step);
        if (od->qmd.time_step < 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "time_step", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.options &= (~QMD_KEEP_INITIAL);
      if ((parameter = get_args(od, "keep_initial"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->qmd.options |= QMD_KEEP_INITIAL;
        }
      }
      od->qmd.options &= (~QMD_GBSA);
      if ((parameter = get_args(od, "gbsa"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->qmd.options |= QMD_GBSA;
        }
      }
      od->qmd.options &= (~QMD_REMOVE_FOLDER);
      if ((parameter = get_args(od, "remove_qmd_folder"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->qmd.options |= QMD_REMOVE_FOLDER;
        }
      }
      od->qmd.diel_const = 1.0;
      if ((!(od->qmd.options & QMD_GBSA)) && (parameter = get_args(od, "diel_const"))) {
        sscanf(parameter, "%lf", &(od->qmd.diel_const));
        if (od->qmd.diel_const <= 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "dielectric constant", QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "QMD", line_orig);
        tee_flush(od);
        memset(file_basename, 0, BUF_LEN);
        if ((parameter = get_args(od, "qmd_dir"))) {
          strcpy(od->qmd.qmd_dir, parameter);
          absolute_path(od->qmd.qmd_dir);
          result = 0;
          if (!dexist(od->qmd.qmd_dir)) {
            #ifndef WIN32
            result = mkdir(od->qmd.qmd_dir, S_IRWXU | S_IRGRP | S_IROTH);
            #else
            result = mkdir(od->qmd.qmd_dir);
            #endif
          }
        }
        else {
          strcpy(file_basename, od->file[MOLFILE_IN]->name);
          get_dirname(file_basename);
          result = open_perm_dir(od, file_basename, "qmd_dir", od->qmd.qmd_dir);
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_DIR_CANNOT_BE_CREATED, od->qmd.qmd_dir,
            QMD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        tee_printf(od, M_INPUT_OUTPUT_LOG_DIR, "qmd_dir", od->qmd.qmd_dir);
        result = call_obenergy(od, O3_MMFF94);
        switch (result) {
          case FL_CANNOT_CREATE_CHANNELS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PIPE, QMD_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CHDIR:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, QMD_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_PROCESS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PROCESS, "OpenBabel", QMD_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, QMD_FAILED);
          return PARSE_INPUT_ERROR;

          case OPENBABEL_ERROR:
          tee_error(od, run_type, overall_line_num,
            E_PROGRAM_ERROR, "OpenBabel");
          if ((od->file[TEMP_LOG]->handle = fopen
            (od->file[TEMP_LOG]->name, "rb"))) {
            while (fgets(buffer, BUF_LEN,
              od->file[TEMP_LOG]->handle)) {
              buffer[BUF_LEN - 1] = '\0';
              tee_printf(od, "%s", buffer);
            }
            fclose(od->file[TEMP_LOG]->handle);
            od->file[TEMP_LOG]->handle = NULL;
            tee_printf(od, "\n%s", QMD_FAILED);
          }
          else {
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", QMD_FAILED);
          }
          return PARSE_INPUT_ERROR;
        }
        result = qmd(od);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        free_threads(od);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, QMD_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create",
            od->error_code, QMD_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join",
            od->error_code, QMD_FAILED);
          return PARSE_INPUT_ERROR;
        }
        for (i = 0, result = 0; i < od->grid.object_num; ++i) {
          if (od->al.task_list[i]->code) {
            result = 1;
            tee_printf(od, "Object %4d:\n", i + 1);
            switch (od->al.task_list[i]->code) {
              case FL_OUT_OF_MEMORY:
              tee_printf(od, E_OUT_OF_MEMORY, "");
              break;
              
              case FL_CANNOT_READ_OB_OUTPUT:
              tee_printf(od, E_ERROR_IN_READING_OB_OUTPUT,
                od->al.task_list[i]->string, "");
              break;
              
              case FL_CANNOT_READ_OUT_FILE:
              tee_printf(od, E_ERROR_IN_READING_TINKER_OUTPUT,
                od->al.task_list[i]->string, "");
              tee_printf(od, E_PROGRAM_ERROR, "TINKER");
              sprintf(log_fd.name, "%s%c%04d%c%04d_%06d.log",
                od->qmd.qmd_dir, SEPARATOR,
                od->al.mol_info[od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM]]->object_id,
                SEPARATOR,
                od->al.mol_info[od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM]]->object_id,
                od->al.task_list[i]->data[TEMPLATE_CONF_NUM]);
              if ((log_fd.handle = fopen(log_fd.name, "rb"))) {
                while (fgets(buffer, BUF_LEN, log_fd.handle)) {
                  buffer[BUF_LEN - 1] = '\0';
                  tee_printf(od, "%s", buffer);
                }
                fclose(log_fd.handle);
                log_fd.handle = NULL;
              }
              else {
                tee_printf(od, E_CANNOT_READ_PROGRAM_LOG, "TINKER", "");
              }
              break;
              
              case FL_CANNOT_WRITE_INP_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_TINKER_INPUT,
                od->al.task_list[i]->string, "");
              break;

              case FL_ABNORMAL_TERMINATION:
              tee_printf(od, E_CALCULATION_ERROR, "TINKER computations", "");
              break;

              case FL_CANNOT_READ_MOL_FILE:
              tee_printf(od, E_ERROR_IN_READING_MOL_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_WRITE_TEMP_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_TEMP_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_CREATE_CHANNELS:
              tee_printf(od, E_CANNOT_CREATE_PIPE, "");
              break;

              case FL_CANNOT_CREATE_PROCESS:
              tee_printf(od, E_CANNOT_CREATE_PROCESS, "TINKER", "");
              break;

              case FL_UNKNOWN_ATOM_TYPE:
              tee_printf(od, E_UNKNOWN_ATOM_TYPE, "atom", "");
              break;

              case FL_CANNOT_CREATE_SCRDIR:
              tee_printf(od, E_TINKER_DIR_CANNOT_BE_CREATED, "working", "");
              break;
            }
            O3_ERROR_PRINT(od->al.task_list[i]);
          }
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_CALCULATION_ERROR, "QMD procedures", QMD_FAILED);
          return PARSE_INPUT_ERROR;
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "QMD");
        tee_flush(od);
        free_array(od->al.task_list);
        od->al.task_list = NULL;
      }
    }
    else if (!strcasecmp(arg->me[0], "energy")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, ENERGY_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, ENERGY_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->qmd.tinker_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_TINKER_PATH, ENERGY_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      sprintf(od->qmd.tinker_prm_path, "%s%c..%cshare%ctinker",
        od->qmd.tinker_exe_path, SEPARATOR, SEPARATOR, SEPARATOR);
      found = 0;
      if (dexist(od->qmd.tinker_prm_path)) {
        sprintf(buffer, "%s%c%s", od->qmd.tinker_prm_path,
          SEPARATOR, TINKER_MMFF94_PRM_FILE);
        if (fexist(buffer)) {
          sprintf(buffer, "%s%c%s", od->qmd.tinker_prm_path,
            SEPARATOR, TINKER_MMFF94S_PRM_FILE);
          found = fexist(buffer);
        }
      }
      if (!found) {
        memset(od->qmd.tinker_prm_path, 0, BUF_LEN);
      }
      if (!(parameter = get_args(od, "prm_dir"))) {
        if (!(od->qmd.tinker_prm_path[0])) {
          tee_error(od, run_type, overall_line_num,
            E_MISSING_DIR, "TINKER MMFF parameter",
            "prm_dir", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else {
        strncpy(od->qmd.tinker_prm_path, parameter, BUF_LEN - 1);
        od->qmd.tinker_prm_path[BUF_LEN - 1] = '\0';
        absolute_path(od->qmd.tinker_prm_path);
      }
      strcpy(od->qmd.minimizer, TINKER_ANALYZE_EXE);
      if ((parameter = get_args(od, "tool"))) {
        if (!strncasecmp(parameter, "min", 3)) {
          strcpy(od->qmd.minimizer, TINKER_MINIMIZE_EXE);
        }
        else if (!strncasecmp(parameter, "opt", 3)) {
          strcpy(od->qmd.minimizer, TINKER_OPTIMIZE_EXE);
        }
        else if (strncasecmp(parameter, "ana", 3)) {
          tee_error(od, run_type, overall_line_num,
            "The only allowed values for the \"tool\" "
            "parameter are \"ANALYZE\", \"MINIMIZE\" "
            "and \"OPTIMIZE\".\n%s", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.min_grad = 0.001;
      if ((parameter = get_args(od, "min_grad"))) {
        sscanf(parameter, "%lf", &(od->qmd.min_grad));
        if (od->qmd.min_grad <= 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "per-atom gradient", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.min_maxiter = 1000;
      if ((parameter = get_args(od, "min_maxiter"))) {
        sscanf(parameter, "%d", &(od->qmd.min_maxiter));
        if (od->qmd.min_maxiter <= 0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "minimization maximum iterations", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.rmsd = 0.2;
      if ((parameter = get_args(od, "rmsd"))) {
        sscanf(parameter, "%lf", &(od->qmd.rmsd));
        if (od->qmd.rmsd < 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "rmsd", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.range = 3.0;
      if ((parameter = get_args(od, "range"))) {
        sscanf(parameter, "%lf", &(od->qmd.range));
        if (od->qmd.range < 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "range", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->qmd.options &= (~QMD_GBSA);
      if ((parameter = get_args(od, "gbsa"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->qmd.options |= QMD_GBSA;
        }
      }
      od->qmd.options &= (~QMD_ALIGN);
      if ((parameter = get_args(od, "align"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->qmd.options |= QMD_ALIGN;
        }
      }
      od->qmd.options &= (~QMD_REMOVE_DUPLICATES);
      if ((parameter = get_args(od, "remove_duplicates"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->qmd.options |= QMD_REMOVE_DUPLICATES;
        }
      }
      od->qmd.options &= (~QMD_DONT_SUPERPOSE);
      if ((parameter = get_args(od, "superpose"))) {
        if (!strncasecmp(parameter, "n", 1)) {
          od->qmd.options |= QMD_DONT_SUPERPOSE;
        }
      }
      od->qmd.diel_const = 1.0;
      if ((!(od->qmd.options & QMD_GBSA)) && (parameter = get_args(od, "diel_const"))) {
        sscanf(parameter, "%lf", &(od->qmd.diel_const));
        if (od->qmd.diel_const <= 0.0) {
          tee_error(od, run_type, overall_line_num,
            E_POSITIVE_NUMBER, "dielectric constant", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->align.type = 0;
      if ((parameter = get_args(od, "candidate"))) {
        if (!strncasecmp(parameter, "multi", 5)) {
          od->align.type |= ALIGN_MULTICONF_CANDIDATE_BIT;
        }
        else if (strncasecmp(parameter, "sing", 4)) {
          tee_error(od, run_type, overall_line_num,
            E_ONLY_SINGLE_MULTIPLE, ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
        if ((parameter = get_args(od, "src_dir"))) {
          strcpy(od->qmd.src, parameter);
          absolute_path(od->qmd.src);
          if (!dexist(od->qmd.src)) {
            tee_error(od, run_type, overall_line_num,
              E_DIR_NOT_EXISTING, od->qmd.src, ENERGY_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            E_MISSING_DIR, "source",
            "src_dir", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(parameter = get_args(od, "dest_dir"))) {
          tee_error(od, run_type, overall_line_num,
            E_MISSING_DIR, "destination",
            "dest_dir", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(od->qmd.dest, parameter);
        absolute_path(od->qmd.dest);
      }
      else {
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_FOR_MODIFIED_COORDINATES,
            "optimized", ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(od->qmd.dest, parameter);
        absolute_path(od->qmd.dest);
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "ENERGY", line_orig);
        tee_flush(od);
        if (od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
          if (!dexist(od->qmd.dest)) {
            #ifndef WIN32
            result = mkdir(od->qmd.dest, S_IRWXU | S_IRGRP | S_IROTH);
            #else
            result = mkdir(od->qmd.dest);
            #endif
          }
        }
        result = call_obenergy(od, O3_MMFF94);
        switch (result) {
          case FL_CANNOT_CREATE_CHANNELS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PIPE, ENERGY_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CHDIR:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, ENERGY_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_PROCESS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PROCESS, "OpenBabel", ENERGY_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ENERGY_FAILED);
          return PARSE_INPUT_ERROR;

          case OPENBABEL_ERROR:
          tee_error(od, run_type, overall_line_num,
            E_PROGRAM_ERROR, "OpenBabel");
          if ((od->file[TEMP_LOG]->handle = fopen
            (od->file[TEMP_LOG]->name, "rb"))) {
            while (fgets(buffer, BUF_LEN,
              od->file[TEMP_LOG]->handle)) {
              buffer[BUF_LEN - 1] = '\0';
              tee_printf(od, "%s", buffer);
            }
            fclose(od->file[TEMP_LOG]->handle);
            od->file[TEMP_LOG]->handle = NULL;
            tee_printf(od, "\n%s", ENERGY_FAILED);
          }
          else {
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", ENERGY_FAILED);
          }
          return PARSE_INPUT_ERROR;
        }
        if (od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT) {
          od->pel.conf_population[ANY_DB] = int_perm_resize
            (od->pel.conf_population[ANY_DB], od->grid.object_num);
          if (!(od->pel.conf_population[ANY_DB])) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, ENERGY_FAILED);
            return PARSE_INPUT_ERROR;
          }
          result = check_conf_db(od, od->qmd.src,
            ANY_DB, &wrong_object_num, &wrong_conf_num);
          if (result) {
            sprintf(buffer, "%s%c%04d.sdf", od->align.template_conf_dir,
              SEPARATOR, od->al.mol_info[wrong_object_num]->object_id);
            switch (result) {
              case PREMATURE_EOF:
              tee_error(od, run_type, overall_line_num, 
                E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "SDF",
                buffer, ENERGY_FAILED);
              return PARSE_INPUT_ERROR;

              case N_ATOM_BOND_MISMATCH:
              tee_error(od, run_type, overall_line_num,
                E_CONF_DB_ATOM_BONDS_NOT_MATCHING,
                wrong_object_num + 1, buffer, wrong_conf_num,
                ENERGY_FAILED);
              return PARSE_INPUT_ERROR;
            }
          }
        }
        result = open_temp_dir(od, od->temp_dir, "energy_scratch", od->align.align_scratch);
        if (result) {
          tee_error(od, run_type, overall_line_num, E_TEMP_DIR_CANNOT_BE_CREATED,
            od->align.align_scratch, ENERGY_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        tee_printf(od, "The energy_scratch directory is:\n%s\n\n", od->align.align_scratch);
        result = energy(od);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        free_threads(od);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ENERGY_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create",
            od->error_code, ENERGY_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join",
            od->error_code, ENERGY_FAILED);
          return PARSE_INPUT_ERROR;
          
          case ERROR_MERGING_FILES:
          tee_error(od, run_type, overall_line_num,
            "Error while merging temporary files in %s",
            od->align.align_scratch, ENERGY_FAILED);
          return PARSE_INPUT_ERROR;
        }
        for (i = 0, result = 0; i < od->grid.object_num; ++i) {
          if (od->al.task_list[i]->code) {
            result = 1;
            tee_printf(od, "Object %4d:\n", i + 1);
            switch (od->al.task_list[i]->code) {
              case FL_OUT_OF_MEMORY:
              tee_printf(od, E_OUT_OF_MEMORY, "");
              break;
              
              case FL_CANNOT_READ_OB_OUTPUT:
              tee_printf(od, E_ERROR_IN_READING_OB_OUTPUT,
                od->al.task_list[i]->string, "");
              break;
              
              case FL_CANNOT_READ_OUT_FILE:
              tee_printf(od, E_ERROR_IN_READING_TINKER_OUTPUT,
                od->al.task_list[i]->string, "");
              tee_printf(od, E_PROGRAM_ERROR, "TINKER");
              sprintf(log_fd.name, "%s%c%04d%c%04d_%06d.log",
                od->qmd.qmd_dir, SEPARATOR,
                od->al.mol_info[od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM]]->object_id,
                SEPARATOR,
                od->al.mol_info[od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM]]->object_id,
                od->al.task_list[i]->data[TEMPLATE_CONF_NUM]);
              if ((log_fd.handle = fopen(log_fd.name, "rb"))) {
                while (fgets(buffer, BUF_LEN, log_fd.handle)) {
                  buffer[BUF_LEN - 1] = '\0';
                  tee_printf(od, "%s", buffer);
                }
                fclose(log_fd.handle);
                log_fd.handle = NULL;
              }
              else {
                tee_printf(od, E_CANNOT_READ_PROGRAM_LOG, "TINKER", "");
              }
              break;
              
              case FL_CANNOT_WRITE_INP_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_TINKER_INPUT,
                od->al.task_list[i]->string, "");
              break;

              case FL_ABNORMAL_TERMINATION:
              tee_printf(od, E_CALCULATION_ERROR, "TINKER computations", "");
              break;

              case FL_CANNOT_READ_MOL_FILE:
              tee_printf(od, E_ERROR_IN_READING_MOL_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_WRITE_TEMP_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_TEMP_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_CREATE_CHANNELS:
              tee_printf(od, E_CANNOT_CREATE_PIPE, "");
              break;

              case FL_CANNOT_CREATE_PROCESS:
              tee_printf(od, E_CANNOT_CREATE_PROCESS, "TINKER", "");
              break;

              case FL_UNKNOWN_ATOM_TYPE:
              tee_printf(od, E_UNKNOWN_ATOM_TYPE, "atom", "");
              break;

              case FL_CANNOT_CREATE_SCRDIR:
              tee_printf(od, E_TINKER_DIR_CANNOT_BE_CREATED, "working", "");
              break;
            }
            O3_ERROR_PRINT(od->al.task_list[i]);
          }
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_CALCULATION_ERROR, "ENERGY calculations", ENERGY_FAILED);
          return PARSE_INPUT_ERROR;
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "ENERGY");
        tee_flush(od);
        free_array(od->al.task_list);
        od->al.task_list = NULL;
      }
    }
    else if (!strcasecmp(arg->me[0], "filter")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, FILTER_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->align.pharao_exe[0])) {
        tee_error(od, run_type, overall_line_num,
          E_PHARAO_PATH, FILTER_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      strcpy(od->align.pharao_exe_path, od->align.pharao_exe);
      get_dirname(od->align.pharao_exe_path);
      if ((parameter = get_args(od, "hybrid"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->align.type |= ALIGN_TOGGLE_HYBRID_BIT;
        }
      }
      if ((parameter = get_args(od, "merge"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->align.type |= ALIGN_TOGGLE_MERGE_BIT;
        }
      }
      od->align.filter_type = FILTER_INTRA_CONF_DB_BIT;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "inter", 5)) {
          od->align.filter_type = FILTER_INTER_CONF_DB_BIT;
        }
        else if (strncasecmp(parameter, "intra", 5)) {
          tee_error(od, run_type, overall_line_num,
            "The only allowed values for the \"type\" "
            "parameter are \"INTRA\" and \"INTER\".\n%s",
            FILTER_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->align.level = 0.7;
      if ((parameter = get_args(od, "level"))) {
        sscanf(parameter, "%lf", &(od->align.level));
      }
      if (!(run_type & DRY_RUN)) {
        memset(od->align.template_conf_dir, 0, BUF_LEN);
        if (od->qmd.qmd_dir[0]) {
          strcpy(od->align.template_conf_dir, od->qmd.qmd_dir);
        }
        if ((parameter = get_args(od, "template_conf_dir"))) {
          strcpy(od->align.template_conf_dir, parameter);
        }
        if (!(od->align.template_conf_dir[0])) {
          tee_error(od, run_type, overall_line_num,
            E_MISSING_TEMPLATE_CONF_DIR, FILTER_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        absolute_path(od->align.template_conf_dir);
        strcpy(file_basename, od->align.pharao_exe);
        get_dirname(file_basename);
        if ((parameter = get_args(od, "filter_conf_dir"))) {
          strcpy(od->align.filter_conf_dir, parameter);
          absolute_path(od->align.filter_conf_dir);
          result = 0;
          if (!dexist(od->align.filter_conf_dir)) {
            #ifndef WIN32
            result = mkdir(od->align.filter_conf_dir, S_IRWXU | S_IRGRP | S_IROTH);
            #else
            result = mkdir(od->align.filter_conf_dir);
            #endif
          }
        }
        else {
          result = open_perm_dir(od, file_basename, "filter_conf_dir", od->align.filter_conf_dir);
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_DIR_CANNOT_BE_CREATED, od->align.filter_conf_dir,
            FILTER_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        tee_printf(od, "The filter_conf_dir where filtered SDF conformational "
          "databases will be put is:\n%s\n\n", od->align.filter_conf_dir);
        od->pel.numberlist[OBJECT_LIST] = int_perm_resize
          (od->pel.numberlist[OBJECT_LIST], od->grid.object_num);
        if (!(od->pel.numberlist[OBJECT_LIST])) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, FILTER_FAILED);
          return PARSE_INPUT_ERROR;
        }
        for (i = 0; i < od->grid.object_num; ++i) {
          od->pel.numberlist[OBJECT_LIST]->pe[i] = i + 1;
        }
        if ((parameter = get_args(od, "pool"))) {
          if (!strncasecmp(parameter, "exclude_test_set", 12)) {
            od->pel.numberlist[OBJECT_LIST] = int_perm_resize
              (od->pel.numberlist[OBJECT_LIST], od->active_object_num);
            if (!(od->pel.numberlist[OBJECT_LIST])) {
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, FILTER_FAILED);
              return PARSE_INPUT_ERROR;
            }
            for (i = 0, j = 0; i < od->grid.object_num; ++i) {
              if (get_object_attr(od, i, ACTIVE_BIT)) {
                od->pel.numberlist[OBJECT_LIST]->pe[j] = i + 1;
                ++j;
              }
            }
          }
          else if (!strncasecmp(parameter, "include_list", 7)) {
            result = parse_synonym_lists(od, "FILTER", FILTER_FAILED,
              (1 << OBJECT_LIST) | (1 << ID_LIST), &list_type,
              OBJECT_LIST, run_type, overall_line_num);
            switch (result) {
              case PARSE_INPUT_RECOVERABLE_ERROR:
              fail = !(run_type & INTERACTIVE_RUN);
              continue;

              case PARSE_INPUT_ERROR:
              return PARSE_INPUT_ERROR;
            }
          }
        }
        od->pel.conf_population[TEMPLATE_DB] = int_perm_resize
          (od->pel.conf_population[TEMPLATE_DB], od->grid.object_num);
        if (!(od->pel.conf_population[TEMPLATE_DB])) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, FILTER_FAILED);
          return PARSE_INPUT_ERROR;
        }
        if (!dexist(od->align.template_conf_dir)) {
          tee_error(od, run_type, overall_line_num, E_DIR_NOT_EXISTING,
          od->align.template_conf_dir, FILTER_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = check_conf_db(od, od->align.template_conf_dir,
          TEMPLATE_DB, &wrong_object_num, &wrong_conf_num);
        if (result) {
          sprintf(buffer, "%s%c%04d.sdf", od->align.template_conf_dir,
            SEPARATOR, od->al.mol_info[wrong_object_num]->object_id);
          switch (result) {
            case CANNOT_READ_ORIGINAL_SDF:
            tee_error(od, run_type, overall_line_num, 
              E_ERROR_IN_READING_CONF_FILE, FILTER_FAILED);
            return PARSE_INPUT_ERROR;

            case PREMATURE_EOF:
            tee_error(od, run_type, overall_line_num, 
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "SDF",
              buffer, FILTER_FAILED);
            return PARSE_INPUT_ERROR;

            case N_ATOM_BOND_MISMATCH:
            tee_error(od, run_type, overall_line_num,
              E_CONF_DB_ATOM_BONDS_NOT_MATCHING,
              wrong_object_num + 1, buffer, wrong_conf_num,
              FILTER_FAILED);
            return PARSE_INPUT_ERROR;
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "FILTER", line_orig);
        tee_flush(od);
        result = open_temp_dir(od, od->temp_dir, "filter_scratch", od->align.align_scratch);
        if (result) {
          tee_error(od, run_type, overall_line_num, E_TEMP_DIR_CANNOT_BE_CREATED,
            od->align.align_scratch, FILTER_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        tee_printf(od, "The filter_scratch directory is:\n%s\n\n", od->align.align_scratch);
        result = filter(od);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        len = 0;
        free_threads(od);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, FILTER_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create",
            od->error_code, FILTER_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join",
            od->error_code, FILTER_FAILED);
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_READ_SDF_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_SDF_FILE, od->task.string, FILTER_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_WRITE_SDF_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_SDF_FILE, od->task.string, FILTER_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, od->task.string, FILTER_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
        }
        for (i = 0; (result && (result != NOTHING_TO_DO_FILTER) && (i < od->align.n_tasks)); ++i) {
          if (od->al.task_list[i]->code) {
            template_object_num = od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM];
            template_conf_num = od->al.task_list[i]->data[TEMPLATE_CONF_NUM];
            sprintf(log_fd.name, "%s%c%04d_phar_%s.log",
              od->align.align_scratch, SEPARATOR,
              od->al.mol_info[template_object_num]->object_id,
              (result == ERROR_IN_FILTER_EXTRACT_SPLIT) ? "extract" : "compare");
            switch (result) {
              case ERROR_IN_FILTER_EXTRACT_SPLIT:
              tee_printf(od, "Error during pharmacophore extraction "
                "from template object %d (ID %d)\n", template_object_num + 1,
                od->al.mol_info[template_object_num]->object_id);
              break;
              
              case ERROR_IN_FILTER_INTRA:
              tee_printf(od, E_INTRA_INTER_FILTRATION, "intra", template_object_num + 1,
                od->al.mol_info[template_object_num]->object_id);
              break;
              
              case ERROR_IN_FILTER_INTER:
              tee_printf(od, E_INTRA_INTER_FILTRATION, "inter", template_object_num + 1,
                od->al.mol_info[template_object_num]->object_id);
              break;
            }
            switch (od->al.task_list[i]->code) {
              case FL_OUT_OF_MEMORY:
              tee_printf(od, E_OUT_OF_MEMORY, "");
              break;

              case FL_CANNOT_READ_PHARAO_OUTPUT:
              tee_printf(od, E_ERROR_IN_READING_PHARAO_OUTPUT,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_CREATE_SCRDIR:
              tee_printf(od, E_SCRDIR_CANNOT_BE_CREATED,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_FIND_CONF:
              tee_printf(od, E_ERROR_IN_FINDING_CONF,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_WRITE_TEMP_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_TEMP_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_READ_SDF_FILE:
              tee_printf(od, E_ERROR_IN_READING_SDF_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_WRITE_SDF_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_SDF_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_CREATE_CHANNELS:
              tee_printf(od, E_CANNOT_CREATE_PIPE, "");
              break;

              case FL_CANNOT_CREATE_PROCESS:
              tee_printf(od, E_CANNOT_CREATE_PROCESS, "PHARAO", "");
              break;

              case FL_PHARAO_ERROR:
              tee_error(od, run_type, overall_line_num,
                E_PROGRAM_ERROR, "PHARAO");
              if ((log_fd.handle = fopen(log_fd.name, "rb"))) {
                while (fgets(buffer, BUF_LEN, log_fd.handle)) {
                  buffer[BUF_LEN - 1] = '\0';
                  tee_printf(od, "%s", buffer);
                }
                tee_printf(od, "\n");
                fclose(log_fd.handle);
              }
              else {
                tee_error(od, run_type, overall_line_num,
                  E_CANNOT_READ_PROGRAM_LOG, "PHARAO", "");
              }
              break;
            }
            O3_ERROR_PRINT(od->al.task_list[i]);
          }
        }
        if (result) {
          if (result == NOTHING_TO_DO_FILTER) {
            tee_printf(od, "SDF databases for all selected template objects "
              "are already present in the %s folder; therefore, "
              "nothing was done.\n\n", od->align.filter_conf_dir);
          }
          else {
            tee_error(od, run_type, overall_line_num, E_CALCULATION_ERROR,
              (result == ERROR_IN_FILTER_EXTRACT_SPLIT)
              ? "pharmacophore extractions with PHARAO"
              : "pharmacophore comparisons with PHARAO", FILTER_FAILED);
            return PARSE_INPUT_ERROR;
          }
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "FILTER");
        tee_flush(od);
        free_array(od->al.task_list);
        od->al.task_list = NULL;
      }
    }
    else if (!strcasecmp(arg->me[0], "align")) {
      gettimeofday(&start, NULL);
      od->align.n_tasks = 0;
      memset(od->align.template_file, 0, BUF_LEN);
      memset(od->align.candidate_file, 0, BUF_LEN);
      memset(od->align.template_dir, 0, BUF_LEN);
      memset(od->align.candidate_dir, 0, BUF_LEN);
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, ALIGN_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, ALIGN_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->align.type = ALIGN_ATOMBASED_BIT;
      if ((parameter = get_args(od, "type"))) {
        if ((!strncasecmp(parameter, "phar", 4))
          || (!strncasecmp(parameter, "mix", 3))) {
          if (!(od->align.pharao_exe[0])) {
            tee_error(od, run_type, overall_line_num,
              E_PHARAO_PATH, ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          strcpy(od->align.pharao_exe_path, od->align.pharao_exe);
          get_dirname(od->align.pharao_exe_path);
        }
        if (!strncasecmp(parameter, "phar", 4)) {
          od->align.type = ALIGN_PHARAO_BIT;
        }
        else if (!strncasecmp(parameter, "rand", 4)) {
          od->align.type = ALIGN_RANDOM_BIT;
        }
        else if (!strncasecmp(parameter, "atom", 4)) {
          od->align.type = ALIGN_ATOMBASED_BIT;
        }
        else if (!strncasecmp(parameter, "mix", 3)) {
          od->align.type = ALIGN_ATOMBASED_BIT | ALIGN_MIXED_BIT;
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "The only allowed values for the \"type\" "
            "parameter are \"ATOM\", \"PHAR\", \"MIXED\" "
            "and RANDOM.\n%s", ALIGN_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (od->align.type & (ALIGN_PHARAO_BIT | ALIGN_MIXED_BIT)) {
        if ((parameter = get_args(od, "hybrid"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            od->align.type |= ALIGN_TOGGLE_HYBRID_BIT;
          }
        }
        if ((parameter = get_args(od, "merge"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            od->align.type |= ALIGN_TOGGLE_MERGE_BIT;
          }
        }
      }
      if (!(od->align.type & ALIGN_RANDOM_BIT)) {
        if ((parameter = get_args(od, "template"))) {
          if (!strncasecmp(parameter, "multi", 5)) {
            od->align.type |= ALIGN_MULTICONF_TEMPLATE_BIT;
          }
          else if (!strncasecmp(parameter, "iter", 4)) {
            od->align.type |= ALIGN_ITERATIVE_TEMPLATE_BIT;
          }
          else if (strncasecmp(parameter, "sing", 4)) {
            tee_error(od, run_type, overall_line_num,
              "The only allowed values for the \"template\" "
              "parameter are \"SINGLE\", \"MULTIPLE\" "
              "and \"ITERATIVE\".\n%s", ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        if ((parameter = get_args(od, "keep_best_template"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            od->align.type |= ALIGN_KEEP_BEST_TEMPLATE_BIT;
          }
        }
        if ((parameter = get_args(od, "fast"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            od->align.type |= ALIGN_TOGGLE_LOOP_BIT;
          }
        }
        od->align.gold = ALIGN_GOLD_COEFFICIENT;
        if ((parameter = get_args(od, "gold"))) {
          sscanf(parameter, "%lf", &(od->align.gold));
          if (od->align.gold < 0.0) {
            tee_error(od, run_type, overall_line_num,
              E_POSITIVE_NUMBER, "gold parameter", ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        od->align.max_iter = DEFAULT_MAX_ITER_ALIGN;
        if ((parameter = get_args(od, "max_iter"))) {
          sscanf(parameter, "%d", &(od->align.max_iter));
          if (od->align.max_iter < 0) {
            tee_error(od, run_type, overall_line_num,
              E_POSITIVE_NUMBER, "max_iter parameter", ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        od->align.max_fail = DEFAULT_MAX_FAIL_ALIGN;
        if ((parameter = get_args(od, "max_fail"))) {
          sscanf(parameter, "%d", &(od->align.max_fail));
          if (od->align.max_fail < 0) {
            tee_error(od, run_type, overall_line_num,
              E_POSITIVE_NUMBER, "max_fail parameter", ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        if ((parameter = get_args(od, "print_rmsd"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            od->align.type |= ALIGN_PRINT_RMSD_BIT;
          }
        }
        if ((parameter = get_args(od, "candidate"))) {
          if (!strncasecmp(parameter, "multi", 5)) {
            od->align.type |= ALIGN_MULTICONF_CANDIDATE_BIT;
          }
          else if (strncasecmp(parameter, "sing", 4)) {
            tee_error(od, run_type, overall_line_num,
              E_ONLY_SINGLE_MULTIPLE, ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        if ((parameter = get_args(od, "template_file"))) {
          strcpy(od->align.template_file, parameter);
          absolute_path(od->align.template_file);
        }
        if ((parameter = get_args(od, "candidate_file"))) {
          strcpy(od->align.candidate_file, parameter);
          absolute_path(od->align.candidate_file);
        }
      }
      result = parse_synonym_lists(od, "ALIGN", ALIGN_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST), &list_type,
        MOST_ACTIVE_PERCENT, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      most_active_percent = 0;
      if ((!list_type) && (!(od->align.type & ALIGN_RANDOM_BIT))) {
        most_active_percent = 25;
        if ((parameter = get_args(od, "most_active_percent"))) {
          sscanf(parameter, "%d", &most_active_percent);
          if ((most_active_percent < 1) || (most_active_percent > 100)) {
            tee_error(od, run_type, overall_line_num,
              "The allowed range for the percent of most "
              "active training set objects used as templates "
              "for the alignment is 1 - 100.\n%s", ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
      }
      if (!(run_type & DRY_RUN)) {
        if (od->align.template_file[0] && (!fexist(od->align.template_file))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            od->align.template_file, ALIGN_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (od->align.candidate_file[0] && (!fexist(od->align.candidate_file))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            od->align.candidate_file, ALIGN_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (most_active_percent) {
          ref_y_var = 1;
          if (!(od->y_vars)) {
            tee_error(od, run_type, overall_line_num,
              "Either an object_list/id_list is indicated, or "
              "at least one dependent variable is "
              "imported so that the most active %d%% "
              "may be selected from the training set.\n%s",
              most_active_percent, ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          else if (od->y_vars > 1) {
            if ((parameter = get_args(od, "ref_y_var"))) {
              sscanf(parameter, "%d", &ref_y_var);
              if ((ref_y_var < 1) || (ref_y_var > od->y_vars)) {
                tee_error(od, run_type, overall_line_num,
                  E_Y_VAR_ALLOWED_RANGE,
                  od->y_vars, ALIGN_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            else {
              tee_error(od, run_type, overall_line_num,
                "Please indicate through the \"ref_y_var\" "
                "parameter which Y variable should ",
                "be used to select the most active %d%% "
                "from the training set.\n%s",
                most_active_percent, ALIGN_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
          }
          criterion = HIGHEST_VALUE;
          if ((parameter = get_args(od, "criterion"))) {
            if (!strncasecmp(parameter, "low", 3)) {
              criterion = LOWEST_VALUE;
            }
            else if (strncasecmp(parameter, "high", 4)) {
              tee_error(od, run_type, overall_line_num,
                "Only \"LOWEST\" and \"HIGHEST\" "
                "types are allowed for the \"criterion\" "
                "parameter.\n%s", ALIGN_FAILED);
            }
          }
          od->pel.activity_rank = int_perm_resize
            (od->pel.activity_rank, od->active_object_num);
          if (!(od->pel.activity_rank)) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, ALIGN_FAILED);
            return PARSE_INPUT_ERROR;
          }
          od->vel.activity_list = double_vec_resize
            (od->vel.activity_list, od->active_object_num);
          if (!(od->vel.activity_list)) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, ALIGN_FAILED);
            return PARSE_INPUT_ERROR;
          }
          for (i = 0, j = 0; i < od->grid.object_num; ++i) {
            if (get_object_attr(od, i, ACTIVE_BIT)) {
              od->vel.activity_list->ve[j] =
                get_y_value(od, i, ref_y_var - 1, 0);
              ++j;
            }
          }
          double_vec_sort(od->vel.activity_list, od->pel.activity_rank);
          len = (int)((double)most_active_percent / 100.0 * (double)(od->active_object_num));
          begin = ((criterion == HIGHEST_VALUE) ? (od->active_object_num - len) : 0);
          od->pel.numberlist[OBJECT_LIST] = int_perm_resize
            (od->pel.numberlist[OBJECT_LIST], len);
          if (!(od->pel.numberlist[OBJECT_LIST])) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, ALIGN_FAILED);
            return PARSE_INPUT_ERROR;
          }
          for (i = 0; i < len; ++i) {
            j = 0;
            active_object_num = 0;
            while ((active_object_num < od->pel.activity_rank->pe[i + begin])
              && (j < od->grid.object_num)) {
              if (get_object_attr(od, j, ACTIVE_BIT)) {
                ++active_object_num;
              }
              ++j;
            }
            od->pel.numberlist[OBJECT_LIST]->pe[(begin ? len - i - 1 : i)] = j + 1;
          }
          qsort(od->pel.numberlist[OBJECT_LIST]->pe, len, sizeof(int), compare_integers);
        }
        if (od->align.type & (ALIGN_MULTICONF_TEMPLATE_BIT | ALIGN_MULTICONF_CANDIDATE_BIT)) {
          memset(od->align.template_conf_dir, 0, BUF_LEN);
          memset(od->align.candidate_conf_dir, 0, BUF_LEN);
          if (od->qmd.qmd_dir[0]) {
            strcpy(od->align.template_conf_dir, od->qmd.qmd_dir);
            strcpy(od->align.candidate_conf_dir, od->qmd.qmd_dir);
          }
          j = 0;
          if ((parameter = get_args(od, "conf_dir"))) {
            strcpy(od->align.template_conf_dir, parameter);
            strcpy(od->align.candidate_conf_dir, parameter);
            j += 2;
          }
          if ((parameter = get_args(od, "template_conf_dir"))) {
            strcpy(od->align.template_conf_dir, parameter);
            ++j;
          }
          if ((parameter = get_args(od, "candidate_conf_dir"))) {
            strcpy(od->align.candidate_conf_dir, parameter);
            ++j;
          }
          if ((j < 2) && (((!(od->align.template_conf_dir[0]))
            && (od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT))
            || ((!(od->align.candidate_conf_dir[0]))
            && (od->align.type & ALIGN_MULTICONF_CANDIDATE_BIT)))) {
            tee_error(od, run_type, overall_line_num,
              E_MISSING_CONF_DIR, ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          else if (j > 2) {
            tee_error(od, run_type, overall_line_num,
              "Please supply either the \"conf_dir\" "
              "keyword alone, or the \"template_conf_dir\" "
              "and \"candidate_conf_dir\" keywords.\n%s",
              ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          absolute_path(od->align.template_conf_dir);
          absolute_path(od->align.candidate_conf_dir);
          for (i = 0; i <= 1; ++i) {
            if (!(od->align.type
              & (i ? ALIGN_MULTICONF_CANDIDATE_BIT : ALIGN_MULTICONF_TEMPLATE_BIT))) {
              continue;
            }
            od->pel.conf_population[i] = int_perm_resize
              (od->pel.conf_population[i], od->grid.object_num);
            if (!(od->pel.conf_population[i])) {
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, ALIGN_FAILED);
              return PARSE_INPUT_ERROR;
            }
            conf_dir = (i ? od->align.candidate_conf_dir : od->align.template_conf_dir);
            if ((!conf_dir[0]) || (!dexist(conf_dir))) {
              fail = 1;
              break;
            }
            result = check_conf_db(od, conf_dir, (i ? CANDIDATE_DB : TEMPLATE_DB),
              &wrong_object_num, &wrong_conf_num);
            if (result) {
              switch (result) {
                case CANNOT_READ_ORIGINAL_SDF:
                tee_error(od, run_type, overall_line_num, 
                  E_ERROR_IN_READING_CONF_FILE,
                  od->task.string, ALIGN_FAILED);
                return PARSE_INPUT_ERROR;

                case PREMATURE_EOF:
                tee_error(od, run_type, overall_line_num, 
                  E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "SDF",
                  od->task.string, ALIGN_FAILED);
                return PARSE_INPUT_ERROR;

                case N_ATOM_BOND_MISMATCH:
                tee_error(od, run_type, overall_line_num,
                  E_CONF_DB_ATOM_BONDS_NOT_MATCHING,
                  wrong_object_num + 1, od->task.string,
                  wrong_conf_num, ALIGN_FAILED);
                return PARSE_INPUT_ERROR;
              }
            }
          }
          if (fail) {
            tee_error(od, run_type, overall_line_num,
              E_DIR_NOT_EXISTING, conf_dir, ALIGN_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "ALIGN", line_orig);
        tee_flush(od);
        if (!(od->align.type & ALIGN_RANDOM_BIT)) {
          tee_printf(od, "The objects which will be used as "
            "alignment templates are the following:\n");
          for (i = 0, template_num = 0; i < od->pel.numberlist[OBJECT_LIST]->size; ++i) {
            template_object_num = od->pel.numberlist[OBJECT_LIST]->pe[i] - 1;
            tee_printf(od, "%d%s", template_object_num + 1,
              ((i != (od->pel.numberlist[OBJECT_LIST]->size - 1)) ? ", " : "\n\n"));
            template_num += ((od->align.type & ALIGN_MULTICONF_TEMPLATE_BIT)
              ? od->pel.conf_population[TEMPLATE_DB]->pe[template_object_num] : 1);
          }
        }
        memset(file_basename, 0, BUF_LEN);
        if ((parameter = get_args(od, "align_dir"))) {
          strcpy(od->align.align_dir, parameter);
          absolute_path(od->align.align_dir);
          result = 0;
          if (!dexist(od->align.align_dir)) {
            #ifndef WIN32
            result = mkdir(od->align.align_dir, S_IRWXU | S_IRGRP | S_IROTH);
            #else
            result = mkdir(od->align.align_dir);
            #endif
          }
        }
        else {
          result = open_perm_dir(od, od->default_folder, "align_dir", od->align.align_dir);
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_DIR_CANNOT_BE_CREATED, od->align.align_dir,
            ALIGN_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        tee_printf(od, M_INPUT_OUTPUT_LOG_DIR, "align_dir", od->align.align_dir);
        result = open_temp_dir(od, od->temp_dir, "align_scratch", od->align.align_scratch);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_DIR_CANNOT_BE_CREATED, od->align.align_scratch,
            ALIGN_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        tee_printf(od, "The align_scratch directory is:\n%s\n\n", od->align.align_scratch);
        result = call_obenergy(od, O3_MMFF94);
        switch (result) {
          case FL_CANNOT_CREATE_CHANNELS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PIPE, ALIGN_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CHDIR:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, ALIGN_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_PROCESS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PROCESS, "OpenBabel", ALIGN_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ALIGN_FAILED);
          return PARSE_INPUT_ERROR;

          case OPENBABEL_ERROR:
          tee_error(od, run_type, overall_line_num,
            E_PROGRAM_ERROR, "OpenBabel");
          if ((od->file[TEMP_LOG]->handle = fopen
            (od->file[TEMP_LOG]->name, "rb"))) {
            while (fgets(buffer, BUF_LEN,
              od->file[TEMP_LOG]->handle)) {
              buffer[BUF_LEN - 1] = '\0';
              tee_printf(od, "%s", buffer);
            }
            fclose(od->file[TEMP_LOG]->handle);
            od->file[TEMP_LOG]->handle = NULL;
            tee_printf(od, "\n%s", ALIGN_FAILED);
          }
          else {
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", ALIGN_FAILED);
          }
          return PARSE_INPUT_ERROR;
        }
        for (i = 0; i < od->grid.object_num ; ++i) {
          od->al.mol_info[i]->done = 0;
        }
        if (od->align.type & ALIGN_RANDOM_BIT) {
          result = align_random(od);
        }
        else {
          if (alloc_threads(od)) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, ALIGN_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
          }
          result = ((od->align.type & ALIGN_ITERATIVE_TEMPLATE_BIT) ? align_iterative(od) : align(od));
        }
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        len = 0;
        free_threads(od);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_DIRECTORY:
          tee_error(od, run_type, overall_line_num,
            E_TEMP_DIR_CANNOT_BE_CREATED,
            od->align.align_dir, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create",
            od->error_code, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join",
            od->error_code, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE,
            od->task.string, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          case FL_CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE,
            od->task.string, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_ALIGNED_SDF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_SDF_FILE,
            od->task.string, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_ORIGINAL_SDF:
          case FL_CANNOT_READ_SDF_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_SDF_FILE,
            od->task.string, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_READ_MOL_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_MOL_FILE,
            od->task.string, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_READ_OB_OUTPUT:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_OB_OUTPUT,
            od->task.string, ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_UNKNOWN_ATOM_TYPE:
          tee_error(od, run_type, overall_line_num,
            E_UNKNOWN_ATOM_TYPE, "atom", ALIGN_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
        }
        for (i = 0; result && (i < od->align.n_tasks); ++i) {
          if (od->al.task_list[i]->code) {
            template_object_num = od->al.task_list[i]->data[TEMPLATE_OBJECT_NUM];
            template_conf_num = od->al.task_list[i]->data[TEMPLATE_CONF_NUM];
            moved_object_num = od->al.task_list[i]->data[MOVED_OBJECT_NUM];
            moved_conf_num = od->al.task_list[i]->data[MOVED_CONF_NUM];
            switch (result) {
              case ERROR_EXTRACTING_PHAR:
              tee_printf(od, "Error during pharmacophore extraction "
                "from object %d (ID %d)", template_object_num + 1,
                od->al.mol_info[template_object_num]->object_id);
              if (template_conf_num != -1) {
                tee_printf(od, ", conformation %d", template_conf_num + 1);
                sprintf(log_fd.name, "%s%c%04d_%06d_phar.log", od->align.align_scratch,
                  SEPARATOR, od->al.mol_info[template_object_num]->object_id,
                  template_conf_num + 1);
              }
              else {
                sprintf(log_fd.name, "%s%c%04d_phar.log", od->align.align_scratch,
                  SEPARATOR, od->al.mol_info[template_object_num]->object_id);
              }
              tee_printf(od, "\n");
              break;
              
              case ERROR_IN_ALIGNMENT:
              tee_printf(od, "Error during alignment ");
              if (moved_object_num != -1) {
                tee_printf(od, "of object %d (ID %d)", moved_object_num + 1,
                  od->al.mol_info[moved_object_num]->object_id);
              }
              else {
                tee_printf(od, "of the dataset");
              }
              if (moved_conf_num != -1) {
                tee_printf(od, ", conformation %d", moved_conf_num + 1);
              }
              if (template_object_num != -1) {
                tee_printf(od, " on template object %d (ID %d)", template_object_num + 1,
                    od->al.mol_info[template_object_num]->object_id);
                if (template_conf_num != -1) {
                  tee_printf(od, ", conformation %d", template_conf_num + 1);
                  sprintf(log_fd.name, "%s%c%04d_on_%04d_%06d.log",
                    od->align.align_scratch, SEPARATOR, od->al.mol_info[i]->object_id,
                    od->al.mol_info[template_object_num]->object_id, template_conf_num + 1);
                }
                else {
                  sprintf(log_fd.name, "%s%c%04d-%04d_on_%04d.log",
                    od->align.align_scratch, SEPARATOR,
                    od->al.mol_info[0]->object_id,
                    od->al.mol_info[od->grid.object_num - 1]->object_id,
                    od->al.mol_info[template_object_num]->object_id);
                }
              }
              tee_printf(od, "\n");
              break;
            }
            switch (od->al.task_list[i]->code) {
              case FL_OUT_OF_MEMORY:
              tee_printf(od, E_OUT_OF_MEMORY, "");
              break;

              case FL_CANNOT_READ_PHARAO_OUTPUT:
              tee_printf(od, E_ERROR_IN_READING_PHARAO_OUTPUT,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_READ_MOL_FILE:
              tee_printf(od, E_ERROR_IN_READING_MOL_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_READ_TEMP_FILE:
              tee_printf(od, E_ERROR_IN_READING_TEMP_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_WRITE_TEMP_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_TEMP_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_READ_SDF_FILE:
              tee_printf(od, E_ERROR_IN_READING_SDF_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_WRITE_SDF_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_SDF_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_FIND_CONF:
              tee_printf(od, E_ERROR_IN_FINDING_CONF,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_READ_CONF_FILE:
              tee_printf(od, E_ERROR_IN_READING_CONF_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_CREATE_CHANNELS:
              tee_printf(od, E_CANNOT_CREATE_PIPE, "");
              break;

              case FL_CANNOT_CREATE_PROCESS:
              tee_printf(od, E_CANNOT_CREATE_PROCESS, "PHARAO", "");
              break;

              case FL_PHARAO_ERROR:
              tee_error(od, run_type, overall_line_num,
                E_PROGRAM_ERROR, "PHARAO");
              if ((log_fd.handle = fopen(log_fd.name, "rb"))) {
                while (fgets(buffer, BUF_LEN, log_fd.handle)) {
                  buffer[BUF_LEN - 1] = '\0';
                  tee_printf(od, "%s", buffer);
                }
                tee_printf(od, "\n");
                fclose(log_fd.handle);
              }
              else {
                tee_error(od, run_type, overall_line_num,
                  E_CANNOT_READ_PROGRAM_LOG, "PHARAO", "");
              }
              break;
            }
            O3_ERROR_PRINT(od->al.task_list[i]);
          }
        }
        if (result == ERROR_EXTRACTING_PHAR) {
          tee_error(od, run_type, overall_line_num, E_CALCULATION_ERROR,
            "pharmacophore extractions with PHARAO", ALIGN_FAILED);
          return PARSE_INPUT_ERROR;
        }
        else if (result == ERROR_IN_ALIGNMENT) {
          tee_error(od, run_type, overall_line_num, E_CALCULATION_ERROR,
            "ligand alignments", ALIGN_FAILED);
          return PARSE_INPUT_ERROR;
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "ALIGN");
        tee_flush(od);
        free_array(od->al.task_list);
        od->al.task_list = NULL;
      }
    }
    else if (!strcasecmp(arg->me[0], "compare")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, COMPARE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, COMPARE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num,
          E_SPECIFY_STRUCTURE_FILE, "compare", COMPARE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      memset(od->file[MOLFILE_IN]->name, 0, BUF_LEN);
      strcpy(od->file[MOLFILE_IN]->name, parameter);
      absolute_path(od->file[MOLFILE_IN]->name);
      memset(od->file[ASCII_IN]->name, 0, BUF_LEN);
      type = 0;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "block", 5)) {
          type |= BLOCK_COMPARE;
        }
        else if (strncasecmp(parameter, "pair", 4)) {
          tee_error(od, run_type, overall_line_num,
            "Only \"PAIRWISE\", and \"BLOCK\" "
            "types are allowed "
            "for the COMPARE keyword.\n%s",
            COMPARE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if ((parameter = get_args(od, "aligned"))) {
        strcpy(od->file[ASCII_IN]->name, parameter);
        absolute_path(od->file[ASCII_IN]->name);
      }
      if (!(run_type & DRY_RUN)) {
        if (!fexist(od->file[MOLFILE_IN]->name)) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            parameter, COMPARE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "COMPARE", line_orig);
        tee_flush(od);
        if (!(od->file[MOLFILE_IN]->handle =
          fopen(od->file[MOLFILE_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            parameter, COMPARE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        memcpy(&od_comp, od, sizeof(O3Data));
        od_comp.al.mol_info = NULL;
        od_comp.pymol.use_pymol = 0;
        result = open_temp_dir(&od_comp, NULL, "comp_mol_dir", od_comp.field.mol_dir);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_DIR_CANNOT_BE_CREATED, od_comp.field.mol_dir,
            COMPARE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = parse_sdf(&od_comp, ALLOC_MOL_INFO_BIT, NULL);
        if (od->file[MOLFILE_IN]->handle) {
          fclose(od->file[MOLFILE_IN]->handle);
          od->file[MOLFILE_IN]->handle = NULL;
        }
        if (result) {
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
        }
        switch (result) {
          case PREMATURE_EOF:
          tee_error(od, run_type, overall_line_num, 
            E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "sdf",
            od->file[MOLFILE_IN]->name, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;
        }
        if (od_comp.grid.object_num != od->grid.object_num) {
          tee_error(od, run_type, overall_line_num,
            E_OBJECT_NUMBER_NOT_MATCHING,
            od->file[MOLFILE_IN]->name, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;
        }
        for (i = 0; i <= 1; ++i) {
          result = call_obenergy(i ? &od_comp : od, O3_MMFF94);
          switch (result) {
            case FL_CANNOT_CREATE_CHANNELS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PIPE, COMPARE_FAILED);
            return PARSE_INPUT_ERROR;

            case FL_CANNOT_CHDIR:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, COMPARE_FAILED);
            return PARSE_INPUT_ERROR;

            case FL_CANNOT_CREATE_PROCESS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PROCESS, "OpenBabel", COMPARE_FAILED);
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, COMPARE_FAILED);
            return PARSE_INPUT_ERROR;

            case OPENBABEL_ERROR:
            tee_error(od, run_type, overall_line_num,
              E_PROGRAM_ERROR, "OpenBabel");
            if ((od->file[TEMP_LOG]->handle = fopen
              (od->file[TEMP_LOG]->name, "rb"))) {
              while (fgets(buffer, BUF_LEN,
                od->file[TEMP_LOG]->handle)) {
                buffer[BUF_LEN - 1] = '\0';
                tee_printf(od, "%s", buffer);
              }
              fclose(od->file[TEMP_LOG]->handle);
              od->file[TEMP_LOG]->handle = NULL;
              tee_printf(od, "\n%s", COMPARE_FAILED);
            }
            else {
              tee_error(od, run_type, overall_line_num,
                E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", COMPARE_FAILED);
            }
            return PARSE_INPUT_ERROR;
          }
        }
        for (i = 0, len = 0; i < od->grid.object_num; ++i) {
          if ((od->al.mol_info[i]->n_atoms != od_comp.al.mol_info[i]->n_atoms)
            || (od->al.mol_info[i]->n_bonds != od_comp.al.mol_info[i]->n_bonds)) {
            ++len;
            tee_printf(od, "Object %4d (%s):\n", i + 1, od->al.mol_info[i]->object_name);
            tee_printf(od, "Reference atoms, bonds: %d, %d\n",
              od->al.mol_info[i]->n_atoms, od->al.mol_info[i]->n_bonds);
            tee_printf(od, "In %s: %d, %d\n", od->file[MOLFILE_IN]->name,
              od_comp.al.mol_info[i]->n_atoms, od_comp.al.mol_info[i]->n_bonds);
          }
        }
        if (len) {
          tee_error(od, run_type, overall_line_num,
            E_OBJECT_ATOMS_BONDS_NOT_MATCHING,
            od->file[MOLFILE_IN]->name, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;
        }
        if (update_mol(&od_comp)) {
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = compare(od, &od_comp, type, VERBOSE_BIT);
        free_threads(od);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create",
            od->error_code, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join",
            od->error_code, COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_MOL_FILE, od->task.string,
            COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_ORIGINAL_SDF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_SDF_FILE, od->task.string,
            COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_ALIGNED_SDF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_SDF_FILE, od->task.string,
            COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_READ_MOL_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_MOL_FILE, od->task.string,
            COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_WRITE_SDF_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_SDF_FILE, od->task.string,
            COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_READ_OB_OUTPUT:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_OB_OUTPUT, od->task.string,
            COMPARE_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_UNKNOWN_ATOM_TYPE:
          tee_error(od, run_type, overall_line_num,
            E_UNKNOWN_ATOM_TYPE, "atom", COMPARE_FAILED);
          return PARSE_INPUT_ERROR;
        }
        for (i = 0, result = 0; i < od->grid.object_num; ++i) {
          if (od->al.task_list[i]->code) {
            result = 1;
            tee_printf(od, "Object %4d:\n", i + 1);
            switch (od->al.task_list[i]->code) {
              case FL_OUT_OF_MEMORY:
              tee_printf(od, E_OUT_OF_MEMORY, "");
              break;

              case FL_CANNOT_READ_MOL_FILE:
              tee_printf(od, E_ERROR_IN_READING_MOL_FILE,
                od->al.task_list[i]->string, "");
              break;

              case FL_CANNOT_WRITE_TEMP_FILE:
              tee_printf(od, E_ERROR_IN_WRITING_TEMP_FILE,
                od->al.task_list[i]->string, "");
              break;
            }
            O3_ERROR_PRINT(od->al.task_list[i]);
          }
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_CALCULATION_ERROR, "RMSD computations", COMPARE_FAILED);
          return PARSE_INPUT_ERROR;
        }
        free_array(od->al.task_list);
        od->al.task_list = NULL;
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "COMPARE");
        tee_flush(od);
        free_array(od_comp.al.mol_info);
        od_comp.al.mol_info = NULL;
        remove_recursive(od_comp.field.mol_dir);
      }
    }
    else if ((!strcasecmp(arg->me[0], "cd"))
      || (!strcasecmp(arg->me[0], "chdir"))) {
      memset(file_basename, 0, BUF_LEN);
      if ((parameter = get_args(od, "dir"))) {
        strncpy(file_basename, parameter, BUF_LEN - 1);
        file_basename[BUF_LEN - 1] = '\0';
      }
      ++command;
      tee_printf(od, M_TOOL_INVOKE, nesting, command, "CHDIR", line_orig);
      tee_flush(od);
      result = 0;
      if (file_basename[0]) {
        result = chdir(file_basename);
      }
      if (result) {
        tee_error(od, run_type, overall_line_num,
          E_CANNOT_CHANGE_DIR, file_basename, CHANGE_DIR_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      else {
        if (getcwd(file_basename, BUF_LEN - 1)) {
          tee_printf(od, "Current directory: %s\n\n", file_basename);
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "CHDIR");
          tee_flush(od);
        }
      }

    }
    else if (!strcasecmp(arg->me[0], "source")) {
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the input file you wish "
          "to source.\n%s",
          SOURCE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      source = &(od->mel.source);
      while (*source) {
        source = &((*source)->next);
      }
      *source = (FileDescriptor *)malloc(sizeof(FileDescriptor));
      if (!(*source)) {
        tee_error(od, run_type, overall_line_num, E_OUT_OF_MEMORY,
          SOURCE_FAILED);
        return PARSE_INPUT_ERROR;
      }
      memset(*source, 0, sizeof(FileDescriptor));
      strcpy((*source)->name, parameter);
      (*source)->handle = fopen((*source)->name, "rb");
      if (!((*source)->handle)) {
        tee_error(od, run_type, overall_line_num,
          E_FILE_CANNOT_BE_OPENED_FOR_READING,
          parameter, SOURCE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      ++command;
      tee_printf(od, M_TOOL_INVOKE, nesting, command, "SOURCE", line_orig);
      tee_flush(od);
      if (run_type & DRY_RUN) {
        result = parse_input(od, (*source)->handle, DRY_RUN);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_WHILE_SOURCING, (*source)->name,
            SOURCE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      result = parse_input(od, (*source)->handle, 0);
      if (result) {
        tee_error(od, run_type, overall_line_num,
          E_WHILE_SOURCING, (*source)->name,
          SOURCE_FAILED);
      }
      fclose((*source)->handle);
      (*source)->handle = NULL;
      free(*source);
      *source = NULL;
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        break;

        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;

        default:
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SOURCE");
        tee_flush(od);
      }
      continue;
    }
    else if (!strcasecmp(arg->me[0], "dataset")) {
      ++command;
      tee_printf(od, M_TOOL_INVOKE, nesting, command, "DATASET", line_orig);
      tee_flush(od);
      update_field_object_attr(od, VERBOSE_BIT);
      calc_active_vars(od, FULL_MODEL);
      tee_printf(od, M_TOOL_SUCCESS, nesting, command, "DATASET");
      tee_flush(od);
    }
    else if ((!strcasecmp(arg->me[0], "exit"))
      || (!strcasecmp(arg->me[0], "quit"))
      || (!strcasecmp(arg->me[0], "stop"))) {
      break;
    }
    else {
      tee_error(od, run_type, overall_line_num,
        "%s: unknown command.\n\n", arg->me[0]);
      if (!(run_type & INTERACTIVE_RUN)) {
        tee_printf(od, PACKAGE_NAME" failed.\n");
        fail = 1;
        continue;
      }
    }
  }
  if (fail) {
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  if (run_type & DRY_RUN) {
    for (i = 0; i < od->file_num; ++i) {
      if ((i != MAIN_INPUT) && (i != MAIN_OUTPUT)) {
        if (od->file[i]->handle) {
          fclose(od->file[i]->handle);
          od->file[i]->handle = NULL;
        }
      }
    }
    od->valid = 0;
    od->object_num = 0;
    rewind(input_stream);
  }
    
  return 0;
}


#define NUM_SYNONYM_LISTS    4

int parse_synonym_lists(O3Data *od, char *tool_name, char *tool_msg,
  int synonym_list, int *list_type, int default_list,
  int run_type, int overall_line_num)
{
  char *list_names[] = {
    "field_list",
    "object_list",
    "id_list",
    "struct_list"
  };
  char *list_elem[] = {
    "fields",
    "objects",
    "object IDs",
    "structure IDs"
  };
  char *parameter = NULL;
  char comma_hyphen_list[BUF_LEN];
  int i;
  int synonym = 0;
  int result = 0;


  *list_type = 0;
  memset(comma_hyphen_list, 0, BUF_LEN);
  if (default_list && (default_list != MOST_ACTIVE_PERCENT)) {
    *list_type = (1 << default_list);
    strcpy(comma_hyphen_list, "all");
  }
  for (i = 0; i < NUM_SYNONYM_LISTS; ++i) {
    if ((synonym_list & (1 << i))
      && (parameter = get_args(od, list_names[i]))) {
      *list_type = (1 << i);
      ++synonym;
      strcpy(comma_hyphen_list, parameter);
    }
  }
  if (((!default_list) && (!synonym)) || (synonym > 1)) {
    tee_error(od, run_type, overall_line_num,
      E_SUPPLY_OBJECT_ID_STRUCT_LIST,
      tool_name, tool_msg);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  if ((*list_type & ((1 << OBJECT_LIST) | (1 << ID_LIST)
    | (1 << STRUCT_LIST))) && (!(od->grid.object_num))) {
    tee_error(od, run_type, overall_line_num,
      E_NO_OBJECTS_PRESENT, tool_msg);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  if ((synonym_list & (1 << FROM_FILE)) && (parameter = get_args(od, "file"))) {
    strcpy(od->file[ASCII_IN]->name, parameter);
    if (!(od->file[ASCII_IN]->handle =
      fopen(od->file[ASCII_IN]->name, "rb"))) {
      tee_error(od, run_type, overall_line_num,
        E_FILE_CANNOT_BE_OPENED_FOR_READING,
        od->file[ASCII_IN]->name, tool_msg);
      return PARSE_INPUT_RECOVERABLE_ERROR;
    }
    *list_type |= (1 << FROM_FILE);
  }
  else if (*list_type) {
    result = parse_comma_hyphen_list_to_array(od,
      comma_hyphen_list, intlog2(*list_type));
    switch (result) {
      case OUT_OF_MEMORY:
      tee_error(od, run_type, overall_line_num,
        E_OUT_OF_MEMORY, tool_msg);
      return PARSE_INPUT_ERROR;

      case INVALID_LIST_RANGE:
      for (i = 0; i < NUM_SYNONYM_LISTS; ++i) {
        if (*list_type & (1 << i)) {
          tee_error(od, run_type, overall_line_num, E_LIST_PARSING,
            list_elem[i], tool_name, tool_msg);
          break;
        }
      }
      return PARSE_INPUT_RECOVERABLE_ERROR;
    }
    if (!(run_type & DRY_RUN)) {
      result = set(od, intlog2(*list_type), OPERATE_BIT, 1, SILENT);
      if (result) {
        for (i = 0; i < NUM_SYNONYM_LISTS; ++i) {
          if (*list_type & (1 << i)) {
            tee_error(od, run_type, overall_line_num,
              "The specified %s are out of range.\n%s",
              list_elem[i], tool_msg);
            break;
          }
        }
        return PARSE_INPUT_RECOVERABLE_ERROR;
      }
    }
  }
  
  return 0;
}
