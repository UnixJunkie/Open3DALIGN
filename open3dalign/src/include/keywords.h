/*

keywords.h

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


typedef struct O3ParameterData O3ParameterData;
typedef struct O3KeywordData O3KeywordData;

struct O3ParameterData {
  unsigned short type;
  char *parameter;
  char *choice[MAX_ARG];
};

struct O3KeywordData {
  char *keyword;
  O3ParameterData parameter_data[MAX_ARG];
};

O3KeywordData keyword_data[] =
{
  {
    "align",
    {
      {
        O3_PARAM_STRING, "type", {
          "ATOM",
          "MIXED",
          "PHAR",
          "RANDOM",
          NULL
        }
      }, {
        O3_PARAM_STRING, "object_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "most_active_percent", {
          "25",
          NULL
        }
      }, {
        O3_PARAM_REF_Y_VAR, "ref_y_var", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "criterion", {
          "HIGHEST",
          "LOWEST",
          NULL
        }
      }, {
        O3_PARAM_STRING, "hybrid", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "merge", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "template", {
          "SINGLE",
          "MULTI",
          NULL
        }
      }, {
        O3_PARAM_STRING, "keep_best_template", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "candidate", {
          "SINGLE",
          "MULTI",
          NULL
        }
      }, {
        O3_PARAM_STRING, "print_rmsd", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "conf_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "template_conf_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "candidate_conf_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "align_dir", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "box",
    {
      {
        O3_PARAM_STRING, "mode", {
          "SET",
          "GET",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "outgap", {
          "5.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "step", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "x_start", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "x_end", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_start", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_end", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_start", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_end", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "cd",
    {
      {
        O3_PARAM_DIRECTORY, "dir", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "chdir",
    {
      {
        O3_PARAM_DIRECTORY, "dir", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "compare",
    {
      {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "type", {
          "PAIRWISE",
          "BLOCK",
          NULL
        }
      }, {
        O3_PARAM_FILE, "aligned", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "dataset",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "env",
    {
      {
        O3_PARAM_NUMERIC, "random_seed", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "temp_dir", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "nice", {
          #ifndef WIN32
          "20",
          #else
          "BELOW_NORMAL",
          "ABOVE_NORMAL",
          "HIGH",
          "IDLE",
          "NORMAL",
          "REALTIME",
          #endif
          NULL
        }
      }, {
        O3_PARAM_N_CPUS, "n_cpus", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "babel_path", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "pharao", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "pymol", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "tinker_path", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "exit",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "stop",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "quit",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "energy",
    {
      {
        O3_PARAM_DIRECTORY, "prm_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "src_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "dest_dir", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "tool", {
          "ANALYZE",
          "OPTIMIZE",
          "MINIMIZE",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "min_maxiter", {
          "1000",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "min_grad", {
          "0.001",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "rmsd", {
          "0.2",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "range", {
          "3.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "diel_const", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "gbsa", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "align", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "remove_duplicates", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "superpose", {
          "YES",
          "NO",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "filter",
    {
      {
        O3_PARAM_STRING, "type", {
          "INTRA",
          "INTER",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "level", {
          "0.7",
          NULL
        }
      }, {
        O3_PARAM_STRING, "pool", {
          "EXCLUDE_TEST_SET",
          "INCLUDE_LIST",
          NULL
        }
      }, {
        O3_PARAM_STRING, "object_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "hybrid", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "merge", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "template_conf_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "filter_conf_dir", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "import",
    {
      {
        O3_PARAM_STRING, "type", {
          "SDF",
          "DEPENDENT",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "y_var_name", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "load",
    {
      {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "mode", {
          "APPEND",
          "NORMAL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "qmd",
    {
      {
        O3_PARAM_DIRECTORY, "prm_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "qmd_dir", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "keep_initial", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "minimizer", {
          "OPTIMIZE",
          "MINIMIZE",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "min_maxiter", {
          "1000",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "min_grad", {
          "0.001",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "rmsd", {
          "0.2",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "range", {
          "3.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "temperature", {
          "1000",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "runs", {
          "200",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "window", {
          "10.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "time_step", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "diel_const", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "gbsa", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "remove_qmd_folder", {
          "NO",
          "YES",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_box",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_object",
    {
      {
        O3_PARAM_STRING, "object_list", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "id_list", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_y_vars",
    {
      {
        O3_PARAM_STRING, "y_var_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "rototrans",
    {
      {
        O3_PARAM_NUMERIC, "x_trans", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_trans", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_trans", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "x_rot", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_rot", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_rot", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "save",
    {
      {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "set",
    {
      {
        O3_PARAM_STRING, "object_list", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "id_list", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "attribute", {
          "EXCLUDED",
          "TRAININGSET",
          "TESTSET",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "source",
    {
      {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {  // this is the terminator
    NULL,
    {
      {
        O3_PARAM_STRING, NULL, {
          NULL
        }
      }
    }
  }
};
