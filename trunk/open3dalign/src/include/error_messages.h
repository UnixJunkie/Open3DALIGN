/*

error_messages.h

is part of

Open3DALIGN
-----------

An open-source software aimed at unsupervised molecular alignment

Copyright (C) 2010-2012 Paolo Tosco, Thomas Balle

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


char M_NUMBER_OF_CPUS[] =
  "The number of CPUs used by %s has been set to %d.\n\n";
char M_TOOL_INVOKE[] =
  "BGN COMMAND #%02d.%04d - %s tool was invoked as follows:\n"
  "> %s\n\n";
char M_TOOL_SUCCESS[] =
  "END COMMAND #%02d.%04d - %s tool succeeded.\n";
char M_EXECUTABLE[] =
  "The %s executable has been set to \"%s\".\n\n";
char M_EXECUTABLE_PATH[] =
  "The path to %s tools has been set to \"%s\".\n\n";
char M_INPUT_OUTPUT_LOG_DIR[] =
  "The %s directory where input, output and log "
  "files will be put is:\n\"%s\"\n\n";
char E_CANNOT_CHANGE_DIR[] =
  "Cannot change directory to \"%s\".\n%s";
char E_DIR_NOT_EXISTING[] =
  "The directory \"%s\" does not exist.\n%s";
char E_IMPORT_MOLFILE_FIRST[] =
  "Please import some molecule structures first.\n%s";
char E_SPECIFY_FILE_TYPE[] =
  "Please specify a type for "
  "the file you wish to import.\n%s";
char E_SPECIFY_STRUCTURE_FILE[] =
  "Please specify a SDF file containing the "
  "structure(s) you wish to %s.\n%s";
char E_CANNOT_ALTER_GRID_BOX[] =
  "The grid box cannot be changed/removed "
  "once fields are present.\n%s\n\n";
char E_CANNOT_ALTER_OBJECT_COORDINATES[] =
  "The object coordinates cannot be altered "
  "once fields are present.\n%s\n\n";
char E_NO_GRID_BOX_PRESENT[] =
  "No grid box is currently present.\n%s";
char E_POSITIVE_NUMBER[] =
  "The %s must "
  "be a positive number.\n%s";
char E_OPENBABEL_NOT_WORKING[] =
  "A simple OpenBabel test failed with the following error:\n";
char E_OPENBABEL_DATA_PLUGINS[] =
  "OpenBabel data/plugins could not be found.\n"
  "This issue can be fixed setting the "BABEL_DATADIR_ENV" and "
  BABEL_LIBDIR_ENV" environment "
  "variables to suitable values.\n%s\n";
char E_OPENBABEL_MISSING_OR_TOO_OLD[] =
  "OpenBabel binaries could not be found.\n"
  "Please consider downloading and installing the "
  "openbabel_for_open3dtools package "
  "for your system from the "PACKAGE_NAME" website.\n%s\n";
char E_OPENBABEL_PATH[] =
  "Please set the O3_BABEL_PATH environment variable "
  "or use the \"env babel_path\" keyword to indicate the path "
  "to OpenBabel executables.\n%s";
char E_PHARAO_NOT_WORKING[] =
  "A simple PHARAO test failed with the following error:\n";
char E_PHARAO_MISSING_OR_TOO_OLD[] =
  "PHARAO binary could not be found.\n"
  "Please consider downloading and installing the "
  "openbabel_for_open3dtools package "
  "for your system from the "PACKAGE_NAME" website.\n%s\n";
char E_PROGRAM_ERROR[] =
  "%s reported the error message which follows:\n";
char E_CANNOT_READ_PROGRAM_LOG[] =
  "Cannot read %s log file.\n%s";
char E_TINKER_PATH[] =
  "Please set the O3_TINKER_PATH environment variable "
  "or use the \"env tinker_path\" keyword to indicate the path "
  "to TINKER executables.\n%s";
char E_MISSING_DIR[] =
  "Please indicate the %s directory "
  "through the %s parameter.\n%s";
char E_PHARAO_PATH[] =
  "Please set the O3_PHARAO_PATH environment variable "
  "or use the \"env pharao_path\" keyword to indicate the path "
  "to the PHARAO executable.\n%s";
char E_ERROR_IN_WRITING_TINKER_INPUT[] =
  "Cannot write TINKER .key file \"%s\".\n%s";
char E_ERROR_IN_READING_TINKER_OUTPUT[] =
  "Cannot read TINKER output file \"%s\".\n%s";
char E_ERROR_IN_READING_OB_OUTPUT[] =
  "Cannot read OpenBabel output file \"%s\".\n%s";
char E_ERROR_IN_READING_MOL_FILE[] =
  "Cannot read MOL file \"%s\".\n%s";
char E_ERROR_IN_READING_PHARAO_OUTPUT[] =
  "Cannot read PHARAO output file \"%s\".\n%s";
char E_ERROR_IN_READING_SDF_FILE[] =
  "Cannot read SDF file \"%s\".\n%s";
char E_ERROR_IN_WRITING_SDF_FILE[] =
  "Cannot write SDF file \"%s\".\n%s";
char E_ERROR_IN_READING_CONF_FILE[] =
  "Cannot read SDF conformational database \"%s\".\n%s";
char E_ERROR_IN_FINDING_CONF[] =
  "Cannot find required conformation in the "
  "SDF conformational database \"%s\".\n%s";
char E_MISSING_CONF_DIR[] =
  "Please supply either the \"conf_dir\" "
  "parameter or the \"template_conf_dir\", "
  "\"candidate_conf_dir\" pair.\n%s";
char E_MISSING_TEMPLATE_CONF_DIR[] =
  "Please supply the \"template_conf_dir\" parameter.\n%s";
char E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING[] =
  "The temporary file \"%s\" cannot be opened for writing.\n%s";
char E_TEMP_DIR_CANNOT_BE_CREATED[] =
  "The folder \"%s\" cannot be opened for writing.\n%s";
char E_FILE_CANNOT_BE_OPENED_FOR_WRITING[] =
  "The file \"%s\" cannot be opened for writing.\n%s";
char E_FILE_CANNOT_BE_OPENED_FOR_READING[] =
  "The file \"%s\" cannot be opened for reading.\n%s";
char E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT[] =
  "The %s file \"%s\" appears to be corrupted "
  "or in the wrong format.\n%s";
char E_CANNOT_CREATE_PIPE[] =
  "Error in creating a pipe.\n%s";
char E_CANNOT_CREATE_PROCESS[] =
  "Error in spawning a %s process.\n%s";
char E_ERROR_IN_WRITING_TEMP_FILE[] =
  "Error in writing temporary file \"%s\".\n%s";
char E_ERROR_IN_READING_TEMP_FILE[] =
  "Error in reading temporary file \"%s\".\n%s";
char E_CALCULATION_ERROR[] =
  "An error occurred during one or more %s; "
  "please check your input/output/log files.\n%s";
char E_INTRA_INTER_FILTRATION[] =
  "Error during %s-conformational database filtration "
  "of template object %d (ID %d)\n";
char E_TINKER_DIR_CANNOT_BE_CREATED[] =
  "A %s folder for TINKER calculations cannot be created.\n%s";
char E_SCRDIR_CANNOT_BE_CREATED[] =
  "The scratch folder \"%s\" cannot be created.\n%s";
char E_SUPPLY_OBJECT_ID_STRUCT_LIST[] =
  "Please supply a single keyword among "
  "\"object_list\", \"id_list\" and eventually "
  "\"struct_list\" on which %s should operate.\n%s";
char E_STRUCT_ATTRIBUTE_ONLY[] =
  "The \"%s\" attribute can "
  "only be assigned to a structure_list, "
  "not to an object_list, when conformers "
  "are present in the dataset.\n%s";
char E_UNKNOWN_ATOM_TYPE[] =
  "Unknown %s type.\n%s";
char E_Y_VAR_LOW_SD[] =
  "The SD associated with the y variable(s) is too low.\n%s";
char E_OUT_OF_MEMORY[] =
  "Out of memory.\n%s";
char E_THREAD_ERROR[] =
  "Thread error; the code returned "
  "from pthread_%s() is %d.\n%s";
char E_LIST_PARSING[] =
  "Error while parsing the list of %s "
  "on which %s should operate.\n%s";
char E_ALLOWED_OBJECT_RANGE[] =
  "The allowed object range is 1 - %d.\n%s";
char E_CHECK_OBJECT_RANGE[] =
  "Please check your object ID range.\n%s";
char E_OBJECT_NUMBER_NOT_MATCHING[] =
  "The number of objects in %s does not match "
  "the number of objects currently loaded.\n%s";
char E_OBJECT_ATOMS_BONDS_NOT_MATCHING[] =
  "The number of %s in some objects in %s does not match "
  "the number of %s in the currently loaded objects.\n%s";
char E_CONF_DB_ATOM_BONDS_NOT_MATCHING[] =
  "There is a atom/bond number mismatch "
  "between the currently loaded object %d "
  "and its conformational database \"%s\" "
  "(conformer %d); please check your input.\n%s";
char E_ONLY_SINGLE_MULTIPLE[] =
  "The only allowed values for the \"candidate\" "
  "parameter are \"SINGLE\" and \"MULTIPLE\".\n%s";
char E_FILE_FOR_MODIFIED_COORDINATES[] =
  "Please specify a file where "
  "%s coordinates should be saved.\n%s";
char E_NOT_ENOUGH_Y_VARS[] =
  "In the file \"%s\" there must be at "
  "least %d biological activit%s.\n%s";
char E_WRONG_NUMBER_OF_Y_VARS[] =
  "In the file \"%s\" for all objects there "
  "must be the same number of biological "
  "activities.\n%s";
char E_CANNOT_FIND_Y_VAR_NAME[] =
  "In the file \"%s\" one of the requested y variable "
  "names does not exist.\n%s";
char E_Y_VAR_ALLOWED_RANGE[] =
  "Please choose biological activity number 1 - %d.\n%s";
char E_NO_OBJECTS_PRESENT[] =
  "No objects are currently present.\n%s";
char E_WHILE_SOURCING[] =
  "Error while sourcing \"%s\".\n%s";
char E_PROGRAM_EXIT[] =
  "\nPress ENTER to leave "PACKAGE_NAME".\n";
char O3_FAILED[] =
  PACKAGE_NAME" failed.\n";
char PARSE_INPUT_FAILED[] =
  "Input file parsing failed.\n";
char IMPORT_FAILED[] =
  "IMPORT failed.\n";
char ROTOTRANS_FAILED[] =
  "ROTOTRANS failed.\n";
char BOX_FAILED[] =
  "BOX failed.\n";
char REMOVE_BOX_FAILED[] =
  "REMOVE_BOX failed.\n";
char LOAD_FAILED[] =
  "LOAD failed.\n";
char SAVE_FAILED[] =
  "SAVE failed.\n";
char REMOVE_OBJECT_FAILED[] =
  "REMOVE_OBJECT failed.\n";
char REMOVE_Y_VARS_FAILED[] =
  "REMOVE_Y_VARS failed.\n";
char QMD_FAILED[] =
  "QMD failed.\n";
char ENERGY_FAILED[] =
  "ENERGY failed.\n";
char ALIGN_FAILED[] =
  "ALIGN failed.\n";
char COMPARE_FAILED[] =
  "COMPARE failed.\n";
char FILTER_FAILED[] =
  "FILTER failed.\n";
char SET_FAILED[] =
  "SET failed.\n";
char ENV_FAILED[] =
  "ENV failed.\n";
char CHANGE_DIR_FAILED[] =
  "CHANGE_DIR failed.\n";
char SOURCE_FAILED[] =
  "SOURCE failed.\n";
