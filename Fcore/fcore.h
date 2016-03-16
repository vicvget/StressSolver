#ifndef FCORE_H

#define FCORE_H


// ограничение длины строковых буферов
#define STRLIMIT 1024
#define IS_ZERO(x) (fabs(x)<1e-6)

// enviroments
#define ENV_FRUND "FSHELL"
#define ENV_PATH  "PATH"
#define ENV_LIB	  "LIB"
#define ENV_INCLUDE	"INCLUDE"
#define ENV_LIC   "INTEL_LICENSE_FILE"

// directories
#define DIR_DATA "Data"

#ifdef GNUCPP
#define DIR_TOOLSET "Intel"
#else
#define DIR_TOOLSET "GNU"
#endif

#define NAME_DEFAULT "default"

#define MAX_LINE_LENGTH 1024

// extensions
#define EXT_FRM				".frm"
#define EXT_MODEL_ELEMENT   ".elf"
#define EXT_BODY_GEOMETRY   ".elg"
#define EXT_MULTIPHYSICS	".mph"
#define EXT_BODY_GRID	    ".grid"
#define EXT_STATUS		    ".sts"
#define EXT_MPH_RESULTS		".mpr"
#define EXT_MACRO			".mac"
#define EXT_SYNC			".sync"
#define EXT_DECOMPOSITION	".dec"
#define EXT_XML				".xml"
#define EXT_GV				".gv"
#define EXT_LBA				".lba"
#define EXT_ICO				".ico"

// files
#define FILE_TIRE_DEFAULT_PARAMS_1		"tfy2123.dat"
#define FILE_TIRE_DEFAULT_PARAMS_2		"tmz2123.dat"
#define FILE_TIRES_PATHS				"tire.cnt"
#define FILE_ROAD_PROFILES_PATHS		"way.cnt"
#define FILE_EXTRA_CONTROL_PARAMS		"epsilon.dat"
#define FILE_CONTROL_PARAMS				"uprf"
#define FILE_INITIAL_CONDITIONS			"initcond.inp"
#define FILE_DEFAULT_INITIAL_CONDITIONS "default.ico"

//#define FILE_CONFIG						"fcore.cfg"
#define FILE_CONFIG						"fcore.xml"
#define FILE_RESOURCES					"resources.xml"

#define FILE_REZR						"rezr"
#define FILE_INITDEFA					"initdefa"
#define FILE_MBS_RESULTS				"results.mbr"
#define FILE_MBS_ANIMATION_DUMP			"results.mba"
#define FILE_FADRES						"fadres.dat"
#define FILE_DEBUG_PARAMS				"rdebug.sys"

#define FILE_DEBUG_BC					"boundary.txt"
#define FILE_ADDRESSES_BC				"genexite.dat"


// Macros preprocessor files
#define FILE_RHL_SPRINGS				"springs.rhl"	// code 3, shell
#define FILE_MHL_SPRINGS				"springs.mhl"	// code 2, shell
#define FILE_MHL_DEFAULT				"defaul.mhl"	// default.mhl перезаписывается в CRMODEL, поэтому без последней буквы
#define FILE_MHL						"modhelp"		//FILE_MHL_DEFAULT + FILE_MHL_SPRINGS
#define FILE_RHL						"rashelp"		//FILE_RHL_DEFAULT + FILE_RHL_SPRINGS
#define FILE_PMODEL_SPRINGM				"springm.dat"
#define FILE_RHL_DEFAULT				"default.rhl"
#define FILE_PMODEL_DEFAULMUPR          "defaulm.upr"
#define FILE_PMODEL_PREPSYS				"preprsys.dat"
#define FILE_PMODEL_MCHAR				"mchar.dat"
#define FILE_PMODEL_RASHELPP            "rashelpp.rhl"  // TODO: не реализовано, выяснить, когда этот файл нужен и как с ним работать
#define FILE_SYNC_TABLE					"table.sync"
// Control code file
#define FILE_CPP_CODE					"control.cpp"
#define FILE_CPP_F_INTERFACE			"CppInterface.f"

#define FILE_INC_GCONTROL				"gcontrol.fi"
#define FILE_INC_GPRCONTROL				"gprcontrol.fi"
#define FILE_CNT_GCONTROL				"gcontrol.cnt"
#define FILE_CNT_GPRCONTROL				"gprcontrol.cnt"

#define FILE_CAD_MATRIXES				"CAD_Matrixes.dat"
#define FILE_WAYGEO						"waygeo.dat"

// VIV files
#define FILE_DEFAULT_CHART				"benafi.bnf"
#define ID_TIME_DEPENDENT				"tfunc"

//Task map files
#define FILE_TASK_MAP					"taskmap"

// Начальные условия для вспомогательного решателя
#define FILE_ICS					"default.ics"

#ifdef NIX
#	define FILE_SOLVER_DLL  "./solver32.so"
#else

#	define FILE_SOLVER "solver"

#	define EXT_DLL ".dll"
#	define EXT_LIB ".lib"
#	define EXT_EXP ".exp"

#	define FILE_SOLVER_DLL FILE_SOLVER EXT_DLL
#	define FILE_SOLVER_LIB FILE_SOLVER EXT_LIB
#	define FILE_SOLVER_EXP FILE_SOLVER EXT_EXP

#endif


#define MSG_SOLVER_START  "============================== SOLVER STARTED ================================="
#define MSG_SOLVER_FINISH "============================== SOLVER FINISHED ================================"

#define MSG_GEN_START     "============================ GENERATOR STARTED ================================"
#define MSG_GEN_FINISH    "============================ GENERATOR FINISHED ==============================="

#define MSG_PREP_START    "============================== PREPROC STARTED ================================"
#define MSG_PREP_FINISH   "============================== PREPROC FINISHED ==============================="

#define MSG_MAKE_START    "================================ MAKE STARTED ================================="
#define MSG_MAKE_FINISH   "================================ MAKE FINISHED ================================"

#define MSG_POST_START    "============================ POSTPROCESSING STARTED ============================"
#define MSG_POST_FINISH   "============================ POSTPROCESSING FINISHED ==========================="

#define MSG_OUT_CONTROL   "============================= OUTPUT CONTROL PARAMS ============================"

#define MSG_MPHCFG_START  "============================= MPH CONFIGURE STARTED ============================"
#define MSG_MPHCFG_FINISH "============================= MPH CONFIGURE FINISHED ==========================="
#define MSG_UNDERLINE     "--------------------------------------------------------------------------------"

#define FLAG_PREVENT_OVERHEAD 0
#define FLAG_FIXED_THREADS 2	// NOT USED IF ZERO
#define FLAG_FIXED_SUBSOLVER_THREADS 1 // NOT USED IF ZERO

#define FLAG_USE_PARALLEL_CONFIG 1
#define FILE_PARALLEL_CONFIG "parallel.cfg"

//#define GFORTRAN


#endif // FCORE_H