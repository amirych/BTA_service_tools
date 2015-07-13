#include<iostream>
#include<fstream>
#include<libgen.h>
//#include<cstdlib>
#include<cstdio>
#include<list>

#include<boost/program_options.hpp>
#include<boost/filesystem.hpp>

using namespace std;

namespace po = boost::program_options;

#define ROTCEN_ERROR_OK 0
#define ROTCEN_ERROR_HELP 1
#define ROTCEN_ERROR_CMD 10
#define ROTCEN_ERROR_INPUT_LIST 20
#define ROTCEN_ERROR_UNKNOWN_OPT 30
#define ROTCEN_ERROR_INVALID_OPT_VALUE 40
#define ROTCEN_ERROR_INVALID_FILENAME 50
#define ROTCEN_ERROR_UNAVAILABLE_CMD 60
#define ROTCEN_ERROR_NOT_ENOUGH_FILES 70
#define ROTCEN_ERROR_APP_FAILED 80
#define ROTCEN_ERROR_CANNOT_CREATE_FILE 90


#define ROTCEN_SEX_PARAM_FILE "sex.param"

static string join_vector(vector<string> &vec)
{
    string str;
    for ( long i = 0; i < vec.size(); ++i) str += vec[i] + " ";
    return str;
}

int main(int argc, char* argv[])
{

    vector<float> sex_thresh(1,5.0); // sextractor's THRESH default value
    vector<float> match_tol(1,0.5); // default radius of coordinate matching. It is in arcsecs for astrometrical solution or
                                    // dimensionless value for 'match application' solution ('matchrad' parameter)

    vector<string> sex_pars = {"-DETECT_TYPE CCD","-SATUR_LEVEL 40000","-CHECKIMAGE_TYPE NONE", "-FILTER N"};

    vector<string> solve_field_pars = {"--no-plots"};


    vector<string> sex_cat_prefix = {"detect"};

    string input_list_filename;

    // commandline options and arguments

    po::options_description visible_opts("Allowed options");
    visible_opts.add_options()
        ("help,h", "produce help message")
        ("threshold,t", po::value<vector<float> >(&sex_thresh), "set threshold level for object detection (sextractor's THRESH keyword)")
        ("radius,r", po::value<vector<float> >(&match_tol), "radius of coordinate matching [arcsecs for astrometrical solution]")
        ("use-match,m","use of 'match' application instead of astrometry (explicitly set '-s' option)")
        ("use-sex,s","use of sextractor to detect objects (in case of astrometrical solution)")
        ("sex-pars",po::value<vector<string> >(), "sextractor's parameters")
        ("solve-filed-pars",po::value<vector<string> >(), "'solve-field' parameters");


    po::options_description hidden_opts("");
    hidden_opts.add_options()("input-file",po::value<string>()->required());

    po::options_description cmd_opts("");
    cmd_opts.add(visible_opts).add(hidden_opts);

    po::positional_options_description pos_arg;
    pos_arg.add("input-file", 1);

    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv).options(cmd_opts).positional(pos_arg).run(), vm);

        if ( vm.count("help") ) {
//            cout << "Usage: " << basename(argv[0]) << " [-t num] [-h] input_list\n\n";
            cout << "Usage: " << boost::filesystem::basename(argv[0]) << "[-h] [-t num] [-r num]    \
                                                                          [--sex-pars] [--solve-field-pars]        \
                                                                          [--use-match] [--use-sex] input_list\n\n";
            cout << visible_opts << "\n";
            return ROTCEN_ERROR_HELP;
        }

        po::notify(vm);
    } catch (boost::program_options::required_option& e) {
        cerr << "The input list of files is missed! Try '-h' option!\n";
        return ROTCEN_ERROR_INPUT_LIST;
    } catch (boost::program_options::unknown_option& e) {
        cerr << "Unknown commandline options! Try '-h' option!\n";
        return ROTCEN_ERROR_UNKNOWN_OPT;
    } catch (boost::program_options::invalid_option_value& e) {
        cerr << "Invalid option value! Try '-h' option!\n";
        return ROTCEN_ERROR_INVALID_OPT_VALUE;
    } catch(boost::program_options::error& e) {
        cerr << "Error in commandline options! Try '-h' option!\n";
        return ROTCEN_ERROR_CMD;
    }


    // parse commandline options

    if ( vm.count("threshold") ) {
        sex_thresh = vm["threshold"].as<vector<float> >();
        for ( int i = 0; i < sex_thresh.size(); ++i ) cout << "THRESH = " << sex_thresh[i] << "\n";
    } else {
        cout << "THRESH = " << sex_thresh.front() << ".\n";
    }

    if ( vm.count("radius") ) {
        match_tol = vm["radius"].as<vector<float> >();
    }

    if ( vm.count("sex-pars") ) {
        sex_pars = vm["sex-pars"].as<vector<string> >();
    }

    if ( vm.count("solve-field-pars") ) {
        solve_field_pars = vm["solve-filed-pars"].as<vector<string> >();
    }

    if (vm.count("input-file")) {
        cout << "Input files are: "
             << vm["input-file"].as<string>() << "\n";
        input_list_filename = vm["input-file"].as<string>();
    }


    bool use_match = false;
    bool use_sex = false;

    if ( vm.count("use-match") ) {
        int ret = system("match --help  >/dev/null 2>&1"); // try to run command 'match'
        int exit_code = WEXITSTATUS(ret);
        if ( ret == -1 || exit_code == 127 ) {
            cerr << "Application 'match' is not available!\n";
            return ROTCEN_ERROR_UNAVAILABLE_CMD;
        }

        use_match = true;

        sex_pars.push_back("-GAIN 1.0");
        sex_pars.push_back("-PARAMETERS_NAME " + ROTCEN_SEX_PARAM_FILE);
        sex_pars.push_back("-CATALOG_TYPE FITS_1.0");

        // create SExtractor's parameter file
        ofstream sex_param_file;
        sex_param_file.open(ROTCEN_SEX_PARAM_FILE);
        if ( !sex_param_file.good() ) {
            return ROTCEN_ERROR_CANNOT_CREATE_FILE;
        }


        sex_param_file << "NUMBER\n";
        sex_param_file << "X_IMAGE\n";
        sex_param_file << "Y_IMAGE\n";
        sex_param_file << "MAG_BEST\n";

        sex_param_file.close();

    } else { // use of 'solve-field' from astrometry.net
        int ret = system("solve-field --help  >/dev/null 2>&1"); // try to run command 'solve-field'
        int exit_code = WEXITSTATUS(ret);
        if ( ret == -1 || exit_code == 127 ) {
            cerr << "Application 'solve-field' is not available!\n";
            return ROTCEN_ERROR_UNAVAILABLE_CMD;
        }
    }

    if ( vm.count("use-sex") ) {
        int ret = system("sex  >/dev/null 2>&1"); // try to run command 'sex' (Bertin's sextractor)
        int exit_code = WEXITSTATUS(ret);
        if ( ret == -1 || exit_code == 127 ) {
            cerr << "Application 'sex' is not available!\n";
            return ROTCEN_ERROR_UNAVAILABLE_CMD;
        }
        use_sex = true;
        solve_field_pars.push_back("--use-sextractor");
    }


    list<string> input_files;
    string str;

    ifstream input_list_file;
    input_list_file.open(input_list_filename.c_str());

    if ( !input_list_file.good() ) {
        cerr << "Cannot find file of input frames list!\n";
        return ROTCEN_ERROR_INVALID_FILENAME;
    }

    try {
        // check input files
        while ( input_list_file.good() ) {
            input_list_file >> str;
            if ( !ifstream(str.c_str()).good() ) throw (int)ROTCEN_ERROR_INVALID_FILENAME;
            input_files.push_back(str);
        }
    } catch (int err) {
        input_list_file.close();
        cerr << "Cannot open " << str.c_str() << " input file!\n";
        return err;
    }

    input_list_file.close();

    if ( input_files.size() < 3 ) {
        cerr << "At least 3 files must be given in the input list!\n";
        return ROTCEN_ERROR_NOT_ENOUGH_FILES;
    }

    // run object detection and astrometry
    for ( auto it_file = input_files.begin(); it_file != input_files.end(); ++it_file ) {
        if ( use_match ) { // skip astrometry, just detect objects using sextractor
            boost::filesystem::path pp = *it_file;
            string path = pp.parent_path().string();
            string file = boost::filesystem::basename(*it_file);

            file = path + sex_cat_prefix.back() + file + ".fits";
            sex_pars.push_back("-CATALOG_NAME " + file);
            string cmd_str = "sex " + join_vector(sex_pars) + *it_file + " >/dev/null 2>&1";

            int ret = system(cmd_str.c_str()); // try to run command 'sex' (Bertin's sextractor)
            int exit_code = WEXITSTATUS(ret);
            if ( ret == -1 || exit_code != 0 ) {
                cerr << "Something wrong while run application 'sex'!\n";
                return ROTCEN_ERROR_APP_FAILED;
            }
        } else {

        }
    }


    return ROTCEN_ERROR_OK;
}
