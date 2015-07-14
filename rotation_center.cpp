#include<iostream>
#include<fstream>
//#include<libgen.h>
//#include<cstdlib>
#include<cstdio>
#include<list>

#define BOOST_NO_CXX11_SCOPED_ENUMS // special definition to fix Boost's copy_file and -std=c++11 linking error
#include<boost/program_options.hpp>
#include<boost/filesystem.hpp>


#include"ascii_file.h"

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
#define ROTCEN_ERROR_BAD_DATA 100


static string ROTCEN_SEX_PARAM_FILE = "sex.param";
static string ROTCEN_MATCH_REF_CAT = "ref.cat";

static string ROTCEN_SEX_EXE = "sex";
static string ROTCEN_AST_EXE = "solve-field";
static string ROTCEN_MATCH_EXE = "match";


static string join_vector(vector<string> &vec)
{
    string str;
    for ( long i = 0; i < vec.size(); ++i) str += vec[i] + " ";
    return str;
}


static int run_external(string &cmd_str)
{
    int ret = system(cmd_str.c_str()); // try to run external command

    int exit_code = WEXITSTATUS(ret);

    if ( ret == -1 ) {
        return ret;
    }
    if ( exit_code != 0 ) {
        return exit_code;
    }

    return 0;
}


static int read_catalog(string &filename, size_t N_items, vector<vector<double> > &data)
{
    AsciiFile cat(filename.c_str());

    if ( !cat.good() ) return ROTCEN_ERROR_INVALID_FILENAME;

    data.clear();

    data = vector<vector<double> >(N_items);

    vector<double> vec;
    AsciiFile::AsciiFileFlag line_flag;


    while ( (line_flag = cat.ReadLine(N_items,&vec)) != AsciiFile::Eof ) {
        if ( line_flag == AsciiFile::InvalidData ) {
            cat.close();
            return ROTCEN_ERROR_BAD_DATA;
        }
        if ( line_flag == AsciiFile::DataString ) {
            for ( size_t i = 0; i < N_items; ++i ) data[i].push_back(vec[i]);
        }
    }

    cat.close();
    return ROTCEN_ERROR_OK;
}


int main(int argc, char* argv[])
{


    vector<float> sex_thresh(1,5.0); // sextractor's THRESH default value
    vector<float> match_tol;        // default radius of coordinate matching. It is in arcsecs for astrometrical solution or
                                    // dimensionless value for 'match application' solution ('matchrad' parameter)
                                    // The actuall default value is set below according to '--use-match'

    vector<string> sex_pars = {"-DETECT_TYPE CCD","-SATUR_LEVEL 40000","-CHECKIMAGE_TYPE NONE", "-FILTER N"};

    vector<string> solve_field_pars = {"--no-plots"};

    vector<string> match_pars = {"id1=0 id2=0 min_scale=0.9 max_scale=1.1 linear"}; // default 'match' commandline parameters

    vector<string> sex_cat_prefix = {"obj_"};

    vector<string> ast_prefix = {"wcs_"}; // astrometry-calibrated file prefix

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
        ("solve-filed-pars",po::value<vector<string> >(), "'solve-field' parameters")
        ("match-pars",po::value<vector<string> >(), "'match' parameters")
        ("dont-delete,d","do not delete temporary files");


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
            string head_str = "Usage: " + boost::filesystem::basename(argv[0]);
            string skip_str(head_str.length()+1,' ');

            cout << head_str << " [-h] [-t num] [-r num] [-d] [--solve-field-pars]\n" << skip_str <<
                                "[--use-match] [--match-pars] \n" << skip_str <<
                                "[--use-sex] [--sex-pars] input_list\n\n";

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

    bool dont_delete = false;
    if ( vm.count("dont-delete") ) {
        dont_delete = true;
    }

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
        sex_pars.erase(sex_pars.begin(),sex_pars.end());
        sex_pars.push_back(vm["sex-pars"].as<vector<string> >().back());
    }

    if ( vm.count("solve-field-pars") ) {
        solve_field_pars.erase(solve_field_pars.begin(),solve_field_pars.end());
        solve_field_pars.push_back(vm["solve-filed-pars"].as<vector<string> >().back());
    }

    if ( vm.count("match-pars") ) {
        match_pars.erase(match_pars.begin(),match_pars.end());
        match_pars.push_back(vm["match-pars"].as<vector<string> >().back());
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

        if ( !vm.count("radius") ) { // use default value
            match_tol = {4.0};
        }

        sex_pars.push_back("-GAIN 1.0");
        sex_pars.push_back("-PARAMETERS_NAME " + ROTCEN_SEX_PARAM_FILE);
        sex_pars.push_back("-CATALOG_TYPE ASCII");

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

        if ( !vm.count("radius") ) { // use default value
            match_tol = {0.5};
        }
    }


    list<string> input_files;
    list<string> sex_cats;
    string str;

    ifstream input_list_file;

    try {
        input_list_file.open(input_list_filename.c_str());
        if ( !input_list_file.good() ) {
            cerr << "Cannot find file of input frames list!\n";
            return ROTCEN_ERROR_INVALID_FILENAME;
        }

        // check input files

        while ( (input_list_file >> str).good() ) {
            if ( !boost::filesystem::exists(str) ) {
                cerr << "Cannot find " << str.c_str() << " input file!\n";
                throw (int)ROTCEN_ERROR_INVALID_FILENAME;
            }
            input_files.push_back(str);
        }

        input_list_file.close();

        if ( input_files.size() < 3 ) {
            cerr << "At least 3 files must be given in the input list!\n";
            throw (int)ROTCEN_ERROR_NOT_ENOUGH_FILES;
        }

        // run object detection and astrometry

        string cmd_pars = join_vector(sex_pars);

        cout << "\nObjects detection:\n";

        for ( auto it_file = input_files.begin(); it_file != input_files.end(); ++it_file ) {
            if ( use_match ) { // skip astrometry, just detect objects using sextractor
                boost::filesystem::path pp = *it_file;
                string path = pp.parent_path().string();
                string file = boost::filesystem::basename(*it_file);

                file = path + boost::filesystem::path::preferred_separator + sex_cat_prefix.back() + file + ".cat";
                string cmd_str = ROTCEN_SEX_EXE + " " + cmd_pars + " -CATALOG_NAME " +
                                 file + " " + *it_file + " >/dev/null 2>&1";

                cout << "  Run SExtractor for " + *it_file + " ... ";

                int ret = run_external(cmd_str); // try to run command 'sex' (Bertin's sextractor)
                if ( ret ) {
                    cout << "Failed!\n";
                    cerr << "Something wrong while run application 'sex'!\n";
                    throw (int)ROTCEN_ERROR_APP_FAILED;
                }

                cout << "OK!\n";

                sex_cats.push_back(file);
            } else { // perform astrometry

            }
        }

        // matching objects

        vector<vector<double> > obj_cat(3*input_files.size()); // NUMBER, X_IMAGE and Y_IMAGE columns
        vector<vector<double> > obj_id(input_files.size());
        vector<vector<double> > current_cat;


        if ( use_match ) { // use of 'match' application
            cout << "\nMatching objects (use of 'match' application):\n";

            auto it_file = sex_cats.begin();
            ++it_file; // point to the second catalog


            boost::filesystem::copy_file(sex_cats.front(),ROTCEN_MATCH_REF_CAT,boost::filesystem::copy_option::overwrite_if_exists);

            // read the first catalog
            int ret = read_catalog(sex_cats.front(), 3, current_cat);
            if ( ret != ROTCEN_ERROR_OK ) {
                cerr << "Something wrong while reading " << sex_cats.front() << " file!\n";
                throw ret;
            }
            obj_cat[0] = current_cat[0]; // NUMBER
            obj_cat[1] = current_cat[1]; // X_IMAGE
            obj_cat[2] = current_cat[2]; // Y_IMAGE

            string cmd_str;
            string pp = join_vector(match_pars);

            string matchedA = "matched.mtA";
            string matchedB = "matched.mtB";

            for ( long i_cat = 1; it_file != sex_cats.end(); ++it_file, ++i_cat ) {

                // read current catalog

                ret = read_catalog(*it_file, 3, current_cat);
                if ( ret != ROTCEN_ERROR_OK ) {
                    cerr << "Something wrong while reading " << *it_file << " file!\n";
                    throw ret;
                }
                obj_cat[i_cat*3] = current_cat[0]; // NUMBER
                obj_cat[i_cat*3+1] = current_cat[1]; // X_IMAGE
                obj_cat[i_cat*3+2] = current_cat[2]; // Y_IMAGE

                cmd_str = ROTCEN_MATCH_EXE + " " + ROTCEN_MATCH_REF_CAT + " 1 2 3 " + *it_file + " 1 2 3 " +
                          pp + " matchrad=" + to_string(match_tol.back()) + " >/dev/null 2>&1";

                cout << "  Run match for " + *it_file + " ... ";

                int ret = run_external(cmd_str);
                if ( ret ) {
                    cout << "Failed!\n";
                    cerr << "Something wrong while run application 'match'!\n";
                    throw (int)ROTCEN_ERROR_APP_FAILED;
                }

                cout << "OK!\n";

                // read result matched catalogs (only ID(NUMBER) columns)

                ret = read_catalog(matchedA, 1, current_cat);
                if ( ret != ROTCEN_ERROR_OK ) {
                    cerr << "Something wrong while reading matched.mtA file!\n";
                    throw ret;
                }
                // REARRANGE HERE!!!
                obj_id[0] = current_cat[0];

                ret = read_catalog(matchedB, 1, current_cat);
                if ( ret != ROTCEN_ERROR_OK ) {
                    cerr << "Something wrong while reading matched.mtB file!\n";
                    throw ret;
                }
                obj_id[i_cat] = current_cat[0];

            }

        } else { // use of astrometrical solution

        }

    } catch (int err) {
        input_list_file.close();
        if ( use_match && !dont_delete ) { // delete temporary files
            boost::filesystem::remove(ROTCEN_SEX_PARAM_FILE);
            for (auto filename = sex_cats.begin(); filename != sex_cats.end(); ++filename) {
                boost::filesystem::remove(*filename);
            }
        }
        return err;
    }

    // delete temporary files
    if ( use_match && !dont_delete ) {
        boost::filesystem::remove(ROTCEN_SEX_PARAM_FILE);
        for (auto filename = sex_cats.begin(); filename != sex_cats.end(); ++filename) {
            boost::filesystem::remove(*filename);
        }
    }

    return ROTCEN_ERROR_OK;
}
