#include<iostream>
#include<fstream>
#include <iomanip>
//#include<libgen.h>
//#include<cstdlib>
#include<cstdio>
#include<list>
#include<regex>
#include<ctime>

#define BOOST_NO_CXX11_SCOPED_ENUMS // special definition to fix Boost's copy_file and -std=c++11 linking error
#include<boost/program_options.hpp>
#include<boost/filesystem.hpp>
#include<boost/algorithm/string.hpp>

#include<fitsio.h>
#include<gsl/gsl_linalg.h>

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
#define ROTCEN_ERROR_EMPTY_CAT 110
#define ROTCEN_ERROR_BAD_ALLOC 120
#define ROTCEN_ERROR_BAD_MATCH 130
#define ROTCEN_ERROR_CANNOT_SOLVE 140
#define ROTCEN_ERROR_CANNOT_CREATE_RESULT_FILE 150

#define ROTCEN_ERROR_CFITSIO 1000 // displacement for CFITSIO error code

static string ROTCEN_SEX_PARAM_FILE = "sex.param";
static string ROTCEN_MATCH_REF_CAT = "ref.cat";

static string ROTCEN_SEX_EXE = "sex";
static string ROTCEN_AST_EXE = "solve-field";
static string ROTCEN_MATCH_EXE = "match";



/*
    The function tries to execute external application given by 'cmd_str' string.
    It checks exit code of the application.
*/
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


/*
    The routine reads data from FITS binary table. It assumes the binary table format
    is according to RDLS-files of 'solve-field' application
*/
static int read_fits_catalog(string &filename, vector<vector<double> > &data)
{
    int fits_status = 0;
    fitsfile *file;
    int ret_code = ROTCEN_ERROR_OK;
    long opt_nrows,nrows;
    double *ra = nullptr;
    double *dec = nullptr;

    try {
        fits_open_table(&file,filename.c_str(),READONLY,&fits_status);
        if ( fits_status ) throw fits_status;

        fits_get_num_rows(file, &nrows, &fits_status);
        if ( fits_status ) throw fits_status;

        fits_get_rowsize(file,&opt_nrows,&fits_status);
        if ( fits_status ) throw fits_status;

        if ( opt_nrows > nrows ) opt_nrows = nrows;

        ra = new double[opt_nrows];
        dec = new double[opt_nrows];

        data = vector<vector<double> >(3);

        long n_chunk = nrows/opt_nrows;
        long N_rest = nrows - n_chunk*opt_nrows;


        for ( long i = 0; i < n_chunk; ++i ) {
            fits_read_col(file,TDOUBLE,1,i*opt_nrows+1,1,opt_nrows,NULL,(void*)ra,NULL,&fits_status);
            if ( fits_status ) throw fits_status;

            fits_read_col(file,TDOUBLE,2,i*opt_nrows+1,1,opt_nrows,NULL,(void*)dec,NULL,&fits_status);
            if ( fits_status ) throw fits_status;

            auto it = data[1].end();
            data[1].insert(it,ra,ra+opt_nrows);

            it = data[2].end();
            data[2].insert(it,dec,dec+opt_nrows);
        }
        if ( N_rest ) {
            fits_read_col(file,TDOUBLE,1,n_chunk*opt_nrows+1,1,N_rest,NULL,(void*)ra,NULL,&fits_status);
            if ( fits_status ) throw fits_status;

            fits_read_col(file,TDOUBLE,2,n_chunk*opt_nrows+1,1,N_rest,NULL,(void*)dec,NULL,&fits_status);
            if ( fits_status ) throw fits_status;

            auto it = data[0].end();
            data[0].insert(it,ra,ra+N_rest);

            it = data[1].end();
            data[1].insert(it,dec,dec+N_rest);
        }
        // generate IDs column (just from 1 to size(RAcol))
        for ( double i = 1; i <= data[1].size(); ++i ) data[0].push_back(i);
    } catch (int err) {
        ret_code =  err + ROTCEN_ERROR_CFITSIO;
    } catch (bad_alloc &ex) {
        ret_code = ROTCEN_ERROR_BAD_ALLOC;
    }

    delete[] ra;
    delete[] dec;

    fits_close_file(file,&fits_status);

    return ret_code;
}

/*
    The function rearranges table of object IDs according to new vector of IDs for the first column.
    The table rows will be permutted according to new order of ID numbers in the new_id vector.
    NOTE: the algorithm assumes:
               1) new_id contains of unique numbers
               2) all numbers from new_id are members of table[0] vector
*/
static void rearrange_table(vector<vector<double> > &table, size_t last_col, vector<double> &new_id)
{
    size_t j;
    double tmp;

    for ( size_t i = 0; i < new_id.size(); ++i ) {
        j = i;
        while ( table[0].at(j) != new_id[i]) ++j;
        for ( size_t k = 0; k <= last_col; ++k ) { // swap elements
            tmp = table[k].at(i);
            table[k].at(i) = table[k].at(j);
            table[k].at(j) = tmp;
        }
    }
    for ( size_t k = 0; k <= last_col; ++k ) { // resize table columns
        table[k].resize(new_id.size());
    }
}


/*
    The function convert string sexagesimal presentation of the celestial coordinates
    to degrees value. The sexagesimal presentation is "[+/-]dd:mm:ss[.sss]" for DEC and
    "hh::mm:ss[.sss]" for RA with possible space characters before and after the string.

    The function does not check validness of the input string!
*/
static float sex2deg(string str, bool hours = false)
{
    string s = str;
    float val = 0.0;
    float p = 1.0;
    size_t idx;

    for ( int i = 0; i < 3; ++i, p *= 60.0 ) {
        val += stof(s,&idx)/p;
        if (idx < s.length() ) s = s.substr(idx+1);
    }

    if ( hours ) val *= 15.0;

    return val;
}


int main(int argc, char* argv[])
{

    // some defaults

    vector<float> sex_thresh(1,5.0); // sextractor's THRESH default value
    vector<float> match_tol;        // default radius of coordinate matching. It is in arcsecs for astrometrical solution or
                                    // dimensionless value for 'match application' solution ('matchrad' parameter)
                                    // The actuall default value is set below according to '--use-match'

    vector<string> sex_pars = {"-DETECT_TYPE CCD -SATUR_LEVEL 40000 -CHECKIMAGE_TYPE NONE -FILTER N -GAIN 1.0 -VERBOSE_TYPE QUIET"};

//    vector<string> solve_field_pars = {"--no-plots -M none --no-fits2fits -O -B none -W none -U none -y -L 0.1 -H 0.7 -u arcsecperpix"};
    vector<string> solve_field_pars = {"--no-plots -M none --no-fits2fits -O -B none -W none -y -L 0.1 -H 0.7 -u arcsecperpix"};

    vector<string> match_pars = {"id1=0 id2=0 min_scale=0.9 max_scale=1.1 linear"}; // default 'match' commandline parameters

    vector<string> sex_cat_prefix = {"obj_"}; // SExtractor output catalog  filename prefix

    vector<string> ast_prefix = {"wcs_"}; // astrometry-calibrated filename prefix

    vector<float> ra_deg, dec_deg; // guess value for RA and DEC for astrometrical solution
    ra_deg.push_back(0.0);
    dec_deg.push_back(0.0);

    vector<float> ast_radius = {0.5}; // search radius aroung guess RA and DEC for astrometry solution
    vector<string> solve_field_config = {"/usr/local/astrometry/etc/astrometry.cfg"};

    vector<string> ra_keyword = {"RA"};
    vector<string> dec_keyword = {"DEC"};

    string input_list_filename;
    string result_file;

    int ret_status = ROTCEN_ERROR_OK;

    // commandline options and arguments definitions

    po::options_description visible_opts("Allowed options");
    visible_opts.add_options()
        ("help,h", "produce help message")
        ("threshold,t", po::value<vector<float> >(&sex_thresh), "set threshold level for object detection (sextractor's DETECT_THRESH keyword)")
        ("radius,r", po::value<vector<float> >(&match_tol), "radius of coordinate matching [arcsecs for astrometrical solution]")
        ("use-match,m","use of 'match' application instead of astrometry (explicitly set '-s' option)")
        ("use-sex,s","use of sextractor to detect objects (in case of astrometrical solution)")
        ("sex-pars",po::value<vector<string> >(), "sextractor's parameters")
        ("solve-field-pars",po::value<vector<string> >(), "'solve-field' parameters")
        ("match-pars",po::value<vector<string> >(), "'match' parameters")
        ("dont-delete,d","do not delete temporary files")
        ("solve-field-config,c",po::value<vector<string> >(), "filename with full path of 'solve-field' config")
        ("ra",po::value<vector<float> >(), "Guess RA for the field (in degrees)")
        ("dec",po::value<vector<float> >(), "Guess DEC for the field (in degrees)")
        ("ra-key",po::value<vector<string> >(), "FITS-keyword name with RA guess value")
        ("dec-key",po::value<vector<string> >(), "FITS-keyword name with DEC guess value")
        ("ra-in-hours","RA value in FITS-keyword is given in hours")
        ("ra-dec-str","RA and DEC values in FITS-keywords are given in form of sexagesimal string (RA: hh:mm:ss.ss, DEC: dd:mm:ss.ss)")
        ("search-radius",po::value<vector<float> >(), "search radius for astrometrical solution (in degrees)")
        ("save-wcs", "Save WCS-calibrated FITS-files (assuming 'solve-field' type of matching)");


    po::options_description hidden_opts("");
    hidden_opts.add_options()
            ("input-file", po::value<string>()->required())
            ("result-file", po::value<string>());

    po::options_description cmd_opts("");
    cmd_opts.add(visible_opts).add(hidden_opts);

    po::positional_options_description pos_arg;
    pos_arg.add("input-file", 1);
    pos_arg.add("result-file", 2);

    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv).options(cmd_opts).positional(pos_arg).run(), vm);

        if ( vm.count("help") ) {
            string head_str = "Usage: " + boost::filesystem::basename(argv[0]);
            string skip_str(head_str.length()+1,' ');

            cout << head_str << " [-h] [-t num] [-r num] [-d] [--solve-field-pars]\n" << skip_str <<
                                "[--use-match] [--match-pars str] \n" << skip_str <<
                                "[--use-sex] [--sex-pars str]\n" << skip_str <<
                                "[--ra num] [--deg num] [--search-radius num]\n" << skip_str <<
                                "[--ra-key str] [--dec-key str] [--ra-in-hours] [--ra-dec-str]\n" << skip_str <<
                                "[--solve-field-pars str] [--solve-field-config str] [--save-wcs] input_list [result_file]\n\n";

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
        solve_field_pars.push_back(vm["solve-field-pars"].as<vector<string> >().back());
    }

    if ( vm.count("solve-field-config") ) {
        solve_field_config.erase(solve_field_config.begin(),solve_field_config.end());
        solve_field_config.push_back(vm["solve-field-config"].as<vector<string> >().back());
    }

    if ( vm.count("match-pars") ) {
        match_pars.erase(match_pars.begin(),match_pars.end());
        match_pars.push_back(vm["match-pars"].as<vector<string> >().back());
    }

    if (vm.count("input-file")) {
        input_list_filename = vm["input-file"].as<string>();
    }

    if ( vm.count("result-file") ) {
        result_file = vm["result-file"].as<string>();
    }


    bool use_match = false;
    bool use_sex = false;
    bool use_guess_radec = false;
    bool save_wcs = false;

    if ( vm.count("use-match") ) {
        int ret = system("match --help  >/dev/null 2>&1"); // try to run command 'match'
        int exit_code = WEXITSTATUS(ret);
        if ( ret == -1 || exit_code == 127 ) {
            cerr << "Application 'match' is not available!\n";
            return ROTCEN_ERROR_UNAVAILABLE_CMD;
        }

        use_match = true;

        if ( !vm.count("radius") ) { // use default value
            match_tol = {1.0};
        }


        sex_pars.back() += " -PARAMETERS_NAME " + ROTCEN_SEX_PARAM_FILE + " -CATALOG_TYPE ASCII" +
                           " -DETECT_THRESH " + to_string(sex_thresh.back()) +
                           " -ANALYSIS_THRESH " + to_string(sex_thresh.back());


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

        if ( vm.count("solve-field-config") ) {
            solve_field_config.back() = vm["solve-field-config"].as<vector<float> >().back();
        }

//        solve_field_pars.back() += " --config " + solve_field_config.back();
        solve_field_pars.back() += " -b " + solve_field_config.back();

        if ( vm.count("ra") && vm.count("dec") ) { // it makes sense only if the both are given
            use_guess_radec = true;

            ra_deg = vm["ra"].as<vector<float> >();
            dec_deg = vm["dec"].as<vector<float> >();

            solve_field_pars.back() += " --ra " + to_string(ra_deg.back()) + " --dec " + to_string(dec_deg.back());
        }

        if ( vm.count("search-radius") ) {
            ast_radius = vm["search-radius"].as<vector<float> >();
//            solve_field_pars.back() += " --radius " + to_string(ast_radius.back());
        }

        solve_field_pars.back() += " --radius " + to_string(ast_radius.back());

        if ( vm.count("ra-key") ) {
            ra_keyword = vm["ra-key"].as<vector<string> >();
        }

        if ( vm.count("dec-key") ) {
            dec_keyword = vm["dec-key"].as<vector<string> >();
        }

        if ( !vm.count("radius") ) { // use default value
            match_tol = {0.3};
        }

        if ( vm.count("save-wcs") ) {
            save_wcs = true;
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

        solve_field_pars.back() += " --use-sextractor";
        sex_pars.back() += " -DETECT_THRESH " + to_string(sex_thresh.back()) +
                           " -ANALYSIS_THRESH " + to_string(sex_thresh.back());
        solve_field_pars.back() += " --sextractor-path \"" + ROTCEN_SEX_EXE + " " + sex_pars.back() + "\" ";

        if ( !vm.count("radius") ) { // use default value
            match_tol = {0.5};
        }
    } else {
        solve_field_pars.back() += " --sigma " + to_string(sex_thresh.back());
    }


    list<string> input_files;
    list<string> sex_cats, ast_cat;
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
            boost::algorithm::trim(str);
            if ( str[0] == '#' ) continue; // comment
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

        cout << "\nObjects detection:\n";

        for ( auto it_file = input_files.begin(); it_file != input_files.end(); ++it_file ) {
            boost::filesystem::path pp = *it_file;
            string path = pp.parent_path().string();
            if ( path.empty() ) path = ".";
            string file = boost::filesystem::basename(*it_file);

            if ( use_match ) { // skip astrometry, just detect objects using sextractor

                file = path + boost::filesystem::path::preferred_separator + sex_cat_prefix.back() + file + ".cat";

                string cmd_str = ROTCEN_SEX_EXE + " " + sex_pars.back() + " -CATALOG_NAME " +
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
//                file = path + boost::filesystem::path::preferred_separator + ast_prefix.back() + file + ".fits";

//                string cmd_str = ROTCEN_AST_EXE + " " + solve_field_pars.back() +
//                                 " -R " + file;

                string solved_file = path + boost::filesystem::path::preferred_separator + file + ".solved";
                string rdls_file = path + boost::filesystem::path::preferred_separator + file + ".rdls";

//                cmd_str += " -S " + solved_file + " -N none";

                string cmd_str = ROTCEN_AST_EXE + " " + solve_field_pars.back();

                if ( save_wcs ) { // save WCS-calibrated FITS-file
                    cmd_str += " -N " + path + boost::filesystem::path::preferred_separator + ast_prefix.back() + file + ".fits";
                } else {
                    cmd_str += " -N none";
                }

                if ( !use_guess_radec ) { // no user's guess RA and DEC in commandline
                    int fits_status = 0;      // try to read RA and DEC from FITS header
                    fitsfile* file;

                    char key_value[81];

                    fits_open_image(&file,(*it_file).c_str(),READONLY,&fits_status);
                    if ( fits_status ) {
                        cerr << "Something wrong while opening " << *it_file << " file!\n";
                        throw ROTCEN_ERROR_CFITSIO + fits_status;
                    }

                    fits_read_keyword(file,ra_keyword.back().c_str(),key_value,NULL,&fits_status);
                    if ( fits_status == KEY_NO_EXIST ) goto jmp; // No RA keyword. Just skip
                    if ( fits_status ) {
                        cerr << "Something wrong while reading '" << ra_keyword.back() << "' keyword in " << *it_file << " file!\n";
                        throw ROTCEN_ERROR_CFITSIO + fits_status;
                    }
                    if ( vm.count("ra-dec-str") ) { // RA should be in sexagesimal form
                        if (regex_match(key_value,regex("^ *'.+' *$"))) { // for SCORPIO format of the FITS keyword
                            string ss = key_value;
                            boost::algorithm::trim(ss);
                            size_t len = ss.copy(key_value,ss.length()-2,1);
                            key_value[len] = '\0';
                        }
                        if (!regex_match(key_value,regex("^ *\\+?\\d\\d:\\d\\d:\\d\\d(\\.{1}\\d*)? *$")) ) { // is it in sexagesimal?
                            cerr << "Invalid RA value in " << ra_keyword.back() << " FITS-keyword of " << *it_file << " file!\n";
                            fits_close_file(file,&fits_status);
                            throw (int)ROTCEN_ERROR_BAD_DATA;
                        }
                        // convert to degrees
                        ra_deg.back() = sex2deg(key_value,true);
                    } else { // RA should be non-negative numeric values
                        if (!regex_match(key_value,regex("^ *\\+?\\d+(\\.{1}\\d*)?([EeDd][+-]?\\d+)? *$")) ) { // is it a number?
                            cerr << "Invalid RA value in " << ra_keyword.back() << " FITS-keyword of " << *it_file << " file!\n";
                            fits_close_file(file,&fits_status);
                            throw (int)ROTCEN_ERROR_BAD_DATA;
                        }

                        ra_deg.back() = stof(key_value);

                        if ( vm.count("ra-in-hours") ) {
                            ra_deg.back() *= 15.0;
                        }
                    }

                    fits_read_keyword(file,dec_keyword.back().c_str(),key_value,NULL,&fits_status);
                    if ( fits_status == KEY_NO_EXIST ) goto jmp; // No DEC keyword. just skip
                    if ( fits_status ) {
                        cerr << "Something wrong while reading '" << ra_keyword.back() << "' keyword in " << *it_file << " file!\n";
                        throw ROTCEN_ERROR_CFITSIO + fits_status;
                    }
                    if ( vm.count("ra-dec-str") ) { // DEC should be in sexagesimal form
                        if (regex_match(key_value,regex("^ *'.+' *$"))) { // for SCORPIO format of the FITS keyword
                            string ss = key_value;
                            boost::algorithm::trim(ss);
                            size_t len = ss.copy(key_value,ss.length()-2,1);
                            key_value[len] = '\0';
                        }
                        if (!regex_match(key_value,regex("^ *[+-]?\\d\\d:\\d\\d:\\d\\d(\\.{1}\\d*)? *$")) ) { // is it in sexagesimal?
                            cerr << "Invalid DEC value in " << dec_keyword.back() << " FITS-keyword of " << *it_file << " file!\n";
                            fits_close_file(file,&fits_status);
                            throw (int)ROTCEN_ERROR_BAD_DATA;
                        }
                        // convert to degrees
                        dec_deg.back() = sex2deg(key_value);
                    } else { // DEC should be non-negative numeric values
                        if (!regex_match(key_value,regex("^ *[+-]?\\d+(\\.{1}\\d*)?([EeDd][+-]?\\d+)? *$")) ) { // is it a number?
                            cerr << "Invalid DEC value in " << dec_keyword.back() << " FITS-keyword of " << *it_file << " file!\n";
                            fits_close_file(file,&fits_status);
                            throw (int)ROTCEN_ERROR_BAD_DATA;
                        }
                        dec_deg.back() = stof(key_value);
                    }

                    cmd_str += " --ra " + to_string(ra_deg.back()) + " --dec " + to_string(dec_deg.back());

                    fits_close_file(file,&fits_status);
                }
                jmp:

                cmd_str +=  " " + *it_file + "  >/dev/null 2>&1";
//                cmd_str +=  " " + *it_file

                cout << "  Run solve-field for " + *it_file + " ... ";

//                cout << cmd_str << endl;

                int ret = run_external(cmd_str); // try to run command 'solve-field'
                bool ok = boost::filesystem::exists(solved_file);
                if ( ret || !ok ) {
                    cerr << "ret=" << ret << endl;
                    cout << "Failed!\n";
                    cerr << "Something wrong while run application 'solve-field'!\n";
                    throw (int)ROTCEN_ERROR_APP_FAILED;
                }

                cout << "OK!\n";

//                ast_cat.push_back(file);
                ast_cat.push_back(rdls_file);
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

            string matchedA = "matched.mtA";
            string matchedB = "matched.mtB";

            for ( long i_cat = 1; it_file != sex_cats.end(); ++it_file, ++i_cat ) {

                // read current catalog

                ret = read_catalog(*it_file, 3, current_cat);
                if ( ret != ROTCEN_ERROR_OK ) {
                    cerr << "Something wrong while reading " << *it_file << " file!\n";
                    throw ret;
                }
                if ( current_cat[0].empty() ) {
                    cerr << "Empty catalog in file " << *it_file << " file!\n";
                    throw (int)ROTCEN_ERROR_EMPTY_CAT;
                }
                obj_cat[i_cat*3] = current_cat[0]; // NUMBER
                obj_cat[i_cat*3+1] = current_cat[1]; // X_IMAGE
                obj_cat[i_cat*3+2] = current_cat[2]; // Y_IMAGE

                cmd_str = ROTCEN_MATCH_EXE + " " + ROTCEN_MATCH_REF_CAT + " 1 2 3 " + *it_file + " 1 2 3 " +
                match_pars.back() + " matchrad=" + to_string(match_tol.back()) + " >/dev/null 2>&1";

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
                if ( current_cat[0].empty() ) {
                    cerr << "Empty catalog in file " << matchedA << " file!\n";
                    throw (int)ROTCEN_ERROR_EMPTY_CAT;
                }

                cout << "    Matched " << current_cat[0].size() << " objects\n";

                if ( i_cat > 1 ) rearrange_table(obj_id,i_cat-1,current_cat[0]); else obj_id[0] = current_cat[0];

                ret = read_catalog(matchedB, 1, current_cat);
                if ( ret != ROTCEN_ERROR_OK ) {
                    cerr << "Something wrong while reading matched.mtB file!\n";
                    throw ret;
                }
                if ( current_cat[0].empty() ) {
                    cerr << "Empty catalog in file " << matchedB << " file!\n";
                    throw (int)ROTCEN_ERROR_EMPTY_CAT;
                }
                obj_id[i_cat] = current_cat[0];

                boost::filesystem::copy_file(matchedA,ROTCEN_MATCH_REF_CAT,boost::filesystem::copy_option::overwrite_if_exists);
            }

//            cout << endl << endl;
//            for (size_t k = 0; k < input_files.size(); ++k ) {
//                cout << obj_cat[3*k+1].at(obj_id[k].at(0)-1) << ", " << obj_cat[3*k+2].at(obj_id[k].at(0)-1) << ", ";
//            }
//            cout << endl << endl;

        } else { // use of astrometrical solution
            cout << "\nMatching objects using astrometrical solution:\n";

            auto it_file = ast_cat.begin();
            ++it_file; // point to the second catalog

            // read the first catalog
            int ret = read_fits_catalog(ast_cat.front(),current_cat);
            if ( ret != ROTCEN_ERROR_OK ) {
                cerr << "Something is wrong while reading " << ast_cat.front() << " file!\n";
                throw ret;
            }
            obj_cat[0] = current_cat[0]; // ID
            obj_cat[1] = current_cat[1]; // RA values
            obj_cat[2] = current_cat[2]; // DEC values

            obj_id[0] = current_cat[0];

            double min_dist = match_tol.back()*match_tol.back()/3600.0/3600.0; // in square degrees
            double dist;

//            vector<double> dist(current_cat[0].size());

            for ( long i_cat = 1; it_file != ast_cat.end(); ++it_file, ++i_cat ) {
                // read current catalog
                ret = read_fits_catalog(*it_file,current_cat);
                if ( ret != ROTCEN_ERROR_OK ) {
                    cerr << "Something is wrong while reading " << *it_file << " file!\n";
                    throw ret;
                }

                long cat_col = 3*i_cat;

                obj_cat[cat_col] = current_cat[0]; // ID
                obj_cat[cat_col+1] = current_cat[1]; // RA
                obj_cat[cat_col+2] = current_cat[2]; // DEC

                // matching

                cout << "  0 <--> " << i_cat << ", ";
                long N_matched = 0;
                for ( long idx = 0; idx < obj_id[0].size(); ++idx ) {
//                    double min_d = min_dist;
                    long min_ind = 0;
                    for ( long j = 0; j < obj_cat[cat_col].size(); ++j ) {

                        /*
                            Actually, here one can use strict equality for RA and DEC because of
                            solve-field RDLS-files consist of coordinates from index-file. By other
                            words the RA and DEC in 'obj_cat' catalogs are not computed.
                        */


                        if ( (obj_cat[1].at(obj_id[0].at(idx)-1) == obj_cat[cat_col+1].at(j)) &&
                             (obj_cat[2].at(obj_id[0].at(idx)-1) == obj_cat[cat_col+2].at(j)) ) {

                            current_cat[0].at(N_matched) = obj_id[0].at(idx);
                            min_ind = obj_cat[cat_col].at(j);
                            break;
                        }

                        /*
                         * The closest points realization ...
                         *
                        // compute distance ("-1" in the index computation since ID starts from 1!)
                        double rd = obj_cat[1].at(obj_id[0].at(idx)-1) - obj_cat[cat_col+1].at(j);
                        double dd = obj_cat[2].at(obj_id[0].at(idx)-1) - obj_cat[cat_col+2].at(j);

                        dist = rd*rd + dd*dd;

                        if ( dist < min_d ) { // store ID of matched objects
                            if ( N_matched == current_cat[0].size() ) {
                                cerr << "Something wrong while matching objects!\n";
                                cerr << "It seems there was to big matching radius!\n";
                                throw ROTCEN_ERROR_BAD_MATCH;
                            }
                            current_cat[0].at(N_matched) = obj_id[0].at(idx);
                            min_ind = obj_cat[cat_col].at(j);
                            min_d = dist;     // search for the closest objects
                        }

                        */
                    }
                    if ( min_ind ) {
                        obj_id[i_cat].push_back(min_ind);
                        ++N_matched;
                    }
                }

                if ( !N_matched ) {
                    cerr << "No matching objects in the input catalogs!\n";
                    throw (int)ROTCEN_ERROR_EMPTY_CAT;
                }

                cout << N_matched << " objects were matched\n";
                if ( N_matched != current_cat[0].size() ) current_cat[0].resize(N_matched);

                if ( i_cat > 1 ) {
                    rearrange_table(obj_id,i_cat-1,current_cat[0]);
                }
            }

            // read catalogs with pixel coordinates
            size_t i_cat = 0;
            for ( it_file = ast_cat.begin(); it_file != ast_cat.end(); ++it_file, ++i_cat ) {
                boost::filesystem::path pp = *it_file;
                string path = pp.parent_path().string();
                if ( path.empty() ) path = ".";
                string file = boost::filesystem::basename(*it_file);

                file = path + boost::filesystem::path::preferred_separator + file + "-indx.xyls";

                int ret = read_fits_catalog(file,current_cat);
                if ( ret != ROTCEN_ERROR_OK ) {
                    cerr << "Something is wrong while reading " << file << " file!\n";
                    throw ret;
                }

                obj_cat[i_cat*3] = current_cat[0]; // ID
                obj_cat[i_cat*3+1] = current_cat[1]; // X values
                obj_cat[i_cat*3+2] = current_cat[2]; // Y values
            }


//            cout << endl << endl;
//            for ( size_t j = 0; j < obj_id[0].size(); ++j ) {
//                for (size_t k = 0; k < input_files.size(); ++k ) {
////                    cout << obj_cat[3*k+1].at(obj_id[k].at(0)-1) << ", " << obj_cat[3*k+2].at(obj_id[k].at(0)-1) << ", ";
//                    cout << obj_id[k].at(j) << " ";
//                }
//                cout << endl;
//            }
//            cout << endl << endl;

//            for (size_t k = 0; k < input_files.size(); ++k ) {
//                cout << obj_cat[3*k+1].at(obj_id[k].at(0)-1) << ", " << obj_cat[3*k+2].at(obj_id[k].at(0)-1) << ", ";
//            }
//            cout << endl << endl;
        }

        // compute rotation center


        cout << "\nSolving ... ";

        gsl_matrix* sys_mat = NULL;
        gsl_vector *b = NULL;
        gsl_vector *x = NULL;
        gsl_vector *tau = NULL;
        gsl_vector *res = NULL;

        size_t N_circles = obj_id[0].size();
        size_t N_objs = obj_id.size();

        size_t N_eq = N_circles*(N_objs-1)*N_objs/2; // number of linear equations

        try {
            gsl_error_handler_t *gsl_eh =  gsl_set_error_handler_off(); // turn off GSL default error handler

            sys_mat = gsl_matrix_alloc(N_eq,2);
            if ( sys_mat == NULL ) {
                cerr << "Cannot allocate memory for system matrix!\n";
                throw (int)ROTCEN_ERROR_BAD_ALLOC;
            }

            b = gsl_vector_alloc(N_eq);
            if ( b == NULL ) {
                cerr << "Cannot allocate memory for system right-hand part!\n";
                throw (int)ROTCEN_ERROR_BAD_ALLOC;
            }

            x = gsl_vector_alloc(2);
            if ( x == NULL ) {
                cerr << "Cannot allocate memory for system solution vector!\n";
                throw (int)ROTCEN_ERROR_BAD_ALLOC;
            }

            tau = gsl_vector_alloc(2);
            if ( tau == NULL ) {
                cerr << "Cannot allocate memory for QR decomposition tau vector!\n";
                throw (int)ROTCEN_ERROR_BAD_ALLOC;
            }

            res = gsl_vector_alloc(N_eq);
            if ( res == NULL ) {
                cerr << "Cannot allocate memory for system residual vector!\n";
                throw (int)ROTCEN_ERROR_BAD_ALLOC;
            }


            // fill system matrix and right-hand part:
            size_t i = 0;
            for ( size_t i_circ = 0; i_circ < N_circles; ++i_circ ) {
                for ( size_t i_obj = 0; i_obj < (N_objs-1); ++i_obj ) {
                    size_t i_col1 = i_obj*3;
                    // "-1" in the index computation since ID starts from 1!
                    double x1 = obj_cat[i_col1 + 1].at(obj_id[i_obj].at(i_circ)-1);
                    double y1 = obj_cat[i_col1 + 2].at(obj_id[i_obj].at(i_circ)-1);

                    for ( size_t j = i_obj+1; j < N_objs; ++j ) {
                        size_t i_col2 = j*3;

                        // "-1" in the index computation since ID starts from 1!
                        double x2 = obj_cat[i_col2 + 1].at(obj_id[j].at(i_circ)-1);
                        double y2 = obj_cat[i_col2 + 2].at(obj_id[j].at(i_circ)-1);


                        gsl_matrix_set(sys_mat,i,0,2.0*(x2-x1));
                        gsl_matrix_set(sys_mat,i,1,2.0*(y2-y1));

                        gsl_vector_set(b,i,x2*x2+y2*y2-x1*x1-y1*y1);
                        ++i;
                    }

                }
            }


            int ret = gsl_linalg_QR_decomp(sys_mat,tau);
            if ( ret ) {
                cerr << "Something wrong while QR decomposition!\n";
                throw ROTCEN_ERROR_CANNOT_SOLVE;
            }

            ret = gsl_linalg_QR_lssolve(sys_mat,tau,b,x,res);
            if ( ret ) {
                cerr << "Something wrong while system solving!\n";
                throw ROTCEN_ERROR_CANNOT_SOLVE;
            }

            // compute residual
            double residual = 0.0;
            for ( size_t i = 0; i < N_eq; ++i ) {
                residual += gsl_vector_get(res,i)*gsl_vector_get(res,i);
            }

            cout << "OK!\n\n";
            cout << "Solution: " << endl;
            cout << "  rotation center: [" << gsl_vector_get(x,0) << ", " << gsl_vector_get(x,1) << "]" <<
                    " (residual: " << sqrt(residual)/(N_eq-1) << ")\n";

            // save result file if given

            if ( !result_file.empty() ) {
                ofstream rfile(result_file);
                if ( !rfile.good() ) {
                    cerr << "Cannot open result file!\n";
                    throw (int)ROTCEN_ERROR_CANNOT_CREATE_RESULT_FILE;
                }


                rfile << "# \n";
                rfile << "# Computation of rotation center ('" << boost::filesystem::basename(argv[0]) << "' application, ";

                time_t t = time(nullptr);
                string tt = asctime(localtime(&t));
                rfile <<  tt.substr(0,tt.length()-1) << ")\n"; // delete trailng '\n' from asctime-string

                rfile << "# \n";
                rfile << "# Input file: " << input_list_filename << "\n";
                rfile << "# Method: ";
                if ( use_match ) {
                    rfile << "match application (pixel coordinates matching using triangles)\n";
                } else {
                    rfile << "astrometrical solution (astrometry.net 'solve-field' application)\n";
                }
                rfile << "# \n";
                rfile << "# Number of points per circle: " << N_objs << endl;
                rfile << "# Number of circles: " << N_circles << endl;
                rfile << "# \n";
                rfile << "# Rotation center in pixel coordinates: \n";

                rfile << std::fixed << std::setprecision(1) << gsl_vector_get(x,0) << " " <<
                         std::fixed << std::setprecision(1) << gsl_vector_get(x,1) << endl;

                rfile.close();
            }

        } catch (int err) {
            gsl_matrix_free(sys_mat);
            gsl_vector_free(b);
            gsl_vector_free(x);
            gsl_vector_free(tau);
            gsl_vector_free(res);
            throw err;
        }

        gsl_matrix_free(sys_mat);
        gsl_vector_free(b);
        gsl_vector_free(x);
        gsl_vector_free(tau);
        gsl_vector_free(res);


    } catch (int err) {
        input_list_file.close();
        ret_status = err;
    }

    // delete temporary files
    if ( use_match ) {
        if ( !dont_delete ) {
            boost::filesystem::remove(ROTCEN_SEX_PARAM_FILE);
            for (auto filename = sex_cats.begin(); filename != sex_cats.end(); ++filename) {
                boost::filesystem::remove(*filename);
            }
            boost::filesystem::remove("matched.mtA");
            boost::filesystem::remove("matched.mtB");
            boost::filesystem::remove("matched.unA");
            boost::filesystem::remove("matched.unB");
            boost::filesystem::remove(ROTCEN_MATCH_REF_CAT);
        }
    } else {
        if ( !dont_delete ) {
            for (auto filename = ast_cat.begin(); filename != ast_cat.end(); ++filename) {
                boost::filesystem::remove(*filename);

                boost::filesystem::path pp = *filename;
                string path = pp.parent_path().string();
                string file = boost::filesystem::basename(*filename);

                boost::filesystem::remove(path + boost::filesystem::path::preferred_separator + file + "-indx.xyls");
                boost::filesystem::remove(path + boost::filesystem::path::preferred_separator + file + ".rdls");
                boost::filesystem::remove(path + boost::filesystem::path::preferred_separator + file + ".axy");
                boost::filesystem::remove(path + boost::filesystem::path::preferred_separator + file + ".solved");
            }
        }
    }


    return ret_status;
}
