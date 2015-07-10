#include<iostream>
#include<fstream>

#include<boost/program_options.hpp>

using namespace std;

namespace po = boost::program_options;

#define CMDOPT_


int main(int argc, char* argv[])
{
    if ( argc == 1 ) {
        cout << "Usage: rotation_center file_list\n";
        return 1;
    }

    vector<float> sex_thresh(1,5.0);

    // commandline options and arguments

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("threshold,t", po::value<vector<float> >(&sex_thresh), "set threshold level for object detection (sextractor THRESH keyword)")
        ("input-file",po::value<string>()->required());

    po::positional_options_description pos_arg;
    pos_arg.add("input-file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_arg).run(), vm);
//    po::notify(vm);

    // check files

    if ( vm.count("help") ) {
        cout << desc << "\n";
        return 1;
    }

    if ( vm.count("threshold") ) {
        sex_thresh = vm["threshold"].as<vector<float> >();
//        cout << "THRESH = " << vm["threshold"].as<float>() << ".\n";
        for ( int i = 0; i < sex_thresh.size(); ++i ) cout << "THRESH = " << sex_thresh[i] << "\n";
    } else {
        cout << "THRESH = " << sex_thresh.front() << ".\n";
    }

    if (vm.count("input-file"))
    {
        cout << "Input files are: "
             << vm["input-file"].as<string>() << "\n";
    }

//    if ( !ifstream().good() ) {
//        cerr << "Cannot find file of input frames list!\n";
//        return 10;
//    }
}
