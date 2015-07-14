#ifndef ASCII_FILE_H
#define ASCII_FILE_H

#include <iostream>
#include <fstream>
#include <vector>


using namespace std;


#define MAX_STR_LEN 512

//
// The class inherited from ifstream
// It adds the method ReadLine.
//
//
class AsciiFile: public ifstream
{
public:
    enum AsciiFileFlag {Eof,InvalidData,CommentString,EmptyString,DataString};

    AsciiFile(const char comment_symbol = '#');
    AsciiFile(const char* filename, const char comment_symbol = '#');

    // inputs:
    //   N_elems - number of elements to be read;
    //   outputs - N_elems pointers of double* types
    AsciiFileFlag ReadLine(size_t N_elems, ...);
    AsciiFileFlag ReadLine(size_t N_elems, vector<double> *data);
private:
    char CommentSymbol; // a char to be interpretated as comment symbol
                        // if inputs string starts by CommentSymbol it assumed to be a comment string
};

#endif // ASCII_FILE_H
