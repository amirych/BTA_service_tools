#include "ascii_file.h"

#include <cstring>
#include <errno.h>
#include <stdarg.h>
#include <cstdlib>

AsciiFile::AsciiFile(const char comment_symbol): ifstream(), CommentSymbol(comment_symbol)
{
}


AsciiFile::AsciiFile(const char* filename, const char comment_symbol):
    ifstream(filename), CommentSymbol(comment_symbol)
{
//    open(filename);
}


//void AsciiFile::Open(const char *filename)
//{
//    fst.open(filename);
//    if ( !fst ) {
//        throw bad_ascii_file(AsciiFile::OpenFailure);
//    }
//}


AsciiFile::AsciiFileFlag AsciiFile::ReadLine(size_t N_elems, ...)
{
    char s[MAX_STR_LEN];
    char *subs, *start, *end;
    double val,*pval;
    size_t i;


//    fst.getline(s,MAX_STR_LEN);
//    if ( fst.eof() ) return AsciiFile::Eof;
    getline(s,MAX_STR_LEN);
    if ( eof() ) return AsciiFile::Eof;

    // delete leading spaces
    for (i = 0; i < MAX_STR_LEN; ++i ) {
        if ( s[i] != ' ' ) break;
    }
    if ( i == strlen(s) ) return AsciiFile::EmptyString;
//    subs = strndup(s+i,MAX_STR_LEN-i);
    subs = strdup(s+i);
    if ( errno == ENOMEM ) throw bad_alloc();

    if ( subs[0] == CommentSymbol ) return AsciiFile::CommentString;

    va_list args;
    va_start(args, N_elems);

    start = subs;
    for (i = 0; i < N_elems; ++i ) {
        val = strtod(start,&end);
        if ( end != start ) {
            pval = va_arg(args,double*);
            *pval = val;
        } else break;
        start = end;
    }
    free(subs);
    va_end(args);

    if ( i != N_elems ) {
        return AsciiFile::InvalidData;
    }

    return AsciiFile::DataString;
}


AsciiFile::AsciiFileFlag AsciiFile::ReadLine(size_t N_elems, vector<double> *data)
{
    char s[MAX_STR_LEN];
    char *subs, *start, *end;
    double val,*pval;
    size_t i;


    getline(s,MAX_STR_LEN);
    if ( eof() ) return AsciiFile::Eof;

    // delete leading spaces
    for (i = 0; i < MAX_STR_LEN; ++i ) {
        if ( s[i] != ' ' ) break;
    }
    if ( i == strlen(s) ) return AsciiFile::EmptyString;
    subs = strdup(s+i);
//    subs = strndup(s+i,MAX_STR_LEN-i);
    if ( errno == ENOMEM ) throw bad_alloc();

    if ( subs[0] == CommentSymbol ) return AsciiFile::CommentString;

    data->clear();
    start = subs;
    for (i = 0; i < N_elems; ++i ) {
        val = strtod(start,&end);
        if ( end != start ) {
            data->push_back(val);
        } else break;
        start = end;
    }
    free(subs);

    if ( i != N_elems ) {
        return AsciiFile::InvalidData;
    }

    return AsciiFile::DataString;
}




//void AsciiFile::Close()
//{
//    fst.close();
//}
