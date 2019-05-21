/*
 * re2.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: xzl
 *
 *  Purdue University, 2016
 *
 *  to build:
 *  g++ re2.cpp -o re2.bin -lre2 -std=c++11
 *
 *  ref:
 *  https://gist.github.com/chezou/1395527
 *
 */


#include <iostream>
#include <re2/re2.h>

using namespace std;
using namespace re2;

int main (int argc, char **argv)
{
    cout << "hello world" << endl;

    int matchResult;

    matchResult = RE2::FullMatch("hello", "h.*o");
    cout << "matchResult = " << matchResult << endl;

    string year;
    int yeari;

    /* match */
    matchResult = RE2::PartialMatch("this is a test 1983 hello",
//    		R"(\d\d\d\d)",  /* does not work */
    		"(\\d\\d\\d\\d)",
    		&year);
    cout << "matchResult = " << matchResult << endl;
    cout << year << endl;

    /* match */
    matchResult = RE2::PartialMatch("this is a test 1983 hello",
    		"(\\d\\d\\d\\d)",
    		&yeari);
    cout << "matchResult = " << matchResult << endl;
    cout << yeari << endl;

    RE2 re("(\\d\\d\\d\\d)");
    re2::StringPiece input("this is 1983 and 1984 and 1985");
    while(RE2::FindAndConsume(&input, re, &year)) {
    	cout << "matched: " << year << endl;
    }

    return 0;
}



