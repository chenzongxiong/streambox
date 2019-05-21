#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

using namespace std;

int main(int argc, char *argv[])
{
	unordered_map<string, int> map;
	
	string s1("hello");
	string s2("world");
	string s3("hello");
	
	map[s1] = 1;
	map[s2] = 2;
	map[s3] = 3;
		
	for (auto && kv : map) {
		cout << kv.first << " "  << kv.second << "\n";
	}
	return 0;
}
