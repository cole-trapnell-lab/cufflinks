#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>
#include <cassert>

using namespace std;

/**
 * Split string s according to given delimiters.  Mostly borrowed
 * from C++ Programming HOWTO 7.3.
 */
void tokenize(const string& s, const string& delims, vector<string>& ss) {
	string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type pos = s.find_first_of(delims, lastPos);
	while (string::npos != pos || string::npos != lastPos) {
		ss.push_back(s.substr(lastPos, pos - lastPos));
        lastPos = s.find_first_not_of(delims, pos);
        pos = s.find_first_of(delims, lastPos);
    }
}

/**
 * Split string s according to given delimiters. If two delimiters occur in
 * succession, sticks an empty string token at that position in ss
 */
void tokenize_strict(const string& s, const string& delims, vector<string>& ss) {
	string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type pos = s.find_first_of(delims, lastPos);
	while (lastPos < s.length() || pos < s.length()) {
		ss.push_back(s.substr(lastPos, pos - lastPos));
		if (pos == string::npos)
			break;
        lastPos = pos + 1;
        pos = s.find_first_of(delims, lastPos);
    }
}
