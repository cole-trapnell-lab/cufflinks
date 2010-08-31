#ifndef TOKENIZE_H_
#define TOKENIZE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

void tokenize(const std::string& s,
              const std::string& delims,
              std::vector<std::string>& ss);

void tokenize_strict(const std::string& s,
                     const std::string& delims, 
                     std::vector<std::string>& ss);
#endif /*TOKENIZE_H_*/
