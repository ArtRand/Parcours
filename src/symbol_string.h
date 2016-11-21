#ifndef PARCOURS_SYMBOL_STRING_H
#define PARCOURS_SYMBOL_STRING_H

#include "stl_includes.h"

typedef enum {
    a=0,
    c=1,
    g=2,
    t=3,
    n=4,
} Symbol;

typedef std::vector<Symbol> SymbolString;

SymbolString SymbolStringFromString(std::string& seq);
Symbol CharToSymbol(char b);
std::string StringFromSymbolString(SymbolString seq);
char SymbolToChar(Symbol s);

#endif // PARCOURS_SYMBOL_STRING_H
