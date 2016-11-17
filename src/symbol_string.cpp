#include "symbol_string.h"
#include "parcours_exceptions.h"


Symbol CharToSymbol(char b) {
    switch (b) {
        case 'A':
        case 'a':
            return a;
        case 'C':
        case 'c':
            return c;
        case 'G':
        case 'g':
            return g;
        case 'T':
        case 't':
            return t;
        case 'N':
        case 'n':
            return n;
        default:
            throw ParcoursException("[SymbolStringFromString] illegal base %c\n", b);
    }
}

SymbolString SymbolStringFromString(std::string& seq) {
    SymbolString S;

    auto base_to_symbol = [] (char base) -> Symbol {
        switch (base) {
            case 'A':
            case 'a':
                return a;
            case 'C':
            case 'c':
                return c;
            case 'G':
            case 'g':
                return g;
            case 'T':
            case 't':
                return t;
            case 'N':
            case 'n':
                return n;
            default:
                throw ParcoursException("[SymbolStringFromString] illegal base %c\n", base);
        }
    };

    for (char s : seq) {
        S.push_back(base_to_symbol(s)); 
    }

    return S;
}
