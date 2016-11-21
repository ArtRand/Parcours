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

char SymbolToChar(Symbol s) {
    switch (s) {
        case a:
            return 'A';
        case c:
            return 'C';
        case g:
            return 'G';
        case t:
            return 'T';
        case n:
            return 'N';
        default:
            throw ParcoursException("[SymbolStringFromString] illegal base %c\n", s);
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

std::string StringFromSymbolString(SymbolString seq) {
    std::string seq_string;
    seq_string.reserve(seq.size());
    for (Symbol s : seq) {
        seq_string += SymbolToChar(s);
    }
    return seq_string;
}
