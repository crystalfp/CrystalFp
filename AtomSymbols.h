
#ifndef ATOM_SYMBOLS_H
#define ATOM_SYMBOLS_H

#include <vector>

extern const char *elementZToSymbol(unsigned int aAtomZ);
extern unsigned int elementSymbolToZ(const char* aElementSymbol, unsigned int aDefaultX=1);
extern void convertAtomSymbolsToZ(const char *aAtomTypesStr, std::vector<unsigned int>& aAtomZ);

#endif

