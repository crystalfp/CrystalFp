
#include <string>
#include <cctype>
#include "AtomSymbols.h"

static const unsigned int MAX_Z = 111;

static const char *gElementSymbols[MAX_Z+1] = {
/*   0 */        "Xx",
/*   1 */        "H",
/*   2 */        "He",
/*   3 */        "Li",
/*   4 */        "Be",
/*   5 */        "B",
/*   6 */        "C",
/*   7 */        "N",
/*   8 */        "O",
/*   9 */        "F",
/*  10 */        "Ne",
/*  11 */        "Na",
/*  12 */        "Mg",
/*  13 */        "Al",
/*  14 */        "Si",
/*  15 */        "P",
/*  16 */        "S",
/*  17 */        "Cl",
/*  18 */        "Ar",
/*  19 */        "K",
/*  20 */        "Ca",
/*  21 */        "Sc",
/*  22 */        "Ti",
/*  23 */        "V",
/*  24 */        "Cr",
/*  25 */        "Mn",
/*  26 */        "Fe",
/*  27 */        "Co",
/*  28 */        "Ni",
/*  29 */        "Cu",
/*  30 */        "Zn",
/*  31 */        "Ga",
/*  32 */        "Ge",
/*  33 */        "As",
/*  34 */        "Se",
/*  35 */        "Br",
/*  36 */        "Kr",
/*  37 */        "Rb",
/*  38 */        "Sr",
/*  39 */        "Y",
/*  40 */        "Zr",
/*  41 */        "Nb",
/*  42 */        "Mo",
/*  43 */        "Tc",
/*  44 */        "Ru",
/*  45 */        "Rh",
/*  46 */        "Pd",
/*  47 */        "Ag",
/*  48 */        "Cd",
/*  49 */        "In",
/*  50 */        "Sn",
/*  51 */        "Sb",
/*  52 */        "Te",
/*  53 */        "I",
/*  54 */        "Xe",
/*  55 */        "Cs",
/*  56 */        "Ba",
/*  57 */        "La",
/*  58 */        "Ce",
/*  59 */        "Pr",
/*  60 */        "Nd",
/*  61 */        "Pm",
/*  62 */        "Sm",
/*  63 */        "Eu",
/*  64 */        "Gd",
/*  65 */        "Tb",
/*  66 */        "Dy",
/*  67 */        "Ho",
/*  68 */        "Er",
/*  69 */        "Tm",
/*  70 */        "Yb",
/*  71 */        "Lu",
/*  72 */        "Hf",
/*  73 */        "Ta",
/*  74 */        "W",
/*  75 */        "Re",
/*  76 */        "Os",
/*  77 */        "Ir",
/*  78 */        "Pt",
/*  79 */        "Au",
/*  80 */        "Hg",
/*  81 */        "Tl",
/*  82 */        "Pb",
/*  83 */        "Bi",
/*  84 */        "Po",
/*  85 */        "At",
/*  86 */        "Rn",
/*  87 */        "Fr",
/*  88 */        "Ra",
/*  89 */        "Ac",
/*  90 */        "Th",
/*  91 */        "Pa",
/*  92 */        "U",
/*  93 */        "Np",
/*  94 */        "Pu",
/*  95 */        "Am",
/*  96 */        "Cm",
/*  97 */        "Bk",
/*  98 */        "Cf",
/*  99 */        "Es",
/* 100 */        "Fm",
/* 101 */        "Md",
/* 102 */        "No",
/* 103 */        "Lr",
/* 104 */        "Rf",
/* 105 */        "Db",
/* 106 */        "Sg",
/* 107 */        "Bh",
/* 108 */        "Hs",
/* 109 */        "Mt",
/* 110 */        "Ds",
/* 111 */        "Rg"
};

const char *elementZToSymbol(unsigned int aAtomZ)
{
	if(aAtomZ < 1 || aAtomZ > MAX_Z) return "";
	return gElementSymbols[aAtomZ];
}

unsigned int elementSymbolToZ(const char* aElementSymbol, unsigned int aDefaultX)
{
	unsigned int z;
	char  n[2];
	const char *p;

	// Transform first uppercase second lowercase
	if(aElementSymbol[0] == '\0') return aDefaultX;
	n[0] = toupper(aElementSymbol[0]);
	n[1] = (aElementSymbol[1] == '\0' || !isalpha(aElementSymbol[1])) ? '\0' : tolower(aElementSymbol[1]);

	int zblank = 0;
	for(z=1; z <= MAX_Z; ++z)
	{
		p = gElementSymbols[z];
		if(n[0] == p[0])
		{
			if(n[1] == p[1]) break;
			if(p[1] == '\0' && !zblank) zblank = z;
		}
	}

	// If not found try to check only the first letter
	// else atom symbol not found so forced to the one passed as argument
	if(z > MAX_Z)
	{
		// Special case: deuterium is converted to hydrogen
		if(n[0] == 'D') return 1;

		// Special case for lone pair (Lp) atoms in MOL
		if(n[0] == 'L' && n[1] == 'p') return 0;

		// Not found, return the default Z value
		z = (zblank) ? zblank : aDefaultX;
	}

	return z;
}


void convertAtomSymbolsToZ(const char *aAtomTypesStr, std::vector<unsigned int>& aAtomZ)
{
	std::string atom_types(aAtomTypesStr);
	aAtomZ.clear();
	if(atom_types.size() > 0)
	{
		const char* separator = " ,;\t\n";
		std::string::size_type start;
		std::string::size_type end = 0;

		std::vector<std::string> names;
		for(;;)
		{
			start = atom_types.find_first_not_of(separator, end);
			if(start == std::string::npos) break;
			end = atom_types.find_first_of(separator, start);

			names.push_back(std::string(atom_types, start, end - start));
		}

		aAtomZ.reserve(names.size());
		for(std::vector<std::string>::const_iterator ia = names.begin(); ia != names.end(); ++ia)
		{
			aAtomZ.push_back(elementSymbolToZ(ia->c_str()));
		}
	}
}
