from re import findall

# Constants
AA_MASS = dict(
    G=57.05132, A=71.0779, S=87.0773, P=97.11518, V=99.13106,
    T=101.10388, C=103.1429, L=113.15764, I=113.15764, N=114.10264,
    D=115.0874, Q=128.12922, K=128.17228, E=129.11398, M=131.19606,
    H=137.13928, F=147.17386, R=156.18568, Y=163.17326, W=186.2099
)

WATER_MASS = 18.01528
H_MASS = 1.00784


MASSES = dict(
    H=1.008, He=4.0026, Li=6.94, Be=9.0122, B=10.81,
    C=12.011, N=14.007, O=15.999, F=18.998, Ne=20.18,
    Na=22.99, Mg=24.305, Al=26.982, Si=28.085, P=30.974,
    S=32.06, Cl=35.45, Ar=39.948, K=39.098, Ca=40.078
)

def read_mass(table):
    m = {}
    with open(table, "r") as f:
        table_lines = f.readlines()
        
        for line in table_lines:
            # Don't forget to strip potential \r\n or \n from end of line and cast str to float
            key, value = line.split(",")
            m[key] = float(value.rstrip())

    return m

def calculate_mass(formula):
    # Chemical compounds form a space separated list in which each entry 
    # starts with exactly one capital letter
    # followed by at most one small letter
    # followed by 0 or more integers 

    # The regex has two match groups.
    # The first pattern looks for the capital letter small letter combination and is index 0 of the resulting tuple
    # The second pattern looks for the corresponding integer and is index 1 of the resulting tuple
    compounds = findall(r"([A-Z][a-z]?)([0-9]*)", formula)
    
    # Compute mass
    mass = 0.0
    for (element, amount) in compounds:
        # In case no amount was specified default to 1
        if amount == '':
            amount = 1
        mass += MASSES[element]*float(amount)

    return mass



class Protein:
    def __init__(self, name, uniprot_id, sequence):
        self.name = name
        self.uniprot_id = uniprot_id
        self.sequence = sequence

    def get_length(self):
        return len(self.sequence)

    def contains(self, peptide):
        return peptide in self.sequence

    def get_mw(self, disulfides = 0):
        # sum of the masses of the amino acid residues
        # + mass of a water molecule
        # - 2 * mass of hydrogen * number of disulfide bridges
        mw = WATER_MASS
        for acid in self.sequence:
            mw += AA_MASS[acid]

        mw -= 2*H_MASS*disulfides

        return mw


# Running the python file executes this section
if __name__ == "__main__":
    ## 
    galanin = Protein("Galanin", "P22466", "GWTLNSAGYLLGPHAVGNHRSFSDKNGLTS")
    print(type(galanin))
    print(galanin.get_length())
    print(galanin.contains("CGSHLV"))
    print(galanin.get_mw())


    insulin_B = Protein("Insulin B chain", "P01308", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    print(insulin_B.get_length())
    print(insulin_B.contains("CGSHLV"))
    print(insulin_B.get_mw(disulfides=1))
    
    
    ##
    m = read_mass("average_mass.csv")
    print(2 * m["H"] + m["He"])

    ## 
    print(calculate_mass("C6 H12 O6"))
    print(calculate_mass("H2 O"))
    print(calculate_mass("C34 H46 Cl N3 O10"))