# References

Every experimental and methodological source used in the project, with DOIs, for building the
paper's reference list. All DOIs were checked against Crossref on 2026-11-09. When a card gains a new
source, add it here as well.

**Status column.** "To confirm" means the value on the card has not yet been checked against the
primary paper by J. Picione. While that is so, the card says "Value to be confirmed against the
primary paper by J. Picione" (`experimental_confirmed = false` in `00_inputs/card_notes.toml`).

## Experimental values on the cards

### SN2: HO⁻ + CH₃Br, water

| Source | Used for | Status |
|---|---|---|
| W. Mabey and T. Mill, "Critical review of hydrolysis of organic compounds in water under environmental conditions", *J. Phys. Chem. Ref. Data* **1978**, *7*, 383–415. doi:[10.1063/1.555572](https://doi.org/10.1063/1.555572). Free reprint at srd.nist.gov. | Table 4.2: k(B) for MeBr at 373 K = 0.352 M⁻¹ s⁻¹, Ea = 96.3 kJ/mol (23.0 kcal/mol), log A = 13.017. Table 5.1: k(B)[OH⁻] at pH 7 and 298 K = 1.41 × 10⁻¹¹ s⁻¹, so k(B) = 1.4 × 10⁻⁴ M⁻¹ s⁻¹. Eyring then gives ΔG‡ ≈ 95.0 kJ/mol (22.7 kcal/mol). | To confirm |
| I. Fells and E. A. Moelwyn-Hughes, "The kinetics of the hydrolysis of the chlorinated methanes", *J. Chem. Soc.* **1959**, 398. doi:[10.1039/JR9590000398](https://doi.org/10.1039/JR9590000398) | The primary source the review gives for Table 4.2 (its ref. 7, in the references for Section 4.1). | To confirm |
| R. H. Bathgate and E. A. Moelwyn-Hughes, "The kinetics of certain ionic exchange reactions of the four methyl halides in aqueous solution", *J. Chem. Soc.* **1959**, 2642. doi:[10.1039/JR9590002642](https://doi.org/10.1039/JR9590002642) | Possible true origin of the CH₃Br + OH⁻ coefficients. Ref. 7's title covers only the chlorinated methanes, yet Table 4.2 lists all four methyl halides under it. | To check |
| E. A. Moelwyn-Hughes, "The kinetics of certain reactions between methyl halides and anions in water", *Proc. R. Soc. Lond. A* **1949**, *196*, 540–553. doi:[10.1098/rspa.1949.0044](https://doi.org/10.1098/rspa.1949.0044) | Earlier measurements of the same reactions; a second candidate. | To check |

### E2: EtO⁻ + 2-bromopropane, ethanol

| Source | Used for | Status |
|---|---|---|
| M. L. Dhar, E. D. Hughes and C. K. Ingold, "Mechanism of elimination reactions. Part X. Kinetics of olefin elimination from isopropyl, sec.-butyl, 2-n-amyl, and 3-n-amyl bromides in acidic and alkaline alcoholic media", *J. Chem. Soc.* **1948**, 2058–2065. doi:[10.1039/JR9480002058](https://doi.org/10.1039/JR9480002058) | Tables III and VII: overall k ≈ 2.88 × 10⁻⁶ M⁻¹ s⁻¹ at 25 °C (19 points, 2.76–2.96 × 10⁻⁶ M⁻¹ s⁻¹) and 80.3 % propene. The E2 share then gives ΔG‡(E2) ≈ 105.2 kJ/mol (25.1 kcal/mol) and ΔG‡(SN2) ≈ 108.7 kJ/mol (26.0 kcal/mol). | To confirm |
| E. D. Hughes, C. K. Ingold, S. Masterman and B. J. McNulty, "Mechanism of elimination reactions. Part V. …", *J. Chem. Soc.* **1940**, 899–912. doi:[10.1039/JR9400000899](https://doi.org/10.1039/JR9400000899) | Companion paper: earlier isopropyl bromide kinetics in alkaline ethanol. | Companion |
| M. L. Dhar, E. D. Hughes, C. K. Ingold, A. M. M. Mandour, G. A. Maw and L. I. Woolf, "Mechanism of elimination reactions. Part XVI. Constitutional influences in elimination. A general discussion", *J. Chem. Soc.* **1948**, 2093. doi:[10.1039/JR9480002093](https://doi.org/10.1039/JR9480002093) | Companion paper: general discussion. | Companion |

### Diels–Alder: butadiene + ethylene, gas phase

| Source | Used for | Status |
|---|---|---|
| D. Rowley and H. Steiner, "Kinetics of diene reactions at high temperatures", *Discuss. Faraday Soc.* **1951**, *10*, 198–213. doi:[10.1039/DF9511000198](https://doi.org/10.1039/DF9511000198), as compiled in the NIST Chemical Kinetics Database (NIST Standard Reference Database 17, [kinetics.nist.gov](https://kinetics.nist.gov/)), squib 1951ROW/STE198-213. | C₂H₄ + 1,3-butadiene → cyclohexene, gas phase, second order, 760–921 K; Ea = 115 kJ/mol, A = 4.98 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹ (3.0 × 10⁷ L mol⁻¹ s⁻¹), n = 0. Pressure 3.72 × 10⁴ to 1.01 × 10⁵ Pa; rate constant reported as pressure independent; static thermal experiment followed by pressure measurement; absolute values. NIST lists k(825 K) = 2.59 × 10⁻²¹ and k(850 K) = 4.23 × 10⁻²¹ cm³ molecule⁻¹ s⁻¹ (A and Ea above give 2.61 × 10⁻²¹ and 4.27 × 10⁻²¹, within 1 %, rounding of the listed parameters). Used for the gas-phase Ea comparison at 840 K (midpoint of the range) and for k(840 K) = 2.12 L mol⁻¹ s⁻¹ (3.52 × 10⁻²¹ cm³ molecule⁻¹ s⁻¹) against transition-state theory. | Confirmed by J. Picione from the NIST entry, 2026-11-10 (Ea and temperature range); original paper not yet read |
| K. N. Houk, Y.-T. Lin and F. K. Brown, "Evidence for the concerted mechanism of the Diels–Alder reaction of butadiene with ethylene", *J. Am. Chem. Soc.* **1986**, *108*, 554–556. doi:[10.1021/ja00263a059](https://doi.org/10.1021/ja00263a059) | Concerted transition state. | Comparison |
| E. Goldstein, B. Beno and K. N. Houk, "Density functional theory prediction of the relative energies and isotope effects for the concerted and stepwise mechanisms of the Diels–Alder reaction of butadiene and ethylene", *J. Am. Chem. Soc.* **1996**, *118*, 6036–6043. doi:[10.1021/ja9601494](https://doi.org/10.1021/ja9601494) | Concerted transition state. | Comparison |

### Nucleophilic addition: cyanide + acetone, water

Leads supplied by J. Picione (2026-11-10); the values on the card were read by J. Picione from
Guthrie 1998, Table 1 (2026-11-11). Water, 25 °C. The measured rate refers to free cyanide (CN⁻, not
HCN), which is the species that adds to the carbonyl, so observed rates depend on pH.

| Source | Used for | Status |
|---|---|---|
| J. P. Guthrie, "Prediction of rate constants for cyanohydrin formation using equilibrium constants and distortion energies", *J. Am. Chem. Soc.* **1998**, *120*, 1688–1694. doi:[10.1021/ja970506q](https://doi.org/10.1021/ja970506q) | Compilation (Table 1) for acetone: log K = −1.84 for K = [R₂C(O⁻)CN]/([R₂C=O][CN⁻]), so K = 0.0145 M⁻¹ and ΔG° = +10.5 kJ/mol (+2.5 kcal/mol); rate = k[CN⁻][R₂C=O], log k = 0.60, so k = 4.0 M⁻¹ s⁻¹ and ΔG‡ = 69.6 kJ/mol (16.6 kcal/mol) (298.15 K, 1 M standard state, Eyring with no transmission correction; conversions checked by Claude). Both values on the card, cited as the compilation. | Confirmed by J. Picione from Table 1, 2026-11-11 |
| G. Schlesinger and S. L. Miller, "Equilibrium and kinetics of glyconitrile formation in aqueous solution", *J. Am. Chem. Soc.* **1973**, *95*, 3729–3735. doi:[10.1021/ja00792a043](https://doi.org/10.1021/ja00792a043) | Original source of the measurements as cited by Guthrie 1998; cited on the card as such. Not yet read; its title names glyconitrile (the formaldehyde cyanohydrin), so whether the acetone values come from this paper is still to be checked. | Not yet read |
| W. J. Svirbely and J. F. Roth, "Carbonyl reactions. I. The kinetics of cyanohydrin formation in aqueous solution", *J. Am. Chem. Soc.* **1953**, *75*, 3106–3111. doi:[10.1021/ja01109a020](https://doi.org/10.1021/ja01109a020) | Lead: measured kinetics in aqueous solution. | Lead, to be confirmed by J. Picione |

### Carbocation 1,2-hydride shift: 2,3-dimethyl-2-butyl cation (degenerate), low-temperature NMR

Leads found by Claude (2026-11-10) at J. Picione's request; the value on the card was read from
the paper by J. Picione (2026-11-11). The measurements are in superacid solution, which SMD has no
parameters for; the calculations use SMD dichloromethane as a polar, non-nucleophilic stand-in, with
gas phase alongside. The comparison is made at 135 K.

| Source | Used for | Status |
|---|---|---|
| M. Saunders and M. R. Kates, "Rates of degenerate 1,2-hydride and 1,2-methide shifts from the carbon-13 nuclear magnetic resonance spectra of tertiary alkyl cations", *J. Am. Chem. Soc.* **1978**, *100*, 7082–7083. doi:[10.1021/ja00490a060](https://doi.org/10.1021/ja00490a060) | Compound I, degenerate 1,2-hydride shift in the 2,3-dimethyl-2-butyl cation: ΔG‡ = 13.0 ± 0.4 kJ/mol (3.1 ± 0.1 kcal/mol) at −138 °C (135 K); k = 3.1 × 10⁷ s⁻¹ (135 K), 4.6 × 10⁷ s⁻¹ (145 K), 1.4 × 10⁸ s⁻¹ (154 K). 50:50 SO₂ClF/SO₂F₂; ¹³C NMR line broadening at 67.9 MHz, fast-exchange limit; temperatures ±2 °C. No Arrhenius fit, so no measured Ea. The authors infer weaker hyperconjugation at the transition state than in the cation. Used for the ΔG‡ and Eyring rate comparison at 135 K. | Confirmed by J. Picione from the original paper, 2026-11-11 |
| M. Saunders and H. A. Jiménez-Vázquez, "Recent studies of carbocations", *Chem. Rev.* **1991**, *91*, 375–397. doi:[10.1021/cr00003a005](https://doi.org/10.1021/cr00003a005) | Lead: review that collects the hydride-shift barriers; for cross-checking the value, not for citing it. | Lead, to be confirmed by J. Picione |

### Secondary source consulted (not cited on the cards)

| Source | Note |
|---|---|
| C. G. Jesudason, "I. Determination of chemical reaction rate constants by numerical nonlinear analysis: differential methods", arXiv:1101.5060 (2011). | The E2 numbers were first read from this paper's reproduction of Dhar, Hughes and Ingold's tables. The card now cites the primary paper directly. Not for the paper's reference list. |

## Methods

### Program

| Source | Used for |
|---|---|
| F. Neese, "Software update: the ORCA program system—version 6.0", *WIREs Comput. Mol. Sci.* **2025**, *15*, e70019. doi:[10.1002/wcms.70019](https://doi.org/10.1002/wcms.70019) | ORCA 6.1.1, all calculations |
| F. Neese, F. Wennmohs, U. Becker and C. Riplinger, "The ORCA quantum chemistry program package", *J. Chem. Phys.* **2020**, *152*, 224108. doi:[10.1063/5.0004608](https://doi.org/10.1063/5.0004608) | ORCA |
| F. Neese, "The SHARK integral generation and digestion system", *J. Comput. Chem.* **2023**, *44*, 381–396. doi:[10.1002/jcc.26942](https://doi.org/10.1002/jcc.26942) | Integrals (ORCA lists it for every job) |

### Density functionals and dispersion

| Source | Used for |
|---|---|
| A. D. Becke, "Density-functional thermochemistry. III. The role of exact exchange", *J. Chem. Phys.* **1993**, *98*, 5648–5652. doi:[10.1063/1.464913](https://doi.org/10.1063/1.464913) | B3LYP |
| C. Lee, W. Yang and R. G. Parr, "Development of the Colle–Salvetti correlation-energy formula into a functional of the electron density", *Phys. Rev. B* **1988**, *37*, 785–789. doi:[10.1103/PhysRevB.37.785](https://doi.org/10.1103/PhysRevB.37.785) | B3LYP |
| P. J. Stephens, F. J. Devlin, C. F. Chabalowski and M. J. Frisch, "Ab initio calculation of vibrational absorption and circular dichroism spectra using density functional force fields", *J. Phys. Chem.* **1994**, *98*, 11623–11627. doi:[10.1021/j100096a001](https://doi.org/10.1021/j100096a001) | B3LYP |
| Y. Zhao and D. G. Truhlar, "The M06 suite of density functionals …", *Theor. Chem. Acc.* **2008**, *120*, 215–241. doi:[10.1007/s00214-007-0310-x](https://doi.org/10.1007/s00214-007-0310-x) | M06-2X comparison curve |
| S. Grimme, J. Antony, S. Ehrlich and H. Krieg, "A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H–Pu", *J. Chem. Phys.* **2010**, *132*, 154104. doi:[10.1063/1.3382344](https://doi.org/10.1063/1.3382344) | D3 (B3LYP and M06-2X) |
| S. Grimme, S. Ehrlich and L. Goerigk, "Effect of the damping function in dispersion corrected density functional theory", *J. Comput. Chem.* **2011**, *32*, 1456–1465. doi:[10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759) | Becke–Johnson damping, D3(BJ) |

### Basis sets and integral approximations

| Source | Used for |
|---|---|
| F. Weigend and R. Ahlrichs, "Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence quality for H to Rn", *Phys. Chem. Chem. Phys.* **2005**, *7*, 3297–3305. doi:[10.1039/b508541a](https://doi.org/10.1039/b508541a) | def2-SVP, def2-TZVP |
| J. Zheng, X. Xu and D. G. Truhlar, "Minimally augmented Karlsruhe basis sets", *Theor. Chem. Acc.* **2011**, *128*, 295–305. doi:[10.1007/s00214-010-0846-z](https://doi.org/10.1007/s00214-010-0846-z) | ma-def2-SVP, ma-def2-TZVP |
| F. Weigend, "Accurate Coulomb-fitting basis sets for H to Rn", *Phys. Chem. Chem. Phys.* **2006**, *8*, 1057–1065. doi:[10.1039/b515623h](https://doi.org/10.1039/b515623h) | def2/J auxiliary basis |
| G. L. Stoychev, A. A. Auer and F. Neese, "Automatic generation of auxiliary basis sets", *J. Chem. Theory Comput.* **2017**, *13*, 554–562. doi:[10.1021/acs.jctc.6b01041](https://doi.org/10.1021/acs.jctc.6b01041) | AutoAux (DLPNO-CCSD(T)) |
| F. Neese, "An improvement of the resolution of the identity approximation for the formation of the Coulomb matrix", *J. Comput. Chem.* **2003**, *24*, 1740–1747. doi:[10.1002/jcc.10318](https://doi.org/10.1002/jcc.10318) | RI-J |
| F. Neese, F. Wennmohs, A. Hansen and U. Becker, "Efficient, approximate and parallel Hartree–Fock and hybrid DFT calculations. A 'chain-of-spheres' algorithm for the Hartree–Fock exchange", *Chem. Phys.* **2009**, *356*, 98–109. doi:[10.1016/j.chemphys.2008.10.036](https://doi.org/10.1016/j.chemphys.2008.10.036) | COSX |
| B. Helmich-Paris, B. de Souza, F. Neese and R. Izsák, "An improved chain of spheres for exchange algorithm", *J. Chem. Phys.* **2021**, *155*, 104109. doi:[10.1063/5.0058766](https://doi.org/10.1063/5.0058766) | COSX (ORCA 6 grids) |

### Solvation

| Source | Used for |
|---|---|
| V. Barone and M. Cossi, "Quantum calculation of molecular energies and energy gradients in solution by a conductor solvent model", *J. Phys. Chem. A* **1998**, *102*, 1995–2001. doi:[10.1021/jp9716997](https://doi.org/10.1021/jp9716997) | C-PCM |
| A. V. Marenich, C. J. Cramer and D. G. Truhlar, "Universal solvation model based on solute electron density and on a continuum model of the solvent defined by the bulk dielectric constant and atomic surface tensions", *J. Phys. Chem. B* **2009**, *113*, 6378–6396. doi:[10.1021/jp810292n](https://doi.org/10.1021/jp810292n) | SMD (water, ethanol) |
| M. Garcia-Ratés and F. Neese, "Effect of the solute cavity on the solvation energy and its derivatives within the framework of the Gaussian charge scheme", *J. Comput. Chem.* **2020**, *41*, 922–939. doi:[10.1002/jcc.26139](https://doi.org/10.1002/jcc.26139) | ORCA C-PCM cavity |
| M. Garcia-Ratés and F. Neese, "Efficient implementation of the analytical second derivatives of Hartree–Fock and hybrid DFT energies within the framework of the conductor-like polarizable continuum model", *J. Comput. Chem.* **2019**, *40*, 1816–1828. doi:[10.1002/jcc.25833](https://doi.org/10.1002/jcc.25833) | Frequencies in C-PCM |
| M. Garcia-Ratés, U. Becker and F. Neese, "Implicit solvation in domain based pair natural orbital coupled cluster (DLPNO-CCSD) theory", *J. Comput. Chem.* **2021**, *42*, 1959–1973. doi:[10.1002/jcc.26726](https://doi.org/10.1002/jcc.26726) | C-PCM with DLPNO-CCSD(T) |

### Coupled cluster

| Source | Used for |
|---|---|
| C. Riplinger and F. Neese, "An efficient and near linear scaling pair natural orbital based local coupled cluster method", *J. Chem. Phys.* **2013**, *138*, 034106. doi:[10.1063/1.4773581](https://doi.org/10.1063/1.4773581) | DLPNO-CCSD |
| C. Riplinger, B. Sandhoefer, A. Hansen and F. Neese, "Natural triple excitations in local coupled cluster calculations with pair natural orbitals", *J. Chem. Phys.* **2013**, *139*, 134101. doi:[10.1063/1.4821834](https://doi.org/10.1063/1.4821834) | (T) triples |
| C. Riplinger, P. Pinski, U. Becker, E. F. Valeev and F. Neese, "Sparse maps—A systematic infrastructure for reduced-scaling electronic structure methods. II. Linear scaling domain based pair natural orbital coupled cluster theory", *J. Chem. Phys.* **2016**, *144*, 024109. doi:[10.1063/1.4939030](https://doi.org/10.1063/1.4939030) | DLPNO-CCSD(T) implementation |
| D. G. Liakos, M. Sparta, M. K. Kesharwani, J. M. L. Martin and F. Neese, "Exploring the accuracy limits of local pair natural orbital coupled-cluster theory", *J. Chem. Theory Comput.* **2015**, *11*, 1525–1539. doi:[10.1021/ct501129s](https://doi.org/10.1021/ct501129s) | TightPNO settings |

ORCA's own citation list for the DLPNO jobs also names: F. Neese, F. Wennmohs and A. Hansen, *J. Chem. Phys.* **2009**, *130*, 114108, doi:[10.1063/1.3086717](https://doi.org/10.1063/1.3086717); F. Neese, A. Hansen and D. G. Liakos, *J. Chem. Phys.* **2009**, *131*, 064103, doi:[10.1063/1.3173827](https://doi.org/10.1063/1.3173827); G. Bistoni et al., *J. Chem. Theory Comput.* **2017**, *13*, 3220–3227, doi:[10.1021/acs.jctc.7b00352](https://doi.org/10.1021/acs.jctc.7b00352); R. Izsák, A. Hansen and F. Neese, *Mol. Phys.* **2012**, *110*, 2413–2417, doi:[10.1080/00268976.2012.687466](https://doi.org/10.1080/00268976.2012.687466).

### Reaction path, thermochemistry and rates

| Source | Used for |
|---|---|
| K. Fukui, "The path of chemical reactions — the IRC approach", *Acc. Chem. Res.* **1981**, *14*, 363–368. doi:[10.1021/ar00072a001](https://doi.org/10.1021/ar00072a001) | Intrinsic reaction coordinate |
| D. Bykov, T. Petrenko, R. Izsák, S. Kossmann, U. Becker, E. Valeev and F. Neese, "Efficient implementation of the analytic second derivatives of Hartree–Fock and hybrid DFT energies: a detailed analysis of different approximations", *Mol. Phys.* **2015**, *113*, 1961–1977. doi:[10.1080/00268976.2015.1025114](https://doi.org/10.1080/00268976.2015.1025114) | Analytic Hessians (TS search, frequencies) |
| S. Grimme, "Supramolecular binding thermodynamics by dispersion-corrected density functional theory", *Chem. Eur. J.* **2012**, *18*, 9955–9964. doi:[10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497) | Quasi-RRHO vibrational entropy (ORCA default, 100 cm⁻¹) |
| H. Eyring, "The activated complex in chemical reactions", *J. Chem. Phys.* **1935**, *3*, 107–115. doi:[10.1063/1.1749604](https://doi.org/10.1063/1.1749604) | ΔG‡ from measured rate constants and the transition-state-theory rate constant (Diels–Alder and hydride-shift cards) |
| E. Wigner, "Über das Überschreiten von Potentialschwellen bei chemischen Reaktionen", *Z. Phys. Chem.* **1932**, *19B*, 203–216. doi:[10.1515/zpch-1932-1920](https://doi.org/10.1515/zpch-1932-1920) | Wigner tunnelling factor, given as a rough indicator only (hydride-shift card) |
| A. Fernández-Ramos, B. A. Ellingson, R. Meana-Pañeda, J. M. C. Marques and D. G. Truhlar, "Symmetry numbers and chemical reaction rates", *Theor. Chem. Acc.* **2007**, *118*, 813–826. doi:[10.1007/s00214-007-0328-0](https://doi.org/10.1007/s00214-007-0328-0) | Rotational symmetry numbers in transition-state theory, including symmetric reactions and chiral reactants and transition states: the treatment used for the degenerate hydride shift (σ = 1 for the C1 cation; the C2 axis of the TS reverses the reaction coordinate and is not counted; the factor 2 for mirror-image conformers cancels). Checked on Crossref; only the abstract read so far. |

### Orbitals and tracking

| Source | Used for |
|---|---|
| G. Knizia, "Intrinsic atomic orbitals: an unbiased bridge between quantum theory and chemical concepts", *J. Chem. Theory Comput.* **2013**, *9*, 4834–4843. doi:[10.1021/ct400687b](https://doi.org/10.1021/ct400687b) | IAO/IBO localization |
| R. S. Mulliken, "Electronic population analysis on LCAO–MO molecular wave functions. I", *J. Chem. Phys.* **1955**, *23*, 1833–1840. doi:[10.1063/1.1740588](https://doi.org/10.1063/1.1740588) | Orbital labels from atomic populations |
| H. W. Kuhn, "The Hungarian method for the assignment problem", *Nav. Res. Logist. Q.* **1955**, *2*, 83–97. doi:[10.1002/nav.3800020109](https://doi.org/10.1002/nav.3800020109) | Frame-to-frame orbital matching |

### Starting structures

| Source | Used for |
|---|---|
| S. Riniker and G. A. Landrum, "Better informed distance geometry: using what we know to improve conformation generation", *J. Chem. Inf. Model.* **2015**, *55*, 2562–2574. doi:[10.1021/acs.jcim.5b00654](https://doi.org/10.1021/acs.jcim.5b00654) | RDKit ETKDG |
| T. A. Halgren, "Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94", *J. Comput. Chem.* **1996**, *17*, 490–519. doi:[10.1002/(SICI)1096-987X(199604)17:5/6<490::AID-JCC1>3.0.CO;2-P](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6%3C490::AID-JCC1%3E3.0.CO;2-P) | MMFF pre-optimization |

### Graphics and software

| Source | Used for |
|---|---|
| W. E. Lorensen and H. E. Cline, "Marching cubes: a high resolution 3D surface construction algorithm", *Comput. Graph.* **1987**, *21*, 163–169. doi:[10.1145/37402.37422](https://doi.org/10.1145/37402.37422) | Isosurfaces |
| T. Lewiner, H. Lopes, A. W. Vieira and G. Tavares, "Efficient implementation of marching cubes' cases with topological guarantees", *J. Graph. Tools* **2003**, *8*, 1–15. doi:[10.1080/10867651.2003.10487582](https://doi.org/10.1080/10867651.2003.10487582) | Isosurfaces (scikit-image algorithm) |
| M. Garland and P. S. Heckbert, "Surface simplification using quadric error metrics", *Proc. SIGGRAPH '97* **1997**, 209–216. doi:[10.1145/258734.258849](https://doi.org/10.1145/258734.258849) | Mesh decimation (fast-simplification 0.1.13) |
| S. van der Walt et al., "scikit-image: image processing in Python", *PeerJ* **2014**, *2*, e453. doi:[10.7717/peerj.453](https://doi.org/10.7717/peerj.453) | scikit-image |
| P. Virtanen et al., "SciPy 1.0: fundamental algorithms for scientific computing in Python", *Nat. Methods* **2020**, *17*, 261–272. doi:[10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) | Hungarian assignment, cubic spline between coupled-cluster frames |
| C. R. Harris et al., "Array programming with NumPy", *Nature* **2020**, *585*, 357–362. doi:[10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) | NumPy |
| J. D. Hunter, "Matplotlib: a 2D graphics environment", *Comput. Sci. Eng.* **2007**, *9*, 90–95. doi:[10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55) | Energy profiles, video inset |

Software without a paper DOI (cite by name, version and URL): RDKit 2026.09.1 (rdkit.org);
three.js r160 (threejs.org); Blender 5.2 (blender.org); FFmpeg 9.0.2 (ffmpeg.org); Pillow
(python-pillow.org).
